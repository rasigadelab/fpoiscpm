#' @rdname msd
#' @title Detect the most significant deviation of a fuzzy event rate in a time series
#' @description
#' The \code{msd} routine tests a series of candidate change points for a deviation
#' (increase or decrease) of a fuzzy event rate, identifies the most significant
#' deviation in the series, then performs a permutation test to estimate the
#' overall significance of the deviation.#'
#' @param s time series of fuzzy event counts, ordered backward in time (more recent bins come first)
#' @param tmax number of candidate change points, defaults to half the series length
#' @param alternative the test direction ("greater" detects an increase of rate toward present time)
#' @param n_perm no. of permutations to estimate the null distribution in P-value estimation
#' @param alpha_perm error rate for the confidence interval of the P-value
#' @param log.p return the logarithm of P-values
#' @param threshold significance threshold (type 1 error rate) for sequential change point detection
#' @return 
#' Subroutine \code{fpr.test} returns a vector of P-values of the same length as \code{s}.
#'
#' Subroutine \code{msd} returns a numeric vector with named elements
#' \item{p.value}{approximate P-value estimated using a Gamma fit (see Details)}
#' \item{p.lo}{lower bound of the P-value confidence interval}
#' \item{p.hi}{upper bound of the P-value confidence interval}
#' \item{tau}{the estimated change point, in units of time before present}
#' \item{rate0}{the estimated fuzzy rate before the change point}
#' \item{rate1}{the estimated fuzzy rate after the change point}
#'
#' Subroutine \code{msd.seq} returns a list of detected change points (whose P-value is below \code{threashold})
#' and behaves as if the time series was discovered sequentially. When a significant change point is 
#' found, the data before this change point is discarded and the analysis resumes at the last change point.
#'@examples
#'\dontrun{
#' s <- c(runif(10) + 0, runif(10) +  1, runif(10) + 0)
#' plot(rev(s))
#' print(msd.seq(s, 0.05))
#'}
#' @export
fpr.test <- function(s, tmax = floor(length(s) / 2), alternative = c("two.sided", "less", "greater"), log.p = FALSE) {

  n <- length(s)
  # if(missing(tmax)) tmax <- floor(n/2)
  X <- sum(s)

  # Change points
  t  <- 1:tmax
  # Count before t (present)
  x1 <- cumsum(s[t])
  # Count after t (past)
  x0 <- X - x1
  # Test statistic
  ratio <- (x1 + 1) * (n - t) / ((x0 + 1) * t)

  alternative <- match.arg(alternative)
  pval <- NULL

  if(alternative == "greater") {
    pval  <- pf(ratio, 2*(x0 + 1), 2*(x1 + 1), lower.tail = FALSE, log.p = log.p)
  }
  else if(alternative == "less") {
    pval  <- pf(ratio, 2*(x0 + 1), 2*(x1 + 1), lower.tail = TRUE, log.p = log.p)
  }
  else if(alternative == "two.sided") {
    pval <- pmin(
      pf(ratio, 2*(x0 + 1), 2*(x1 + 1), lower.tail = FALSE, log.p = log.p),
      pf(ratio, 2*(x0 + 1), 2*(x1 + 1), lower.tail = TRUE, log.p = log.p)
    )
    if(log.p == TRUE) {
      pval <- pval + log(2)
    } else {
      pval <- pval * 2
    }
  }

  return(pval)
}

#' @rdname msd
#' @export
msd <- function(s, n_perm = 1e3, alpha_perm = 0.05, ...) {

  n <- length(s)


  pval_log <- fpr.test(s, log.p = TRUE, ...)
  change.point <- which.min(pval_log)

  # Estimate null distribution of min P-value under permutations
  # REMARK: replicate() not used because of the ellipsis
  perm_series <- sapply(1:n_perm, function(i, series, ...) {
    min(fpr.test(sample(series), log.p = TRUE, ...))
    }, series = s, ...)


  # Permutation P-value with Clopper-Pearson interval
  n_success <- sum(perm_series <= pval_log[change.point])
  pval_perm <- n_success / n_perm
  pval_ci_lo <- qbeta(alpha_perm / 2, n_success, n_perm - n_success + 1)
  pval_ci_hi <- qbeta(1 - alpha_perm / 2, n_success + 1, n_perm - n_success)

  # Permutation P-value with Gamma approximation, method of moments
  # Assume that -log(P) has Gamma distribution
  pl <- -perm_series
  pm <- mean(pl)
  pv <- var(pl)
  gamma_shape <- pm * pm / pv
  gamma_rate  <- pm / pv

  if(is.finite(gamma_shape + gamma_rate)) {
    pval_gamma <- pgamma(-pval_log[change.point], gamma_shape, gamma_rate, lower.tail = FALSE)
  } else {
    pval_gamma <- 1.0
  }

  # Consistency: coerce pval_gamma within Clopper-Pearson interval
  pval_gamma <- min(pval_ci_hi, pval_gamma)
  pval_gamma <- max(pval_ci_lo, pval_gamma)

  return(c(
    p.value = pval_gamma,
    p.lo = pval_ci_lo,
    p.hi = pval_ci_hi,
    change.point = change.point,
    rate0 = sum(tail(s, n - change.point)) / (n - change.point),
    rate1 = sum(head(s, change.point)) / change.point
  ))
}

#' @rdname msd
#' @export
msd.seq <- function(s, threshold = 0.05, ...) {

  change.point.list <- list()

  last.change.point <- 0

  for(i in 1:length(s)) {

    s_present <- length(s) - i + 1
    s_past    <- length(s) - last.change.point

    si <- s[ s_present:s_past ]
    r  <- msd(si)
    if(r["p.value"] < threshold) {

      change.point.list[[length(change.point.list) + 1]] <- data.table(t(c(
        index = i,
        period.start = last.change.point,
        period.end = setNames(i - r["change.point"], NULL),
        r[-4]
      )))

      last.change.point <- setNames(i - r["change.point"], NULL)
    }
  }

  return(rbindlist( change.point.list ))
}
