#' @rdname msd
#' @title Detect the most significant deviation of a fuzzy event rate in a time series
#' @description
#' The \code{msd} routine tests a series of candidate change points for a deviation
#' (increase or decrease) of a fuzzy event rate, identifies the most significant
#' deviation in the series, then performs a permutation test to estimate the
#' overall significance of the deviation.#'
#' @param s time series of fuzzy event counts, ordered backward in time (more recent bins come first)
#' @param tmax number of candidate change points, defaults to half the series length
#' @param alternative either \code{greater} (default) or \code{less}, the test direction
#' @param n_perm no. of permutations to estimate the null distribution
#' @param alpha_perm error rate for the confidence interval of the P-value
#' @param ... additional arguments passed to \code{pf()}, notably \code{log} and \code{lower.tail}
#' @return Subroutine \code{fpr.test} returns a vector of P-values of the same length as \code{s}.
#'
#' Subroutine \code{msd} returns a numeric vector with named elements
#' \item{p.value}{approximate P-value estimated using a Gamma fit (see Details)}
#' \item{p.lo}{lower bound of the P-value confidence interval}
#' \item{p.hi}{upper bound of the P-value confidence interval}
#' \item{tau}{the estimated change point, in units of time before present}
#' \item{rate0}{the estimated fuzzy rate before the change point}
#' \item{rate1}{the estimated fuzzy rate after the change point}
#' @export
fpr.test <- function(s, tmax, ...) {

  n <- length(s)
  X <- sum(s)

  # Change points
  t  <- 1:tmax
  # Count before t (present)
  x1 <- cumsum(s[t])
  # Count after t (past)
  x0 <- X - x1
  # Test statistic
  ratio <- (x1 + 1) * (n - t) / ((x0 + 1) * t)
  pval  <- pf(ratio, 2*(x0 + 1), 2*(x1 + 1), ...)

  return(pval)
}

#' @rdname msd
#' @export
msd <- function(s, tmax, alternative = "greater", n_perm = 1e3, alpha_perm = 0.05) {

  n <- length(s)
  if(missing(tmax)) tmax <- ceiling(n / 2)
  X <- sum(s)

  # Return lower or upper tail of cumulative distribution depending on H0
  if(alternative == "greater") {
    # One-sided log P-value (H0: r1 <= r0)
    lower <- FALSE
  } else if(alternative == "less") {
    # One-sided log P-value (H0: r1 >= r0)
    lower <- TRUE
  } else stop("Unimplemented alternative.")


  pval_log <- fpr.test(s, tmax, log = TRUE, lower.tail = lower)

  t_msd <- which.min(pval_log)

  # Estimate null distribution of min P-value under permutations
  perm_series <- replicate(n_perm, {
    min(fpr.test(sample(s), tmax, log = TRUE, lower.tail = lower))
  })

  # Permutation P-value with Clopper-Pearson interval
  n_success <- sum(perm_series <= pval_log[t_msd])
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
    pval_gamma <- pgamma(-pval_log[t_msd], gamma_shape, gamma_rate, lower.tail = FALSE)
  } else {
    pval_gamma <- 1.0
  }

  # Consistency: coerce pval_gamma within Clopper-Pearson interval
  pval_gamma <- min(pval_ci_hi, pval_gamma)
  pval_gamma <- max(pval_ci_lo, pval_gamma)

  return(c(
    p.value = pval_gamma,
    p.lo = pval_ci_lo,
    p.center = pval_perm,
    p.hi = pval_ci_hi,
    tau = t_msd,
    rate0 = sum(tail(s, n - t_msd)) / (n - t_msd),
    rate1 = sum(head(s, t_msd)) / t_msd
  ))
}




msd_DEPREC <- function(s, tmax, alternative = "greater", n_perm = 1e3, alpha_perm = 0.05) {

  if(missing(tmax)) tmax <- ceiling(length(s) / 2)
  n <- length(s)
  X <- sum(s)
  n <- length(s)

  # Change points
  t  <- 1:tmax
  # Count before t (present)
  x1 <- cumsum(s[t])
  # Count after t (past)
  x0 <- X - x1

  # Test statistic
  ratio <- (x1 + 1) * (n - t) / ((x0 + 1) * t)
  # One-sided log P-value (H0: r1 <= r0)
  pval_log <- pf(ratio, 2*(x0 + 1), 2*(x1 + 1), log = TRUE, lower.tail = FALSE)
  # One-sided log P-value (H0: r1 >= r0)
  pval_log_alt <- pf(1/ratio, 2*(x0 + 1), 2*(x1 + 1), log = TRUE, lower.tail = FALSE)
  # Two-sided log P-value
  pval_log_2sided <- log(2) + pmin(pval_log, pval_log_alt)

  if(alternative == "greater") {
    # Do nothing, keep pval_log
  } else if(alternative == "less") {
    pval_log <- pval_log_alt
  } else if(alternative == "two.sided") {
    pval_log <- pval_log_2sided
  } else stop("Unimplemented alternative.")

  # Most significant deviation index
  # p_min <- min(pval)
  t_msd <- which.min(pval_log)
  ratio_msd <- ratio[t_msd]
  rate_before <- x0[t_msd] / (n - t_msd)
  rate_after  <- x1[t_msd] / t_msd

  perm_series <- replicate(n_perm, {
    s_perm <- sample(s)
    # Count before t (present)
    x1_perm <- cumsum(s_perm[t])
    # Count after t (past)
    x0_perm <- X - x1_perm

    # Test statistic
    ratio_perm <- (x1_perm + 1) * (n - t) / ((x0_perm + 1) * t)

    # One-sided log P-value (H0: r1 <= r0)
    pval_log_perm <- pf(ratio_perm, 2*(x0_perm + 1), 2*(x1_perm + 1), log = TRUE, lower.tail = FALSE)
    # One-sided log P-value (H0: r1 >= r0)
    pval_log_alt_perm <- pf(1/ratio_perm, 2*(x0_perm + 1), 2*(x1_perm + 1), log = TRUE, lower.tail = FALSE)
    # Two-sided log P-value
    pval_log_2sided_perm <- log(2) + pmin(pval_log_perm, pval_log_alt_perm)

    if(alternative == "greater") {
      # Do nothing, keep pval_log_perm
    } else if(alternative == "less") {
      pval_log_perm <- pval_log_alt_perm
    } else if(alternative == "two.sided") {
      pval_log_perm <- pval_log_2sided_perm
    } else stop("Unimplemented alternative.")

    min(pval_log_perm)
  })

  # Permutation P-value with Clopper-Pearson interval
  n_success <- sum(perm_series <= pval_log[t_msd])
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

  pval_gamma <- pgamma(-pval_log[t_msd], gamma_shape, gamma_rate, lower.tail = FALSE)

  # DEBUG
  # hist(pl, freq = FALSE); curve(dgamma(x, gamma_shape, gamma_rate), lwd = 2, add = TRUE)

  return(c(
    pval_lo = pval_ci_lo,
    pval_center = pval_perm,
    pval_hi = pval_ci_hi,
    t = t_msd,
    tmax = tmax,
    pval_log_raw = pval_log[t_msd],
    pval_gamma = pval_gamma,
    stat = ratio_msd,
    rate_before = rate_before,
    rate_after = rate_after
  ))
}