pfpr <- function(t0, x0, t1, x1) pf(
  (x0 + 1) * t1 / ((x1 + 1) * t0),
  2*(x1 + 1), 2*(x0 + 1)
)

# Generate a series of exponential waiting times
# series_length <- 20
# event_rate  <- 1

#' @title Simulate a time series of exponential waiting times
#' @param pfun function returning the fuzzy event probability, taking
#' the series length as parameter
#' @export
series_sim <- function(series_length, event_rate, pfun) {
  if(missing(pfun)) pfun <- function(n) rep(1, n)

  # Generate much more events than needed
  event_date <- cumsum(rexp(4 * series_length * event_rate, event_rate))
  if(length(event_date) == 0) return(numeric(0))

  stopifnot(tail(event_date, 1) > series_length)

  # Last index within time series
  end <- max(which(event_date < series_length))
  if(end < 1) return(numeric(0))

  # Limit time series
  event_date <- head(event_date, end)

  # Interval series along with event dates
  interval <- factor(findInterval(event_date, 0:series_length), levels = 1:series_length)

  # Fuzzy event value count function (or anything that must be summed)
  # p <- runif(length(interval))
  # p <- rep(1, length(event_date))
  p <- pfun(length(interval))

  event_count <- aggregate(p, list(t = interval), sum, na.rm = TRUE, drop = FALSE)$x
  event_count[is.na(event_count)] <- 0.

  stopifnot(length(event_count) == series_length)
  stopifnot(sum(event_count) == sum(p))

  return(event_count)
}

# s <- series_sim(series_length, event_rate)

# Now that we have a series, apply F-test for rate equality

# n_perm = 1e3
# alpha_perm = 0.05
#

msd <- function(s, n_perm = 1e3, alpha_perm = 0.05) {

  n <- length(s)
  X <- sum(s)
  n <- length(s)

  # Change points
  t  <- 1:n
  # Count before t (present)
  x1 <- cumsum(s)
  # Count after t (past)
  x0 <- X - x1

  # Test statistic
  ratio <- (x1 + 1) * (n - t) / ((x0 + 1) * t)
  # One-sided P-value
  pval <- 1. - pf(ratio, 2*(x0 + 1), 2*(x1 + 1))
  # Two-sided P-value
  # pval <- pf(1/ratio, 2*(x0 + 1), 2*(x1 + 1)) + 1. - pf(ratio, 2*(x0 + 1), 2*(x1 + 1))

  # Most significant deviation index
  # p_min <- min(pval)
  t_msd <- which.min(pval)
  ratio_msd <- ratio[t_msd]
  rate_before <- x0[t_msd] / (n - t_msd)
  rate_after  <- x1[t_msd] / t_msd

  perm_series <- replicate(n_perm, {
    s_perm <- sample(s)
    # Change points
    # t  <- 1:n
    # Count before t (present)
    x1_perm <- cumsum(s_perm)
    # Count after t (past)
    x0_perm <- X - x1_perm

    # Test statistic
    ratio_perm <- (x1_perm + 1) * (n - t) / ((x0_perm + 1) * t)

    # One-sided P-value
    pval_perm <- 1. - pf(ratio_perm, 2*(x0_perm + 1), 2*(x1_perm + 1))

    # Two-sided P-value
    # pval_perm <- pf(1/ratio_perm, 2*(x0_perm + 1), 2*(x1_perm + 1)) + 1. - pf(ratio_perm, 2*(x0_perm + 1), 2*(x1_perm + 1))

    min(pval_perm)
  })

  # Permutation P-value with Clopper-Pearson interval
  n_success <- sum(perm_series < min(pval))
  pval_perm <- n_success / n_perm
  pval_ci_lo <- qbeta(alpha_perm / 2, n_success, n_perm - n_success + 1)
  pval_ci_hi <- qbeta(1 - alpha_perm / 2, n_success + 1, n_perm - n_success)

  return(c(
    pval_lo = pval_ci_lo,
    pval_center = pval_perm,
    pval_hi = pval_ci_hi,
    t = t_msd,
    stat = ratio_msd,
    rate_before = rate_before,
    rate_after = rate_after
  ))
}
