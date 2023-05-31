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

#' @title Simulate a time series of exponential waiting times with a change point
#' @param t1 no. of bins after the change point
#' @param r1 event rate after the change point
#' @param t0 no. of bins before the change point
#' @param r0 event rate before the change point
#' @param ... arguments passed to \code{series_sim}, typically the probability-generating
#' function \code{pfun}
#' @export
twoseries_sim <- function(t1, r1, t0, r0, ...) {

  s1 <- series_sim(t1, r1, ...)
  s0 <- series_sim(t0, r0, ...)

  return(c(s1, s0))

}