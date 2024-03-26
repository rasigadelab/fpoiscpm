# DEMO - EVENT DETECTION USING THE FUZZY POISSON MSD METHOD

library(fpoiscpm)
library(data.table)

rm(list = objects())


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

if(FALSE) {
  s <- c(runif(20) +  0, runif(20) + 10)
  s <- c(runif(20) + 10, runif(20) +  0)

  fpr.test(s, alternative = "g")
  fpr.test(s, alternative = "l")
  fpr.test(s, alternative = "t")
  fpr.test(s)
}

s <- c(runif(10) + 10, runif(10) +  0)

fpr.test(s)
# print(msd(s, alternative = "two", tmax = 15))

# stop()

### Runner version: run msd repeatedly and manage series to find change points


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


s <- c(runif(10) + 0, runif(10) +  1, runif(10) + 0)


print(msd.seq(s, 0.05))


stop()


# learning.phase.length <- 5

p.cut <- 0.05

reslist <- list()

s_active <- s

last.change.point <- 0

for(i in 1:length(s)) {

  s_present <- length(s) - i + 1
  s_past    <- length(s) - last.change.point

  si <- s[ s_present:s_past ]
  r  <- msd(si)
  if(r["p.value"] < p.cut) {
    reslist[[length(reslist) + 1]] <- c(r, i = i, change.point.new = i - r["change.point"], change.point.last = last.change.point)
    last.change.point <- i - r["change.point"]
  }
}

print(reslist)

plot(s)

stop()








### TEST ZONE ###################################################


# Parameters ####################################################

# Length of learning phase, during which testing is not performed
lphase <- 60
# Sliding window length
tmax_ref <- 180
# Alternative hypothesis (greater or less): detect an increase or
# decrease of the event rate ?
alt <- "greater"
# alt <- "less"
# Number of permutations used to estimate the test P-value
n_perm <- 200
#  Alarm threshold
warning_threshold <- 0.05

# source("R/msd.R")

# Load data #####################################################

load("../monitoring_data/nmetric_indicator_covid_dataset.Rdata")

d <- data.table(data.tr_id)
rm(data.tr_id)
setnames(d, "tr_id", "x")

d[, .N, by = .(ward, species)]

unique(d$ward)

# Select ward and species ---
ward_filt <- unique(d$ward)[11]
species_filt <- "KLEPNE"

# Time series (/!\ forward in time)
s <- d[ward == ward_filt & species == species_filt][order(date)]$x

s <- rev(s)

# Scan the time series, leave 2 months (60 days) of learning phase
n <- length(s)

res <- data.table(t(sapply((lphase+1):n, function(i) {
  if(i %% 10 == 0) cat(".")
  # Suppose we're day i. Can test at most i-lphase-1 change points
  s_i <- rev(head(s, i))
  tmax <- pmin(i - lphase - 1, tmax_ref)
  c(day = i, msd(s_i, tmax = tmax, alternative = alt, n_perm = n_perm))
})))

res

{
  par(mfrow = c(2,1))

  plot(s, type = "l", xlim = c(0, n), ylab = "fuzzy event count", xlab = "time (days)")
  plot(res$day, (res$p.value), xlim = c(0, n), col = (1 + (res$p.value < warning_threshold)),
       log = "y", type = "b", pch = 19, cex = 0.5,
       ylab = "P-value", xlab = "time (days)")
  abline(h=warning_threshold, lty = 2, col = "gray")
}




