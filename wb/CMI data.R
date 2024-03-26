# DEMO - EVENT DETECTION USING THE FUZZY POISSON MSD METHOD

library(fpoiscpm)
library(data.table)

rm(list = objects())

# Parameters ####################################################

# Length of learning phase, during which testing is not performed
lphase <- 60
# Sliding window length
tmax_ref <- 180
# Alternative hypothesis (greater or less): detect an increase or
# decrease of the event rate ?
# alt <- "greater"
alt <- "less"
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




