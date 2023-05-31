##################################################################
# FUZZY POISSON CHANGE POINT MODEL

library(parallel)
library(data.table)
# library(rateratio.test)

if(!exists("cl")) cl <- makeCluster(64)

rm(list = setdiff(objects(), "cl"))

# source("wb/FUNC.R")

# Generate two series



n_perm = 1e3
alpha_perm = 0.05



n  <- 100
t1 <- 10
t0 <- n - t1

r0 <- 1
r1 <- r0 * 2
tmax <- 20

# s <- twoseries_sim(t1, r1, t0, r0)
# msd(s)

# data.table(t1, r1, t1, r0, t(msd(s, tmax, n_perm = 100)))





# stop()

repfun <- function(r1_ratio = 1, nrep = 100) {
  rbindlist(replicate(nrep, {
    r1 <- r0 * r1_ratio
    s <- twoseries_sim(t1, r1, t0, r0)
    data.table(t1, r1, t1, r0, t(msd(s, tmax, n_perm = 100)))
  }, simplify = FALSE) )
}

QUIET <- clusterEvalQ(cl, {
  rm(list = setdiff(objects(), ""))
  library(data.table)
  library(fpoiscpm)
})

QUIET <- clusterExport(cl, setdiff(objects(), "cl"))

res <- rbindlist(parSapply(cl, seq(1, 4, by = 0.25), repfun, nrep = 100, simplify = FALSE))

xx <- res[ , .(rejection = mean(pval_center < 0.05)), by = r1]

# jpeg("power curve t10.jpeg", 3000, 4000, pointsize = 96, quality = 95)
{
  par(mfrow = c(2,1))

  plot(xx, xlab = "Rate fold-change after change point",
       ylab = "Rejection probability", type = "b", col = "blue", lwd = 2)
  hist(res[r1 > 2]$t, breaks = n, col = "lightgreen", freq = FALSE,
       xlab = "Change point", main = "Change point estimation")
  hist(res[r1 <= 2]$t, breaks = n, add = TRUE, freq = FALSE, col = rgb(0.8,0,0,0.3))
  legend("topleft", bty = "n", fill = c("lightgreen", rgb(0.8,0,0,0.3)),
         title = "Rate fold-change", legend = c("<= 2", "> 2"))
}
# dev.off()

stop()


# Generate a series of exponential waiting times
series_length <- 40
event_rate  <- 10
n_perm = 1e2
alpha_perm = 0.05
n_rep = 1000

repfun <- function() {
  s <- series_sim(series_length, event_rate)

  # t <- ceiling(series_length / 2)
  t <- 1

  X  <- sum(s)
  x1 <- sum(head(s, t))
  x0 <- X - x1

  # stopifnot(x0 + x1 == sum(s))

  data.table(Ftest = pfpr(series_length - t, x0, t, x1),
             RRtest = rateratio.test(c(x0, x1), c(series_length - t, t), alternative = "less")$p.value)
}

# res <- rbindlist(replicate(10000, repfun(), simplify = FALSE))

# colMeans(as.matrix(res) < 0.05)

# stop()


QUIET <- clusterEvalQ(cl, {
  rm(list = setdiff(objects(), ""))
  library(data.table)
  library(rateratio.test)
})

QUIET <- clusterExport(cl, setdiff(objects(), "cl"))

res <- rbindlist(clusterEvalQ(cl, {
  rbindlist( replicate(n_rep, repfun(), simplify = FALSE) )
}))

colMeans(as.matrix(res) < 0.05)
colMeans(as.matrix(res) < 0.01)
colMeans(as.matrix(res) < 0.001)

# TEST CALIBRATION

jpeg("single t calibration t1.jpeg", 3000, 3000, pointsize = 96, quality = 95)
{
  par(mar = c(4,4,2,2))
  plot(ecdf(res$Ftest), col = "red", lwd = 4,
       xlab = "Theoretical P-value", ylab = "Observed P-value",
       main = "Test calibration")
  lines(ecdf(res$RRtest), col = "blue", lwd = 4)
  abline(0, 1, lty = 2)
  legend("top", bty = "n", fill = c("blue", "red"), legend = c("Poisson", "Fuzzy Poisson"))
  grid()
}
dev.off()


stop()

# jpeg("calibration .jpeg", 4000, 4000, pointsize = 96, quality = 95)
{
  par(mfrow = c(2,2))
  par(mar = c(4,4,2,2))

  # TEST CALIBRATION
  plot(ecdf(res$pval_center), col = "red", pch = 19, cex = 0.2,
       xlab = "Theoretical P-value", ylab = "Observed P-value",
       main = "Test calibration")
  abline(0, 1, lty = 2)
  grid()


  hist(res$t, freq = FALSE, xlab = "MSD change point before present", col = "lightblue",
       main = "Change point distribution")

  boxplot(pval_center ~ t, res, col = "lightblue", outline = FALSE,
          xlab = "MSD change point before present", ylab = "P-value")

  boxplot(stat ~ t, res, col = "lightblue", outline = FALSE,
          xlab = "MSD change point before present", ylab = "Test statistic (rate ratio)")


  par(mfrow = c(1,1))
}
# dev.off()

mean(res$pval_center < 0.05)
mean(res$pval_center < 0.01)
mean(res$pval_center < 0.001)

# jpeg("change point histogram rate01.jpeg", 2000, 2000, pointsize = 96, quality = 95)
{
  hist(res[pval_center < 0.05]$t, breaks = series_length/2, freq = FALSE, col = rgb(0.8,0,0,0.3),
       main = "Change point distribution", ylim = c(0,0.1), xlab = "MSD change point before present")
  hist(res[pval_center > 0.05]$t, breaks = series_length/2, freq = FALSE, add = TRUE, col = rgb(0,0,0.8,0.3))
  legend("top", bty = "n", fill = c(rgb(0.8,0,0,0.3), rgb(0,0,0.8,0.3)), legend = c("P < 0.05", "P >= 0.05"))

}
# dev.off()

stop()

################################################
# COMPARISON WITH PREVIOUS TEST rateratio.test

library(rateratio.test)

pfpr <- function(t0, x0, t1, x1) pf(
  (x0 + 1) * t1 / ((x1 + 1) * t0),
  2*(x1 + 1), 2*(x0 + 1)
)

{
  x0 <- 0
  t0 <- 100
  x1 <- 0.15
  t1 <- 1


  # print(pfpr(t0, x0, t1, x1) + 1 - pfpr(t1, x1, t0, x0))
  print(pfpr(t0, x0, t1, x1))

  print(rateratio.test(c(x0, x1), c(t0, t1), alternative = "l")$p.value)
}


stop()



boxplot(t ~ I(pval_center < 0.01), res)

# CHANGE POINT CALIBRATION
plot(ecdf(res$t), xlab = "MSD change point before present",
     ylab = "Cumulative proportion", main = "Change point calibration")
abline(0, 1/series_length, lty = 2)
grid()

# Check distribution of counts in series: looks OK
# res2 <- rbindlist(clusterEvalQ(cl, {data.frame(t(replicate(1000,
#     series_sim(series_length, event_rate)
#   )))}))
#
# plot(colSums(res2))

stop()


z <- data.table(t(replicate(1000, {
  s <- series_sim(series_length, event_rate)
  msd(s, n_perm, alpha_perm)
})))


z

hist(z$pval_center, breaks = 20)

qqnorm(qnorm(z$pval_center[z$pval_center > 0]))
qqline(qnorm(z$pval_center[z$pval_center > 0]))

mean(z$pval_center < 0.10)
mean(z$pval_center < 0.05)
mean(z$pval_center < 0.01)

plot(ecdf(z$pval_center), col = "red", pch = 19, cex = 0.2)
abline(0, 1, lty = 2)
grid()

# qqline(z$pval_center, distribution = runif)

# s

msd(s)



# Compute permutations




print(mean(perm_series < p_min))

# hist(perm_series)

# Check calibration
# => nominal error rate should not depend on rate or length
# => distribution of most significant t should be uniform

hist(z$t)

stop()




# cut(cw)

cbind(cw, findInterval(cw, 0:ts))


t1 <- 1
t0 <- ts - t1







stop()



# Rethink a little... ;-)

# Suppose true events at rate r_t and false events at rate r_f
# There're independent so total rate is r_t + r_f
# Each event is true with probability p or false with probability
# (1-p), measure rate may be p r_t + (1 - p) r_f ?


library(data.table)

# Simulate a fuzzy Poisson over a time series S

Ts <- 100 # Time series length
tau <- 10 # Change point before present
r <- 1    # Event rate, constant over S

# Avoid the normalization integral here.

# Split the process into a time series with exact exponential waiting times

m <- 1e3

w <- rexp(m, r)

# Then assume that each true event has some beta distributed probability
# with mean p

wprob <- rbeta(m, 1, 1)

# Then the fuzzy poisson count is the cumulative sum over the period

d <- data.table(w = w, t = cumsum(w), p = wprob, fuzzycount = cumsum(wprob))

# ratio <- (x0 + 1) * t1 / ((x1 + 1) * t0)
# pfpr <- function(t0, x0, t1, x1) pf(
#   (x0 + 1) * t1 / ((x1 + 1) * t0),
#   2*(x1 + 1), 2*(x0 + 1)
# )

# For some length Ts and change point before present tau, compute x1 and x0 from
# teh time series, then do the same for all tau's

ts <- 365
t1 <- 1
t0 <- ts - t1

# Check series is long enough
stopifnot(tail(d$t, 1) > ts)

# Find index in the series
t1_index <- max(which(d$t < t1)) # Returns -Inf if no event
ts_index <- max(which(d$t < ts))

# Count over the whole period
xs <- d[ts_index]$fuzzycount

# Count over the x1 period
x1 <- d[t1_index]$fuzzycount
x0 <- xs - x1

ratio <- (x0 + 1) * t1 / ((x1 + 1) * t0)
pfpr <- function(t0, x0, t1, x1) pf(
  (x0 + 1) * t1 / ((x1 + 1) * t0),
  2*(x1 + 1), 2*(x0 + 1)
)

pfpr(t0, x0, t1, x1)

# Steps for MSD test: given a series, compute ratios and P-values,
# accumulate permutations and store P-value distribution.
# Bonus : given a threshold and confidence, stop permutations when
# the P-values is above or below threshold with enough confidence based on
# Wald confint or Clopper-Pearson,
# https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval#Clopper%E2%80%93Pearson_interval


stop()




mean(wprob)
hist(wprob)



# Draw random deviates from the fuzzy Poisson ?

# Looks like we need the normalizing constant first to obtain the deviates ?





# Note on incomplete Gamma here: https://search.r-project.org/CRAN/refmans/expint/html/gammainc.html
# gammainc(a, x)
# gamma(a) * pgamma(x, a, 1, lower = FALSE) # same

gammainc <- function(a, b) gamma(a) * pgamma(b, a, 1, lower = FALSE)

# Fuzzy Poisson CDF
pfpois <- function(q, r) pgamma(r, q + 1, 1, lower = FALSE)

curve(ppois(x, 1), 0, 10)
curve(pfpois(x, 1), add = TRUE, col = "red")


# Fuzzy Poisson quantile
# qfpois <-



# Fuzzy Poisson quantile function ?

# Fuzzy Poisson density (unnormalized)
dfpois <- function(x, r) exp(-r) * r^x / gamma(x+1)

curve(dfpois(x, 1), 0, 10)

r <- 1

integrate(function(x, r) r^x / gamma(x+1), 0, 0.01, r = r)$value * 0.01

exp(-r) / gammainc(1, r)

r <- 1/2

dfpois(0, r)

dfpois(0, r) / integrate(dfpois, 0, Inf, r = r)$value

h <- 0.000001

((pfpois(h, r) - pfpois(0, r)) / h)

exp(-r)*(log(r) - digamma(1) )

###################################

x0 <- 10
t0 <- 100
x1 <- 20
t1 <- 10

ratio <- (x0 + 1) * t1 / ((x1 + 1) * t0)

curve(df(x, 2*(x1 + 1), 2*(x0 + 1)), 0, 4)

pf(ratio, 2*(x1 + 1), 2*(x0 + 1))

# Fuzzy Poisson rate distribution (FPR)

dfpr <- function(t0, x0, t1, x1) df(
  (x0 + 1) * t1 / ((x1 + 1) * t0),
  2*(x1 + 1), 2*(x0 + 1)
)

pfpr <- function(t0, x0, t1, x1) pf(
  (x0 + 1) * t1 / ((x1 + 1) * t0),
  2*(x1 + 1), 2*(x0 + 1)
)

pfpr(10, 1, 10, 10)

# pdf("CPM F-test curves.pdf", 8, 8)
{
  par(mfrow = c(2,2))

  # DENSITY OF STATISTIC UNDER EQUAL OBSERVATION INTERVALS
  xseq <- c(0, 1, 10, 50)
  cramp <- colorRampPalette(c("blue", "green", "orange", "red"))(length(xseq))
  curve(df(x, 2*(xseq[1] + 1), 2*(xseq[1] + 1)), 0, 4,
        xlab = "CPM statistic", ylab = "F density", lwd = 2, ylim = c(0, 2),
        col = cramp[1], main = "Equal length observation intervals")
  for(i in 2:length(xseq)) {
    curve(df(x, 2*(xseq[i] + 1), 2*(xseq[i] + 1)), add = TRUE, lwd = 2, col = cramp[i])
  }
  legend("right", bty = "n", fill = cramp, legend = xseq, title = "Fuzzy event count\nper interval")

  # P-VALUE OVER SINGLE DAY FUZZY COUNT
  xseq <- c(1, 2, 3, 4)
  curve(pfpr(365, 365, xseq[1], x), 0, 10, lwd = 2, col = cramp[1], ylim = c(0, 1),
        xlab = "Fuzzy event count", ylab = "Equal rate P-value",
        main = "1 event/day, 1 year left interval")
  for(i in 2:length(xseq)) {
    curve(pfpr(365, 365, xseq[i], x), add = TRUE, lwd = 2, col = cramp[i])
  }
  legend("right", bty = "n", fill = cramp, legend = xseq, title = "Right interval duration (days)")
  abline(h = 0.05, lty = 2)

  # IMPACT OF LEARNING PHASE LENGTH
  xseq <- c(1, 7, 30, 365)
  curve(pfpr(xseq[1], xseq[1], 1, x), 0, 10, lwd = 2, col = cramp[1], ylim = c(0, 1),
        xlab = "Fuzzy event count", ylab = "Equal rate P-value",
        main = "1 event/day left interval\n1-day right interval")
  for(i in 2:length(xseq)) {
    curve(pfpr(xseq[i], xseq[i], 1, x), add = TRUE, lwd = 2, col = cramp[i])
  }
  legend("right", bty = "n", fill = cramp, legend = xseq, title = "Left interval duration (days)")
  abline(h = 0.05, lty = 2)

  # IMPACT OF LEARNING PHASE LENGTH
  xseq <- c(1, 7, 30, 365)
  curve(pfpr(xseq[1], 0, 1, x), 0, 4, lwd = 2, col = cramp[1], ylim = c(0, 0.5),
        xlab = "Fuzzy event count", ylab = "Equal rate P-value",
        main = "Event-free left interval\n1-day right interval")
  for(i in 2:length(xseq)) {
    curve(pfpr(xseq[i], 0, 1, x), add = TRUE, lwd = 2, col = cramp[i])
  }
  legend("right", bty = "n", fill = cramp, legend = xseq, title = "Left interval duration (days)")
  abline(h = 0.05, lty = 2)

  par(mfrow = c(1,1))
}
# dev.off()

