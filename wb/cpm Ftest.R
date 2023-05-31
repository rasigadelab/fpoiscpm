##################################################################
# EPITRACK - NOTES ON PROCESS CONTROL

# This script generates figures of the density
# of test statistic.

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

