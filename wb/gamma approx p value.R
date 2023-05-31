# A simple Gamma approximation to permutation P-value

m <- 10000
n <- 10

p <- replicate(m, -log(min(runif(n)^0.5)))

hist(p)


# px <- -log(min(runif(n)))

opt <- optim(c(0,0), function(par) {
  -sum(dgamma(p, exp(par[1]), exp(par[2]), log = TRUE))
})

hist(p, freq = FALSE)
curve(dgamma(x, exp(opt$par[1]), exp(opt$par[2])), add = TRUE, lwd = 2, col = "red")

# Moments-based estimator of Gamma parameters
gm <- mean(p)
gv <- var(p)

a <- gm*gm/gv
b <- gm/gv

curve(dgamma(x, a, b), add = TRUE, lwd = 2, col = "blue")
