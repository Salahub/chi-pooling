## simple helper
klDivU <- function(tab) {
    tot <- sum(tab)
    unif <- tot/length(tab)
    sum((unif/tot)*log(unif/tab))
}

## chi-square pool function
poolChi <- function(p, k) {
    M <- length(p) # dimension
    pchisq(sum(qchisq(p, df = k, lower.tail = FALSE)), df = M*k,
           lower.tail = FALSE)
}

## looking at pooled p-values paths for beta distributions
## first assay to get a sense
step <- 0.5
a <- seq(0.5, 5, by = step)
b <- seq(0.5, 5, by = step)
params <- expand.grid(b = b, a = a)
betas <- lapply(seq_len(nrow(params)),
                function(ii) function(n) rbeta(n, params$a[ii],
                                               params$b[ii]))
## generate random p-values
nsim <- 1e4
samps <- lapply(betas, function(f) f(nsim))
## sweep with pooled p kappas
kseq <- exp(seq(-8, 8, by = 0.5))
pools <- lapply(samps, function(smp) sapply(kseq, poolChi, p = smp))

## plot one example
ii <- 13
par(mfrow = c(1,2))
plot(seq(0, 1, 0.01), xlab = "x", ylab = "Density", type = "l",
     main = paste("a =", params$a[ii], ",", "b =", params$b[ii]),
     dbeta(seq(0, 1, 0.01), params$a[ii], params$b[ii]),
     ylim = c(0,5))
plot(log(kseq), log(pools[[ii]]), main = "Chi p-value path",
     xlab = expression(paste0("log(", kappa, ")")),
     ylab = "log(p)", type = "l", ylim = c(-10, 1))
abline(h = log(0.05), lty = 2)

## make this more robust with repeated samples...
simBetaPath <- function(a = 1, b = 1, n = 1e3, nsim = 100,
                        kseq = exp(seq(-8, 8, by = 0.25))) {
    betaSamps <- matrix(rbeta(n*nsim, a, b), ncol = n) # sim matrix
    betaTrans <- lapply(kseq, qchisq, p = betaSamps,
                        lower.tail = FALSE) # quantile values
    paths <- mapply(function(mt, k) pchisq(rowSums(mt), df = k*n,
                                           lower.tail = FALSE),
                    betaTrans, kseq)
    t(paths)
}

## and construct a similar plot
plotBetaPath <- function(a, b, paths = simBetaPath(a, b),
                         kseq = exp(seq(-8, 8, by = 0.25)),
                         ylim1 = c(0,5), ylim2 = c(-10,1)) {
    par(mfrow = c(1,2))
    plot(seq(0, 1, 0.01), xlab = "x", ylab = "Density", type = "l",
         main = paste("Beta(a =", a, ",", "b =", b, ")"),
         dbeta(seq(0, 1, 0.01), a, b), ylim = ylim1)
    plot(NA, xlim = range(log(kseq, base = 10)),
         main = "Chi p-value path",
         xlab = expression(paste(log[10], "(", kappa, ")")),
         ylab = expression(paste(log[10], "(p)")), type = "n",
         ylim = ylim2)
    for (jj in seq_len(ncol(paths))) {
        lines(log(kseq, base = 10), log(paths[,jj], base = 10),
              col = adjustcolor("black", 0.6))
    }
    abline(h = log(0.05, base = 10), lty = 2,
           col = adjustcolor("firebrick", 0.5))
}

## what about a quantile plot?
plotBetaQuants <- function(a, b, qnts = c(0.005, 0.025, 0.25),
                           paths = simBetaPath(a, b),
                           kseq = exp(seq(-8, 8, by = 0.25)),
                           ylim1 = c(0,5), ylim2 = c(-10,1)) {
    par(mfrow = c(1,2))
    plot(seq(0, 1, 0.005), xlab = "x", ylab = "Density", type = "l",
         main = paste("Beta(a =", a, ",", "b =", b, ")"),
         dbeta(seq(0, 1, 0.005), a, b), ylim = ylim1)
    abline(v = seq(0, 1, by = 0.2), h = seq(0, 5, by = 1),
           col = "gray90")
    plot(NA, xlim = range(log(kseq, base = 10)),
         main = "Chi p-value path quantiles",
         xlab = expression(paste(log[10], "(", kappa, ")")),
         ylab = expression(paste(log[10], "(p)")), type = "n",
         ylim = ylim2)
    abline(h = seq(-10, 0, by = 2), v = seq(-3, 3, by = 1),
           col = "gray90")
    for (qn in qnts) {
        polygon(log(c(kseq, rev(kseq)), base = 10),
                c(apply(log(paths, base = 10), 1, quantile, probs = qn),
                  rev(apply(log(paths, base = 10), 1, quantile,
                            probs = 1-qn))),
                col = adjustcolor("gray", 0.25), border = NA)
    }
    lines(x = log(kseq, 10),
          y = apply(log(paths, base = 10), 1, median))
    abline(h = log(0.05, base = 10), lty = 2,
           col = adjustcolor("firebrick", 0.5))
}

## try a mixture of two betas as well
simBetaMix <- function(p = 0.5, b1 = list(a = 1, b = 1),
                       b2 = list(a = 1, b = 1),
                       n = 100, nsim = 1e3,
                       kseq = exp(seq(-8, 8, by = 0.25))) {
    beta1 <- runif(n*nsim) <= p # sampled from first beta
    betaSeq <- matrix(ncol = n, nrow = nsim) # pre-allocate
    betaSeq[beta1] <- rbeta(sum(beta1), b1$a, b1$b)
    betaSeq[!beta1] <- rbeta(sum(!beta1), b2$a, b2$b)
    betaTrans <- lapply(kseq, qchisq, p = betaSeq,
                        lower.tail = FALSE) # quantiles
    paths <- mapply(function(mt, k) pchisq(rowSums(mt), df = k*n,
                                           lower.tail = FALSE),
                    betaTrans, kseq)
    t(paths)
}

## what about a quantile plot?
plotBetaMix <- function(p, b1, b2, qnts = c(0.005, 0.025, 0.25),
                        paths = simBetaMix(p, b1, b2),
                        kseq = exp(seq(-8, 8, by = 0.25)),
                        ylim1 = c(0,5), ylim2 = c(-10,1)) {
    x <- seq(0, 1, 0.005)
    par(mfrow = c(1,2))
    plot(x = x, xlab = "x", ylab = "Density", type = "l",
         main = paste0("Mixture of ", b1$a, ", ", b1$b, " and ",
                       b2$a, ", ", b2$b, " (", p, "/", 1-p, ")"),
         y = p*dbeta(x, b1$a, b1$b) + (1-p)*dbeta(x, b2$a, b2$b),
         ylim = ylim1)
    abline(v = seq(0, 1, by = 0.2), h = seq(0, 5, by = 1),
           col = "gray90")
    plot(NA, xlim = range(log(kseq, base = 10)),
         main = "Chi p-value path quantiles",
         xlab = expression(paste(log[10], "(", kappa, ")")),
         ylab = expression(paste(log[10], "(p)")), type = "n",
         ylim = ylim2)
    for (qn in qnts) {
        polygon(log(c(kseq, rev(kseq)), base = 10),
                c(apply(log(paths, base = 10), 1, quantile, probs = qn),
                  rev(apply(log(paths, base = 10), 1, quantile,
                            probs = 1-qn))),
                col = adjustcolor("gray", 0.25), border = NA)
    }
    lines(x = log(kseq, 10),
          y = apply(log(paths, base = 10), 1, median))
    abline(h = log(0.05, base = 10), lty = 2,
           col = adjustcolor("firebrick", 0.5))
}

## call based on a, b, simulation settings
nsim <- 1e3
n <- 1e3
a <- 1.4
b <- 1.5
kseq <- exp(seq(-8, 8, by = 0.1))
sims <- simBetaPath(a = a, b = b, n = n, nsim = nsim)
plotBetaQuants(a = a, b = b, paths = sims)
##' some interesting settings:
##' a = 0.4, b = 0.2
##' a = 0.82, b = 0.85
##' a = 3, b = 4
##' a = 4, b = 3
##' a = 0.85, b = 0.82
##' a = 0.2, b = 0.4
##' a = 1, b = 2
##' a = 0.9, b = 1.1

## now try a mixture
nsim <- 1e3
n <- 1e3
p <- 0.2
b1 <- list(a = 0.4, b = 0.2)
b2 <- list(a = 3, b = 5)
kseq <- exp(seq(-8, 8, by = 0.1))
mixsims <- simBetaMix(p = p, b1 = b1, b2 = b2, n = n,
                      nsim = nsim, kseq = kseq)
plotBetaMix(p = p, b1 = b1, b2 = b2, paths = mixsims, kseq = kseq,
            ylim2 = c(-25, 0))
abline(h = seq(-25, 0, by = 5), v = seq(-3, 3, by = 1),
       col = adjustcolor("gray", 0.4))
##' settings to include:
##' (1, 1), (0.4, 1), range of p
##' (3, 4), (0.4, 1), p = 0.5, 0.8
##' (1, 0.4), (0.4, 1), p = 0.5
##' (1, 1), (3, 5), p = 0.3
##' (0.4, 0.2), (3, 5), p = 0.3
