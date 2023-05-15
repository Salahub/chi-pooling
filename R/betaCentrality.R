## custom plotting function with narrow margins
narrowPlot <- function(xgrid, ygrid, main = "", xlab = "", ylab = "",
                       xticks = xgrid, yticks = ygrid,
                       mars = c(2.1, 2.1, 1.1, 1.1),
                       xlim = range(xgrid), ylim = range(ygrid),
                       addGrid = TRUE, ...) {
    par(mar = mars) # set narrow margins
    plot(NA, ylim = ylim, xlim = xlim, xaxt = 'n', xlab = "",
         yaxt = 'n', ylab = "", main = "", ...)
    ## add labels
    mtext(main, side = 3, line = 0, cex = 0.8) # main
    mtext(ylab, side = 2, line = 1, cex = 0.8) # ylab
    mtext(xlab, side = 1, line = 1, padj = 0, cex = 0.8) # xlab
    ## add grid lines
    abline(h = ygrid, v = xgrid, lty = 1,
           col = adjustcolor("gray", alpha.f = 0.4))
    ## and ticks
    mtext(side = 1, at = xgrid, text = "|", line = 0, cex = 0.5,
          padj = -2)
    mtext(text = xticks, at = xgrid, side = 1, cex = 0.8)
    mtext(side = 2, at = ygrid, text = "|", line = 0, cex = 0.5,
          padj = 1)
    mtext(text = yticks, at = ygrid, side = 2, cex = 0.8)
}

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

## make this more robust with repeated samples...
simBetaPath <- function(a = 1, b = 1, n = 1e3, nsim = 100,
                        kseq = exp(seq(-8, 8, by = 0.25))) {
    betaSamps <- matrix(rbeta(n*nsim, a, b), ncol = n) # sim matrix
    betaTrans <- lapply(kseq, qchisq,
                        p = betaSamps,
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
addQuantPoly <- function(qnts = c(0.005, 0.025, 0.25),
                         paths = simBetaPath(a, b),
                         kseq = exp(seq(-8, 8, by = 0.25)),
                         ...) {
    for (qn in qnts) {
        polygon(log(c(kseq, rev(kseq)), base = 10),
                c(apply(log(paths, base = 10), 1, quantile, probs = qn),
                  rev(apply(log(paths, base = 10), 1, quantile,
                            probs = 1-qn))),
                col = adjustcolor("gray", 0.25), border = NA)
        text(x = log(kseq[1], 10), labels = 1-2*qn,
             y = log(quantile(paths[1,], 1-qn), 10),
             adj = c(0, 0.5), ...)
    }
    lines(x = log(kseq, 10),
          y = apply(log(paths, base = 10), 1, median))
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
n <- 1e2
a <- 0.4
b <- 0.4
kseq <- exp(seq(-8, 8, by = 0.1))
sims <- simBetaPath(a = a, b = b, n = n, nsim = nsim)

## plot these results
xpos <- seq(0, 1, 0.005)
## the density
narrowPlot(xgrid = seq(0, 1, by = 0.2), ygrid = seq(0, 5, by = 1),
           main = paste("Beta(", a, ",", "", b, ")"),
           xlab = "x", ylab = "Density")
lines(xpos, dbeta(xpos, a, b))

narrowPlot(xgrid = seq(-3, 3, by = 1), ygrid = seq(-15, 0, by = 3),
           xlim = c(-3.5, 3.5),
           main = "Pooled p-value central quantiles",
           xlab = expression(paste(log[10], "(", kappa, ")")),
           ylab = expression(paste(log[10], "(p)")))
addQuantPoly(paths = sims, cex = 0.8)
abline(h = log(0.05, base = 10), lty = 2,
       col = adjustcolor("firebrick", 0.5))


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
p <- 0.3
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

##' get min values across curves for the null case
tol <- 1e6 # tolerance for load balancing
nsim <- 1e4 # number of simulations
M <- c(2, 10, 100, 500, 1000, 10000) # sample sizes
total <- nsim*M # total random numbers
batches <- vector(mode = "list", length = length(M))
for (ii in seq_along(M)){ # split into batches
    if (total[ii] > tol) {
        batches[[ii]] <- rep(tol, total[ii]/tol)
    } else batches[[ii]] <- nsim
}
Ms <- rep(M, times = sapply(batches, length)) # repeat M for batches
batches <- unlist(batches)
kseq <- exp(seq(-8, 8, by = 0.1))
library(parallel)
ncores <- detectCores()/2
clust <- makeCluster(ncores)
clusterExport(clust, varlist = c("simBetaPath", "Ms", "batches",
                                 "kseq", "nsim"))
nullCurves <- lapply(M, # store object
                     function(m) simBetaPath(a = 1, b = 1,
                                             n = m, nsim = nsim,
                                             kseq = kseq))
stopCluster(clust)
saveRDS(nullCurves, file = "nullCurves.Rds")
