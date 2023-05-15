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

## quantile plot of central range
addQuantPoly <- function(qnts = c(0.005, 0.025, 0.25),
                         paths = simBetaPath(a, b),
                         kseq = exp(seq(-8, 8, by = 0.25)),
                         labpos = c("left", "top"), ...) {
    for (qn in qnts) {
        polygon(log(c(kseq, rev(kseq)), base = 10),
                c(apply(log(paths, base = 10), 1, quantile, 
                        probs = qn),
                  rev(apply(log(paths, base = 10), 1, quantile,
                            probs = 1-qn))),
                col = adjustcolor("gray", 0.25), border = NA)
        if (labpos[1] == "left") {
            xind <- 1
        } else if (labpos[1] == "right") {
            xind <- length(kseq)
        }
        if (labpos[2] == "top") {
            yq <- 1-qn
        } else if (labpos[2] == "bottom") {
            yq <- qn
        }
        text(x = log(kseq[xind], 10),
             y = quantile(log(paths[xind,], base = 10), yq),
             labels = 1-2*qn, adj = c(0.5, 0.5), ...)
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

## load the null curve min quantiles (100,000 reps)
nullQuants <- readRDS("curveMinQuantiles.Rds")

## SINGLE DENSITY ####################################################
##' cases for chapter
##' a = 1, b = 1, n = 1e2
##' a = 1, b = 0.5, n = 1e2
##' a = 0.5, b = 1, n = 1e2
##' a = 0.4, b = 0.4, n = 1e2
##' a = 0.4, b = 0.2, n = 1e2
##' a = 4, b = 4, n = 1e2
##' a = 2, b = 4, n = 1e2

## call based on a, b, simulation settings
nsim <- 1e3
n <- 1e2
a <- 1
b <- 1
kseq <- exp(seq(-8, 8, by = 0.25))
sims <- simBetaPath(a = a, b = b, n = n, nsim = nsim,
                    kseq = kseq)

## plot these results
xpos <- seq(0, 1, 0.005)
png(paste0("beta", a, b, "dens.png"), width = 360, height = 360)
## the density
narrowPlot(xgrid = seq(0, 1, by = 0.2), ygrid = seq(0, 5, by = 1),
           main = paste0("Beta(", a, ", ", "", b, ")"),
           xlab = "x", ylab = "Density")
lines(xpos, dbeta(xpos, a, b))
dev.off()

png(paste0("beta", a, b, "quants.png"), width = 360, height = 360)
## the central quantiles
narrowPlot(xgrid = seq(-3, 3, by = 1), ygrid = seq(-15, 0, by = 3),
           xlim = c(-3.5, 3.5),
           main = "Pooled p-value central quantiles",
           xlab = expression(paste(log[10], "(", kappa, ")")),
           ylab = expression(paste(log[10], "(p)")))
addQuantPoly(paths = sims, cex = 0.8, labpos = c("right", "bottom"),
             kseq = kseq)
quantLevs <- c("5%", "1%", "0.1%")
abline(h = log(nullQuants[quantLevs, as.character(n)], 10),
       lty = 2, col = "firebrick")
text(x = rep(-3.2, 3), labels = quantLevs, cex = 0.8,
     y = log(nullQuants[quantLevs, as.character(n)], 10),
     adj = c(0.5, -0.2))
dev.off()

## MIXTURE ###########################################################
##' settings to include:
##' (0.1, 1), (1, 1), n = 1e2, p = 0.05
##' (3, 4), (0.4, 1), n = 1e2, p = 0.8

## now try a mixture
nsim <- 1e3
n <- 1e2
p <- 0.05
b1 <- list(a = 0.1, b = 1)
b2 <- list(a = 1, b = 1)
kseq <- exp(seq(-8, 8, by = 0.25))
mixsims <- simBetaMix(p = p, b1 = b1, b2 = b2, n = n,
                      nsim = nsim, kseq = kseq)

## plot these results
xpos <- seq(0, 1, 0.005)
png(paste0("betamix", b1$a, b1$b,"-", b2$a, b2$b, p, "dens.png"),
    width = 360, height = 360)
## the density
narrowPlot(xgrid = seq(0, 1, by = 0.2), ygrid = seq(0, 5, by = 1),
           main = paste0("Beta(", b1$a, ", ", b1$b, ") and Beta(",
                         b2$a, ", ", b2$b, ") mixture (", p, "/", 1-p,
                         ")"),
           xlab = "x", ylab = "Density")
lines(xpos, p*dbeta(xpos, b1$a, b1$b) + (1-p)*dbeta(xpos, b2$a, b2$b))
dev.off()

png(paste0("betamix", b1$a, b1$b,"-", b2$a, b2$b, p, "quants.png"),
    width = 360, height = 360)
## the central quantiles
narrowPlot(xgrid = seq(-3, 3, by = 1), ygrid = seq(-15, 0, by = 3),
           xlim = c(-3.5, 3.5),
           main = "Pooled p-value central quantiles",
           xlab = expression(paste(log[10], "(", kappa, ")")),
           ylab = expression(paste(log[10], "(p)")))
addQuantPoly(paths = mixsims, cex = 0.8, labpos = c("right", "bottom"))
quantLevs <- c("5%", "1%", "0.1%")
abline(h = log(nullQuants[quantLevs, as.character(n)], 10),
       lty = 2, col = "firebrick")
text(x = rep(-3.2, 3), labels = quantLevs, cex = 0.8,
     y = log(nullQuants[quantLevs, as.character(n)], 10),
     adj = c(0.5, -0.2))
dev.off()

## NULL QUANTILES ####################################################
##' get min values across curves for the null case
##' meant to be run in a large computing environment with many cores
tol <- 1e6 # tolerance for load balancing
nsim <- 1e5 # number of simulations
M <- c(2, 10, 100, 500, 1000, 10000) # sample sizes
total <- nsim*M # total random numbers
batches <- vector(mode = "list", length = length(M))
for (ii in seq_along(M)){ # split into batches to control computation
    if (total[ii] > tol) {
        batches[[ii]] <- rep(tol, total[ii]/tol)/M[ii]
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
nullCurves <- parLapply(clust, # run batches across cluster
                        seq_along(Ms),
                        function(ii) simBetaPath(a = 1, b = 1,
                                                 n = Ms[ii],
                                                 nsim = batches[ii],
                                                 kseq = kseq))
stopCluster(clust)
curvesByM <- split(nullCurves, Ms)
nullCurves <- lapply(testM, function(el) do.call(cbind, el))
saveRDS(list(Ms, nullCurves), file = "nullCurves.Rds")
