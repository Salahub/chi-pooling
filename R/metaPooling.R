library(metadat)

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
    if (addGrid) {
        abline(h = ygrid, v = xgrid, lty = 1,
               col = adjustcolor("gray", alpha.f = 0.4))
    }
    ## and ticks
    mtext(side = 1, at = xgrid, text = "|", line = 0, cex = 0.5,
          padj = -2)
    mtext(text = xticks, at = xgrid, side = 1, cex = 0.8)
    mtext(side = 2, at = ygrid, text = "|", line = 0, cex = 0.5,
          padj = 1)
    mtext(text = yticks, at = ygrid, side = 2, cex = 0.8)
}

## quantile plot of central range
addQuantPoly <- function(mat, qnts = c(0.005, 0.025, 0.25),
                         kseq = exp(seq(-8, 8, by = 0.25)),
                         labpos = c("left", "top"), ...) {
    for (qn in qnts) {
        polygon(log(c(kseq, rev(kseq)), base = 10),
                c(apply(mat, 1, quantile,
                        probs = qn),
                  rev(apply(mat, 1, quantile,
                            probs = 1-qn))),
                col = adjustcolor("gray", 0.25), border = NA)
        if (labpos[1] == "left") {
            xind <- 1
            adj <- c(0, 0.5)
        } else if (labpos[1] == "right") {
            xind <- length(kseq)
            adj <- c(1, 0.5)
        }
        if (labpos[2] == "top") {
            yq <- 1-qn
        } else if (labpos[2] == "bottom") {
            yq <- qn
        }
        text(x = log(kseq[xind], 10),
             y = quantile(mat[xind,], yq),
             labels = 1-2*qn, adj = adj, ...)
    }
    lines(x = log(kseq, 10),
          y = apply(mat, 1, median))
}

## parameter sweep function
paramSweep <- function(mus, sds, thetas) {
    mus2 <- rep(mus, each = length(thetas))
    sds2 <- rep(sds, each = length(thetas))
    matrix(2*pnorm(abs((mus2 - rep(thetas, times = length(mus)))/sds2),
                   lower.tail = FALSE),
           nrow = length(thetas))
}

## sweeping with kappa afterwards
kappaSweep <- function(ps, kseq = exp(seq(-8, 8, by = 0.1))) {
    sapply(kseq,
           function(kap) {
               apply(ps, 2,
                     function(p) pchisq(sum(qchisq(p, kap,
                                               lower.tail = FALSE)),
                                        length(p)*kap,
                                        lower.tail = FALSE))
           })
}

## repeat this process many times to compare bias...
mn <- 0
std <- 2
nsim <- 1000
ngroup <- 8
npop <- 104
minQuants <- readRDS("curveMinQuantiles.Rds")
kseq <- exp(seq(-8, 8, by = 0.2))
thetas <- seq(-2.5, 2.5, by = 0.01)
muMat <- sdMat <- matrix(nrow = nsim, ncol = ngroup)
nMat <- matrix(nrow = nsim, ncol = ngroup)
pMat <- array(dim = c(nsim, ngroup, length(thetas)))
poolMat <- array(dim = c(nsim, length(kseq), length(thetas)))
set.seed(9381278)
for (ii in 1:nsim) {
    ## start with simulated data
    obs <- rnorm(npop, mean = mn, sd = std)
    ## randomly split
    group <- sample(rep(1:ngroup, times = npop/ngroup))
    #group <- sample(c(sample(1:ngroup),
    #                  sample(1:ngroup, npop - ngroup, replace = TRUE)))
    nMat[ii,] <- table(group)
    ## get means and sds from split groups
    muMat[ii,] <- mus <- tapply(obs, group, mean)
    ## get sds from split groups
    sdMat[ii,] <- sds <- tapply(obs, group, sd)/sqrt(table(group))
    ## simulated p-values
    pMat[ii,,] <- simps <- t(paramSweep(mus, sds, thetas))
    ## the pooled p-value curves
    poolMat[ii,,] <- simpool <- t(kappaSweep(simps, kseq = kseq))
    if ((ii %% 10) == 0) cat("\r Done ", ii, " of ", nsim)
}
## compute the mini-max kappa on each repetition
minimax <- apply(apply(poolMat, c(1,2), max), 1, which.min)
minpool <- apply(poolMat, c(1,3), min)

## plot a realization
cols <- RColorBrewer::brewer.pal(4, "Dark2")
real <- 80 # sample(1:nsim, 1)
png("metaPoolCurves.png", width = 5, height = 3, res = 480,
    units = "in")
narrowPlot(xgrid = seq(-2, 2, by = 1), ygrid = seq(0, 1, by = 0.2),
           xlab = expression(hat(theta)),
           ylab = expression(paste("chi(", bold(p), ";", kappa, ")")),
           addGrid = FALSE)
lines(thetas,  minpool[real, ], type = 'l',
      col = adjustcolor(cols[4], 0.8))
      ##poolMat[real, minimax[real], ], type = 'l',
for (ii in c(1, 41, 81)) lines(thetas, poolMat[real,ii,], type = 'l',
                         col = if (ii <= 40) adjustcolor(cols[1], 0.8) else if (ii == 41) adjustcolor(cols[2], 0.8) else adjustcolor(cols[3], 0.8))
                                        #lwd = if (log(kseq)[ii] == 0) 2 else 1)
#lines(thetas, poolMat[real,81,], col = adjustcolor(cols[3], 0.8))
abline(v = 0, col = adjustcolor("black", 0.8), lty = 2, lwd = 1)
abline(v = mean(muMat[real,]), col = adjustcolor("black", 0.6), lwd = 1)
abline(v = muMat[real,], col = adjustcolor("gray50", 0.4))
abline(h = 0.05, lty = 3)
legend(x = "topleft", legend = c(round(log(kseq[c(1, 41, 81)], 10),
                                       1), "Min"),
       lty = 1, col = cols, cex = 0.8,
       title = expression(paste(log[10], "(", kappa, ")")))
dev.off()

## 80, 120, 240, 280, 329, 973, 709

## coverage probabilities, estimates for a range of cutoffs
safeRange <- function(x) {
    if (length(x) == 0) {
        c(NA, NA)
    } else range(x)
}
##cutoffs <- minQuants[c("0.1%", "1%", "5%", "10%"), "100"]
cutoffs <- c(0.1, 0.05, 0.02, 0.01, 0.005)
poolTheta <- matrix(thetas[apply(poolMat, c(1,2), which.max)],
                    ncol = length(kseq))
poolIntervals <- lapply(cutoffs,
                        function(ct) apply(poolMat >= ct, c(1,2),
                                           function(x) {
                                               safeRange(thetas[which(x)])
                                           }))
poolInclude <- lapply(poolIntervals,
                      function(mat) mat[1,,] <= mn & mat[2,,] >= mn)
poolCovP <- lapply(poolInclude,
                   function(mat) apply(mat, 2, mean, na.rm = TRUE))
## using the mean
meanTheta <- apply(muMat, 1, mean)
meansd <- sqrt(apply(sdMat^2*nMat*(nMat-1), 1, sum)/(npop-1))/sqrt(npop)
meanInt <- cbind(meanTheta - qnorm(0.975)*meansd,
                 meanTheta + qnorm(0.975)*meansd)
meanCovP <- mean(meanInt[,1] <= mn & meanInt[,2] >= mn)
## using the minimum
minTheta <- sapply(1:nsim, function(ii) poolTheta[ii, minimax[ii]])
minInterval <- lapply(poolIntervals,
                      function(ints) {
                          sapply(1:nsim,
                                 function(ii) {
                                     ints[, ii, minimax[ii]]
                                 })})
minInclude <- lapply(minInterval,
                     function(mat) mat[1,] <= mn & mat[2,] >= mn)
minCovP <- lapply(minInclude, mean, na.rm = TRUE)

## plot the coverage probabilities by kappa
ctind <- 2
cutoff <- cutoffs[ctind]
png(paste0("meta", 100*cutoff, "pctCovP.png"), width = 3, height = 3,
    res = 480, units = "in")
tempCovP <- lowess(poolCovP[[ctind]], f = 1/6)$y
narrowPlot(xgrid = seq(-3, 3, by = 3), ygrid = seq(0.8, 1, by = 0.05),
           ylab = "Coverage probability", xlim = c(-3.5, 3.5),
           xlab = expression(paste(log[10], "(", kappa, ")")))
lines(log(kseq, 10), tempCovP)
polygon(c(log(kseq, 10), rev(log(kseq, 10))),
        c(tempCovP + qnorm(0.975)*sqrt(tempCovP*(1 - tempCovP)/nsim),
          rev(tempCovP - qnorm(0.975)*sqrt(tempCovP*(1 - tempCovP)/nsim))),
        col = adjustcolor("gray80", 0.5), border = NA)
abline(h = 1 - cutoff, lty = 2)
dev.off()

## plot estimates by kappa
png("metaEstimates.png", width = 3, height = 3, res = 480,
    units = "in")
narrowPlot(xgrid = seq(-3, 3, by = 1.5), ygrid = seq(-1, 1, by = 0.5),
           xlab = expression(paste(log[10], "(", kappa, ")")),
           xlim = c(-3.5, 3.5),
           ylab = "Error")
addQuantPoly(t(poolTheta), kseq = kseq, cex = 0.6)
dev.off()

## plot all intervals for a particular kappa
kapInd <- 41
pltMat <- poolIntervals[[ctind]][,,kapInd] #t(meanInt[order(meanInt[,1]),])
pltMat <- pltMat[, order(pltMat[1, ])]
plot(NA, xlim = range(pltMat, na.rm = TRUE), ylim = c(1, nsim),
     ylab = "Interval", xlab = "Bounds")
abline(h = seq(0, 1000, by = 200), v = seq(-1, 1, by = 0.5),
       col = "gray80")
for (ii in 1:nsim) lines(pltMat[,ii], rep(ii, 2),
                         col = adjustcolor("black", 0.5))

## middle line (kappa = 1) looks normal...
## can we fit a normal curve to it?
## suggested: nls() with c phi((theta - b)/a), a,b,c params

## look at the covid dataset
data(dat.axfors2021) ## requires estimating odds ratios
## or the difference in means data
data(dat.normand1999) ## differences in means but classic otherwise
