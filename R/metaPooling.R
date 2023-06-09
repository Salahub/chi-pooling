## FUNCTIONS #########################################################
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
paramSweep <- function(mus, sds, thetas, pfun = pnorm) {
    mus2 <- rep(mus, each = length(thetas))
    sds2 <- rep(sds, each = length(thetas))
    stdDs <- abs((mus2 - rep(thetas, times = length(mus)))/sds2)
    matrix(2*pfun(q = stdDs, lower.tail = FALSE),
           nrow = length(thetas))
}

## sweeping with kappa afterwards
kappaSweep <- function(ps, np, kseq = exp(seq(-8, 8, by = 0.1))) {
    sapply(kseq,
           function(kap) {
               apply(ps, 2,
                     function(p) pchisq(sum(qchisq(p, kap,
                                               lower.tail = FALSE)),
                                        np*kap,
                                        lower.tail = FALSE))
           })
}

## the case of inhomogeneous means (i.e. biased means)
inhomNormal <- function(npop, ngroup, mns = runif(-2, 2, ngroup),
                        sd = 1) {
    obs <- rep(mns, npop/ngroup) + rnorm(npop, mean = 0, sd = sd)
    grp <- rep(1:ngroup, npop/ngroup)
    cbind(x = obs, group = grp)
}
## the fixed effects model generation function
fixedNormal <- function(npop, ngroup, mn = 0, sd = 1) {
    obs <- rnorm(npop, mean = mn, sd = sd)
    grp <- sample(rep(1:ngroup, npop/ngroup))
    cbind(x = obs, group = grp)
}
## fixed normal with unbalanced data
fixedNormalUB <- function(npop, groups = sample(rep(1:8, npop/8)),
                          mn = 0, sd = 1) {
    obs <- rnorm(npop, mean = mn, sd = sd)
    grp <- groups
    cbind(x = obs, group = grp)
}
## the random effects model generation function
randomNormal <- function(npop, ngroup, mn = 0, mnsd = 1,
                         sd = 1) {
    mns <- rnorm(ngroup, mean = mn, sd = mnsd)
    obs <- rep(mns, npop/ngroup) + rnorm(npop, mean = 0, sd = sd)
    grp <- rep(1:ngroup, npop/ngroup)
    cbind(x = obs, group = grp)
}

## wrapper to generate meta-analysis results
simMetaStudies <- function(genFn, nsim = 1e3, npop = 240,
                           ngroup = 8) {
    muMat <- sdMat <- matrix(nrow = nsim, ncol = ngroup)
    nMat <- matrix(nrow = nsim, ncol = ngroup)
    for (ii in 1:nsim) {
        sim <- genFn(npop, ngroup) # generate data using genFun
        group <- sim[, "group"]
        obs <- sim[, "x"]
        nMat[ii,] <- table(group) # counts by group
        muMat[ii,] <- tapply(obs, group, mean) # means
        sdMat[ii,] <- tapply(obs, group, sd)/sqrt(table(group)) # sds
        if ((ii %% 10) == 0) cat("\r Generated ", ii, " of ", nsim)
    }
    list(means = muMat, sds = sdMat, ns = nMat)
}

## convert meta-analyses to p-values based on a probability function
metaToP <- function(pfun, studies, thetas) {
    simplify2array(lapply(thetas,
                          function(x) {
                              pfun(x, studies$means, studies$sds,
                                   studies$ns)
                          }))
}

## the t and normal pfuns
pfunNorm <- function(x, mus, sds, ns) {
    2*pnorm(-abs(mus - x)/sds)
}
pfunT <- function(x, mus, sds, ns) {
    2*pt(-abs(mus - x)/sds, df = ns - 1)
}

## convert p-values to chi pool estimates for a range of kappa
chiMetaSweep <- function(ps, kseq = exp(seq(-8, 8, by = 0.1))) {
    K <- dim(ps)[2]
    arr <- simplify2array(lapply(
        kseq,
        function(k) {
            pchisq(apply(qchisq(ps, k, lower.tail = FALSE),
                         c(1,3), sum),
                   df = K*k, lower.tail = FALSE)
        }))
    aperm(arr, c(1,3,2))
}

## weighted chi computation
chiWeighted <- function(ps, studies) {
    vs <- studies$sds
    K <- nrow(vs)
    wgts <- 1/vs^2
    t(sapply(1:nrow(wgts),
           function(ii) {
               pmat <- ps[ii,,]
               pchisq(apply(pmat, 2,
                            function(col) {
                                sum(qchisq(col,
                                           df = K*(wgts[ii,]/
                                                     sum(wgts[ii,])),
                                           lower.tail = FALSE))
                            }),
                      df = K,
                      lower.tail = FALSE)
           }))
}

## safe estimates for a range of cutoffs
safeRange <- function(x) {
    if (length(x) == 0) {
        c(NA, NA)
    } else range(x)
}

## compute coverage probabilities for a given sequence of chis
getChiEsts <- function(chiseqs, thetas, kseq, mn = 0,
                       cutoffs = c(0.1, 0.05, 0.02, 0.01)) {
    ests <- matrix(thetas[apply(chiseqs, c(1,2),
                                     which.max)],
                     ncol = length(kseq))
    intervals <- lapply(cutoffs,
                        function(ct) {
                            apply(chiseqs >= ct,
                                  c(1,2),
                                  function(x) {
                                      safeRange(thetas[which(x)])
                                  })
                        })
    include <- lapply(intervals,
                      function(mat) mat[1,,] <= mn & mat[2,,] >= mn)
    covP <- lapply(include,
                   function(mat) apply(mat, 2, mean, na.rm = TRUE))
    list(ests = ests, intervals = intervals, covP = covP)
}

## and something for univariate estimates
getUniEsts <- function(seqs, thetas, mn = 0,
                       cutoffs = c(0.1, 0.05, 0.02, 0.01)) {
    ests <- thetas[apply(seqs, 1, which.max)]
    intervals <- lapply(cutoffs,
                           function(ct) {
                               apply(seqs >= ct,
                                     1,
                                     function(x) {
                                         safeRange(thetas[which(x)])
                                     })
                           })
    include <- lapply(intervals,
                         function(mat) mat[1,] <= mn & mat[2,] >= mn)
    covP <- lapply(include, mean, na.rm = TRUE)
    list(ests = ests, intervals = intervals, covP = covP)
}

## the mean estimate function
getMeanEsts <- function(sim, cutoffs = c(0.1, 0.05, 0.02, 0.01)) {
    mns <- sim$means
    wgts <- 1/sim$sds^2
    wgted <- mns*(wgts)
    ests <- rowSums(wgted)/rowSums(wgts)
    se <- sqrt(1/rowSums(wgts))
    intervals <- lapply(cutoffs,
                        function(a) {
                            matrix(rep(ests, each = 2) +
                                   c(-1, 1)*qnorm(1-a)*rep(se, each = 2),
                                   nrow = 2)
                        })
    covP <- lapply(intervals,
                   function(mat) mean(mat[1,] <= mn & mat[2,] >= mn,
                                      na.rm = TRUE))
    list(ests = ests, intervals = intervals, covP = covP)
}

## plot a realization based on the simulated data
plotRealization <- function(chi, minchi, means, meanEst, thetas, kseq,
                            ind = sample(1:nsim, 1),
                            cols = 1:4, kaps = c(1, 41, 81),
                            refKap = 41) {
    lines(thetas,  minchi[ind, ], type = 'l',
          col = adjustcolor(cols[4], 0.8))
    for (ii in kaps) lines(thetas, chi[ind,ii,], type = 'l',
                           col = if (ii <= refKap-1) {
                                     adjustcolor(cols[1], 0.8)
                                 } else if (ii == refKap) {
                                     adjustcolor(cols[2], 0.8)
                                 } else adjustcolor(cols[3], 0.8))
    abline(v = 0, col = adjustcolor("black", 0.8), lty = 2, lwd = 1)
    abline(v = meanEst[ind], col = adjustcolor("black", 0.6), lwd = 1)
    abline(v = means[ind,], col = adjustcolor("gray50", 0.4))
    abline(h = 0.05, lty = 3)
    legend(x = "topleft", legend = c(round(log(kseq[kaps], 10),
                                           1), "Min"),
           lty = 1, col = cols, cex = 0.8,
           title = expression(paste(log[10], "(", kappa, ")")))
}


## SIMULATIONS #######################################################
## settings
mn <- 0
std <- 2
nsim <- 1000
ngroup <- 8
npop <- 240
minQuants <- readRDS("curveMinQuantiles.Rds")
kseq <- exp(seq(-8, 8, by = 0.1))
thetas <- seq(-1.5, 1.5, by = 0.01)
cols <- RColorBrewer::brewer.pal(4, "Dark2")
cutoffs <- seq(0.20, 0.01, by = -0.01)

## FIXED EFFECTS ##
## compute all the features
fixedGen <- function(np, ng) fixedNormal(np, ng,
                                         sd = std)
set.seed(9381278)
fixedSim <- simMetaStudies(fixedGen, nsim = nsim, npop = npop,
                           ngroup = ngroup)
fixedps <- metaToP(pfunT, fixedSim, thetas = thetas)
fixedChis <- chiMetaSweep(fixedps, kseq = kseq)
fixedMnChi <- apply(fixedChis, c(1,3), min)
fixedwgtChi <- chiWeighted(fixedps, fixedSim)
## get intervals
fixedChiInt <- getChiEsts(fixedChis, thetas = thetas, kseq = kseq,
                          cutoffs = cutoffs)
fixedMinInt <- getUniEsts(fixedMnChi, thetas = thetas,
                          cutoffs = cutoffs)
fixedWgtInt <- getUniEsts(fixedwgtChi, thetas = thetas,
                          cutoffs = cutoffs)
fixedMnInt <- getMeanEsts(fixedSim, cutoffs = cutoffs)

## plot estimates by kappa
png("metaEstFixed.png", width = 3, height = 3, res = 480,
    units = "in")
narrowPlot(xgrid = seq(-3, 3, by = 1.5),
           ygrid = seq(-0.5, 0.5, by = 0.25),
           xlab = expression(paste(log[10], "(", kappa, ")")),
           xlim = c(-3.5, 3.5), ylim = c(-0.55, 0.55),
           ylab = expression(hat(theta)^{(E)}),
           mars = c(2.1, 2.2, 1.1, 1.1))
addQuantPoly(t(fixedChiInt$ests), kseq = kseq, cex = 0.6)
lines(log(kseq, 10), colMeans(fixedChiInt$ests, na.rm = TRUE),
      lty = 2, col = "firebrick")
abline(h = quantile(fixedMnInt$ests, c(0.025, 0.25, 0.75, 0.975)),
       lty = 3)
mtext(c("0.95", "0.5"), at = quantile(fixedMnInt$ests,
                                      c(0.025, 0.25)),
      side = 4, cex = 0.6, las = 1)
dev.off()

## compare mean estimates to kappa estimates

## plot estimates by kappa
kapInd <- 81
png("metaEstvsMeanFixed.png", width = 3, height = 3, res = 480,
    units = "in")
narrowPlot(xgrid = seq(-0.5, 0.5, by = 0.25),
           ygrid = seq(-0.5, 0.5, by = 0.25),
           xlab = expression(hat(theta)^{"("~E~")"}),
           ylab = expression(hat(theta)))
points(fixedChiInt$est[, kapInd], fixedMnInt$ests)
dev.off()

## plot the coverage probabilities by kappa
ctind <- 20
cutoff <- cutoffs[ctind]
png(paste0("meta", 100*cutoff, "pctCovPFix.png"), width = 2.5,
    height = 2.5, res = 480, units = "in")
tempCovP <- lowess(fixedChiInt$covP[[ctind]], f = 1/6)$y
tempCovP <- fixedChiInt$covP[[ctind]]
narrowPlot(xgrid = seq(-3, 3, by = 3), ygrid = seq(0.8, 1, by = 0.05),
           ylab = "Coverage probability", xlim = c(-3.5, 3.5),
           xlab = expression(paste(log[10], "(", kappa, ")")))
lines(log(kseq, 10), tempCovP)
points(log(kseq, 10), tempCovP, cex = 0.5)
polygon(c(log(kseq, 10), rev(log(kseq, 10))),
        c(tempCovP + qnorm(0.975)*sqrt(tempCovP*(1 - tempCovP)/nsim),
          rev(tempCovP - qnorm(0.975)*sqrt(tempCovP*(1 - tempCovP)/nsim))),
        col = adjustcolor("gray80", 0.5), border = NA)
abline(h = 1 - cutoff, lty = 2)
abline(v = log(2,10))
abline(h = 1 - cutoff/2, lty = 2, col = "firebrick")
dev.off()

## plot all intervals for a particular kappa
ctind <- 20
cutoff <- cutoffs[ctind]
kapInd <- 88 # 88  for Fisher's
pltMat <- fixedChiInt$intervals[[ctind]][,,kapInd]
##pltMat <- t(fixedMnInt[order(fixedMnInt[,1]),])
pltMat <- pltMat[, order(pltMat[1, ])]
png(paste0("poolInts", 100*cutoff, "pctFixed.png"), height = 2.5,
    width = 2.5, res = 480, units = "in")
narrowPlot(xgrid = seq(-1, 1, by = 0.5), ygrid = seq(0, 1000, by = 200),
           ylab = "Interval", xlab = "Bounds")
abline(v = 0)
for (ii in 1:nsim) lines(pltMat[,ii], rep(ii, 2),
                         col = adjustcolor("black", 0.2))
abline(h = nsim*(1 - cutoff), lty = 2)
abline(h = nsim*(1 - cutoff/2), lty = 2, col = "firebrick")
dev.off()

## check level of implicit test
kapInd <- 88
levels <- lapply(fixedChiInt$intervals,
                 function(el) apply(el, 3,
                                    function(mat) mean(is.na(mat))))
kapLevs <- sapply(levels, function(vec) vec[kapInd])
png("metaRejectPFixed.png", width = 3, height = 3, res = 480,
    units = "in")
narrowPlot(xgrid = seq(0, 0.2, by = 0.05),
           ygrid = seq(0, 0.2, by = 0.05),
           xlab = "a", ylab = expression(alpha))
points(cutoffs, kapLevs)
for (ii in 1:length(cutoffs)) {
    lines(rep(cutoffs[ii], 2),
          (kapLevs[ii] +
           c(-1,1)*1.96/sqrt(nsim)*sqrt(kapLevs[ii]*(1 - kapLevs[ii]))),
          col = "gray50")
}
abline(a = 0, b = 1)
abline(a = 0, b = 0.5, lty = 3)
dev.off()

## check coverage probabilities for a particular kappa
kapInd <- 88
png("metaCovPFixed.png", width = 3, height = 3, res = 480,
    units = "in")
covPs <- sapply(fixedChiInt$covP, function(el) el[kapInd])
narrowPlot(xgrid = seq(0, 0.2, by = 0.05),
           ygrid = seq(0, 0.2, by = 0.05),
           xlab = "a", ylab = expression("2"~{(1-pi)}))
points(cutoffs, 2*(1 - covPs))
for (ii in 1:length(cutoffs)) {
    lines(rep(cutoffs[ii], 2),
          2*(1 - covPs[ii] +
             c(-1,1)*1.96/sqrt(nsim)*sqrt(covPs[ii]*(1 - covPs[ii]))),
          col = "gray50")
}
abline(a = 0, b = 1)
dev.off()

## plot a realization
real <- 55 # sample(1:nsim, 1)
png("metaPoolCurvesFixed.png", width = 5, height = 3, res = 480,
    units = "in")
narrowPlot(xgrid = seq(-1, 1, by = 0.5),
           ygrid = seq(0, 0.6, by = 0.15),
           xlab = expression(x),
           ylab = expression(paste("chi(", bold(p), "(x)",
                                   ";", kappa, ")")),
           addGrid = FALSE)
plotRealization(chi = fixedChis, minchi = fixedMnChi,
                means = fixedSim$means, meanEst = fixedMnInt$ests,
                thetas = thetas, kseq = kseq, cols = cols,
                kaps = c(1, 88, 161), ind = real, refKap = 88)
dev.off()

## reject homogeneity: 19
## separation of curves: 55, 79, 86


## FIXED EFFECTS UNBALANCED ##
fixedGenUB <- function(np, ng) {
    grps <- sample(1:ng, size = np, replace = TRUE,
                   prob = c(3, 3, 6, 2, 2, 2, 3, 4))
    fixedNormalUB(np, groups = grps, sd = std)
}
## generate and compute on studies
set.seed(73183)
ubSim <- simMetaStudies(fixedGenUB, nsim = nsim, npop = npop,
                        ngroup = ngroup)
ubps <- metaToP(pfunT, ubSim, thetas = thetas)
ubChis <- chiMetaSweep(ubps, kseq = kseq)
ubMnChi <- apply(ubChis, c(1,3), min)
ubwgtChi <- chiWeighted(ubps, ubSim)
## get intervals
ubChiInt <- getChiEsts(ubChis, thetas = thetas, kseq = kseq,
                       cutoffs = cutoffs)
ubMinInt <- getUniEsts(ubMnChi, thetas = thetas,
                       cutoffs = cutoffs)
ubWgtInt <- getUniEsts(ubwgtChi, thetas = thetas,
                       cutoffs = cutoffs)
ubMnInt <- getMeanEsts(ubSim, cutoffs = cutoffs)

## group size histogram
png("metaGroupHist.png", width = 3, height = 3, res = 480,
    unit = "in")
par(mar = c(3.1, 3, 1.1, 0.1))
hist(ubSim$ns, main = "", cex.axis = 0.8, xlab = "Group size",
     ylab = "Frequency", cex.lab = 0.8)
mtext("Group size", side =1, line = 2, cex = 0.8)
mtext("Frequency", side =2, line = 2, cex = 0.8)
dev.off()

## plot estimates by kappa
png("metaEstUB.png", width = 3, height = 3, res = 480,
    units = "in")
narrowPlot(xgrid = seq(-3, 3, by = 1.5),
           ygrid = seq(-0.5, 0.5, by = 0.25),
           xlab = expression(paste(log[10], "(", kappa, ")")),
           xlim = c(-3.5, 3.5), ylim = c(-0.55, 0.55),
           ylab = expression(hat(theta)^{(E)}),
           mar = c(2.1, 2.2, 1.1, 1.1))
addQuantPoly(t(ubChiInt$ests), kseq = kseq, cex = 0.6)
lines(log(kseq, 10), colMeans(ubChiInt$ests, na.rm = TRUE),
      lty = 2, col = "firebrick")
abline(h = quantile(ubMnInt$ests, c(0.025, 0.25, 0.75, 0.975)),
       lty = 3)
mtext(c("0.95", "0.5"), at = quantile(ubMnInt$ests,
                                      c(0.025, 0.25)),
      side = 4, cex = 0.6, las = 1)
#abline(h = quantile(ubWgtInt$ests, c(0.025, 0.25, 0.75, 0.975)),
#       lty = 3, col = "firebrick")
dev.off()

## plot estimates by kappa
kapInd <- 81
png("metaEstvsMeanUB.png", width = 3, height = 3, res = 480,
    units = "in")
narrowPlot(xgrid = seq(-0.5, 0.5, by = 0.25),
           ygrid = seq(-0.5, 0.5, by = 0.25),
           xlab = expression(hat(theta)^{"("~E~")"}),
           ylab = expression(hat(theta)))
points(ubChiInt$est[, kapInd], ubMnInt$ests)
dev.off()

## plot the coverage probabilities by kappa
ctind <- 20
cutoff <- cutoffs[ctind]
png(paste0("meta", 100*cutoff, "pctCovPUB.png"), width = 2.5,
    height = 2.5, res = 480, units = "in")
tempCovP <- ubChiInt$covP[[ctind]]
narrowPlot(xgrid = seq(-3, 3, by = 3), ygrid = seq(0.8, 1, by = 0.05),
           ylab = "Coverage probability", xlim = c(-3.5, 3.5),
           xlab = expression(paste(log[10], "(", kappa, ")")))
lines(log(kseq, 10), tempCovP)
points(log(kseq, 10), tempCovP, cex = 0.5)
polygon(c(log(kseq, 10), rev(log(kseq, 10))),
        c(tempCovP + qnorm(0.975)*sqrt(tempCovP*(1 - tempCovP)/nsim),
          rev(tempCovP - qnorm(0.975)*sqrt(tempCovP*(1 - tempCovP)/nsim))),
        col = adjustcolor("gray80", 0.5), border = NA)
abline(h = 1 - cutoff, lty = 2)
abline(v = log(2,10))
abline(h = 1 - cutoff/2, lty = 2, col = "firebrick")
dev.off()

## plot all intervals for a particular kappa
ctind <- 20
cutoff <- cutoffs[ctind]
kapInd <- 88 # 88  for Fisher's
pltMat <- ubChiInt$intervals[[ctind]][,,kapInd]
##pltMat <- ubWgtInt$intervals[[ctind]]
##pltMat <- t(fixedInts$meanInt[order(fixedInts$meanInt[,1]),])
pltMat <- pltMat[, order(pltMat[1, ])]
png(paste0("poolInts", 100*cutoff, "pctUB.png"), height = 2.5,
    width = 2.5, res = 480, units = "in")
narrowPlot(xgrid = seq(-1, 1, by = 0.5), ygrid = seq(0, 1000, by = 200),
           ylab = "Interval", xlab = "Bounds")
abline(v = 0)
for (ii in 1:nsim) lines(pltMat[,ii], rep(ii, 2),
                         col = adjustcolor("black", 0.2))
abline(h = nsim*(1 - cutoff), lty = 2)
abline(h = nsim*(1 - cutoff/2), lty = 2, col = "firebrick")
dev.off()

## check level of implicit test
kapInd <- 88
levels <- lapply(ubChiInt$intervals,
                 function(el) apply(el, 3,
                                    function(mat) mean(is.na(mat))))
kapLevs <- sapply(levels, function(vec) vec[kapInd])
png("metaRejectPUB.png", width = 3, height = 3, res = 480,
    units = "in")
narrowPlot(xgrid = seq(0, 0.2, by = 0.05),
           ygrid = seq(0, 0.2, by = 0.05),
           xlab = "a", ylab = expression(alpha))
points(cutoffs, kapLevs)
for (ii in 1:length(cutoffs)) {
    lines(rep(cutoffs[ii], 2),
          (kapLevs[ii] +
           c(-1,1)*1.96/sqrt(nsim)*sqrt(kapLevs[ii]*(1 - kapLevs[ii]))),
          col = "gray50")
}
abline(a = 0, b = 1)
abline(a = 0, b = 0.5, lty = 3)
dev.off()

## check coverage probabilities for a particular kappa
kapInd <- 88
png("metaCovPUB.png", width = 3, height = 3, res = 480,
    units = "in")
covPs <- sapply(ubChiInt$covP, function(el) el[kapInd])
narrowPlot(xgrid = seq(0, 0.2, by = 0.05),
           ygrid = seq(0, 0.2, by = 0.05),
           xlab = "a", ylab = expression("2"~{(1-pi)}))
points(cutoffs, 2*(1 - covPs))
for (ii in 1:length(cutoffs)) {
    lines(rep(cutoffs[ii], 2),
          2*(1 - covPs[ii] +
             c(-1,1)*1.96/sqrt(nsim)*sqrt(covPs[ii]*(1 - covPs[ii]))),
          col = "gray50")
}
abline(a = 0, b = 1)
dev.off()

## plot a realization
real <- 86 # sample(1:nsim, 1)
png("metaPoolCurvesUB.png", width = 5, height = 3, res = 480,
    units = "in")
narrowPlot(xgrid = seq(-1.5, 1.5, by = 0.5),
           ygrid = seq(0, 1, by = 0.2),
           xlab = expression(x),
           ylab = expression(paste("chi(", bold(p), ";", kappa, ")")),
           addGrid = FALSE)
plotRealization(chi = ubChis, minchi = ubMnChi,
                means = ubSim$means, meanEst = ubMnInt$ests,
                thetas = thetas, kseq = kseq, cols = cols,
                kaps = c(1, 88, 161), ind = real, refKap = 88)
dev.off()

## 15


## RANDOM EFFECTS ##
mnsd <- 1
randGen <- function(np, ng) randomNormal(np, ng, mnsd = mnsd,
                                         sd = sqrt(std^2 - mnsd^2))
set.seed(8251506)
randSim <- simMetaStudies(randGen, nsim = nsim, npop = npop,
                          ngroup = ngroup)
randps <- metaToP(pfunT, randSim, thetas = thetas)
randChis <- chiMetaSweep(randps, kseq = kseq)
## get intervals
randChiInt <- getChiEsts(randChis, thetas = thetas, kseq = kseq,
                         cutoffs = cutoffs)
randMnInt <- getMeanEsts(randSim, cutoffs = cutoffs)

## plot estimates by kappa
png("metaEstRandom.png", width = 3, height = 3, res = 480,
    units = "in")
narrowPlot(xgrid = seq(-3, 3, by = 1.5),
           ygrid = seq(-1.5, 1.5, by = 0.75),
           xlab = expression(paste(log[10], "(", kappa, ")")),
           xlim = c(-3.5, 3.5), #ylim = c(-0.55, 0.55),
           ylab = "Error")
addQuantPoly(t(randChiInt$ests), kseq = kseq, cex = 0.6)
lines(log(kseq, 10), colMeans(randChiInt$ests, na.rm = TRUE),
      lty = 2)
abline(h = quantile(randMnInt$ests, c(0.025, 0.25, 0.75, 0.975)),
       lty = 3)
mtext(c("0.95", "0.5"),
      at = quantile(randMnInt$ests, c(0.025, 0.25)),
      side = 4, cex = 0.6, las = 1)
dev.off()

## plot all intervals for a particular kappa
kapInd <- 88 # 88  for Fisher's
ctind <- 16
cutoff <- cutoffs[ctind]
pltMat <- randChiInt$intervals[[ctind]][,,kapInd]
##pltMat <- t(fixedInts$meanInt[order(fixedInts$meanInt[,1]),])
pltMat <- pltMat[, order(pltMat[1, ])]
png(paste0("poolInts", 100*cutoff, "pctRand.png"), height = 3,
    width = 3, res = 480, units = "in")
narrowPlot(xgrid = seq(-1.5, 1.5, by = 0.75),
           ygrid = seq(0, 1000, by = 200),
           ylab = "Interval", xlab = "Bounds")
abline(v = 0)
for (ii in 1:nsim) lines(pltMat[,ii], rep(ii, 2),
                         col = adjustcolor("black", 0.2))
abline(h = nsim*(1 - cutoff), lty = 2)
dev.off()

## 0.975 is the proportion of empty intervals
## check for a range of sds
mnsdSeq <- seq(0.1, 1.9, by = 0.1)
## a list of closures with the correct variances
randGens <- lapply(mnsdSeq,
                   function(mnsd) {
                       function(np, ng) {
                           randomNormal(np, ng,
                                        mnsd = mnsd,
                                        sd = sqrt(std^2 - mnsd^2))
                       }})
## simulate for each
set.seed(54910913)
randSims <- lapply(randGens, simMetaStudies, nsim = nsim, npop = npop,
                   ngroup = ngroup)
randps <- lapply(randSims, metaToP, pfun = pfunT, thetas = thetas)
randChis <- lapply(randps, chiMetaSweep, kseq = kseq[c(1,81,88,161)])
## get intervals
randChiInt <- lapply(randChis, getChiEsts, thetas = thetas,
                     kseq = kseq[c(1,81,88,161)], cutoffs = cutoffs)
## and rejection probabilities
randRejectP <- lapply(randChiInt,
                      function(ints) {
                          lapply(ints$intervals,
                                 function(mat) {
                                     colMeans(is.na(mat[1,,]))
                                 })})

## for rejection rule >= 0.05, what does the power look like in tau?
powers5cv <- sapply(randRejectP, function(lst) lst[[16]][2])
powers5fis <- sapply(randRejectP, function(lst) lst[[16]][3])
powers5sto <- sapply(randRejectP, function(lst) lst[[16]][4])
powers5tip <- sapply(randRejectP, function(lst) lst[[16]][1])
png(paste0("poolPowersRand.png"), height = 3.5,
    width = 3, res = 480, units = "in")
narrowPlot(xgrid = seq(0, 1, by = 0.25), ygrid = seq(0, 1, by = 0.25),
           ylab = "Power", mars = c(3.5, 2.1, 1.1, 1.1))
mtext(expression(frac(tau^2, tau^2 + s^2)), side = 1, line = 2.5,
      cex = 0.8)
points(mnsdSeq^2/4, powers5cv, col = cols[4], cex = 0.8)
lines(mnsdSeq^2/4, powers5cv, col = cols[4])
points(mnsdSeq^2/4, powers5fis, col = cols[2], cex = 0.8)
lines(mnsdSeq^2/4, powers5fis, col = cols[2])
points(mnsdSeq^2/4, powers5sto, col = cols[3], cex = 0.8)
lines(mnsdSeq^2/4, powers5sto, col = cols[3])
points(mnsdSeq^2/4, powers5tip, col = cols[1], cex = 0.8)
lines(mnsdSeq^2/4, powers5tip, col = cols[1])
legend(x = "bottomright", legend = c("-3.5", "0", "0.3", "3.5"),
       lty = 1, pch = 1, col = cols[c(1,4,2,3)], cex = 0.8,
       pt.cex = 0.8,
       title = expression(paste(log[10], "(", kappa, ")")))
dev.off()

## REAL DATA #########################################################
library(metadat)
library(metafor)

## look at the covid dataset
data(dat.axfors2021) ## requires estimating odds ratios
## get log odds ratios
ors <- escalc(measure="OR", ai=hcq_arm_event, n1i=hcq_arm_total,
              ci=control_arm_event, n2i=control_arm_total,
              data=dat.axfors2021)
## filter out to consider only small doses of HCQ
ors <- ors[ors$hcq_cq == "hcq" & ors$high_dose == "no", ]
## sweep log odds ratio values with kappa = 2
xseq <- seq(-5, 5, by = 0.05)
curve <- apply(qchisq(paramSweep(ors$yi, sqrt(ors$vi), thetas = xseq),
                      2, lower.tail = FALSE), 1,
               function(row) pchisq(sum(row), 2*26,
                                    lower.tail = FALSE))
curveWgtd <- apply(paramSweep(ors$yi, sqrt(ors$vi), thetas = xseq),
                   1, function(row) {
                       pchisq(sum(qchisq(row, 1/ors$vi,
                                     lower.tail = FALSE)),
                              df = sum(1/ors$vi),
                              lower.tail = FALSE)
                       })

plot(xseq, curveWgtd, type = 'l')
abline(v = ors$yi, col = adjustcolor("gray50", 0.5))
meanEst <- sum(ors$yi*1/ors$vi)/sum(1/ors$vi)
meanSE <- sqrt(1/sum(1/ors$vi))

## school calendar data
data(dat.konstantopoulos2011)
schoolDat <- dat.konstantopoulos2011
schoolDat$district <- factor(schoolDat$district)
rm(dat.konstantopoulos2011)

## plot the data
schoolCol <- RColorBrewer::brewer.pal(11, "Set3")
png("metaSchoolMeans.png", width = 3, height = 3, units = "in",
    res = 480)
narrowPlot(ygrid = seq(-1, 1, by = 0.5),
           xgrid = seq(0, 12, by = 3),
           ylab = "Mean difference", xlab = "District")
points(unclass(as.factor(schoolDat$district)),
       y = schoolDat$yi, pch = 20,
       col = schoolCol[unclass(schoolDat$district)])
for (ii in 1:nrow(schoolDat)) {
    lines(rep(unclass(schoolDat$district)[ii], 2),
          rep(schoolDat$yi[ii],2) +
          c(-1.96, 1.96)*sqrt(schoolDat$vi[ii]),
          col = schoolCol[unclass(schoolDat$district)[ii]])
}
dev.off()

## sweep of values
xseq <- seq(-2, 2, by = 0.01)
pvals <- apply(qchisq(paramSweep(schoolDat$yi,
                                 sqrt(schoolDat$vi),
                                 thetas = xseq),
                      2, lower.tail = FALSE), 1,
               function(row) pchisq(sum(row), 2*nrow(schoolDat),
                                    lower.tail = FALSE))
## strong evidence of heterogeneity
## pool by district
schoolSplit <- split(schoolDat[, c("yi", "vi")],
                     schoolDat$district)
## pool for each district
districtPools <- lapply(schoolSplit,
                        function(sdat) {
                            apply(qchisq(paramSweep(sdat$yi,
                                 sqrt(sdat$vi),
                                 thetas = xseq),
                      2, lower.tail = FALSE), 1,
               function(row) pchisq(sum(row), 2*nrow(sdat),
                                    lower.tail = FALSE))
                            })
## do this one...super easy analysis

## plot the data
schoolCol <- RColorBrewer::brewer.pal(11, "Set3")
png("metaSchoolIntervals.png", width = 3, height = 3, units = "in",
    res = 480)
narrowPlot(ygrid = seq(-1, 1, by = 0.5),
           xgrid = seq(0, 12, by = 3),
           ylab = "Mean difference", xlab = "District")
points(unclass(as.factor(schoolDat$district)),
       y = schoolDat$yi, pch = 20,
       col = schoolCol[unclass(schoolDat$district)])
for (ii in 1:nrow(schoolDat)) {
    lines(rep(unclass(schoolDat$district)[ii], 2),
          rep(schoolDat$yi[ii],2) +
          c(-1.96, 1.96)*sqrt(schoolDat$vi[ii]),
          col = schoolCol[unclass(schoolDat$district)[ii]])
}
for (ii in 1:length(districtPools)) {
    temp <-  safeRange(xseq[(districtPools[[ii]] >= 0.05)])
    tempch <- if (is.na(temp)[1]) 4 else 0
    points(y = xseq[which.max(districtPools[[ii]])],
           x = ii, pch = tempch)
    lines(x = rep(ii, 2), y = temp)
}
dev.off()
