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

## clean up code by defining the simulation loop here
simEvidentialRegion <- function(genFn, nsim = 1e3, npop = 104,
                                ngroup = 8,
                                kseq = seq(-8, 8, by = 0.2),
                                thetas = seq(-2.5, 2.5, by = 0.01)) {
    muMat <- sdMat <- matrix(nrow = nsim, ncol = ngroup)
    nMat <- matrix(nrow = nsim, ncol = ngroup)
    pMat <- array(dim = c(nsim, ngroup, length(thetas)))
    poolMat <- array(dim = c(nsim, length(kseq), length(thetas)))
    pwgtMat <- matrix(nrow = nsim, ncol = length(thetas))
    pfun <- function(q, ...) pnorm(q, ...)
    for (ii in 1:nsim) {
        sim <- genFn(npop, ngroup) # generate data using genFun
        group <- sim[, "group"]
        obs <- sim[, "x"]
        nMat[ii,] <- table(group) # counts by group
        muMat[ii,] <- tapply(obs, group, mean) # means
        sdMat[ii,] <- tapply(obs, group, sd)/sqrt(table(group)) # sds
        pMat[ii,,] <- simps <- t(paramSweep(muMat[ii,], sdMat[ii,],
                                            thetas,
                                            pfun = pfun)) # p-vals
        poolMat[ii,,] <- t(kappaSweep(simps, np = ngroup,
                                      kseq = kseq)) # pooled vals
        pwgtMat[ii,] <- apply(simps,
                              2,
                              function(ps) {
                                  pchisq(sum(qchisq(ps,
                                                    df = 1/sdMat[ii,]^2,
                                      lower.tail = FALSE)),
                                      df = sum(1/sdMat[ii,]^2),
                                      lower.tail = FALSE)
                              })
        if ((ii %% 10) == 0) cat("\r Done ", ii, " of ", nsim)
    }
    ## compute the mini-max kappa for each repetition
    minimax <- apply(apply(poolMat, c(1,2), max), 1, which.min)
    ## and the minimum over kappa for each repetition
    minpool <- apply(poolMat, c(1,3), min)
    ## return everything
    list(muMat = muMat, sdMat = sdMat,
         nMat = nMat, pMat = pMat, poolMat = poolMat,
         minimax = minimax, minKappa = minpool,
         wgtdPool= pwgtMat)
}

## safe estimates for a range of cutoffs
safeRange <- function(x) {
    if (length(x) == 0) {
        c(NA, NA)
    } else range(x)
}

## use a simulation and compute coverage probabilities, etc.
getCoverage <- function(sim, thetas, kseq,
                        cutoffs = c(0.1, 0.05, 0.02, 0.01, 0.005)) {
    poolTheta <- matrix(thetas[apply(sim$poolMat, c(1,2),
                                     which.max)],
                        ncol = length(kseq))
    poolIntervals <- lapply(cutoffs,
                            function(ct) {
                                apply(sim$poolMat >= ct,
                                      c(1,2),
                                      function(x) {
                                          safeRange(thetas[which(x)])
                                      })
                            })
    poolInclude <- lapply(poolIntervals,
                          function(mat) mat[1,,] <= mn & mat[2,,] >= mn)
    poolCovP <- lapply(poolInclude,
                       function(mat) apply(mat, 2, mean, na.rm = TRUE))
    ## for the weighted mean
    wgtTheta <- thetas[apply(sim$wgtdPool, 1, which.max)]
    wgtIntervals <- lapply(cutoffs,
                           function(ct) {
                               apply(sim$wgtdPool >= ct,
                                     1,
                                     function(x) {
                                         safeRange(thetas[which(x)])
                                     })
                           })
    wgtInclude <- lapply(wgtIntervals,
                         function(mat) mat[1,] <= mn & mat[2,] >= mn)
    wgtCovP <- lapply(wgtInclude, mean)
    ## using the mean and the classic combination function
    meanTheta <- apply(sim$muMat*sim$sdMat^(-2), 1, sum)/
        rowSums(sim$sdMat^(-2))
    meansd <- sqrt(1/rowSums(sim$sdMat^(-2)))
    #meansd <- sqrt(apply(sdMat^2*nMat*(nMat-1), 1,
    #                     sum)/(npop-1))/sqrt(npop)
    meanInt <- cbind(meanTheta - qnorm(0.975)*meansd,
                     meanTheta + qnorm(0.975)*meansd)
    meanCovP <- mean(meanInt[,1] <= mn & meanInt[,2] >= mn)
    ## the minimum theta
    minTheta <- thetas[apply(sim$minKappa, 1, which.max)]
    minInt <-  lapply(cutoffs,
                      function(ct) {
                          apply(sim$minKappa >= ct,
                                1,
                                function(x) {
                                    safeRange(thetas[which(x)])
                                })
                      })
    minCovP <- lapply(minInt,
                      function(mat) mean(mat[1,] <= mn &
                                         mat[2,] >= mn,
                                         na.rm = TRUE))
    ## return everything
    list(poolTheta = poolTheta, poolIntervals = poolIntervals,
         poolInclude = poolInclude, poolCovP = poolCovP,
         meanTheta = meanTheta, meansd = meansd, meanInt = meanInt,
         meanCovP = meanCovP, minTheta = minTheta, minInt = minInt,
         minCovP = minCovP, wgtTheta = wgtTheta,
         wgtIntervals = wgtIntervals, wgtInclude = wgtInclude,
         wgtCovP = wgtCovP)
}

## plot a realization based on the simulated data
plotRealization <- function(sims, thetas, kseq,
                            ind = sample(1:nsim, 1),
                            cols = 1:4, kaps = c(1, 41, 81),
                            refKap = 41) {
    lines(thetas,  sims$minKappa[ind, ], type = 'l',
          col = adjustcolor(cols[4], 0.8))
    for (ii in kaps) lines(thetas, sims$poolMat[ind,ii,], type = 'l',
                           col = if (ii <= refKap-1) {
                                     adjustcolor(cols[1], 0.8)
                                 } else if (ii == refKap) {
                                     adjustcolor(cols[2], 0.8)
                                 } else adjustcolor(cols[3], 0.8))
    abline(v = 0, col = adjustcolor("black", 0.8), lty = 2, lwd = 1)
    abline(v = sum(sims$muMat[ind,]*sims$sdMat[ind,]^(-2))/
           sum(sims$sdMat[ind,]^(-2)),
           col = adjustcolor("black", 0.6), lwd = 1)
    abline(v = sims$muMat[ind,], col = adjustcolor("gray50", 0.4))
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
fixedGen <- function(np, ng) fixedNormal(np, ng,
                                         sd = std)
set.seed(9381278)
fixedSim <- simEvidentialRegion(fixedGen, nsim = nsim, npop = npop,
                                ngroup = ngroup, kseq = kseq,
                                thetas = thetas)
## compute intervals and coverage probabilities
fixedInts <- getCoverage(fixedSim, thetas = thetas, kseq = kseq,
                         cutoffs = cutoffs)

## plot estimates by kappa
png("metaEstFixed.png", width = 3, height = 3, res = 480,
    units = "in")
narrowPlot(xgrid = seq(-3, 3, by = 1.5),
           ygrid = seq(-0.5, 0.5, by = 0.25),
           xlab = expression(paste(log[10], "(", kappa, ")")),
           xlim = c(-3.5, 3.5), ylim = c(-0.55, 0.55),
           ylab = "Error")
addQuantPoly(t(fixedInts$poolTheta), kseq = kseq, cex = 0.6)
lines(log(kseq, 10), colMeans(fixedInts$poolTheta, na.rm = TRUE),
      lty = 2, col = "firebrick")
abline(h = quantile(fixedInts$meanTheta, c(0.025, 0.25, 0.75, 0.975)),
       lty = 3)
dev.off()

## plot the coverage probabilities by kappa
ctind <- 20
cutoff <- cutoffs[ctind]
png(paste0("meta", 100*cutoff, "pctCovPFix.png"), width = 2.5,
    height = 2.5, res = 480, units = "in")
tempCovP <- lowess(fixedInts$poolCovP[[ctind]], f = 1/6)$y
tempCovP <- fixedInts$poolCovP[[ctind]]
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
ctind <- 16
cutoff <- cutoffs[ctind]
kapInd <- 88 # 88  for Fisher's
pltMat <- fixedInts$poolIntervals[[ctind]][,,kapInd]
##pltMat <- t(fixedInts$meanInt[order(fixedInts$meanInt[,1]),])
pltMat <- pltMat[, order(pltMat[1, ])]
png(paste0("poolInts", 100*cutoff, "pctFixed.png"), height = 2.5,
    width = 2.5, res = 480, units = "in")
narrowPlot(xgrid = seq(-1, 1, by = 0.5), ygrid = seq(0, 1000, by = 200),
           ylab = "Interval", xlab = "Bounds")
abline(v = 0)
for (ii in 1:nsim) lines(pltMat[,ii], rep(ii, 2),
                         col = adjustcolor("black", 0.2))
abline(h = nsim*(1 - cutoff), lty = 2)
dev.off()

## check level of implicit test
kapInd <- 81
levels <- lapply(fixedInts$poolIntervals,
                 function(el) apply(el, 3,
                                    function(mat) mean(is.na(mat))))
kapLevs <- sapply(levels, function(vec) vec[kapInd])
plot(cutoffs, kapLevs)
for (ii in 1:length(cutoffs)) {
    lines(rep(cutoffs[ii], 2),
          kapLevs[ii] +
          c(-1,1)*1.96/sqrt(nsim)*sqrt(kapLevs[ii]*(1 - kapLevs[ii])))
}
abline(a = 0, b = 1)

## plot a realization
real <- 86 # sample(1:nsim, 1)
png("metaPoolCurvesFixed.png", width = 5, height = 3, res = 480,
    units = "in")
narrowPlot(xgrid = seq(-1, 1, by = 0.5),
           ygrid = seq(0, 1, by = 0.2),
           xlab = expression(x),
           ylab = expression(paste("chi(", bold(p), ";", kappa, ")")),
           addGrid = FALSE)
plotRealization(fixedSim, thetas = thetas, kseq = kseq, cols = cols,
                kaps = c(1, 88, 161), ind = real, refKap = 88)
dev.off()

## reject homogeneity: 19
## separation of curves: 55, 79, 86


## FIXED EFFECTS UNBALANCED ##
fixedGenUB <- function(np, ng) fixedNormalUB(np, ,
                                             sd = std)
set.seed(9381278)
fixedSim <- simEvidentialRegion(fixedGen, nsim = nsim, npop = npop,
                                ngroup = ngroup, kseq = kseq,
                                thetas = thetas)
## compute intervals and coverage probabilities
fixedInts <- getCoverage(fixedSim, thetas = thetas, kseq = kseq,
                         cutoffs = cutoffs)



## RANDOM EFFECTS ##
mnsd <- 1
randGen <- function(np, ng) randomNormal(np, ng, mnsd = mnsd,
                                         sd = sqrt(std^2 - mnsd^2))
set.seed(8251506)
randSim <- simEvidentialRegion(randGen, nsim = nsim, kseq = kseq,
                               thetas = thetas)
## compute intervals and coverage probabilities
randInts <- getCoverage(randSim, thetas = thetas, kseq = kseq,
                        cutoffs = cutoffs)

## plot estimates by kappa
png("metaEstRandom.png", width = 3, height = 3, res = 480,
    units = "in")
narrowPlot(xgrid = seq(-3, 3, by = 1.5),
           ygrid = seq(-1.5, 1.5, by = 0.75),
           xlab = expression(paste(log[10], "(", kappa, ")")),
           xlim = c(-3.5, 3.5), #ylim = c(-0.55, 0.55),
           ylab = "Error")
addQuantPoly(t(randInts$poolTheta), kseq = kseq, cex = 0.6)
lines(log(kseq, 10), colMeans(randInts$poolTheta, na.rm = TRUE),
      lty = 2)
abline(h = quantile(randInts$meanTheta, c(0.025, 0.25, 0.75, 0.975)),
       lty = 3)
dev.off()

## plot all intervals for a particular kappa
kapInd <- 88 # 88  for Fisher's
pltMat <- randInts$poolIntervals[[ctind]][,,kapInd]
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

## which samples give the greatest disagreement between classic and
## evidential intervals?
classLens <- meanInt[,2] - meanInt[,1]
poolLens <- apply(poolIntervals[[ctind]][,,kapInd], 2, diff)
plot(classLens, poolLens,
     pch = as.numeric(poolInclude[[ctind]][,kapInd]),
     col = as.numeric(meanInt[,1] <= mn & meanInt[,2] >= mn) + 1)
abline(a = 0, b = 1)


## RATIOS OF MEAN TO POOL ##
## ratios
lenRat <- poolLens/classLens
## get the max and min examples
lenRatSrt <- order(lenRat, na.last = FALSE)
real <- lenRatSrt[15]
## plot some of the NA cases
narrowPlot(xgrid = seq(-2, 2, by = 1), ygrid = seq(0, 1, by = 0.2),
           xlab = expression(hat(theta)),
           ylab = expression(paste("chi(", bold(p), ";", kappa, ")")),
           addGrid = FALSE)
lines(thetas, poolMat[real, 45,], type = 'l',
      col = adjustcolor(cols[2], 0.8))
lines(thetas, poolMat[real, 10,], type = 'l',
      col = adjustcolor(cols[1], 0.8))
abline(v = 0, col = adjustcolor("black", 0.8), lty = 2, lwd = 1)
abline(v = muMat[real,], col = adjustcolor("gray50", 0.4))
included <- matrix(nrow = ngroup, ncol = length(thetas))
for (ii in 1:ngroup) {
    included[ii,] <- muMat[real, ii] + 1.96*sdMat[real, ii] >= thetas &
        muMat[real, ii] - 1.96*sdMat[real, ii] <= thetas
}
allInclude <- apply(included, 2, sum)
lines(thetas, allInclude/ngroup)
abline(h = 0.05, lty = 3)
legend(x = "topleft", legend = c(round(log(kseq[c(1, 41, 81)], 10),
                                       1), "Min"),
       lty = 1, col = cols, cex = 0.8,
       title = expression(paste(log[10], "(", kappa, ")")))
abline(v = c(meanInt[real, 1], meanInt[real, 2]))


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
