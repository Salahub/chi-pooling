## FUNCTIONS #########################################################
## plot marginal histograms on a scatterplot
addMarHists <- function(x, y, xcuts, ycuts) {
    bds <- par()$usr
    rowDist <- table(cut(x, xcuts))
    colDist <- table(cut(y, ycuts)) # marginal distributions
    vboxBds <- c(bds[2], bds[2] + 0.1*(bds[2] - bds[1]), bds[3:4])
    hboxBds <- c(bds[1:2], bds[4], bds[4] + 0.1*(bds[4] - bds[3]))
    ## density boxes
    rect(vboxBds[1], vboxBds[3], vboxBds[2], vboxBds[4], xpd = NA)
    rect(hboxBds[1], hboxBds[3], hboxBds[2], hboxBds[4], xpd = NA)
    ## add marginal histograms
    vseq <- xcuts
    rect(vboxBds[1], vseq[1:length(colDist)],
         vboxBds[1] + 0.9*diff(vboxBds[1:2])*(colDist/max(colDist)),
         vseq[2:(length(colDist) + 1)], xpd = NA,
         col = adjustcolor("firebrick", 0.5))
    hseq <- ycuts
    rect(hseq[1:length(rowDist)], hboxBds[3],
         hseq[2:(length(rowDist) + 1)], xpd = NA,
         hboxBds[3] + 0.9*diff(hboxBds[3:4])*(rowDist/max(rowDist)),
         col = adjustcolor("firebrick", 0.5))
}

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
        if (labpos[1] == "left") { # change label positions...
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
        } # ...these can be modified to make them legible
        text(x = log(kseq[xind], 10), # labels
             y = quantile(mat[xind,], yq),
             labels = 1-2*qn, adj = adj, ...)
    }
    lines(x = log(kseq, 10), # median line
          y = apply(mat, 1, median))
}

## parameter sweep function, generates p-values given estimates,
## their standard deviations, candidate combined estimates, and a
## p-value function
paramSweep <- function(mus, sds, thetas, pfun = pnorm) {
    mus2 <- rep(mus, each = length(thetas))
    sds2 <- rep(sds, each = length(thetas)) # vectorize operation
    stdDs <- abs((mus2 - rep(thetas, times = length(mus)))/sds2)
    matrix(2*pfun(q = stdDs, lower.tail = FALSE),
           nrow = length(thetas)) # convert to matrix
}

## assuming the output of a paramSweep call, this function sweeps
## through each collection of p-values for proposed combined
## estimates and then computes the pooled p-value over a range of
## kappa values in the chi-squared pooled p-value
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

## the following functions generate particular meta-analysis settings
## 1. the case of inhomogeneous means (i.e. no shared underlying mean)
inhomNormal <- function(npop, ngroup, mns = runif(-2, 2, ngroup),
                        sd = 1) {
    obs <- rep(mns, npop/ngroup) + rnorm(npop, mean = 0, sd = sd)
    grp <- rep(1:ngroup, npop/ngroup)
    cbind(x = obs, group = grp)
}
## 2. a fixed effects model generation function with equal sample
## sizes
fixedNormal <- function(npop, ngroup, mn = 0, sd = 1) {
    obs <- rnorm(npop, mean = mn, sd = sd)
    grp <- sample(rep(1:ngroup, npop/ngroup))
    cbind(x = obs, group = grp)
}
## 3. fixed effects with different sample sizes
fixedNormalUB <- function(npop, groups = sample(rep(1:8, npop/8)),
                          mn = 0, sd = 1) {
    obs <- rnorm(npop, mean = mn, sd = sd)
    grp <- groups
    cbind(x = obs, group = grp)
}
## 4. random effects model generation function
randomNormal <- function(npop, ngroup, mn = 0, mnsd = 1,
                         sd = 1) {
    mns <- rnorm(ngroup, mean = mn, sd = mnsd)
    obs <- rep(mns, npop/ngroup) + rnorm(npop, mean = 0, sd = sd)
    grp <- rep(1:ngroup, npop/ngroup)
    cbind(x = obs, group = grp)
}

## this wrapper generates nsim simulated meta-analysis results using
## repeated calls to genFn based on the settings provided in further
## arguments
simMetaStudies <- function(genFn, nsim = 1e3, npop = 240,
                           ngroup = 8) {
    muMat <- sdMat <- matrix(nrow = nsim, ncol = ngroup)
    nMat <- matrix(nrow = nsim, ncol = ngroup)
    for (ii in 1:nsim) { # repeat generation nsim times
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

## this function takes a probability generating function which
## gives p-values (pfun) for estimates for each in a range of
## candidate combined estimate values (thetas) and applies this to
## every collection of estimates in studies
metaToP <- function(pfun, studies, thetas) {
    simplify2array(lapply(thetas,
                          function(x) {
                              pfun(x, studies$means, studies$sds,
                                   studies$ns)
                          }))
}

## the following functions accept x (a candidate combined estimate),
## mus (observed estimates), sds (corresponding standard deviations),
## and ns (corresponding sample sizes) and return p-values for each
## estimate presuming x is the true population value
## 1. a function based on normal p-values
pfunNorm <- function(x, mus, sds, ns) {
    2*pnorm(-abs(mus - x)/sds)
}
## 2. a function based on t_{n-1} distribution p-values
pfunT <- function(x, mus, sds, ns) {
    2*pt(-abs(mus - x)/sds, df = ns - 1)
}

## given p-values resulting from a sweep of some pfun over a range of
## candidate thetas over a bunch of experimental results, this
## converts every group of p-values to chi pool estimates for every
## kappa in kseq
chiMetaSweep <- function(ps, kseq = exp(seq(-8, 8, by = 0.1))) {
    M <- dim(ps)[2] # number of p-values
    arr <- simplify2array(lapply(
        kseq, # for every kappa
        function(k) {
            pchisq(apply(qchisq(ps, k, lower.tail = FALSE), # pool ps
                         c(1,3), sum),
                   df = M*k, lower.tail = FALSE)
        }))
    aperm(arr, c(1,3,2)) # rear
}

## weighted chi computation: this is an extension not covered in the
## thesis document, but applies Lancaster's weighting to the chi
## pooling function where each p-value is converted by a different
## chi-squared quantile with degrees of freedom inversely proportional
## to the variance
chiWeighted <- function(ps, studies) {
    vs <- studies$sds # get standard deviations
    K <- nrow(vs)
    wgts <- 1/vs^2 # inverse variance
    t(sapply(1:nrow(wgts), # applies the chi pooled p-value with
           function(ii) {  # different dfs for each p to the data
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

## safe estimates for a range of cutoffs, produces NAs rather than
## empty elements to make summaries sensible
safeRange <- function(x) {
    if (length(x) == 0) {
        c(NA, NA)
    } else range(x)
}

## compute coverage probabilities for a many sequences of chi pooled
## p-values generated by a sweep of kappa values on some samples
getChiEsts <- function(chiseqs, thetas, kseq, mn = 0,
                       cutoffs = c(0.1, 0.05, 0.02, 0.01)) {
    ests <- matrix(thetas[apply(chiseqs, c(1,2), # get EME (the max)
                                     which.max)],
                     ncol = length(kseq))
    intervals <- lapply(cutoffs, # get evidential regions (> cutoff)
                        function(ct) {
                            apply(chiseqs >= ct,
                                  c(1,2),
                                  function(x) {
                                      safeRange(thetas[which(x)])
                                  })
                        })
    include <- lapply(intervals, # does this include the true value?
                      function(mat) mat[1,,] <= mn & mat[2,,] >= mn)
    covP <- lapply(include, # compute coverage probabilities
                   function(mat) apply(mat, 2, mean, na.rm = TRUE))
    list(ests = ests, intervals = intervals, covP = covP)
}

## a version of the above which is accepts only a single sequence
## rather than many sequences
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

## the mean estimate function, performs classical combination of
## estimates by weighting estimates by the inverse of their variances
getMeanEsts <- function(sim, cutoffs = c(0.1, 0.05, 0.02, 0.01)) {
    mns <- sim$means # estimates
    wgts <- 1/sim$sds^2 # weights
    wgted <- mns*(wgts)
    ests <- rowSums(wgted)/rowSums(wgts)
    se <- sqrt(1/rowSums(wgts)) # resulting standard error
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

## plot a realization based on simulated data
plotRealization <- function(chi, means, sds, thetas, kseq,
                            cols = 1:4, kaps = c(1, 41, 81),
                            refKap = 41, legend = TRUE) {
    ## lines for chosen kappa indices
    for (ii in kaps) lines(thetas, chi[ii,], type = 'l',
                           col = if (ii <= refKap-1) {
                                     adjustcolor(cols[1], 0.8)
                                 } else if (ii == refKap) {
                                     adjustcolor(cols[2], 0.8)
                                 } else adjustcolor(cols[3], 0.8))
    ## classical combined estimate
    meanEst <- sum(means/sds^2)/sum(1/sds^2)
    lines(x = c(meanEst - qnorm(0.975)/sqrt(sum(1/sds^2)),
                meanEst + qnorm(0.975)/sqrt(sum(1/sds^2))),
          y = rep(-0.1, 2))
    points(meanEst, y = -0.1, cex = 0.8)
    ## plot each estimate with a reference line
    for (ii in seq_along(means)) {
        lines(x = c(means[ii] - 2*sds[ii], means[ii] + 2*sds[ii]),
              y = rep(-0.025 - 0.05*(ii-1)/(length(means)-1), 2),
              lwd = 1, col = adjustcolor("black", 1/2))
    }
    abline(v = meanEst, col = adjustcolor("black", 0.6), lwd = 1)
    abline(v = means, col = adjustcolor("gray50", 0.4))
    abline(h = 0.05, lty = 3)
    if (legend) { # add a legend (maybe)
        temp <- legend("topleft", legend = rep("", length(kaps)),
                       text.width = strwidth("-3.5"), col = cols,
                       cex = 0.8, lty = 1, xjust = 1, yjust = 1,
                       title = expression(paste(log[10], "(",
                                                kappa, ")")))
        text(temp$rect$left + temp$rect$w, temp$text$y,
             round(log(kseq[kaps], 10), 1), pos = 2, cex = 0.8)
    }
}

## CANONICAL EXAMPLES ################################################
## equal variance: different arrangements of estimates where all have
## the same variance
## parameters
mn <- 0
mnstd <- 1 # sqrt(0.5)
kseq <- exp(seq(-8, 8, by = 0.1)) # kappa values
xseq <- seq(-4, 4, by = 0.01) # estimate range
cols <- RColorBrewer::brewer.pal(3, "Dark2") # palette
## symmetric settings: the mean of estimates is zero and the pattern
## is symmetric about zero
syms <- list(c(-3.001, -1.501, -0.501, -0.051, 0.051, 0.501,
               1.501, 3.001),
             c(-2.001, -1.501, -0.501, -0.051, 0.051, 0.501,
               1.501, 2.001),
             seq(-3.001, 3.001, length.out = 8),
             seq(-2.001, 2.001, length.out = 8),
             seq(-0.501, 0.501, length.out = 8))
## p values
syms.p <- simplify2array(lapply(syms, function(mns) {
    2*pnorm(-abs(outer(xseq, mns, `-`)), sd = mnstd)
}))
## sweep kappa values
syms.pool <- chiMetaSweep(aperm(syms.p, c(3, 2, 1)), kseq = kseq)

## plot realizations
ind <- 2 # choose a value in 1-5
wid <- if (ind == 2) 2.3 else 2.1 # graphical parameters
hei <- 2.3
mars <- if (ind == 2) {
            c(1.1, 2.1, 0.1, 0.1)
        } else c(1.1, 1.1, 0.1, 0.1)
suff <- if (mnstd == 1) "" else "sd0_5" # naming of files
leg <- ind == 2 & mnstd == 1 # add a legend?
png(paste0("syms", ind, suff, ".png"), width = wid, height = hei,
    res = 480, units = "in")
narrowPlot(xgrid = seq(-3, 3, by = 1.5), xlab = "x",
           ylab = expression(paste("chi(", bold(p), "(x)",
                                   ";", kappa, ")")),
           ygrid = seq(0, 1, by = 0.25),
           addGrid = FALSE, ylim = c(-0.1, 1),
           mars = mars)
abline(h = 0)
plotRealization(syms.pool[ind,,], means = syms[[ind]],
                sds = rep(mnstd, 8), thetas = xseq,
                kseq = kseq, kaps = c(1, 88, 161),
                refKap = 88, cols = cols,
                legend = leg)
dev.off()

## unbalanced settings: the mean is no longer constant across settings
## and the patterns of points are no longer symmetric
unbs <- list(-c(-2.501, -1.501, -0.751, -0.301, -0.101, -0.051,
               0.051, 2.501),
             -c(-2.001,-1.901, -1.751, -1.501, -1.351, -1.121,
               -1.001, 2.001),
             -c(-2.001, -1.901, -1.751, -1.501, -1.351, -1.121,
               -1.001, 1.501),
             -c(-2.001, -1.901, -1.751, -1.501, -1.351, -1.121,
               -1.001, 3.001),
             -c(-2.001, -1.901, -1.751, -1.501, -1.351, -1.121,
               -1.001, 4.001))
unbs.p <- simplify2array(lapply(unbs, function(mns) {
    2*pnorm(-abs(outer(xseq, mns, `-`)), sd = mnstd)
}))
## sweep kappa values
unbs.pool <- chiMetaSweep(aperm(unbs.p, c(3, 2, 1)), kseq = kseq)

## plot realizations: same as before
ind <- 2 # choose a value in 1-5
wid <- if (ind == 2) 2.3 else 2.1
hei <- 2.3
mars <- if (ind == 2) {
            c(1.1, 2.1, 0.1, 0.1)
        } else c(1.1, 1.1, 0.1, 0.1)
suff <- if (mnstd == 1) "" else "sd0_5"
leg <- ind == 2 & mnstd == 1
png(paste0("unbs", ind, suff, ".png"), width = wid, height = hei,
    res = 480, units = "in")
narrowPlot(xgrid = seq(-3, 3, by = 1.5), xlab = "x",
           ylab = expression(paste("chi(", bold(p), "(x)",
                                   ";", kappa, ")")),
           ygrid = seq(0, 1, by = 0.25),
           addGrid = FALSE, ylim = c(-0.1, 1),
           mars = mars)
abline(h = 0)
plotRealization(unbs.pool[ind,,], means = unbs[[ind]],
                sds = rep(mnstd, 8), thetas = xseq,
                kseq = kseq, kaps = c(1, 88, 161),
                refKap = 88, cols = cols,
                legend = leg)
dev.off()

## unequal variance: the same symmetric and unbalanced settings are
## repeated, but this time the estimates have different variances
mnstds <- c(2, rep(1, 6), 1/2)
## symmetric settings
uvsyms <- list(c(-3.001, -1.501, -0.501, -0.051, 0.051, 0.501,
               1.501, 3.001),
             c(-2.001, -1.501, -0.501, -0.051, 0.051, 0.501,
               1.501, 2.001),
             seq(-3.001, 3.001, length.out = 8),
             seq(-2.001, 2.001, length.out = 8),
             seq(-0.501, 0.501, length.out = 8))
## p values
uvsyms.p <- simplify2array(lapply(uvsyms, function(mns) {
    2*pnorm(-abs(sweep(outer(xseq, mns, `-`), 2, mnstds, `/`)))
}))
## sweep kappa values
uvsyms.pool <- chiMetaSweep(aperm(uvsyms.p, c(3, 2, 1)), kseq = kseq)

## plot realizations: same as last two times
ind <- 2 # choose a value in 1-5
wid <- if (ind == 2) 2.3 else 2.1
hei <- 2.5
mars<- if (ind == 2) {
           c(2.1, 2.1, 0.1, 0.1)
       } else c(2.1, 1.1, 0.1, 0.1)
png(paste0("unevsym", ind, ".png"), width = wid, height = hei,
    res = 480, units = "in")
narrowPlot(xgrid = seq(-3, 3, by = 1.5), xlab = "x",
           ylab = expression(paste("chi(", bold(p), "(x)",
                                   ";", kappa, ")")),
           ygrid = seq(0, 1, by = 0.25),
           addGrid = FALSE, ylim = c(-0.1, 1),
           mars = mars)
abline(h = 0)
plotRealization(uvsyms.pool[ind,,], means = uvsyms[[ind]],
                sds = mnstds, thetas = xseq,
                kseq = kseq, kaps = c(1, 88, 161),
                refKap = 88, cols = cols,
                legend = FALSE)
dev.off()
## compared to the equal variance case, the evidential regions and
## EMEs are pulled towards the less variable estimates

## unbalanced settings
uvunbs <- list(-c(-2.501, -1.501, -0.751, -0.301, -0.101, -0.051,
                  0.051, 2.501),
               -c(-2.001, -1.901, -1.751, -1.501, -1.351, -1.121,
                  -1.001, 2.001),
               -c(-2.001, -1.901, -1.751, -1.501, -1.351, -1.121,
                  -1.001, 1.501),
               -c(-2.001, -1.901, -1.751, -1.501, -1.351, -1.121,
                  -1.001, 3.001),
               -c(-2.001, -1.901, -1.751, -1.501, -1.351, -1.121,
                  -1.001, 4.001))
uvunbs.p <- simplify2array(lapply(uvunbs, function(mns) {
    2*pnorm(-abs(sweep(outer(xseq, mns, `-`), 2, mnstds, `/`)))
}))
## sweep kappa values
uvunbs.pool <- chiMetaSweep(aperm(uvunbs.p, c(3, 2, 1)), kseq = kseq)

## plot realizations
ind <- 5 # choose a value in 1-5
wid <- if (ind == 2) 2.3 else 2.1
hei <- 2.5
mars <- if (ind == 2) {
            c(2.1, 2.1, 0.1, 0.1)
        } else c(2.1, 1.1, 0.1, 0.1)
png(paste0("unevunb", ind, ".png"), width = wid, height = hei,
    res = 480, units = "in")
narrowPlot(xgrid = seq(-3, 3, by = 1.5), xlab = "x",
           ylab = expression(paste("chi(", bold(p), "(x)",
                                   ";", kappa, ")")),
           ygrid = seq(0, 1, by = 0.25),
           addGrid = FALSE, ylim = c(-0.1, 1),
           mars = mars)
abline(h = 0)
plotRealization(uvunbs.pool[ind,,], means = uvunbs[[ind]],
                sds = mnstds, thetas = xseq,
                kseq = kseq, kaps = c(1, 88, 161),
                refKap = 88, cols = cols,
                legend = leg)
dev.off()

## final case: adding internal points between a pair of distant
## estimates, what happens?
addm <- list(c(-3.001, 3.001),
             c(-3.001, rep(1, 5), 3.001),
             c(-3.001, rep(1, 10), 3.001),
             c(-3.001, rep(1, 15), 3.001),
             c(-3.001, rep(1, 25), 3.001),
             c(-3.001, rep(1, 1000), 3.001))
addm.p <- simplify2array(lapply(addm, function(mns) {
    2*pnorm(-abs(sweep(outer(xseq, mns, `-`), 2, mnstd, `/`)))
}))
## sweep kappa values
addm.pool <- lapply(addm.p, function(ps) {
    simplify2array(lapply(
        kseq,
        function(k) {
            M <- ncol(ps)
            pchisq(apply(qchisq(ps, k, lower.tail = FALSE),
                         1, sum),
                   df = M*k, lower.tail = FALSE)
        }))})

## plot realizations
ind <- 1 # choose a value in 1-6
wid <- if (ind == 1) 2.3 else 2.1
mars <- if (ind == 1) {
            c(2.1, 2.1, 1.1, 0.1)
        } else c(2.1, 1.1, 1.1, 0.1)
leg <- ind == 1
png(paste0("adm", ind, ".png"), width = wid, height = 2.7,
    res = 480, units = "in")
narrowPlot(xgrid = seq(-3, 3, by = 1.5), xlab = "x",
           ylab = expression(paste("chi(", bold(p), "(x)",
                                   ";", kappa, ")")),
           ygrid = seq(0, 1, by = 0.25),
           addGrid = FALSE, ylim = c(-0.1, 1),
           mars = mars)
abline(h = 0)
plotRealization(t(addm.pool[[ind]]), means = addm[[ind]],
                sds = rep(mnstd, length(addm[[ind]])),
                thetas = xseq, kseq = kseq, kaps = c(1, 88, 161),
                refKap = 88, cols = cols,
                legend = leg)
dev.off()


## SIMULATIONS #######################################################
## the following perform more in depth simulations of different meta
## analysis settings in order to evaluate the power and performance
## of evidential methods
## settings
mn <- 0 # true mean
std <- 2 # total standard deviation
nsim <- 1000 # number of simulations
ngroup <- 8 # number of groups
npop <- 240 # total population size
minQuants <- readRDS("./results/curveMinQuantiles.Rds") # null dist
kseq <- exp(seq(-8, 8, by = 0.1)) # kappa sequence
thetas <- seq(-1.5, 1.5, by = 0.01) # proposed combined estimates
cols <- RColorBrewer::brewer.pal(4, "Dark2") # palette
cutoffs <- seq(0.20, 0.01, by = -0.01) # evidential threshold values

## FIXED EFFECTS ##
## consider first the case of fixed effects with equal variances
## for every estimate
## compute all the features
fixedGen <- function(np, ng) fixedNormal(np, ng,
                                         sd = std)
set.seed(9381278) # reproducibility
## simulate the data
fixedSim <- simMetaStudies(fixedGen, nsim = nsim, npop = npop,
                           ngroup = ngroup)
## compute p-values
fixedps <- metaToP(pfunT, fixedSim, thetas = thetas)
## compute evidential estimates
fixedChis <- chiMetaSweep(fixedps, kseq = kseq) ## ~ 20 mins
fixedMnChi <- apply(fixedChis, c(1,3), min) # classical estimate
fixedwgtChi <- chiWeighted(fixedps, fixedSim) # lancaster wgts
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

## plot estimates by kappa agaisnt the classic estimate
kapInd <- 81
png("metaEstvsMeanFixed.png", width = 3, height = 3, res = 480,
    units = "in")
narrowPlot(xgrid = seq(-0.5, 0.5, by = 0.25),
           ygrid = seq(-0.5, 0.5, by = 0.25),
           xlab = expression(hat(theta)^{"("~E~")"}),
           ylab = expression(hat(theta)))
points(fixedChiInt$est[, kapInd], fixedMnInt$ests, cex = 0.5, pch = 19,
       col = adjustcolor("black", 0.3))
dev.off()

## plot the coverage probabilities by kappa
ctind <- 20
cutoff <- cutoffs[ctind]
png(paste0("meta", 100*cutoff, "pctCovPFix.png"), width = 2.5,
    height = 2.5, res = 480, units = "in")
## optionally: fit a smooth to these
tempCovP <- lowess(fixedChiInt$covP[[ctind]], f = 1/6)$y
tempCovP <- fixedChiInt$covP[[ctind]]
narrowPlot(xgrid = seq(-3, 3, by = 3), ygrid = seq(0.8, 1, by = 0.05),
           ylab = expression(pi), xlim = c(-3.5, 3.5),
           xlab = expression(paste(log[10], "(", kappa, ")")))
lines(log(kseq, 10), tempCovP)
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
pltMat <- fixedChiInt$intervals[[ctind]][,,kapInd] # chi case
##pltMat <- t(fixedMnInt[order(fixedMnInt[,1]),]) # classical case
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
kapInd <- 88
levels <- lapply(fixedChiInt$intervals,
                 function(el) apply(el, 3,
                                    function(mat) mean(is.na(mat))))
kapLevs <- sapply(levels, function(vec) vec[kapInd])
png("metaRejectPFixed.png", width = 2.5, height = 2.5, res = 480,
    units = "in")
narrowPlot(xgrid = seq(0, 0.2, by = 0.05),
           ygrid = seq(0, 0.2, by = 0.05),
           xlab = "a", ylab = expression(alpha))
points(cutoffs, kapLevs, cex = 0.8)
abline(a = 0, b = 1)
## fit and add a linear model
alphaA <- lm(alpha ~ a + I(a^2) - 1,
             data = data.frame(alpha = kapLevs, a = cutoffs))
lines(seq(0, 0.2, by = 0.01),
      predict(alphaA,
              newdata = data.frame(a = seq(0, 0.2, by = 0.01))),
      col = "firebrick")
dev.off()

## check coverage probabilities for a particular kappa
kapInd <- 88
png("metaCovPFixed.png", width = 2.5, height = 2.5, res = 480,
    units = "in")
covPs <- sapply(fixedChiInt$covP, function(el) el[kapInd])
narrowPlot(xgrid = seq(0, 0.2, by = 0.05),
           ygrid = seq(0, 0.2, by = 0.05),
           xlab = "a", ylab = expression({1-pi}))
points(cutoffs, 1 - covPs, cex = 0.8)
## fit and add a linear model
piA <- lm(pi ~ a + I(a^2) - 1,
          data = data.frame(pi = 1 - covPs, a = cutoffs))
lines(seq(0, 0.2, by = 0.01),
      predict(piA, newdata = data.frame(a = seq(0, 0.2, by = 0.01))),
      col = "firebrick")
abline(a = 0, b = 1)
dev.off()

## compare widths of kappa = 2 to mean CI
ctind <- 16
cutoff <- cutoffs[ctind]
kapInd <- 88
kapWid <- abs(apply(fixedChiInt$intervals[[ctind]][,,kapInd], 2, diff))
CIwid <- abs(apply(fixedMnInt$intervals[[ctind]], 2, diff))
png("metaCIWidToEvWidDetl.png", width = 3, height = 3, units = "in",
    res = 480)
narrowPlot(xgrid = seq(0.35, 0.5, by = 0.05),
           ygrid = seq(0, 1, by = 0.25),
           xlab = "Width of confidence interval",
           ylab = "Width of evidential interval",
           mars = c(2.1, 2.1, 1.2, 1.2))
points(CIwid, kapWid, pch = 20, cex = 0.8,
       col = adjustcolor("black", 0.2))
abline(0, 1)
points(CIwid[c(671, 208, 61)],
       kapWid[c(671, 208, 61)],
       pch = c(15, 17, 19), col = "firebrick", cex = 0.8)
addMarHists(CIwid, kapWid, xcuts = seq(0, 1, by = 0.05),
            ycuts = seq(0, 1, by = 0.05))
bds <- par()$usr
abline(h = bds[4], xpd = NA, col = "white")
lines(y = c(bds[4], bds[4]), x = c(bds[1], bds[2] + 0.1*0.162),
      xpd = NA)
dev.off()

## plot a realization
real <- 61 # sample(1:nsim, 1)
wid <- if (real == 61) 2.3 else 2.1
mars <- if (real == 61) {
            c(2.1, 2.1, 1.1, 0.1)
        } else c(2.1, 1.1, 1.1, 0.1)
leg <- real == 61
png(paste0("metaPoolCurveFixed", real, ".png"),
    width = wid, height = 2.7, res = 480, units = "in")
narrowPlot(xgrid = seq(-1, 1, by = 0.5), xlab = "x",
           ylab = expression(paste("chi(", bold(p), "(x)",
                                   ";", kappa, ")")),
           ygrid = seq(0, 1, by = 0.25),
           addGrid = FALSE, ylim = c(-0.1, 1),
           mars = mars)
abline(h = 0)
plotRealization(fixedChis[real,,], means = fixedSim$means[real,],
                sds = fixedSim$sds[real,],
                thetas = thetas, kseq = kseq, kaps = c(1, 88, 161),
                refKap = 88, cols = cols,
                legend = leg)
dev.off()


## FIXED EFFECTS UNBALANCED ##
## next: the case of fixed effects where not every estimate has
## equal variance
fixedGenUB <- function(np, ng) {
    grps <- sample(1:ng, size = np, replace = TRUE,
                   prob = c(3, 3, 6, 2, 2, 2, 3, 4))
    fixedNormalUB(np, groups = grps, sd = std)
}
## generate and compute on studies
set.seed(73183)
## generate data
ubSim <- simMetaStudies(fixedGenUB, nsim = nsim, npop = npop,
                        ngroup = ngroup)
## compute p-values
ubps <- metaToP(pfunT, ubSim, thetas = thetas)
## compute chi pooled p-values
ubChis <- chiMetaSweep(ubps, kseq = kseq) ## ~ 20 mins
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
dev.off()

## plot evidential estimates against classical estimates
kapInd <- 81 # 88 is fishers, 81 matches closest
png("metaEstvsMeanUB.png", width = 3, height = 3, res = 480,
    units = "in")
narrowPlot(xgrid = seq(-0.5, 0.5, by = 0.25),
           ygrid = seq(-0.5, 0.5, by = 0.25),
           xlab = expression(hat(theta)^{"("~E~")"}),
           ylab = expression(hat(theta)))
points(ubChiInt$est[, kapInd], ubMnInt$ests, cex = 0.5, pch = 19,
       col = adjustcolor("black", 0.3))
dev.off()

## plot the coverage probabilities for a given kappa
ctind <- 16
cutoff <- cutoffs[ctind]
png(paste0("meta", 100*cutoff, "pctCovPUB.png"), width = 2.5,
    height = 2.5, res = 480, units = "in")
tempCovP <- ubChiInt$covP[[ctind]]
narrowPlot(xgrid = seq(-3, 3, by = 3), ygrid = seq(0.8, 1, by = 0.05),
           ylab = expression(pi), xlim = c(-3.5, 3.5),
           xlab = expression(paste(log[10], "(", kappa, ")")))
lines(log(kseq, 10), tempCovP)
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
pltMat <- ubChiInt$intervals[[ctind]][,,kapInd] # chi method
##pltMat <- ubWgtInt$intervals[[ctind]] # another method
##pltMat <- t(fixedInts$meanInt[order(fixedInts$meanInt[,1]),]) #classic
pltMat <- pltMat[, order(pltMat[1, ])]
png(paste0("poolInts", 100*cutoff, "pctUB.png"), height = 2.5,
    width = 2.5, res = 480, units = "in")
narrowPlot(xgrid = seq(-1, 1, by = 0.5), ygrid = seq(0, 1000, by = 200),
           ylab = "Interval", xlab = "Bounds")
abline(v = 0)
for (ii in 1:nsim) lines(pltMat[,ii], rep(ii, 2),
                         col = adjustcolor("black", 0.2))
abline(h = nsim*(1 - cutoff), lty = 2)
dev.off()

## check level of implicit test
kapInd <- 88
levels <- lapply(ubChiInt$intervals,
                 function(el) apply(el, 3,
                                    function(mat) mean(is.na(mat))))
kapLevs <- sapply(levels, function(vec) vec[kapInd])
png("metaRejectPUB.png", width = 2.5, height = 2.5, res = 480,
    units = "in")
narrowPlot(xgrid = seq(0, 0.2, by = 0.05),
           ygrid = seq(0, 0.2, by = 0.05),
           xlab = "a", ylab = expression(alpha))
points(cutoffs, kapLevs, cex = 0.8)
abline(a = 0, b = 1)
## fit a linear model and add this to the plot
alphaA <- lm(alpha ~ a + I(a^2) - 1,
             data = data.frame(alpha = kapLevs, a = cutoffs))
lines(seq(0, 0.2, by = 0.01),
      predict(alphaA,
              newdata = data.frame(a = seq(0, 0.2, by = 0.01))),
      col = "firebrick")
dev.off()

## check coverage probabilities for a particular kappa
kapInd <- 88
png("metaCovPUB.png", width = 2.5, height = 2.5, res = 480,
    units = "in")
covPs <- sapply(ubChiInt$covP, function(el) el[kapInd])
narrowPlot(xgrid = seq(0, 0.2, by = 0.05),
           ygrid = seq(0, 0.2, by = 0.05),
           xlab = "a", ylab = expression({1-pi}))
points(cutoffs, 1 - covPs, cex = 0.8)
## fit a linear model and add this to the plot
piA <- lm(pi ~ a - 1,
          data = data.frame(pi = 1-covPs, a = cutoffs))
lines(seq(0, 0.2, by = 0.01),
      predict(piA,
              newdata = data.frame(a = seq(0, 0.2, by = 0.01))),
      col = "firebrick")
abline(a = 0, b = 1)
dev.off()

## compare widths of kappa = 2 to mean CI
ctind <- 16
cutoff <- cutoffs[ctind]
kapInd <- 88
kapWid <- abs(apply(ubChiInt$intervals[[ctind]][,,kapInd], 2, diff))
CIwid <- abs(apply(ubMnInt$intervals[[ctind]], 2, diff))
png("metaCIWidToEvWidDetlub.png", width = 3, height = 3, units = "in",
    res = 480)
narrowPlot(xgrid = seq(0.3, 0.5, by = 0.05),
           ygrid = seq(0, 1.2, by = 0.3),
           xlab = "Width of confidence interval",
           ylab = "Width of evidential interval",
           mars = c(2.1, 2.1, 1.2, 1.2))
points(CIwid, kapWid, pch = 20, cex = 0.8,
       col = adjustcolor("black", 0.2))
abline(0, 1)
addMarHists(CIwid, kapWid, xcuts = seq(0, 1.2, by = 0.05),
            ycuts = seq(0, 1.2, by = 0.05))
bds <- par()$usr
abline(h = bds[4], xpd = NA, col = "white")
lines(y = c(bds[4], bds[4]), x = c(bds[1], bds[2] + 0.1*0.4),
      xpd = NA)
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
plotRealization(ubChis[real,,], means = ubSim$means[real,],
                sds = ubSim$sds[real,],
                thetas = thetas, kseq = kseq, kaps = c(1, 88, 161),
                refKap = 88, cols = cols,
                legend = leg)
dev.off()


## RANDOM EFFECTS ##
## this final case should commonly lead to a rejection of
## homogeneity, as it explicitly simulates the estimates to have
## different underlying means
mnsd <- sqrt(0.5)
## generation function
randGen <- function(np, ng) randomNormal(np, ng, mnsd = mnsd,
                                         sd = sqrt(std^2 - mnsd^2))
set.seed(53611707)
## generate data
randSim <- simMetaStudies(randGen, nsim = nsim, npop = npop,
                          ngroup = ngroup)
## p-values
randps <- metaToP(pfunT, randSim, thetas = thetas)
## pooled p-values
randChis <- chiMetaSweep(randps, kseq = kseq) ## ~ 20 mins
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
pltMat <- randChiInt$intervals[[ctind]][,,kapInd] # evidential
## classical
##pltMat <- t(fixedInts$meanInt[order(fixedInts$meanInt[,1]),])
pltMat <- pltMat[, order(pltMat[1, ])]
png(paste0("poolInts", 100*cutoff, "pctRand.png"), height = 2.5,
    width = 2.5, res = 480, units = "in")
narrowPlot(xgrid = seq(-1.5, 1.5, by = 0.75),
           ygrid = seq(0, 1000, by = 200),
           ylab = "Interval", xlab = "Bounds")
abline(v = 0)
for (ii in 1:nsim) lines(pltMat[,ii], rep(ii, 2),
                         col = adjustcolor("black", 0.2))
abline(h = nsim*(1 - cutoff), lty = 2)
dev.off()
## indeed, most intervals are empty!

## check for a range of sds
mnVSeq <- seq(0.05, 3.95, by = 0.15)
## a list of closures with the correct variances
randGens <- lapply(mnVSeq,
                   function(mnV) {
                       function(np, ng) {
                           randomNormal(np, ng,
                                        mnsd = mnV,
                                        sd = sqrt(std^2 - mnV))
                       }})
## simulate for each (lower resolution to speed it up)
kseq2 <- exp(seq(-8, 8, by = 0.5))
kseq2 <- c(kseq2[1:17], 2, kseq2[18:33])
set.seed(54910913)
## simulate data
randSims <- lapply(randGens, simMetaStudies, nsim = nsim, npop = npop,
                   ngroup = ngroup)
## second of these still takes ~ 20 mins
randps <- lapply(randSims, metaToP, pfun = pfunT, thetas = thetas)
randChis <- lapply(randps, chiMetaSweep, kseq = kseq2)
## get intervals
randChiInt <- lapply(randChis, getChiEsts, thetas = thetas,
                     kseq = kseq2, cutoffs = cutoffs)
## and rejection probabilities
randRejectP <- lapply(randChiInt,
                      function(ints) {
                          sapply(ints$intervals,
                                 function(mat) {
                                     colMeans(is.na(mat[1,,]))
                                 })})
randRejectP <- simplify2array(randRejectP)
## do the same for the MLE
randMn <- lapply(randSims, getMeanEsts, cutoffs = cutoffs)
randTests <- mapply(function(sim, mns) {
    rowSums((sweep(sim$means, 1, mns$ests, `-`)/sim$sds)^2)},
    randSims, randMn)
randMnPow <- colMeans(pchisq(randTests, df = 7, lower.tail = FALSE) <=
                      #0.05) # naive threshold
                      0.514*0.05 + 0.891*0.05^2) # corrected threshold
## corrected threshold is based on the linear model, gives comparable
## rejection probability as the classic test

## set these powers up as data frame for further analysis
powerdf <- data.frame(pow = c(randRejectP),
                      expand.grid(lgk = log(kseq2, 10), a = cutoffs,
                                  taup = mnVSeq/4))

## the next section plots some power contours
## some plotting packages
library(ggplot2)
library(grid)

## define a helper that makes nice labels
ggLabs <- function(breaks){
    paste0(c("[", rep("(", length(breaks)-2)), # bottom brackets
           breaks[-length(breaks)], # break limits
           ", ", breaks[-1], "]") # top brackets
}

## set up some graphical parameters
nbr <- 12 # number of breaks for simple power
powBreaks <- c(seq(0, 1 - 1/nbr, length.out = nbr), 1.01)
powLabs <- ggLabs(round(powBreaks, 1)) # nice labels
powPal <- colorRampPalette(c("floralwhite", "firebrick")) # palette

## set a new theme based on the "lined raw" template
myTheme <- theme_linedraw()
myTheme$text$size <- 8 # text size
myTheme$plot.title$size <- rel(1.1) # scale down title
myTheme$strip.text$colour <- "black" # facet title colour
myTheme$strip.text$size <- 8 # facet title size
myTheme$strip.background$fill <- "white"
myTheme$panel.grid$colour <- adjustcolor("gray", 0.4)
myTheme$plot.title$hjust <- 0.5 # centre title

## plot the contours
inds <- abs(powerdf$a - 0.05) < 0.001 & powerdf$taup < 0.5
powerContourBase <- ggplot(data = powerdf[inds,],
                           aes(taup, lgk, z = pow)) # base ggplot
powGrob <- ggplotGrob(powerContourBase + # save as a grid grob
                      xlab(expression({tau^2/4})) +
                      ylab(expression(log[10]~{(kappa)})) +
                      geom_contour_filled(breaks = powBreaks) +
                      scale_fill_manual(values = powPal(nbr+1),
                                        name = "Power",
                                        labels = powLabs,
                                        drop = FALSE,
                                        guide = "none") +
                      geom_hline(yintercept = 0.3) +
                      myTheme) # set theme
png("poolPowersRand.png", width = 3, height = 3, res = 480,
    units = "in")
grid.newpage()
pushViewport(viewport(x = unit(0.45, "npc"), # facet plot viewport
                      width = unit(0.9, "npc")))
grid.draw(powGrob) # facet plot
popViewport()
pushViewport(viewport(x = unit(0.93, "npc"), # legend viewport
                      height = unit(0.4, "npc"),
                      width = unit(0.05, "npc")))
grid.text("Power", y = unit(1, "npc") + unit(1, "lines"), # title
          x = 0.7, gp = gpar(fontsize = 8))
pushViewport(viewport(x = unit(0.25, "npc"), # scale viewport
                      height = unit(1, "npc"),
                      width = unit(0.5, "npc")))
grid.rect(gp = gpar(bg = NA)) # scale border
grid.raster(matrix(rev(powPal(25)), ncol = 1), width = unit(1, "npc"),
            height = unit(1, "npc")) # scale fill (raster)
grid.segments(x0 = unit(rep(1, 3), "npc"), # scale ticks
              y0 = unit(seq(0, 1, by = 0.5), "npc"),
              x1 = unit(rep(1, 3), "npc") + unit(rep(0.25, 3), "lines"),
              y1 = unit(seq(0, 1, by = 0.5), "npc"),
              gp = gpar(lwd = 0.5))
grid.text(label = c(0, 0.5, 1), hjust = 0,
          y = unit(seq(0, 1, by = 0.5), "npc"),
          x = unit(rep(1, 3), "npc") + unit(rep(0.5, 3), "lines"),
          gp = gpar(fontsize = 8)) # tick labels
dev.off()

## power comparison to the previous test
png("poolPowersRandVsMn.png", width = 3, height = 3, res = 480,
    units = "in")
narrowPlot(xgrid = seq(0, 0.5, by = 0.25), ygrid = seq(0, 1, by = 0.5),
           xlab = expression({tau^2/4}), ylab = "Power")
lines(powerdf$taup[abs(powerdf$lgk - log(2,10)) < 0.01 &
                   abs(powerdf$a - 0.05) < 0.001],
      powerdf$pow[abs(powerdf$lgk - log(2,10)) < 0.01 &
                  abs(powerdf$a - 0.05) < 0.001],
      col = "firebrick")
lines(mnVSeq/4, randMnPow, col = "steelblue")
legend(x = "bottomright",
       legend = c("Classical", "Evidential"), title = "Test",
       cex = 0.8, lty = 1, col = c("steelblue", "firebrick"))
dev.off()


## REAL DATA #########################################################
library(metadat)
library(metafor)

## school calendar data
data(dat.konstantopoulos2011) # load
schoolDat <- dat.konstantopoulos2011 # rename
schoolDat$district <- factor(schoolDat$district) # convert to a factor
rm(dat.konstantopoulos2011) # clean up workspace
schoolDat <- schoolDat[order(schoolDat$district,
                             schoolDat$yi),] # sort by district

## set up parameters to plot the data
schoolCol <- RColorBrewer::brewer.pal(11, "Set3")
district <- unclass(schoolDat$district)
## determine grouped x values by district
distRle <- rle(c(district)) # run length encoding
xvals <- unlist(sapply(1:11, function(ii) { # determine spacing from
    seq(distRle$values[ii] - 0.2,           # run lengths
        distRle$values[ii] + 0.2,
        length.out = distRle$lengths[ii])
}))

## plot the data
png("metaSchoolMeans.png", width = 6, height = 3, units = "in",
    res = 480)
narrowPlot(ygrid = seq(-1.5, 1.5, by = 0.75),
           xgrid = seq(0, 12, by = 3),
           ylab = "Difference in means", xlab = "District")
for (ii in 1:nrow(schoolDat)) {
    lines(rep(xvals[ii], 2),
          rep(schoolDat$yi[ii],2) +
          c(-1.96, 1.96)*sqrt(schoolDat$vi[ii]))
          #col = schoolCol[unclass(schoolDat$district)[ii]])
}
points(xvals,
       y = schoolDat$yi, pch = 21,
       bg = schoolCol[unclass(schoolDat$district)])
dev.off()

## sweep potential combined estimate values
xseq <- seq(-2, 2, by = 0.01)
## generate p-values and combine them
pvals <- apply(qchisq(paramSweep(schoolDat$yi,
                                 sqrt(schoolDat$vi),
                                 thetas = xseq),
                      2, lower.tail = FALSE), 1,
               function(row) pchisq(sum(row), 2*nrow(schoolDat),
                                    lower.tail = FALSE))
## strong evidence of heterogeneity...

## maybe pool by district separately
schoolSplit <- split(schoolDat[, c("yi", "vi")],
                     schoolDat$district)
## apply the process to each district
districtPools <- lapply(schoolSplit,
                        function(sdat) {
                            apply(qchisq(paramSweep(sdat$yi,
                                 sqrt(sdat$vi),
                                 thetas = xseq),
                      2, lower.tail = FALSE), 1,
               function(row) pchisq(sum(row), 2*nrow(sdat),
                                    lower.tail = FALSE))
                            })

## plot the data
schoolCol <- RColorBrewer::brewer.pal(11, "Set3")
png("metaSchoolIntervals.png", width = 6, height = 3, units = "in",
    res = 480)
narrowPlot(ygrid = seq(-1.5, 1.5, by = 0.75),
           xgrid = seq(0, 12, by = 3),
           ylab = "Difference in means", xlab = "District")
for (ii in 1:nrow(schoolDat)) {
    lines(rep(xvals[ii], 2),
          rep(schoolDat$yi[ii],2) +
          c(-1.96, 1.96)*sqrt(schoolDat$vi[ii]),
          col = schoolCol[unclass(schoolDat$district)[ii]])
}
points(xvals,
       y = schoolDat$yi, pch = 20,
       col = schoolCol[unclass(schoolDat$district)])
for (ii in 1:length(districtPools)) {
    temp <-  safeRange(xseq[(districtPools[[ii]] >= 0.05)])
    tempch <- if (is.na(temp)[1]) 4 else 0
    points(y = xseq[which.max(districtPools[[ii]])],
           x = ii + 0.3, pch = tempch, col = "black")
    lines(x = rep(ii + 0.3, 2), y = temp, col = "black")
}
dev.off()
