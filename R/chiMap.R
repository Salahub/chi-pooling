source("poolFunctions.R")

## CONSTANTS #########################################################
## necessary number of trials for a given standard error
ster <- 0.005
nsim <- (0.5/ster)^2

## given a w seq, the corresponding a values to span D(a,w)
wDsets <- expand.grid(w = exp(seq(-6, 0, by = 1)), # w values
                      logD = seq(-5, 5, by = 0.25)) # target logDivs
aseq <- apply(wDsets, 1, # find a values given w, logDiv
              function(row) findA(row[1], row[2],
                                  tol = .Machine$double.eps))
logD <- log(betaDiv(aseq, wDsets[,1])) # compute actual logDiv

## a sequence of M values, k values, a and w values to parameterize
## data generation
kseq <- exp(seq(-8, 8, by = 0.5)) # k values for chi
Mval <- 40 # total number of tests
m1seq <- seq(0, Mval, by = 1)
preparams <- data.frame(a = rep(aseq, times = length(m1seq)),
                        w = rep(wDsets[,"w"], times = length(m1seq)),
                        m1 = rep(m1seq, each = length(aseq)),
                        logD = rep(logD, times = length(m1seq)))
params <- data.frame(preparams, k = rep(kseq, each = nrow(preparams)))
altSeq <- lapply(seq_len(nrow(params)), # generate alternatives
                 function(ii) {
                     with(params, pgenMix(a[ii], w[ii], m1[ii]))
                 })
poolSeq <- lapply(seq_len(nrow(params)), # generate pool functions
                  function(ii) {
                      function(p) poolChi(p, k = params$k[ii])
                  })
## control data set generation with random seeds
seeds <- round(runif(nrow(params), min = 1e5, max = 1e9))
## helper to modify powersim
chiSimPower <- function(altLst, poolLst, nsim, seeds, M) {
    mapply(powerSim, pcomb = poolLst, pgen = altLst,
           nsim = nsim, seed = seeds, M = M)
}

## get the map for the chi method (requires parallel computing)
library(parallel)
ncores <- detectCores()/2
spltInds <- rep(1:ncores, # indices for splitting functions
                each = ceiling(nrow(params)/ncores))[1:nrow(params)]
sdLst <- split(seeds, spltInds) # split seeds
altLst <- split(altSeq, spltInds)
poolLst <- split(poolSeq, spltInds) # functions separated by indices
powerschi <- mcmapply(chiSimPower, altLst = altLst, poolLst = poolLst,
                      nsim = nsim, seeds = sdLst, M = Mval,
                      mc.cores = ncores)

## store the simulation results
chiPowers <- list(pars = params,
                  chi = unlist(powerschi))
saveRDS(chiPowers, "chiPowersMap.Rds")

## read in this result
chiPowers <- readRDS("chiPowersMap.Rds")
## make into single data frame
chiPowdf <- cbind(chiPowers$pars, power = chiPowers$chi)
## simplest map: average the power by prop and log divergence
chiPowAve <- tapply(chiPowdf$power,
                    list(logD = round(chiPowdf$logD, 2),
                         m1 = chiPowdf$m1,
                         logk = log(chiPowdf$k)),
                    mean)

## get the max for each, compute z-scores, take the smallest not
## significantly different from the max
leastMax <- function(p) {
    mx <- max(p)
    pc <- (p + mx)/2
    eqmax <- (p - mx)/sqrt(pc*(1-pc)*2/nsim) > -1.96
    eqmax[is.na(eqmax)] <- TRUE # fix 1, 1 case
    first <- which(eqmax)[1]
    if (all(eqmax)) NA else first
}

largestMax <- function(p) {
    mx <- max(p)
    pc <- (p + mx)/2
    eqmax <- (p - mx)/sqrt(pc*(1-pc)*2/nsim) > -1.96
    eqmax[is.na(eqmax)] <- TRUE # fix 1, 1 case
    last <- which(eqmax)[sum(eqmax)]
    if (all(eqmax)) NA else last
}

chiPowLstMax <- apply(chiPowAve, c(1,2), leastMax)
chiPowLgstMax <- apply(chiPowAve, c(1,2), largestMax)
chiPowRng <- apply(chiPowAve, c(1,2), function(p) diff(range(p)))
## clean this up...
logDNames <- as.character(seq(-5, 5, by = 0.25))
## store in a clean matrix
ks <- dimnames(chiPowAve)$logk
mapMat <- matrix(as.numeric(ks[chiPowLstMax]),
                 nrow = nrow(chiPowLstMax),
                 dimnames = dimnames(chiPowLstMax))
mapMat[chiPowRng < 0.01] <- NA

## plot as a heat map
mapPal <- hcl.colors(length(ks)/2)
par(mar = c(2.1, 2.1, 1.1, 2.2))
image(mapMat[logDNames,], xaxt = "n", yaxt = "n", col = mapPal,
      ylab = "", xlab = "", main = "")
mtext(expression(paste("Most powerful ", kappa, " by region")),
      side = 3, line = 0, cex = 0.8) # main
mtext(expression(rho), side = 2, line = 1, cex = 0.8) # ylab
mtext("logD(a,w)", side = 1, line = 1, padj = 0, cex = 0.8) # xlab
## add ticks
mtext(side = 1, at = seq(0, 1, by = 0.25), text = "|", line = 0,
      cex = 0.5, padj = -2)
mtext(text = seq(-5, 5, by = 2.5), at = seq(0, 1, by = 0.25),
      side = 1, cex = 0.8)
mtext(side = 2, at = seq(0, 1, by = 0.2), text = "|", line = 0,
      cex = 0.5, padj = 1)
mtext(text = seq(0, 1, by = 0.2), at = seq(0, 1, by = 0.2), side = 2,
      cex = 0.8)
bds <- par()$usr
rasterImage(as.raster(matrix(rev(mapPal), ncol = 1)),
            bds[2] + 0.02, 0.2, bds[2] + 0.05, 0.8, xpd = NA)
rect(bds[2] + 0.02, 0.2, bds[2] + 0.05, 0.8, xpd = NA)
text(x = bds[2] + 0.035, y = 0.82, xpd = NA, cex = 0.8,
     labels = expression(paste("log", kappa)), adj = c(0.5, 0.5))
text(x = rep(bds[2] + 0.05, 5), y = c(0.2, 0.35, 0.5, 0.65, 0.8),
     labels = c(-8, -4, 0, 4, 8), cex = 0.8, xpd = NA,
     adj = c(-0.3,0.5))
## filled contours instead?
filled.contour(mapMat[logDNames,],
               color.palette = function(n) hcl.colors(n))
contour(mapMat[logDNames,], xaxs = "i", yaxs = "i")
.filled.contour(x = seq(0,1, length.out = 41),
                y = seq(0,1, length.out = 41),
                z = mapMat[logDNames,],
                levels = seq(-8,8,by = 1), col = mapPal)
contour(mapMat[logDNames,], add = TRUE, col = "gray90")
