library(pooledCentR)

## FUNCTIONS #########################################################
## binomial difference function comparing an observed power to the
## maximum value observed
binDiff <- function(p, mx, crit = -1.96) {
    pc <- (p + mx)/2 # pooled power
    notless <- (p - mx)/sqrt(pc*(1-pc)*2/nsim) > crit
    notless[is.na(notless)] <- TRUE # handle zeros
    notless
}

## function to convolve two matrices, the default uses a 5x5 gaussian
## filter to reduce noise
smoothPower <- function(x, filt = gFilt) {
    dims <- dimnames(x) # before changing dimensions
    nullCol <- x[,1]
    x <- x[, -1] # insert null values later
    nr <- nrow(x)
    nc <- ncol(x) # new dimensions
    xpd <- rbind(x[c(2,1), ], x, x[c(nr, nr-1),])
    xpd <- cbind(xpd[, c(2,1)], xpd, xpd[, c(nc, nc-1)]) # x padded
    conv <- matrix(NA, nrow = nr, ncol = nc)
    for (ii in (1:nr)) { # convolution steps
        for (jj in (1:nc)) {
            ix <- ii + 2
            jx <- jj + 2
            rinds <- (ix-2):(ix+2)
            cinds <- (jx-2):(jx+2)
            conv[ii,jj] <- sum(xpd[rinds, cinds]*filt)
        }
    }
    conv <- cbind(nullCol, conv) # null case
    dimnames(conv) <- dims
    conv
}

## wrapper to apply above to array
smoothParray <- function(arr) {
    dims <- dimnames(arr) # keep dimnames for later
    smth <- simplify2array(
        lapply(dims[[3]], function(k) smoothPower(arr[,,k])))
    dimnames(smth) <- dimnames(arr)
    smth
}

## mask a matrix of values using a logical matrix
masker <- function(x, mask) {
    x[!mask] <- 0
    x[is.na(x)] <- 0
    x
}

## function to aggregate matrices safely
accumMat <- function(m1, m2) {
    m1[rownames(m2), colnames(m2)] <- m1[rownames(m2),
                                         colnames(m2)] + m2
    m1
}

## and 3d arrays
accumArr <- function(a1, a2) {
    dims2 <- dimnames(a2)
    a1[dims2[[1]], dims2[[2]], dims2[[3]]] <- a1[dims2[[1]],
                                                 dims2[[2]],
                                                 dims2[[3]]] + a2
    a1
}

## mask arrays in our setting
maskArr <- function(arr, mask) {
    for (kk in seq_len(dim(arr)[[3]])) {
        tempMat <- arr[,,kk]
        tempMat[!mask] <- 0
        arr[,,kk] <- tempMat
    }
    arr
}

## PRELIMINARIES #####################################################
## define convolutional gaussian filter for smoothing
library(mvtnorm)
## set grid
lowy <- rep(seq(-2.5, 1.5, by = 1), times = 5)
lowx <- rep(seq(-2.5, 1.5, by = 1), each = 5)
upy <- rep(seq(-1.5, 2.5, by = 1), times = 5)
upx <- rep(seq(-1.5, 2.5, by = 1), each = 5)
## compute integrals over appropriate regions
gaussFilt <- matrix(sapply(1:length(lowx),
                           function(ii) pmvnorm(lower = c(lowx[ii],
                                                          lowy[ii]),
                                                upper = c(upx[ii],
                                                          upy[ii]))),
                    ncol = 5)
gFilt <- gaussFilt/sum(gaussFilt)


## LOADING ###########################################################
## read in the power data from a large simulation
#dataFile <- "./results/chiPowersMap80.Rds" # log(w)
dataFile <- "./results/chiPowersMap80_unifW.Rds"
chiPowers <- readRDS(dataFile)
## make into single data frame
powdf <- cbind(chiPowers$pars, power = chiPowers$chi)


## CONSTANTS #########################################################
## necessary number of trials for a given standard error
ster <- 0.005
nsim <- (0.5/ster)^2

## logD and kappa sequences in the data
logDseq <- seq(-5, 5, by = 0.125)
kapSeq <- seq(-8, 8, by = 0.25)

## null matrix and array for aggregation (ensure constant plots)
nullMat <- matrix(0, ncol = 81, nrow = 81,
                  dimnames = list("logD" = as.character(logDseq),
                                  "m1" = 0:80))
nullArr <- array(0, dim = c(81, 81, length(kapSeq)),
                 dimnames = list("logD" = as.character(logDseq),
                                 "m1" = 0:80,
                                 "logk" = as.character(kapSeq)))

## PROCESSING: USED TO GENERATE DATA MASKED MATRICES FOR PACKAGE #####
## clean up the log divergences
powdf$logD <- round(powdf$logD, 3)
powRegD <- powdf[powdf$logD %in% as.character(logDseq),]

## split the power data by logw and compute each case separately
pow_byW <- split(powRegD[, c("m1", "logD", "k", "power")],
                 if (grepl("unifW", dataFile)) {
                     powRegD$w
                 }else log(powRegD$w))
## arrange each case as an array
pow_byW <- lapply(pow_byW,
                  function(df) {
                      tapply(df$power,
                             list(logD = round(df$logD, 3),
                                  m1 = df$m1,
                                  logk = log(df$k)),
                             mean) # identity, one value
                  })
## smooth the powers
smth_byW <- lapply(pow_byW, smoothParray)
## get the max and min powers
max_byW <- lapply(smth_byW,
                  function(arr) apply(arr, c(1,2), max))
min_byW <-  lapply(smth_byW,
                   function(arr) apply(arr, c(1,2), min))
## check the max against each value
sameMax_byW <- mapply(function(smth, mx) sweep(smth, c(1,2),
                                               mx, binDiff),
                      smth_byW, max_byW)

## reduce this
sameMax <- Reduce(accumArr, sameMax_byW, init = nullArr)

## two possible masks seem reasonable
## the range/sd/change in power
sd_byW <- lapply(pow_byW,
                 function(arr) apply(arr, c(1,2), sd))
## the number of matching methods at max power
match_byW <- lapply(sameMax_byW,
                    function(arr) apply(arr, c(1,2), sum))
## get individual w masks based on thresholds
sdLim <- ster*2
matchLim <- 20
sdmask_byW <- lapply(sd_byW, function(mat) mat > sdLim)
mtmask_byW <- lapply(match_byW, function(mat) mat < matchLim)
## count cases for all
caseMat <- Reduce(accumMat, mtmask_byW, init = nullMat)
maskMat <- Reduce(accumMat, sdmask_byW, init = nullMat) > 0

## use the masks on the sameMax arrays
sameMaxMask_byW <- mapply(maskArr, sameMax_byW, mtmask_byW)
## reduce this
sameMaxMask <- Reduce(accumArr, sameMaxMask_byW, nullArr)
## scale this
maxProp <- sweep(sameMaxMask, c(1,2), caseMat, `/`)

## PLOTS #############################################################
##' everything up until here is effectively used to generate the data
##' used in the marbR package for the altFrequencyMat function, and
##' the display and final tidying of the data has been completed by
##' functions in the package

## plots for the paper (Figs 4.15 and 4.16)
## -8 -8
## -4 -4
## -1 -1
## 0 0
## log(2) log(2)
## 1 1
## 4 4
## 8 8
## select indices by kappa
## clean up name variables
kpnm <- as.character(kaps)
if (length(unique(kpnm)) == 1) kpnm <- kpnm[1]
kpnm <- gsub("\\.", "_", kpnm)
suffix <- if (grepl("unifW", dataFile)) "unifW" else ""
## plot the regions
png(paste0("regionPlot", paste(kpnm, collapse = "-"), suffix, ".png"),
    width = 3, height = 3, units = "in", res = 480)
par(mar = c(2.1, 2.1, 1.5, 1.5))
marHistHeatMap(altFrequencyMat(c(log(2), log(2))))
abline(h = seq(0, 1, by = 0.2), v = seq(0, 1, by = 0.25),
       col = adjustcolor("grey50", 0.5), lty = 2)
dev.off()
