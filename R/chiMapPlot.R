source("poolFunctions.R")

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

## helper to make plotting easier
powerHeatMap <- function(mat, inds, main,
                         pal = hcl.colors(length(ks)/2,
                                          palette = "Temps"),
                         legendLabs = c(-8, -4, 0, 4, 8),
                         legendTitle = expression(paste("log",
                                                        kappa))) {
    image(mat[inds,], xaxt = "n", yaxt = "n", col = pal,
          ylab = "", xlab = "", main = "")
    mtext(main, side = 3, line = 0, cex = 0.8) # main
    mtext(expression(rho), side = 2, line = 1, cex = 0.8) # ylab
    mtext("logD(a,w)", side = 1, line = 1, padj = 0, cex = 0.8) # xlab
    ## add ticks
    mtext(side = 1, at = seq(0, 1, by = 0.25), text = "|", line = 0,
          cex = 0.5, padj = -2)
    mtext(text = seq(-5, 5, by = 2.5), at = seq(0, 1, by = 0.25),
          side = 1, cex = 0.8)
    mtext(side = 2, at = seq(0, 1, by = 0.2), text = "|", line = 0,
          cex = 0.5, padj = 1)
    mtext(text = seq(0, 1, by = 0.2), at = seq(0, 1, by = 0.2),
          side = 2, cex = 0.8)
    bds <- par()$usr
    rasterImage(as.raster(matrix(rev(pal), ncol = 1)),
                bds[2] + 0.02, 0.2, bds[2] + 0.05, 0.8, xpd = NA)
    rect(bds[2] + 0.02, 0.2, bds[2] + 0.05, 0.8, xpd = NA)
    text(x = bds[2] + 0.035, y = 0.82, xpd = NA, cex = 0.8,
         labels = legendTitle, adj = c(0.5, 0.5))
    text(x = rep(bds[2] + 0.05, length(legendLabs)), xpd = NA,
         y = seq(0.2, 0.8, length.out = length(legendLabs)),
         labels = legendLabs, cex = 0.8,
         adj = c(-0.3,0.5))    
}

## CONSTANTS #########################################################
## necessary number of trials for a given standard error
ster <- 0.005
nsim <- (0.5/ster)^2

## read in this result
chiPowers <- readRDS("chiPowersMap.Rds")
## make into single data frame
chiPowdf <- cbind(chiPowers$pars, power = chiPowers$chi)
## simplest map: average the power by prop and log divergence
powAve <- tapply(chiPowdf$power,
                 list(logD = round(chiPowdf$logD, 2),
                      m1 = chiPowdf$m1,
                      logk = log(chiPowdf$k)),
                 mean)
## other maps: by each value of w
powByW <- lapply(split(chiPowdf, chiPowdf$w),
                 function(df) {
                     tapply(df$power,
                            list(logD = round(df$logD, 2),
                                 m1 = df$m1,
                                 logk = log(df$k)),
                            mean)
                 })

## store together in one list
fullList <- c(powByW, list(powAve))
names(fullList) <- c(log(as.numeric(names(powByW))),
                     "Mean")

## plot each case
for (logw in names(fullList)) {
    targMat <- fullList[[logw]]

    ## compute some different maxima
    minMax <- apply(targMat, c(1,2), leastMax)
    maxMax <- apply(targMat, c(1,2), largestMax)
    powRngs <- apply(targMat, c(1,2), function(p) diff(range(p)))
    ## clean this up...
    logDNames <- as.character(seq(-5, 5, by = 0.25))
    logDNames <- logDNames[logDNames %in% dimnames(minMax)$logD]

    ## store in a clean matrix
    ks <- dimnames(targMat)$logk
    minMat <- matrix(as.numeric(ks[minMax]),
                     nrow = nrow(minMax),
                     dimnames = dimnames(minMax))
    maxMat <- matrix(as.numeric(ks[maxMax]),
                     nrow = nrow(maxMax),
                     dimnames = dimnames(maxMax))
    maxMat[powRngs < 0.02] <- NA
    minMat[powRngs < 0.02] <- NA
    meanMat <- (maxMat + minMat)/2

    ## plot the heatmaps of these
    par(mar = c(2.1, 2.1, 1.1, 3.1), mfrow = c(1,2))
    powerHeatMap(minMat, inds = logDNames,
                 main = expression(paste("Smallest ", kappa,
                                         " giving maximum power")))
    powerHeatMap(meanMat, inds = logDNames,
                 main = expression(paste("Mean ", kappa,
                                         " giving maximum power")))
    #powerHeatMap(powRngs, inds = logDNames,
    #             main = expression(paste("Range of powers in ", kappa)),
    #             legendLabs = c(0, 1), legendTitle = "",
    #             pal = colorRampPalette(c("white",
    #                                      "firebrick"))(length(ks)/2))
}

## filled contours instead?
filled.contour(minMat[logDNames,],
               color.palette = function(n) hcl.colors(n))
contour(minMat[logDNames,], xaxs = "i", yaxs = "i")
.filled.contour(x = seq(0,1, length.out = 41),
                y = seq(0,1, length.out = 41),
                z = mapMat[logDNames,],
                levels = seq(-8,8,by = 1), col = mapPal)
contour(mapMat[logDNames,], add = TRUE, col = "gray90")
