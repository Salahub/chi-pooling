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
powerHeatMap <- function(mat, main,
                         pal = hcl.colors(17,
                                          palette = "Temps"),
                         legendLabs = c(-8, -4, 0, 4, 8),
                         legendTitle = expression(paste("log",
                                                        kappa))) {
    image(mat, xaxt = "n", yaxt = "n", col = pal,
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
powdf <- cbind(chiPowers$pars, power = chiPowers$chi)

## clean up
powdf$logD <- round(chiPowdf$logD, 2)
powRegD <- powdf[powdf$logD %in% as.character(seq(-5, 5, by = 0.25)),]

## split by logw and compute each case separately
pow_byW <- split(powRegD[, c("m1", "logD", "k", "power")],
                 log(powRegD$w))
## arrange as arrays
pow_byW <- lapply(pow_byW,
                  function(df) {
                      tapply(df$power,
                             list(logD = round(df$logD, 2),
                                  m1 = df$m1,
                                  logk = log(df$k)),
                             mean) # identity, one value
                  })

## get max and min powers
pow_max <- lapply(pow_byW,
                  function(mt) apply(mt, c(1,2), max))
pow_min <- lapply(pow_byW,
                  function(mt) apply(mt, c(1,2), min))
pow_rng <- mapply(`-`, pow_max, pow_min)
## indices of matches and gap giivng number of max matches
pow_minMax <- lapply(pow_byW,
                     function(mt) apply(mt, c(1,2), leastMax))
pow_maxMax <- lapply(pow_byW,
                     function(mt) apply(mt, c(1,2), largestMax))
pow_gaps <- mapply(`-`, pow_maxMax, pow_minMax)

## convert back to parameter values
ks <- lapply(pow_byW, function(ar) dimnames(ar)$logk)
pow_minMats <- mapply(function(ks, mat) {
    matrix(as.numeric(ks[mat]), nrow = nrow(mat),
           dimnames = dimnames(mat)) },
    ks, pow_minMax)
pow_maxMats <- mapply(function(ks, mat) {
    matrix(as.numeric(ks[mat]), nrow = nrow(mat),
           dimnames = dimnames(mat)) },
    ks, pow_maxMax)

## mask the "uninteresting cases"
pow_minMask <- mapply(function(x, y, tol) {x[y > tol] <- NA; x},
                      pow_minMats, pow_gaps, 8)
pow_maxMask <- mapply(function(x, y, tol) {x[y > tol] <- NA; x},
                      pow_maxMats, pow_gaps, 8)

## take a kappa
kap <- 3
## see where it fits
kapMax <- mapply(function(m1, m2, k) m1 <= k & m2 >= k,
                 pow_minMask, pow_maxMask, kap)
kapMax[gapSize > 15] <- NA
#kapMax[powMin > 0.9 | powMax < 0.1] <- NA
image(kapMax[[1]])

## plot the heatmaps of these
ind <- 7
par(mar = c(2.1, 2.1, 1.1, 3.1), mfrow = c(1,2))
powerHeatMap(pow_minMask[[ind]],
             main = expression(paste("Smallest ", kappa,
                                     " giving maximum power")))
powerHeatMap(pow_maxMask[[ind]],
             main = expression(paste("Largest ", kappa,
                                     " giving maximum power")))
powerHeatMap(pow_rng[[ind]],
             main = expression(paste("Range of powers in ", kappa)),
             legendLabs = c(0, 1), legendTitle = "",
             pal = colorRampPalette(c("white",
                                      "firebrick"))(17))

## new idea: take a kappa/range of kappas and produce the region where
## it isn't different than the maximum power across all kappas, plot
## this as a suggestive visual of likely alternatives

## in support of that:
## - prove that increasing M only increases the resolution of the plot
## - increase the resolution and see what it looks like

## filled contours instead?
filled.contour(minMat[logDNames,],
               color.palette = function(n) hcl.colors(n))
contour(minMat[logDNames,], xaxs = "i", yaxs = "i")
.filled.contour(x = seq(0,1, length.out = 41),
                y = seq(0,1, length.out = 41),
                z = mapMat[logDNames,],
                levels = seq(-8,8,by = 1), col = mapPal)
contour(mapMat[logDNames,], add = TRUE, col = "gray90")
