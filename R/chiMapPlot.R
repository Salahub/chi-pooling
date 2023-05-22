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
powerHeatMap <- function(mat, main = "",
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

## helper to make plotting easier
alternativeHeatMap <- function(mat, main = "", pal = NULL) {
    if (is.null(pal)) {
        pal <- colorRampPalette(c("white",
                                  "firebrick"))(max(mat,
                                                    na.rm = TRUE) + 1)
    }
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
    rowDist <- rowSums(mat, na.rm = TRUE)
    colDist <- colSums(mat, na.rm = TRUE) # marginal distributions
    vboxBds <- c(bds[2], bds[2] + 0.1, bds[3:4])
    hboxBds <- c(bds[1:2], bds[4], bds[4] + 0.1)
    ## density boxes
    rect(vboxBds[1], vboxBds[3], vboxBds[2], vboxBds[4], xpd = NA)
    rect(hboxBds[1], hboxBds[3], hboxBds[2], hboxBds[4], xpd = NA)
    ## add marginal histograms
    vseq <- seq(vboxBds[3], vboxBds[4],
                length.out = length(colDist) + 1)
    rect(vboxBds[1], vseq[1:length(colDist)],
         vboxBds[1] + 0.9*diff(vboxBds[1:2])*(colDist/max(colDist)),
         vseq[2:(length(colDist) + 1)], xpd = NA,
         col = adjustcolor("firebrick", 0.5))
    hseq <- seq(hboxBds[1], hboxBds[2],
                length.out = length(rowDist) + 1)
    rect(hseq[1:length(rowDist)], hboxBds[3],
         hseq[2:(length(rowDist) + 1)], xpd = NA,
         hboxBds[3] + 0.9*diff(hboxBds[3:4])*(rowDist/max(rowDist)),
         col = adjustcolor("firebrick", 0.5))
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
powdf$logD <- round(powdf$logD, 2)
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
tol <- 8
masker <- function(x, y, tol = 8) {
    x[y > tol] <- 0
    x[is.na(x)] <- 0
    x
}
## apply to data
pow_minMask <- mapply(masker, x = pow_minMats, y = pow_gaps,
                      tol = tol)
pow_maxMask <- mapply(masker, x = pow_maxMats, y = pow_gaps,
                      tol = tol)

## aggregate a total mask: cases where it is never interesting
accum <- function(m1, m2) {
    m1[rownames(m2), colnames(m2)] <- m1[rownames(m2),
                                         colnames(m2)] + m2
    m1
}
maskMat <- Reduce(accum, lapply(pow_gaps,
                                function(mat) {
                                    mat[is.na(mat)] <- tol + 1
                                    mat <= tol
                                }))
maskMat <- maskMat > 0

## take a kappa
kap <- 7
## see where it fits
kapMax <- mapply(function(m1, m2, k) m1 <= k & m2 > k,
                 pow_minMask, pow_maxMask, kap)
## reduce to a single matrix
kapMaxDist <- Reduce(accum, kapMax,
                     init = matrix(0, nrow = nrow(kapMax[[1]]),
                                   ncol = ncol(kapMax[[1]]),
                                   dimnames = dimnames(kapMax[[1]])))
kapMaxMask <- kapMaxDist
kapMaxMask[!maskMat] <- NA
## plot it
png(paste0("regionPlot", kap, ".png"), width = 3.5, height = 3.5,
    units = "in", res = 240)
par(mar = c(2.1, 2.1, 2.8, 1.5))
alternativeHeatMap(kapMaxMask, main = "")
mtext(bquote("Alternatives for log"*kappa==.(kap)),
      line = 1.5, cex = 0.8)
abline(h = seq(0, 1, by = 0.2), v = seq(0, 1, by = 0.25),
       col = adjustcolor("grey50", 0.5), lty = 2)
dev.off()

## plot the heatmaps of these
ind <- 1
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

## report tests most contributing: max proportion*M largest statistics
## for the minimum kappa
