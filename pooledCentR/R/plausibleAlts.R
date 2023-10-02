## get plausible alternative regions
altFrequencyMat <- function(kappaRange, logW = FALSE) {
    if (logW) {
        maskMat <- kappaPowersMaskedLogW
    } else {
        maskMat <- kappaPowersMasked
    }
    kapSeq <- as.numeric(dimnames(maskMat)$logk)
    kapInd <- kapSeq <= kappaRange[2] & kapSeq >= kappaRange[1]
    if (kappaRange[2] > max(kapSeq) | kappaRange[1] < min(kapSeq)) {
        warning(cat("kappaRange extends beyond the simulation range",
                    "(exp(-8) to exp(8)), map reported is only",
                    "accurate for values within this range"))
    }
    if (!any(kapInd)) { # if interval is between simulation steps
        kapInd <- logical(length(kapSeq))
        kapInd[which.min(abs(kapSeq - mean(kaps)))] <- TRUE
    }
    ## aggregate
    apply(maskMat[,,kapInd], c(1,2), sum)
}

## helper to make plotting easier
marHistHeatMap <- function(mat, main = "", pal = NULL, ...) {
    if (is.null(pal)) {
        pal <- colorRampPalette(c("white",
                                  "firebrick"))(20)
    }
    image(mat, xaxt = "n", yaxt = "n", col = pal, ylab = "", xlab = "",
          main = main, ...)
    mtext(main, side = 3, line = 0, cex = 0.8) # main
    mtext(expression(eta), side = 2, line = 1, cex = 0.8) # ylab
    mtext("lnD(a,w)", side = 1, line = 1, padj = 0, cex = 0.8) # xlab
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
