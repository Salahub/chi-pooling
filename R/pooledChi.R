source("poolFunctions.R")

## SETUP :: visualization parameters
library(RColorBrewer)
kpal <- colorRamp( # define a simple palette for k values
    brewer.pal(3, "Dark2")[c(2,3)])
chiCol <- "slategray" # a colour for chi curves

## SETUP :: function to make a switch for temporary parameters
makeSwitch <- function(input, output) {
    function(x) {
        inds <- apply(outer(x, input,
                            function(x,y) { # fuzzy matching
                                abs(x-y) < .Machine$double.eps
                            }),
                      1, function(row) {
                          tr <- which(row) # match inds
                          if (length(tr) == 0) tr <- NA # NA if empty
                          tr
                      })
        output[inds] # match to output table
    }
}


## SIMULATION ########################################################
## compute the powers when H4 is true under correct
## specification, misspecification, or the chi method
## COMPUTE :: necessary number of trials for a given standard error
ster <- 0.005
nsim <- (0.5/ster)^2

## LOAD :: simulated data or else run simulation
readData <- TRUE
datafile <- "constMPowers_fullgrid.Rds"
## datafile <- "constMPowers.Rds"
if (readData) { # load specified file
    tempDat <- readRDS(datafile)
    kseq <- tempDat$ks # chi k values in simulation
    params <- tempDat$pars
    misSpec <- tempDat$misw
    powershr <- tempDat$tw
    powerschi <- tempDat$chi
    powersmis <- tempDat$mis
    rm(tempDat)
} else {
    if (grepl("fullgrid", datafile)) {
        source("simulateConstM_fullgrid.R")
    } else {
        source("simulateConstM.R")
    }
}

######################################################################


## PLOTS CURRENTLY IN THE PAPER ######################################
## SETUP :: compute the KL divergence of each setting
klDivs <- betaDiv(params$a, params$w)
## the temporary lty, pch, and color switches
acol <- makeSwitch(c(0.05, 0.5, 0.95),
                   c("firebrick", "seagreen", "steelblue"))
apch <- makeSwitch(c(0.05, 0.5, 0.95), c(21, 24, 23))
wlty <- makeSwitch(misSpec, c(1,2,3,6))
Mlwd <- makeSwitch(c(2, 5, 10, 20, 50), c(1, 2, 3, 4, 5))
wcol <- makeSwitch(c(unique(params$w), 0.5),
                   c(brewer.pal(length(unique(params$w)), "Set2"),
                     "cadetblue2"))
wpch <- makeSwitch(unique(params$w), c(15,12,13,16,18,14,17))
## establish unique values
if (!grepl("fullgrid", datafile)) {
    uniqueAW <- expand.grid(a = unique(params$a),
                            w = unique(params$w))
    uniqueKL <- betaDiv(uniqueAW$a, uniqueAW$w)
} else {
    uniqueAW <- params[params$M == 2, c("a", "w", "logD")]
    uniqueKL <- exp(uniqueAW$logD)
}
## get the indices of M = 2 and M = 20 for later plots
M2inds <- which(params$M == 2)
M20inds <- which(params$M == 20)

## IMAGE :: KL divergence with inset densities by w
dispLog <- rep(rep(c(TRUE, FALSE), each = 7), # display only a subset
               length.out = nrow(uniqueAW))
png("klDivs.png", width = 4, height = 4, units = "in", res = 400)
narrowPlot(xgrid = seq(-6, 0, by = 2), # plot area
           mars = c(2.1, 2.1, 1.1, 1.1), ylim = c(-5.1,5.1),
           xlim = c(-6.1,0.1), ygrid = seq(-5, 5, length.out = 5),
           xlab = "ln(w)", ylab = "ln(D(a,w))",
           main ="Beta densities by w and D(a,w)")
tempRat <- aspectRatio() # capture device aspect ratio
for (ii in seq_along(uniqueAW$w)) { # plot glyphs, aspect ratio = 1
    if (dispLog[ii]) { # & log(uniqueAW$w[ii]) %in% c(-6, -3, -1, 0)) {
        densityGlyph(function(x) dbeta(x, uniqueAW$a[ii],
                                       bwa(uniqueAW$w[ii],
                                           uniqueAW$a[ii])),
                     xseq = seq(0, 1, length.out = 100),
                     x = log(uniqueAW$w[ii]),
                     y = log(uniqueKL[ii]),
                     devRat = tempRat, size = 0.22,
                     col = adjustcolor(wcol(uniqueAW$w[ii]), 0.5))
    }
}
dev.off()

## IMAGE :: the distribution of lw for different w values
set.seed(1782491)
wVals <- c(exp(-6), exp(-3), 0.5, 1)
lwDat <- lapply(wVals,
                function(w) matrix(runif(1e5*10), ncol = 10))
lwDists <- mapply(function(dt, w) apply(dt, 1, hrBeta, w = w),
                  lwDat, wVals)
## plot densities of these
reprs <- c(1, 2, 3, 4) # select densities
lwmax <-  max(abs(lwDists)) # limits
png("lhrDens.png", width = 3.5, height = 2.5, units = "in",
    res = 400)
narrowPlot(xgrid = seq(-20, 20, by = 20), xlim = c(-lwmax, lwmax),
           ygrid = seq(0, 0.15, by = 0.05), ylab = "Density",
           xlab = expression(l[HR]),
           main = expression(paste("Distribution of ", l[HR],
                                   " by ", w)))
## add the densities over top
for (ii in reprs) {
    tempdens <- density(lwDists[,ii])
    polygon(tempdens, col = adjustcolor("gray", 0.3),
            border = "black", lty = wlty(wVals[ii]))
}
## add a legend
legend(x = "topright", lty = wlty(wVals),
       legend = c(expression(e^-6), expression(e^-3), 1/2, 1),
       title = "w", bg = NA, bty = "n", cex = 0.8)
dev.off()

## IMAGE :: show the impact of M on power for every KL div
png("klDivPowByM.png", width = 3, height = 2.6, units = "in", res = 400)
narrowPlot(ygrid = seq(0.2, 1, by = 0.2), xlab = "ln(D(a, w))",
           xgrid = seq(-5, 5, length.out = 5), ylab = "Power",
           ylim = c(0.05, 1.02), xlim = c(-5, 5),
           mars = c(2.1, 2.1, 1.1, 3.1),
           main = expression(paste("Power of HR(", bold(p), ", w)",
                                   " by D(a, w) and M")))
for (ii in seq_along(M2inds)) { # lines connecting M = 2 to M = 20
    if (params$w[ii] %in% exp(c(-6, 0))) {
        lines(x = c(log(klDivs[M2inds[ii]]), log(klDivs[M20inds[ii]])),
              y = c(powershr[M2inds[ii]], powershr[M20inds[ii]]),
              col = adjustcolor(wcol(params$w[M2inds[ii]]), 0.5))
        points(x = rep(log(klDivs[M2inds[ii]]), 2), pch = 16,
               y = c(powershr[M2inds[ii]], powershr[M20inds[ii]]),
               col = adjustcolor(wcol(params$w[M2inds[ii]]), 0.5),
               cex = c(0.5, 1.5))
    }
}
## cex = log(params$M[ii], base = 4), # for dynamic point size
##legend(x = 5.6, y = 0.35, legend = unique(params$M), xpd = NA,
##       col = adjustcolor("gray", 0.6), pch = 16, bty = "n",
##       cex = 0.8, pt.cex = log(unique(params$M), base = 4),
##       title = "M") # legend for dynamic point size
legend(x = 5.6, y = 0.35, legend = c(2, 20), xpd = NA,
       col = adjustcolor("gray", 0.6), pch = 16, bty = "n",
       cex = 0.8, pt.cex = c(0.5, 1.5),
       title = "M") # legend for point size
legend(x = 5.5, y = 1.15, xpd = NA, bty = "n", cex = 0.8,
       legend = c(log(unique(params$w))[-length(unique(params$w))],
                  " 0"),
       fill = adjustcolor(wcol(unique(params$w)), 0.6),
       title = "ln(w)", bg = "white") # legend for w
dev.off()

## IMAGE :: UMP power curves plotted vs KL divergence by w
png("klDivPowByw.png", width = 2.6, height = 2.6, units = "in", res = 400)
narrowPlot(ygrid = seq(0.2, 1, by = 0.2), xlab = "ln(D(a, w))",
           xgrid = seq(-5, 5, length.out = 5), ylab = "Power",
           ylim = c(0.05, 1), xlim = c(-5, 5),
           main = expression(paste("HR(", bold(p), ", w)",
                                   " power curves by D(a,w)")))
for (w in unique(params$w)) { # curve for each w
    inds <- abs(params$w - w) < .Machine$double.eps &
        params$M == 2 # select cases from design mat
    points(log(klDivs[inds]), powershr[inds], cex = log(2, base = 3),
           col = adjustcolor(wcol(params$w[inds]), 0.5), pch = 20)
    lines(log(klDivs[inds]), powershr[inds],
          col = adjustcolor(wcol(w), 0.5)) # curve
}
dev.off()

## IMAGE :: mis-specified KL plots by w, shows clear ordering of lines
png("klDivPowMisSpec.png", width = 4.5, height = 3, units = "in",
    res = 400)
narrowPlot(ygrid = seq(0.2, 1, by = 0.2),
           xlab = expression(paste("ln(D(a, ", omega, "))")),
           xgrid = seq(-5, 5, length.out = 5), ylab = "Power",
           ylim = c(0.05, 1), xlim = c(-5, 5),
           mars = c(2.1, 2.1, 1.1, 4.1),
           main = expression(paste("Power of ", g[w], " by D(a, ",
                                   omega, ")")))
for (w in exp(c(-6, 0))) { # subset of w cases
    inds <- abs(params$w - w) < .Machine$double.eps &
        params$M == 2 # select from design mat
    for (jj in seq_len(ncol(powersmis))) { # compare to each misspec
        lines(log(klDivs[inds]), powersmis[inds, jj],
              col = adjustcolor(wcol(misSpec[jj]), 0.7))
        points(log(klDivs[inds]), powersmis[inds, jj],
               pch = wpch(w), cex = 0.7,
               col = adjustcolor(wcol(misSpec[jj]), 0.7))
    }
}
legend(x = 5.5, y = 1, title = "ln(w)", lty = 1, bty = "n", # lty by w
       legend = c("    -6", "    -3", "-ln(2)", "     0"),
       col = wcol(exp(c(-6, -3, -log(2), 0))), xpd = NA, cex = 0.8)
legend(x = 6.2, y = 0.4, pch = wpch(exp(c(-6,0))), xpd = NA,
       bty = "n", legend = c("-6", " 0"),
       col = adjustcolor("black", 0.5),  cex = 0.8,
       title = expression(paste("ln(", omega, ")"))) # pch by omega
dev.off()

## IMAGE :: mis-specified plots by KL divs including chi square line
png("klDivPowMisChi.png", width = 4.5, height = 3, units = "in",
    res = 400)
narrowPlot(ygrid = seq(0.2, 1, by = 0.2),
           xlab = expression(paste("ln(D(a, ", omega, "))")),
           xgrid = seq(-5, 5, length.out = 5), ylab = "Power",
           ylim = c(0.05, 1), xlim = c(-5, 5),
           mars = c(2.1, 2.1, 1.1, 4.1),
           main = expression(paste("Power of ", g[chi],
                                   " and ", g[w],
                                   " by D(a,", omega, ")")))
for (w in exp(c(-6, 0))) {
    inds <- abs(params$w - w) < .Machine$double.eps &
        params$M == 2
    for (jj in seq_len(ncol(powersmis))) {
        lines(log(klDivs[inds]), powersmis[inds, jj],
              col = adjustcolor(wcol(misSpec[jj]), 0.7))
        points(log(klDivs[inds]), powersmis[inds, jj],
               pch = wpch(w), cex = 0.7,
               col = adjustcolor(wcol(misSpec[jj]), 0.7))
    }
    lines(log(klDivs[inds]), powerschi[inds, 6], lwd = 1,
          col = adjustcolor(chiCol)) # add chi method line
}
legend(x = 5.5, y = 1, title = "ln(w)", lty = 1, bty = "n",
       legend = c("    -6", "    -3", "-ln(2)", "     0"),
       col = wcol(exp(c(-6, -3, -log(2), 0))), xpd = NA, cex = 0.8)
legend(x = 6.2, y = 0.4, pch = wpch(exp(c(-6,0))), xpd = NA,
       bty = "n", legend = c("-6", " 0"),
       col = adjustcolor("black", 0.5),  cex = 0.8,
       title = expression(paste("ln(", omega, ")")))
dev.off()

######################################################################


## OLD PLOTS FOR REFERENCE, POTENTIAL USE ############################

## IMAGE :: KL divergence with inset densities by w and a
dispLog <- rep(rep(c(TRUE, FALSE), each = 7), # display only a subset
               length.out = nrow(uniqueAW))
png("klDivs_wa.png", width = 5, height = 4, units = "in", res = 400)
narrowPlot(xgrid = seq(0, 1, length.out = 5), # plot area
           mars = c(2.1, 2.1, 1.1, 3.1), ylim = c(-5.1,5.1),
           ygrid = seq(-5, 5, length.out = 5), xlab = "a",
           ylab = "ln(D(a,w))", main ="D(a, w) by a and w")
tempRat <- aspectRatio() # capture device aspect ratio
for (ww in unique(uniqueAW$w)) {
    tempinds <- uniqueAW$w == ww
    lines(uniqueAW$a[tempinds], uniqueAW$logD[tempinds],
          col = adjustcolor(wcol(ww), 0.5))
}
for (ii in seq_along(uniqueAW$w)) { # plot glyphs, aspect ratio = 1
    if (dispLog[ii]) { # & log(uniqueAW$w[ii]) %in% c(-6, -3, -1, 0)) {
        densityGlyph(function(x) dbeta(x, uniqueAW$a[ii],
                                       bwa(uniqueAW$w[ii],
                                           uniqueAW$a[ii])),
                     xseq = seq(0, 1, length.out = 100),
                     x = uniqueAW$a[ii], y = log(uniqueKL[ii]),
                     devRat = tempRat, size = 0.22,
                     col = adjustcolor(wcol(uniqueAW$w[ii]), 0.5))
    }
}
legend(x = 1.04, y = 2, xpd = NA, bty = "n", cex = 0.8,
       legend = c(log(unique(params$w))[-length(unique(params$w))],
                  " 0"), # align -x, 0
       fill = adjustcolor(wcol(unique(params$w)), 0.6),
       title = "ln(w)", bg = "white")
dev.off()

## IMAGE :: a contour plot of the KL div by a and w
contour(x = seq(0, 1, length.out = 101), xlab = "a",
        y = seq(0, 1, length.out = 101), ylab = "w",
        z = log(outer(seq(0, 1, length.out = 101),
                      seq(0, 1, length.out = 101),
                      betaDiv)),
        xaxs = "i", yaxs = "i", main = "Beta KL divergence by a, w")
klDivMat <- log(outer(seq(0, 1, length.out = 101),
                      seq(0, 1, length.out = 101),
                      betaDiv)) # compute for image
image(klDivMat # plot as image
      breaks = exp(seq(-10, 10, by = 1)),
      col = colorRampPalette(c("floralwhite", "firebrick"))(20))

## IMAGE :: beta densities by a and w
x <- seq(0, 1, length.out = 401)
## add density lines
for (a in c(0.05, 0.5, 0.95)) {
    png(paste0("repBetaDens", a*100, ".png"), width = 2,
        height = 2, units = "in", res = 400) # start png
    narrowPlot(xgrid = seq(0, 1, by = 0.5),
               ygrid = seq(0, 2, by = 1),
               main = paste("a =", a), ylab = "Density",
               xlab = expression(x)) # plot area
    abline(h = c(1.5, 0.5), v = c(0.25, 0.75), lty = 1,
           col = adjustcolor("gray", alpha.f = 0.4)) # extra gridding
    for(w in exp(c(-6, -3, -log(2), 0))) { # density lines
        lines(x, y = dbeta(x, a, bwa(w, a)),
              lty = wlty(w), col = acol(a))
    }
    if (a == 0.05) { # legend for chosen case
        legend(x = "topright", lty = c(6, 3, 2, 1), title = "w",
               legend = c(1, "1/2", expression(e^-3),
                          expression(e^-6)),
               bty = "n", cex = 0.8, bg = NA)
    }
    dev.off()
}

## IMAGE :: UMP plotted by KL divergence, linked by a, coloured by w
png("klDivPowBya.png", width = 2.9, height = 2.6, units = "in", res = 400)
narrowPlot(ygrid = seq(0.2, 1, by = 0.2), xlab = "ln(D(a, w))",
           xgrid = seq(-6, 6, length.out = 5), ylab = "Power",
           ylim = c(0.05, 1), xlim = c(-7, 6),
           mars = c(2.1, 2.1, 1.1, 3.1),
           main = expression(paste(l[w], " power curves joined by a")))
for (a in unique(params$a)) {
    inds <- abs(params$a - a) < .Machine$double.eps &
        params$M == 2
    lines(log(klDivs[inds]), powershr[inds],
          col = adjustcolor("gray50", 0.6))
    points(log(klDivs[inds]), powershr[inds], cex = log(2, base = 3),
           col = adjustcolor(wcol(params$w[inds]), 0.8), pch = 20)
}
legend(x = 6.6, y = 1, xpd = NA, bty = "n", cex = 0.8, pch = 20,
       legend = c(log(unique(params$w))[-length(unique(params$w))],
                  " 0"),
       col = adjustcolor(wcol(unique(params$w)), 0.8), lty = 1,
       title = "ln(w)")
dev.off()

## IMAGE :: show the impact of M on power for a = 0.95
png("lwPowersByM.png", width = 3, height = 3, units = "in", res = 400)
narrowPlot(xgrid = seq(from = min(log(params$w)),
                       to = max(log(params$w)),
                       length.out = 4),
           ygrid = seq(0, 1, by = 0.25),
           main = expression(paste("Power of ", l[w], " by ", w,
                                   " and M")),
           ylab = "Power", xlab = "ln(w)")
## add density lines
a <- 0.95
for (m in unique(params$M)) {
    inds <- abs(params$a - a) < .Machine$double.eps &
        params$M == m & params$w %in% exp(seq(-6, 0, by = 1))
    lines(x = seq(-6, 0, by = 1), powershr[inds], lwd = Mlwd(m),
          col = "steelblue")
    points(x = seq(-6, 0, by = 1), powershr[inds], pch = 23,
           cex = 0.7, col = "steelblue")
}
## add a legend
legend(x = "bottomleft", legend = c(2, 5, 10, 20), lty = 1,
       lwd = c(1, 2, 3, 4), col = "steelblue", title = "M",
       bty = 'n', bg = NA, cex = 0.8)
dev.off()

## IMAGE :: plot showing the impact of a on power
png("lwPowersBya.png", width = 3, height = 3, units = "in", res = 400)
narrowPlot(xgrid = seq(from = min(log(params$w)),
                       to = max(log(params$w)),
                       length.out = 4),
           ygrid = seq(0, 1, by = 0.25),
           main = expression(paste("Power of ", l[w], " by ", w,
                                   " and ", a)),
           ylab = "Power", xlab = "ln(w)")
m <- 2 # set M value
for (a in c(0.05, 0.5, 0.95)) {
    inds <- abs(params$a - a) < .Machine$double.eps &
        params$M == m & params$w %in% exp(seq(-6, 0, by = 1))
    lines(x = seq(-6, 0, by = 1), powershr[inds], col = acol(a),
          lwd = Mlwd(m))
    points(x = seq(-6, 0, by = 1), powershr[inds], pch = apch(a),
           col = acol(a), cex = 0.7)
}
## add a legend
legend(x = "bottomleft", legend = c(0.05, 0.5, 0.95), lty = 1,
       pch = c(21, 24, 23), title = "a", bty = 'n', bg = NA,
       col = c("firebrick", "seagreen", "steelblue"),
       cex = 0.8)
dev.off()

## IMAGE :: misspeced final plot for a = 0.05, 0.5, 0.95
png("lwPowersMis.png", width = 3, height = 3, units = "in", res = 400)
narrowPlot(xgrid = seq(from = min(log(params$w)),
                       to = max(log(params$w)),
                       length.out = 4),
           ygrid = seq(0, 1, by = 0.25),
           main = expression(paste("Power of ", l[w], " by ", omega)),
           ylab = "Power",
           xlab = expression(paste("ln(", omega, ")")))
## add faint power curves
for (a in c(0.05, 0.50, 0.95)) {
    inds <- abs(params$a - a) < .Machine$double.eps &
        params$M == m & params$w %in% exp(seq(-6, 0, by = 1))
    #lines(x = log(params$w[inds]), y = powershr[inds],
    #      col = adjustcolor(acol(a), alpha.f = 0.5))
    tempster <- 1.96*sqrt(powershr[inds]*(1 - powershr[inds])/nsim)
    polygon(x = c(log(params$w[inds]), rev(log(params$w[inds]))),
            y = c(powershr[inds] + tempster,
                  rev(powershr[inds] - tempster)),
            col = adjustcolor(acol(a), alpha.f = 0.25), border = NA)
    for (ii in seq_len(ncol(powersmis))) {
        lines(x = log(params$w[inds]),
              y = powersmis[inds, ii],
              col = acol(a), lty = wlty(misSpec[ii]))
        points(x = seq(-6, 0, by = 1),
               y = powersmis[inds & params$w %in%
                            exp(seq(-6, 0, by = 1)), ii],
               pch = apch(a), col = acol(a), cex = 0.7)
    }
}
## add a legend
legend(x = "bottomleft", lty = wlty(misSpec), title = "ln(w)",
       legend = c("-6", "-3", "-ln(2)", "0"), bg = NA, bty = "n",
       cex = 0.8)
dev.off()

## IMAGE :: misspeced final plot with chi method lines added
png("chiComparedMis.png", width = 3, height = 3, units = "in",
    res = 400)
narrowPlot(xgrid = seq(from = min(log(params$w)),
                       to = max(log(params$w)),
                       length.out = 4),
           ygrid = seq(0, 1, by = 0.25),
           main = expression(paste("Power of ", l[w], " and ",
                                   g[chi], "(;2981)",
                                   " by ", omega)),
           ylab = "Power",
           xlab = expression(paste("ln(", omega, ")")))
## add faint power curves
for (a in c(0.05, 0.50, 0.95)) {
    inds <- abs(params$a - a) < .Machine$double.eps &
        params$M == m & params$w %in% exp(seq(-6, 0, by = 1))
    lines(x = log(params$w[inds]), y = powershr[inds],
          col = adjustcolor(acol(a), alpha.f = 0.5))
    tempster <- 1.96*sqrt(powershr[inds]*(1 - powershr[inds])/nsim)
    polygon(x = c(log(params$w[inds]), rev(log(params$w[inds]))),
            y = c(powershr[inds] + tempster,
                  rev(powershr[inds] - tempster)),
            col = adjustcolor(acol(a), alpha.f = 0.25), border = NA)
    for (ii in seq_len(ncol(powersmis))) {
        lines(x = log(params$w[inds]),
              y = powersmis[inds, ii],
              col = acol(a), lty = wlty(misSpec[ii]))
        points(x = seq(-6, 0, by = 1),
               y = powersmis[inds & params$w %in%
                            exp(seq(-6, 0, by = 1)), ii],
               pch = apch(a), col = acol(a), cex = 0.7)
        lines(x = log(params$w[inds]),
              y = powerschi[inds, length(kseq)])
    }
}
## add a legend
legend(x = "bottomleft", lty = wlty(misSpec), title = "ln(w)",
       legend = c("-6", "-3", "-ln(2)", "0"), bg = NA, bty = "n",
       cex = 0.8)
dev.off()

######################################################################


## SIMULATION ########################################################
## compute simulated powers when H3 is true under correct
## specification, misspecification, and using the chi method
## LOAD :: simulated data for variable M or else run simulation
readData <- TRUE
datafile <- "varMPowers_fullgrid.Rds"
## datafile <- "varMPowers.Rds"
if (readData) { # load specified file
    tempDat <- readRDS(datafile)
    kseq <- tempDat$ks # chi k values in simulation
    params <- tempDat$pars
    misSpec <- tempDat$misw
    powershr <- tempDat$tw
    powerschi <- tempDat$chi
    powersmis <- tempDat$mis
    rm(tempDat)
} else {
    if (grepl("fullgrid", datafile)) {
        source("simulateVarM_fullgrid.R")
    } else {
        source("simulateVarM.R")
    }
}

######################################################################


## PLOTS CURRENTLY IN THE PAPER ######################################
## some plotting packages
library(ggplot2)
library(grid)

## SETUP :: packages and some palettes
propPal <- colorRampPalette(c("white", "indianred4"))(100)
powPal <- colorRampPalette(c("white", "dodgerblue4"))(100)
difPal <-  colorRampPalette(c("steelblue", "floralwhite",
                              "firebrick"))(100)

## COMPUTE :: the KL divergence for each beta setting
klDivs <- betaDiv(params$a, params$w)

## DEFINE :: a helper that makes nice labels
ggLabs <- function(breaks){
    paste0(c("[", rep("(", length(breaks)-2)), # bottom brackets
           breaks[-length(breaks)], # break limits
           ", ", breaks[-1], "]") # top brackets
}

## SETUP :: some graphical parameters
nbr <- 10 # number of breaks for simple power
powBreaks <- c(seq(0, 1 - 1/nbr, length.out = nbr), 1.01)
powLabs <- ggLabs(round(powBreaks, 1)) # nice labels
powPal <- colorRampPalette(c("floralwhite", "firebrick")) # palette
## parameters for the difference plots
difBreaks <- c(-1.01, seq(-0.95, 0.95, by = 0.1), 1.01)
difLabs <- ggLabs(round(difBreaks, 1))
difPal <- colorRampPalette(c("steelblue", "floralwhite",
                             "firebrick")) # diverging palette
## set a new theme based on the "lined raw" template
myTheme <- theme_linedraw()
myTheme$text$size <- 8 # text size
myTheme$plot.title$size <- rel(1.1) # scale down title
myTheme$strip.text$colour <- "black" # facet title colour
myTheme$strip.text$size <- 8 # facet title size
myTheme$strip.background$fill <- "white"
myTheme$panel.grid$colour <- adjustcolor("gray", 0.4)
myTheme$plot.title$hjust <- 0.5 # centre title

## IMAGE :: facetted power contours of the UMP using ggplot
## change ump = powershr -> fis = powersmis[,4]
##        umpPower -> fisPower
##        z = ump -> z = fis
##        l[w] -> l[1]
##        D(a,w) -> D(a, omega)
##        ln(w) -> ln( omega )
## and set width = height = 3 to obtain analogous plot of fisher's
## method
inds <- log(params$w) %in% c(-6, -3, 0) & -5 <= log(klDivs) &
    5 >= log(klDivs)
powerContourBase <- ggplot(data = data.frame(ump = powershr,
                                      lnkl = log(klDivs),
                                      rho = params$m1/10,
                                      a = params$a,
                                      w = params$w)[inds,], # temporary df
                           aes(lnkl, rho, z = ump)) # base ggplot
png("umpPower.png", width = 6, height = 2.5, units = "in", res = 400)
powGrob <- ggplotGrob(powerContourBase + # save as a grid grob
                      xlab(expression(paste("ln(D(a, ", "w", "))"))) +
                      ylab(expression(rho)) +
                      ggtitle(expression(paste("Power of ", g[w],
                                               " by ln(D(a,", "w",
                                               ")) and ", rho,
                                               " facetted by ln(",
                                               "w", ")"))) +
                      geom_contour_filled(breaks = powBreaks) +
                      facet_wrap(~log(w)) + # facet by w
                      scale_fill_manual(values = powPal(nbr+1),
                                        name = "Power",
                                        labels = powLabs,
                                        drop = FALSE,
                                        guide = "none") +
                      myTheme) # set theme
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

## COMPUTE :: power differences from fis to others
fisMin <- apply(cbind(ump = powershr, chi = powerschi,
                      mis = powersmis[,1:3]), 2,
                function(col) powersmis[,4] - col) # all in one matrix
colnames(fisMin) <- c("ump", paste0("chi", round(log(kseq), 1)),
                      paste0("mis", round(log(misSpec[1:3]), 1)))
difTitles <- c("ump" = quote(g[omega]), # nice names for plotting
               "minchi" = quote(paste(g[chi], " (", kappa, " = 0.0003)")),
               "smchi" = quote(paste(g[chi], " (", kappa, " = 0.018)")),
               "cvchi" = quote(paste(g[chi], " (", kappa, " = 1)")),
               "fischi" = quote(paste(g[chi], " (", kappa, " = 2)")),
               "modchi" = quote(paste(g[chi], " (", kappa, " = 55)")),
               "bigchi" = quote(paste(g[chi], " (", kappa, " = 2981)")),
               "le-6" = quote(g[0.002]),
               "le-3" = quote(g[0.05]),
               "l05" = quote(g[0.5]))

## IMAGE :: facetting power difference contours by w
inds <- log(params$w) %in% c(-6, -3, 0) & -5 <= log(klDivs) &
    5 >= log(klDivs)
for (ii in seq_len(ncol(fisMin))) {
    # temporary base for the particular case
    diffContourBase <- ggplot(data = data.frame(lnkl = log(klDivs),
                                                rho = params$m1/10,
                                                a = params$a,
                                                w = params$w,
                                                dif = fisMin[,ii])[inds,],
                              mapping = aes(lnkl, rho, z = dif))
    diffGrob <- ggplotGrob( # save in a grid grob as before
        diffContourBase +
        xlab(expression(paste("ln(D(a, ", omega, "))"))) +
        ylab(expression(rho)) +
        ggtitle(bquote(paste("Power of ", g[1], " minus power of ",
                       .(difTitles[[ii]])))) + # dynamic title
        geom_contour_filled(breaks = difBreaks) +
        facet_wrap(~log(w)) + # facet by w
        scale_fill_manual(values = difPal(length(difBreaks)-1),
                          labels = difLabs, name = "Difference",
                          drop = FALSE, guide = "none") +
        myTheme) # set theme
    png(paste0("diff", names(difTitles)[ii], ".png"), width = 6,
        height = 2.5, units = "in", res = 400) # open named png
    grid.newpage()
    pushViewport(viewport(x = unit(0.45, "npc"), # facet plot vp
                          width = unit(0.9, "npc")))
    grid.draw(diffGrob)
    popViewport()
    pushViewport(viewport(x = unit(0.93, "npc"), # legend vp
                          height = unit(0.4, "npc"),
                          width = unit(0.05, "npc")))
    grid.text("", y = unit(1, "npc") + unit(1, "lines"), # title
              x = 0.7, gp = gpar(fontsize = 8))
    pushViewport(viewport(x = unit(0.25, "npc"), # scale vp
                          height = unit(1, "npc"),
                          width = unit(0.5, "npc")))
    grid.rect(gp = gpar(bg = NA)) # scale border
    grid.raster(matrix(rev(difPal(50)), ncol = 1), # scale raster
                width = unit(1, "npc"),
                height = unit(1, "npc"))
    grid.segments(x0 = unit(rep(1, 5), "npc"), # scale ticks
                  y0 = unit(seq(0, 1, by = 0.25), "npc"),
                  x1 = unit(rep(1, 5), "npc") +
                      unit(rep(0.25, 5), "lines"),
                  y1 = unit(seq(0, 1, by = 0.25), "npc"),
                  gp = gpar(lwd = 0.5))
    grid.text(label = c(-1, -0.5, 0, 0.5, 1), hjust = 0, # tick labels
              y = unit(seq(0, 1, by = 0.25), "npc"),
              x = unit(rep(1, 5), "npc") + unit(rep(0.5, 5), "lines"),
              gp = gpar(fontsize = 8))
   dev.off()
}

######################################################################


## OLD PLOTS FOR REFERENCE, POTENTIAL USE ############################
## IMAGE :: the power by point size on klDiv by rho
plot(log(klDivs), params$m1/10, cex = powershr*2 + 1, pch = 20,
     col = adjustcolor(wcol(params$w), 0.4))

## IMAGE :: plot the power by klDiv, size by M1
plot(log(klDivs), powershr, cex = params$m1/5 + 0.5, pch = 20,
     col = adjustcolor(wcol(params$w), 0.4))
## alt version
narrowPlot(xgrid = seq(-6, 6, length.out = 5),
           ygrid = seq(0, 1, by = 0.5), xlim = c(-7, 6))
for (m1 in unique(params$m1)) {
    for (w in unique(params$w)) {
        tempinds <- abs(params$w - w) <= .Machine$double.eps &
            params$m1 == m1
        lines(log(klDivs[tempinds]), powershr[tempinds],
              col = adjustcolor(wcol(w), 0.8))
        #points(log(klDivs[tempinds]), logit(powershr[tempinds]), pch = 19,
        #       col = adjustcolor(wcol(w), 0.4), cex = m1/5 + 0.5)
    }
}

## COMPUTE :: interpolate divergences
## library(mgcv) # for gams with smoothing

## COMPUTE :: the most powerful test for different w, omega, settings
## when not all tests are non-null and w may be misspecified
allpowers <- cbind(true = powershr, mis = powersmis) # table of powers
bestPow <- apply(allpowers, 1, max) # identify the maximum power
bPowByPar <- tapply(bestPow, # aggregate by parameter combination
                    list(w = log(params$w), Mrho = params$m1),
                    mean)
nearMax <- sweep(allpowers, 1, bestPow, # binomial test of whether
                 function(p1, p2) {     # a row matches this max
                     p <- (p1 + p2)/2
                     abs(p1 - p2)/sqrt(p*(1-p)*2/nsim) <= 1.96
                 })
nearMax[is.na(nearMax)] <- TRUE # cases where both powers = 1
isMax <- sweep(allpowers, 1, bestPow, `==`) # simpler: equal to max
table(omega = isMax[,1], zero = isMax[,5], half = isMax[,4])
## for what parameter settings does knowledge beat w = 1?
omgaBts1 <- !nearMax[,5] & nearMax[,1] & params$w != 1
cbind(params[omgaBts1,], allpowers[omgaBts1,c(1,5)])
## mabybe plot some of these
tapply(allpowers[,5], params$m1, mean)
tapply(allpowers[,4], params$m1, mean)
tapply(allpowers[,1], params$m1, mean)
## or the differences by m1 and omega
propFisBest <- tapply(nearMax[,5],
                      list(w = log(params$w), Mrho = params$m1),
                      mean)
proplwBest <-  tapply(nearMax[,1],
                      list(w = log(params$w), Mrho = params$m1),
                      mean)

## IMAGE :: a sliced plot of the maximum power
bestPower <- array(bestPow,
                   dim = lapply(params,
                                function(a) length(unique(a))),
                   dimnames = lapply(params, unique))
par(mfrow = c(4,2))
for (a in unique(params$a)) {
    par(mar = c(2.1, 2.1, 1.1, 1.1)) # set narrow margins
    narrowImage(t(bestPower[as.character(a), , -1]),
                col = powPal, ynames = -6:0, main = paste("a =", a),
                ylab = "ln(w)", xlab = expression(paste("M",rho)),
                breaks = seq(0, 1, length.out = 16))
}

## IMAGE :: sliced plot of which method is best
bestSet <- apply(allpowers, 1, which.max)
bestSet[bestSet == 1 & params$w == 1] <- 5
bestSet[bestSet == 1 & log(params$w) == -6] <- 2
bestSet[bestSet == 1 & log(params$w) == -3] <- 3
par(mfrow = c(4,2))
for (a in unique(params$a)) {
    narrowPlot(xgrid = seq(0, 10, by = 2), ygrid = seq(-6, 0, by = 1),
               main = paste("a =", a))
    text(x = params$m1[params$a == a],
         y = log(params$w)[params$a == a],
         labels = bestSet[params$a == a])
}

## IMAGE :: sliced plot of where Fisher's method is equal to the max
par(mfrow = c(4,2))
for (a in unique(params$a)) {
    par(mar = c(2.1, 2.1, 1.1, 1.1))
    narrowImage(t(matrix(nearMax[params$a == a, 5], nrow = 7))[-1,],
                ynames = log(unique(params$w)),
                xnames = (unique(params$m1))[-1],
                col = c("white", "indianred4"),
                breaks = c(-0.1, 0.5, 1.1),
                main = paste("a =", a))
}

## IMAGE :: plot the proportion of cases where Fisher's method has
## the same power as the best as determined by a binomial
## difference test
png(paste0("propFishBest.png"), width = 3, height = 2,
    units = "in", res = 400)
par(mar = c(2.1, 2.1, 1.1, 2.1)) # set narrow margins
narrowImage(t(propFisBest[,-1]), xlab = expression(paste("M",rho)),
            ylab = expression(paste("ln(", omega, ")")), las = 1,
            col = propPal,
            breaks = seq(0, 1, length.out = 16))
bds <- par()$usr
rasterImage(as.raster(matrix(rev(propPal), ncol = 1)),
            bds[2] + 0.05, 0.2, bds[2] + 0.1, 0.8, xpd = NA)
rect(bds[2] + 0.05, 0.2, bds[2] + 0.1, 0.8, xpd = NA)
text(x = rep(bds[2] + 0.1, 3), y = c(0.2, 0.5, 0.8),
     labels = c(0, "0.5", 1), cex = 0.8, xpd = NA, adj = c(-0.15,0.5))
dev.off()

## IMAGE :: plot the proportion of cases where lw is best
png(paste0("proplwBest.png"), width = 3, height = 2,
    units = "in", res = 400)
par(mar = c(2.1, 2.1, 1.1, 2.1)) # set narrow margins
narrowImage(t(proplwBest[,-1]), xlab = expression(paste("M",rho)),
            ylab = expression(paste("ln(", omega, ")")), las = 1,
            col = propPal,
            breaks = seq(0, 1, length.out = 16))
dev.off()

## COMPUTE :: the z-scores of the chi method against the most powerful
## lw, use these to see when it is best by k
chiZ <- sweep(powerschi, 1, bestPow,
              function(p1, p2) {
                  p <- (p1 + p2)/2
                  (p1 - p2)/sqrt(p*(1-p)*2/nsim)
              })
chiZ[is.nan(chiZ)] <- 0
chiBest <- chiZ >= -1.96
chiPropBest <- apply(chiBest, 2,
                     function(col) {
                         tapply(col,
                                list(w = log(params$w),
                                     Mrho = params$m1),
                                mean)
                     }, simplify = FALSE)

## EXPLORE :: plot these side-by-side
par(mfrow = c(3,2), mar = c(2.1, 2.1, 1.1, 2.1))
for (ii in seq_along(chiPropBest)) {
    narrowImage(t(chiPropBest[[ii]][,-1]),
                main = bquote(kappa == .(kseq[ii])),
                col = propPal, breaks = seq(0, 1, length.out = 16),
                xlab = expression(paste("M", rho)),
                ylab = expression(paste("ln(", omega, ")")))
}

## IMAGE :: output heatmap for k = e^-8
for (ii in match(c(exp(-8), 2, exp(8)), kseq)) {
    if (ii == length(kseq)) {
        png(paste0("propChiBest", round(log(kseq[ii]), 2), ".png"),
        width = 2.2, height = 1.5, units = "in", res = 400)
        par(mar = c(2.1, 2.1, 1.1, 2.1))
    } else  {
        png(paste0("propChiBest", round(log(kseq[ii]), 2), ".png"),
            width = 1.8, height = 1.5, units = "in", res = 400)
        par(mar = c(2.1, 2.1, 1.1, 0.1))
    }
    narrowImage(t(chiPropBest[[ii]][,-1]),
                xlab = expression(paste("M",rho)),
                ylab = expression(paste("ln(", omega, ")")),
                col = propPal, las = 1,
                main = bquote(kappa == .(kseq[ii])),
                breaks = seq(0, 1, length.out = 16))
    if (ii == length(kseq)) {
        bds <- par()$usr
        rasterImage(as.raster(matrix(rev(propPal), ncol = 1)),
                    bds[2] + 0.05, 0.2, bds[2] + 0.1, 0.8, xpd = NA)
        rect(bds[2] + 0.05, 0.2, bds[2] + 0.1, 0.8, xpd = NA)
        text(x = rep(bds[2] + 0.1, 3), y = c(0.2, 0.5, 0.8),
             labels = c(0, "0.5", 1), cex = 0.8, xpd = NA,
             adj = c(-0.15,0.5))
    }
    dev.off()
}

## IMAGE :: split by a and plot each separately
par(mfrow = c(4,2), mar = c(2.1, 2.1, 1.1, 2.1))
ki <- 4
for (a in unique(params$a)) {
    narrowImage(t(matrix(chiBest[params$a == a, ki], nrow = 7))[-1,],
                ynames = log(unique(params$w)),
                xnames = (unique(params$m1))[-1],
                col = c("white", "indianred4"),
                breaks = c(-0.1, 0.5, 1.1),
                main = paste("a =", a))
}

## IMAGE ::  a plot of the power of l_w by rho
for (w in c(-6, -3, 0)) { # split plots by w
    png(paste0("powerByM1w", w, ".png"), width = 2,
        height = 2, units = "in", res = 400)
    narrowPlot(xgrid = seq(0, 10, by = 5)/10,
               ygrid = seq(0, 1, by = 0.5),
               main = paste("ln(w) = ", w),
               xlab = expression(rho), ylab = "Power")
    abline(v = c(0.25, 0.75), h = c(0.25, 0.75), lty = 1,
           col = adjustcolor("gray", alpha.f = 0.4))
    for (a in c(0.05, 0.5, 0.95)) { # split plots by a
        inds <- abs(params$a - a) < .Machine$double.eps &
            log(params$w) == w
        lines(params$m1[inds]/10, powershr[inds], col = acol(a),
              lty = wlty(exp(w)))
        ptinds <- abs(params$a - a) < .Machine$double.eps &
            log(params$w) == w & params$m1 %in% c(2,4,6,8)
        points(params$m1[ptinds]/10, powershr[ptinds],
               pch = apch(a), col = acol(a))
    }
    if (w == -6) { # add a legend
       legend(x = "bottomright", legend = c(0.05, 0.5, 0.95), lty = 1,
              pch = c(21, 24, 23), title = "a", bty = 'n', bg = NA,
              col = c("firebrick", "seagreen", "steelblue"),
              cex = 0.8)
    }
    dev.off()
}

## IMAGE :: add lines for the chi squared method, compare to the
## best of the lw lines
for (w in c(-6, -3, 0)) { # split plots by w
    png(paste0("powerChiByM1w", w, ".png"),
        width = 2, height = 2, units = "in", res = 400)
    narrowPlot(xgrid = seq(0, 10, by = 5)/10,
               ygrid = seq(0, 1, by = 0.5),
               main = paste("ln(w) = ", w),
               xlab = expression(rho), ylab = "Power")
    abline(v = c(0.25, 0.75), h = c(0.25, 0.75), lty = 1,
           col = adjustcolor("gray", alpha.f = 0.4))
    for (a in c(0.05, 0.5, 0.95)) { # split plots by a
        inds <- abs(params$a - a) < .Machine$double.eps &
            log(params$w) == w
        lines(params$m1[inds]/10, powershr[inds], col = acol(a),
              lty = wlty(exp(w)))
        ## chi lines
        lines(params$m1[inds]/10, powerschi[inds, 4])
        ptinds <- abs(params$a - a) < .Machine$double.eps &
            log(params$w) == w & params$m1 %in% c(2,4,6,8)
        points(params$m1[ptinds]/10, powershr[ptinds],
               pch = apch(a), col = acol(a))
        ## chi points
        points(params$m1[ptinds]/10, powerschi[ptinds, 4],
               pch = apch(a))
    }
    if (w == -6) { # add a legend
        legend(x = "bottomright", legend = c(0.05, 0.5, 0.95),
               lty = 1, pch = c(21, 24, 23), title = "a",
               bty = 'n', bg = NA,
               col = c("firebrick", "seagreen", "steelblue"),
               cex = 0.8)
    }
    dev.off()
}

## IMAGE :: mis-specified w under H3
for (omga in c(-6, -3, 0)) { # split plots by omega
    for (a in c(0.05, 0.5, 0.95)) { # and a
        png(paste0("misByM1w", omga, "a", a*100, ".png"), width = 2,
            height = 2, units = "in", res = 400)
        narrowPlot(xgrid = seq(0, 10, by = 5)/10,
                   ygrid = seq(0, 1, by = 0.5),
                   main = bquote(paste("ln(", omega, ") = ", .(omga))),
                   xlab = expression(rho), ylab = "Power")
        abline(v = c(0.25, 0.75), h = c(0.25, 0.75), lty = 1,
               col = adjustcolor("gray", alpha.f = 0.4))
        inds <- abs(params$a - a) < .Machine$double.eps &
            log(params$w) == omga
        for (w in log(misSpec)) { # misspecified w
            tempmiscol <- which(log(misSpec) == w)
            lines(params$m1[inds]/10, powersmis[inds, tempmiscol],
                  col = acol(a), lty = wlty(exp(w)))
            ptinds <- abs(params$a - a) < .Machine$double.eps &
                log(params$w) == omga & params$m1 %in% c(2,4,6,8)
            points(params$m1[ptinds]/10, powersmis[ptinds, tempmiscol],
                   pch = apch(a), col = acol(a))
        }
        if (omga == -6 & a %in% c(0.5, 0.95)) { # add a legend
            legend(x = "bottomright", lty = wlty(misSpec),
                   title = "ln(w)", bg = NA, bty = "n",
                   legend = c("-6", "-3", "-ln(2)", "0"),
                   cex = 0.8)
        }
        dev.off()
    }
}

######################################################################


## CENTRAL/MARGINAL REJECTION ########################################
## COMPUTE :: the central/marginal rejection levels of l_w for three
## w: ln(w) = -6, -3, 0
Mseq <- c(2, 5, 10, 20)
wseq <- exp(c(-6, -3, 0))
params <- expand.grid(w = wseq, M = Mseq)
lwPcPool <- mapply(pcHR, w = params$w,  M = params$M)
lwPrPool <- mapply(prHR, w = params$w,  M = params$M)
## determine central rejection
lwPc <- sapply(lwPcPool,
               function(f) {
                   uniroot(function(x) f(x) - 0.05,
                           interval = c(0, 0.999))$root
               })
## and marginal rejection
lwPr <- sapply(lwPrPool,
               function(f) {
                   uniroot(function(x) f(x) - 0.05,
                           interval = c(0, 0.999))$root
               })
## compute g
lwg <- (lwPc - lwPr)/lwPc

## IMAGE :: the centrality quotient for the chi method
Mseq <- c(2, 5, 20, 100, 500, 2000, 10000) # M sequence
a <- 0.05 # alpha value
k <- exp(seq(-16, 16, by = 0.2)) # k on log scale
png("cgChiPool.png", width = 4, height = 4, units = "in",
    res = 400)
narrowPlot(xgrid = seq(-15, 15, by = 5), ygrid = seq(0, 1, by = 0.2),
           main = expression(paste("c(",g[chi],") by M and ln(",
                                   kappa, ")")),
           xlim = c(-16, 16),
           xlab = expression(paste("ln(", kappa, ")")),
           ylab = "Centrality quotient")
for (m in Mseq) {
    tempQuot <- centQuot(k, M = m, alpha = a)
    lines(log(k), tempQuot, #lwd = templwd,
          col = adjustcolor(chiCol, 1))
    text(log(k)[which.min(abs(tempQuot - 0.5))] - 0.6,
         y = 0.5, labels = m, srt = 80, cex = 0.7)
}
dev.off()

## TABLE :: the choice of kappa to given a particular centrality
## quotient given a particular M
centSeq <- seq(0.1, 0.9, by = 0.1)
centTab <- log(sapply(Mseq,
                      function(m) sapply(centSeq, chiKappa,
                                         M = m,
                                         interval = c(0, 1e4),
                                         tol = .Machine$double.eps^0.75)))
paste(apply(cbind(c("", Mseq), t(cbind(centSeq, round(centTab, 1)))),
            1, paste, collapse = " & "), collapse = "\\")

## IMAGE :: the central rejection level for the chi method by M, kappa
png("pcChiPool.png", width = 2, height = 2, units = "in",
    res = 400)
narrowPlot(xgrid = c(-10, 0, 10), ygrid = seq(0, 0.4, by = 0.2),
           ylim = c(0, 0.45),
           main = expression(paste(p[c], " by M and ln(", kappa, ")")),
           ylab = expression(p[c]),
           xlab = expression(paste("ln(", kappa, ")")))
abline(h = c(0.1, 0.3), v = c(-5, 5),
       col = adjustcolor("gray", alpha.f = 0.4))
points(x = rep(log(min(k)), length(Mseq)), y = 1 - (1 - a)^(1/Mseq),
       col = "firebrick", pch = 7)
points(x = rep(log(max(k)), length(Mseq)),
       y = 1 - pnorm(qnorm(1 - a), sd = sqrt(Mseq)),
       col = "steelblue", pch = 9)
for (m in Mseq) {
        templwd <- if (m == 2) {
                   1 } else if (m == 5) {
                         2 } else if (m == 10) {
                               3 } else if (m == 20) {
                                     4 } else 5
        lines(log(k), chiPc(k, alpha = a, M = m), lwd = templwd,
              col = adjustcolor(chiCol, 0.8))
}
dev.off()

## IMAGE :: pr by kappa, M for the chi method
png("prChiPool.png", width = 2, height = 2, units = "in",
    res = 400)
narrowPlot(xgrid = c(-10, 0, 10), ygrid = seq(0, 0.03, 0.01),
           ylim = c(0, 0.03),
           main = expression(paste(p[r], " by M and ln(", kappa, ")")),
           ylab = expression(p[r]),
           xlab = expression(paste("ln(", kappa, ")")))
abline(v = c(-5, 5), col = adjustcolor("gray", alpha.f = 0.4))
points(x = rep(log(min(k)), length(Mseq)), y = 1 - (1 - a)^(1/Mseq),
       col = "firebrick", pch = 7)
points(x = rep(log(max(k)), length(Mseq)), y = rep(0, length(Mseq)),
       col = "steelblue", pch = 9)
for (m in Mseq) {
    templwd <- if (m == 2) {
                   1 } else if (m == 5) {
                         2 } else if (m == 10) {
                               3 } else if (m == 20) {
                                     4 } else 5
    lines(log(k), chiPr(k, alpha = a, M = m), lwd = templwd,
          col = adjustcolor(chiCol, 0.8))
}
legend(x = "topright", legend = c(2, 5, 10, 20, 50), lty = 1,
       lwd = c(1, 2, 3, 4, 5), col = adjustcolor(chiCol, 0.8),
       title = "M", bty = 'n', bg = NA, cex = 0.8)
dev.off()

## SIMULATE :: random uniform data to see how chi, tippett, and
## normal pooling functions compare
set.seed(1209149)
M <- 5
nsim <- 1000
randDat <- matrix(runif(M*nsim), ncol = M)
## k settings
kseq <- exp(seq(-8, 8, length.out = 200))
kseq[100] <- 1
kcurves <- apply(randDat, 1,
                 function(col) sapply(kseq, poolChi, p = col))
## compute the p-values for other methods
tipPs <- apply(randDat, 1, poolTip)
normPs <- apply(randDat, 1, poolNorm)

## IMAGE :: compare given chi method to normal, tippett methods
kind <- 30
chikPs <- kcurves[kind,] # select k
## plot in one display
png(paste0("chiTipSto", kind, ".png"), width = 2.5, height = 2.5,
    units = "in", res = 400)
par(mar = c(2.1, 2.1, 1.1, 1.1))
plot(c(0.48*tipPs, 0.48*normPs + 0.52), c(chikPs, chikPs),
     xaxs = 'i', yaxs = 'i', xlim = c(-0.02, 1.02),
     ylim = c(-0.02, 1.02), xaxt = 'n', xlab = "", yaxt = 'n',
     ylab = "", main = "", pch = 20,
     col = adjustcolor(chiCol, 0.4))
mtext(expression(paste(g[chi], "(", bold(p)[i], ")")),
      side = 2, line = 1, cex = 0.8)
tempk <- round(kseq[kind], 4)
mtext(bquote(paste(kappa, " = ", .(tempk))), side = 3, line = 0,
      cex = 0.8)
## grids
abline(v = 0.5)
abline(h = seq(0, 1, by = 0.25), col = adjustcolor("gray", 0.4),
       lty = 1, v = c(seq(0, 0.48, length.out = 5),
                      seq(0.52, 1, length.out = 5)))
## custom y axis
mtext(side = 2, at = seq(0, 1, length.out = 5),
      text = "|", line = 0, cex = 0.5, padj = 1)
mtext(text = seq(0, 1, length.out = 5), side = 2, cex = 0.8,
      at = seq(0, 1, length.out = 5), padj = -0.5)
## custom x axes
mtext(side = 1, at = seq(0, 0.48, length.out = 5),
      text = "|", line = 0, cex = 0.5, padj = -2)
mtext(text = seq(0, 1, length.out = 3), side = 1, cex = 0.8,
      at = seq(0, 0.48, length.out = 3))
mtext(text = expression(g[Tip](bold(p)[i])), at = 0.25, side = 1,
      line = 1, padj = 0, cex = 0.8)
mtext(side = 1, at = seq(0.52, 1, length.out = 5),
      text = "|", line = 0, cex = 0.5, padj = -2)
mtext(text = seq(0, 1, length.out = 3), side = 1, cex = 0.8,
      at = seq(0.52, 1, length.out = 3))
mtext(text = expression(g[Sto](bold(p)[i])), at = 0.75, side = 1,
      line = 1, padj = 0, cex = 0.8)
dev.off()

######################################################################


## OLD PLOTS #########################################################
## compute correlations with other method for all k
corrs <- apply(kcurves, 1,
               function(row) apply(otherPs, 1, cor, y = row,
                                   method = "spearman"))
## plot these
krng <- 30:200
png("chiStoTipCorrk.png", width = 5, height = 5, units = "in",
    res = 400)
plot(x = log(kseq)[krng], type = 'n',
     y = corrs["stouffer", krng], xlab = expression(paste(log[e],k)),
     ylab = expression(paste("Correlation with ", ~g[chi], "(k)")))
abline(h = seq(0.65, 1, by = 0.05), lty = 2, col = 'gray50',
       v = seq(-6, 8, by = 2))
for (mth in rownames(corrs)) {
    lines(log(kseq)[krng], y = corrs[mth,krng], col = pal[mth])
}
dev.off()

## extension idea: plot k curves, coloured by minimum value
minVals <- apply(randDat, 1, min)
minCol <- as.numeric(cut(minVals, breaks = quantile(minVals),
                         include.lowest = TRUE))
pal <- hcl.colors(4)
plot(x = log(kseq), y = rep(range(kcurves), length.out = length(kseq)),
     type = 'n', xlab = expression(logk), ylab = "p")
for (ii in 1:ncol(kcurves)) {
    lines(log(kseq), kcurves[, ii],
          col = adjustcolor(pal[minCol[ii]], alpha.f = 0.4))
}

## which data changed values the most along the range -5, 5
changes <- apply(kcurves, 2, function(col) sum(diff(col)))
changOrd <- order(changes, decreasing = TRUE)

## the image showing the range of rejection regions
x <- exp(seq(-7, 0, by = 0.01)) # x values (log scale)
kseq <- exp(seq(-8, 8, by = 1)) # k values (log scale)
axislogk <- seq(-7, 0, by = 1) # for nice labels
chipal <- colorRampPalette( # define a palette
    brewer.pal(3, "Dark2")[c(2,3)] )(length(kseq)+2)
png("ChiRejectSeq.png", width = 5, height = 5, units = "in",
    res = 400) # create image file
plot(x = x, y = x, type = 'n', xlab = expression(paste(log[e],p[1])),
     xlim = c(0.001, 1.2), ylim = c(0.001, 1.2),
     xaxs = "i", yaxs = "i",
     ylab = expression(paste(log[e],p[2])), log = "xy",
     main = expression(paste("Rejection boundaries for ", g[chi],
                             " by k")),
     xaxt = 'n', yaxt = 'n')
abline(h = exp(axislogk), v = exp(axislogk), col = "gray50", lty = 2)
axis(side = 1, at = exp(axislogk), labels = axislogk)
axis(side = 2, at = exp(axislogk), labels = axislogk)
for (k in seq_along(kseq)) lines(x, chiCrv(x, k = kseq[k]),
                                 col = chipal[k+1])
lines(c(0, rep(1-0.95^0.5, 2), 1), c(1, 1, rep(1 - 0.95^0.5, 2)),
      lwd = 2, lty = 2)
text(1-0.95^0.5, 1-0.95^0.5, expression(g[Tip]), adj = c(1,1))
lines(x, normCrv(x), lwd = 2, lty = 2)
stoInd <- which.min(abs(x - normCrv(x)))
text(x[stoInd], normCrv(x[stoInd]), expression(g[Sto]), adj = c(0,0))
rasterImage(as.raster(matrix(rev(chipal), ncol = 1)), 0.002, 0.002,
            0.003, 0.02)
text(x = 0.0025, y = 0.022, expression(paste(log[e],k)),
     adj = c(0.5,0))
ticks <- exp(seq(log(0.002), log(0.02), length.out = 5))
for (ii in seq_along(ticks)) {
    lines(c(0.003, 0.0032), rep(ticks[ii], 2))
    text(x = 0.0035, y = ticks[ii], log(kseq[4*ii-3]))
}
dev.off()

## other methods in the same plot region
source("~/Research/Core/experiments/3. Pooled pvals/poolFunctions.R")
png("ChiOtherMethods.png", width = 5, height = 5, units = "in",
    res = 400)
plot(x = x, y = x, type = 'n', xlab = expression(paste(log[e],p[1])),
     xlim = c(0.001, 1.2), ylim = c(0.001, 1.2),
     xaxs = "i", yaxs = "i",
     ylab = expression(paste(log[e],p[2])), log = "xy",
     main = expression("Rejection boundaries for previous methods"),
     xaxt = 'n', yaxt = 'n')
abline(h = exp(axislogk), v = exp(axislogk), col = "gray50", lty = 2)
axis(side = 1, at = exp(axislogk), labels = axislogk)
axis(side = 2, at = exp(axislogk), labels = axislogk)
for (mth in names(methsNiceName)) lines(x, methsNiceName[[mth]](x),
                                        lty = 1, col = pal[mth])
lines(c(0, rep(1-0.95^0.5, 2), 1), c(1, 1, rep(1 - 0.95^0.5, 2)),
      lwd = 2, lty = 2)
text(1-0.95^0.5, 1-0.95^0.5, expression(g[Tip]), adj = c(1,1))
lines(x, Stocrv(x), lwd = 2, lty = 2)
text(x[stoInd], Stocrv(x[stoInd]), expression(g[Sto]), adj = c(0,0))
legend(x = "bottomleft", title = "Method", lty = 1,
       legend = c(expression(g[Bon]), expression(g[Tip]),
                  expression(g[Fis]), expression(g[Sto]),
                  expression(g[CV]), expression(g[Wil]),
                  expression(g[Edg]), expression(g[MG]),
                  expression(g[HMP]),
                  expression(paste(g[HMP], " (exact)"))),
       col = pal[names(methsNiceName)], bg = "white", cex = 0.8)
dev.off()

## some power curve plots: generate curves
powerTip <- powerCrv(crv = tipCrv, crvArgs = list(),
                     medianScale = TRUE)
power0.005 <- powerCrv(crvArgs = list(k = 0.005), medianScale = TRUE)
power0.05 <- powerCrv(crvArgs = list(k = 0.5), medianScale = TRUE)
power0.5 <- powerCrv(crvArgs = list(k = 0.5), medianScale = TRUE)
power1 <- powerCrv(crvArgs = list(k = 1), medianScale = TRUE)
power2 <- powerCrv(crvArgs = list(k = 2), medianScale = TRUE)
power20 <- powerCrv(crvArgs = list(k = 20), medianScale = TRUE)
power200 <- powerCrv(crvArgs = list(k = 200), medianScale = TRUE)
powerSto <- powerCrv(crv = Stocrv, crvArgs = list(),
                     medianScale = TRUE)

## plot the curves
x <- seq(-7, -1, length.out = 7)
plot(x = x, y = x, type = 'n', # xaxs = 'i', yaxs = 'i',
     xlab = expression(paste(log[2], "(Median ",b[1], ")")), xaxt = 'n',
     ylab = expression(paste(log[2], "(Median ",b[2], ")")), yaxt = 'n')
axis(1, x, labels = x)
axis(2, x, labels = x)
abline(h = x, v = x, lty = 2, col = "gray50")
with(powerSto,
     contour(x = log(medb1, base = 2), y = log(medb2, base = 2),
             z = matrix(power, ncol = 100), add = TRUE,
             levels = c(0.5, 0.65, 0.8, 0.9), method = "simple",
             col = "black"))
