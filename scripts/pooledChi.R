source("poolFunctions.R")

## SETUP #############################################################
## visualization parameters for later
library(RColorBrewer)
kpal <- colorRamp( # define a simple palette for k values
    brewer.pal(3, "Dark2")[c(2,3)])
chiCol <- "slategray" # a colour for chi curves

## function to make a switch for temporary parameters
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
datafile <- "./results/constMPowers_fullgrid.Rds"
## datafile <- "./results/constMPowers.Rds"
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


## PLOTS #############################################################
## compute the KL divergence of each setting
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

## KL divergence with inset densities by w (Figure 4.1)
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

## IMAGE :: the distribution of lw for different w values (Figure 4.2)
set.seed(1782491)
wVals <- c(exp(-6), exp(-3), 0.5, 1)
lwDat <- lapply(wVals, # uniform samples
                function(w) matrix(runif(1e5*10), ncol = 10))
lwDists <- mapply(function(dt, w) apply(dt, 1, hrBeta, w = w),
                  lwDat, wVals) # statistic values
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

## show the impact of M on power for every KL div (Figure 4.3(a))
png("klDivPowByM.png", width = 3, height = 2.6, units = "in", res = 400)
narrowPlot(ygrid = seq(0.2, 1, by = 0.2), xlab = "ln(D(a, w))",
           xgrid = seq(-5, 5, length.out = 5), ylab = "Power",
           ylim = c(0.05, 1.02), xlim = c(-5, 5),
           mars = c(2.1, 2.1, 1.1, 3.1),
           main = expression(paste("Power of HR(", bold(p), "; w)",
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

## UMP power curves plotted vs KL divergence by w (Figure 4.3(b))
png("klDivPowByw.png", width = 2.6, height = 2.6, units = "in", res = 400)
narrowPlot(ygrid = seq(0.2, 1, by = 0.2), xlab = "ln(D(a, w))",
           xgrid = seq(-5, 5, length.out = 5), ylab = "Power",
           ylim = c(0.05, 1), xlim = c(-5, 5),
           main = expression(paste("HR(", bold(p), "; w)",
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

## mis-specified KL plots by w (Figure 4.4)
## the lines have a clear ordering based on how close they are to the
## correct specification
png("klDivPowMisSpec.png", width = 4.5, height = 3, units = "in",
    res = 400)
narrowPlot(ygrid = seq(0.2, 1, by = 0.2),
           xlab = expression(paste("ln(D(a, ", omega, "))")),
           xgrid = seq(-5, 5, length.out = 5), ylab = "Power",
           ylim = c(0.05, 1), xlim = c(-5, 5),
           mars = c(2.1, 2.1, 1.1, 4.1),
           main = expression(paste("Power of HR(", bold(p), "; w)",
                                   " by D(a, ", omega, ")")))
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

## mis-specified plots by KL divs including chi square line
## (Figure 4.10)
## note how close the chi-square test is to the UMP curve
png("klDivPowMisChi.png", width = 4.5, height = 3, units = "in",
    res = 400)
narrowPlot(ygrid = seq(0.2, 1, by = 0.2),
           xlab = expression(paste("ln(D(a, ", omega, "))")),
           xgrid = seq(-5, 5, length.out = 5), ylab = "Power",
           ylim = c(0.05, 1), xlim = c(-5, 5),
           mars = c(2.1, 2.1, 1.1, 4.1),
           main = expression(paste("Power of chi(;2981)",
                                   " and HR(;w)",
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


## SIMULATION ########################################################
## compute simulated powers when H3 is true under correct
## specification, misspecification, and using the chi method
## load simulated data for variable M or else run simulation
readData <- TRUE
datafile <- "./results/varMPowers_fullgrid.Rds"
## datafile <- "./results/varMPowers.Rds"
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


## PLOTS #############################################################
## some plotting packages
library(ggplot2)
library(grid)

## set up packages and some palettes
propPal <- colorRampPalette(c("white", "indianred4"))(100)
powPal <- colorRampPalette(c("white", "dodgerblue4"))(100)
difPal <-  colorRampPalette(c("steelblue", "floralwhite",
                              "firebrick"))(100)

## compute the KL divergence for each beta setting
klDivs <- betaDiv(params$a, params$w)

## define a helper that makes nice labels
ggLabs <- function(breaks){
    paste0(c("[", rep("(", length(breaks)-2)), # bottom brackets
           breaks[-length(breaks)], # break limits
           ", ", breaks[-1], "]") # top brackets
}

## set up further graphical parameters
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

## facetted power contours of the UMP using ggplot (Figs 4.5, 4.6)
## change ump = powershr -> fis = powersmis[,4]
##        umpPower -> fisPower
##        z = ump -> z = fis
##        HR(..., w) -> HR(..., 1)
##        D(a,w) -> D(a, omega)
##        ln(w) -> ln( omega )
## to obtain analogous plot of fisher's method
inds <- log(params$w) %in% c(-6, -3, 0) & -5 <= log(klDivs) &
    5 >= log(klDivs)
powerContourBase <- ggplot(data = data.frame(fis = powersmis[,4],
                                      lnkl = log(klDivs),
                                      eta = params$m1/10,
                                      a = params$a,
                                      w = params$w)[inds,], # temporary df
                           aes(lnkl, eta, z = fis)) # base ggplot
png("fisPower.png", width = 6, height = 2.5, units = "in", res = 400)
powGrob <- ggplotGrob(powerContourBase + # save as a grid grob
                      xlab(expression(paste("ln(D(a, ", omega, "))"))) +
                      ylab(expression(eta)) +
                      ggtitle(expression(paste("Power of HR(", bold(p),
                                               "; 1) by ln(D(a,", omega,
                                               ")) and ", eta,
                                               " facetted by ln(",
                                               omega, ")"))) +
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

## compute power differences from fis to other pooling methods
fisMin <- apply(cbind(ump = powershr, chi = powerschi,
                      mis = powersmis[,1:3]), 2,
                function(col) powersmis[,4] - col) # all in one matrix
colnames(fisMin) <- c("ump", paste0("chi", round(log(kseq), 1)),
                      paste0("mis", round(log(misSpec[1:3]), 1)))
## nice names for plotting
difTitles <- c("ump" = quote(paste("HR(", bold(p), ";", omega, ")")),
               "minchi" = quote(paste("chi", "(", bold(p), "; 0.003)")),
               "smchi" = quote(paste("chi", "(", bold(p), "; 0.018)")),
               "cvchi" = quote(paste("chi", "(", bold(p), "; 1)")),
               "fischi" = quote(paste("chi", "(", bold(p), "; 2)")),
               "modchi" = quote(paste("chi", "(", bold(p), "; 55)")),
               "bigchi" = quote(paste("chi", "(", bold(p), "; 2981)")),
               "le-6" = quote(paste("HR(", bold(p), "; 0.002)")),
               "le-3" = quote(paste("HR(", bold(p), "; 0.05)")),
               "l05" = quote(paste("HR(", bold(p), "; 0.5)")))

## display power difference contours by w (Figs 4.11, 4.7)
inds <- log(params$w) %in% c(-6, -3, 0) & -5 <= log(klDivs) &
    5 >= log(klDivs)
for (ii in seq_len(ncol(fisMin))) {
    # temporary base for the particular case
    diffContourBase <- ggplot(data = data.frame(lnkl = log(klDivs),
                                                eta = params$m1/10,
                                                a = params$a,
                                                w = params$w,
                                                dif = fisMin[,ii])[inds,],
                              mapping = aes(lnkl, eta, z = dif))
    diffGrob <- ggplotGrob( # save in a grid grob as before
        diffContourBase +
        xlab(expression(paste("ln(D(a, ", omega, "))"))) +
        ylab(expression(eta)) +
        ggtitle(bquote(paste("Power of HR(", bold(p), ";1)",
                             " minus power of ",
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

## CENTRAL/MARGINAL REJECTION ########################################
## compute the central/marginal rejection levels of l_w for three w:
## ln(w) = -6, -3, 0
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

## plots the centrality quotient for the chi method (Figure 4.9)
Mseq <- c(2, 5, 20, 100, 500, 2000, 10000) # M sequence
a <- 0.05 # alpha value
k <- exp(seq(-16, 16, by = 0.2)) # k on log scale
png("cgChiPool.png", width = 4, height = 4, units = "in",
    res = 400)
narrowPlot(xgrid = seq(-6, 6, by = 3), # set up plot area
           ygrid = seq(0, 1, by = 0.2),
           main = expression(paste("q(chi(;", kappa, ")) by M and ",
                                   log[10], "(", kappa, ")")),
           xlim = c(-7, 7),
           xlab = expression(paste(log[10], "(", kappa, ")")),
           ylab = "Centrality quotient")
for (m in Mseq) { # add lines for each
    tempQuot <- centQuot(k, M = m, alpha = a)
    lines(log10(k), tempQuot,
          col = adjustcolor(chiCol, 1))
    text(log10(k)[which.min(abs(tempQuot - 0.5))] - 0.3,
         y = 0.5, labels = m, srt = 80, cex = 0.7)
}
dev.off()

## tabulate choices of kappa to generate a chi-squared pooled p-value
## with a particular centrality quotient given a particular M
## (Table 4.2)
centSeq <- seq(0.1, 0.9, by = 0.1)
centTab <- log10(sapply(Mseq,
                      function(m) sapply(centSeq, chiKappa,
                                         M = m,
                                         interval = c(0, 1e4),
                                         tol = .Machine$double.eps^0.75)))
## automate the conversion of this to text for the thesis
paste(apply(cbind(c("", Mseq), t(cbind(centSeq, round(centTab, 1)))),
            1, paste, collapse = " & "), collapse = "\\")

## simulate random uniform data to see how chi, tippett, and
## normal pooling functions compare empirically (Figure 4.8)
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

## compare given chi method to normal, tippett methods (Fig 4.8)
kind <- 30
chikPs <- kcurves[kind,] # select k
## plot in one display
png(paste0("chiTipSto", kind, ".png"), width = 4, height = 2,
    units = "in", res = 400)
par(mar = c(2.1, 2.1, 1.1, 1.1))
plot(c(0.48*tipPs, 0.48*normPs + 0.52), c(chikPs, chikPs),
     xaxs = 'i', yaxs = 'i', xlim = c(-0.02, 1.02),
     ylim = c(-0.02, 1.02), xaxt = 'n', xlab = "", yaxt = 'n',
     ylab = "", main = "", pch = 20,
     col = adjustcolor(chiCol, 0.4))
mtext(expression(paste("chi(", bold(p)[i], ";", kappa, ")")),
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
mtext(text = seq(0, 1, length.out = 3), side = 2, cex = 0.8,
      at = seq(0, 1, length.out = 3), padj = -0.5)
## custom x axes
mtext(side = 1, at = seq(0, 0.48, length.out = 5),
      text = "|", line = 0, cex = 0.5, padj = -2)
mtext(text = seq(0, 1, length.out = 3), side = 1, cex = 0.8,
      at = seq(0, 0.48, length.out = 3))
mtext(text = expression(Tip(bold(p)[i])), at = 0.25, side = 1,
      line = 1, padj = 0, cex = 0.8)
mtext(side = 1, at = seq(0.52, 1, length.out = 5),
      text = "|", line = 0, cex = 0.5, padj = -2)
mtext(text = seq(0, 1, length.out = 3), side = 1, cex = 0.8,
      at = seq(0.52, 1, length.out = 3))
mtext(text = expression(Sto(bold(p)[i])), at = 0.75, side = 1,
      line = 1, padj = 0, cex = 0.8)
dev.off()



## OLD INVESTIGATIONS ################################################
## another aspect of this not explored is that of rejection boundaries
## and regions, which are visualized in this section for the chi
## and many others

## the following plots follow those in Loughin (2004)
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
powerSto <- powerCrv(crv = normCrv, crvArgs = list(),
                     medianScale = TRUE)

## plot the curves
x <- seq(-7, -1, length.out = 7)
plot(x = x, y = x, type = 'n', # xaxs = 'i', yaxs = 'i',
     xlab = expression(paste(log[2], "(Median ",b[1], ")")), xaxt = 'n',
     ylab = expression(paste(log[2], "(Median ",b[2], ")")), yaxt = 'n')
axis(1, x, labels = x)
axis(2, x, labels = x)
abline(h = x, v = x, lty = 2, col = "gray50")
with(power20,
     contour(x = log(medb1, base = 2), y = log(medb2, base = 2),
             z = matrix(power, ncol = 100), add = TRUE,
             levels = c(0.5, 0.65, 0.8, 0.9), method = "simple",
             col = "black"))

## this image shows the range of rejection regions
x <- exp(seq(-7, 0, by = 0.01)) # x values (log scale)
kseq <- exp(seq(-8, 8, by = 1)) # k values (log scale)
axislogk <- seq(-7, 0, by = 1) # for nice labels
chipal <- colorRampPalette( # define a palette
    brewer.pal(3, "Dark2")[c(2,3)] )(length(kseq)+2)
png("ChiRejectSeq.png", width = 5, height = 5, units = "in",
    res = 400) # create image file
## plot area
plot(x = x, y = x, type = 'n', xlab = expression(paste(log[e],p[1])),
     xlim = c(0.001, 1.2), ylim = c(0.001, 1.2),
     xaxs = "i", yaxs = "i",
     ylab = expression(paste(log[e],p[2])), log = "xy",
     main = expression(paste("Rejection boundaries for ", g[chi],
                             " by k")),
     xaxt = 'n', yaxt = 'n')
## axes and grid lines
abline(h = exp(axislogk), v = exp(axislogk), col = "gray50", lty = 2)
axis(side = 1, at = exp(axislogk), labels = axislogk)
axis(side = 2, at = exp(axislogk), labels = axislogk)
## rejection boundary lines by kappa
for (k in seq_along(kseq)) lines(x, chiCrv(x, k = kseq[k]),
                                 col = chipal[k+1])
## tippett rejection lines
lines(c(0, rep(1-0.95^0.5, 2), 1), c(1, 1, rep(1 - 0.95^0.5, 2)),
      lwd = 2, lty = 2)
text(1-0.95^0.5, 1-0.95^0.5, expression(g[Tip]), adj = c(1,1))
## stouffer normal rejection lines
lines(x, normCrv(x), lwd = 2, lty = 2)
stoInd <- which.min(abs(x - normCrv(x)))
text(x[stoInd], normCrv(x[stoInd]), expression(g[Sto]), adj = c(0,0))
## raster scale as a legend
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
source("otherMethods.R")
x <- exp(seq(-7, 0, by = 0.01)) # x values (log scale)
kseq <- exp(seq(-8, 8, by = 1)) # k values (log scale)
axislogk <- seq(-7, 0, by = 1) # for nice labels
chipal <- colorRampPalette( # define a palette
    brewer.pal(3, "Dark2")[c(2,3)] )(length(kseq)+2)
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
text(1-0.95^0.5, 1-0.95^0.5, expression(Tip), adj = c(1,1))
lines(x, Stocrv(x), lwd = 2, lty = 2)
text(x[stoInd], Stocrv(x[stoInd]), expression(Sto), adj = c(0,0))
legend(x = "bottomleft", title = "Method", lty = 1,
       legend = c(expression(Bon), expression(Tip),
                  expression(Fis), expression(Sto),
                  expression(CV), expression(Wil),
                  expression(Edg), expression(MG),
                  expression(HMP),
                  expression(paste(HMP, " (exact)"))),
       col = pal[names(methsNiceName)], bg = "white", cex = 0.8)
dev.off()
