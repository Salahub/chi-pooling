## some pooling functions
## 1. the chi-squared quantile pooled p-value
poolChi <- function(p, k) {
    M <- length(p) # dimension
    pchisq(sum(qchisq(p, df = k, lower.tail = FALSE)), df = M*k,
           lower.tail = FALSE)
}
## 2. stouffer's normal pooled p-value
poolNorm <- function(p, mu = 0, sd = 1) {
    M <- length(p) # dimension
    pnorm(sum(qnorm(p, mean = mu, sd = sd, lower.tail = FALSE)),
          mean = M*mu, sd = sqrt(M)*sd, lower.tail = FALSE)
}
## 3. a general gamma pooled p-value
poolGamma <- function(p, shape = 1, rate = 1/2) {
    M <- length(p) # dimension
    pgamma(sum(qgamma(p, shape = shape, rate = rate,
                      lower.tail = FALSE)),
           shape = M*shape, rate = rate,
           lower.tail = FALSE) # sum shapes
}
## 4. tippett's minimum pooled p-value
poolTip <- function(p) {
    1 - (1 - min(p))^(length(p))
}

## central and marginal rejection levels for the chi-squared
## pooled p-value
## central rejection level
chiPc <- function(k, M, alpha = 0.05) {
    pchisq(qchisq(alpha, df = M*k, lower.tail = FALSE)/M, df = k,
           lower.tail = FALSE)
}
## marginal rejection level
chiPr <- function(k, M, alpha = 0.05) {
    pchisq(qchisq(alpha, df = M*k, lower.tail = FALSE), df = k,
           lower.tail = FALSE)
}

## two ways of computing the centrality quotient:
## conditional distribution
centQuot <- function(k, M, alpha = 0.05) {
    chiCrit <- qchisq(alpha, df = M*k, lower.tail = FALSE)
    1 - pchisq(chiCrit, df = k, lower.tail = FALSE)/
        pchisq(chiCrit/M, df = k, lower.tail = FALSE)
}
## directly using previous functions (seems more stable)
centQuotRat <- function(k, M, alpha = 0.05) {
    1 - chiPr(k, M, alpha)/chiPc(k, M, alpha)
}

## use uniroot to compute the kappa value that gives a particular
## centrality quotient in the chi-squared pooled p-value
chiKappa <- function(cq, M, alpha = 0.05, interval = c(0,100),
                     tol = .Machine$double.eps^0.5) {
    uniroot(function(x) centQuotRat(x, M, alpha) - cq,
            interval = interval, tol = tol)$root
}

## functions for the UMP pooled p-value
## 1. the statistic given p and w
hrBeta <- function(p, w = 1) {
    w*sum(log(p)) - (1 - w)*sum(log(1 - p))
}
## 2. estimated rejection bounds for the UMP
hrBnd <- function(w, alpha = 0.05, M = 2, nsim = 1e4) {
    dat <- matrix(runif(M*nsim), ncol = M) # simulate data
    vals <- apply(dat, 1, hrBeta, w = w)
    quantile(vals, 1 - alpha) # get simulated quantiles
}
## 3. estimated pooled p-value, produces a closure to avoid repeated
## simulation (which is time instensive)
poolHR <- function(w = 1, M = 2, nsim = 1e5) { # simulated pooled p
    dat <- matrix(runif(M*nsim), ncol = M)
    pools <- apply(dat, 1, hrBeta, w = w)
    function(p) { # closure to limit repeated computation
        mean(hrBeta(p, w) >= pools) # obs quant
    }
}
## 4. a univariate version of the UMP for central rejection
pcHR <- function(w = 1, M = 2, nsim = 1e5) {
    dat <- matrix(runif(M*nsim), ncol = M)
    pools <- apply(dat, 1, hrBeta, w = w)
    function(pc) {
        mean(hrBeta(rep(pc, M), w) >= pools)
    }
}
## 5. univariate version of the UMP for marginal rejection
prHR <- function(w = 1, M = 2, nsim = 1e5) {
    dat <- matrix(runif(M*nsim), ncol = M)
    pools <- apply(dat, 1, hrBeta, w = w)
    function(pr) {
        mean(hrBeta(c(pr, rep(0.999, M-1)), w) >= pools)
    }
}

## H4 alternative generation function assuming a beta distribution
pgenHR <- function(a, w) { # accept parameters
    b <- bwa(w, a)
    function(M, nsim) { # closure to generate values
        ps <- rbeta(M*nsim, shape1 = a, shape2 = b)
        matrix(ps, ncol = M) # organize output
    }
}

## H3 alternative generation function assuming a beta distribution
pgenMix <- function(a, w, m1) {
    b <- bwa(w, a)
    function(M, nsim) {
        h1 <- rbeta(m1*nsim, shape1 = a, shape2 = b) # true non-nulls
        h0 <- runif((M - m1)*nsim) # true nulls
        cbind(matrix(h1, ncol = m1, nrow = nsim),
              matrix(h0, ncol = M - m1, nrow = nsim))
    }
}

## the next functions define bivariate curves that show the
## rejection boundaries of the different methods
## 1. the curve for the chi-squared pooled p-value
chiCrv <- function(x, k = 2, alpha = 0.05) {
    pchisq(qchisq(alpha, df = 2*k, lower.tail = FALSE) -
           qchisq(x, df = k, lower.tail = FALSE),
           df = k, lower.tail = FALSE)
}
## 2. the curve for stouffer's normal pooled p-value
normCrv <- function(x, alpha = 0.05, mu = 0, sd = 1) {
    pnorm(qnorm(alpha, mean = mu, sd = sqrt(2*sd),
                lower.tail = FALSE) -
          qnorm(x, mean = mu, sd = sd, lower.tail = FALSE),
          mean = mu, sd = sd, lower.tail = FALSE)
}
## 3. tippett's curve
tipCrv <- function(x, alpha = 0.05) {
    bnd <- 1 - (1 - alpha)^(1/2) # the boundary
    below <- x <= bnd
    x[below] <- 1
    x[!below] <- bnd
    x
}
## 4. curve for the gamma pooled p-value
gamCrv <- function(x, shape = 1, rate = 1/2, alpha = 0.05) {
    pgamma(qgamma(alpha, shape = 2*shape, rate = rate,
                  lower.tail = FALSE) -
           qgamma(x, shape = shape, rate = rate, lower.tail = FALSE),
           shape = shape, rate = rate, lower.tail = FALSE)
}

## compute the Kullback-Leibler divergence
klDiv <- function(f1, f2, lower = 0, upper = 1) {
    div <- function(x) { # define a divergence function pointwise
        vals <- f1(x)*log(f1(x)/f2(x))
        vals[f1(x) == 0] <- 0 # by limiting arguments
        vals
    }
    integrate(div, lower = lower, upper = upper) # integrate
}

## compute the KL div for the beta case in particular
betaDiv <- function(a, w) {
    b <- bwa(w, a) # second parameter
    lbeta(a, b) + a + b - 2 # compute
}

## a helper to compute b given the UMP parameter w and a
bwa <- function(w, a) 1/w + (1 - 1/w)*a

## use uniroot to find the a valu that gives a particular log beta
## divergence given w
findA <- function(w, logd = 0, ...) {
    scr <- function(a) log(betaDiv(a, w)) - logd
    uniroot(scr, c(0,1), ...)$root
}

## another helper to repeat loughin's work
medianP <- function(b) {
    1 - 0.5^(1/b)
}

## beta power function: computes the power for a chi square pool
## over a grid of locations assuming independent beta variables are
## the true underlying distribution of p1, p2
## following Loughin (2004), a = 1 for all betas
powerCrv <- function(crv = chiCrv, crvArgs = list(k = 2),
                     alpha = 0.05, medianScale = FALSE,
                     b1seq = seq(1, 100, length.out = 100),
                     b2seq = seq(1, 100, length.out = 100)) {
    ## define the integrating power function
    pow <- function(bs) {
        b1 <- bs[1]
        b2 <- bs[2] # bs defined within a row
        integrate(function(x) {
            sapply(x, function(y) { # handle bivariate itegration
                dbeta(y, shape1 = 1, shape2 = b1)*
                    integrate(dbeta, lower = 0,
                              upper = do.call(crv, c(x = y,
                                                     alpha = alpha,
                                                     crvArgs)),
                              shape1 = 1, shape2 = b2)$value
            })
        }, lower = 0, upper = 1)$value
    }
    if (!medianScale) { # on the scale of the parameters
        ## grid based on the given specifications
        bcomb <- expand.grid(b1seq, b2seq)
        list(b1 = b1seq, b2 = b2seq,
             power = apply(bcomb, 1, pow)) # get the powers
    } else { # median p-value scale
        # make sure order is still increasing
        bcomb <- expand.grid(rev(b1seq), rev(b2seq))
        power <- apply(bcomb, 1, pow) # compute powers
        list(medb1 = medianP(rev(b1seq)), medb2 = medianP(rev(b2seq)),
             power = power)
    }
}

## a general function to simulate power based on a p generation
## function, a combination function, and a prescribed alpha level the
## combination function should be calibrated to, i.e. when it is less
## than or equal to alpha that should lead to rejection
powerSim <- function(pcomb, pgen, alpha = 0.05, M = 2, nsim = 1e4,
                     seed = as.numeric(Sys.Date())*Sys.getpid()) {
    set.seed(seed)
    ps <- pgen(M, nsim)
    pooled <- apply(ps, 1, pcomb)
    mean(pooled <= alpha) # proportion of rejected tests
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
    abline(h = ygrid, v = xgrid, lty = 1,
           col = adjustcolor("gray", alpha.f = 0.4))
    ## and ticks
    mtext(side = 1, at = xgrid, text = "|", line = 0, cex = 0.5,
          padj = -2)
    mtext(text = xticks, at = xgrid, side = 1, cex = 0.8)
    mtext(side = 2, at = ygrid, text = "|", line = 0, cex = 0.5,
          padj = 1)
    mtext(text = yticks, at = ygrid, side = 2, cex = 0.8)
}

## custom image plot with narrower margins
narrowImage <- function(mat, col, xnames = rownames(mat),
                        ynames = colnames(mat), main = "",
                        xlab = names(dimnames(mat))[1],
                        ylab = names(dimnames(mat))[2],
                        las = 1, ...) {
    if (missing(col)) {
        col <- colorRampPalette(c("white", "indianred4"))(15)
    }
    image(mat, col = col, xaxt = "n", yaxt = "n", ...) # add image
    bds <- par()$usr # user parameters
    mtext(main, side = 3, line = 0, cex = 0.8) # main
    mtext(ylab, side = 2, line = 1, cex = 0.8) # ylab
    mtext(xlab, side = 1, line = 1, padj = 0, cex = 0.8) # xlab
    ## add ticks
    xtcks <- seq(0, 1, length.out = length(xnames))
    ytcks <- seq(0, 1, length.out = length(ynames))
    mtext(text = xnames, at = xtcks, side = 1, cex = 0.8,
          las = las)
    mtext(text = ynames, at = ytcks, side = 2, cex = 0.8,
          las = las, adj = 1.2)
    ## separating lines
    abline(v = c(filter(xtcks, c(0.5, 0.5)), bds[c(1,2)]),
           h = c(filter(ytcks, c(0.5, 0.5)), bds[c(3,4)]))
}

## helper that provides multipliers to keep a given aspect ratio for
## an active plot on an active device
aspectRatio <- function(asp = 1, units = "in") {
    devln <- dev.size(units) # capture the current device dimensions
    pltlm <- diff(par()$plt)[c(1,3)] # plot limits on device
    pltd <- diff(par()$usr)[c(1,3)] # plot dimensions in plot units
    pltd/(devln*pltlm)*c(1, asp) # plot unit/length adjusted for asp
}

## a glyph that displays a density
densityGlyph <- function(f, xseq, x, y, size = 0.5, units = "in",
                         asp = 1, devRat = aspectRatio(asp, units),
                         flim = 2, ...) {
    dims <- size*devRat # get dim in plot units
    xleft <- x - dims[1]/2
    xright <- x + dims[1]/2
    ybot <- y - dims[2]/2
    ytop <- y + dims[2]/2 # convert to bounds
    fseq <- f(xseq) # apply function to x sequence
    fseq[fseq > flim] <- flim
    rect(xleft, ybot, xright, ytop) # bounding box
    polygon(c(xleft, xleft + (xright - xleft)*xseq/max(xseq), xright),
            c(ybot, ybot + (ytop - ybot)*fseq/max(fseq), ybot), ...)
}
