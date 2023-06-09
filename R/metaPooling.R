library(metadat)

## parameter sweep function
paramSweep <- function(mus, sds, thetas) {
    mus2 <- rep(mus, each = length(thetas))
    sds2 <- rep(sds, each = length(thetas))
    matrix(2*pnorm(abs((mus2 - rep(thetas, times = length(mus)))/sds2),
                   lower.tail = FALSE),
           nrow = length(thetas))
}

## sweeping with kappa afterwards
kappaSweep <- function(ps, kseq = exp(seq(-8, 8, by = 0.1))) {
    sapply(kseq,
           function(kap) {
               apply(ps, 2,
                     function(p) pchisq(sum(qchisq(p, kap,
                                               lower.tail = FALSE)),
                                        length(p)*kap,
                                        lower.tail = FALSE))
           })
}

## repeat this process many times to compare bias...
nsim <- 1000
ngroup <- 8
npop <- 200
minQuants <- readRDS("curveMinQuantiles.Rds")
kseq <- exp(seq(-8, 8, by = 0.2))
thetas <- seq(60, 70, by = 0.02)
muMat <- sdMat <- matrix(nrow = nsim, ncol = ngroup)
pMat <- array(dim = c(nsim, ngroup, length(thetas)))
poolMat <- array(dim = c(nsim, length(kseq), length(thetas)))
for (ii in 1:nsim) {
    ## start with simulated data
    obs <- rnorm(npop, mean = 65, sd = 5)
    ## randomly split
    group <- sample(1:8, 200, replace = TRUE,
                    prob = c(2, 1, 4, 3, 1, 1, 2, 3))
    ## get means and sds from split groups
    muMat[ii,] <- mus <- tapply(obs, group, mean)
    ## get sds from split groups
    sdMat[ii,] <- sds <- tapply(obs, group, sd)/sqrt(table(group))
    ## simulated p-values
    pMat[ii,,] <- simps <- t(paramSweep(mus, sds, thetas))
    ## the pooled p-value curves
    poolMat[ii,,] <- simpool <- t(kappaSweep(simps, kseq = kseq))
    if ((ii %% 10) == 0) cat("\r Done ", ii, " of ", nsim)
}

## plot a realization
real <- sample(1:nsim, 1)
plot(thetas, poolMat[real,81,], type = 'l', ylim = c(0,1),
     col = adjustcolor("firebrick", 0.25), xlab = expression(theta),
     ylab = expression(g[chi]))
for (ii in 10:80) lines(thetas, poolMat[real,ii,], type = 'l',
                         col = if (ii <= 40) adjustcolor("steelblue", 0.25) else if (ii == 41) adjustcolor("black", 0.25) else adjustcolor("firebrick", 0.25),
                         lwd = if (log(kseq)[ii] == 0) 2 else 1)
abline(v = 65, col = "black", lty = 2, lwd = 2)
abline(v = mean(muMat[real,]), col = "black")
abline(v = muMat[real,], col = "gray50")
abline(h = 0.05, lty = 3)

## coverage probabilities, estimates
poolTheta <- apply(poolMat, c(1,2), max)
poolAccept <- apply(poolMat >= 0.05, c(1,2), which)
poolInclude <- apply(poolAccept, c(1,2), function(inds) 251 %in% unlist(inds))
poolCovP <- apply(poolInclude, 2, mean)

## middle line (kappa = 1) looks normal...
## can we fit a normal curve to it?
## suggested: nls() with c phi((theta - b)/a), a,b,c params

## look at the covid dataset
data(dat.axfors2021) ## requires estimating odds ratios
data(dat.normand1999) ## differences in means but classic otherwise
data(dat.lim2014) ## offspring size data (counts)
data(dat.franchini2012) ## parkinsons, means and sds
data(dat.bakdash2021) ## big but multi-leveled, pretty involved

## try this simple idea out
normandps <- paramSweep(dat.normand1999$m1i, dat.normand1999$sd1i,
                        range = c()
