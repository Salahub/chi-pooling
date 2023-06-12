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
nsim <- 2000
ngroup <- 8
npop <- 104
minQuants <- readRDS("curveMinQuantiles.Rds")
kseq <- exp(seq(-8, 8, by = 0.2))
thetas <- seq(60, 70, by = 0.01)
muMat <- sdMat <- matrix(nrow = nsim, ncol = ngroup)
nMat <- matrix(nrow = nsim, ncol = ngroup)
pMat <- array(dim = c(nsim, ngroup, length(thetas)))
poolMat <- array(dim = c(nsim, length(kseq), length(thetas)))
set.seed(9381278)
for (ii in 1:nsim) {
    ## start with simulated data
    obs <- rnorm(npop, mean = 65, sd = 5)
    ## randomly split
    group <- sample(rep(1:ngroup, times = npop/ngroup))
    #group <- sample(c(sample(1:ngroup),
    #                  sample(1:ngroup, npop - ngroup, replace = TRUE)))
    nMat[ii,] <- table(group)
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
## compute the minimum by theta
minpool <- apply(poolMat, c(1,3), min)

## plot a realization
cols <- RColorBrewer::brewer.pal(4, "Dark2")
real <- sample(1:nsim, 1)
plot(thetas,  minpool[real,], type = 'l', ylim = c(0,1),
     col = adjustcolor(cols[4], 0.8), xlab = expression(theta),
     ylab = expression(g[chi]), xlim = c(60, 70))
for (ii in c(10, 41)) lines(thetas, poolMat[real,ii,], type = 'l',
                         col = if (ii <= 40) adjustcolor(cols[1], 0.8) else if (ii == 41) adjustcolor(cols[2], 0.8) else adjustcolor(cols[3], 0.8))
                                        #lwd = if (log(kseq)[ii] == 0) 2 else 1)
lines(thetas, poolMat[real,81,], col = adjustcolor(cols[3], 0.8))
abline(v = 65, col = adjustcolor("black", 0.4), lty = 2, lwd = 2)
abline(v = mean(muMat[real,]), col = adjustcolor("black", 0.4), lwd = 2)
abline(v = muMat[real,], col = adjustcolor("gray50", 0.4))
abline(h = minQuants[c("5%"), "100"], lty = 3)

## 789, 1061, 1192, 1921

## coverage probabilities, estimates for a range of cutoffs
safeRange <- function(x) {
    if (length(x) == 0) {
        c(NA, NA)
    } else range(x)
}
cutoffs <- minQuants[c("0.1%", "1%", "5%", "10%"), "100"]
poolTheta <- matrix(thetas[apply(poolMat, c(1,2), which.max)],
                    ncol = length(kseq))
poolIntervals <- lapply(cutoffs,
                        function(ct) apply(poolMat >= ct, c(1,2),
                                           function(x) {
                                               safeRange(thetas[which(x)])
                                           }))
poolInclude <- lapply(poolIntervals,
                      function(mat) mat[1,,] <= 65 & mat[2,,] >= 65)
poolCovP <- lapply(poolInclude,
                   function(mat) apply(mat, 2, mean, na.rm = TRUE))
## using the mean
meanTheta <- apply(muMat, 1, mean)
meansd <- sqrt(apply(sdMat^2*nMat*(nMat-1), 1, sum)/(npop-1))/sqrt(npop)
meanInt <- cbind(meanTheta - qnorm(0.975)*meansd,
                 meanTheta + qnorm(0.975)*meansd)
meanCovP <- mean(meanInt[,1] <= 65 & meanInt[,2] >= 65)

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
