##' An introduction to the use of the PoolBal package and the concepts
##' of marginal and central dependence
library(PoolBal)


## GENERATING P-VALUES ###############################################
##' start by setting  the number of repeated samples (N) and size of
##' each (M)
M <- 20
N <- 100
##' the PoolBal package provides two helpful functions to simulate
##' samples which are mixtures of beta and uniform distributions
##' as well as some parameter settings to feed into rUnifBeta
aSeq <- c(1, rep(0.05, 3), rep(0.8, 3))
bSeq <- c(1, rep(5, 3), rep(1, 3))
propSeq <- c(0, rep(c(0.05, 0.5, 0.95), 2))
##' these seven settings correspond to:
##' 1. completely uniform random data
##' 2. strong non-null evidence in one test
##' 3. strong non-null evidence in half the tests
##' 4. strong non-null evidence in all but one test
##' 5. weak non-null evidence in one test
##' 6. weak non-null evidence in half the tests
##' 7. weak non-null evidence in all but one test
##' "strong" and "weak" are here defined based on the Kullback-Leibler
##' divergence from uniform of the non-null distribution, which
##' can be computed for the beta case using betaDiv from PoolBal
betaDiv(a = aSeq, b = bSeq)
## all of this can be combined into a data.frame of cases/settings
cases <- data.frame(case = 1:7, a = aSeq, b = bSeq, prop = propSeq,
                    klDiv = betaDiv(a = aSeq, b = bSeq))
## we can use this to generate corresponding simulated samples
simPs <- do.call(rbind, # store in one matrix for ease later
                 lapply(cases$case,
                        function(ind) rUnifBeta(N, M,
                                                prop = propSeq[ind],
                                                a = aSeq[ind],
                                                b = bSeq[ind])))
## and add a case number to each sample
simPs <- cbind(simPs, case = rep(1:7, each = N))


## POOLING THE P-VALUES ##############################################
##' the two pooled p-value functions provided in PoolBal can be used
##' to process these pooled p-values
##' chiPool uses default R functions to do this without simulation,
##' while hrPool is a closure which runs and stores simulations in
##' advance to save time and so must be defined before use
##' hrPool requires the sample size (M), the number of simulations to
##' compute the empirical p-value (nsim), and a w parameter explained
##' in Heard & Rubin-Delanchy (2018)
hr1 <- hrPool(w = 1, M = M, nsim = 1e3)
hr.5 <- hrPool(w = 0.5, M = M, nsim = 1e3)
hr.01 <- hrPool(w = 0.01, M = M, nsim = 1e3)
##' each of these gives a slightly different distribution of pooled
##' p-values when applied to each data set
hr1Dists <- apply(simPs[, 1:M], 1, hr1)
hr.5Dists <- apply(simPs[, 1:M], 1, hr.5)
hr.01Dists <- apply(simPs[, 1:M], 1, hr.01)

##' we can do the same thing with chiPool, parameterized by the
##' degrees of freedom kappa which can range from 0 to infinity
chi0.01 <- apply(simPs[, 1:M], 1, chiPool, kappa = 0.01)
chi2 <- apply(simPs[, 1:M], 1, chiPool, kappa = 2)
chi2000 <- apply(simPs[, 1:M], 1, chiPool, kappa = 2000)

##' let's compare for a given case using overlaid quantile plots
pal <- hcl.colors(6, palette = "Dark2")
caseInds <- simPs[, "case"] == 7
par(mfrow = c(1,2)) # separate chi and hr functions
plot(NA, xlim = c(0,1), ylim = c(0,1), xlab = "Quantile",
     ylab = "Pooled p-value", main = "hrPool")
abline(h = seq(0, 1, by = 0.2), v = seq(0, 1, by = 0.2),
       col = "gray50", lty = 2)
points(ppoints(N), sort(hr1Dists[caseInds]), col = pal[1])
points(ppoints(N), sort(hr.5Dists[caseInds]), col = pal[2])
points(ppoints(N), sort(hr.01Dists[caseInds]), col = pal[3])
abline(h = 0.05, col = "firebrick") # example rejection bound
legend(x = "topleft", legend = c(1, 0.5, 0.01), title = "w",
       pch = 21, col = pal[1:3], bg = "white")
## now the chi quantile plot
plot(NA, xlim = c(0,1), ylim = c(0,1), xlab = "Quantile",
     ylab = "Pooled p-value", main = "chiPool")
abline(h = seq(0, 1, by = 0.2), v = seq(0, 1, by = 0.2),
       col = "gray50", lty = 2)
points(ppoints(N), sort(chi0.01[caseInds]), col = pal[4])
points(ppoints(N), sort(chi2[caseInds]), col = pal[5])
points(ppoints(N), sort(chi2000[caseInds]), col = pal[6])
abline(h = 0.05, col = "firebrick") # example rejection bound
legend(x = "topleft", legend = c(0.01, 2, 2000),
       title = expression(kappa), bg = "white",
       pch = 21, col = pal[4:6])


##' this section defines some previous pooled p-value functions for
##' reference
## Bonferroni's method
bonPool <- function(p) min(1, length(p)*min(p))
## Tippett's method (first order statistic)
tipPool <- function(p) 1 - (1 - min(p))^length(p)
## Fisher's method
fisPool <- function(p) pchisq(-2*sum(log(p)), 2*length(p),
                              lower.tail = FALSE)
## Pearson's method
peaPool <- function(p) pchisq(-2*sum(log(1 - p)), 2*length(p))
## Stouffer's method (normal transform)
Sto <- function(p) 1 - pnorm(1/sqrt(length(p))*sum(qnorm(1 - p)))
## Wilkinson's binomial method
Wil <- function(p, each = 0.05) 1 - pbinom(sum(p <= each),
                                           length(p), each)
## the general order statistic pooled p-value
Ord <- function(p, k = 5) {
    pk <- sort(p)[k]
    1 - pbinom(k - 1, length(p), pk)
}
## Mudholkar and George's logit function
Log <- function(p) {
    M <- length(p)
    scl <- pi*sqrt((M*(5*M+2))/(3*(5*M+4)))
    pt(sum(log(p/(1-p)))/scl, df = 5*M+4)
}
## the harmonic mean p-value
HMP <- function(p) length(p)/sum(1/p)

