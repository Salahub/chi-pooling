##' An introduction to the use of the PoolBal package and the concepts
##' of marginal and central dependence
library(PoolBal)


## GENERATING P-VALUES ###############################################
##' start by setting  the number of repeated samples (N) and size of
##' each (M)
M <- 20
N <- 200
##' as well as some parameter settings to define simulations
aSeq <- c(1, rep(0.33, 4), rep(0.99, 4))
bSeq <- c(1, rep(13.7, 4), rep(1.11, 4))
propSeq <- c(0, rep(c(0.05, 0.25, 0.5, 0.95), 2))
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
cases <- data.frame(case = 1:9, a = aSeq, b = bSeq, prop = propSeq,
                    klDiv = betaDiv(a = aSeq, b = bSeq))

##' the PoolBal package contains helpful functions which can be used
##' to generate data based on settings like these: rBetaH3 and rBetaH4
##' rBetaH3 works for this example, where a mixture of uniform and
##' beta distributions is generated
## we can use this to generate corresponding simulated samples
simPs <- do.call(rbind, # store in one matrix for ease later
                 lapply(cases$case,
                        function(ind) rBetaH3(a = aSeq[ind],
                                              b = bSeq[ind],
                                              eta = propSeq[ind],
                                              M = M,
                                              N = N)))
## and add a case number to each sample
simPs <- cbind(simPs, case = rep(1:9, each = N))


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

##' let's compare for a given case using overlaid quantile plots,
##' by adding a candidate rejection threshold at 0.05, we can get a
##' sense of the power of these pooled p-values
##' additionally, by seeing the location of the distribution of
##' pooled p-values, we get a sense of the power of each method
##' for *any* rejection threshold: the lower the distribution for
##' a pooled p-value and case the more powerful the pooled p-value
##' against that particular case
pal <- hcl.colors(6, palette = "Dark2")
caseInds <- simPs[, "case"] == 3
par(mfrow = c(1,2)) # separate chi and hr functions
plot(NA, xlim = c(0,1), ylim = c(0,1), xlab = "Quantile",
     ylab = "Pooled p-value", main = "hrPool")
abline(h = seq(0, 1, by = 0.2), v = seq(0, 1, by = 0.2),
       col = "gray50", lty = 2)
abline(a = 0, b = 1, col = "gray70")
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
abline(a = 0, b = 1, col = "gray70")
points(ppoints(N), sort(chi0.01[caseInds]), col = pal[4])
points(ppoints(N), sort(chi2[caseInds]), col = pal[5])
points(ppoints(N), sort(chi2000[caseInds]), col = pal[6])
abline(h = 0.05, col = "firebrick") # example rejection bound
legend(x = "topleft", legend = c(0.01, 2, 2000),
       title = expression(kappa), bg = "white",
       pch = 21, col = pal[4:6])
##' case 2 and case 9 are of particular interest, as they show
##' the distribution for the case of strong evidence concentrated
##' in a few tests and weak evidence spread among almost all tests,
##' respectively
##' note how w and kappa display a similar pattern: the order of
##' lines reverses between case 2 and 9, suggesting pooled p-values
##' that are powerful against concentrated strong evidence are weak
##' against diffuse weak evidence
##' this is made precise by the concepts of central and marginal
##' rejection and the centrality quotient


## ESTIMATING CENTRAL AND MARGINAL REJECTION AND CENTRALITY ##########
##' there are functions in PoolBal to estimate these quantities
##' for the hrPool functions
hrPc(1, alpha = 0.05, M = M)
hrPc(0.5, alpha = 0.05, M = M)
hrPc(0.01, alpha = 0.05, M = M)
## and chi functions
chiPc(0.01, alpha = 0.05, M = M)
chiPc(2, alpha = 0.05, M = M)
chiPc(2000, alpha = 0.05, M = M)
##' Pc (the central rejection level), is the largest p-value shared
##' among all tests that leads to rejection at alpha
##' Pr (the marginal rejection level), is the largest p-value in a
##' single test leading to rejection when all other p-values are 1
hrPr(1, alpha = 0.05, M = M)
hrPr(0.5, alpha = 0.05, M = M)
hrPr(0.01, alpha = 0.05, M = M) # hard to estimate empirically
## and chi functions
chiPr(0.01, alpha = 0.05, M = M)
chiPr(2, alpha = 0.05, M = M)
chiPr(2000, alpha = 0.05, M = M)
##' it can be proven that Pc >= Pr for any pooled p-value, which
##' justifies the centrality quotient (Pc - Pr)/Pc that measures the
##' balance between these oppositional tendencies
hrQ(1, alpha = 0.05, M = M)
hrQ(0.5, alpha = 0.05, M = M)
hrQ(0.01, alpha = 0.05, M = M) # hard to estimate empirically
## and chi functions
chiQ(0.01, alpha = 0.05, M = M)
chiQ(2, alpha = 0.05, M = M)
chiQ(2000, alpha = 0.05, M = M)


## OTHER POOLING METHODS #############################################
##' estimation of this quotient is also provided for user-defined
##' functions and those from other packages, for example, consider:
## Tippett's method (first order statistic)
tipPool <- function(p) 1 - (1 - min(p))^length(p)
## Fisher's method
fisPool <- function(p) pchisq(-2*sum(log(p)), 2*length(p),
                              lower.tail = FALSE)
## Pearson's method
peaPool <- function(p) pchisq(-2*sum(log(1 - p)), 2*length(p))
## Stouffer's method (normal transform)
stoPool <- function(p) 1 - pnorm(1/sqrt(length(p))*sum(qnorm(1 - p)))
## Wilkinson's binomial method
wilPool <- function(p, each = 0.05) 1 - pbinom(sum(p <= each),
                                           length(p), each)
## the general order statistic pooled p-value
ordPool <- function(p, k = 5) {
    pk <- sort(p)[k]
    1 - pbinom(k - 1, length(p), pk)
}
## Mudholkar and George's logit function
logPool <- function(p) {
    M <- length(p)
    scl <- pi*sqrt((M*(5*M+2))/(3*(5*M+4)))
    pt(sum(log(p/(1-p)))/scl, df = 5*M+4)
}
## the harmonic mean p-value
HMPpool <- function(p) length(p)/sum(1/p)

##' it can be proven that the centrality quotient is 0 uniquely for
##' Tippett's p-value (the minimum pooled p-value), every other
##' pooled p-value is greater than 0 when M >= 2
estimatePrb(HMPpool, M = M)
estimatePc(HMPpool, M = M)
estimateQ(HMPpool, M = M) ## harmonic mean p-value
estimatePrb(tipPool, M = M)
estimatePc(tipPool, M = M)
estimateQ(tipPool, M = M) ## minimum/tippett's
estimatePrb(fisPool, M = M)
estimatePc(fisPool, M = M)
estimateQ(fisPool, M = M) ## fisher's method

## chiPool controls this quotient
plot(seq(-10, 10, by = 0.1), chiQ(exp(seq(-10, 10, by = 0.1)), M = 20),
     type = "l", xlab = expression(log~kappa),
     ylab = "Centrality quotient")

##' as M increases, basically every pooled p-value approaches a
##' centrality of 1: can be understood geometrically or using the CLT
plot(1:20, sapply(1:20, estimateQ, poolFun = fisPool, alpha = 0.05),
     xlab = "M", ylab = "Centrality quotient", type = "b")
plot(1:20, sapply(1:20, estimateQ, poolFun = HMPpool, alpha = 0.05),
     xlab = "M", ylab = "Centrality quotient", type = "b")
