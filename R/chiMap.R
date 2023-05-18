source("poolFunctions.R")

## CONSTANTS #########################################################
## necessary number of trials for a given standard error
ster <- 0.005
nsim <- (0.5/ster)^2

## given a w seq, the corresponding a values to span D(a,w)
wDsets <- expand.grid(w = exp(seq(-6, 0, by = 1)), # w values
                      logD = seq(-5, 5, by = 0.25)) # target logDivs
aseq <- apply(wDsets, 1, # find a values given w, logDiv
              function(row) findA(row[1], row[2],
                                  tol = .Machine$double.eps))
logD <- log(betaDiv(aseq, wDsets[,1])) # compute actual logDiv

## a sequence of M values, k values, a and w values to parameterize
## data generation
kseq <- exp(seq(-8, 8, by = 0.5)) # k values for chi
Mval <- 40 # total number of tests
m1seq <- seq(0, Mval, by = 1)
preparams <- data.frame(a = rep(aseq, times = length(m1seq)),
                        w = rep(wDsets[,"w"], times = length(m1seq)),
                        m1 = rep(m1seq, each = length(aseq)),
                        logD = rep(logD, times = length(m1seq)))
params <- data.frame(preparams, k = rep(kseq, each = nrow(preparams)))
altSeq <- lapply(seq_len(nrow(params)), # generate alternatives
                 function(ii) {
                     with(params, pgenMix(a[ii], w[ii], m1[ii]))
                 })
poolSeq <- lapply(seq_len(nrow(params)), # generate pool functions
                  function(ii) {
                      function(p) poolChi(p, k = params$k[ii])
                  })
## control data set generation with random seeds
seeds <- round(runif(nrow(params), min = 1e5, max = 1e9))
## helper to modify powersim
chiSimPower <- function(altLst, poolLst, nsim, seeds, M) {
    mapply(powerSim, pcomb = poolLst, pgen = altLst,
           nsim = nsim, seed = seeds, M = M)
}

## get the map for the chi method (requires parallel computing)
library(parallel)
ncores <- detectCores()/2
spltInds <- rep(1:ncores, # indices for splitting functions
                each = ceiling(nrow(params)/ncores))[1:nrow(params)]
sdLst <- split(seeds, spltInds) # split seeds
altLst <- split(altSeq, spltInds)
poolLst <- split(poolSeq, spltInds) # functions separated by indices
powerschi <- mcmapply(chiSimPower, altLst = altLst, poolLst = poolLst,
                      nsim = nsim, seeds = sdLst, M = Mval,
                      mc.cores = ncores)

## store the simulation results
chiPowers <- list(pars = params,
                  chi = unlist(powerschi))
saveRDS(chiPowers, "chiPowersMap.Rds")
