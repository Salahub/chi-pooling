source("poolFunctions.R")

## CONSTANTS #########################################################
## necessary number of trials for a given standard error
ster <- 0.005
nsim <- (0.5/ster)^2

## given a w seq, the corresponding a values to span D(a,w)
logscaleW <- FALSE
if (logscaleW) {
    w <- exp(seq(-6, 0, by = 0.5))
} else {w <- seq(1/12, 1, length.out = 12) }
wDsets <- expand.grid(w = w, # w values
                      logD = seq(-5, 5, by = 0.125)) # target logDivs
aseq <- apply(wDsets, 1, # find a values given w, logDiv
              function(row) findA(row[1], row[2],
                                  tol = .Machine$double.eps))
logD <- log(betaDiv(aseq, wDsets[,1])) # compute actual logDiv

## a sequence of M values, k values, a and w values to parameterize
## data generation
kseq <- exp(seq(-8, 8, by = 0.25)) # k values for chi
Mval <- 80 # total number of tests
m1seq <- seq(0, Mval, by = 1)
preparams <- data.frame(a = rep(aseq, times = length(m1seq)),
                        w = rep(wDsets[,"w"], times = length(m1seq)),
                        m1 = rep(m1seq, each = length(aseq)),
                        logD = rep(logD, times = length(m1seq)))
params <- data.frame(preparams, k = rep(kseq, each = nrow(preparams)))
#altSeq <- lapply(seq_len(nrow(params)), # generate alternatives
#                 function(ii) {
#                     with(params, pgenMix(a[ii], w[ii], m1[ii]))
#                 })
#poolSeq <- lapply(seq_len(nrow(params)), # generate pool functions
#                  function(ii) {
#                      function(p) poolChi(p, k = params$k[ii])
#                  })
## control data set generation with random seeds
##seeds <- round(runif(nrow(params), min = 1e5, max = 1e9))

## helper to modify powersim
chiSimPower <- function(parRow, alpha = 0.05) {
    b <- bwa(parRow[2], parRow[1])
    rands <- cbind(matrix(rbeta(parRow[3]*nsim, parRow[1], b),
                          ncol = parRow[3], nrow = nsim),
                   matrix(runif((Mval-parRow[3])*nsim),
                          ncol = Mval - parRow[3], nrow = nsim))
    pcomb <- apply(rands, 1, poolChi, k = parRow[5])
    pow <- mean(pcomb <= alpha)
    pow
}

## get the map for the chi method (requires parallel computing)
library(parallel)
ncores <- ceiling(detectCores()*0.6)
clst <- makeCluster(ncores)
clusterExport(clst, c("bwa", "Mval", "nsim", "poolChi"))
powerschi <- parRapply(clst, params, chiSimPower)
stopCluster(clst)
#spltInds <- rep(1:ncores, # indices for splitting functions
#                each = ceiling(nrow(params)/ncores))[1:nrow(params)]
#ncores <- spltInds[length(spltInds)] # fix number of cores
#sdLst <- split(seeds, spltInds) # split seeds
#altLst <- split(altSeq, spltInds)
#poolLst <- split(poolSeq, spltInds) # functions separated by indices
#powerschi <- mcmapply(powerSim, pcomb = poolSeq, pgen = altSeq,
#                      nsim = nsim, seed = seeds, M = Mval,
#                      mc.cores = ncores) #nm = 1:ncores
                      ##mc.preschedule = FALSE, affinity.list = 1:ncores)

## store the simulation results
chiPowers <- list(pars = params,
                  chi = unlist(powerschi))
saveRDS(chiPowers, "chiPowersMap80_unifW.Rds")
