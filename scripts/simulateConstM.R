source("poolFunctions.R")

## simulating t_w under H_4: all data are from the same alternative
## distribution

## COMPUTE :: necessary number of trials for a given standard error
ster <- 0.005
nsim <- (0.5/ster)^2
Mval <- 2

## SETUP :: power parameters for the HR setting
kseq <- c(exp(c(-8,-4,0)), 2, exp(c(4, 8))) # k values for chi
Mseq <- c(2, 5, 10, 20) # M values to simulate
params <- expand.grid(a = seq(0.05, 0.95, by = 0.15), # params for beta
                      w = exp(seq(-6, 0, by = 1)),
                      M = c(2, 5, 10, 20))
betaSeq <- lapply(seq_len(nrow(params)), # beta gen functions
                  function(ii) pgenHR(params$a[ii], params$w[ii]))
misSpec <- exp(c(-6, -3, -log(2), 0))
## control data set generation with random seeds
seeds <- round(runif(nrow(params), min = 1e5, max = 1e9))

## SIMULATE :: run the simulations
## function giving the hr pooled value for each w
poolHRs <- mapply(poolHR, w = params$w, nsim = 1e5, M = params$M,
                  SIMPLIFY = FALSE)
poolMis <- mapply(poolHR, w = rep(misSpec, each = nrow(params)),
                  nsim = 1e5,
                  M = rep(params$M, times = length(misSpec)),
                  SIMPLIFY = FALSE) # misspecified w
## get powers by w
powershr <- mapply(powerSim, pcomb = poolHRs, pgen = betaSeq,
                   nsim = nsim, M = params$M, seed = seeds)
powersmis <- mapply(powerSim, pcomb = poolMis, # misspecified w
                    pgen = rep(betaSeq, times = length(misSpec)),
                    nsim = nsim, seed = seeds,
                    M = rep(params$M, times = length(misSpec)))
powerschi <- sapply(kseq, # for the chi method
                    function(k) {
                        mapply(powerSim,
                               pcomb = list(function(p) poolChi(p, k)),
                               pgen = betaSeq,
                               nsim = nsim,
                               seed = seeds,
                               M = params$M)
                    })

## OUTPUT :: output the simulation results
constMPowers <- list(pars = params, misw = misSpec, ks = kseq,
                     tw = powershr, chi = powerschi,
                     mis = matrix(powersmis, ncol = length(misSpec),
                                  byrow = FALSE))
saveRDS(constMPowers, file = "constMPowers.Rds")
