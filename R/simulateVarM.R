source("poolFunctions.R")

## simulating t_w under H_3: not all data come from the null
## distribution, but all departures are of the same form

## COMPUTE :: necessary number of trials for a given standard error
ster <- 0.005
nsim <- (0.5/ster)^2

## SETUP: a sequence of M values, k values (for the chi method), and
## a and w values to parameterize data generation
kseq <- c(exp(c(-8,-4,0)), 2, exp(c(4, 8))) # k values for chi
Mval <- 10 # total number of tests
params <- expand.grid(a = seq(0.05, 0.95, by = 0.15),
                      w = exp(seq(-6, 0, by = 1)),
                      m1 = seq(0, Mval, by = 1))
altSeq <- lapply(seq_len(nrow(params)), # generate alteratives
                 function(ii) {
                     with(params, pgenMix(a[ii], w[ii], m1[ii]))
                 })
misSpec <- exp(c(-6, -3, -log(2), 0))
## control data set generation with random seeds
seeds <- round(runif(nrow(params), min = 1e5, max = 1e9))

## SIMULATE: run the simulation using the above settings
poolHRs <- mapply(poolHR, w = params$w, nsim = 1e5, M = Mval,
                  SIMPLIFY = FALSE)
poolMis <- mapply(poolHR, w = rep(misSpec, each = nrow(params)),
                  nsim = 1e5, M = Mval, SIMPLIFY = FALSE)
## compute powers
powershr <- mapply(powerSim, pcomb = poolHRs, pgen = altSeq,
                   nsim = nsim, M = Mval, seed = seeds)
powersmis <- mapply(powerSim, pcomb = poolMis,
                    pgen = rep(altSeq, times = length(misSpec)),
                    nsim = nsim, M = Mval, seed = seeds)
powerschi <- sapply(kseq, ## for the chi method
                    function(k) {
                        mapply(powerSim,
                               pcomb = list(function(p) poolChi(p, k)),
                               pgen = altSeq,
                               nsim = nsim,
                               seed = seeds,
                               M = Mval)
                    })

## OUTPUT :: output the simulation results
varMPowers <- list(pars = params, misw = misSpec, ks = kseq,
                   tw = powershr, chi = powerschi,
                   mis = matrix(powersmis, ncol = length(misSpec),
                                byrow = FALSE))
saveRDS(varMPowers, "varMPowers.Rds")
