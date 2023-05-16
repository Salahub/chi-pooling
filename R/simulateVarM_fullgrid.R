source("poolFunctions.R")

## simulating t_w under H_3: not all data come from the null
## distribution, but all departures are of the same form

## COMPUTE :: necessary number of trials for a given standard error
ster <- 0.005
nsim <- (0.5/ster)^2

## COMPUTE :: given w seq, the corresponding a values to span D(a,w)
wDsets <- expand.grid(w = exp(seq(-6, 0, by = 0.1)), # w values
                      logD = seq(-5, 5, by = 0.1)) # target logDivs
aseq <- apply(wDsets, 1, # find a values given w, logDiv
              function(row) findA(row[1], row[2],
                                  tol = .Machine$double.eps))
logD <- log(betaDiv(aseq, wDsets[,1])) # compute actual logDiv

## SETUP: a sequence of M values, k values (for the chi method), and
## a and w values to parameterize data generation
kseq <- exp(seq(-8, 8, by = 0.5)) # k values for chi
Mval <- 20 # total number of tests
m1seq <- seq(0, Mval, by = 1)
params <- data.frame(a = rep(aseq, times = length(m1seq)),
                     w = rep(wDsets[,"w"], times = length(m1seq)),
                     m1 = rep(m1seq, each = length(aseq)),
                     logD = rep(logD, times = length(m1seq)))
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
saveRDS(varMPowers, "varMPowers_fullgrid.Rds")
