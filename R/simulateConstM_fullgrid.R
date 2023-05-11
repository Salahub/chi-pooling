source("poolFunctions.R")

## simulating t_w under H_4: all data are from the same alternative
## distribution

## COMPUTE :: necessary number of trials for a given standard error
ster <- 0.005
nsim <- (0.5/ster)^2
Mval <- 2

## COMPUTE :: given w seq, the corresponding a values to span D(a,w)
wDsets <- expand.grid(w = exp(seq(-6, 0, by = 1)), # w values
                      logD = seq(-5, 5, by = 0.5)) # target logDivs
aseq <- apply(wDsets, 1, # find a values given w, logDiv
              function(row) findA(row[1], row[2],
                                  tol = .Machine$double.eps))
logD <- log(betaDiv(aseq, wDsets[,1])) # compute actual logDiv

## SETUP :: power parameters
kseq <- c(exp(c(-8,-4,0)), 2, exp(c(4, 8))) # k values for chi
Mseq <- c(2, 5, 10, 20) # M values to simulate
params <- data.frame(a = rep(aseq, times = length(Mseq)), # par matrix
                     w = rep(wDsets[,"w"], times = length(Mseq)),
                     M = rep(Mseq, each = length(aseq)),
                     logD = rep(logD, times = length(Mseq)))
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
saveRDS(constMPowers, file = "constMPowers_fullgrid.Rds")
