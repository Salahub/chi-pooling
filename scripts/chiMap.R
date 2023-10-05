source("poolFunctions.R")

## CONSTANTS #########################################################
## necessary number of trials for a given standard error
ster <- 0.005
nsim <- (0.5/ster)^2

## given a w seq, the corresponding a values to span D(a,w)
logscaleW <- TRUE
if (logscaleW) {
    w <- exp(seq(-6, 0, by = 0.25))
} else {w <- seq(1/20, 1, length.out = 20) }
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

## modify powerSim for the particular case of chi pooling
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

## SIMULATION ########################################################
## get the map for the chi method (requires parallel computing, takes
## a very long time to run)
library(parallel)
ncores <- ceiling(detectCores()*0.6)
clst <- makeCluster(ncores)
clusterExport(clst, c("bwa", "Mval", "nsim", "poolChi"))
powerschi <- parRapply(clst, params, chiSimPower)
stopCluster(clst)

## store the simulation results
chiPowers <- list(pars = params,
                  chi = unlist(powerschi))
filnm <- if (logscaleW) {
             "./results/chiPowerMaps80.Rds"
         } else "./results/chiPowerMaps80_unifW.Rds"
saveRDS(chiPowers, filnm)
