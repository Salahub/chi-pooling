## get the central rejection level
estimatePc <- function(poolFun, alpha = 0.05, M = 2,
                       poolArgs = list(), ...) {
    pcf <- function(x) do.call(poolFun,
                               args = c(p = list(rep(x, M)),
                                        poolArgs)) - alpha
    uniroot(pcf, ...)
}

## get the marginal rejection level at b
estimatePrb <- function(poolFun, alpha = 0.05, b = 1, M = 2,
                       poolArgs = list(), ...) {
    pcr <- function(x) do.call(poolFun,
                               args = c(p = list(c(x, rep(b, M-1))),
                                        poolArgs)) - alpha
    uniroot(pcr, ...)
}

## the centrality quotient for a pooled p-value
estimateQ <- function(poolFun, alpha = 0.05, M = 2,
                      poolArgs = list(), ...) {
    pc <- estimatePc(poolFun, alpha = alpha, M = M,
                     poolArgs = poolArgs, ...)$root
    pr <- estimatePc(poolFun, alpha = alpha, M = M,
                     poolArgs = poolArgs, ...)$root
    (pc - pr)/pc
}
