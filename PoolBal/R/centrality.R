##' @title Compute the central rejection level
##' @description Estimates the central rejection level for an
##' arbitrary pooled p-value function.
##' @details The central rejection level is the maximum p-value
##' shared among all tests which still results in rejection of the
##' null using a pooled p-value.
##'
##' This function is essentially a wrapper for uniroot, and accepts
##' a pooling function which takes a numeric vector as its first
##' argument and potentially other arguments given in poolArgs and
##' returns a single value. Using this pooling function, a specified
##' dimension M and a rejection level alpha, uniroot searches for the
##' root to poolFun - alpha along the line where all p-values are
##' equal.
##' @param poolFun function accepting a vector of p-values
##' @param alpha numeric between 0 and 1
##' @param M integer, how many p-values are there?
##' @param interval two numerics giving the bounds of root-searching
##' @param poolArgs (optional) additional named arguments for poolFun
##' @param ... additional arguments to uniroot
##' @return The uniroot output.
##' @examples
##' tippool <- function(p) 1 - (1 - p)^(length(p))
##' estimatePc(tippool, 0.05, M = 10, interval = c(0, 1))
##' @author Chris Salahub
estimatePc <- function(poolFun, alpha = 0.05, M = 2,
                       interval = c(0, 1), poolArgs = list(), ...) {
    pcf <- function(x) do.call(poolFun,
                               args = c(list(rep(x, M)),
                                        poolArgs)) - alpha
    rt <- try(uniroot(pcf, interval = interval, ...)$root,
              silent = TRUE)
    if ((is(rt, "try-error")|is(rt, "error")) &
        interval[1] < 2*.Machine$double.eps &
        interval[2] > 2*.Machine$double.eps) { # basically 0 and 1
        if (grepl("values at end points not of opposite sign", rt)) {
            0 # handle floating point issues without errors
        } else {
            stop(rt)
        }
    } else {
        rt
    }
}

##' @title Compute the marginal rejection level
##' @description Estimates the marginal rejection level for an
##' arbitrary pooled p-value function.
##' @details The marginal rejection level is the maximum p-value
##' in a single test less than b which still results in rejection of
##' the null when all other tests have a p-value of b.
##'
##' This function is essentially a wrapper for uniroot, and accepts
##' a pooling function which takes a numeric vector as its first
##' argument and potentially other arguments given in poolArgs and
##' returns a single value. Using this pooling function, a specified
##' dimension M and a rejection level alpha, uniroot searches for the
##' root to poolFun - alpha along one margin when all other p-values
##' are equal to b.
##' @param poolFun function accepting a vector of p-values
##' @param alpha numeric between 0 and 1
##' @param M integer, how many p-values are there?
##' @param interval two numerics giving the bounds of root-searching
##' @param poolArgs (optional) additional named arguments for poolFun
##' @param ... additional arguments to uniroot
##' @return The uniroot output.
##' @examples
##' stopool <- function(p) pnorm(sum(qnorm(p, lower.tail = FALSE))/
##' sqrt(length(p), lower.tail = FALSE)
##' estimatePrb(stopool, 0.05, M = 10)
##' estimatePrb(stopool, 0.05, M = 10, b = 0.5)
##' @author Chris Salahub
estimatePrb <- function(poolFun, alpha = 0.05, b = 1, M = 2,
                        interval = c(0, b), poolArgs = list(), ...) {
    ## check upper limit
    atb <- do.call(poolFun, args = c(list(rep(b, M)), poolArgs))
    if (atb < alpha) {
        return(b)
    } else { # if it's greater than alpha, a root exists
        pcr <- function(x) do.call(poolFun,
                                   args = c(list(c(x, rep(b, M-1))),
                                            poolArgs)) - alpha
        rt <- try(uniroot(pcr, interval = interval, ...)$root,
                  silent = TRUE) # error handling
        if ((is(rt, "try-error")|is(rt, "error")) &
            interval[1] < 2*.Machine$double.eps &
            interval[2] > 1 - 2*.Machine$double.eps) { # basically 0 and 1
            if (grepl("values at end points not of opposite sign",
                      rt)) {
                0 # handle floating point issues without errors
            } else {
                stop(rt)
            }
        } else {
            rt
        }
    }
}

##' @title Compute the centrality quotient
##' @description Estimates the centrality quotient for an arbitrary
##' pooled p-value function.
##' @details The centrality quotient communicates the tendency for a
##' test to favour evidence shared among all tests over strong
##' evidence in a single test.
##'
##' This function uses the individual estimation functions for
##' central and marginal rejection levels to compute the centrality
##' quotient for an arbitrary pooled p-value function. The option to
##' specify b for marginal rejection is included in case the pooled p
##' -value has strange behaviour when p-values are equal to 1.
##' @param poolFun function accepting a vector of p-values
##' @param alpha numeric between 0 and 1
##' @param M integer, how many p-values are there?
##' @param interval two numerics giving the bounds of root-searching
##' @param poolArgs (optional) additional named arguments for poolFun
##' @param ... additional arguments to uniroot
##' @return The uniroot output.
##' @examples
##' stopool <- function(p) pnorm(sum(qnorm(p, lower.tail = FALSE))/
##' sqrt(length(p), lower.tail = FALSE)
##' estimatePrb(stopool, 0.05, M = 10)
##' estimatePrb(stopool, 0.05, M = 10, b = 0.5)
##' @author Chris Salahub
estimateQ <- function(poolFun, alpha = 0.05, M = 2,
                      interval = c(0,1), poolArgs = list(), ...) {
    pc <- estimatePc(poolFun, alpha = alpha, M = M,
                     poolArgs = poolArgs, interval = interval, ...)
    pr <- estimatePrb(poolFun, alpha = alpha, M = M, b = 1,
                     poolArgs = poolArgs, interval = interval, ...)
    (pc - pr)/pc
}
