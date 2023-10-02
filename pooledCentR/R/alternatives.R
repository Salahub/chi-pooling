## the Kullback-Leibler divergence
klDiv <- function(f1, f2, lower = 0, upper = 1) {
    div <- function(x) {
        vals <- f1(x)*log(f1(x)/f2(x))
        vals[f1(x) == 0] <- 0
        vals
    }
    integrate(div, lower = lower, upper = upper)
}

## the KL div for the beta case in particular
betaDiv <- function(a, w) {
    b <- bwa(w, a) # second parameter
    lbeta(a, b) + a + b - 2 # compute
}

## find the a that gives a particular log beta divergence given w
findA <- function(w, logd = 0, ...) {
    scr <- function(a) log(betaDiv(a, w)) - logd
    uniroot(scr, c(0,1), ...)$root
}

## compute b given the ratio w and a
bwa <- function(w, a) 1/w + (1 - 1/w)*a

## HR beta generation functions
pBetaH4 <- function(a, w, M, nsim) {
    b <- bwa(w, a)
    ps <- rbeta(M*nsim, shape1 = a, shape2 = b)
    matrix(ps, ncol = M)
}

## m1 are not uniform, M - m1 are
pBetaH3 <- function(a, w, m1) {
    b <- bwa(w, a)
    function(M, nsim) {
        h1 <- rbeta(m1*nsim, shape1 = a, shape2 = b)
        h0 <- runif((M - m1)*nsim) # true nulls
        cbind(matrix(h1, ncol = m1, nrow = nsim),
              matrix(h0, ncol = M - m1, nrow = nsim))
    }
}
