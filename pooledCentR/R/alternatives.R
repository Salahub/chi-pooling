##' @title Compute the Kullback-Leibler divergence
##' @description Computes the Kullback-Leibler divergence for two
##' arbitrary densities f1 and f2.
##' @details Given lower and upper bounds, this function integrates
##' the expression for the Kullback-Leibler divergence KL(f1|f2).
##' @param f1 density function of a real-valued random variable
##' @param f2 density function of a real-values random variable
##' @param lower real value, the lower bound of integration
##' @param upper real value, the upper bound of integration
##' @return A real value.
##' @examples
##' klDiv(dunif, function(x) dbeta(x, 0.5, 1))
##' @author Chris Salahub
klDiv <- function(f1, f2, lower = 0, upper = 1) {
    div <- function(x) {
        vals <- f1(x)*log(f1(x)/f2(x))
        vals[f1(x) == 0] <- 0
        vals
    }
    integrate(div, lower = lower, upper = upper)
}

##' @title Compute the Kullback-Leibler divergence between the beta
##' and uniform distributions
##' @description Computes the Kullback-Leibler divergence for the
##' special case of the uniform density against the beta density.
##' @details This function accepts either the a/b parameterization
##' (equivalent to shape1/shape2 in R), or the a/w parameterization
##' which links the divergence to the UMP test.
##' @param a first shape parameter between 0 and infinity
##' @param w UMP parameter between 0 and 1
##' @param b second shape parameter between 0 and infinity
##' @return A real value.
##' @examples
##' betaDiv(a = 0.5, w = 0.5)
##' betaDiv(a = 0.1, b = 1)
##' @author Chris Salahub
betaDiv <- function(a, w = (1 - a)/(b - a),
                    b = 1/w + a*(1 - 1/w)) {
    lbeta(a, b) + a + b - 2 # compute
}

##' @title Estimate parameter for a given beta KL divergence and UMP
##' test
##' @description Computes the first parameter value for a given KL
##' divergence and UMP test.
##' @details This function uses uniroot to invert the beta divergence
##' for a given w and return the a value which gives that beta
##' divergence given the UMP parameter w. The search interval is
##' specified internally, so should not be passed in using additional
##' argument.
##' @param w UMP parameter between 0 and 1
##' @param logd numeric value, the log KL divergence
##' @param ... additional arguments to uniroot
##' @return A real value.
##' @examples
##' findA(0.5, logd = 0)
##' @author Chris Salahub
findA <- function(w, logd = 0, ...) {
    scr <- function(a) log(betaDiv(a, w = w)) - logd
    uniroot(scr, interval = c(0,1), ...)$root
}

##' Alternatives
##' @title Generate realizations of beta alternative distributions
##' @description These functions can be used to generate samples of
##' p-values all following a beta distribution (H4) or following
##' either uniform or beta distributions according to proportion
##' eta (H3).
##' @details These functions are provided as a convenience, and
##' support a/b (shape1/shape2) or a/w specification of beta
##' parameters.
##' @param a first beta parameter, numeric between 0 and infinity
##' @param b second beta parameter, numeric  between 0 and infinity
##' @param w UMP parameter between 0 and 1
##' @param eta numeric between 0 and 1, proportion of non-null
##' tests per sample
##' @param M number of p-values per realization
##' @param N number of realizations
##' @return An N by M matrix of simulated p-values.
##' @examples
##' pBetaH4(a = 0.5, b = 1.5, M = 10, N = 100)
##' pBetaH3(a = 0.5, b = 1.5, eta = 0.5, M = 10, N = 100)
##' @author Chris Salahub
##' @describeIn alternatives iid Beta(a,w) p-values
pBetaH4 <- function(a, b = 1/w + a*(1 - 1/w), w = (1 - a)/(b - a),
                    M = 2, N = 10) {
    ps <- rbeta(M*N, shape1 = a, shape2 = b)
    matrix(ps, ncol = M)
}
##' @describeIn alternatives M*eta iid Beta(a,w) p-values, others uniform
pBetaH3 <- function(a, b = 1/w + a*(1 - 1/w), w = (1 - a)/(b - a),
                    eta = 0.5, M = 2, N = 10) {
    m1 <- floor(M*eta) # number of non-null tests to generate
    h1 <- rbeta(m1*N, shape1 = a, shape2 = b)
    h0 <- runif((M - m1)*N) # true nulls
    cbind(matrix(h1, ncol = m1, nrow = nsim),
          matrix(h0, ncol = M - m1, nrow = nsim))
}
