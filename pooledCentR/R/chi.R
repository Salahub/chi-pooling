##' @title Chi-squared p-value pooling
##' @description This implements the chi-squared pooled p-value
##' which can be used to control the centrality quotient when pooling
##' p-values.
##' @details The chi-squared pooled p-value is a quantile
##' transformation pooled p-value based on the chi-squared
##' distribution with degrees of freedom kappa. By setting kappa
##' between 0 and infinity, smooth interpolation is achieved between
##' Tippett's minimum pooled p-value and Stouffer's normal quantile
##' pooled p-value respectively. Choosing a kappa value of 2, Fisher's
##' pooling function is obtained. Tippett's pooled p-value is
##' maximally non-central and Stouffer's is maximally central, while
##' Fisher's presents a balance between marginal and central
##' rejection.
##' @param p numeric vector of p-values between 0 and 1
##' @param kappa numeric value between 0 and infinity
##' @return A pooled p-value between 0 and 1.
##' @examples
##' p <- c(0.1, 0.5, 0.9)
##' chipool(p, exp(-4))
##' chipool(p, 2)
##' chipool(p, exp(4))
##' @author Chris Salahub
chiPool <- function(p, kappa) {
    M <- length(p) # dimension
    pchisq(sum(qchisq(p, df = kappa, lower.tail = FALSE)),
           df = M*kappa, lower.tail = FALSE)
}

##' @title Chi-squared central rejection level
##' @description Computes the central rejection level for the
##' chi-squared pooled p-value.
##' @details The central rejection level is the maximum p-value
##' shared among all tests which still results in rejection of the
##' null using a pooled p-value. For the chi-squared pooled p-value,
##' this is an upper tail probability of the chi-squared distribution.
##' This function computes the upper tail probability for a given
##' sample size M, degrees of freedom kappa, and rejection level
##' alpha.
##' @param kappa numeric between 0 and infinity
##' @param M integer sample size greater than 0
##' @param alpha numeric between 0 and 1
##' @return A numeric between 0 and 1.
##' @examples
##' chiPc(2, 10, 0.05)
##' chiPc(2, 20, 0.05) # increases in sample size
##' @author Chris Salahub
chiPc <- function(kappa, M, alpha = 0.05) {
    pchisq(qchisq(alpha, df = M*kappa, lower.tail = FALSE)/M,
           df = kappa, lower.tail = FALSE)
}

##' @title Chi-squared marginal rejection level
##' @description Computes the marginal rejection level for the
##' chi-squared pooled p-value.
##' @details The marginal rejection level is the maximum p-value
##' in a single test which results in rejection when all other tests
##' produce p-values of one. For the chi-squared pooled p-value,
##' this is an upper tail probability of the chi-squared distribution.
##' This function computes the upper tail probability for a given
##' sample size M, degrees of freedom kappa, and rejection level
##' alpha.
##' @param kappa numeric between 0 and infinity
##' @param M integer sample size greater than 0
##' @param alpha numeric between 0 and 1
##' @return A numeric between 0 and 1.
##' @examples
##' chiPr(2, 10, 0.05)
##' chiPr(2, 20, 0.05)
##' @author Chris Salahub
chiPr <- function(kappa, M, alpha = 0.05) {
    pchisq(qchisq(alpha, df = M*kappa, lower.tail = FALSE),
           df = kappa, lower.tail = FALSE)
}

##' @title Chi-squared centrality quotient
##' @description Computes the centrality quotient of the chi-square
##' pooled p-value.
##' @details The centrality quotient of a pooled p-value measures the
##' relative preference it gives to p-values all sharing the same
##' level of evidence over a single test with strong evidence relative
##' to others. For the chi-square pooled p-value, this is a
##' conditional probability which this function computes.
##' @param kappa numeric between 0 and infinity
##' @param M integer sample size greater than 0
##' @param alpha numeric between 0 and 1
##' @return A numeric between 0 and 1.
##' @examples
##' chiCentQuot(2, 10, 0.05)
##' chiCentQuot(2, 20, 0.05)
##' chiCentQuot(0.5, 20, 0.05)
##' @author Chris Salahub
chiCentQuot <- function(kappa, M, alpha = 0.05) {
    1 - chiPr(kappa, M, alpha)/chiPc(kappa, M, alpha)
}

##' @title Chi-squared kappa for a given centrality quotient
##' @description Computes the kappa (degrees of freedom) required
##' to obtain a given centrality quotient using the chi-square pooled
##' p-value.
##' @details This function is essentially a wrapper for uniroot which
##' finds where chiCentQuot gives an output equal to the given
##' centrality quotient to provide an approximate kappa giving that
##' quotient.
##' @param cq numeric between 0 and 1
##' @param M integer sample size greater than 0
##' @param alpha numeric between 0 and 1
##' @param interval numeric of length 2, where should roots be sought?
##' @param tol numeric, how close do values need to be for equality?
##' @return A numeric within interval.
##' @examples
##' chiKappa(0.5, 10, 0.05)
##' chiKappa(0.5, 20, 0.05)
##' chiKappa(0.5, 100, 0.05, interval = c(10, 1000))
##' @author Chris Salahub
chiKappa <- function(cq, M, alpha = 0.05, interval = c(0,100),
                     tol = .Machine$double.eps^0.5) {
    uniroot(function(x) chiCentQuot(x, M, alpha) - cq,
            interval = interval, tol = tol)$root
}
