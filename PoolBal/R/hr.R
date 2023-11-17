##' @title UMP beta p-value pooled statistic
##' @description Computes the UMP p-value pooling statistic for a
##' restricted beta family.
##' @details To test the null hypotheses that all p-values are uniform
##' against a restricted beta family 0 < a <= 1 <= b, the most
##' powerful pooled p-value linearly combines upper and lower tail
##' probabilities of the chi-squared distribution with two degrees
##' of freedom with weights w and (1 - w) where w = (1 - a)/(b - a).
##'
##' This function computes the statistic given by this combination
##' for a collection of p-values, simulation or approximation is
##' required to convert this to a p-value.
##' @param p numeric vector of p-values between 0 and 1
##' @param w numeric value between 0 and 1
##' @return A numeric value giving the pooled statistic.
##' @examples
##' p <- c(0.1, 0.5, 0.9)
##' hrStat(p, 0.2)
##' hrStat(p, 0.5)
##' hrStat(p, 0.9)
##' @author Chris Salahub
hrStat <- function(p, w = 1) {
    w*sum(log(p)) - (1 - w)*sum(log(1 - p))
}

##' @title Empirical UMP beta pooled p-value
##' @description Uses simulation under the null to approximate the UMP
##' pooled p-value for a restricted beta family.
##' @details To test the null hypotheses that all p-values are uniform
##' against a restricted beta family 0 < a <= 1 <= b, the most
##' powerful pooled p-value linearly combines upper and lower tail
##' probabilities of the chi-squared distribution with two degrees
##' of freedom with weights w and (1 - w) where w = (1 - a)/(b - a).
##'
##' This function computes the statistic given by this combination
##' for a collection of p-values, and then simulates a specified
##' number of null cases to give an empirical pooled p-value. It
##' produces a closure so that the time-intensive simulation step
##' doesn't need to be repeated.
##' @param w numeric value between 0 and 1
##' @param nsim integer, the number of simulated null cases generated
##' @param M integer, the number of tests to pool
##' @return A closure which accepts a vector of values between 0 and 1
##' and returns a single numeric between 0 and 1
##' @examples
##' p <- c(0.1, 0.5, 0.9)
##' hr2 <- hrpool(w = 0.2, M = 3)
##' hr2(p)
##' hr5 <- hrpool(w = 0.5, M = 3, nsim = 100)
##' hr5(p)
##' @author Chris Salahub
hrPool <- function(w = 1, M = 10, nsim = 1e5) {
    dat <- matrix(runif(M*nsim), ncol = M)
    pools <- apply(dat, 1, hrStat, w = w)
    fun <- function(p) {
        mean(hrStat(p, w) >= pools) # observed quantile
    }
    fun
}

##' @title Empirical UMP beta central rejection level
##' @description Uses simulation to estimate the central rejection
##' level for the UMP pooled p-value of a restricted beta family
##' @details The central rejection level is the maximum p-value
##' shared among all tests which still results in rejection of the
##' null using a pooled p-value.
##'
##' To test the null hypotheses that all p-values are uniform
##' against a restricted beta family 0 < a <= 1 <= b, the most
##' powerful pooled p-value linearly combines upper and lower tail
##' probabilities of the chi-squared distribution with two degrees
##' of freedom with weights w and (1 - w) where w = (1 - a)/(b - a).
##'
##' This function estimates the central rejection level empirically
##' by simulating a specified number of null cases to give an empirical
##' pooled p-value for the rejection level alpha.
##' @param w numeric between 0 and 1
##' @param alpha numeric between 0 and 1
##' @param M integer sample size greater than 0
##' @param nsim integer, the number of simulated null cases generated
##' @return A numeric between 0 and 1.
##' @examples
##' hrPc(2, 10, 0.05)
##' hrPc(2, 20, 0.05)
##' @author Chris Salahub
hrPc <- function(w, alpha = 0.05, M = 2, nsim = 1e5) {
    dat <- matrix(runif(M*nsim), ncol = M)
    pools <- apply(dat, 1, hrStat, w = w)
    tempPc <- function(pc) {
        mean(hrStat(rep(pc, M), w) >= pools) - alpha
    }
    uniroot(tempPc, interval = c(0, 1-1/(2*nsim)),
            tol = 1/(2*nsim))$root
}

##' @title Empirical UMP beta marginal rejection level
##' @description Uses simulation to estimate the marginal rejection
##' level for the UMP pooled p-value of a restricted beta family
##' @details The marginal rejection level is the maximum p-value
##' in a single tests which still results in rejection of the null
##' when all other tests have a p-value of 1.
##'
##' To test the null hypotheses that all p-values are uniform
##' against a restricted beta family 0 < a <= 1 <= b, the most
##' powerful pooled p-value linearly combines upper and lower tail
##' probabilities of the chi-squared distribution with two degrees
##' of freedom with weights w and (1 - w) where w = (1 - a)/(b - a).
##'
##' This function estimates the marginal rejection level empirically
##' by simulating a specified number of null cases to give an empirical
##' pooled p-value for the rejection level alpha.
##' @param w numeric between 0 and 1
##' @param alpha numeric between 0 and 1
##' @param M integer sample size greater than 0
##' @param nsim integer, the number of simulated null cases generated
##' @return A numeric between 0 and 1.
##' @examples
##' hrPr(2, 10, 0.05)
##' hrPr(2, 20, 0.05) # decreases in sample size
##' @author Chris Salahub
hrPr <- function(w, alpha = 0.05, M = 2, nsim = 1e5) {
    dat <- matrix(runif(M*nsim), ncol = M)
    pools <- apply(dat, 1, hrStat, w = w)
    tempPr <- function(pr) {
        mean(hrStat(c(pr, rep(0.999, M-1)), w) >= pools) -
            alpha
    }
    uniroot(tempPr, interval = c(0, 1-1/(2*nsim)),
            tol = 1/(2*nsim))$root
}

##' @title Empirical UMP beta centrality quotient
##' @description Estimates the centrality quotient for the UMP
##' pooled p-value of a restricted beta family.
##' @details The centrality quotient communicates the tendency for a
##' test to favour evidence shared among all tests over strong
##' evidence in a single test.
##'
##' To test the null hypotheses that all p-values are uniform
##' against a restricted beta family 0 < a <= 1 <= b, the most
##' powerful pooled p-value linearly combines upper and lower tail
##' probabilities of the chi-squared distribution with two degrees
##' of freedom with weights w and (1 - w) where w = (1 - a)/(b - a).
##'
##' This function uses the individual estimation functions for
##' central and marginal rejection levels to compute the centrality
##' quotient for the UMP pooled p-value.
##' @param w numeric between 0 and 1
##' @param alpha numeric between 0 and 1
##' @param M integer sample size greater than 0
##' @param nsim integer, the number of simulated null cases generated
##' @return An empirical estimate of the centrality quotient.
##' @examples
##' estimatePrb(stopool, 0.05, M = 10)
##' estimatePrb(stopool, 0.05, M = 10, b = 0.5)
##' @author Chris Salahub
hrQ <- function(w, alpha = 0.05, M = 2, nsim = 1e5) {
    pc <- hrPc(w, alpha, M, nsim)
    pr <- hrPr(w, alpha, M, nsim)
    1 - pr/pc
}
