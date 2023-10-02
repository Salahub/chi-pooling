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
##' number of null cases to give an empirical pooled p-value.
##' @param p numeric vector of p-values between 0 and 1
##' @param w numeric value between 0 and 1
##' @param nsim integer, the number of simulated null cases generated
##' @return A numeric value between 0 and 1.
##' @examples
##' p <- c(0.1, 0.5, 0.9)
##' hrPool(p, 0.2)
##' hrPool(p, 0.5)
##' hrPool(p, 0.9)
##' @author Chris Salahub
hrPool <- function(p, w = 1, nsim = 1e5) {
    M <- length(p) # get dimension
    dat <- matrix(runif(M*nsim), ncol = M)
    pools <- apply(dat, 1, hrStat, w = w)
    mean(hrStat(p, w) >= pools) # observed quantile
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
##' @param kappa numeric between 0 and infinity
##' @param M integer sample size greater than 0
##' @param alpha numeric between 0 and 1
##' @param nsim integer, the number of simulated null cases generated
##' @return A numeric between 0 and 1.
##' @examples
##' hrPc(2, 10, 0.05)
##' hrPc(2, 20, 0.05)
##' @author Chris Salahub
hrPc <- function(w, alpha = 0.05, M = 2, nsim = 1e5) {
    dat <- matrix(runif(M*nsim), ncol = M)
    pools <- apply(dat, 1, hrStat, w = w)
    tempPc <- function(pc, alpha = alpha) {
        mean(hrStat(rep(pc, M), w) >= pools) - alpha
    }
    uniroot(tempPc, interval = c(0, 1-1/(2*nsim)),
            tol = 1/(2*nsim))
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
##' @param kappa numeric between 0 and infinity
##' @param M integer sample size greater than 0
##' @param alpha numeric between 0 and 1
##' @param nsim integer, the number of simulated null cases generated
##' @return A numeric between 0 and 1.
##' @examples
##' hrPr(2, 10, 0.05)
##' hrPr(2, 20, 0.05) # decreases in sample size
##' @author Chris Salahub
hrPr <- function(w, alpha = 0.05, M = 2, nsim = 1e5) {
    dat <- matrix(runif(M*nsim), ncol = M)
    pools <- apply(dat, 1, hrStat, w = w)
    tempPr <- function(pr, alpha = alpha) {
        mean(hrStat(c(pr, rep(0.999, M-1)), w) >= pools) -
            alpha
    }
    uniroot(tempPr, interval = c(0, 1-1/(2*nsim)),
            tol = 1/(2*nsim))
}
