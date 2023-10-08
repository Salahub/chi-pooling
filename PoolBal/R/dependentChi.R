##' @title Convert p-value correlation to chi-squared covariance
##' @description Convert a matrix of correlations between p-values to
##' a matrix of covariances between their chi-squared transforms.
##' @details This function uses models fit to large simulated data
##' sets to convert a matrix of correlations between genetic markers
##' the covariance matrix of chi-squared random variables gained from
##' transforming p-values on these markers. The simulations used to
##' create data for these models assume the p-values for each marker
##' arise from tests of association with a common, normally
##' distributed trait independent of all markers. As a result, this
##' conversion function should be used only in analogous settings.
##'
##' Models were fit for degrees of freedom at increments of 0.1
##' between -8 and 8 on the log scale, and interpolation is applied
##' if the degrees of freedom given to the function does not fall
##' exactly on this grid (with a warning provided to the user).
##'
##' If a user wants to generalize this setting, the option to
##' provide a custom list of models which predict based on a named
##' argument `zcor` is supported. Each model must have a name in the
##' list that can be converted to a numeric, and these are assumed
##' to be on the natural log scale.
##' @param sigma M by M correlation matrix between markers
##' @param kappa numeric degrees of freedom
##' @param models model object with a predict method
##' @return M by M matrix of chi-squared covariances
##' @author Chirs Salahub
convertGeneticSigma <- function(sigma, kappa, models = chiCorMods) {
    kaps <- as.numeric(names(models)) # kappa values of models
    kapDifs <- abs(log(kappa) - kaps)
    kapInds <- order(kapDifs)[1:2] # best indices
    mod1 <- matrix(predict(models[[kapInds[1]]], # predict 1st
                           newdata = data.frame(zcor = c(sigma))),
                   ncol = ncol(sigma))
    mod2 <- matrix(predict(models[[kapInd[2]]], # predict 2nd
                           newdata = data.frame(zcor = c(sigma))),
                   ncol = ncol(sigma))
    ## perform linear interpolation on the log scale
    (mod1*kapDifs[kapInds[2]] + mod2*kapDifs[kapInds[1]])/
        (sum(kapDifs[kapInds]))
}

##' @title Satterthwaite p-values
##' @description p-value of the sum of dependent chi-squared using
##' the Satterthwaite approximation for the degrees of freedom.
##' @details Computes the p-value of an observed vector of chi-squared
##' variables using the Satterthwaite approximation. This approximates
##' the sum of dependent chi-squared variables with a scaled
##' chi-squared distribution with degrees of freedom chosen to match
##' the first two moments of the dependent sum.
##' @param qs M numeric values (observed chi-squared values)
##' @param covmat M by M covariance matrix of qs
##' @param kappa degrees of freedom of qs
##' @return a numeric in [0,1], the p-value of the sum
##' @author Chris Salahub
satterApproxP <- function(qs, covmat, kappa) {
    diag(covmat) <- 0 # makes later sums easier
    M <- length(qs) # dimension
    ex <- M*kappa # expectation of sum
    varx <- 2*M*kappa + 2*kappa*sum(covmat) # variance
    f <- 2*ex^2/varx # adjustment factors
    c <- varx/(2*ex)
    pchisq(sum(qs)/c, df = f, lower.tail = FALSE) # p-value
}

##' @title Pool p-values using the Satterthwaite approximation
##' @description Compute the pooled p-value of dependent p-values
##' based on the dependence present when they are all converted
##' to chi-squared random variables by the same chi-squared
##' quantile function.
##' @details Care must be taken in the arguments for this function, as
##' the covmat argument accepts the covariance of the transformed
##' variables rather than the covariance of the p-values, and so
##' passes the argument covmat directly to the function that computes
##' the Satterthwaite approximation. For the case of genetic markers,
##' the `convertGeneticSigma` function provides the appropriate matrix
##' given a genetic correlation matrix.
##' @param ps numeric vector of M p-values
##' @param covmat M by M covariance matrix of chi-squared random
##' variables arising from quantile transformations of ps
##' @param kappa numeric degrees of freedom
##' @return A pooled p-value between 0 and 1.
##' @author Chris Salahub
satterChiPool <- function(ps, covmat, kappa) {
    chis <- qchisq(ps, df = kappa, lower.tail = FALSE)
    satterApproxP(chis, covmat, kappa)
}
