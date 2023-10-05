library(harmonicmeanp)
library(geometry)

## this script provides basic implementations of previous pooled
## p-value functions and functions which can be used to plot power
## curves for the bivariate case for each

## pooling methods
Bon <- function(p) min(1, length(p)*min(p))

Tip <- function(p) 1 - (1 - min(p))^length(p)

Fis <- function(p) 1 - pchisq(-2*sum(log(p)), 2*length(p))

Cin <- function(p) 1 - pchisq(sum(qchisq(1 - p, df = 1)), length(p))

Chi <- function(p, shape = 1) {
    1 - pchisq(sum(qchisq(1 - p, df = shape)), shape*length(p))
}

Sto <- function(p) 1 - pnorm(1/sqrt(length(p))*sum(qnorm(1 - p)))

Wil <- function(p, each = 0.05) 1 - pbinom(sum(p <= each),
                                           length(p), each)

Ord <- function(p, k = 5) {
    pk <- sort(p)[k]
    1 - pbinom(k - 1, length(p), pk)
}

Exp <- function(p, lambda = 1) {
    1 - pgamma(sum(qexp(1 - p, rate = lambda)),
               shape = length(p), rate = lambda)
}

Log <- function(p) {
    M <- length(p)
    scl <- pi*sqrt((M*(5*M+2))/(3*(5*M+4)))
    pt(sum(log(p/(1-p)))/scl, df = 5*M+4)
}

Edg <- function(p) {
    M <- length(p)
    S <- sum(p)
    inds <- 0:(floor(S))
    sum((-1)^inds*choose(M, inds)*(S - inds)^M/(factorial(M)))
}

HMP <- function(p) length(p)/sum(1/p)

HMP2 <- function(p) p.hmp(p, L = length(p))

## try all methods out on the test grid
poolMeths <- list(bonferroni = Bon, tippett = Tip,
                  ord2 = function(p) Ord(p, k = 2),
                  ord10 = function(p) Ord(p, k = 10),
                  fisher = Fis, cinvie = Cin, stouffer = Sto,
                  wilk = Wil, mudgeo = Log, edgington = Edg,
                  hmp = HMP, hmpex = HMP2)

## set up a palette for each
pal <- hcl.colors(length(poolMeths), palette = "Dark3")
names(pal) <- c("fisher", "wilk", "cinvie", "bonferroni", "hmp",
                "hmpex", "tippett", "ord2", "ord10", "stouffer",
                "edgington", "mudgeo")

## plotting these as parametric curves for different alpha instead
Boncrv <- function(x, alpha = 0.05) {
    bd <- alpha/2
    ret <- x
    ret[x <= bd] <- 1
    ret[x > bd] <- bd
    ret
}

Tipcrv <- function(x, alpha = 0.05) {
    bd <- 1 - (1 - alpha)^0.5
    ret <- x
    ret[x <= bd] <- 1
    ret[x > bd] <- bd
    ret
}

Fiscrv <- function(x, alpha = 0.05) {
    vals <- (1/x)*exp(-1/2*qchisq(1 - alpha, df = 4))
    vals[vals > 1] <- 1
    vals
}

Stocrv <- function(x, alpha = 0.05) {
    1 - pnorm(sqrt(2)*qnorm(1 - alpha) - qnorm(1 - x))
}

Cincrv <- function(x, alpha = 0.05) {
    1 - pchisq(qchisq(1 - alpha, df = 2) -
               qchisq(1 - x, df = 1), df = 1)
}

Chicrv <- function(x, alpha = 0.05, df = 1) {
    1 - pchisq(qchisq(1 - alpha, df = 2*df) -
               qchisq(1 - x, df = df), df = df)
}

Wilcrv <- function(x, alpha = 0.05, ae = 0.05) {
    ps <- dbinom(c(0,1,2), 2, ae)
    needed <- sum(1 - ps < 1 - alpha)
    xtest <- x <= ae
    ret <- rep(1, length(x))
    if (needed == 3) {
        ret <- rep(0, length(x))
    } else if (needed == 2) {
        ret[xtest] <- ae
        ret[!xtest] <- 0
    } else if (needed == 1) {
        ret[xtest] <- 1
        ret[!xtest] <- ae
    }
    ret
}

MGcrv <- function(x, alpha = 0.05) {
    c <- pi*sqrt(24)/sqrt(42)
    a <- c*qt(alpha, df = 14) - log(x/(1-x))
    vals <- exp(a)/(1 + exp(a))
    vals[x == 0] <- 1 # correct for x = 0
    vals
}

Edgcrv <- function(x, alpha = 0.05) {
    vals <- sqrt(2*alpha) - x
    vals[vals < 0] <- 0
    vals
}

HMPcrv <- function(x, alpha = 0.05) {
    vals <- 1/(2/alpha - 1/x)
    vals[x <= alpha/2] <- 1
    vals[vals > 1] <- 1
    vals
}

HMPexcrv <- function(x, alpha = 0.05) {
    qtl <- qLandau(1 - alpha, mu = log(2) + 0.874)
    vals <- 1/(2*qtl - 1/x)
    vals[1/x >= 2*qtl] <- 1
    vals[vals > 1] <- 1
    vals
}

## some nicer names
methsNiceName <- list(bonferroni = Boncrv, tippett = Tipcrv,
                      fisher = Fiscrv, stouffer = Stocrv,
                      cinvie = Cincrv, wilk = Wilcrv,
                      edgington = Edgcrv, mudgeo = MGcrv,
                      hmp = HMPcrv, hmpex = HMPexcrv)
