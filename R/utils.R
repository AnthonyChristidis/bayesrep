# Internal helpers shared by the Gaussian and binomial MCMC engines.
# None of these are exported.

#' Initialise a streaming pointwise log-likelihood accumulator
#'
#' WAIC needs, for every observation, the posterior mean of the likelihood and
#' the posterior variance of the log-likelihood. Both are accumulated online so
#' that no n-by-draws matrix is ever materialised.
#'
#' @param n Integer. Number of observations.
#' @return A list of accumulator state.
#' @noRd
newLogLikAcc <- function(n) {
    list(lse = NULL, sum = numeric(n), sumSq = numeric(n), count = 0L)
}

#' Fold one draw's pointwise log-likelihood into the accumulator
#'
#' @param acc Accumulator from \code{newLogLikAcc}.
#' @param ll Numeric vector of pointwise log-likelihoods for the current draw.
#' @return The updated accumulator.
#' @noRd
accumulateLogLik <- function(acc, ll) {
    if (is.null(acc$lse)) {
        acc$lse <- ll
    } else {
        mx <- pmax(acc$lse, ll)
        acc$lse <- mx + log(exp(acc$lse - mx) + exp(ll - mx))
    }
    acc$sum   <- acc$sum + ll
    acc$sumSq <- acc$sumSq + ll^2
    acc$count <- acc$count + 1L
    acc
}

#' Convert a log-likelihood accumulator into WAIC
#'
#' @param acc Accumulator from \code{newLogLikAcc}.
#' @return A list with \code{waic}, \code{lppd}, \code{pWaic} and \code{elpd}.
#' @noRd
finalizeWaic <- function(acc) {
    tDraws <- acc$count
    if (tDraws < 2L) {
        return(list(waic = NA_real_, lppd = NA_real_, pWaic = NA_real_, elpd = NA_real_))
    }
    lppd <- sum(acc$lse - log(tDraws))
    meanLl <- acc$sum / tDraws
    varLl <- (acc$sumSq / tDraws - meanLl^2) * tDraws / (tDraws - 1L)
    pWaic <- sum(pmax(varLl, 0))
    list(waic = -2 * (lppd - pWaic), lppd = lppd, pWaic = pWaic, elpd = lppd - pWaic)
}

#' Draw Polya-Gamma latent weights
#'
#' At the default learning rate of one the required shape is exactly one, where
#' the exact Devroye sampler in \pkg{pgdraw} applies and is much faster. Any
#' other learning rate needs a non-integer shape, which falls back to the
#' truncated series representation.
#'
#' @param shape Numeric. The Polya-Gamma shape, equal to the learning rate.
#' @param psi Numeric vector of linear predictors.
#' @param trunc Integer. Series truncation for the non-integer fallback.
#' @return A numeric vector of latent weights.
#' @noRd
drawPolyaGamma <- function(shape, psi, trunc) {
    if (isTRUE(all.equal(shape, 1))) {
        pgdraw::pgdraw(1, psi)
    } else {
        rpgVecCpp(shape, psi, trunc)
    }
}

#' Numerically stable Bernoulli log-likelihood from a linear predictor
#'
#' @param y Numeric 0/1 vector.
#' @param psi Numeric vector of linear predictors.
#' @return Pointwise log-likelihoods.
#' @noRd
bernoulliLogLik <- function(y, psi) {
    y * psi - (pmax(psi, 0) + log1p(exp(-abs(psi))))
}

#' Robbins-Monro adaptation of a random-walk proposal scale
#'
#' Applied during burn-in only, so the post-burn-in chain is a valid
#' Metropolis-Hastings chain with a fixed proposal.
#'
#' @param logSd Current log proposal standard deviation.
#' @param accepted 0/1 indicator for the latest proposal.
#' @param iteration Current iteration index.
#' @param target Target acceptance rate. Default 0.44, optimal for scalar
#'   random-walk Metropolis.
#' @return The updated log proposal standard deviation.
#' @noRd
adaptProposalSd <- function(logSd, accepted, iteration, target = 0.44) {
    logSd + (accepted - target) / sqrt(iteration)
}

#' Allocate the retained-draw storage shared by both MCMC engines
#'
#' @param p,numModels,nKeep Dimensions of the problem.
#' @param storeChains Logical. Allocate the per-sub-model chains as well.
#' @return A list of pre-allocated storage objects.
#' @noRd
newChainStorage <- function(p, numModels, nKeep, storeChains) {
    list(
        ensembleBetaChain = matrix(0, nrow = p, ncol = nKeep),
        gammaChain        = raw(p * numModels * nKeep),
        interceptChain    = numeric(nKeep),
        lambda2Chain      = numeric(nKeep),
        modelSize         = matrix(0L, nrow = numModels, ncol = nKeep),
        pipSum            = matrix(0, nrow = p, ncol = numModels),
        emipCount         = numeric(p),
        thetaSum          = numeric(p),
        thetaSqSum        = numeric(p),
        betaChain         = if (storeChains) array(0, dim = c(p, numModels, nKeep)) else NULL,
        thetaChain        = if (storeChains) matrix(0, nrow = p, ncol = nKeep) else NULL
    )
}

#' Report MCMC progress at roughly ten evenly spaced points
#'
#' @param i Current iteration.
#' @param iter Total iterations.
#' @param burnin Burn-in length.
#' @param label Character tag identifying the engine.
#' @return Invisibly \code{NULL}; called for the message side effect.
#' @noRd
reportProgress <- function(i, iter, burnin, label) {
    step <- max(1L, floor(iter / 10))
    if (i %% step == 0L) {
        message(sprintf("  %s: iteration %d / %d%s", label, i, iter,
                        if (i <= burnin) " (burn-in)" else ""))
    }
    invisible(NULL)
}
