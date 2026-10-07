#' Gibbs sampler for Gaussian Bayesian Repulsive Ensembles
#'
#' Internal MCMC engine behind \code{\link{bre}} for continuous responses.
#'
#' Under the composite likelihood every sub-model is fitted to the full response
#' rather than to a partial residual, and the sub-models are conditionally
#' independent given \eqn{\sigma^2}, coupled only through the Ising prior. The
#' \eqn{1/G} tempering enters as an inflated error variance \eqn{G\sigma^2} in
#' each sub-model's conditional, which keeps the total information in the
#' posterior equal to \eqn{n} observations rather than \eqn{Gn}.
#'
#' @param y Numeric vector (length n). Centred continuous response.
#' @param X Numeric matrix (n x p). Standardised design matrix.
#' @param numModels Integer. Number of sub-models G.
#' @param iter,burnin,thin Integers controlling chain length and storage.
#' @param tauSq Numeric. Slab variance for a sub-model coefficient.
#' @param lambda2 Numeric. Starting (or fixed) repulsion strength.
#' @param learnLambda2 Logical. Sample \code{lambda2} rather than hold it fixed.
#' @param lambda2Shape,lambda2Rate Numeric. Gamma prior for \code{lambda2}.
#' @param thetaShape1,thetaShape2 Numeric. Beta hyperprior for theta.
#' @param nu0,sigma0Sq Numeric. Inverse-Gamma prior for \eqn{\sigma^2}.
#' @param learningRate Numeric. Composite-likelihood learning rate.
#' @param transferMoves Integer. Feature-relocation proposals per iteration.
#' @param thetaSweeps Integer. Metropolis sweeps per iteration for theta.
#' @param storeChains Logical. Retain per-sub-model chains.
#' @param verbose Logical. Report progress.
#'
#' @return A list of retained chains, posterior summaries and diagnostics.
#'
#' @keywords internal
runMcmcContinuous <- function(y, X, numModels, iter, burnin, thin,
                              tauSq, lambda2, learnLambda2,
                              lambda2Shape, lambda2Rate,
                              thetaShape1, thetaShape2,
                              nu0, sigma0Sq, learningRate,
                              transferMoves, thetaSweeps,
                              storeChains, verbose) {

    n <- nrow(X)
    p <- ncol(X)

    xTx <- colSums(X^2)

    betaMatrix  <- matrix(0, nrow = p, ncol = numModels)
    gammaMatrix <- matrix(0, nrow = p, ncol = numModels)
    thetaVec    <- rep(thetaShape1 / (thetaShape1 + thetaShape2), p)

    sigmaSq <- stats::var(y)
    if (!is.finite(sigmaSq) || sigmaSq <= 0) sigmaSq <- 1
    intercept <- mean(y)

    logProposalSd <- log(0.5)
    thetaAccepted <- 0; thetaProposed <- 0
    lambdaAccepted <- 0; lambdaProposed <- 0
    transferAccepted <- 0; transferProposed <- 0

    keepAt <- seq.int(burnin + 1L, iter, by = thin)
    isKept <- logical(iter); isKept[keepAt] <- TRUE
    nKeep <- length(keepAt)

    st <- newChainStorage(p, numModels, nKeep, storeChains)
    sigmaSqChain <- numeric(nKeep)
    llAcc <- newLogLikAcc(n)
    blockSize <- p * numModels
    unitW <- matrix(1, nrow = n, ncol = numModels)
    k <- 0L

    for (i in seq_len(iter)) {

        betaBar <- rowSums(betaMatrix) / numModels

        # 1. Error variance, from the mean sub-model residual sum of squares.
        sigmaSq <- updateSigmaSqCpp(y = y, X = X, betaMatrix = betaMatrix,
                                    intercept = intercept,
                                    learningRate = learningRate, nu0 = nu0,
                                    sigma0Sq = sigma0Sq)

        # 2. Intercept. It is shared by all sub-models, and its full conditional
        #    under the tempered composite likelihood is the familiar one.
        intercept <- stats::rnorm(1, mean(y - as.numeric(X %*% betaBar)),
                                  sqrt(sigmaSq / n))

        yAdj <- y - intercept
        sigmaSqScaled <- sigmaSq / learningRate

        # 3. Every sub-model is fitted to the whole response; no partial
        #    residual is formed. Order is randomised so that no sub-model is
        #    systematically served first.
        for (g in sample.int(numModels)) {

            upd <- updateBetaGammaCpp(yScaled = yAdj,
                                      X = X, xTx = xTx,
                                      betaGroup = betaMatrix[, g],
                                      gammaMatrix = gammaMatrix,
                                      targetGroup = g,
                                      sigmaSqScaled = sigmaSqScaled,
                                      tauSqScaled = tauSq,
                                      thetaVec = thetaVec,
                                      lambda2 = lambda2)

            betaMatrix[, g] <- upd$betaGroup
            gammaMatrix <- upd$gammaMatrix
        }

        # 4. Relocate features between sub-models in single accepted steps.
        if (transferMoves > 0L && numModels > 1L) {
            tr <- transferFeatureCpp(resp = matrix(yAdj, nrow = n, ncol = numModels),
                                     sqrtW = unitW,
                                     X = X,
                                     betaMatrix = betaMatrix,
                                     gammaMatrix = gammaMatrix,
                                     sigmaSqScaled = sigmaSqScaled,
                                     tauSq = tauSq,
                                     numMoves = transferMoves)
            betaMatrix  <- tr$betaMatrix
            gammaMatrix <- tr$gammaMatrix
            transferAccepted <- transferAccepted + tr$numAccepted
            transferProposed <- transferProposed + tr$numAttempted
        }

        mCounts <- rowSums(gammaMatrix)

        # 5. Baseline inclusion probabilities, exact up to the Ising constant.
        thUpd <- updateThetaCpp(thetaVec = thetaVec, mCounts = mCounts,
                                aPrior = thetaShape1, bPrior = thetaShape2,
                                lambda2 = lambda2, numModels = numModels,
                                numSweeps = thetaSweeps)
        thetaVec <- as.numeric(thUpd$thetaVec)
        thetaAccepted <- thetaAccepted + thUpd$numAccepted
        thetaProposed <- thetaProposed + thUpd$numProposed

        # 6. Repulsion strength.
        if (learnLambda2) {
            lamUpd <- updateLambda2Cpp(lambda2 = lambda2, thetaVec = thetaVec,
                                       mCounts = mCounts, numModels = numModels,
                                       priorShape = lambda2Shape,
                                       priorRate = lambda2Rate,
                                       proposalSd = exp(logProposalSd))
            lambda2 <- lamUpd$lambda2
            lambdaAccepted <- lambdaAccepted + lamUpd$accepted
            lambdaProposed <- lambdaProposed + 1L
            if (i <= burnin) {
                logProposalSd <- adaptProposalSd(logProposalSd, lamUpd$accepted, i)
            }
        }

        # 7. Storage.
        if (isKept[i]) {
            k <- k + 1L
            betaBar <- rowSums(betaMatrix) / numModels

            st$ensembleBetaChain[, k] <- betaBar
            st$gammaChain[((k - 1L) * blockSize + 1L):(k * blockSize)] <- as.raw(gammaMatrix)
            st$interceptChain[k] <- intercept
            st$lambda2Chain[k]   <- lambda2
            st$modelSize[, k]    <- as.integer(colSums(gammaMatrix))
            st$pipSum    <- st$pipSum + gammaMatrix
            st$emipCount <- st$emipCount + (mCounts > 0)
            st$thetaSum  <- st$thetaSum + thetaVec
            st$thetaSqSum <- st$thetaSqSum + thetaVec^2
            sigmaSqChain[k] <- sigmaSq

            if (storeChains) {
                st$betaChain[, , k] <- betaMatrix
                st$thetaChain[, k]  <- thetaVec
            }

            # Pointwise composite log-density: the 1/G-weighted average of the
            # sub-model log-densities, which sums to the composite objective.
            resid <- yAdj - X %*% betaMatrix
            ll <- numModels * learningRate * (-0.5 * log(2 * pi * sigmaSq)) -
                learningRate * rowSums(resid^2) / (2 * sigmaSq)
            llAcc <- accumulateLogLik(llAcc, as.numeric(ll))
        }

        if (verbose) reportProgress(i, iter, burnin, "gaussian")
    }

    thetaMean <- st$thetaSum / nKeep
    thetaVar  <- pmax(st$thetaSqSum / nKeep - thetaMean^2, 0)

    list(
        ensembleBetaChain = st$ensembleBetaChain,
        gammaChain        = st$gammaChain,
        betaChain         = st$betaChain,
        thetaChain        = st$thetaChain,
        interceptChain    = st$interceptChain,
        sigmaSqChain      = sigmaSqChain,
        lambda2Chain      = st$lambda2Chain,
        modelSize         = st$modelSize,
        pip               = st$pipSum / nKeep,
        emip              = st$emipCount / nKeep,
        thetaMean         = thetaMean,
        thetaSd           = sqrt(thetaVar),
        waic              = finalizeWaic(llAcc),
        acceptance        = list(
            theta    = thetaAccepted / max(thetaProposed, 1L),
            lambda2  = if (learnLambda2) lambdaAccepted / max(lambdaProposed, 1L) else NA_real_,
            transfer = if (transferProposed > 0L) transferAccepted / transferProposed else NA_real_,
            lambda2ProposalSd = exp(logProposalSd)
        )
    )
}
