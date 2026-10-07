#' Gibbs sampler for binary Bayesian Repulsive Ensembles
#'
#' Internal MCMC engine behind \code{\link{bre}} for binary responses.
#'
#' Each sub-model carries its own Bernoulli likelihood raised to the power
#' \eqn{1/G}, so the Polya-Gamma augmentation needed is \eqn{PG(1/G, \psi)}
#' rather than \eqn{PG(1, \psi)}. Latent weights are drawn separately for every
#' sub-model, which gives each one its own heteroskedastic pseudo-data, and the
#' intercept is drawn from the full conditional pooled across sub-models.
#'
#' @param y Numeric vector (length n) of zeros and ones.
#' @param X Numeric matrix (n x p). Standardised design matrix.
#' @param numModels Integer. Number of sub-models G.
#' @param iter,burnin,thin Integers controlling chain length and storage.
#' @param tauSq Numeric. Slab variance for a sub-model coefficient.
#' @param lambda2 Numeric. Starting (or fixed) repulsion strength.
#' @param learnLambda2 Logical. Sample \code{lambda2} rather than hold it fixed.
#' @param lambda2Shape,lambda2Rate Numeric. Gamma prior for \code{lambda2}.
#' @param thetaShape1,thetaShape2 Numeric. Beta hyperprior for theta.
#' @param transferMoves Integer. Feature-relocation proposals per iteration.
#' @param thetaSweeps Integer. Metropolis sweeps per iteration for theta.
#' @param pgTrunc Integer. Series truncation for the Polya-Gamma sampler.
#' @param learningRate Numeric. Composite-likelihood learning rate.
#' @param storeChains Logical. Retain per-sub-model chains.
#' @param verbose Logical. Report progress.
#'
#' @return A list of retained chains, posterior summaries and diagnostics.
#'
#' @keywords internal
runMcmcBinary <- function(y, X, numModels, iter, burnin, thin,
                          tauSq, lambda2, learnLambda2,
                          lambda2Shape, lambda2Rate,
                          thetaShape1, thetaShape2,
                          transferMoves, thetaSweeps, pgTrunc, learningRate,
                          storeChains, verbose) {

    n <- nrow(X)
    p <- ncol(X)

    # The 1/G tempering is carried inside the Polya-Gamma shape and the kappa
    # offset, so the working Gaussian likelihood has unit variance.
    sigmaSqScaled <- 1
    pgShape <- learningRate
    kappa <- (y - 0.5) * learningRate

    betaMatrix  <- matrix(0, nrow = p, ncol = numModels)
    gammaMatrix <- matrix(0, nrow = p, ncol = numModels)
    thetaVec    <- rep(thetaShape1 / (thetaShape1 + thetaShape2), p)

    # Start the intercept at the empirical log-odds so the sampler does not have
    # to walk there from zero on an imbalanced sample.
    pBar <- min(max(mean(y), 1 / (n + 1)), 1 - 1 / (n + 1))
    intercept <- log(pBar / (1 - pBar))

    logProposalSd <- log(0.5)
    thetaAccepted <- 0; thetaProposed <- 0
    lambdaAccepted <- 0; lambdaProposed <- 0
    transferAccepted <- 0; transferProposed <- 0

    keepAt <- seq.int(burnin + 1L, iter, by = thin)
    isKept <- logical(iter); isKept[keepAt] <- TRUE
    nKeep <- length(keepAt)

    st <- newChainStorage(p, numModels, nKeep, storeChains)
    llAcc <- newLogLikAcc(n)
    blockSize <- p * numModels
    k <- 0L

    sqrtOmega <- matrix(0, nrow = n, ncol = numModels)
    respMat   <- matrix(0, nrow = n, ncol = numModels)

    for (i in seq_len(iter)) {

        linPred <- X %*% betaMatrix                 # n x G, unweighted
        psi <- intercept + linPred

        # 1. Polya-Gamma weights, one set per sub-model.
        omega <- matrix(drawPolyaGamma(pgShape, as.numeric(psi), pgTrunc),
                        nrow = n, ncol = numModels)
        omega <- pmax(omega, 1e-10)
        sqrtOmega[] <- sqrt(omega)

        # 2. Intercept, pooled across sub-models.
        sumOmega <- sum(omega)
        interceptMean <- sum(kappa - omega * linPred) / sumOmega
        intercept <- stats::rnorm(1, interceptMean, sqrt(1 / sumOmega))

        # 3. Sub-models, each against its own pseudo-data. Forming the weighted
        #    response as kappa / sqrt(omega) avoids a 1 / omega blow-up.
        for (g in sample.int(numModels)) {

            sw <- sqrtOmega[, g]
            zWeighted <- kappa / sw - intercept * sw
            XWeighted <- X * sw
            xTxWeighted <- colSums(XWeighted^2)

            upd <- updateBetaGammaCpp(yScaled = zWeighted,
                                      X = XWeighted, xTx = xTxWeighted,
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

        # 4. Relocate features between sub-models, in each sub-model's own
        #    weighted space.
        if (transferMoves > 0L && numModels > 1L) {
            respMat[] <- kappa / omega - intercept
            tr <- transferFeatureCpp(resp = respMat,
                                     sqrtW = sqrtOmega,
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

        # 5. Baseline inclusion probabilities.
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

            st$ensembleBetaChain[, k] <- rowSums(betaMatrix) / numModels
            st$gammaChain[((k - 1L) * blockSize + 1L):(k * blockSize)] <- as.raw(gammaMatrix)
            st$interceptChain[k] <- intercept
            st$lambda2Chain[k]   <- lambda2
            st$modelSize[, k]    <- as.integer(colSums(gammaMatrix))
            st$pipSum    <- st$pipSum + gammaMatrix
            st$emipCount <- st$emipCount + (mCounts > 0)
            st$thetaSum  <- st$thetaSum + thetaVec
            st$thetaSqSum <- st$thetaSqSum + thetaVec^2

            if (storeChains) {
                st$betaChain[, , k] <- betaMatrix
                st$thetaChain[, k]  <- thetaVec
            }

            # Pointwise composite log-density, averaged over sub-models.
            psiNew <- intercept + X %*% betaMatrix
            llMat <- y * psiNew - (pmax(psiNew, 0) + log1p(exp(-abs(psiNew))))
            llAcc <- accumulateLogLik(llAcc,
                numModels * learningRate * as.numeric(rowMeans(llMat)))
        }

        if (verbose) reportProgress(i, iter, burnin, "binomial")
    }

    thetaMean <- st$thetaSum / nKeep
    thetaVar  <- pmax(st$thetaSqSum / nKeep - thetaMean^2, 0)

    list(
        ensembleBetaChain = st$ensembleBetaChain,
        gammaChain        = st$gammaChain,
        betaChain         = st$betaChain,
        thetaChain        = st$thetaChain,
        interceptChain    = st$interceptChain,
        sigmaSqChain      = NULL,
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
