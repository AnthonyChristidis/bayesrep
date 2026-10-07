#' @useDynLib bayesrep, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' Fit a Bayesian Repulsive Ensemble (BRE)
#'
#' Fits an ensemble of \code{G} sparse linear or logistic sub-models whose
#' feature-allocation matrix carries an Ising repulsive prior, so that the
#' sub-models are pushed towards structurally diverse predictive pathways. A
#' Beta hyperprior on the feature-specific baseline inclusion probabilities
#' lets universally predictive covariates ("super-predictors") be shared across
#' sub-models in spite of the repulsion.
#'
#' @section Model:
#' Every sub-model is fitted to the response in full. The posterior is built
#' from the composite likelihood
#' \deqn{p(y \mid B, \sigma^2) \propto \prod_{g=1}^{G}
#'       N\!\left(y \mid \alpha_0 + X\beta^{(g)}, \sigma^2\right)^{1/G},}
#' and the logistic analogue for \code{family = "binomial"}. The \eqn{1/G}
#' tempering keeps the total information in the posterior equal to \eqn{n}
#' observations rather than \eqn{Gn}, and enters each sub-model's conditional as
#' an inflated variance \eqn{G\sigma^2}.
#'
#' Predictions come from the ensemble average
#' \eqn{\bar\beta = G^{-1}\sum_g \beta^{(g)}}. Each cell of the \eqn{p \times G}
#' allocation matrix carries a spike-and-slab prior, and each row carries an
#' Ising prior with pairwise repulsion \eqn{\lambda_2} and external field
#' \eqn{logit(\theta_j)}.
#'
#' @section Why the sub-models become distinct:
#' Two forces act together. Because every sub-model has to predict \eqn{y} on
#' its own, each needs a competent set of features. Because the Ising prior
#' penalises reusing a feature another sub-model already holds, they cannot all
#' reach for the same one. On a collinear block the result is that the block is
#' split across sub-models, each of which reconstructs the signal from its own
#' member: this is the Bayesian analogue of split regularised regression.
#'
#' Unlike an averaged-likelihood formulation, here the assignment of features to
#' sub-models \strong{is} identified: relocating a feature changes each affected
#' sub-model's fit and therefore the posterior density. \code{\link{coallocation}}
#' is consequently informative. Co-allocation below \eqn{1/G} between two
#' features means the ensemble is actively separating them, which is the
#' signature of collinearity being resolved; above \eqn{1/G} means they are
#' complementary and tend to be used together.
#'
#' The model retains the usual label symmetry: permuting all \eqn{G} sub-model
#' labels at once leaves the posterior unchanged. \code{\link{coallocation}} is
#' invariant to that, and \code{\link{allocation}} removes it by relabelling.
#'
#' @param X Numeric matrix (n x p). The design matrix. Column names, if present,
#'   are retained throughout the fitted object.
#' @param y Numeric vector (length n). Continuous response for
#'   \code{family = "gaussian"}, 0/1 response for \code{family = "binomial"}.
#' @param family Character. \code{"gaussian"} for continuous responses,
#'   \code{"binomial"} for binary classification.
#' @param numModels Integer. The number of sub-models \eqn{G} in the ensemble.
#'   Default 5. Compare values with \code{\link{waic}}.
#' @param iter Integer. Total number of MCMC iterations. Default 5000.
#' @param burnin Integer. Iterations discarded before storage. Default
#'   \code{floor(iter / 2)}.
#' @param thin Integer. Keep every \code{thin}-th post-burn-in draw. Default 1.
#' @param tauSq Numeric. Prior variance of a non-zero sub-model coefficient.
#'   Because every sub-model predicts the response in full, this is on the scale
#'   of an ordinary regression coefficient and does not depend on
#'   \code{numModels}. Default 1, suited to standardised data.
#' @param lambda2 Numeric or \code{NULL}. The repulsion strength on the
#'   log-odds scale: each sub-model that already carries a feature multiplies
#'   the prior odds of another one taking it by \eqn{e^{-\lambda_2}}. Default 2.
#'   \code{lambda2 = 0} turns the repulsion off and is the matched control.
#'
#'   \strong{This is a genuine tuning parameter and the default will not suit
#'   every problem.} It competes on the log-odds scale against the evidence for
#'   including a feature, which grows with the sample size, the signal strength
#'   and \code{learningRate}, so the value needed to separate a correlated
#'   group can be an order of magnitude larger on one dataset than another.
#'   Choose it with \code{\link{cv.bre}} rather than by hand, and inspect
#'   \code{\link{coallocation}} to confirm the ensemble is actually separating
#'   anything.
#'
#'   Setting \code{lambda2 = NULL} samples it under the \code{lambda2Prior}, but
#'   this is \strong{not recommended}: under the composite likelihood every
#'   sub-model individually prefers to keep every useful predictor, so the
#'   posterior for \code{lambda2} is pulled towards zero and the ensemble
#'   collapses to \eqn{G} near-copies of the same model. See
#'   \code{\link{cv.bre}}.
#' @param lambda2Prior Numeric vector of length two giving the shape and rate of
#'   the Gamma prior on \code{lambda2}. Default \code{c(2, 1)}. Used only when
#'   \code{lambda2 = NULL}.
#' @param thetaShape1,thetaShape2 Numeric. Shape parameters of the
#'   \eqn{Beta(a, b)} hyperprior on the feature-specific baseline inclusion
#'   probability \eqn{\theta_j}. \code{thetaShape2 = NULL} (the default) uses
#'   \code{p}, giving a prior mean inclusion rate of \eqn{1 / (1 + p)}.
#' @param nu0,sigma0Sq Numeric. Prior degrees of freedom and scale of the
#'   inverse-Gamma prior on \eqn{\sigma^2}. Gaussian family only. Default 1, 1.
#' @param transferMoves Integer. Number of Metropolis proposals per iteration
#'   that relocate a feature from one sub-model to another. Default 10. A
#'   single-site sweep cannot relocate a feature without passing through a
#'   duplicated state that the repulsion penalises, so without these moves the
#'   assignment mixes very slowly. Set to 0 to disable.
#' @param thetaSweeps Integer. Metropolis sweeps of the \eqn{\theta_j} update
#'   per iteration. Default 5; the independence proposal is centred on the
#'   no-repulsion conjugate posterior, which mixes slowly when \code{lambda2}
#'   is large.
#' @param learningRate Numeric in \eqn{(0, 1]}. The composite-likelihood
#'   learning rate \eqn{\eta}: the power to which each sub-model's likelihood is
#'   raised. Default 1, at which every sub-model is fitted to the data in full
#'   and the objective matches the frequentist split-regression loss.
#'
#'   Lower values temper each sub-model's evidence. The value \eqn{1/G} makes
#'   the composite normalising constant that of a single \eqn{n}-observation
#'   Gaussian, which is appealing in principle but costs selection power: a
#'   sub-model then sees only \eqn{1/G} of the evidence for including a feature
#'   while facing an unchanged sparsity threshold, and when \eqn{p \gg n} this
#'   can suppress selection entirely. Tune it with \code{\link{cv.bre}} if in
#'   doubt.
#' @param pgTrunc Integer. Number of series terms used by the Polya-Gamma
#'   sampler for the non-integer shape \eqn{1/G} required by the composite
#'   logistic likelihood. Default 50. The expectation of the discarded tail is
#'   added back analytically, so the sampler is unbiased in the mean; raise this
#'   only if you want to reduce the residual tail variance. Binomial family
#'   only.
#' @param scaleData Logical. Standardise the columns of \code{X} (and centre
#'   \code{y} for the Gaussian family) before fitting. Default \code{TRUE}.
#'   Coefficients are reported on the original scale regardless.
#' @param storeChains Logical. Also retain the per-sub-model coefficient and
#'   \eqn{\theta} chains, which cost \eqn{8pG} and \eqn{8p} bytes per stored
#'   draw. Default \code{FALSE}; the ensemble coefficient chain and the
#'   one-byte-per-entry allocation chain are always retained.
#' @param maxMemoryGb Numeric. Refuse to start if the chains would exceed this
#'   many gigabytes. Default 4. Increase \code{thin} or \code{burnin} to fit
#'   within a budget.
#' @param verbose Logical. Report progress. Default \code{TRUE}.
#' @param seed Integer or \code{NULL}. Optional seed for reproducibility.
#'
#' @return An object of class \code{"bre"}: a list with the retained chains,
#'   posterior summaries (\code{emip}, \code{pip}, \code{waic}), MCMC
#'   diagnostics (\code{acceptance}, \code{modelSize}) and the information
#'   needed to map coefficients back to the original data scale.
#'
#' @seealso \code{\link{summary.bre}} for feature ranking,
#'   \code{\link{coallocation}} for label-invariant allocation summaries,
#'   \code{\link{waic}} for choosing \code{numModels}.
#'
#' @examples
#' set.seed(1)
#' n <- 60; p <- 12
#' X <- matrix(rnorm(n * p), n, p)
#' X[, 2] <- X[, 1] + rnorm(n, sd = 0.1)
#' y <- as.numeric(X %*% c(1, 1, rep(0, p - 2)) + rnorm(n))
#'
#' fit <- bre(X, y, numModels = 3, iter = 400, burnin = 200, verbose = FALSE)
#' fit
#' summary(fit)
#' coallocation(fit, top = 4)
#'
#' @export
bre <- function(X, y,
                family = c("gaussian", "binomial"),
                numModels = 5,
                iter = 5000,
                burnin = floor(iter / 2),
                thin = 1,
                tauSq = 1,
                lambda2 = 2,
                lambda2Prior = c(2, 1),
                thetaShape1 = 1,
                thetaShape2 = NULL,
                nu0 = 1,
                sigma0Sq = 1,
                transferMoves = 10,
                thetaSweeps = 5,
                pgTrunc = 50,
                learningRate = 1,
                scaleData = TRUE,
                storeChains = FALSE,
                maxMemoryGb = 4,
                verbose = TRUE,
                seed = NULL) {

    family <- match.arg(family)

    if (!is.null(seed)) set.seed(seed)

    # ---- Input validation -------------------------------------------------
    if (!is.matrix(X)) {
        if (is.data.frame(X)) X <- as.matrix(X) else stop("'X' must be a matrix.")
    }
    if (!is.numeric(X)) stop("'X' must be numeric.")
    y <- as.numeric(y)
    if (nrow(X) != length(y)) {
        stop("The number of rows of 'X' must match the length of 'y'.")
    }
    if (anyNA(X) || anyNA(y)) stop("'X' and 'y' must not contain missing values.")
    if (family == "binomial" && !all(y %in% c(0, 1))) {
        stop("For family = 'binomial', 'y' must contain only 0s and 1s.")
    }

    n <- nrow(X)
    p <- ncol(X)

    if (numModels < 1 || numModels != round(numModels)) {
        stop("'numModels' must be a positive integer.")
    }
    if (iter < 1 || iter != round(iter)) stop("'iter' must be a positive integer.")
    if (burnin < 0 || burnin != round(burnin)) stop("'burnin' must be a non-negative integer.")
    if (burnin >= iter) stop("'burnin' must be smaller than 'iter'.")
    if (thin < 1 || thin != round(thin)) stop("'thin' must be a positive integer.")
    if (tauSq <= 0) stop("'tauSq' must be positive.")
    if (nu0 <= 0 || sigma0Sq <= 0) stop("'nu0' and 'sigma0Sq' must be positive.")
    if (thetaShape1 <= 0) stop("'thetaShape1' must be positive.")
    if (thetaSweeps < 1 || thetaSweeps != round(thetaSweeps)) {
        stop("'thetaSweeps' must be a positive integer.")
    }

    if (is.null(learningRate)) learningRate <- 1
    if (length(learningRate) != 1L || !is.finite(learningRate) ||
        learningRate <= 0 || learningRate > 1) {
        stop("'learningRate' must be a single number in (0, 1], or NULL.")
    }
    if (pgTrunc < 1 || pgTrunc != round(pgTrunc)) {
        stop("'pgTrunc' must be a positive integer.")
    }
    if (transferMoves < 0 || transferMoves != round(transferMoves)) {
        stop("'transferMoves' must be a non-negative integer.")
    }

    if (is.null(thetaShape2)) thetaShape2 <- p
    if (thetaShape2 <= 0) stop("'thetaShape2' must be positive.")

    learnLambda2 <- is.null(lambda2)
    if (learnLambda2) {
        if (length(lambda2Prior) != 2L || any(lambda2Prior <= 0)) {
            stop("'lambda2Prior' must be two positive numbers (shape, rate).")
        }
        lambda2 <- lambda2Prior[1L] / lambda2Prior[2L]
    } else {
        if (length(lambda2) != 1L || !is.finite(lambda2) || lambda2 < 0) {
            stop("'lambda2' must be a single non-negative number, or NULL to learn it.")
        }
        if (lambda2 == 0 && verbose) {
            message("lambda2 = 0: the sub-models are independent spike-and-slab fits.")
        }
    }

    featureNames <- colnames(X)
    if (is.null(featureNames)) featureNames <- paste0("X", seq_len(p))

    # A constant column is collinear with the intercept and cannot be
    # standardised; fail loudly rather than propagating NaN through the chain.
    colSds <- apply(X, 2L, stats::sd)
    constCols <- which(colSds == 0)
    if (length(constCols) > 0) {
        stop(sprintf("Column(s) %s of 'X' are constant; remove them before fitting.",
                     paste(utils::head(featureNames[constCols], 10L), collapse = ", ")))
    }

    # ---- Memory budget ----------------------------------------------------
    nKeep <- length(seq.int(burnin + 1L, iter, by = thin))
    bytes <- 8 * p * nKeep +                      # ensemble coefficient chain
             1 * p * numModels * nKeep +          # allocation chain (raw)
             8 * numModels * nKeep +              # per sub-model sizes
             (if (storeChains) 8 * p * numModels * nKeep + 8 * p * nKeep else 0)
    gb <- bytes / 1024^3
    if (gb > maxMemoryGb) {
        stop(sprintf(paste0("The requested chains need about %.1f GB, above maxMemoryGb = %.1f. ",
                            "Increase 'thin' or 'burnin', set storeChains = FALSE, ",
                            "or raise 'maxMemoryGb'."), gb, maxMemoryGb))
    }

    # ---- Standardisation --------------------------------------------------
    if (scaleData) {
        Xs <- scale(X)
        xCenter <- attr(Xs, "scaled:center")
        xScale  <- attr(Xs, "scaled:scale")
        attributes(Xs) <- list(dim = dim(Xs))
    } else {
        Xs <- X
        xCenter <- rep(0, p)
        xScale  <- rep(1, p)
    }
    if (scaleData && family == "gaussian") {
        yCenter <- mean(y)
        ys <- y - yCenter
    } else {
        yCenter <- 0
        ys <- y
    }

    if (verbose) {
        message(sprintf("Fitting %s BRE: G = %d, p = %d, n = %d, %d iterations (%d retained, %.2f GB).",
                        family, numModels, p, n, iter, nKeep, gb))
        message(sprintf("Repulsion lambda2 is %s.",
                        if (learnLambda2) sprintf("learned, Gamma(%g, %g) prior",
                                                  lambda2Prior[1L], lambda2Prior[2L])
                        else sprintf("fixed at %g", lambda2)))
    }

    mcmcArgs <- list(y = ys, X = Xs,
                     numModels = numModels, iter = iter, burnin = burnin, thin = thin,
                     tauSq = tauSq,
                     lambda2 = lambda2, learnLambda2 = learnLambda2,
                     lambda2Shape = lambda2Prior[1L], lambda2Rate = lambda2Prior[2L],
                     thetaShape1 = thetaShape1, thetaShape2 = thetaShape2,
                     transferMoves = transferMoves,
                     thetaSweeps = thetaSweeps,
                     learningRate = learningRate,
                     storeChains = storeChains, verbose = verbose)

    if (family == "gaussian") {
        mcmcArgs$nu0 <- nu0
        mcmcArgs$sigma0Sq <- sigma0Sq
        res <- do.call(runMcmcContinuous, mcmcArgs)
    } else {
        mcmcArgs$pgTrunc <- pgTrunc
        res <- do.call(runMcmcBinary, mcmcArgs)
    }

    out <- c(res, list(
        call         = match.call(),
        family       = family,
        numModels    = numModels,
        n            = n,
        p            = p,
        featureNames = featureNames,
        iter         = iter,
        burnin       = burnin,
        thin         = thin,
        nKeep        = nKeep,
        tauSq        = tauSq,
        learnLambda2 = learnLambda2,
        transferMoves = transferMoves,
        learningRate = learningRate,
        scaleData    = scaleData,
        scaleInfo    = list(xCenter = xCenter, xScale = xScale, yCenter = yCenter)
    ))

    rownames(out$pip) <- featureNames
    colnames(out$pip) <- paste0("Model", seq_len(numModels))
    names(out$emip) <- featureNames
    names(out$thetaMean) <- featureNames
    names(out$thetaSd) <- featureNames
    names(out$emip) <- featureNames
    rownames(out$modelSize) <- paste0("Model", seq_len(numModels))
    rownames(out$ensembleBetaChain) <- featureNames

    class(out) <- "bre"
    out
}
