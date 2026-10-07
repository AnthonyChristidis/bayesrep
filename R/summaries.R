#' Widely Applicable Information Criterion
#'
#' Extracts the WAIC accumulated while sampling, computed from the pointwise
#' composite log-density that the posterior targets.
#'
#' @section Do not use this to choose lambda2:
#' The composite objective scores each sub-model on how well it predicts the
#' response \emph{on its own}. Repulsion necessarily makes each individual
#' sub-model worse, even when it makes the ensemble better, so WAIC decreases
#' monotonically as \code{lambda2} goes to zero and would always select no
#' repulsion. The same caution applies to \code{numModels}.
#'
#' Use \code{\link{cv.bre}}, which scores held-out error of the ensemble
#' prediction, to choose either. WAIC is retained as a descriptive measure of
#' in-sample composite fit and for comparing models at a fixed \code{lambda2}
#' and \code{numModels}.
#'
#' @param object An object of class \code{"bre"}.
#' @param ... Unused, for method consistency.
#'
#' @return A list with \code{waic}, \code{lppd} (log pointwise predictive
#'   density), \code{pWaic} (effective number of parameters) and \code{elpd}.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(60 * 8), 60, 8)
#' y <- as.numeric(X %*% c(1.5, rep(0, 7)) + rnorm(60))
#' fit <- bre(X, y, numModels = 3, iter = 300, burnin = 150, verbose = FALSE)
#' waic(fit)$waic
#'
#' @export
waic <- function(object, ...) UseMethod("waic")

#' @rdname waic
#' @export
waic.bre <- function(object, ...) object$waic

#' Select features at a controlled false discovery rate
#'
#' Ranks features by their ensemble marginal inclusion probability and admits
#' them while the estimated Bayesian false discovery rate stays within
#' \code{fdr}.
#'
#' @section Why not just threshold eMIP at 0.5:
#' Because \eqn{\mathrm{eMIP}_j = P(\textrm{feature } j \textrm{ is carried by at
#' least one of the } G \textrm{ sub-models})}, a null feature gets \eqn{G}
#' independent chances to be picked up, and the meaning of any fixed cutoff
#' therefore shifts with \code{numModels}. A feature whose per-sub-model
#' inclusion probability is only \eqn{0.12} has
#' \eqn{1 - (1 - 0.12)^5 \approx 0.47} at \eqn{G = 5}: it crosses a 0.5
#' threshold on arithmetic alone. In a simulation with 150 active predictors
#' and 350 nulls, \code{emip > 0.5} admitted 347 of the 350 nulls, while the
#' rule used here admitted 3.
#'
#' Thresholding the per-sub-model probabilities instead does not work either:
#' the posterior is invariant to permuting the sub-model labels, so a feature
#' reliably carried by exactly one of five sub-models has probability near
#' \eqn{1/5} in each of them and never reaches 0.5 anywhere.
#'
#' The estimated false discovery rate among the top \eqn{k} features is the
#' average of \eqn{1 - \mathrm{eMIP}} over those features, which needs no
#' reference to \eqn{G}.
#'
#' @param object An object of class \code{"bre"}.
#' @param fdr Numeric in (0, 1). Target false discovery rate. Default 0.1.
#'
#' @return A data frame of the selected features, ordered by decreasing eMIP,
#'   with their ensemble coefficients on the original data scale. The attributes
#'   \code{targetFdr} and \code{estimatedFdr} record the requested and achieved
#'   values. Zero rows means nothing could be selected at that rate.
#'
#' @seealso \code{\link{summary.bre}} for the full ranking.
#'
#' @examples
#' set.seed(1)
#' n <- 120; p <- 20
#' X <- matrix(rnorm(n * p), n, p)
#' colnames(X) <- paste0("V", seq_len(p))
#' y <- as.numeric(X %*% c(2, -2, 1.5, rep(0, p - 3)) + rnorm(n))
#'
#' fit <- bre(X, y, numModels = 3, iter = 600, burnin = 300,
#'            verbose = FALSE, seed = 1)
#' selectFeatures(fit, fdr = 0.1)
#'
#' @export
selectFeatures <- function(object, fdr = 0.1) {

    stopifnot(inherits(object, "bre"))
    if (length(fdr) != 1L || !is.finite(fdr) || fdr <= 0 || fdr >= 1) {
        stop("'fdr' must be a single number strictly between 0 and 1.")
    }

    ord <- order(object$emip, decreasing = TRUE)
    running <- cumsum(1 - object$emip[ord]) / seq_along(ord)

    k <- if (any(running <= fdr)) max(which(running <= fdr)) else 0L
    idx <- if (k == 0L) integer(0) else ord[seq_len(k)]

    betaOrig <- rowMeans(object$ensembleBetaChain) / object$scaleInfo$xScale

    out <- data.frame(
        Feature     = object$featureNames[idx],
        eMIP        = unname(object$emip[idx]),
        PostMean    = unname(betaOrig[idx]),
        MeanNModels = unname(rowSums(object$pip)[idx]),
        row.names   = NULL,
        stringsAsFactors = FALSE
    )

    attr(out, "targetFdr") <- fdr
    attr(out, "estimatedFdr") <- if (k == 0L) 0 else unname(running[k])
    out
}

#' Posterior co-allocation matrix
#'
#' Computes \eqn{C_{jk} = P(\textrm{features } j \textrm{ and } k \textrm{ are
#' carried by a common sub-model})}, with the ensemble marginal inclusion
#' probability on the diagonal.
#'
#' @section Interpretation:
#' Co-allocation is invariant to relabelling of the \eqn{G} sub-models, so it is
#' well defined despite the label symmetry of the posterior, and under the
#' composite likelihood used by \code{\link{bre}} it is informative: relocating
#' a feature changes each affected sub-model's fit, so the assignment is
#' identified.
#'
#' Read values against the baseline \eqn{1/G} that would arise if features were
#' assigned at random and each were carried exactly once:
#' \describe{
#'   \item{well below \eqn{1/G}}{the ensemble is actively separating the two
#'     features. This is the signature of a collinear group being split across
#'     sub-models, which is the behaviour the repulsive prior exists to produce.}
#'   \item{near \eqn{1/G}}{no detectable relationship.}
#'   \item{well above \eqn{1/G}}{the features complement each other and tend to
#'     be used together, or at least one of them is carried by several
#'     sub-models at once.}
#' }
#'
#' The exact baseline depends on how often each feature is carried: for features
#' in \eqn{m_j} and \eqn{m_k} sub-models it is
#' \eqn{1 - {G - m_j \choose m_k} / {G \choose m_k}}, which equals \eqn{1/G}
#' only when both are carried exactly once. Check \code{MeanNModels} in
#' \code{\link{summary.bre}} before reading too much into a single number.
#'
#' @param object An object of class \code{"bre"}.
#' @param features Optional character vector of feature names, or integer
#'   vector of column indices, to restrict the matrix to. Defaults to the
#'   \code{top} features by ensemble marginal inclusion probability.
#' @param top Integer. Number of top-ranked features to use when
#'   \code{features} is \code{NULL}. Default 25.
#'
#' @return A square numeric matrix with feature names on both margins.
#'
#' @seealso \code{\link{allocation}}, \code{\link{summary.bre}}.
#'
#' @examples
#' set.seed(1)
#' n <- 80; p <- 10
#' X <- matrix(rnorm(n * p), n, p)
#' X[, 2] <- X[, 1] + rnorm(n, sd = 0.1)
#' y <- as.numeric(X %*% c(1, 1, rep(0, p - 2)) + rnorm(n))
#' fit <- bre(X, y, numModels = 2, iter = 400, burnin = 200, verbose = FALSE)
#' round(coallocation(fit, top = 4), 2)
#'
#' @export
coallocation <- function(object, features = NULL, top = 25) {

    stopifnot(inherits(object, "bre"))

    idx <- resolveFeatures(object, features, top)

    cm <- coallocationCpp(gammaRaw = object$gammaChain,
                          p = object$p,
                          numModels = object$numModels,
                          numDraws = object$nKeep,
                          featureIdx = as.integer(idx - 1L))

    dimnames(cm) <- list(object$featureNames[idx], object$featureNames[idx])
    cm
}

#' Label-aligned sub-model inclusion probabilities
#'
#' Returns the posterior inclusion probability of each feature within each
#' sub-model, after relabelling every draw to best match the running aligned
#' mean.
#'
#' @section Reading this summary:
#' The posterior is invariant to permuting all \eqn{G} sub-model labels at once,
#' so the columns of this matrix are only meaningful after the draws have been
#' put on a common labelling, which is what this function does. Once relabelled
#' the columns describe genuinely different sub-models, since the composite
#' likelihood identifies the assignment.
#'
#' A near-zero \code{switchRate} is not in itself a problem. It means the chain
#' stayed on one labelling, which is the usual outcome when the assignment is
#' sharply identified: the posterior genuinely prefers one arrangement of
#' features across sub-models, and the raw \code{object$pip} was already
#' interpretable. Check \code{acceptance$transfer} in the fitted object to
#' confirm the relocation move is running; if it is zero, raise
#' \code{transferMoves} in \code{\link{bre}}.
#'
#' For a summary that needs no relabelling at all, use
#' \code{\link{coallocation}}.
#'
#' @param object An object of class \code{"bre"}.
#' @param maxExact Integer. Solve the relabelling assignment exactly by
#'   enumeration when \code{numModels} is at most this value, greedily above it.
#'   Default 7.
#'
#' @return A list with \code{pip} (p by G matrix) and \code{switchRate}, the
#'   proportion of draws whose optimal relabelling was not the identity.
#'
#' @seealso \code{\link{coallocation}}.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(60 * 8), 60, 8)
#' y <- as.numeric(X %*% c(1.5, rep(0, 7)) + rnorm(60))
#' fit <- bre(X, y, numModels = 2, iter = 300, burnin = 150, verbose = FALSE)
#' allocation(fit)$switchRate
#'
#' @export
allocation <- function(object, maxExact = 7) {

    stopifnot(inherits(object, "bre"))

    res <- alignedPipCpp(gammaRaw = object$gammaChain,
                         p = object$p,
                         numModels = object$numModels,
                         numDraws = object$nKeep,
                         maxExact = as.integer(maxExact))

    dimnames(res$pip) <- list(object$featureNames,
                              paste0("Model", seq_len(object$numModels)))
    res
}

#' Resolve a user feature selection to column indices
#'
#' @param object An object of class \code{"bre"}.
#' @param features Names, indices, or \code{NULL}.
#' @param top Number of top features to take when \code{features} is NULL.
#' @return An integer vector of column indices.
#' @noRd
resolveFeatures <- function(object, features, top) {

    if (is.null(features)) {
        top <- min(as.integer(top), object$p)
        if (top < 1L) stop("'top' must be at least 1.")
        return(order(object$emip, decreasing = TRUE)[seq_len(top)])
    }

    if (is.character(features)) {
        idx <- match(features, object$featureNames)
        if (anyNA(idx)) {
            stop("Unknown feature name(s): ",
                 paste(features[is.na(idx)], collapse = ", "))
        }
        return(idx)
    }

    idx <- as.integer(features)
    if (any(idx < 1L | idx > object$p)) stop("'features' index out of range.")
    idx
}
