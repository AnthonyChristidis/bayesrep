#' Choose the repulsion strength (and G) by cross-validation
#'
#' Runs K-fold cross-validation over a grid of \code{lambda2} and
#' \code{numModels} values, scoring each by held-out predictive error of the
#' ensemble, and refits at the best setting.
#'
#' @section Why cross-validation and not WAIC:
#' Under the composite likelihood each sub-model is scored on how well it
#' predicts the response \emph{on its own}, and repulsion necessarily makes each
#' individual sub-model worse even when it makes the ensemble better. The
#' posterior for \code{lambda2} therefore drifts towards zero, and WAIC, which
#' measures the same composite in-sample fit, also prefers zero. Neither is a
#' valid criterion for choosing the repulsion strength.
#'
#' Held-out error of the \emph{ensemble} prediction is the right target, which
#' is what this function uses. \code{\link{waic}} remains available but should
#' not be used to select \code{lambda2}.
#'
#' @param X Numeric matrix (n x p). The design matrix.
#' @param y Numeric vector (length n). The response.
#' @param family Character. \code{"gaussian"} or \code{"binomial"}.
#' @param lambda2 Numeric vector of candidate repulsion strengths. Default
#'   \code{c(0, 1, 2, 4, 8, 16, 32)}; include 0 so the no-repulsion control is
#'   always in the comparison. The useful range scales roughly with the inverse
#'   of \code{learningRate}, because the repulsion competes on the log-odds
#'   scale against evidence that the learning rate multiplies, so a grid that
#'   suits \eqn{\eta = 1} will be far too coarse at \eqn{\eta = 1/G}.
#' @param numModels Integer vector of candidate ensemble sizes. Default 5. A
#'   vector here is crossed with \code{lambda2}.
#' @param learningRate Numeric vector of candidate learning rates, or
#'   \code{NULL} to leave it at the \code{\link{bre}} default of 1. Crossed with
#'   the other grids. Worth tuning when \eqn{p \gg n}; see \code{\link{bre}}.
#' @param tauSq Numeric vector of candidate slab variances, or \code{NULL} to
#'   leave it at the \code{\link{bre}} default. Worth crossing with
#'   \code{lambda2}: stronger repulsion concentrates each feature into fewer
#'   sub-models, which changes the scale of the coefficients those sub-models
#'   need, so a slab variance tuned at one repulsion strength can handicap
#'   another.
#' @param nfolds Integer. Number of folds. Default 5.
#' @param foldid Optional integer vector of length n assigning each observation
#'   to a fold, overriding \code{nfolds}.
#' @param measure Character. \code{"mse"} for the Gaussian family;
#'   \code{"deviance"}, \code{"auc"} or \code{"misclass"} for the binomial
#'   family. Defaults to \code{"mse"} and \code{"deviance"} respectively.
#' @param oneSe Logical. Return the most parsimonious setting within one
#'   standard error of the best, as is conventional for cross-validated
#'   selection. Default \code{TRUE}; "most parsimonious" here means the largest
#'   \code{lambda2} and the smallest \code{numModels}, which is the sparsest,
#'   most strongly separated ensemble that the data support.
#' @param refit Logical. Refit on the full data at the chosen setting and return
#'   it as \code{$fit}. Default \code{TRUE}.
#' @param verbose Logical. Report progress. Default \code{TRUE}.
#' @param seed Integer or \code{NULL}. Seed for fold assignment and sampling.
#' @param ... Further arguments passed to \code{\link{bre}}, for example
#'   \code{iter}, \code{burnin}, \code{thetaShape1}. Do not pass grid
#'   parameters this way; use the dedicated arguments above.
#'
#' @return An object of class \code{"cv.bre"} with the scored \code{grid}, the
#'   per-fold errors, the selected \code{best} and \code{bestOneSe} rows, and
#'   (if \code{refit}) the refitted model in \code{$fit}.
#'
#' @seealso \code{\link{bre}}, \code{\link{plot.cv.bre}}.
#'
#' @examples
#' set.seed(1)
#' n <- 80; p <- 10
#' X <- matrix(rnorm(n * p), n, p)
#' X[, 2] <- X[, 1] + rnorm(n, sd = 0.1)
#' y <- as.numeric(X %*% c(1.5, 1.5, -2, rep(0, p - 3)) + rnorm(n))
#'
#' cv <- cv.bre(X, y, lambda2 = c(0, 2), nfolds = 3,
#'              iter = 300, burnin = 150, verbose = FALSE, seed = 1)
#' cv$grid
#'
#' @export
cv.bre <- function(X, y,
                   family = c("gaussian", "binomial"),
                   lambda2 = c(0, 1, 2, 4, 8, 16, 32),
                   numModels = 5,
                   learningRate = NULL,
                   tauSq = NULL,
                   nfolds = 5,
                   foldid = NULL,
                   measure = NULL,
                   oneSe = TRUE,
                   refit = TRUE,
                   verbose = TRUE,
                   seed = NULL,
                   ...) {

    family <- match.arg(family)
    if (!is.null(seed)) set.seed(seed)

    if (!is.matrix(X)) X <- as.matrix(X)
    y <- as.numeric(y)
    n <- nrow(X)
    if (length(y) != n) stop("The number of rows of 'X' must match the length of 'y'.")

    if (is.null(measure)) {
        measure <- if (family == "gaussian") "mse" else "deviance"
    }
    measure <- match.arg(measure, c("mse", "deviance", "auc", "misclass"))
    if (family == "gaussian" && measure != "mse") {
        stop("For family = 'gaussian' the only available measure is 'mse'.")
    }
    if (family == "binomial" && measure == "mse") {
        stop("For family = 'binomial' use 'deviance', 'auc' or 'misclass'.")
    }

    if (any(lambda2 < 0)) stop("'lambda2' values must be non-negative.")
    if (any(numModels < 1)) stop("'numModels' values must be positive integers.")

    if (is.null(foldid)) {
        if (nfolds < 2 || nfolds > n) stop("'nfolds' must be between 2 and n.")
        foldid <- sample(rep(seq_len(nfolds), length.out = n))
    } else {
        foldid <- as.integer(foldid)
        if (length(foldid) != n) stop("'foldid' must have length n.")
        nfolds <- length(unique(foldid))
    }

    tauSqDefault <- eval(formals(bre)$tauSq)
    if (!is.null(tauSq) && any(tauSq <= 0)) {
        stop("'tauSq' values must be positive.")
    }
    if (!is.null(learningRate) &&
            (any(learningRate <= 0) || any(learningRate > 1))) {
        stop("'learningRate' values must lie in (0, 1].")
    }

    grid <- expand.grid(
        lambda2 = lambda2,
        numModels = numModels,
        learningRate = if (is.null(learningRate)) NA_real_ else learningRate,
        tauSq = if (is.null(tauSq)) NA_real_ else tauSq,
        KEEP.OUT.ATTRS = FALSE)
    nGrid <- nrow(grid)
    errors <- matrix(NA_real_, nrow = nGrid, ncol = nfolds)

    folds <- sort(unique(foldid))

    for (k in seq_along(folds)) {

        inTest <- foldid == folds[k]
        if (all(inTest) || !any(inTest)) next

        Xtr <- X[!inTest, , drop = FALSE]; ytr <- y[!inTest]
        Xte <- X[inTest, , drop = FALSE];  yte <- y[inTest]

        # A column that is constant within a training fold cannot be
        # standardised; drop it from this fold for every grid point alike.
        keep <- apply(Xtr, 2L, stats::sd) > 0
        if (!all(keep)) {
            Xtr <- Xtr[, keep, drop = FALSE]
            Xte <- Xte[, keep, drop = FALSE]
        }

        for (r in seq_len(nGrid)) {

            fitK <- bre(Xtr, ytr, family = family,
                        numModels = grid$numModels[r],
                        lambda2 = grid$lambda2[r],
                        tauSq = if (is.na(grid$tauSq[r])) tauSqDefault else grid$tauSq[r],
                        learningRate = if (is.na(grid$learningRate[r])) NULL
                                       else grid$learningRate[r],
                        verbose = FALSE, ...)

            pred <- predict(fitK, Xte, type = "response")
            errors[r, k] <- cvLoss(yte, pred, measure)
        }

        if (verbose) {
            message(sprintf("  fold %d of %d complete", k, length(folds)))
        }
    }

    higherIsBetter <- measure == "auc"

    grid$cvMean <- rowMeans(errors, na.rm = TRUE)
    grid$cvSe <- apply(errors, 1L, function(e) {
        e <- e[is.finite(e)]
        if (length(e) < 2L) NA_real_ else stats::sd(e) / sqrt(length(e))
    })

    bestIdx <- if (higherIsBetter) which.max(grid$cvMean) else which.min(grid$cvMean)

    # One-standard-error rule: among settings statistically indistinguishable
    # from the best, take the most strongly separated ensemble.
    seThreshold <- grid$cvMean[bestIdx] +
        (if (higherIsBetter) -1 else 1) * (grid$cvSe[bestIdx] %||% 0)
    withinSe <- if (higherIsBetter) {
        which(grid$cvMean >= seThreshold)
    } else {
        which(grid$cvMean <= seThreshold)
    }
    if (length(withinSe) == 0L) withinSe <- bestIdx
    oneSeIdx <- withinSe[order(-grid$lambda2[withinSe],
                               grid$numModels[withinSe])][1L]

    chosen <- if (oneSe) oneSeIdx else bestIdx

    out <- list(
        grid      = grid,
        errors    = errors,
        foldid    = foldid,
        measure   = measure,
        family    = family,
        best      = grid[bestIdx, , drop = FALSE],
        bestOneSe = grid[oneSeIdx, , drop = FALSE],
        chosen    = grid[chosen, , drop = FALSE],
        oneSe     = oneSe,
        call      = match.call()
    )

    if (refit) {
        if (verbose) {
            message(sprintf("Refitting at lambda2 = %g, numModels = %d",
                            grid$lambda2[chosen], grid$numModels[chosen]))
        }
        out$fit <- bre(X, y, family = family,
                       numModels = grid$numModels[chosen],
                       lambda2 = grid$lambda2[chosen],
                       tauSq = if (is.na(grid$tauSq[chosen])) tauSqDefault else grid$tauSq[chosen],
                       learningRate = if (is.na(grid$learningRate[chosen])) NULL
                                      else grid$learningRate[chosen],
                       verbose = FALSE, ...)
    }

    class(out) <- "cv.bre"
    out
}

#' Held-out loss for cross-validation
#'
#' @param y Observed responses.
#' @param pred Predicted means or probabilities.
#' @param measure One of "mse", "deviance", "auc", "misclass".
#' @return A single numeric loss.
#' @noRd
cvLoss <- function(y, pred, measure) {
    switch(measure,
        mse = mean((y - pred)^2),
        deviance = {
            pr <- pmin(pmax(pred, 1e-10), 1 - 1e-10)
            -2 * mean(y * log(pr) + (1 - y) * log(1 - pr))
        },
        misclass = mean((pred > 0.5) != (y == 1)),
        auc = {
            n1 <- sum(y == 1); n0 <- sum(y == 0)
            if (n1 == 0 || n0 == 0) return(NA_real_)
            r <- rank(pred)
            (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
        },
        stop("Unknown measure.")
    )
}

`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a

#' Print a cross-validated BRE selection
#'
#' @param x An object of class \code{"cv.bre"}.
#' @param ... Unused, for method consistency.
#' @return \code{x}, invisibly.
#' @export
print.cv.bre <- function(x, ...) {

    cat("\nCross-validated BRE selection\n")
    cat("-----------------------------\n")
    cat(sprintf("Family  : %s\n", x$family))
    cat(sprintf("Measure : %s (%s is better)\n", x$measure,
                if (x$measure == "auc") "higher" else "lower"))
    cat(sprintf("Folds   : %d\n\n", length(unique(x$foldid))))

    g <- x$grid
    g$cvMean <- round(g$cvMean, 4)
    g$cvSe <- round(g$cvSe, 4)
    print(g, row.names = FALSE)

    cat(sprintf("\nBest            : lambda2 = %g, numModels = %d\n",
                x$best$lambda2, x$best$numModels))
    cat(sprintf("Within one SE   : lambda2 = %g, numModels = %d\n",
                x$bestOneSe$lambda2, x$bestOneSe$numModels))
    cat(sprintf("Used for refit  : lambda2 = %g, numModels = %d%s\n",
                x$chosen$lambda2, x$chosen$numModels,
                if (x$oneSe) " (one-SE rule)" else ""))

    invisible(x)
}

#' Plot a cross-validation curve
#'
#' @param x An object of class \code{"cv.bre"}.
#' @param ... Unused, for method consistency.
#' @return Invisibly \code{NULL}; called for the plot side effect.
#' @export
plot.cv.bre <- function(x, ...) {

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar), add = TRUE)

    g <- x$grid
    sizes <- sort(unique(g$numModels))
    cols <- grDevices::hcl.colors(max(length(sizes), 2), "Dark 3")

    graphics::par(mar = c(5, 5, 4, 2))
    graphics::plot(range(g$lambda2), range(c(g$cvMean - g$cvSe, g$cvMean + g$cvSe),
                                           na.rm = TRUE),
                   type = "n", xlab = expression(lambda[2]),
                   ylab = sprintf("CV %s", x$measure),
                   main = "Cross-validated repulsion strength")

    for (i in seq_along(sizes)) {
        sub <- g[g$numModels == sizes[i], , drop = FALSE]
        sub <- sub[order(sub$lambda2), ]
        graphics::segments(sub$lambda2, sub$cvMean - sub$cvSe,
                           sub$lambda2, sub$cvMean + sub$cvSe, col = cols[i])
        graphics::lines(sub$lambda2, sub$cvMean, col = cols[i], lwd = 2)
        graphics::points(sub$lambda2, sub$cvMean, col = cols[i], pch = 19)
    }

    graphics::abline(v = x$chosen$lambda2, lty = 2, col = "darkred")
    if (length(sizes) > 1L) {
        graphics::legend("topright", legend = paste("G =", sizes),
                         col = cols[seq_along(sizes)], lwd = 2, bty = "n")
    }

    invisible(NULL)
}
