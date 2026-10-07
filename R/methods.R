#' Print a Bayesian Repulsive Ensemble fit
#'
#' @param x An object of class \code{"bre"}.
#' @param ... Unused, for method consistency.
#'
#' @return \code{x}, invisibly.
#'
#' @export
print.bre <- function(x, ...) {

    cat("\nBayesian Repulsive Ensemble (BRE)\n")
    cat("---------------------------------\n")
    cat(sprintf("Family        : %s\n", x$family))
    cat(sprintf("Sub-models (G): %d\n", x$numModels))
    cat(sprintf("Data          : n = %d, p = %d\n", x$n, x$p))
    cat(sprintf("MCMC          : %d iterations, %d burn-in, thin %d, %d retained\n",
                x$iter, x$burnin, x$thin, x$nKeep))

    if (x$learnLambda2) {
        cat(sprintf("Repulsion     : learned, posterior mean lambda2 = %.3f (95%% CI %.3f, %.3f)\n",
                    mean(x$lambda2Chain),
                    stats::quantile(x$lambda2Chain, 0.025),
                    stats::quantile(x$lambda2Chain, 0.975)))
        cat(sprintf("                MH acceptance %.1f%%\n", 100 * x$acceptance$lambda2))
    } else {
        cat(sprintf("Repulsion     : fixed, lambda2 = %.3f\n", x$lambda2Chain[1L]))
    }

    cat(sprintf("Slab variance : tauSq = %g (sub-model coefficient scale)\n", x$tauSq))
    cat(sprintf("theta MH acc. : %.1f%%\n", 100 * x$acceptance$theta))
    if (!is.na(x$acceptance$transfer)) {
        cat(sprintf("Transfer acc. : %.1f%% of %d proposals per iteration\n",
                    100 * x$acceptance$transfer, x$transferMoves))
    }
    cat(sprintf("Sub-model size: mean %.1f active features (range %d to %d)\n",
                mean(x$modelSize), min(x$modelSize), max(x$modelSize)))
    cat(sprintf("WAIC          : %.1f (pWAIC %.1f)\n", x$waic$waic, x$waic$pWaic))
    sel <- selectFeatures(x, fdr = 0.1)
    cat(sprintf("Selected (FDR 10%%): %d features\n", nrow(sel)))

    cat("\nUse summary() for the feature ranking, selectFeatures() to select at a\n")
    cat("controlled false discovery rate, coallocation() for the allocation\n")
    cat("structure, and plot() for diagnostics.\n")

    invisible(x)
}

#' Summarise a Bayesian Repulsive Ensemble fit
#'
#' Ranks features by their ensemble marginal inclusion probability (eMIP), the
#' posterior probability that a feature is carried by at least one sub-model.
#' eMIP is estimated directly from the retained draws as the proportion in
#' which \eqn{\sum_g \gamma_j^{(g)} > 0}; it is not approximated by a product
#' over sub-models, which would assume an independence the repulsive prior is
#' designed to violate.
#'
#' eMIP is a sound \emph{ranking} score and a correct per-feature probability,
#' but it is \strong{not} a selection rule. Because a null feature gets \eqn{G}
#' chances to be picked up by some sub-model, the meaning of a fixed cutoff such
#' as 0.5 changes with \code{numModels}; in a simulation with 350 null
#' predictors and \eqn{G = 5}, \code{emip > 0.5} admitted 347 of them. Use
#' \code{\link{selectFeatures}} to select at a controlled false discovery rate.
#'
#' @param object An object of class \code{"bre"}.
#' @param level Numeric. Credible interval level. Default 0.95.
#' @param ... Unused, for method consistency.
#'
#' @return A data frame of class \code{"summary.bre"} with one row per feature,
#'   ordered by decreasing eMIP, holding the posterior mean and credible
#'   interval of the ensemble coefficient on the original data scale, the
#'   posterior mean of \eqn{\theta_j}, and the posterior mean number of
#'   sub-models carrying the feature.
#'
#' @export
summary.bre <- function(object, level = 0.95, ...) {

    if (level <= 0 || level >= 1) stop("'level' must lie strictly between 0 and 1.")
    probs <- c((1 - level) / 2, 1 - (1 - level) / 2)

    # Ensemble coefficient chain is on the standardised scale; map back.
    scaleVec <- object$scaleInfo$xScale
    betaChainOrig <- object$ensembleBetaChain / scaleVec

    postMean <- rowMeans(betaChainOrig)
    ci <- apply(betaChainOrig, 1L, stats::quantile, probs = probs)

    res <- data.frame(
        Feature     = object$featureNames,
        eMIP        = object$emip,
        PostMean    = postMean,
        Lower       = ci[1L, ],
        Upper       = ci[2L, ],
        Theta       = object$thetaMean,
        MeanNModels = rowSums(object$pip),
        row.names   = NULL,
        stringsAsFactors = FALSE
    )

    res <- res[order(res$eMIP, decreasing = TRUE), ]
    rownames(res) <- NULL

    attr(res, "level") <- level
    attr(res, "numModels") <- object$numModels
    attr(res, "burnin") <- object$burnin
    class(res) <- c("summary.bre", "data.frame")
    res
}

#' Print a Bayesian Repulsive Ensemble summary
#'
#' @param x An object of class \code{"summary.bre"}.
#' @param top Integer. Number of leading rows to display. Default 15.
#' @param digits Integer. Significant digits. Default 4.
#' @param ... Unused, for method consistency.
#'
#' @return \code{x}, invisibly.
#'
#' @export
print.summary.bre <- function(x, top = 15, digits = 4, ...) {

    cat("\nBRE feature ranking by ensemble marginal inclusion probability\n")
    cat(sprintf("G = %d sub-models, %g%% credible intervals, %d burn-in draws discarded\n",
                attr(x, "numModels"), 100 * attr(x, "level"), attr(x, "burnin")))
    cat(sprintf("Showing %d of %d features\n\n", min(top, nrow(x)), nrow(x)))

    body <- as.data.frame(x)[seq_len(min(top, nrow(x))), , drop = FALSE]
    print(body, row.names = FALSE, digits = digits)

    cat("\nMeanNModels is the posterior mean number of sub-models carrying the feature;\n")
    cat("values near G indicate a shared super-predictor, values near 1 a split pathway.\n")

    invisible(x)
}

#' Extract ensemble coefficients from a Bayesian Repulsive Ensemble
#'
#' @param object An object of class \code{"bre"}.
#' @param scaled Logical. Return coefficients on the internal standardised
#'   scale rather than the original data scale. Default \code{FALSE}.
#' @param intercept Logical. Prepend the intercept. Default \code{TRUE}.
#' @param ... Unused, for method consistency.
#'
#' @return A named numeric vector of posterior mean ensemble coefficients.
#'
#' @export
coef.bre <- function(object, scaled = FALSE, intercept = TRUE, ...) {

    betaStd <- rowMeans(object$ensembleBetaChain)
    a0Std <- mean(object$interceptChain)

    if (scaled) {
        out <- betaStd
        names(out) <- object$featureNames
        if (intercept) out <- c("(Intercept)" = a0Std, out)
        return(out)
    }

    betaOrig <- betaStd / object$scaleInfo$xScale
    a0 <- a0Std + object$scaleInfo$yCenter -
        sum(betaOrig * object$scaleInfo$xCenter)

    names(betaOrig) <- object$featureNames
    if (intercept) betaOrig <- c("(Intercept)" = a0, betaOrig)
    betaOrig
}

#' Predict from a Bayesian Repulsive Ensemble
#'
#' Predictions are averaged over the retained posterior draws rather than
#' computed from a single plug-in coefficient vector. For
#' \code{family = "binomial"} and \code{type = "response"} this matters: the
#' inverse logit is applied within each draw before averaging, so the result is
#' a posterior mean probability rather than the probability implied by the mean
#' coefficient.
#'
#' @param object An object of class \code{"bre"}.
#' @param newdata Numeric matrix or data frame of new covariates with \code{p}
#'   columns in the same order as the training matrix.
#' @param type Character. \code{"response"} returns the mean response
#'   (probabilities for the binomial family); \code{"link"} returns the linear
#'   predictor.
#' @param interval Character. \code{"none"} for point predictions,
#'   \code{"credible"} for intervals on the mean response, \code{"prediction"}
#'   for intervals that also include observation noise (Gaussian family only).
#' @param level Numeric. Interval level. Default 0.95.
#' @param ... Unused, for method consistency.
#'
#' @return A numeric vector when \code{interval = "none"}, otherwise a matrix
#'   with columns \code{fit}, \code{lwr} and \code{upr}.
#'
#' @export
predict.bre <- function(object, newdata,
                        type = c("response", "link"),
                        interval = c("none", "credible", "prediction"),
                        level = 0.95, ...) {

    type <- match.arg(type)
    interval <- match.arg(interval)

    if (missing(newdata)) stop("'newdata' is required.")
    if (is.data.frame(newdata)) newdata <- as.matrix(newdata)
    if (!is.matrix(newdata)) newdata <- as.matrix(newdata)
    if (!is.numeric(newdata)) stop("'newdata' must be numeric.")
    if (ncol(newdata) != object$p) {
        stop(sprintf("'newdata' has %d columns but the model was fitted with %d.",
                     ncol(newdata), object$p))
    }
    if (level <= 0 || level >= 1) stop("'level' must lie strictly between 0 and 1.")

    if (interval == "prediction" && object$family != "gaussian") {
        stop("interval = 'prediction' is only available for the Gaussian family.")
    }

    Xs <- scale(newdata,
                center = object$scaleInfo$xCenter,
                scale  = object$scaleInfo$xScale)

    # n_new by nKeep linear predictors, one column per retained draw.
    eta <- Xs %*% object$ensembleBetaChain
    eta <- sweep(eta, 2L, object$interceptChain, "+")

    if (object$family == "gaussian") {
        eta <- eta + object$scaleInfo$yCenter
        if (interval == "prediction") {
            noise <- matrix(stats::rnorm(length(eta), 0,
                                         rep(sqrt(object$sigmaSqChain), each = nrow(eta))),
                            nrow = nrow(eta))
            eta <- eta + noise
        }
        draws <- eta
    } else {
        draws <- if (type == "response") stats::plogis(eta) else eta
    }

    fit <- rowMeans(draws)
    names(fit) <- rownames(newdata)

    if (interval == "none") return(fit)

    probs <- c((1 - level) / 2, 1 - (1 - level) / 2)
    qs <- t(apply(draws, 1L, stats::quantile, probs = probs))

    out <- cbind(fit = fit, lwr = qs[, 1L], upr = qs[, 2L])
    rownames(out) <- rownames(newdata)
    out
}

#' Diagnostic and summary plots for a Bayesian Repulsive Ensemble
#'
#' @param x An object of class \code{"bre"}.
#' @param type Character. \code{"emip"} ranks features by ensemble marginal
#'   inclusion probability; \code{"coallocation"} shows the label-invariant
#'   probability that two features share a sub-model; \code{"allocation"} shows
#'   label-aligned per-sub-model ownership; \code{"intervals"} is a caterpillar
#'   plot of ensemble coefficients on the original scale; \code{"trace"} shows
#'   MCMC traces for the global quantities.
#' @param top Integer. Number of top-ranked features to display. Default 15.
#' @param level Numeric. Credible interval level for \code{type = "intervals"}.
#' @param ... Unused, for method consistency.
#'
#' @return Invisibly \code{NULL}; called for the plot side effect.
#'
#' @export
plot.bre <- function(x, type = c("emip", "coallocation", "allocation", "intervals", "trace"),
                     top = 15, level = 0.95, ...) {

    type <- match.arg(type)

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar), add = TRUE)

    idxTop <- order(x$emip, decreasing = TRUE)[seq_len(min(top, x$p))]
    idx <- rev(idxTop)                   # so the leading feature sits on top
    yPos <- seq_along(idx)
    labels <- x$featureNames[idx]

    if (type == "emip") {

        vals <- x$emip[idx]
        graphics::par(mar = c(5, 7, 4, 2))
        graphics::plot(vals, yPos, pch = 21, bg = "darkred", col = "darkred", cex = 1.4,
                       xlim = c(0, 1), ylim = c(0.5, length(idx) + 0.5), yaxt = "n",
                       xlab = "Ensemble marginal inclusion probability", ylab = "",
                       main = sprintf("Top %d features", length(idx)), cex.lab = 1.1)
        graphics::segments(0, yPos, vals, yPos, col = "darkgray", lwd = 2)
        graphics::axis(2, at = yPos, labels = labels, las = 1)
        graphics::abline(v = 0.5, lty = 2)

    } else if (type == "coallocation") {

        cm <- coallocation(x, top = min(top, x$p))
        ord <- rev(seq_len(nrow(cm)))
        pal <- grDevices::colorRampPalette(c("white", "steelblue", "darkblue"))(100)

        graphics::par(mar = c(7, 7, 4, 2))
        graphics::image(x = seq_len(ncol(cm)), y = seq_len(nrow(cm)),
                        z = t(cm[ord, , drop = FALSE]), col = pal, zlim = c(0, 1),
                        xlab = "", ylab = "", axes = FALSE,
                        main = "Posterior co-allocation")
        graphics::axis(1, at = seq_len(ncol(cm)), labels = colnames(cm), las = 2, cex.axis = 0.8)
        graphics::axis(2, at = seq_len(nrow(cm)), labels = rownames(cm)[ord], las = 1, cex.axis = 0.8)
        graphics::box()
        graphics::title(sub = "P(two features share a sub-model); label-invariant", cex.sub = 0.8)

    } else if (type == "allocation") {

        al <- allocation(x)
        mat <- al$pip[idx, , drop = FALSE]
        rs <- rowSums(mat)
        mat <- mat / ifelse(rs == 0, 1, rs)      # relative ownership per feature

        pal <- grDevices::colorRampPalette(c("white", "steelblue", "darkblue"))(100)

        graphics::par(mar = c(6, 7, 4, 2))
        graphics::image(x = seq_len(x$numModels), y = seq_len(nrow(mat)), z = t(mat),
                        col = pal, zlim = c(0, 1), xlab = "Sub-model", ylab = "",
                        main = "Relative feature allocation", axes = FALSE)
        graphics::axis(1, at = seq_len(x$numModels), labels = colnames(al$pip), las = 2, cex.axis = 0.8)
        graphics::axis(2, at = seq_len(nrow(mat)), labels = labels, las = 1)
        graphics::box()
        graphics::title(sub = sprintf("Label-aligned; relabelling rate %.0f%%", 100 * al$switchRate),
                        cex.sub = 0.8)

    } else if (type == "intervals") {

        betaOrig <- x$ensembleBetaChain / x$scaleInfo$xScale
        probs <- c((1 - level) / 2, 1 - (1 - level) / 2)
        ci <- apply(betaOrig[idx, , drop = FALSE], 1L, stats::quantile, probs = probs)
        means <- rowMeans(betaOrig[idx, , drop = FALSE])

        graphics::par(mar = c(5, 7, 4, 2))
        graphics::plot(means, yPos, pch = 20, col = "darkblue", cex = 1.4,
                       xlim = range(c(ci, 0)), ylim = c(0.5, length(idx) + 0.5), yaxt = "n",
                       xlab = sprintf("Ensemble coefficient (%g%% CI)", 100 * level),
                       ylab = "", main = sprintf("Top %d features", length(idx)), cex.lab = 1.1)
        graphics::segments(ci[1L, ], yPos, ci[2L, ], yPos, col = "darkblue", lwd = 2)
        graphics::axis(2, at = yPos, labels = labels, las = 1)
        graphics::abline(v = 0, lty = 2, col = "darkred", lwd = 1.5)

    } else {

        nPanel <- 2L + as.integer(x$family == "gaussian")
        graphics::par(mfrow = c(nPanel, 1), mar = c(4, 4, 2, 1))

        graphics::plot(x$lambda2Chain, type = "l", col = "darkblue",
                       xlab = "Retained draw", ylab = expression(lambda[2]),
                       main = if (x$learnLambda2) "Repulsion strength" else "Repulsion strength (fixed)")

        graphics::matplot(t(x$modelSize), type = "l", lty = 1,
                          xlab = "Retained draw", ylab = "Active features",
                          main = "Sub-model sizes")

        if (x$family == "gaussian") {
            graphics::plot(x$sigmaSqChain, type = "l", col = "darkred",
                           xlab = "Retained draw", ylab = expression(sigma^2),
                           main = "Error variance")
        }
    }

    invisible(NULL)
}
