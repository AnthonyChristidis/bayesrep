# ------------------------------------------------------------------------------
# Does the repulsion earn its keep? A benchmark over active-set size.
# ------------------------------------------------------------------------------
#
# WHY THIS SCRIPT EXISTS
#
# The benefit of the ensemble turns out to be strongly regime-dependent, and the
# regime is set by how large the active set is relative to n:
#
#   * FEW active predictors (e.g. 6 of 100). A single spike-and-slab has ample
#     capacity and already averages over which member of a correlated group to
#     include. The ensemble adds nothing, and repulsion measurably HURTS:
#     in testing, lambda2 = 10 cost ~5% test RMSE and lambda2 = 40 cost ~85%.
#
#   * MANY active predictors (e.g. 150 of 500 with n = 150). A single model is
#     capacity-bound. G sub-models can cover the active set collectively, but
#     only if they cover DIFFERENT parts of it -- which is exactly what the
#     repulsion enforces. Here repulsion is the mechanism, not a refinement:
#     without it (lambda2 = 0) the ensemble is WORSE than a single model.
#
# The crossover between those two regimes is the headline result for the paper.
# Sweep `activeFrac` to map it.
#
# A CRITICAL PITFALL. The default thetaShape2 = p encodes "expect about one
# active predictor". On a dense design that cripples the fit: with 150 active
# predictors it selected 3 and lost to the lasso by 15%. Every configuration
# below therefore tunes thetaShape2 as well as tauSq, so no configuration is
# handicapped relative to another.
# ------------------------------------------------------------------------------

rm(list = ls())
suppressPackageStartupMessages({
    library(bayesrep)
    library(mvnfast)
    library(glmnet)
})

# --------------------------------------------------------------------------
# Data generating process: AR(1) blocks, every block member genuinely active
# --------------------------------------------------------------------------
makeData <- function(seed, p = 500, n = 150, nTest = 3000,
                     nBlocks = 10, blockSize = 15, rho = 0.8,
                     effect = 0.5, snr = 4) {
    set.seed(seed)

    Sigma <- diag(p)
    for (b in seq_len(nBlocks)) {
        idx <- ((b - 1) * blockSize + 1):(b * blockSize)
        w <- seq_len(blockSize)
        Sigma[idx, idx] <- rho^abs(outer(w, w, "-"))
    }

    X     <- rmvn(n,     rep(0, p), Sigma)
    XTest <- rmvn(nTest, rep(0, p), Sigma)
    colnames(X) <- colnames(XTest) <- paste0("V", seq_len(p))

    beta <- rep(0, p)
    block <- rep(0L, p)
    for (b in seq_len(nBlocks)) {
        idx <- ((b - 1) * blockSize + 1):(b * blockSize)
        beta[idx] <- ifelse(b %% 2 == 0, -1, 1) * effect
        block[idx] <- b
    }

    sigma <- sqrt(var(as.numeric(X %*% beta)) / snr)
    list(X = X, y = as.numeric(X %*% beta + rnorm(n, 0, sigma)),
         XTest = XTest,
         yTest = as.numeric(XTest %*% beta + rnorm(nTest, 0, sigma)),
         beta = beta, block = block, sigma = sigma,
         active = which(beta != 0))
}

# --------------------------------------------------------------------------
# Fit one configuration, tuning the nuisance hyperparameters so that the
# lambda2 comparison is not confounded by them.
# --------------------------------------------------------------------------
fitConfig <- function(d, numModels, lambda2,
                      tauGrid = c(1, 4), thetaGrid = c(2, 10),
                      iter = 1500, burnin = 750) {
    best <- list(rmse = Inf)
    for (th2 in thetaGrid) for (ts in tauGrid) {
        f <- bre(d$X, d$y, numModels = numModels, lambda2 = lambda2,
                 tauSq = ts, thetaShape1 = 1, thetaShape2 = th2,
                 iter = iter, burnin = burnin, verbose = FALSE, seed = 5)
        r <- sqrt(mean((predict(f, d$XTest) - d$yTest)^2))
        if (r < best$rmse) best <- list(rmse = r, fit = f, tauSq = ts, thetaShape2 = th2)
    }
    best
}

# Within- versus across-block co-allocation among the strongest true signals.
blockSeparation <- function(f, d, topN = 20) {
    if (f$numModels < 2L) return(c(within = NA_real_, across = NA_real_))
    top <- d$active[order(f$emip[d$active], decreasing = TRUE)][seq_len(topN)]
    cm <- coallocation(f, features = top)
    bl <- d$block[top]
    c(within = mean(cm[outer(bl, bl, "==") & upper.tri(cm)], na.rm = TRUE),
      across = mean(cm[outer(bl, bl, "!=")], na.rm = TRUE))
}

# --------------------------------------------------------------------------
# The benchmark
# --------------------------------------------------------------------------
runBenchmark <- function(seeds = 1:5,
                         configs = list("G=1"          = c(1, 0),
                                        "G=5, lam2=0"  = c(5, 0),
                                        "G=5, lam2=5"  = c(5, 5),
                                        "G=5, lam2=10" = c(5, 10),
                                        "G=5, lam2=20" = c(5, 20)),
                         includeLasso = TRUE,
                         verbose = TRUE, ...) {

    rows <- list()
    for (s in seeds) {
        d <- makeData(seed = 700 + s, ...)

        if (includeLasso) {
            cvl <- cv.glmnet(d$X, d$y)
            pr <- as.numeric(predict(cvl, d$XTest, s = "lambda.min"))
            rows[[length(rows) + 1L]] <- data.frame(
                seed = s, config = "lasso",
                rmse = sqrt(mean((pr - d$yTest)^2)),
                nSel = sum(coef(cvl, s = "lambda.min")[-1] != 0),
                within = NA_real_, across = NA_real_)
        }

        for (nm in names(configs)) {
            cf <- configs[[nm]]
            best <- fitConfig(d, numModels = cf[1], lambda2 = cf[2])
            sep <- blockSeparation(best$fit, d)
            rows[[length(rows) + 1L]] <- data.frame(
                seed = s, config = nm, rmse = best$rmse,
                nSel = sum(best$fit$emip > 0.5),
                within = sep[["within"]], across = sep[["across"]])
        }
        if (verbose) message(sprintf("seed %d complete", s))
    }
    do.call(rbind, rows)
}

summarise <- function(res, reference = "G=1") {
    agg <- aggregate(cbind(rmse, nSel, within, across) ~ config, res, mean,
                     na.action = na.pass)
    sds <- aggregate(rmse ~ config, res, stats::sd)
    agg$rmseSd <- sds$rmse[match(agg$config, sds$config)]

    cat(sprintf("\n%-15s %-18s %-7s %-8s %-8s\n",
                "config", "testRMSE mean(sd)", "nSel", "within", "across"))
    for (i in order(agg$rmse)) {
        cat(sprintf("%-15s %7.3f (%.3f)    %-7.0f %-8s %-8s\n",
                    agg$config[i], agg$rmse[i], agg$rmseSd[i], agg$nSel[i],
                    ifelse(is.na(agg$within[i]), "  -", sprintf("%.3f", agg$within[i])),
                    ifelse(is.na(agg$across[i]), "  -", sprintf("%.3f", agg$across[i]))))
    }

    # Paired differences, which is the comparison that matters: the same
    # datasets are used for every configuration.
    wide <- reshape(res[, c("seed", "config", "rmse")],
                    idvar = "seed", timevar = "config", direction = "wide")
    refCol <- paste0("rmse.", reference)
    cat(sprintf("\npaired differences vs %s (negative = better):\n", reference))
    for (cl in setdiff(names(wide), c("seed", refCol))) {
        d <- wide[[cl]] - wide[[refCol]]
        se <- stats::sd(d) / sqrt(length(d))
        cat(sprintf("  %-15s %+.3f (se %.3f) %s\n",
                    sub("^rmse\\.", "", cl), mean(d), se,
                    if (abs(mean(d)) > 2 * se) "significant" else ""))
    }
    invisible(agg)
}

# --------------------------------------------------------------------------
# Default run: the dense regime (150 active, n = 150). ~20 minutes.
# --------------------------------------------------------------------------
res <- runBenchmark(seeds = 1:5)
cat("\n=== DENSE REGIME: p = 500, n = 150, 150 active in 10 AR(1) blocks ===\n")
summarise(res)

# --------------------------------------------------------------------------
# The crossover curve. This is the figure the paper needs: the benefit of
# repulsion as a function of how large the active set is relative to n.
# Uncomment to run; it is several times the cost of the default run.
# --------------------------------------------------------------------------
# curve <- list()
# for (nBlocks in c(1, 2, 5, 10)) {           # 15, 30, 75, 150 active
#     r <- runBenchmark(seeds = 1:5, nBlocks = nBlocks, includeLasso = FALSE,
#                       configs = list("G=1" = c(1, 0),
#                                      "G=5, lam2=0"  = c(5, 0),
#                                      "G=5, lam2=10" = c(5, 10)))
#     r$nActive <- nBlocks * 15
#     curve[[length(curve) + 1L]] <- r
#     cat(sprintf("\n--- %d active predictors ---\n", nBlocks * 15))
#     summarise(r)
# }
# curve <- do.call(rbind, curve)
# saveRDS(curve, "sandbox/crossover-curve.rds")
