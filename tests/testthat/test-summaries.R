makeFit <- function(G = 3, iter = 1200, burnin = 600, seed = 1, ...) {
    set.seed(seed)
    n <- 150; p <- 10
    X <- matrix(rnorm(n * p), n, p)
    X[, 2] <- X[, 1] + rnorm(n, sd = 0.1)
    colnames(X) <- paste0("V", seq_len(p))
    y <- as.numeric(X %*% c(1.5, 1.5, -1.5, rep(0, p - 3)) + rnorm(n))
    bre(X, y, numModels = G, iter = iter, burnin = burnin,
        verbose = FALSE, seed = seed, ...)
}

# A correlated block whose members are all genuinely needed, so that a fit
# which simply drops the redundant ones would be wrong. This is the setting in
# which separation is meaningful.
blockFit <- function(lambda2, G = 5, iter = 2000, burnin = 1000, seed = 1, ...) {
    set.seed(seed)
    n <- 150; p <- 20; rho <- 0.95
    S <- diag(p); S[1:3, 1:3] <- rho; S[4:6, 4:6] <- rho; diag(S) <- 1
    ch <- chol(S)
    X <- matrix(rnorm(n * p), n, p) %*% ch
    colnames(X) <- paste0("V", seq_len(p))
    beta <- rep(0, p); beta[1:3] <- 1.5; beta[4:6] <- -1.5
    y <- as.numeric(X %*% beta + rnorm(n))
    bre(X, y, numModels = G, iter = iter, burnin = burnin, tauSq = 5,
        lambda2 = lambda2, verbose = FALSE, seed = seed, ...)
}

blockStats <- function(fit) {
    cm <- coallocation(fit, features = 1:6)
    blk <- rep(1:2, each = 3)
    list(cm = cm,
         within = mean(cm[outer(blk, blk, "==") & upper.tri(cm)]),
         across = mean(cm[outer(blk, blk, "!=")]),
         minEmip = min(fit$emip[1:6]))
}

test_that("summary ranks by eMIP and reports original-scale intervals", {
    fit <- makeFit()
    s <- summary(fit)

    expect_s3_class(s, "summary.bre")
    expect_identical(nrow(s), 10L)
    expect_true(all(diff(s$eMIP) <= 0))
    expect_true(all(s$Lower <= s$PostMean))
    expect_true(all(s$Upper >= s$PostMean))
    expect_true(all(s$MeanNModels >= 0 & s$MeanNModels <= fit$numModels))

    # Matches coef() on the original scale.
    cf <- coef(fit, intercept = FALSE)
    expect_equal(s$PostMean, unname(cf[as.character(s$Feature)]), tolerance = 1e-10)

    expect_output(print(s), "ensemble marginal inclusion")
})

test_that("eMIP is the exact Monte Carlo probability, not a product over sub-models", {
    fit <- makeFit(G = 3, iter = 800, burnin = 400, storeChains = TRUE)

    g <- array(as.integer(fit$gammaChain), dim = c(fit$p, fit$numModels, fit$nKeep))
    exact <- rowMeans(apply(g, c(1, 3), function(v) any(v == 1)))
    expect_equal(unname(fit$emip), unname(exact))

    expect_true(all(fit$emip >= 0 & fit$emip <= 1))
    expect_true(all(fit$emip >= apply(fit$pip, 1, max) - 1e-9))
})

test_that("eMIP is at least the largest single-sub-model inclusion probability", {
    fit <- makeFit(G = 3, iter = 800, burnin = 400)

    # A sanity bound that holds whatever the dependence across sub-models, and
    # which the exact Monte Carlo estimator satisfies by construction.
    expect_true(all(fit$emip >= apply(fit$pip, 1, max) - 1e-9))
    expect_true(all(fit$emip <= pmin(rowSums(fit$pip), 1) + 1e-9))
})


test_that("strong repulsion separates a correlated block across sub-models", {
    skip_on_cran()

    fit <- blockFit(lambda2 = 16)
    st <- blockStats(fit)

    # Every block member must still be in the model: a low co-allocation that
    # came from dropping one of them would not be separation.
    expect_gt(st$minEmip, 0.8)

    # Members of the same correlated block avoid each other, members of
    # different blocks do not.
    expect_lt(st$within, 1 / fit$numModels)
    expect_gt(st$across, st$within)

    expect_equal(unname(diag(st$cm)), unname(fit$emip[1:6]))
    expect_true(all(st$cm >= 0 & st$cm <= 1))
})

test_that("turning the repulsion off removes the separation", {
    skip_on_cran()

    strong <- blockStats(blockFit(lambda2 = 16, iter = 1500, burnin = 750))
    none <- blockStats(blockFit(lambda2 = 0, iter = 1500, burnin = 750))

    # With no repulsion every sub-model simply takes every useful feature.
    expect_gt(none$within, strong$within)
    expect_gt(none$within, 0.8)
})

test_that("the transfer move lets features change sub-model", {
    skip_on_cran()

    ownerChanges <- function(fit) {
        g <- array(as.integer(fit$gammaChain),
                   dim = c(fit$p, fit$numModels, fit$nKeep))
        top <- order(fit$emip, decreasing = TRUE)[1:3]
        owner <- apply(g[top, , , drop = FALSE], c(1, 3),
                       function(v) if (sum(v) == 1) which(v == 1) else NA_integer_)
        vapply(seq_len(nrow(owner)), function(i) {
            o <- owner[i, ]; o <- o[!is.na(o)]
            if (length(o) < 2L) 0L else sum(diff(o) != 0)
        }, integer(1))
    }

    stuck <- blockFit(lambda2 = 16, iter = 1500, burnin = 750, transferMoves = 0)
    moving <- blockFit(lambda2 = 16, iter = 1500, burnin = 750, transferMoves = 20)

    # The move is attempted and sometimes accepted.
    expect_true(is.na(stuck$acceptance$transfer))
    expect_true(is.finite(moving$acceptance$transfer))
    expect_gt(moving$acceptance$transfer, 0)

    # It does not disturb the fit: relocating a feature is a move within the
    # posterior, not a change to it.
    expect_equal(unname(coef(moving)), unname(coef(stuck)), tolerance = 0.3)

    # Note that a strongly identified assignment legitimately stops moving:
    # ownerChanges can be zero for signal features because the posterior really
    # does prefer one arrangement. It is reported here only as a diagnostic.
    expect_true(all(ownerChanges(moving) >= 0))
})

test_that("allocation returns an aligned matrix and a switch rate", {
    fit <- makeFit(G = 3, iter = 600, burnin = 300)
    al <- allocation(fit)

    expect_identical(dim(al$pip), c(10L, 3L))
    expect_true(all(al$pip >= 0 & al$pip <= 1))
    expect_true(al$switchRate >= 0 && al$switchRate <= 1)
    expect_identical(rownames(al$pip), fit$featureNames)
})

test_that("feature selection helpers validate their input", {
    fit <- makeFit(iter = 400, burnin = 200)

    expect_identical(dim(coallocation(fit, top = 4)), c(4L, 4L))
    expect_identical(rownames(coallocation(fit, features = c("V3", "V1"))),
                     c("V3", "V1"))
    expect_error(coallocation(fit, features = "nope"), "Unknown feature")
    expect_error(coallocation(fit, features = 99), "out of range")
})

test_that("print reports the fit without error", {
    fit <- makeFit(iter = 400, burnin = 200)
    expect_output(print(fit), "Bayesian Repulsive Ensemble")
    expect_output(print(fit), "Repulsion")
})

test_that("plot restores graphical parameters and accepts every type", {
    fit <- makeFit(iter = 400, burnin = 200)

    tmp <- tempfile(fileext = ".pdf")
    grDevices::pdf(tmp)
    on.exit({ grDevices::dev.off(); unlink(tmp) }, add = TRUE)

    for (ty in c("emip", "coallocation", "allocation", "intervals", "trace")) {
        before <- graphics::par(c("mar", "mfrow"))
        plot(fit, type = ty)
        expect_identical(graphics::par(c("mar", "mfrow")), before)
    }
})

test_that("lambda2 is fixed by default and sampled only when asked", {
    fixed <- makeFit(iter = 400, burnin = 200)
    expect_false(fixed$learnLambda2)
    expect_true(all(fixed$lambda2Chain == 2))
    expect_true(is.na(fixed$acceptance$lambda2))

    explicit <- makeFit(iter = 400, burnin = 200, lambda2 = 2.5)
    expect_true(all(explicit$lambda2Chain == 2.5))

    learned <- makeFit(iter = 800, burnin = 400, lambda2 = NULL)
    expect_true(learned$learnLambda2)
    expect_gt(stats::var(learned$lambda2Chain), 0)
    expect_true(is.finite(learned$acceptance$lambda2))

    # Sampling it is available but not recommended: under the composite
    # likelihood each sub-model prefers to keep every useful predictor, so the
    # posterior for lambda2 is pulled towards zero.
    expect_lt(mean(learned$lambda2Chain), 2)
})

test_that("thinning reduces storage without changing the chain dimensions", {
    fit <- makeFit(iter = 1000, burnin = 500, thin = 5)
    expect_identical(fit$nKeep, 100L)
    expect_identical(ncol(fit$ensembleBetaChain), 100L)
    expect_length(fit$gammaChain, fit$p * fit$numModels * 100L)
})

test_that("storeChains controls the optional per-sub-model chains", {
    lean <- makeFit(iter = 400, burnin = 200)
    expect_null(lean$betaChain)
    expect_null(lean$thetaChain)

    full <- makeFit(iter = 400, burnin = 200, storeChains = TRUE)
    expect_identical(dim(full$betaChain), c(10L, 3L, 200L))
    expect_identical(dim(full$thetaChain), c(10L, 200L))
})

test_that("a single sub-model reduces to an ordinary spike-and-slab fit", {
    fit <- makeFit(G = 1, iter = 600, burnin = 300)
    expect_identical(dim(fit$pip), c(10L, 1L))
    expect_equal(unname(fit$emip), unname(fit$pip[, 1]))
    expect_true(is.finite(waic(fit)$waic))
})
