makeData <- function(n = 120, p = 12, seed = 1, collinear = TRUE) {
    set.seed(seed)
    X <- matrix(rnorm(n * p), n, p)
    if (collinear) X[, 2] <- X[, 1] + rnorm(n, sd = 0.1)
    colnames(X) <- paste0("V", seq_len(p))
    beta <- c(1.5, 1.5, -1.5, rep(0, p - 3))
    list(X = X, y = as.numeric(5 + X %*% beta + rnorm(n)), beta = beta)
}

test_that("bre is exported and returns a well formed object", {
    expect_true(is.function(bayesrep::bre))

    d <- makeData()
    fit <- bre(d$X, d$y, numModels = 3, iter = 400, burnin = 200,
               verbose = FALSE, seed = 1)

    expect_s3_class(fit, "bre")
    expect_identical(fit$nKeep, 200L)
    expect_identical(dim(fit$ensembleBetaChain), c(12L, 200L))
    expect_identical(dim(fit$pip), c(12L, 3L))
    expect_length(fit$emip, 12L)
    expect_true(all(fit$emip >= 0 & fit$emip <= 1))
    expect_true(is.finite(waic(fit)$waic))
})

test_that("column names are carried through every summary", {
    d <- makeData()
    fit <- bre(d$X, d$y, numModels = 2, iter = 300, burnin = 150,
               verbose = FALSE, seed = 1)

    expect_identical(rownames(fit$pip), colnames(d$X))
    expect_identical(names(fit$emip), colnames(d$X))
    expect_identical(names(coef(fit)), c("(Intercept)", colnames(d$X)))
    expect_identical(as.character(summary(fit)$Feature[1:3]) %in% colnames(d$X),
                     rep(TRUE, 3))
    expect_identical(rownames(coallocation(fit, top = 3)),
                     colnames(d$X)[order(fit$emip, decreasing = TRUE)][1:3])
})

test_that("input validation rejects malformed arguments", {
    d <- makeData()

    expect_error(bre(d$X, d$y[-1], verbose = FALSE), "match the length")
    expect_error(bre(d$X, d$y, iter = 100, burnin = 100, verbose = FALSE),
                 "smaller than")
    expect_error(bre(d$X, round(abs(d$y)), family = "binomial", verbose = FALSE),
                 "only 0s and 1s")
    expect_error(bre(cbind(d$X, 1), d$y, verbose = FALSE), "constant")
    expect_error(bre(d$X, d$y, tauSq = 0, verbose = FALSE), "tauSq")
    expect_error(bre(d$X, d$y, lambda2 = -1, verbose = FALSE), "non-negative")
    expect_error(bre(d$X, d$y, iter = 2000, burnin = 0, maxMemoryGb = 1e-7,
                     verbose = FALSE), "maxMemoryGb")

    Xna <- d$X; Xna[1, 1] <- NA
    expect_error(bre(Xna, d$y, verbose = FALSE), "missing values")
})

test_that("the Gaussian fit recovers the intercept and the signal", {
    d <- makeData(n = 200, seed = 4)
    fit <- bre(d$X, d$y, numModels = 3, iter = 1500, burnin = 750,
               verbose = FALSE, seed = 2)

    cf <- coef(fit)
    expect_equal(unname(cf[1]), 5, tolerance = 0.3)          # true intercept 5
    expect_gt(fit$emip[["V3"]], 0.8)                         # isolated signal found
    expect_gt(max(fit$emip[c("V1", "V2")]), 0.8)             # collinear pair found
    expect_true(max(fit$emip[4:12]) < 0.5)                   # noise rejected

    # The collinear pair 1 and 2 share the signal; their sum is what is
    # identified, so check that rather than the individual coefficients.
    expect_equal(unname(cf["V1"] + cf["V2"]), 3, tolerance = 0.8)
    expect_equal(unname(cf["V3"]), -1.5, tolerance = 0.5)
})

test_that("coefficients are returned on the original data scale", {
    d <- makeData()
    fit <- bre(d$X, d$y, numModels = 2, iter = 400, burnin = 200,
               verbose = FALSE, seed = 1)

    scaled <- coef(fit, scaled = TRUE)
    orig <- coef(fit)

    expect_equal(unname(orig[-1]),
                 unname(unname(scaled[-1]) / unname(fit$scaleInfo$xScale)))
    expect_equal(unname(orig[1]),
                 unname(scaled[1]) + fit$scaleInfo$yCenter -
                     sum(unname(orig[-1]) * fit$scaleInfo$xCenter))

    # Coefficients must reproduce predictions on the original scale.
    expect_equal(as.numeric(orig[1] + d$X %*% orig[-1]),
                 unname(predict(fit, d$X)), tolerance = 1e-8)
})

test_that("the fit is stable as the ensemble size changes", {
    skip_on_cran()
    d <- makeData(n = 200, seed = 4)

    est <- vapply(c(1, 3, 6), function(g) {
        fit <- bre(d$X, d$y, numModels = g, iter = 1200, burnin = 600,
                   lambda2 = 2, verbose = FALSE, seed = 3)
        # V1 and V2 are collinear so only their sum is identified; V3 is not.
        c(unname(coef(fit)["V1"] + coef(fit)["V2"]), unname(coef(fit)["V3"]))
    }, numeric(2))

    # Because every sub-model is fitted to the response in full, tauSq is on the
    # scale of an ordinary regression coefficient and needs no G adjustment.
    expect_true(all(abs(est[1, ] - est[1, 1]) < 0.5))
    expect_true(all(abs(est[2, ] - est[2, 1]) < 0.5))
    expect_true(all(abs(est[2, ]) > 0.8))
})

test_that("the learning rate controls how much evidence each sub-model sees", {
    skip_on_cran()
    d <- makeData(n = 200, seed = 4)

    strong <- bre(d$X, d$y, numModels = 5, iter = 1000, burnin = 500,
                  learningRate = 1, verbose = FALSE, seed = 3)
    weak <- bre(d$X, d$y, numModels = 5, iter = 1000, burnin = 500,
                learningRate = 1 / 5, verbose = FALSE, seed = 3)

    # Tempering shrinks the evidence for inclusion, so fewer features clear the
    # sparsity threshold.
    expect_gte(sum(strong$emip > 0.5), sum(weak$emip > 0.5))
    expect_equal(strong$learningRate, 1)
    expect_equal(weak$learningRate, 1 / 5)

    expect_error(bre(d$X, d$y, learningRate = 0, verbose = FALSE), "learningRate")
    expect_error(bre(d$X, d$y, learningRate = 2, verbose = FALSE), "learningRate")
})

test_that("predict supports posterior-averaged intervals", {
    d <- makeData()
    fit <- bre(d$X, d$y, numModels = 2, iter = 500, burnin = 250,
               verbose = FALSE, seed = 1)

    point <- predict(fit, d$X)
    expect_length(point, nrow(d$X))

    cred <- predict(fit, d$X, interval = "credible")
    expect_identical(colnames(cred), c("fit", "lwr", "upr"))
    expect_equal(unname(cred[, "fit"]), unname(point))
    expect_true(all(cred[, "lwr"] <= cred[, "fit"]))
    expect_true(all(cred[, "upr"] >= cred[, "fit"]))

    pred <- predict(fit, d$X, interval = "prediction")
    expect_true(mean(pred[, "upr"] - pred[, "lwr"]) >
                    mean(cred[, "upr"] - cred[, "lwr"]))

    expect_error(predict(fit, d$X[, 1:3], verbose = FALSE), "columns")
    expect_error(predict(fit), "required")
})

test_that("the binomial family recovers an imbalanced intercept", {
    skip_on_cran()
    set.seed(8)
    n <- 400; p <- 10
    X <- matrix(rnorm(n * p), n, p)
    colnames(X) <- paste0("V", seq_len(p))
    a0 <- -2.2
    y <- rbinom(n, 1, stats::plogis(a0 + as.numeric(X %*% c(1.5, -1.5, rep(0, p - 2)))))

    fit <- bre(X, y, family = "binomial", numModels = 3, iter = 1500,
               burnin = 750, verbose = FALSE, seed = 5)

    expect_equal(unname(coef(fit)[1]), a0, tolerance = 0.5)
    expect_true(all(fit$emip[1:2] > 0.8))

    probs <- predict(fit, X, type = "response")
    expect_true(all(probs >= 0 & probs <= 1))
    expect_equal(mean(probs), mean(y), tolerance = 0.05)

    # Averaging within draws is not the same as the plug-in probability.
    plugin <- stats::plogis(coef(fit)[1] + as.numeric(X %*% coef(fit)[-1]))
    expect_gt(max(abs(probs - plugin)), 1e-4)

    expect_error(predict(fit, X, interval = "prediction"), "Gaussian")
})

test_that("predicted probabilities respond to type", {
    set.seed(9)
    n <- 150; p <- 6
    X <- matrix(rnorm(n * p), n, p)
    y <- rbinom(n, 1, stats::plogis(as.numeric(X %*% c(1.5, rep(0, p - 1)))))
    fit <- bre(X, y, family = "binomial", numModels = 2, iter = 400,
               burnin = 200, verbose = FALSE, seed = 1)

    link <- predict(fit, X, type = "link")
    resp <- predict(fit, X, type = "response")
    expect_true(all(resp >= 0 & resp <= 1))
    expect_false(isTRUE(all.equal(unname(link), unname(resp))))
})
