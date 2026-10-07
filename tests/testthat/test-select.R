selData <- function(n = 120, p = 20, seed = 1) {
    set.seed(seed)
    X <- matrix(rnorm(n * p), n, p)
    colnames(X) <- paste0("V", seq_len(p))
    list(X = X,
         y = as.numeric(X %*% c(2, -2, 1.5, rep(0, p - 3)) + rnorm(n)),
         truth = 1:3)
}

test_that("selectFeatures recovers the true signal at a controlled FDR", {
    d <- selData()
    fit <- bre(d$X, d$y, numModels = 3, iter = 800, burnin = 400,
               verbose = FALSE, seed = 1)

    sel <- selectFeatures(fit, fdr = 0.1)

    expect_s3_class(sel, "data.frame")
    expect_identical(names(sel), c("Feature", "eMIP", "PostMean", "MeanNModels"))
    expect_true(all(diff(sel$eMIP) <= 1e-12))          # ordered by eMIP
    expect_setequal(sel$Feature, colnames(d$X)[d$truth])
    expect_lte(attr(sel, "estimatedFdr"), 0.1)
    expect_identical(attr(sel, "targetFdr"), 0.1)
})

test_that("a stricter FDR target never selects more features", {
    d <- selData()
    fit <- bre(d$X, d$y, numModels = 3, iter = 800, burnin = 400,
               verbose = FALSE, seed = 1)

    sizes <- vapply(c(0.01, 0.05, 0.1, 0.3, 0.5),
                    function(f) nrow(selectFeatures(fit, f)), numeric(1))
    expect_true(all(diff(sizes) >= 0))
})

test_that("the estimated FDR is the mean of 1 - eMIP over the selected set", {
    d <- selData()
    fit <- bre(d$X, d$y, numModels = 3, iter = 600, burnin = 300,
               verbose = FALSE, seed = 2)

    sel <- selectFeatures(fit, fdr = 0.2)
    if (nrow(sel) > 0L) {
        expect_equal(attr(sel, "estimatedFdr"), mean(1 - sel$eMIP))
    }
})

test_that("selectFeatures does not inherit the G-dependence of an eMIP cutoff", {
    skip_on_cran()

    # A null feature gets G chances to be carried by some sub-model, so a fixed
    # eMIP cutoff admits more nulls as G grows. The FDR rule must not.
    set.seed(5)
    n <- 120; p <- 60
    X <- matrix(rnorm(n * p), n, p)
    colnames(X) <- paste0("V", seq_len(p))
    y <- as.numeric(X %*% c(2.5, -2.5, rep(0, p - 2)) + rnorm(n))
    nulls <- 3:p

    cutoffFp <- integer(0); fdrFp <- integer(0)
    for (G in c(2, 8)) {
        f <- bre(X, y, numModels = G, iter = 1000, burnin = 500,
                 thetaShape1 = 1, thetaShape2 = 2, verbose = FALSE, seed = 3)
        cutoffFp <- c(cutoffFp, sum(f$emip[nulls] > 0.5))
        fdrFp <- c(fdrFp, sum(selectFeatures(f, 0.1)$Feature %in% colnames(X)[nulls]))
    }

    # The naive cutoff lets in more nulls at the larger G; the FDR rule holds.
    expect_gt(cutoffFp[2], cutoffFp[1])
    expect_lte(max(fdrFp), 0.25 * length(nulls))
})

test_that("selectFeatures validates its input", {
    d <- selData(n = 60, p = 10)
    fit <- bre(d$X, d$y, numModels = 2, iter = 300, burnin = 150,
               verbose = FALSE, seed = 1)

    expect_error(selectFeatures(fit, fdr = 0), "between 0 and 1")
    expect_error(selectFeatures(fit, fdr = 1), "between 0 and 1")
    expect_error(selectFeatures(fit, fdr = c(0.1, 0.2)), "single number")
    expect_error(selectFeatures(list()), "inherits")
})

test_that("print reports the FDR-controlled count rather than an eMIP cutoff", {
    d <- selData(n = 80, p = 10)
    fit <- bre(d$X, d$y, numModels = 3, iter = 400, burnin = 200,
               verbose = FALSE, seed = 1)
    expect_output(print(fit), "Selected \\(FDR 10%\\)")
})
