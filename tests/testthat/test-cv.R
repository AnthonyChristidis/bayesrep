cvData <- function(n = 100, p = 10, seed = 1) {
    set.seed(seed)
    X <- matrix(rnorm(n * p), n, p)
    X[, 2] <- X[, 1] + rnorm(n, sd = 0.1)
    colnames(X) <- paste0("V", seq_len(p))
    list(X = X,
         y = as.numeric(X %*% c(1.5, 1.5, -2, rep(0, p - 3)) + rnorm(n)))
}

test_that("cv.bre scores a grid and refits at the chosen setting", {
    d <- cvData()
    cv <- cv.bre(d$X, d$y, lambda2 = c(0, 2), nfolds = 3,
                 iter = 300, burnin = 150, verbose = FALSE, seed = 1)

    expect_s3_class(cv, "cv.bre")
    expect_identical(nrow(cv$grid), 2L)
    expect_true(all(c("lambda2", "numModels", "learningRate",
                      "cvMean", "cvSe") %in% names(cv$grid)))
    expect_true(all(is.finite(cv$grid$cvMean)))
    expect_identical(dim(cv$errors), c(2L, 3L))

    expect_s3_class(cv$fit, "bre")
    expect_equal(cv$fit$lambda2Chain[1L], cv$chosen$lambda2)
    expect_identical(cv$fit$numModels, cv$chosen$numModels)

    expect_output(print(cv), "Cross-validated BRE selection")
})

test_that("cv.bre crosses lambda2, numModels and learningRate", {
    d <- cvData()
    cv <- cv.bre(d$X, d$y, lambda2 = c(0, 2), numModels = c(2, 3),
                 learningRate = c(0.5, 1), nfolds = 2,
                 iter = 200, burnin = 100, refit = FALSE,
                 verbose = FALSE, seed = 1)

    expect_identical(nrow(cv$grid), 8L)
    expect_setequal(unique(cv$grid$learningRate), c(0.5, 1))
    expect_null(cv$fit)
})

test_that("the one-standard-error rule prefers the stronger repulsion", {
    d <- cvData()
    cv <- cv.bre(d$X, d$y, lambda2 = c(0, 1, 2), nfolds = 3,
                 iter = 300, burnin = 150, verbose = FALSE, seed = 2)

    expect_gte(cv$bestOneSe$lambda2, cv$best$lambda2)
    expect_identical(cv$chosen$lambda2, cv$bestOneSe$lambda2)

    cvBest <- cv.bre(d$X, d$y, lambda2 = c(0, 1, 2), nfolds = 3, oneSe = FALSE,
                     iter = 300, burnin = 150, refit = FALSE,
                     verbose = FALSE, seed = 2)
    expect_identical(cvBest$chosen$lambda2, cvBest$best$lambda2)
})

test_that("cv.bre handles the binomial family and its measures", {
    skip_on_cran()
    set.seed(5)
    n <- 150; p <- 8
    X <- matrix(rnorm(n * p), n, p)
    colnames(X) <- paste0("V", seq_len(p))
    y <- rbinom(n, 1, stats::plogis(as.numeric(X %*% c(1.5, -1.5, rep(0, p - 2)))))

    for (m in c("deviance", "auc", "misclass")) {
        cv <- cv.bre(X, y, family = "binomial", measure = m,
                     lambda2 = c(0, 2), nfolds = 3, iter = 250, burnin = 125,
                     refit = FALSE, verbose = FALSE, seed = 1)
        expect_true(all(is.finite(cv$grid$cvMean)))
        expect_identical(cv$measure, m)
    }

    expect_error(cv.bre(X, y, family = "binomial", measure = "mse",
                        verbose = FALSE), "deviance")
})

test_that("cv.bre validates its inputs", {
    d <- cvData(n = 60)

    expect_error(cv.bre(d$X, d$y[-1], verbose = FALSE), "match the length")
    expect_error(cv.bre(d$X, d$y, lambda2 = -1, verbose = FALSE), "non-negative")
    expect_error(cv.bre(d$X, d$y, learningRate = 3, verbose = FALSE), "learningRate")
    expect_error(cv.bre(d$X, d$y, nfolds = 1, verbose = FALSE), "nfolds")
    expect_error(cv.bre(d$X, d$y, foldid = 1:3, verbose = FALSE), "length n")
})

test_that("a user supplied foldid is honoured", {
    d <- cvData(n = 60)
    folds <- rep(1:3, length.out = 60)
    cv <- cv.bre(d$X, d$y, lambda2 = 2, foldid = folds,
                 iter = 200, burnin = 100, refit = FALSE,
                 verbose = FALSE, seed = 1)
    expect_identical(cv$foldid, folds)
    expect_identical(ncol(cv$errors), 3L)
})

test_that("the cross-validation curve plots without altering par", {
    d <- cvData(n = 60)
    cv <- cv.bre(d$X, d$y, lambda2 = c(0, 1, 2), nfolds = 2,
                 iter = 200, burnin = 100, refit = FALSE,
                 verbose = FALSE, seed = 1)

    tmp <- tempfile(fileext = ".pdf")
    grDevices::pdf(tmp)
    on.exit({ grDevices::dev.off(); unlink(tmp) }, add = TRUE)

    before <- graphics::par(c("mar", "mfrow"))
    plot(cv)
    expect_identical(graphics::par(c("mar", "mfrow")), before)
})
