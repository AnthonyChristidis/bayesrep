# Correctness of the individual MCMC kernels, checked against quantities that
# can be computed exactly without the sampler.

logIsingZ <- function(theta, lambda2, G) {
    m <- 0:G
    sum(choose(G, m) * theta^m * (1 - theta)^(G - m) * exp(-lambda2 * m * (m - 1) / 2))
}

test_that("the theta update targets the exact Ising-corrected posterior", {
    skip_on_cran()

    G <- 4; lambda2 <- 2; a <- 1; b <- 1; mj <- 2

    # Exact posterior: Beta kernel divided by the Ising normalising constant.
    kernel <- function(t) {
        t^(a + mj - 1) * (1 - t)^(b + G - mj - 1) /
            vapply(t, logIsingZ, numeric(1), lambda2 = lambda2, G = G)
    }
    norm <- stats::integrate(kernel, 0, 1)$value
    exactMean <- stats::integrate(function(t) t * kernel(t) / norm, 0, 1)$value

    set.seed(99)
    nRep <- 50000
    theta <- rep(0.5, nRep)
    for (s in 1:400) {
        theta <- as.numeric(updateThetaCpp(theta, rep(mj, nRep), a, b, lambda2, G)$thetaVec)
    }

    expect_equal(mean(theta), exactMean, tolerance = 0.01)

    # The naive conjugate draw, which the Ising constant corrects, is badly off.
    naiveMean <- (a + mj) / (a + b + G)
    expect_gt(abs(naiveMean - exactMean), 0.2)
})

test_that("the theta update is conjugate and always accepted when lambda2 is 0", {
    set.seed(7)
    res <- updateThetaCpp(rep(0.5, 500), rep(1, 500), 1, 1, 0, 3)
    expect_identical(res$numAccepted, res$numProposed)
    expect_equal(mean(res$thetaVec), 2 / 5, tolerance = 0.05)
})

test_that("the composite likelihood identifies which sub-model carries a feature", {
    # This is the property the averaged-likelihood formulation lacks: there,
    # relocating a feature left the posterior density exactly unchanged.
    set.seed(3)
    n <- 40; p <- 5; G <- 3; tauSq <- 1; lambda2 <- 2.5; sigma <- 0.7
    eta <- 1 / G
    X <- matrix(rnorm(n * p), n, p)
    y <- as.numeric(X %*% c(1, -1, 0, 0, 0) + rnorm(n, sd = sigma))
    theta <- c(0.3, 0.4, 0.2, 0.25, 0.35)

    logPost <- function(B, Gam) {
        m <- rowSums(Gam)
        ssr <- sum(vapply(seq_len(G), function(g) {
            r <- y - as.numeric(X %*% B[, g]); sum(r^2)
        }, numeric(1)))
        composite <- -eta * ssr / (2 * sigma^2)
        slab <- sum(dnorm(B[Gam == 1], 0, sqrt(tauSq), log = TRUE))
        ising <- sum(m * log(theta) + (G - m) * log(1 - theta) -
                         lambda2 * m * (m - 1) / 2 -
                         log(vapply(theta, logIsingZ, numeric(1),
                                    lambda2 = lambda2, G = G)))
        composite + slab + ising
    }

    B <- matrix(0, p, G); Gam <- matrix(0, p, G)
    Gam[1, 1] <- 1; B[1, 1] <- 1.1
    Gam[2, 2] <- 1; B[2, 2] <- -0.9
    reference <- logPost(B, Gam)

    # Moving feature 1 into the sub-model that already carries feature 2 changes
    # that sub-model's fit, so the density moves.
    B2 <- B; G2 <- Gam
    B2[1, 1] <- 0; G2[1, 1] <- 0
    B2[1, 2] <- 1.1; G2[1, 2] <- 1
    expect_false(isTRUE(all.equal(logPost(B2, G2), reference)))

    # Collecting both features into one sub-model likewise changes it.
    B3 <- B; G3 <- Gam
    B3[2, 2] <- 0; G3[2, 2] <- 0
    B3[2, 1] <- -0.9; G3[2, 1] <- 1
    expect_false(isTRUE(all.equal(logPost(B3, G3), reference)))

    # The remaining symmetry is the global one: permuting all G labels at once.
    perm <- c(2, 3, 1)
    expect_equal(logPost(B[, perm], Gam[, perm]), reference)
})

test_that("the transfer move preserves sparsity bookkeeping", {
    set.seed(11)
    n <- 60; p <- 8; G <- 3
    X <- matrix(rnorm(n * p), n, p)
    y <- as.numeric(X %*% c(2, -2, rep(0, p - 2)) + rnorm(n))

    beta <- matrix(0, p, G); gamma <- matrix(0, p, G)
    gamma[1, 1] <- 1; beta[1, 1] <- 2
    gamma[2, 2] <- 1; beta[2, 2] <- -2
    gamma[3, 3] <- 1; beta[3, 3] <- 0.3

    res <- transferFeatureCpp(resp = matrix(y, n, G),
                              sqrtW = matrix(1, n, G),
                              X = X, betaMatrix = beta, gammaMatrix = gamma,
                              sigmaSqScaled = G, tauSq = 1, numMoves = 50)

    # Each feature is still carried exactly as often as before.
    expect_equal(rowSums(res$gammaMatrix), rowSums(gamma))
    # Zero cells stay zero and active cells stay non-zero.
    expect_true(all(res$betaMatrix[res$gammaMatrix == 0] == 0))
    expect_true(all(res$betaMatrix[res$gammaMatrix == 1] != 0))
    expect_true(res$numAttempted > 0)
    expect_true(res$numAccepted <= res$numAttempted)
})

test_that("the transfer move does nothing with a single sub-model", {
    set.seed(2)
    X <- matrix(rnorm(20 * 3), 20, 3)
    beta <- matrix(c(1, 0, 0), 3, 1); gamma <- matrix(c(1, 0, 0), 3, 1)
    res <- transferFeatureCpp(matrix(rnorm(20), 20, 1), matrix(1, 20, 1),
                              X, beta, gamma, 1, 1, 10)
    expect_equal(res$gammaMatrix, gamma)
    expect_identical(res$numAttempted, 0L)
})

test_that("the Polya-Gamma sampler matches its exact mean at non-integer shape", {
    skip_on_cran()

    exactMean <- function(b, cc) if (cc == 0) b / 4 else (b / (2 * cc)) * tanh(cc / 2)

    set.seed(5)
    for (cfg in list(c(1 / 5, 0), c(1 / 3, 1), c(1 / 2, 3), c(1, 2))) {
        b <- cfg[1]; cc <- cfg[2]
        draws <- rpgVecCpp(b, rep(cc, 200000), 50)
        se <- stats::sd(draws) / sqrt(length(draws))
        expect_lt(abs(mean(draws) - exactMean(b, cc)), 4 * se)
        expect_true(all(draws > 0))
    }
})

test_that("the Polya-Gamma sampler agrees with pgdraw at shape one", {
    skip_on_cran()
    skip_if_not_installed("pgdraw")

    set.seed(6)
    mine <- rpgVecCpp(1, rep(2, 100000), 50)
    theirs <- pgdraw::pgdraw(1, rep(2, 100000))

    expect_equal(mean(mine), mean(theirs), tolerance = 0.01)
    expect_equal(stats::sd(mine), stats::sd(theirs), tolerance = 0.03)
})
