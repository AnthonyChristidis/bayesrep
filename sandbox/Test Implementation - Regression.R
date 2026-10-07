# ------------------------------------------------------------------------------
# Test Simulation for Bayesian Repulsive Ensembles (BRE) - Regression Case
# ------------------------------------------------------------------------------

rm(list = ls())

library(MASS)
library(bayesrep)

# 1. Synthetic data with two highly correlated blocks of true signal ----------
# Every member of each block carries a real effect, so a fit that simply kept
# one member per block would be wrong. This is the setting in which separating
# the block across sub-models is the right answer.
set.seed(0)
n <- 150
p <- 20
rho <- 0.95

Sigma <- diag(p)
Sigma[1:3, 1:3] <- rho          # block 1
Sigma[4:6, 4:6] <- rho          # block 2
diag(Sigma) <- 1

X <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
colnames(X) <- paste0("V", seq_len(p))

trueBeta <- rep(0, p)
trueBeta[1:3] <- 1.5
trueBeta[4:6] <- -1.5
trueIntercept <- 7

y <- as.numeric(trueIntercept + X %*% trueBeta + rnorm(n, sd = 1))

# 2. Choose the repulsion strength by cross-validation ------------------------
# lambda2 cannot be learned from the posterior and must not be chosen by WAIC:
# the composite objective scores each sub-model on its own fit, and repulsion
# always makes that worse. Held-out error of the ensemble is the right target.
cv <- cv.bre(X, y,
             lambda2   = c(0, 2, 4, 8, 16, 32),
             numModels = 5,
             nfolds    = 5,
             iter      = 2000,
             burnin    = 1000,
             tauSq     = 5,
             seed      = 1)

print(cv)
fit <- cv$fit
print(fit)

# 3. Feature ranking ----------------------------------------------------------
print(summary(fit))

cat("\nTrue intercept:", trueIntercept,
    " estimated:", round(coef(fit)[1], 3), "\n")

# Within a 0.95-correlated block only the block total is well identified.
cat("Block 1 total  true", sum(trueBeta[1:3]),
    " estimated", round(sum(coef(fit)[2:4]), 3), "\n")
cat("Block 2 total  true", sum(trueBeta[4:6]),
    " estimated", round(sum(coef(fit)[5:7]), 3), "\n")

# 4. Allocation structure -----------------------------------------------------
# Under the composite likelihood the assignment of features to sub-models IS
# identified, so this matrix is informative. Read it against the 1/G baseline
# that random assignment would give.
cm <- coallocation(fit, features = paste0("V", 1:6))
cat("\nCo-allocation among the true signal (baseline 1/G =",
    round(1 / fit$numModels, 2), "):\n")
print(round(cm, 2))

blk <- rep(1:2, each = 3)
cat(sprintf("\nmean within-block : %.3f   (below baseline = block is being split)\n",
            mean(cm[outer(blk, blk, "==") & upper.tri(cm)])))
cat(sprintf("mean across-block : %.3f   (above baseline = used together)\n",
            mean(cm[outer(blk, blk, "!=")])))

# Essential check: low co-allocation only means separation if both features are
# actually in the model. The diagonal holds the inclusion probabilities.
cat("\neMIP of the signal features (want all high):\n")
print(round(diag(cm), 3))
cat("E[number of sub-models] per signal feature:\n")
print(round(rowSums(fit$pip)[1:6], 2))

# 5. Matched ablation ---------------------------------------------------------
noRep <- bre(X, y, numModels = 5, lambda2 = 0, iter = 2000, burnin = 1000,
             tauSq = 5, verbose = FALSE, seed = 1)
cmNo <- coallocation(noRep, features = paste0("V", 1:6))

cat(sprintf("\nwithin-block co-allocation  with repulsion: %.3f   without: %.3f\n",
            mean(cm[outer(blk, blk, "==") & upper.tri(cm)]),
            mean(cmNo[outer(blk, blk, "==") & upper.tri(cmNo)])))

# 6. Prediction with uncertainty ----------------------------------------------
pred <- predict(fit, X, interval = "credible")
cat("\nIn-sample RMSE:", round(sqrt(mean((pred[, "fit"] - y)^2)), 3),
    " (true sigma = 1)\n")

predInt <- predict(fit, X, interval = "prediction")
cat("95% prediction interval coverage of y:",
    round(mean(y >= predInt[, "lwr"] & y <= predInt[, "upr"]), 3),
    " (nominal 0.95)\n")

# 7. Plots --------------------------------------------------------------------
oldpar <- par(no.readonly = TRUE)
par(mfrow = c(2, 2))
plot(fit, type = "emip", top = 10)
plot(fit, type = "coallocation", top = 8)
plot(fit, type = "intervals", top = 10)
plot(cv)
par(oldpar)
