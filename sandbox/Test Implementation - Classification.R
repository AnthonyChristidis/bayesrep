# ------------------------------------------------------------------------------
# Test Simulation for Bayesian Repulsive Ensembles (BRE) - Classification Case
# ------------------------------------------------------------------------------

rm(list = ls())

library(MASS)
library(bayesrep)

# 1. Synthetic binary data with correlated blocks and class imbalance ---------
set.seed(0)
n <- 400
p <- 20
rho <- 0.95

Sigma <- diag(p)
Sigma[1:3, 1:3] <- rho
Sigma[4:6, 4:6] <- rho
diag(Sigma) <- 1

X <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
colnames(X) <- paste0("V", seq_len(p))

trueBeta <- rep(0, p)
trueBeta[1:3] <- 1.0
trueBeta[4:6] <- -1.0
trueIntercept <- -2.2          # deliberately imbalanced

trueProbs <- plogis(trueIntercept + as.numeric(X %*% trueBeta))
y <- rbinom(n, size = 1, prob = trueProbs)

cat("Class balance:", round(mean(y), 3), "\n")

# 2. Fit ----------------------------------------------------------------------
# The binomial engine now carries an explicit intercept. Without one the
# log-odds are forced through zero at x = 0, which biases every slope whenever
# the classes are not balanced.
fit <- bre(X, y,
           family    = "binomial",
           numModels = 3,
           iter      = 4000,
           burnin    = 2000,
           tauSq     = 1,
           lambda2   = 8,          # tune with cv.bre() for your own data
           seed      = 1)

print(fit)
print(summary(fit))

cat("\nTrue intercept:", trueIntercept,
    " estimated:", round(coef(fit)[1], 3), "\n")
cat("Block 1 total  true", sum(trueBeta[1:3]),
    " estimated", round(sum(coef(fit)[2:4]), 3), "\n")
cat("Block 2 total  true", sum(trueBeta[4:6]),
    " estimated", round(sum(coef(fit)[5:7]), 3), "\n")

# 3. Predictive performance ---------------------------------------------------
# Probabilities are averaged within draws, so these are posterior mean
# probabilities rather than the probability implied by the mean coefficient.
probs <- predict(fit, X, type = "response")

auc <- function(p, y) {
    r <- rank(p)
    (sum(r[y == 1]) - sum(y) * (sum(y) + 1) / 2) / (sum(y) * sum(1 - y))
}

cat("\nIn-sample AUC :", round(auc(probs, y), 3), "\n")
cat("Brier score   :", round(mean((probs - y)^2), 3), "\n")
cat("Calibration   : mean predicted", round(mean(probs), 3),
    "vs observed", round(mean(y), 3), "\n")

plugin <- plogis(coef(fit)[1] + as.numeric(X %*% coef(fit)[-1]))
cat("Max |posterior mean - plug-in| probability:",
    round(max(abs(probs - plugin)), 4), "\n")

# Patient-level credible intervals on the probability scale.
ci <- predict(fit, X, type = "response", interval = "credible")
cat("\nFirst five subjects (probability with 95% CI):\n")
print(round(head(ci, 5), 3))

# 4. Allocation structure -----------------------------------------------------
# The assignment of features to sub-models is identified under the composite
# likelihood, so read this against the 1/G baseline: below means the ensemble is
# separating the two features, above means it uses them together. Check the
# diagonal (the eMIPs) as well, because a feature that was simply dropped also
# shows a low co-allocation without that meaning anything.
cm <- coallocation(fit, features = paste0("V", 1:6))
cat("\nCo-allocation among the true signal (baseline 1/G =",
    round(1 / fit$numModels, 2), "):\n")
print(round(cm, 2))

blk <- rep(1:2, each = 3)
cat(sprintf("\nmean within-block : %.3f\nmean across-block : %.3f\n",
            mean(cm[outer(blk, blk, "==") & upper.tri(cm)]),
            mean(cm[outer(blk, blk, "!=")])))
cat("eMIP of the signal features:\n"); print(round(diag(cm), 3))

# 5. Plots --------------------------------------------------------------------
oldpar <- par(no.readonly = TRUE)
par(mfrow = c(2, 2))
plot(fit, type = "emip", top = 10)
plot(fit, type = "coallocation", top = 8)
plot(fit, type = "intervals", top = 10)
plot(fit, type = "trace")
par(oldpar)
