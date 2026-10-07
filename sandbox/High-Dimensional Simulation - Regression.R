# ------------------------------------------------------------------------------
# High-Dimensional Test: 10 AR(1) Pathways, p = 500, n = 200
# Compares a single spike-and-slab (G = 1) against the ensemble (G = 5) on
# prediction, discovery and predictive uncertainty.
#
# NOTE ON INTERPRETATION. Under the composite likelihood the assignment of
# features to sub-models is identified, so the co-allocation structure reported
# below is meaningful. But lambda2 is a genuine tuning parameter: it cannot be
# learned from the posterior and must not be chosen by WAIC, both of which drive
# it to zero. The lambda2 = 0 fit below is the matched control that isolates
# what the repulsion contributes; for your own data choose lambda2 with cv.bre().
# ------------------------------------------------------------------------------

rm(list = ls())
suppressPackageStartupMessages({
    library(mvnfast)
    library(bayesrep)
})

# --------------------------------------------------------
# 1. Data generating process: AR(1) blocks
# --------------------------------------------------------
set.seed(0)
nTrain <- 200
nTest  <- 5000
p      <- 500
numBlocks <- 10
blockSize <- 15
rho <- 0.80

Sigma <- diag(p)
for (b in seq_len(numBlocks)) {
    idx <- ((b - 1) * blockSize + 1):(b * blockSize)
    within <- seq_len(blockSize)
    Sigma[idx, idx] <- rho^(abs(outer(within, within, "-")))
}

XTrain <- rmvn(nTrain, mu = rep(0, p), sigma = Sigma)
XTest  <- rmvn(nTest,  mu = rep(0, p), sigma = Sigma)
colnames(XTrain) <- colnames(XTest) <- paste0("V", seq_len(p))

trueBeta <- rep(0, p)
trueLabels <- rep(0, p)
for (b in seq_len(numBlocks)) {
    idx <- ((b - 1) * blockSize + 1):(b * blockSize)
    trueBeta[idx] <- ifelse(b %% 2 == 0, -1, 1) * 2.0
    trueLabels[idx] <- 1
}

signalTrain <- as.numeric(XTrain %*% trueBeta)
sigmaNoise <- sqrt(var(signalTrain) / 1.0)          # SNR = 1
trueNoiseVar <- sigmaNoise^2

yTrain <- signalTrain + rnorm(nTrain, 0, sigmaNoise)
yTest  <- as.numeric(XTest %*% trueBeta) + rnorm(nTest, 0, sigmaNoise)

cat(sprintf("Simulated n = %d, p = %d, %d AR(1) pathways, SNR = 1\n",
            nTrain, p, numBlocks))

# --------------------------------------------------------
# 2. Fit
# --------------------------------------------------------
# tauSq is the prior variance of a sub-model coefficient. Because every
# sub-model predicts the response in full, it needs no adjustment for G.
common <- list(X = XTrain, y = yTrain, family = "gaussian",
               iter = 4000, burnin = 2000, thin = 2,
               tauSq = 5, thetaShape1 = 1, thetaShape2 = 10,
               verbose = FALSE, seed = 1)

cat("\nFitting three models...\n")
fitBase <- do.call(bre, c(common, list(numModels = 1)))                     # G = 1
fitAbl  <- do.call(bre, c(common, list(numModels = 5, lambda2 = 0)))        # no repulsion
fitEns  <- do.call(bre, c(common, list(numModels = 5, lambda2 = 16)))       # repulsion on

models <- list("Spike-and-slab (G=1)" = fitBase,
               "Ensemble, lambda2=0"  = fitAbl,
               "Ensemble, lambda2=16" = fitEns)

# --------------------------------------------------------
# 3. Predictive accuracy and uncertainty
# --------------------------------------------------------
evalPredictive <- function(fit) {
    pr <- predict(fit, XTest, interval = "credible")
    list(relMse = mean((yTest - pr[, "fit"])^2) / trueNoiseVar,
         uqWidth = mean(pr[, "upr"] - pr[, "lwr"]),
         pred = pr)
}

# --------------------------------------------------------
# 4. Variable discovery, ranked by exact eMIP
# --------------------------------------------------------
evalDiscovery <- function(fit) {
    emip <- fit$emip
    o <- order(emip, decreasing = TRUE)
    lab <- trueLabels[o]
    tpr <- cumsum(lab) / sum(lab)
    fpr <- cumsum(1 - lab) / sum(1 - lab)
    auc <- sum(diff(fpr) * (tpr[-1] + tpr[-length(tpr)]) / 2)

    sel <- emip > 0.5
    list(auc = auc,
         tpr = sum(sel & trueLabels == 1) / sum(trueLabels == 1),
         fpr = sum(sel & trueLabels == 0) / sum(trueLabels == 0),
         selected = sum(sel))
}

pred <- lapply(models, evalPredictive)
disc <- lapply(models, evalDiscovery)

# --------------------------------------------------------
# 5. Report
# --------------------------------------------------------
cat("\n", strrep("-", 78), "\n", sep = "")
cat("  BRE PERFORMANCE REPORT: HIGH-DIMENSIONAL PATHWAYS\n")
cat(strrep("-", 78), "\n", sep = "")
cat(sprintf("True signal: %d features across %d AR(1) blocks\n",
            sum(trueLabels), numBlocks))
cat("Theoretical minimum scaled MSE = 1.000; perfect AUC = 1.000\n")

cat("\n--- A. VARIABLE SELECTION ---\n")
cat(sprintf("%-28s | %-6s | %-5s | %-6s | %-8s\n",
            "Model", "AUC", "TPR", "FPR", "Selected"))
cat(strrep("-", 78), "\n", sep = "")
for (nm in names(models)) {
    d <- disc[[nm]]
    cat(sprintf("%-28s | %-6.4f | %-5.2f | %-6.4f | %-8d\n",
                nm, d$auc, d$tpr, d$fpr, d$selected))
}

cat("\n--- B. PREDICTION AND UNCERTAINTY (test set) ---\n")
cat(sprintf("%-28s | %-12s | %-13s | %-8s\n",
            "Model", "Scaled MSE", "Avg CI width", "WAIC*"))
cat(strrep("-", 78), "\n", sep = "")
for (nm in names(models)) {
    cat(sprintf("%-28s | %-12.4f | %-13.4f | %-8.1f\n",
                nm, pred[[nm]]$relMse, pred[[nm]]$uqWidth, waic(models[[nm]])$waic))
}
cat("* WAIC measures composite in-sample fit and always prefers lambda2 = 0;\n")
cat("  it is shown for reference only and must not be used to select lambda2.\n")

cat("\n--- C. WHAT THE REPULSION IS DOING ---\n")
cat(sprintf("lambda2 held fixed at %.1f (choose it with cv.bre for your own data)\n",
            fitEns$lambda2Chain[1]))
cat(sprintf("E[m_j] over true signal : %.3f  (1 = no duplication, G = fully shared)\n",
            mean(rowSums(fitEns$pip)[trueLabels == 1])))
cat(sprintf("E[m_j] over noise       : %.3f\n",
            mean(rowSums(fitEns$pip)[trueLabels == 0])))
cat(sprintf("Co-allocation among top signal features: %.3f (combinatorial value %.3f)\n",
            {
                cm <- coallocation(fitEns, top = 10)
                mean(cm[upper.tri(cm)])
            }, 1 / fitEns$numModels))

cat("\n--- D. EXAMPLE SUBJECT PREDICTIONS (95% credible intervals) ---\n")
for (i in 1:3) {
    cat(sprintf("Subject %d (true y = %6.2f):\n", i, yTest[i]))
    for (nm in names(models)) {
        pr <- pred[[nm]]$pred
        cat(sprintf("  %-28s [%6.2f, %6.2f]  mean %6.2f\n",
                    nm, pr[i, "lwr"], pr[i, "upr"], pr[i, "fit"]))
    }
}

cat("\n--- E. TOP FEATURES ---\n")
print(summary(fitEns))
