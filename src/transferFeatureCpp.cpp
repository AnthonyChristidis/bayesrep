// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>

// Metropolis-Hastings move that relocates a feature from one sub-model to
// another in a single step.
//
// A single-site Gibbs sweep cannot move feature j from sub-model g to sub-model
// h directly: it must pass through the state in which both carry j, which the
// Ising prior penalises by exp(-lambda2). The barrier grows with the repulsion,
// so the sweep alone leaves the assignment stuck wherever burn-in reached.
//
// This move proposes the relocation as one joint step. Restricting it to
// features carried by exactly one sub-model keeps the proposal symmetric (the
// relocated feature is still carried exactly once, so the reverse move is drawn
// from the same uniform choice) and leaves m_j and the overlap count unchanged,
// so the Ising and Beta factors cancel. What remains is the ratio of the
// collapsed Bayes factors for carrying j in the destination versus the origin,
//
//   accept with probability  min(1, BF_h / BF_g),
//
// after which beta_j is redrawn from its exact conditional in the destination.
//
// Under the composite likelihood every sub-model is fitted to the response in
// full, so BF_g and BF_h genuinely differ: a feature migrates towards whichever
// sub-model can use it best given what that sub-model already carries. This is
// the mechanism that identifies the partition.
//
// Both families are handled through a common weighted representation. Sub-model
// g is fitted to response column resp.col(g) with observation weights
// sqrtW.col(g) squared. For the Gaussian family sqrtW is all ones and every
// column of resp is y - intercept; for the binomial family the columns carry
// that sub-model's own Polya-Gamma weights and pseudo-response.
//
// [[Rcpp::export]]
Rcpp::List transferFeatureCpp(const arma::mat& resp,
                              const arma::mat& sqrtW,
                              const arma::mat& X,
                              arma::mat betaMatrix,
                              arma::mat gammaMatrix,
                              double sigmaSqScaled,
                              double tauSq,
                              int numMoves) {

    const int p = gammaMatrix.n_rows;
    const int numModels = gammaMatrix.n_cols;

    int numAccepted = 0;
    int numAttempted = 0;

    Rcpp::List out = Rcpp::List::create(
        Rcpp::Named("betaMatrix")   = betaMatrix,
        Rcpp::Named("gammaMatrix")  = gammaMatrix,
        Rcpp::Named("numAccepted")  = 0,
        Rcpp::Named("numAttempted") = 0);

    if (numModels < 2 || numMoves < 1) return out;

    // Candidate features: those carried by exactly one sub-model. An accepted
    // move leaves m_j at one, so this set is fixed for the whole call.
    std::vector<int> candidates;
    candidates.reserve(p);
    for (int j = 0; j < p; ++j) {
        double m = 0.0;
        for (int g = 0; g < numModels; ++g) m += gammaMatrix(j, g);
        if (m == 1.0) candidates.push_back(j);
    }
    if (candidates.empty()) return out;

    // Unweighted linear predictors, maintained incrementally so that each
    // attempt costs O(n) rather than O(np).
    arma::mat fitted = X * betaMatrix;

    const int nCand = static_cast<int>(candidates.size());

    for (int move = 0; move < numMoves; ++move) {

        int pick = static_cast<int>(R::unif_rand() * nCand);
        if (pick >= nCand) pick = nCand - 1;
        const int j = candidates[pick];

        int g = -1;
        for (int k = 0; k < numModels; ++k) {
            if (gammaMatrix(j, k) == 1.0) { g = k; break; }
        }
        if (g < 0) continue;

        // Destination: uniform among the other G - 1 sub-models.
        int offset = static_cast<int>(R::unif_rand() * (numModels - 1));
        if (offset >= numModels - 1) offset = numModels - 2;
        const int h = offset < g ? offset : offset + 1;

        ++numAttempted;

        // Origin residual with feature j removed, and destination residual,
        // both in their own weighted spaces.
        const arma::vec xjG = sqrtW.col(g) % X.col(j);
        const arma::vec xjH = sqrtW.col(h) % X.col(j);

        const arma::vec rG = sqrtW.col(g) %
            (resp.col(g) - fitted.col(g) + X.col(j) * betaMatrix(j, g));
        const arma::vec rH = sqrtW.col(h) % (resp.col(h) - fitted.col(h));

        const double precG = (arma::dot(xjG, xjG) / sigmaSqScaled) + (1.0 / tauSq);
        const double precH = (arma::dot(xjH, xjH) / sigmaSqScaled) + (1.0 / tauSq);
        const double varG = 1.0 / precG;
        const double varH = 1.0 / precH;

        const double meanG = varG * (arma::dot(xjG, rG) / sigmaSqScaled);
        const double meanH = varH * (arma::dot(xjH, rH) / sigmaSqScaled);

        const double logBfG = 0.5 * std::log(varG) - 0.5 * std::log(tauSq)
                            + 0.5 * meanG * meanG * precG;
        const double logBfH = 0.5 * std::log(varH) - 0.5 * std::log(tauSq)
                            + 0.5 * meanH * meanH * precH;

        if (std::log(R::unif_rand()) < (logBfH - logBfG)) {

            fitted.col(g) -= X.col(j) * betaMatrix(j, g);
            betaMatrix(j, g)  = 0.0;
            gammaMatrix(j, g) = 0.0;

            const double betaNew = R::rnorm(meanH, std::sqrt(varH));
            betaMatrix(j, h)  = betaNew;
            gammaMatrix(j, h) = 1.0;
            fitted.col(h) += X.col(j) * betaNew;

            ++numAccepted;
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("betaMatrix")   = betaMatrix,
        Rcpp::Named("gammaMatrix")  = gammaMatrix,
        Rcpp::Named("numAccepted")  = numAccepted,
        Rcpp::Named("numAttempted") = numAttempted
    );
}
