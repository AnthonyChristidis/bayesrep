// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "bre_math.h"

// Joint single-site update of (gamma_j^(g), beta_j^(g)) for every feature j
// within one sub-model g.
//
// For each j the coefficient beta_j^(g) is analytically integrated out, the
// indicator gamma_j^(g) is drawn from its collapsed full conditional, and
// beta_j^(g) is then drawn from its exact conditional Normal (or set to zero).
// A running working residual avoids any matrix inversion, so the sweep costs
// O(n p) rather than O(p_active^3).
//
// [[Rcpp::export]]
Rcpp::List updateBetaGammaCpp(const arma::vec& yScaled,
                              const arma::mat& X,
                              const arma::vec& xTx,
                              arma::vec betaGroup,
                              arma::mat gammaMatrix,
                              int targetGroup,
                              double sigmaSqScaled,
                              double tauSqScaled,
                              const arma::vec& thetaVec,
                              double lambda2) {

    const int p = X.n_cols;
    const int g_idx = targetGroup - 1;

    arma::vec currentResidual = yScaled - (X * betaGroup);

    for (int j = 0; j < p; ++j) {

        // Remove feature j's current contribution from the working residual.
        if (betaGroup(j) != 0.0) {
            currentResidual += X.col(j) * betaGroup(j);
        }

        const double precPost = (xTx(j) / sigmaSqScaled) + (1.0 / tauSqScaled);
        const double varPost  = 1.0 / precPost;

        const double xjTrj   = arma::dot(X.col(j), currentResidual);
        const double meanPost = varPost * (xjTrj / sigmaSqScaled);

        // log Bayes factor for gamma_j = 1 vs gamma_j = 0, beta_j collapsed out.
        const double logBf = 0.5 * std::log(varPost)
                           - 0.5 * std::log(tauSqScaled)
                           + 0.5 * (meanPost * meanPost * precPost);

        const double logPriorOdds = calcRepulsiveLogOdds(gammaMatrix, g_idx, j,
                                                         thetaVec(j), lambda2);

        const double probInclusion = R::plogis(logPriorOdds + logBf, 0.0, 1.0, 1, 0);
        const double newGamma = R::rbinom(1, probInclusion);

        gammaMatrix(j, g_idx) = newGamma;

        if (newGamma == 1.0) {
            betaGroup(j) = R::rnorm(meanPost, std::sqrt(varPost));
            currentResidual -= X.col(j) * betaGroup(j);
        } else {
            betaGroup(j) = 0.0;
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("betaGroup")   = betaGroup,
        Rcpp::Named("gammaMatrix") = gammaMatrix
    );
}
