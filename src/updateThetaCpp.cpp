// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "bre_math.h"

// Exact full-conditional update of the feature-specific baseline inclusion
// probabilities theta_j.
//
// The naive conjugate draw Beta(a + m_j, b + G - m_j) is the correct full
// conditional only when lambda2 = 0, because the Ising normalising constant
// Z(theta_j, lambda2) depends on theta_j:
//
//   p(theta_j | gamma_j) propto theta^{a+m-1} (1-theta)^{b+G-m-1} / Z(theta, lambda2).
//
// We therefore use the Beta draw as an independence Metropolis-Hastings
// proposal. Every factor except Z cancels from the acceptance ratio, leaving
//
//   alpha = min(1, Z(theta_current, lambda2) / Z(theta_proposed, lambda2)),
//
// which costs a sum over G + 1 terms and restores exactness of the sampler.
//
// [[Rcpp::export]]
Rcpp::List updateThetaCpp(arma::vec thetaVec,
                          const arma::vec& mCounts,
                          double aPrior,
                          double bPrior,
                          double lambda2,
                          int numModels,
                          int numSweeps = 1) {

    const int p = thetaVec.n_elem;
    int numAccepted = 0;
    int numProposed = 0;

    // With no repulsion the Ising prior collapses to independent Bernoullis and
    // the conjugate Beta draw is exact, so one sweep suffices.
    const bool conjugate = (lambda2 <= 0.0);
    if (conjugate) numSweeps = 1;

    for (int sweep = 0; sweep < numSweeps; ++sweep) {
        for (int j = 0; j < p; ++j) {

            const double shape1 = aPrior + mCounts(j);
            const double shape2 = bPrior + static_cast<double>(numModels) - mCounts(j);

            const double thetaProp = clampProb(R::rbeta(shape1, shape2));
            ++numProposed;

            if (conjugate) {
                thetaVec(j) = thetaProp;
                ++numAccepted;
                continue;
            }

            const double logZCur  = logIsingNormConst(thetaVec(j), lambda2, numModels);
            const double logZProp = logIsingNormConst(thetaProp,   lambda2, numModels);

            if (std::log(R::unif_rand()) < (logZCur - logZProp)) {
                thetaVec(j) = thetaProp;
                ++numAccepted;
            }
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("thetaVec")    = thetaVec,
        Rcpp::Named("numAccepted") = numAccepted,
        Rcpp::Named("numProposed") = numProposed
    );
}
