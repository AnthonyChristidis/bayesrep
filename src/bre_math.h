#ifndef BRE_MATH_H
#define BRE_MATH_H

#include <RcppArmadillo.h>
#include <vector>
#include <limits>

// ---------------------------------------------------------------------------
// Shared numerical helpers for the Bayesian Repulsive Ensemble sampler.
//
// The prior on a single row gamma_j = (gamma_j^(1), ..., gamma_j^(G)) of the
// feature-allocation matrix is the Ising / autologistic distribution
//
//   p(gamma_j | theta_j, lambda2)
//       = theta_j^{m} (1 - theta_j)^{G - m} exp(-lambda2 * m (m - 1) / 2)
//         / Z(theta_j, lambda2),        m = sum_g gamma_j^(g)
//
// whose normalising constant is a sum over the G + 1 possible values of m,
//
//   Z(theta, lambda2)
//       = sum_{m=0}^{G} choose(G, m) theta^m (1-theta)^{G-m} exp(-lambda2 m (m-1) / 2).
//
// Z cancels from the full conditional of gamma_j^(g) but NOT from the full
// conditionals of theta_j or lambda2, both of which therefore need it.
// ---------------------------------------------------------------------------

inline double clampProb(double t) {
    if (t < 1e-10) return 1e-10;
    if (t > 1.0 - 1e-10) return 1.0 - 1e-10;
    return t;
}

// log Z(theta, lambda2), evaluated with a log-sum-exp for numerical stability.
inline double logIsingNormConst(double theta, double lambda2, int numModels) {

    theta = clampProb(theta);

    const double logTheta    = std::log(theta);
    const double log1mTheta  = std::log1p(-theta);

    std::vector<double> logTerms(numModels + 1);
    double maxLogTerm = -std::numeric_limits<double>::infinity();

    for (int m = 0; m <= numModels; ++m) {
        const double dm = static_cast<double>(m);
        const double lt = R::lchoose(static_cast<double>(numModels), dm)
                        + dm * logTheta
                        + (static_cast<double>(numModels) - dm) * log1mTheta
                        - lambda2 * dm * (dm - 1.0) / 2.0;
        logTerms[m] = lt;
        if (lt > maxLogTerm) maxLogTerm = lt;
    }

    double acc = 0.0;
    for (int m = 0; m <= numModels; ++m) acc += std::exp(logTerms[m] - maxLogTerm);

    return maxLogTerm + std::log(acc);
}

// Conditional prior log-odds of gamma_j^(g) = 1 given the other sub-models.
// Derived from the Ising joint above; the normalising constant cancels.
inline double calcRepulsiveLogOdds(const arma::mat& gammaMatrix,
                                   int targetGroupZeroIndexed,
                                   int targetFeatureZeroIndexed,
                                   double theta_j,
                                   double lambda2) {

    theta_j = clampProb(theta_j);

    const double baseLogOdds = std::log(theta_j) - std::log1p(-theta_j);

    const int numModels = gammaMatrix.n_cols;
    double overlapCount = 0.0;

    for (int g = 0; g < numModels; ++g) {
        if (g != targetGroupZeroIndexed) {
            overlapCount += gammaMatrix(targetFeatureZeroIndexed, g);
        }
    }

    return baseLogOdds - (lambda2 * overlapCount);
}

#endif
