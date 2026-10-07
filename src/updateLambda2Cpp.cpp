// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "bre_math.h"

// Metropolis-Hastings update of the repulsion strength lambda2.
//
// Collecting the lambda2-dependent terms of log p(Gamma | theta, lambda2) over
// all p feature rows gives
//
//   log L(lambda2) = -lambda2 * sum_j m_j (m_j - 1) / 2
//                    - sum_j log Z(theta_j, lambda2),
//
// which is O(p G) to evaluate. Combined with a Gamma(shape, rate) prior this
// lets lambda2 be learned from the data instead of fixed by a grid search.
// The proposal is a Gaussian random walk on the log scale; the Jacobian of the
// transform contributes log(lambda2Prop) - log(lambda2Cur) to the ratio.
//
// [[Rcpp::export]]
Rcpp::List updateLambda2Cpp(double lambda2,
                            const arma::vec& thetaVec,
                            const arma::vec& mCounts,
                            int numModels,
                            double priorShape,
                            double priorRate,
                            double proposalSd) {

    const int p = thetaVec.n_elem;

    // Sufficient statistic: total number of overlapping pairs in Gamma.
    double pairCount = 0.0;
    for (int j = 0; j < p; ++j) {
        pairCount += mCounts(j) * (mCounts(j) - 1.0) / 2.0;
    }

    const double logLambdaCur  = std::log(lambda2);
    const double logLambdaProp = logLambdaCur + R::norm_rand() * proposalSd;
    const double lambdaProp    = std::exp(logLambdaProp);

    // log posterior kernel, including the log-scale Jacobian.
    double logPostCur  = (priorShape - 1.0) * logLambdaCur  - priorRate * lambda2
                       - lambda2 * pairCount + logLambdaCur;
    double logPostProp = (priorShape - 1.0) * logLambdaProp - priorRate * lambdaProp
                       - lambdaProp * pairCount + logLambdaProp;

    for (int j = 0; j < p; ++j) {
        logPostCur  -= logIsingNormConst(thetaVec(j), lambda2,    numModels);
        logPostProp -= logIsingNormConst(thetaVec(j), lambdaProp, numModels);
    }

    int accepted = 0;
    if (std::log(R::unif_rand()) < (logPostProp - logPostCur)) {
        lambda2  = lambdaProp;
        accepted = 1;
    }

    return Rcpp::List::create(
        Rcpp::Named("lambda2")  = lambda2,
        Rcpp::Named("accepted") = accepted
    );
}
