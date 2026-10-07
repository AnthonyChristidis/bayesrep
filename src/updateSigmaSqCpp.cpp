// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// Full conditional draw of the error variance under the composite likelihood
//
//   p(y | B, sigma^2) = prod_g N(y | a0 + X beta^(g), sigma^2)^eta
//                     = (2 pi sigma^2)^{-n G eta / 2}
//                       exp( -(eta / (2 sigma^2)) sum_g || y - a0 - X beta^(g) ||^2 ).
//
// The learning rate eta sets how much evidence each sub-model is allowed to
// draw from the data. At eta = 1/G the normalising constant is that of a single
// Gaussian with n observations, so the ensemble posterior carries the
// information of one pass through the data; at eta = 1 each sub-model uses the
// data in full.
//
// Note that sigma^2 absorbs the mean within-sub-model lack of fit, not only the
// observation noise, and is larger than the error variance of a single
// regression. This is intrinsic to the composite formulation.
//
// [[Rcpp::export]]
double updateSigmaSqCpp(const arma::vec& y,
                        const arma::mat& X,
                        const arma::mat& betaMatrix,
                        double intercept = 0.0,
                        double learningRate = 1.0,
                        double nu0 = 1.0,
                        double sigma0Sq = 1.0) {

    const double n = static_cast<double>(y.n_elem);
    const int numModels = betaMatrix.n_cols;

    const arma::vec centred = y - intercept;

    double ssrTotal = 0.0;
    for (int g = 0; g < numModels; ++g) {
        const arma::vec resid = centred - X * betaMatrix.col(g);
        ssrTotal += arma::dot(resid, resid);
    }

    const double effectiveN = n * static_cast<double>(numModels) * learningRate;

    const double shapePost = (nu0 + effectiveN) / 2.0;
    const double scalePost = (nu0 * sigma0Sq + learningRate * ssrTotal) / 2.0;

    // R::rgamma is parameterised by scale, so pass the reciprocal of the rate.
    return 1.0 / R::rgamma(shapePost, 1.0 / scalePost);
}
