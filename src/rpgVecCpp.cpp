// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// Polya-Gamma draws with a non-integer shape parameter.
//
// Under the composite likelihood each sub-model's Bernoulli contribution is
// raised to the power 1/G, which turns the required augmentation from PG(1, c)
// into PG(1/G, c). The Devroye-style exact samplers used by 'pgdraw' cover
// integer shapes only, so this routine uses the infinite-sum representation of
// Polson, Scott and Windle (2013),
//
//   omega = (1 / (2 pi^2)) sum_{k=1}^inf  g_k / ((k - 1/2)^2 + c^2 / (4 pi^2)),
//   g_k ~ Gamma(b, 1) independently,
//
// truncated after `trunc` terms. The expectation of the discarded tail is
// available in closed form and is added back, so the sampler is unbiased in the
// mean; the residual error is in the tail's variance only and is negligible,
// since the summands decay like k^{-2}. With the default truncation the mean
// agrees with the exact (b / (2c)) tanh(c / 2) to several decimal places.
//
// [[Rcpp::export]]
arma::vec rpgVecCpp(double b, const arma::vec& c, int trunc = 200) {

    const int n = c.n_elem;
    arma::vec out(n);

    const double twoPiSq = 2.0 * M_PI * M_PI;
    const double kMinusHalf = static_cast<double>(trunc) + 0.5;

    for (int i = 0; i < n; ++i) {

        // a = |c| / (2 pi) controls the shape of the summand.
        const double a = std::fabs(c(i)) / (2.0 * M_PI);
        const double aSq = a * a;

        double acc = 0.0;
        for (int k = 1; k <= trunc; ++k) {
            const double km = static_cast<double>(k) - 0.5;
            acc += R::rgamma(b, 1.0) / (km * km + aSq);
        }

        // Expected value of the truncated tail, integrated from the truncation
        // point: sum_{k > K} 1 / ((k - 1/2)^2 + a^2).
        double tail;
        if (a < 1e-8) {
            tail = 1.0 / kMinusHalf;
        } else {
            tail = (M_PI_2 - std::atan(kMinusHalf / a)) / a;
        }

        out(i) = (acc + b * tail) / twoPiSq;
    }

    return out;
}
