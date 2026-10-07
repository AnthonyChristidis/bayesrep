// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// Posterior co-allocation matrix for a selected set of features.
//
//   C[j, k] = P( features j and k are carried by a common sub-model )
//
// This summary is invariant to relabelling of the G sub-models, so unlike a
// per-sub-model PIP matrix it is well defined whether or not the chain mixes
// over the label-permutation symmetry of the posterior. The diagonal is the
// ensemble marginal inclusion probability, P(feature j active in some model).
//
// The allocation chain is passed as the raw (one byte per entry) buffer used
// for storage, with layout index = j + p * (g + G * t).
//
// [[Rcpp::export]]
arma::mat coallocationCpp(const Rcpp::RawVector& gammaRaw,
                          int p,
                          int numModels,
                          int numDraws,
                          const Rcpp::IntegerVector& featureIdx) {

    const int pSel = featureIdx.size();
    arma::mat counts(pSel, pSel, arma::fill::zeros);

    std::vector<unsigned char> slice(pSel * numModels);

    for (int t = 0; t < numDraws; ++t) {

        const R_xlen_t drawOffset = static_cast<R_xlen_t>(p) * numModels * t;

        for (int a = 0; a < pSel; ++a) {
            const int j = featureIdx[a];
            for (int g = 0; g < numModels; ++g) {
                slice[a * numModels + g] =
                    gammaRaw[drawOffset + j + static_cast<R_xlen_t>(p) * g];
            }
        }

        for (int a = 0; a < pSel; ++a) {
            for (int b = a; b < pSel; ++b) {
                bool shared = false;
                for (int g = 0; g < numModels; ++g) {
                    if (slice[a * numModels + g] && slice[b * numModels + g]) {
                        shared = true;
                        break;
                    }
                }
                if (shared) {
                    counts(a, b) += 1.0;
                    if (a != b) counts(b, a) += 1.0;
                }
            }
        }
    }

    return counts / static_cast<double>(numDraws);
}
