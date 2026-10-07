// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <numeric>

// Label alignment for the sub-model-specific inclusion probabilities.
//
// The BRE posterior is invariant under permutation of the G sub-model labels,
// so raw per-sub-model PIPs averaged across draws are only interpretable if the
// chain fails to mix over that symmetry. This routine removes the ambiguity by
// aligning each draw's columns to the running aligned mean before accumulating,
// which is the standard relabelling fix used for mixture models.
//
// The assignment problem is tiny (G x G) and is solved exactly by enumeration
// when G <= maxExact, and greedily otherwise. The fraction of draws whose
// optimal permutation is not the identity is returned as an honest diagnostic
// of how much label switching the chain actually exhibits.
//
// [[Rcpp::export]]
Rcpp::List alignedPipCpp(const Rcpp::RawVector& gammaRaw,
                         int p,
                         int numModels,
                         int numDraws,
                         int maxExact = 7) {

    arma::mat alignedSum(p, numModels, arma::fill::zeros);

    const bool useExact = (numModels <= maxExact);
    int numSwitched = 0;

    std::vector<int> perm(numModels), bestPerm(numModels);
    std::vector<double> sim(numModels * numModels);

    for (int t = 0; t < numDraws; ++t) {

        const R_xlen_t drawOffset = static_cast<R_xlen_t>(p) * numModels * t;

        if (t == 0) {
            for (int g = 0; g < numModels; ++g) {
                for (int j = 0; j < p; ++j) {
                    alignedSum(j, g) +=
                        static_cast<double>(gammaRaw[drawOffset + j + static_cast<R_xlen_t>(p) * g]);
                }
            }
            continue;
        }

        // Similarity between draw column a and reference column b, where the
        // reference is the running aligned mean. Only active features
        // contribute, so the inner loop is driven by the sparsity of Gamma.
        std::fill(sim.begin(), sim.end(), 0.0);
        const double refScale = 1.0 / static_cast<double>(t);

        for (int a = 0; a < numModels; ++a) {
            for (int j = 0; j < p; ++j) {
                if (gammaRaw[drawOffset + j + static_cast<R_xlen_t>(p) * a]) {
                    for (int b = 0; b < numModels; ++b) {
                        sim[a * numModels + b] += alignedSum(j, b) * refScale;
                    }
                }
            }
        }

        double bestScore = -1.0;

        if (useExact) {
            std::iota(perm.begin(), perm.end(), 0);
            do {
                double score = 0.0;
                for (int a = 0; a < numModels; ++a) score += sim[a * numModels + perm[a]];
                if (score > bestScore) {
                    bestScore = score;
                    bestPerm = perm;
                }
            } while (std::next_permutation(perm.begin(), perm.end()));
        } else {
            std::vector<bool> usedCol(numModels, false), usedRef(numModels, false);
            std::fill(bestPerm.begin(), bestPerm.end(), -1);
            for (int step = 0; step < numModels; ++step) {
                int bestA = -1, bestB = -1;
                double bestVal = -1.0;
                for (int a = 0; a < numModels; ++a) {
                    if (usedCol[a]) continue;
                    for (int b = 0; b < numModels; ++b) {
                        if (usedRef[b]) continue;
                        if (sim[a * numModels + b] > bestVal) {
                            bestVal = sim[a * numModels + b];
                            bestA = a;
                            bestB = b;
                        }
                    }
                }
                usedCol[bestA] = true;
                usedRef[bestB] = true;
                bestPerm[bestA] = bestB;
            }
        }

        bool isIdentity = true;
        for (int a = 0; a < numModels; ++a) {
            if (bestPerm[a] != a) { isIdentity = false; break; }
        }
        if (!isIdentity) ++numSwitched;

        for (int a = 0; a < numModels; ++a) {
            const int b = bestPerm[a];
            for (int j = 0; j < p; ++j) {
                alignedSum(j, b) +=
                    static_cast<double>(gammaRaw[drawOffset + j + static_cast<R_xlen_t>(p) * a]);
            }
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("pip")        = alignedSum / static_cast<double>(numDraws),
        Rcpp::Named("switchRate") = static_cast<double>(numSwitched) /
                                    static_cast<double>(std::max(numDraws - 1, 1))
    );
}
