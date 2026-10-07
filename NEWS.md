# bayesrep 1.0.0

First release.

## Model

* `bre()` fits an ensemble of `G` sparse sub-models for Gaussian or binomial
  responses. Every sub-model is fitted to the response in full under a
  composite likelihood, and an Ising prior on the feature-allocation matrix
  penalises re-use of a predictor across sub-models, so correlated groups are
  split across the ensemble.
* A Beta hyperprior on the feature-specific inclusion probabilities handles
  sparsity hierarchically, so it is not a tuning parameter.
* `learningRate` exposes the composite-likelihood learning rate. The default of
  1 fits every sub-model to the data in full; lower values temper each
  sub-model's evidence, at a cost in selection power that can be severe when
  `p` is large.
* Both families carry an explicit intercept.

## Inference

* Collapsed single-site Gibbs updates with a running working residual: `O(np)`
  per sub-model sweep, with no matrix inversions anywhere.
* The feature-specific inclusion probabilities are updated exactly. The
  conjugate Beta draw is correct only when `lambda2 = 0`, since the Ising
  normalising constant depends on the parameter being drawn; the correction is
  material rather than cosmetic.
* A Metropolis move relocates features between sub-models in a single step. A
  single-site sweep cannot do this without passing through a duplicated state
  that the repulsion penalises, so without the move the assignment freezes.
* Pólya-Gamma augmentation for binary responses, including at the non-integer
  shape required when `learningRate` is not 1.
* Posterior summaries accumulate online and the allocation chain is stored one
  byte per entry, with `burnin`, `thin` and a memory budget, so the sampler is
  usable at the dimensions it targets.

## Interface

* `cv.bre()` chooses `lambda2`, `numModels`, `learningRate` and `tauSq` by
  held-out predictive error of the ensemble, with a one-standard-error rule.
* `selectFeatures()` selects at a controlled false discovery rate.
* `coallocation()` reports the posterior probability that two features share a
  sub-model, invariant to relabelling of the sub-models.
* `allocation()` reports label-aligned per-sub-model inclusion probabilities.
* `summary()`, `coef()`, `predict()` and `plot()` methods. Coefficients are
  returned on the original data scale; predictions are averaged over the
  posterior rather than evaluated at a plug-in coefficient vector, and offer
  credible and predictive intervals.
* `waic()` extracts the WAIC accumulated during sampling.

## Notes for users

* `lambda2` is a genuine tuning parameter. It cannot be learned from the
  posterior and must not be selected by WAIC: under the composite likelihood
  every sub-model prefers to keep every useful predictor, so both criteria are
  driven towards no repulsion. Use `cv.bre()`.
* The ensemble marginal inclusion probability is a ranking score, not a
  selection rule. Because a null feature has `G` chances to be carried by some
  sub-model, the meaning of a fixed cutoff changes with `numModels`. Use
  `selectFeatures()`.
* The default `thetaShape2 = p` encodes an expectation of roughly one active
  predictor. It suits sparse problems and should be lowered when many
  predictors are expected to be active, or the fit will badly under-select.
