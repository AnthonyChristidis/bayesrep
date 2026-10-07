# bayesrep

**Bayesian Repulsive Ensembles for high-dimensional regression and classification.**

When predictors are correlated, a sparse model picks one from each correlated
group and discards the rest. The choice is close to arbitrary — resample the
data and a different member often wins — which leaves coefficients unstable and
uncertainty understated.

`bayesrep` fits an ensemble of `G` sparse sub-models at once. Two forces act
together:

- every sub-model has to predict the response **on its own**, so each needs a
  competent set of predictors;
- an **Ising repulsive prior** penalises a sub-model for reaching for a
  predictor another sub-model already holds.

Correlated groups therefore get split *across* sub-models, each reconstructing
the signal through its own member. You get full posterior uncertainty, and an
explicit account of which predictors are acting as substitutes for one another.

## Installation

```r
# install.packages("remotes")
remotes::install_github("AnthonyChristidis/bayesrep")
```

Requires a C++ compiler (Rtools on Windows, Xcode command line tools on macOS).

## Quick start

```r
library(bayesrep)

set.seed(1)
n <- 150; p <- 20
Sigma <- diag(p); Sigma[1:5, 1:5] <- 0.9; diag(Sigma) <- 1
X <- matrix(rnorm(n * p), n, p) %*% chol(Sigma)
colnames(X) <- paste0("V", 1:p)

# V1-V5 are a correlated block, all genuinely active; V6 is independent
beta <- rep(0, p); beta[1:5] <- 1.2; beta[6] <- -2.5
y <- as.numeric(10 + X %*% beta + rnorm(n))

fit <- bre(X, y, numModels = 5, lambda2 = 20, tauSq = 2,
           iter = 3000, burnin = 1500, seed = 1)
fit
#> Bayesian Repulsive Ensemble (BRE)
#> ---------------------------------
#> Family        : gaussian
#> Sub-models (G): 5
#> Data          : n = 150, p = 20
#> MCMC          : 3000 iterations, 1500 burn-in, thin 1, 1500 retained
#> Repulsion     : fixed, lambda2 = 20.000
#> Slab variance : tauSq = 2 (sub-model coefficient scale)
#> theta MH acc. : 96.0%
#> Transfer acc. : 14.2% of 10 proposals per iteration
#> Sub-model size: mean 3.1 active features (range 3 to 6)
#> WAIC          : 2836.1 (pWAIC 32.7)
#> Selected (FDR 10%): 6 features
```

Selection at a controlled false discovery rate:

```r
selectFeatures(fit, fdr = 0.1)
#>   Feature eMIP   PostMean MeanNModels
#> 1      V1    1  1.8306801           3
#> 2      V2    1  1.2585074           2
#> 3      V3    1  1.1481101           2
#> 4      V4    1  1.1235864           2
#> 5      V5    1  0.5235829           1
#> 6      V6    1 -2.6080213           5
```

All six true signals, no false positives. The correlated block sums to 5.88
against a truth of 6.0 (only the block total is identified), and `V6` comes back
at −2.61 against −2.5.

Read `MeanNModels` as how widely a predictor is shared: `V6` carries information
nothing else does, so **all five** sub-models take it. The block members sit in
one to three apiece, because they substitute for each other.

## What you get

| Function | Purpose |
|---|---|
| `bre()` | Fit the ensemble. Gaussian or binomial. |
| `cv.bre()` | Choose `lambda2`, `numModels`, `learningRate`, `tauSq` by held-out error. |
| `summary()` | Features ranked by ensemble marginal inclusion probability, with intervals on the original data scale. |
| `selectFeatures()` | Select at a controlled false discovery rate. |
| `coallocation()` | Posterior probability that two features share a sub-model — which predictors are substitutes. |
| `predict()` | Posterior-averaged predictions with credible or predictive intervals. |
| `plot()` | Importance, co-allocation, allocation, coefficient intervals, MCMC traces. |

## Two things to know before you use it

**`lambda2` must be tuned.** It competes on the log-odds scale against evidence
that grows with sample size and signal strength, so the value that separates a
correlated group on one dataset can be an order of magnitude off on another. It
cannot be learned from the posterior and must not be selected by WAIC — both
drive it to zero, because the objective scores each sub-model on its own fit and
repulsion necessarily makes that worse. Use `cv.bre()`, and widen the grid until
the chosen value is interior.

**Do not threshold `eMIP` at 0.5.** It is the probability that a feature is
carried by *at least one* of `G` sub-models, so a null predictor gets `G` chances
and the meaning of a fixed cutoff shifts with `G`. In a simulation with 350 null
predictors and `G = 5`, `eMIP > 0.5` admitted 347 of them, where
`selectFeatures()` admitted 3. Use eMIP to rank; use `selectFeatures()` to
select.

## When the ensemble helps

The benefit is regime-dependent, and it is worth knowing which side you are on.

- **Many active predictors relative to `n`.** A single sparse model is
  capacity-bound; the sub-models can cover the active set collectively, but only
  if the repulsion stops them covering the same parts. Here repulsion is the
  mechanism rather than a refinement — without it the ensemble is *worse* than a
  single model.
- **Few active predictors.** A single spike-and-slab has ample capacity and
  already averages over which member of a correlated group to include. The
  ensemble adds little, and strong repulsion costs predictive accuracy.

`sandbox/Simulation - Dense AR1 Benchmark.R` benchmarks both regimes against a
single spike-and-slab, the no-repulsion ablation and the lasso.

## Documentation

```r
vignette("bayesrep")   # worked example, start here
?bre                   # the model and every argument
?cv.bre                # choosing the tuning parameters
?coallocation          # reading the allocation structure
```

## Status

Research software under active development. The interface may change. Not yet
on CRAN.

## License

GPL-3.
