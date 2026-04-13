# metafrontier 0.2.1

## New features

### Panel SFA and time-varying inefficiency

- `simulate_panel_metafrontier()`: Simulate balanced panel data with
  time-varying inefficiency and technical change for Monte Carlo studies.
- Panel SFA estimation with BC92 (Battese & Coelli, 1992) and BC95
  (Battese & Coelli, 1995) specifications, supporting time-varying
  inefficiency via the `eta` parameter.

### Bootstrap inference

- `boot_tgr()`: Bootstrap confidence intervals for technology gap ratios.
  Supports parametric (residual resampling) and nonparametric (case
  resampling) approaches, with percentile and BCa interval types.
- Parallel bootstrap via `ncores` argument using the `parallel` package.
- S3 methods: `print`, `confint`, and `plot` for `boot_tgr` objects.

### Latent class metafrontier

- `latent_class_metafrontier()`: EM algorithm for estimating latent class
  stochastic frontier models within the metafrontier framework.
- `select_n_classes()`: Automatic selection of the optimal number of
  latent classes via BIC.
- S3 methods: `print`, `summary`, `coef`, and `efficiencies` for
  `lc_metafrontier` objects.

### ggplot2 visualisation

- `autoplot.metafrontier()`: Visualise TGR distributions and efficiency
  decompositions using ggplot2.
- `autoplot.malmquist_meta()`: Plot Malmquist index components over time.
- `autoplot.boot_tgr()`: Visualise bootstrap distributions and
  confidence intervals for TGR.

### Model interoperability

- `as_metafrontier_model()`: Convert pre-fitted model objects from
  `sfaR` (`sfacross`), `frontier` (`sfa`), and `Benchmarking` (`Farrell`)
  packages into metafrontier-compatible format.

### Directional distance functions

- DEA-based metafrontier now supports directional distance functions
  alongside radial efficiency measurement.

## Improvements

- **Murphy-Topel variance correction**: Stochastic metafrontier standard
  errors now account for the generated-regressor problem using the
  Murphy & Topel (1985) correction with PSD enforcement.
- **Improved numerical stability**: Safe Mills ratio computation
  (`.safe_mills()`) prevents NaN in extreme tail regions for both
  panel SFA and latent class models.
- **Memory-efficient latent class**: EM algorithm uses `lm.wfit()` for
  weighted regression instead of explicit diagonal weight matrices.
- **DEA batch solver**: LP object reuse in `dea_batch_fast()` for
  improved performance on large datasets.
- **Input validation**: `boot_tgr()` now validates `R >= 1` with an
  informative error message.
- All `summary()` methods now return invisible S3 objects for
  programmatic access, with `print()` methods for display.
- Expanded `@examples` on all exported S3 methods and autoplot functions.
- Three vignettes added: introduction, Malmquist, and methods.

## Documentation

- Comprehensive vignette: "Introduction to metafrontier" covering
  estimation, inference, bootstrap CIs, panel SFA, latent class,
  and directional distance functions.
- Vignette: "Metafrontier Malmquist Productivity Index" with worked
  panel data examples and identity verification.
- Vignette: "Metafrontier Methods: Theory and Computation" with
  mathematical details and Monte Carlo comparisons.
- Expanded `@param` documentation for `control` argument in
  `metafrontier()` with `optim` options.

---

# metafrontier 0.1.0

Initial CRAN release.

## Estimation

- `metafrontier()`: Main estimation function supporting SFA- and DEA-based
  metafrontiers.
- Deterministic metafrontier (Battese, Rao, and O'Donnell, 2004) via
  constrained optimisation.
- Stochastic metafrontier (Huang, Huang, and Liu, 2014) via second-stage SFA.
- DEA-based metafrontier with CRS, VRS, DRS, and IRS technology assumptions.
- Half-normal, truncated-normal, and exponential inefficiency distributions.
- Heteroscedastic SFA via formula interface (`y ~ x1 + x2 | z1 + z2`).
- Pre-fitted model support via `models=` argument (sfaR, frontier).

## Productivity

- `malmquist_meta()`: Metafrontier Malmquist TFP index with three-way
  decomposition (TEC x TGC x TC*) following O'Donnell, Rao, and Battese (2008).

## Analysis tools

- `technology_gap_ratio()` / `tgr_summary()`: Extract and summarise TGR by group.
- `efficiencies()`: Extract group, metafrontier, and TGR efficiency scores.
- `poolability_test()`: Likelihood ratio test for technology heterogeneity.
- `simulate_metafrontier()`: Data-generating process for Monte Carlo studies.

## Methods

- Full S3 method suite: `print`, `summary`, `coef`, `vcov`, `logLik`,
  `fitted`, `residuals`, `nobs`, `plot`, `confint`, `predict`.
- Four built-in plot types: TGR distributions, efficiency scatter, decomposition,
  and frontier comparison.
