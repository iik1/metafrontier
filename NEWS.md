# metafrontier (development version)

## Breaking changes

- `boot_tgr(type = "nonparametric")` no longer returns
  observation-level intervals: `ci` is `NULL` and `confint()` signals
  an error. Each replicate refits the model on rows resampled within
  groups, so column `i` of `tgr_boot` is a random draw from a group,
  not DMU `i`, and the previous `ci[i, ]` was not specific to DMU `i`
  whatever the row order. Evaluating the original DMUs against each
  bootstrap frontier instead would, for DEA, be the naive bootstrap,
  which is inconsistent (Kneip, Simar and Wilson, 2008). The
  parametric bootstrap keeps the design fixed and still returns
  observation-level intervals.

## Bug fixes

- `boot_tgr(type = "nonparametric")`: group-level intervals mixed
  bootstrap draws across groups whenever the data were not sorted by
  group. Each replicate stacks the resampled rows group by group, so
  column `j` of `tgr_boot` is a draw from the `j`-th position block,
  but `ci_group`, and the distribution panels of `plot()` and
  `autoplot()`, selected columns by the groups of the original rows.
  With unsorted groups of different TGR levels, the intervals were
  pulled towards the pooled mean and could exclude the group estimate.
  Group intervals now use the position blocks, stored in the new
  element `boot_group`, and do not depend on row order
  (regression-tested). Results for data already sorted by group are
  unchanged. Nonparametric `boot_tgr` objects created with earlier
  versions from unsorted data should be recomputed.
- `boot_tgr()` now bootstraps only the rows used in the fit. For an
  SFA fit with rows dropped for missing values, the nonparametric
  bootstrap previously stopped with a replacement-length error, and the
  parametric bootstrap, with only recycling warnings, simulated each
  group's response from fitted frontier values shifted out of line with
  the inputs from the group's first dropped row onwards, giving wrong
  intervals. Both now match the bootstrap on the complete cases
  (regression-tested).
- `boot_tgr()` refitted each replicate by looking up a column named
  `"group"` in the data whenever the fit was given `group` as a vector.
  Without such a column every replicate failed; with an unrelated
  column of that name the replicates were silently regrouped by it. The
  refits now use the fitted group labels directly, and results match
  those of the same fit made with a column name (regression-tested).
- `boot_tgr()` now stops with an informative error for metafrontiers
  built from pre-fitted `models`, whose group models cannot be
  re-estimated on bootstrap samples. Every replicate previously failed,
  leaving only a warning. Fit from `formula`, `data` and `group`
  instead.
- `boot_tgr(type = "parametric")` failed in every replicate for fits
  with `engine = "sfaR"` or `engine = "frontier"`: their group models
  stored no fitted frontier values, and the error scales `sigma_v` and
  `sigma_u` of `sfaR` fits were `NA` because every `sfaR` release
  reports them as log variances (`Zv_(Intercept)`, `Zu_(Intercept)`),
  not under the names the extractor looked for. External-engine group
  models now store `fitted`, and `sfaR` fits yield `sigma_v`, `sigma_u`
  and, for `dist = "tnormal"`, `mu`, matching the internal estimates
  (regression-tested for `"hnormal"`, `"tnormal"` and `"exponential"`).
  The corrected scales also apply to `sfacross` objects passed to
  `metafrontier(models = )` or `as_metafrontier_model()`.

## New features

- `boot_tgr()` returns percentile intervals for the group median TGR
  (`ci_group_median`) alongside those for the group mean
  (`ci_group`), and `print()` shows both.

# metafrontier 0.3.1

## Bug fixes

- BC88/BC92 efficiencies are now evaluated on the log scale
  (`pnorm(log.p = TRUE)`), with an asymptotic closed form in the far
  left tail of the conditional distribution. The previous
  natural-scale ratio `Phi(mu*/sigma* - sigma*) / Phi(mu*/sigma*)`
  could return 0, `Inf`, or values above one on platforms without
  extended-precision long doubles (the CRAN noLD additional check)
  when both `Phi` terms were subnormal, which occurred in degenerate
  fits with `sigma_u` near zero. On standard platforms results are
  unchanged to near machine precision; in the affected tail region the
  new evaluation is also more accurate than the previous fallback to
  the JLMS value.

# metafrontier 0.3.0

## Breaking changes

- The default technical efficiency estimator for SFA models is now the
  Battese-Coelli (1988) conditional expectation `E[exp(-u)|eps]`
  (`estimator = "bc88"`), which provides consistent efficiency
  estimates. The previous JLMS estimator remains available via
  `estimator = "jlms"`, and both are stored on every fit;
  `efficiencies(fit, estimator = )` switches without refitting.
  Deterministic and stochastic TGRs are unaffected (they do not depend
  on the estimator), but TE and TE* values change slightly relative to
  0.2.x.
- `malmquist_meta()` output `id` column now carries the user-supplied
  firm identifier (previously a within-group loop index).
- The `as_metafrontier_model()` method for `frontier::sfa()` fits is
  now registered for the correct class `"frontier"`; the previous
  registration (`"sfa"`) never dispatched.

## Bug fixes

- Panel SFA (BC92/BC95): fixed a row-alignment bug where
  observation-level efficiencies were returned in internally sorted
  (id, time) order but assigned back in input-row order, scrambling
  `te_group`/`te_meta` whenever the input rows were not already in
  string-sorted order. Coefficients, TGRs, and group means were
  unaffected. Results are now row-order invariant (regression-tested).
- BC92 decay is now anchored at the panel-wide final period `T`
  (previously the firm-specific last period), matching Battese and
  Coelli's (1992) unbalanced formulation and the package's own
  simulator. Balanced panels are numerically unchanged.
- Panel SFA now aligns the firm/time index with rows dropped by
  `na.action`, and uses the same BFGS-to-Nelder-Mead fallback cascade
  as the cross-sectional path.
- DEA: removed an erroneous lower bound (`phi >= 1`) in the
  output-oriented LPs that made cross-period evaluation of
  super-efficient DMUs infeasible even under CRS, silently biasing
  Malmquist TC and MPI means. Cross-period scores with `phi < 1` now
  solve correctly.
- `as_metafrontier_model()` is now idempotent (converting an already
  converted object is a no-op), so the previously documented
  pre-conversion workflow works.
- `poolability_test()` now derives `data.name` from the passed
  expression instead of deparsing the stored call.
- `autoplot(boot)` (and the base `plot()` method) now actually draw
  the dashed CI bound lines promised by the documentation.
- `autoplot(malm, which = "mpi_trend")` now plots each transition at
  its end period (axis 2..T), consistent with the caption convention
  "change relative to the previous period".
- Removed a redundant marginal likelihood computation in the latent
  class EM loop.

## New features

- `estimator = c("bc88", "jlms")` on `metafrontier()` and
  `malmquist_meta()` (see Breaking changes).
- `objective = c("lp", "qp")` on `metafrontier()`: both identification
  criteria of Battese, Rao and O'Donnell (2004) for the deterministic
  metafrontier. The default LP minimises the sum of absolute
  deviations (O'Donnell, Rao and Battese, 2008, Eqs. 23-25); the QP
  minimises the sum of squared deviations, solved exactly via
  `quadprog` (new in Suggests) with a `constrOptim()` barrier
  fallback. The bootstrap respects the choice.
- `engine = c("internal", "sfaR", "frontier", "Benchmarking")` on
  `metafrontier()`: delegate group-frontier estimation to external
  packages; `engine = "Benchmarking"` also delegates the pooled DEA
  metafrontier via `XREF`/`YREF`.
- `check_convergence()`: new exported diagnostic reporting one row per
  estimation stage (group frontiers and metafrontier). `print()` and
  `summary()` methods now include convergence status; every stage
  warns on non-zero optimiser codes.
- `malmquist_meta(id = )`: explicit firm matching across periods, with
  errors on duplicated (id, period) pairs, warnings counting dropped
  observations on unbalanced panels, and a message when falling back
  to positional matching. Cross-period infeasible programmes (possible
  under vrs/drs/irs/fdh) are now counted and reported in a
  consolidated warning, stored as `n_infeasible`, and shown by
  `print()`/`summary()`; the SFA path announces its pointwise-maximum
  approximation in a message and in the documentation.
- Unbalanced panels are fully supported in BC92 (firm-specific period
  sets), and `simulate_panel_metafrontier()` gains an `attrition`
  argument for generating unbalanced test panels.
- DEA: `rts = "fdh"` (free disposable hull; exact enumeration for
  radial measures, binary MIP for DDF), `type = "hyperbolic"` (graph
  efficiency; closed form under CRS, bisection otherwise; always
  feasible cross-period), user-supplied direction vectors or
  firm-specific direction matrices for DDF (reported via the additive
  `ddf_gap`), and `slack = TRUE` two-stage slack maximisation.
- `poolability_test()` now dispatches a permutation test for DEA fits
  (group labels exchangeable under the pooled-technology null), with
  `B` and `seed` arguments.
- `coef()`, `vcov()`, and `summary()` expose all estimated parameters:
  `extraPar = TRUE` returns variance parameters (and `eta` for BC92)
  with back-transformed values; `vcov(which = "group")` returns full
  per-group covariance matrices; group summary tables now include the
  variance parameters and `eta` with standard errors.
- `simulate_metafrontier()` gains `beta_groups` (group-specific slope
  coefficients; the true TGR is then computed against the pointwise
  maximum over group frontiers and varies within groups),
  `input_means` (group-specific input distributions), and
  `input_corr` (correlated log inputs).
- Latent class: the EM now records whether the best start met the
  convergence tolerance (`em_converged`) and warns when it did not.

---

# metafrontier 0.2.2

## Bug fixes

- `boot_tgr()`: Fixed orientation/rts not propagating to bootstrap
  replicates (always defaulted to output/CRS).
- `boot_tgr()`: Fixed hardcoded group column name; now respects
  the user's original group variable.
- Latent class `.loglik_to_u_hat()` now respects the `dist` argument
  with correct JLMS formulas for half-normal, truncated-normal, and
  exponential distributions.
- Latent class parameter count now correct for truncated-normal
  (k+3 instead of k+2).
- Fixed two incorrect DOIs: Huang et al. (2014) and O'Donnell et al.
  (2008) now link to the correct papers.

## Improvements

- `autoplot` methods now use proper conditional S3 registration
  (`@exportS3Method ggplot2::autoplot`) instead of direct `export()`.
- `.extract_benchmarking()` now attempts to retrieve XREF/YREF from
  Farrell objects and `.estimate_from_models()` gives a clear error
  when DEA models lack the required X/y/beta.
- `technology_gap_ratio()` documentation now correctly notes that TGR
  can exceed 1 under the stochastic metafrontier.
- Pre-built vignettes included in inst/doc for reliable R CMD check.

---

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
