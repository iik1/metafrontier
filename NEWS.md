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
