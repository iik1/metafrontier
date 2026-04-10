# metafrontier 0.1.0.9000 (development)

## Initial features

- `metafrontier()`: Main estimation function for SFA- and DEA-based metafrontiers
- Deterministic metafrontier (Battese, Rao, and O'Donnell, 2004) via LP
- Stochastic metafrontier (Huang, Huang, and Liu, 2014) via second-stage SFA
- DEA-based metafrontier with CRS/VRS/DRS/IRS
- `technology_gap_ratio()`: Extract and summarise TGR by group
- `efficiencies()`: Extract group, metafrontier, and TGR efficiency scores
- `poolability_test()`: LR test for common vs group-specific frontiers
- `simulate_metafrontier()`: Data-generating process for Monte Carlo studies
- Full S3 method suite: print, summary, coef, vcov, logLik, fitted, residuals, nobs, plot
- Half-normal, truncated-normal, and exponential inefficiency distributions
