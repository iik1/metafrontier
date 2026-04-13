# metafrontier

Analysis of Metafrontier Models for Efficiency and Productivity

## Overview

`metafrontier` provides a unified R implementation of metafrontier production function models for estimating technical efficiencies and technology gaps across groups of firms operating under different technologies.

### Estimation methods

- **Deterministic metafrontier** (Battese, Rao & O'Donnell, 2004) via constrained LP/QP optimisation
- **Stochastic metafrontier** (Huang, Huang & Liu, 2014) via second-stage SFA with Murphy-Topel corrected standard errors
- **DEA-based metafrontier** with CRS, VRS, DRS, and IRS technology assumptions
- **Latent class metafrontier** via EM algorithm with BIC-based class selection

### Productivity analysis

- **Metafrontier Malmquist TFP index** (O'Donnell, Rao & Battese, 2008) with three-way decomposition (TEC x TGC x TC*)
- **Panel SFA** with time-varying inefficiency (BC92 and BC95 specifications)
- **Directional distance functions** for DEA-based efficiency measurement

### Inference and diagnostics

- **Bootstrap confidence intervals** for TGR (parametric and nonparametric; percentile and BCa)
- **Murphy-Topel variance correction** for stochastic metafrontier standard errors
- **Poolability tests** (likelihood ratio) for common vs group-specific frontiers
- **Half-normal, exponential, and truncated-normal** inefficiency distributions

### Visualisation

- Base R `plot()` with four plot types: TGR distributions, efficiency scatter, frontier decomposition, and frontier comparison
- `ggplot2` integration via `autoplot()` methods for metafrontier, Malmquist, and bootstrap results

### Interoperability

- Import pre-fitted models from `sfaR`, `frontier`, and `Benchmarking` via `as_metafrontier_model()`
- Formula interface with heteroscedastic SFA support (`y ~ x1 + x2 | z1 + z2`)
- Full S3 method suite: `print`, `summary`, `coef`, `vcov`, `logLik`, `fitted`, `residuals`, `nobs`, `confint`, `predict`, `plot`

## Installation

Install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("iik1/metafrontier")
```

## Quick start

```r
library(metafrontier)

# Simulate metafrontier data with two technology groups
sim <- simulate_metafrontier(n_groups = 2, n_per_group = 200, seed = 42)

# Estimate a deterministic SFA-based metafrontier
fit_det <- metafrontier(log_y ~ log_x1 + log_x2,
                        data = sim$data, group = "group")

# Estimate a stochastic metafrontier (with Murphy-Topel SEs)
fit_sto <- metafrontier(log_y ~ log_x1 + log_x2,
                        data = sim$data, group = "group",
                        meta_type = "stochastic")

# DEA-based metafrontier (requires level-scale inputs/outputs)
dat_lev <- within(sim$data, { y <- exp(log_y); x1 <- exp(log_x1); x2 <- exp(log_x2) })
fit_dea <- metafrontier(y ~ x1 + x2,
                        data = dat_lev, group = "group",
                        method = "dea", rts = "vrs")

# Inspect results
summary(fit_det)
tgr_summary(fit_det)
confint(fit_det)
```

### Bootstrap confidence intervals

```r
boot <- boot_tgr(fit_det, R = 999, seed = 1)
confint(boot)

# Parallel bootstrap
boot_par <- boot_tgr(fit_det, R = 999, ncores = 4, seed = 1)
```

### Malmquist productivity index

```r
# Simulate panel data
panel <- simulate_panel_metafrontier(n_groups = 3, n_firms_per_group = 50,
                                     n_periods = 5, seed = 42)

# Three-way Malmquist decomposition
malm <- malmquist_meta(log_y ~ log_x1 + log_x2,
                       data = panel$data, group = "group",
                       time = "year")
summary(malm)
```

### Latent class metafrontier

```r
# Automatic class selection via BIC (discovers groups endogenously)
lc <- latent_class_metafrontier(log_y ~ log_x1 + log_x2,
                                data = sim$data,
                                n_classes = 3)
summary(lc)
```

### Visualisation

```r
# Base R
plot(fit_det, which = "tgr")
plot(fit_det, which = "decomposition")

# ggplot2
library(ggplot2)
autoplot(fit_det)
autoplot(boot)
autoplot(malm)
```

### Using pre-fitted models

```r
library(sfaR)

# Fit group-specific SFA models externally
sfa_g1 <- sfacross(log_y ~ log_x1 + log_x2,
                   data = subset(sim$data, group == "G1"))
sfa_g2 <- sfacross(log_y ~ log_x1 + log_x2,
                   data = subset(sim$data, group == "G2"))

# Pass to metafrontier
fit <- metafrontier(models = list(G1 = sfa_g1, G2 = sfa_g2))
```

## Vignettes

The package includes three vignettes:

- **Introduction to metafrontier** -- end-to-end walkthrough covering estimation, inference, bootstrap CIs, panel SFA, latent class, and directional distance functions
- **Metafrontier Malmquist Productivity Index** -- panel data productivity decomposition with worked examples
- **Metafrontier Methods: Theory and Computation** -- mathematical details, comparison of deterministic/stochastic/DEA approaches, and Monte Carlo evidence

```r
browseVignettes("metafrontier")
```

## References

- Battese, G.E., Rao, D.S.P. and O'Donnell, C.J. (2004). A metafrontier production function for estimation of technical efficiencies and technology gaps for firms operating under different technologies. *Journal of Productivity Analysis*, 21(1), 91--103. [doi:10.1023/B:PROD.0000012454.06094.29](https://doi.org/10.1023/B:PROD.0000012454.06094.29)
- Huang, C.J., Huang, T.-H. and Liu, N.-H. (2014). A new approach to estimating the metafrontier production function based on a stochastic frontier framework. *Journal of Productivity Analysis*, 42(3), 241--254. [doi:10.1007/s11123-014-0402-2](https://doi.org/10.1007/s11123-014-0402-2)
- O'Donnell, C.J., Rao, D.S.P. and Battese, G.E. (2008). Metafrontier frameworks for the study of firm-level efficiencies and technology ratios. *Empirical Economics*, 34(2), 231--255. [doi:10.1007/s00181-007-0119-4](https://doi.org/10.1007/s00181-007-0119-4)

## License

GPL (>= 3)
