# metafrontier

Metafrontier Analysis for Efficiency and Productivity Studies

## Overview

`metafrontier` provides a comprehensive R implementation of metafrontier production function models for estimating technical efficiencies and technology gaps for firms operating under different technologies.

The package implements:

- **Deterministic metafrontier** (Battese, Rao, and O'Donnell, 2004) via constrained optimisation
- **Stochastic metafrontier** (Huang, Huang, and Liu, 2014) via second-stage SFA
- **DEA-based metafrontier** using pooled data envelopment analysis
- **Metafrontier Malmquist index** (O'Donnell, Rao, and Battese, 2008) for productivity decomposition
- **Technology gap ratio (TGR)** computation with summary statistics and visualisation
- **Poolability tests** (likelihood ratio) for testing common vs group-specific frontiers

## Installation

Install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("iik1/metafrontier")
```

## Quick Start

```r
library(metafrontier)

# Simulate data with two technology groups
sim <- simulate_metafrontier(n_groups = 2, n_per_group = 200, seed = 42)

# Estimate metafrontier (SFA-based, deterministic)
fit <- metafrontier(log_y ~ log_x1 + log_x2,
                    data = sim$data,
                    group = "group")
summary(fit)

# Technology gap ratios by group
tgr_summary(fit)

# Visualise
plot(fit, which = "decomposition")
```

## Key References

- Battese, G.E., Rao, D.S.P. and O'Donnell, C.J. (2004). A metafrontier production function for estimation of technical efficiencies and technology gaps for firms operating under different technologies. *Journal of Productivity Analysis*, 21(1), 91-103.
- Huang, C.J., Huang, T.-H. and Liu, N.-H. (2014). A new approach to estimating the metafrontier production function based on a stochastic frontier framework. *Journal of Productivity Analysis*, 42(3), 241-254.
- O'Donnell, C.J., Rao, D.S.P. and Battese, G.E. (2008). Metafrontier frameworks for the study of firm-level efficiencies and technology ratios. *Empirical Economics*, 34(2), 231-255.

## License

GPL (>= 3)
