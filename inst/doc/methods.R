## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----setup--------------------------------------------------------------------
library(metafrontier)

## ----det-example--------------------------------------------------------------
sim <- simulate_metafrontier(
  n_groups = 2, n_per_group = 300,
  tech_gap = c(0, 0.4),
  sigma_u = c(0.2, 0.35),
  seed = 123
)

fit_det <- metafrontier(
  log_y ~ log_x1 + log_x2,
  data = sim$data,
  group = "group",
  meta_type = "deterministic"
)

# Metafrontier coefficients (no standard errors)
coef(fit_det, which = "meta")

# Group coefficients for comparison
coef(fit_det, which = "group")

## ----verify-envelop-----------------------------------------------------------
meta_b0 <- coef(fit_det, which = "meta")[1]
group_b0 <- sapply(coef(fit_det, which = "group"), `[`, 1)
meta_b0 >= group_b0

## ----sto-example--------------------------------------------------------------
fit_sto <- metafrontier(
  log_y ~ log_x1 + log_x2,
  data = sim$data,
  group = "group",
  meta_type = "stochastic"
)

summary(fit_sto)

## ----sto-inference------------------------------------------------------------
# Variance-covariance matrix
vcov(fit_sto)

# Log-likelihood of the metafrontier model
logLik(fit_sto)

## ----tgr-range----------------------------------------------------------------
tgr_vals <- efficiencies(fit_sto, type = "tgr")
summary(tgr_vals)

## ----dea-example--------------------------------------------------------------
# CRS metafrontier
fit_crs <- metafrontier(
  log_y ~ log_x1 + log_x2,
  data = sim$data,
  group = "group",
  method = "dea",
  rts = "crs"
)

# VRS metafrontier
fit_vrs <- metafrontier(
  log_y ~ log_x1 + log_x2,
  data = sim$data,
  group = "group",
  method = "dea",
  rts = "vrs"
)

# Compare mean TGR
cbind(
  CRS = tapply(fit_crs$tgr, fit_crs$group_vec, mean),
  VRS = tapply(fit_vrs$tgr, fit_vrs$group_vec, mean)
)

## ----compare-methods----------------------------------------------------------
# Compare TGR estimates across methods
tgr_det <- tapply(fit_det$tgr, fit_det$group_vec, mean)
tgr_sto <- tapply(fit_sto$tgr, fit_sto$group_vec, mean)
tgr_dea <- tapply(fit_crs$tgr, fit_crs$group_vec, mean)
true_tgr <- tapply(sim$data$true_tgr, sim$data$group, mean)

comparison <- data.frame(
  True = true_tgr,
  Deterministic = tgr_det,
  Stochastic = tgr_sto,
  DEA_CRS = tgr_dea
)
round(comparison, 4)

## ----poolability--------------------------------------------------------------
poolability_test(fit_det)

## ----monte-carlo, eval=FALSE--------------------------------------------------
# # Monte Carlo: check parameter recovery over 100 replications
# set.seed(1)
# n_rep <- 100
# beta_hat <- matrix(NA, n_rep, 3)
# 
# for (r in seq_len(n_rep)) {
#   sim_r <- simulate_metafrontier(
#     n_groups = 2, n_per_group = 200,
#     tech_gap = c(0, 0.3),
#     sigma_u = c(0.2, 0.3),
#     sigma_v = 0.15
#   )
#   fit_r <- metafrontier(
#     log_y ~ log_x1 + log_x2,
#     data = sim_r$data,
#     group = "group",
#     meta_type = "deterministic"
#   )
#   beta_hat[r, ] <- coef(fit_r, which = "meta")
# }
# 
# # Bias
# true_beta <- c(1.0, 0.5, 0.3)
# colMeans(beta_hat) - true_beta

