## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----setup--------------------------------------------------------------------
library(metafrontier)

## ----simulate-----------------------------------------------------------------
sim <- simulate_metafrontier(
  n_groups = 3,
  n_per_group = 200,
  beta_meta = c(1.0, 0.5, 0.3),  # intercept, elasticity_1, elasticity_2
  tech_gap = c(0, 0.25, 0.5),    # intercept shifts (0 = best technology)
  sigma_u = c(0.2, 0.3, 0.4),    # inefficiency SD per group
  sigma_v = 0.15,                 # noise SD
  seed = 42
)

str(sim$data[, c("log_y", "log_x1", "log_x2", "group")])
table(sim$data$group)

## ----estimate-----------------------------------------------------------------
fit <- metafrontier(
  log_y ~ log_x1 + log_x2,
  data = sim$data,
  group = "group",
  method = "sfa",
  meta_type = "deterministic"
)

fit

## ----deterministic------------------------------------------------------------
fit_det <- metafrontier(
  log_y ~ log_x1 + log_x2,
  data = sim$data,
  group = "group",
  meta_type = "deterministic"
)

summary(fit_det)

## ----stochastic---------------------------------------------------------------
fit_sto <- metafrontier(
  log_y ~ log_x1 + log_x2,
  data = sim$data,
  group = "group",
  meta_type = "stochastic"
)

summary(fit_sto)

## ----vcov---------------------------------------------------------------------
vcov(fit_sto)

## ----dea----------------------------------------------------------------------
fit_dea <- metafrontier(
  log_y ~ log_x1 + log_x2,
  data = sim$data,
  group = "group",
  method = "dea",
  rts = "vrs"
)

fit_dea

## ----efficiencies-------------------------------------------------------------
te <- efficiencies(fit_det, type = "group")
tgr <- efficiencies(fit_det, type = "tgr")
te_star <- efficiencies(fit_det, type = "meta")

# Verify the fundamental identity: TE* = TE x TGR
all.equal(te_star, te * tgr)

## ----tgr----------------------------------------------------------------------
tgr_by_group <- technology_gap_ratio(fit_det)
lapply(tgr_by_group, summary)

## ----tgr-summary--------------------------------------------------------------
tgr_summary(fit_det)

## ----coefs--------------------------------------------------------------------
# Metafrontier coefficients
coef(fit_det, which = "meta")

# Group-specific coefficients
coef(fit_det, which = "group")

## ----model-info---------------------------------------------------------------
# Log-likelihood (sum of group log-likelihoods for deterministic)
logLik(fit_det)

# Number of observations
nobs(fit_det)

# AIC and BIC (available automatically via logLik method)
AIC(fit_det)

## ----plot-tgr, fig.height=4---------------------------------------------------
plot(fit_det, which = "tgr")

## ----plot-eff, fig.height=4---------------------------------------------------
plot(fit_det, which = "efficiency")

## ----plot-decomp, fig.height=4, fig.width=9-----------------------------------
plot(fit_det, which = "decomposition")

## ----poolability--------------------------------------------------------------
poolability_test(fit_det)

## ----distributions, eval=FALSE------------------------------------------------
# # Half-normal (default): u ~ |N(0, sigma_u^2)|
# fit_hn <- metafrontier(log_y ~ log_x1 + log_x2,
#                        data = sim$data, group = "group",
#                        dist = "hnormal")
# 
# # Truncated normal: u ~ N+(mu, sigma_u^2)
# fit_tn <- metafrontier(log_y ~ log_x1 + log_x2,
#                        data = sim$data, group = "group",
#                        dist = "tnormal")
# 
# # Exponential: u ~ Exp(1/sigma_u)
# fit_exp <- metafrontier(log_y ~ log_x1 + log_x2,
#                         data = sim$data, group = "group",
#                         dist = "exponential")

## ----compare-truth------------------------------------------------------------
# True vs estimated metafrontier coefficients
cbind(
  True = sim$params$beta_meta,
  Estimated = coef(fit_det, which = "meta")
)

# True vs estimated mean TGR by group
true_tgr <- tapply(sim$data$true_tgr, sim$data$group, mean)
est_tgr <- tapply(fit_det$tgr, fit_det$group_vec, mean)
cbind(True = true_tgr, Estimated = est_tgr)

# Correlation between true and estimated efficiency
cor(sim$data$true_te, fit_det$te_group)
cor(sim$data$true_te_star, fit_det$te_meta)

## ----panel-sfa, eval=FALSE----------------------------------------------------
# # Simulate panel data
# panel_sim <- simulate_panel_metafrontier(
#   n_groups = 2, n_firms_per_group = 20, n_periods = 5, seed = 42
# )
# 
# # BC92: time-varying inefficiency u_it = u_i * exp(-eta*(t-T))
# fit_panel <- metafrontier(
#   log_y ~ log_x1 + log_x2,
#   data = panel_sim$data,
#   group = "group",
#   panel = list(id = "firm", time = "year"),
#   panel_dist = "bc92"
# )
# summary(fit_panel)
# 
# # The eta parameter captures time-varying inefficiency
# # eta > 0: inefficiency decreasing over time
# # eta < 0: inefficiency increasing over time

## ----bootstrap, eval=FALSE----------------------------------------------------
# sim <- simulate_metafrontier(n_groups = 2, n_per_group = 100, seed = 42)
# fit <- metafrontier(log_y ~ log_x1 + log_x2, data = sim$data,
#                     group = "group", meta_type = "stochastic")
# 
# # Nonparametric bootstrap (case resampling within groups)
# boot <- boot_tgr(fit, R = 499, type = "nonparametric", seed = 1)
# print(boot)
# 
# # Observation-level CIs
# ci <- confint(boot)
# head(ci)
# 
# # Group-level mean TGR CIs
# boot$ci_group
# 
# # Parametric bootstrap (resample from estimated error distributions)
# boot_par <- boot_tgr(fit, R = 499, type = "parametric", seed = 1)

## ----murphy-topel, eval=FALSE-------------------------------------------------
# fit <- metafrontier(log_y ~ log_x1 + log_x2, data = sim$data,
#                     group = "group", meta_type = "stochastic")
# 
# # Naive (uncorrected) standard errors
# vcov(fit)
# 
# # Murphy-Topel corrected standard errors
# vcov(fit, correction = "murphy-topel")
# 
# # Corrected confidence intervals
# confint(fit, correction = "murphy-topel")

## ----latent-class, eval=FALSE-------------------------------------------------
# sim <- simulate_metafrontier(n_groups = 2, n_per_group = 100, seed = 42)
# 
# # Fit with 2 latent classes
# lc <- latent_class_metafrontier(
#   log_y ~ log_x1 + log_x2,
#   data = sim$data, n_classes = 2,
#   n_starts = 5, seed = 123
# )
# print(lc)
# summary(lc)
# 
# # Select optimal number of classes via BIC
# bic_table <- select_n_classes(
#   log_y ~ log_x1 + log_x2, data = sim$data,
#   n_classes_range = 2:4, n_starts = 3, seed = 42
# )
# print(bic_table)  # choose n_classes with lowest BIC

## ----ddf, eval=FALSE----------------------------------------------------------
# sim <- simulate_metafrontier(n_groups = 2, n_per_group = 50, seed = 42)
# # Use raw (non-log) data for DEA
# sim$data$y <- exp(sim$data$log_y)
# sim$data$x1 <- exp(sim$data$log_x1)
# sim$data$x2 <- exp(sim$data$log_x2)
# 
# fit_ddf <- metafrontier(
#   y ~ x1 + x2, data = sim$data, group = "group",
#   method = "dea", type = "directional", direction = "output"
# )
# summary(fit_ddf)
# 
# # Additive decomposition: beta_meta = beta_group + ddf_tgr
# head(data.frame(
#   beta_meta = fit_ddf$beta_meta,
#   beta_group = fit_ddf$beta_group,
#   ddf_tgr = fit_ddf$ddf_tgr
# ))

