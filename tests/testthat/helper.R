# Shared test setup
# This file is loaded before all tests by testthat

# Tolerances for numerical comparisons
TOLERANCE_COEF <- 1e-4
TOLERANCE_SE   <- 1e-3
TOLERANCE_LL   <- 1e-3
TOLERANCE_EFF  <- 1e-4

# Small simulated dataset for fast tests
set.seed(12345)
test_sim <- simulate_metafrontier(
  n_groups = 2,
  n_per_group = 100,
  n_inputs = 2,
  beta_meta = c(1.0, 0.5, 0.3),
  tech_gap = c(0, 0.3),
  sigma_u = c(0.3, 0.4),
  sigma_v = 0.2,
  seed = 12345
)
test_data <- test_sim$data
test_params <- test_sim$params
