# Tests for panel SFA estimation (P3)

test_that("simulate_panel_metafrontier returns valid data", {
  sim <- simulate_panel_metafrontier(
    n_groups = 2, n_firms_per_group = 20,
    n_periods = 4, eta = 0.05, seed = 42
  )

  expect_true(is.data.frame(sim$data))
  expect_equal(nrow(sim$data), 2 * 20 * 4)
  expect_true("firm" %in% names(sim$data))
  expect_true("year" %in% names(sim$data))
  expect_true("true_te" %in% names(sim$data))
  expect_true(all(sim$data$true_te > 0 & sim$data$true_te <= 1))
})

test_that("panel SFA BC92 converges", {
  sim <- simulate_panel_metafrontier(
    n_groups = 2, n_firms_per_group = 30,
    n_periods = 5, eta = 0.05, seed = 123
  )

  fit <- metafrontier(
    log_y ~ log_x1 + log_x2,
    data = sim$data,
    group = "group",
    meta_type = "deterministic",
    panel = list(id = "firm", time = "year"),
    panel_dist = "bc92"
  )

  expect_s3_class(fit, "metafrontier")
  expect_true(all(fit$te_group > 0 & fit$te_group <= 1))
  expect_true(all(fit$tgr > 0))
})

test_that("BC92 eta recovery is reasonable", {
  sim <- simulate_panel_metafrontier(
    n_groups = 2, n_firms_per_group = 50,
    n_periods = 6, eta = 0.1, sigma_u = 0.4,
    seed = 456
  )

  fit <- metafrontier(
    log_y ~ log_x1 + log_x2,
    data = sim$data,
    group = "group",
    meta_type = "deterministic",
    panel = list(id = "firm", time = "year"),
    panel_dist = "bc92"
  )

  # Check that eta is estimated (exists in group models)
  for (g in fit$groups) {
    gm <- fit$group_models[[g]]
    expect_true("eta" %in% names(gm))
    # Eta should be in a reasonable range (true = 0.1)
    expect_true(abs(gm$eta) < 1.0,
                label = paste("eta for group", g))
  }
})

test_that("TE* = TE x TGR identity holds for panel", {
  sim <- simulate_panel_metafrontier(
    n_groups = 2, n_firms_per_group = 25,
    n_periods = 4, seed = 789
  )

  fit <- metafrontier(
    log_y ~ log_x1 + log_x2,
    data = sim$data,
    group = "group",
    meta_type = "deterministic",
    panel = list(id = "firm", time = "year"),
    panel_dist = "bc92"
  )

  # TE* should equal TE x TGR
  expect_equal(fit$te_meta, fit$te_group * fit$tgr, tolerance = 1e-10)
})

test_that("panel=NULL gives backward-compatible results", {
  # Without panel argument, should work as cross-sectional
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data, group = "group")

  expect_s3_class(fit, "metafrontier")
  expect_true(length(fit$tgr) == nrow(test_data))
})
