# DGP extensions: group-specific slopes, input distributions, attrition
# (v0.3.0). Legacy reference values were generated with
# the pre-0.3.0 implementation.

test_that("simulate_metafrontier reproduces legacy output for seed 42", {
  sim <- simulate_metafrontier(seed = 42)

  expect_equal(nrow(sim$data), 200L)
  expect_equal(ncol(sim$data), 9L)
  expect_equal(sim$data$log_x1[1], 4.574030217481777, tolerance = 1e-12)
  expect_equal(sim$data$log_y[1], 3.553174757064474, tolerance = 1e-12)
  expect_equal(sim$data$true_te[57], 0.801290099145637, tolerance = 1e-12)
  expect_equal(sim$data$log_y[200], 1.470010450855711, tolerance = 1e-12)
  expect_equal(sim$data$true_te_star[150], 0.401389206182210,
               tolerance = 1e-12)
})

test_that("beta_groups with differing slopes gives per-observation TGR in (0, 1]", {
  bg <- rbind(c(1.0, 0.5, 0.2),
              c(0.9, 0.6, 0.1))
  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 200,
                               beta_groups = bg, seed = 1)

  expect_true(all(sim$data$true_tgr > 0 & sim$data$true_tgr <= 1))

  # TGR varies within each group when slopes differ
  sds <- tapply(sim$data$true_tgr, sim$data$group, stats::sd)
  expect_true(all(sds > 0))

  # No single log-linear metafrontier: beta_meta is NULL,
  # beta_groups is returned
  expect_null(sim$params$beta_meta)
  expect_equal(sim$params$beta_groups$G1, c(1.0, 0.5, 0.2))
  expect_equal(sim$params$beta_groups$G2, c(0.9, 0.6, 0.1))

  # Decomposition identity still holds
  expect_equal(sim$data$true_te_star,
               sim$data$true_te * sim$data$true_tgr,
               tolerance = 1e-10)
})

test_that("beta_groups replicating the intercept-only design gives exp(-gap)", {
  bg <- rbind(c(1.0, 0.5, 0.2),
              c(0.75, 0.5, 0.2))
  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 100,
                               beta_groups = bg, seed = 7)

  tgr1 <- sim$data$true_tgr[sim$data$group == "G1"]
  tgr2 <- sim$data$true_tgr[sim$data$group == "G2"]
  expect_equal(tgr1, rep(1, 100), tolerance = 1e-10)
  expect_equal(tgr2, rep(exp(-0.25), 100), tolerance = 1e-10)
})

test_that("beta_groups accepts a list and warns when tech_gap is also given", {
  bg_list <- list(c(1.0, 0.5, 0.2), c(0.9, 0.6, 0.1))
  sim <- simulate_metafrontier(n_groups = 2, beta_groups = bg_list, seed = 2)
  expect_equal(sim$params$beta_groups$G2, c(0.9, 0.6, 0.1))

  expect_warning(
    simulate_metafrontier(n_groups = 2, beta_groups = bg_list,
                          tech_gap = c(0, 0.3), seed = 2),
    "tech_gap"
  )
  expect_error(
    simulate_metafrontier(n_groups = 2,
                          beta_groups = rbind(c(1, 0.5, 0.2)),
                          seed = 2),
    "beta_groups"
  )
})

test_that("input_corr induces the target input correlation", {
  R <- matrix(c(1, 0.6,
                0.6, 1), nrow = 2)
  sim <- simulate_metafrontier(n_groups = 1, n_per_group = 2000,
                               input_corr = R, seed = 3)

  realised <- stats::cor(sim$data$log_x1, sim$data$log_x2)
  expect_lt(abs(realised - 0.6), 0.15)

  # Spread roughly matches the legacy uniform draws
  expect_lt(abs(stats::sd(sim$data$log_x1) - 5 / sqrt(12)), 0.15)
})

test_that("input_corr is validated", {
  not_psd <- matrix(c(1, 2,
                      2, 1), nrow = 2)
  expect_error(simulate_metafrontier(input_corr = not_psd, seed = 1),
               "positive definite")
  not_sym <- matrix(c(1, 0.2,
                      0.6, 1), nrow = 2)
  expect_error(simulate_metafrontier(input_corr = not_sym, seed = 1),
               "symmetric")
})

test_that("input_means centres per-group log-inputs", {
  m <- rbind(c(1, 2),
             c(3, 4))
  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 500,
                               input_means = m, seed = 4)

  for (g in 1:2) {
    d <- sim$data[sim$data$group == paste0("G", g), ]
    expect_lt(abs(mean(d$log_x1) - m[g, 1]), 0.15)
    expect_lt(abs(mean(d$log_x2) - m[g, 2]), 0.15)
  }

  expect_error(simulate_metafrontier(n_groups = 2,
                                     input_means = rbind(c(1, 2)),
                                     seed = 4),
               "input_means")
})

test_that("attrition = 0 reproduces the legacy balanced panel", {
  sim <- simulate_panel_metafrontier(seed = 99)

  expect_equal(nrow(sim$data), 300L)
  expect_equal(sim$data$log_y[1], 2.103142882313839, tolerance = 1e-12)
  expect_equal(sim$data$log_y[300], 0.928893291616277, tolerance = 1e-12)
  expect_equal(sim$data$true_tgr[250], 0.606530659712633, tolerance = 1e-12)

  sim0 <- simulate_panel_metafrontier(seed = 99, attrition = 0)
  expect_identical(sim$data, sim0$data)
  expect_equal(sim0$params$attrition_share, 0)
})

test_that("attrition yields an unbalanced panel with all firms in period 1", {
  sim <- simulate_panel_metafrontier(seed = 123, attrition = 0.3)
  d <- sim$data
  n_firms <- 2 * 30

  expect_lt(nrow(d), 300)
  expect_equal(length(unique(d$firm)), n_firms)
  expect_equal(sum(d$year == 1), n_firms)

  # Unbalanced: firms have differing numbers of periods
  periods_per_firm <- table(d$firm)
  expect_gt(length(unique(as.integer(periods_per_firm))), 1L)

  # Realised share stored and consistent with the rows dropped
  expect_equal(sim$params$attrition_share,
               (300 - nrow(d)) / (300 - n_firms))
  expect_gt(sim$params$attrition_share, 0)

  expect_error(simulate_panel_metafrontier(seed = 1, attrition = 0.9),
               "attrition")
})
