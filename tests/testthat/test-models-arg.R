# Test models= argument for pre-fitted objects

test_that("models= works with internal fit objects", {
  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 60, seed = 60)

  # Fit group models manually using internal function
  f <- Formula::Formula(log_y ~ log_x1 + log_x2)
  g1_data <- sim$data[sim$data$group == "G1", ]
  g2_data <- sim$data[sim$data$group == "G2", ]

  fit_g1 <- metafrontier:::.fit_sfa_group(f, g1_data, "hnormal", list())
  fit_g2 <- metafrontier:::.fit_sfa_group(f, g2_data, "hnormal", list())

  # Pass as models
  result <- metafrontier(models = list(G1 = fit_g1, G2 = fit_g2))

  expect_s3_class(result, "metafrontier")
  expect_length(result$groups, 2)
  expect_equal(unname(result$nobs["total"]),
               nrow(g1_data) + nrow(g2_data))
  expect_true(all(result$tgr > 0 & result$tgr <= 1 + 1e-6))
})


test_that("models= with stochastic metafrontier works", {
  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 60, seed = 61)

  f <- Formula::Formula(log_y ~ log_x1 + log_x2)
  g1_data <- sim$data[sim$data$group == "G1", ]
  g2_data <- sim$data[sim$data$group == "G2", ]

  fit_g1 <- metafrontier:::.fit_sfa_group(f, g1_data, "hnormal", list())
  fit_g2 <- metafrontier:::.fit_sfa_group(f, g2_data, "hnormal", list())

  result <- metafrontier(models = list(G1 = fit_g1, G2 = fit_g2),
                         meta_type = "stochastic")

  expect_s3_class(result, "metafrontier")
  expect_true(!is.null(result$meta_vcov))
})


test_that("models= requires named list", {
  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 30, seed = 62)
  f <- Formula::Formula(log_y ~ log_x1 + log_x2)
  g1 <- metafrontier:::.fit_sfa_group(
    f, sim$data[sim$data$group == "G1", ], "hnormal", list())

  expect_error(metafrontier(models = list(g1)), "named list")
})


test_that("models= requires at least 2 groups", {
  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 30, seed = 63)
  f <- Formula::Formula(log_y ~ log_x1 + log_x2)
  g1 <- metafrontier:::.fit_sfa_group(
    f, sim$data[sim$data$group == "G1", ], "hnormal", list())

  expect_error(metafrontier(models = list(G1 = g1)), "2 groups")
})


test_that("models= with sfaR works", {
  skip_if_not_installed("sfaR")

  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 80, seed = 64)
  g1_data <- sim$data[sim$data$group == "G1", ]
  g2_data <- sim$data[sim$data$group == "G2", ]

  fit_g1 <- sfaR::sfacross(log_y ~ log_x1 + log_x2, data = g1_data)
  fit_g2 <- sfaR::sfacross(log_y ~ log_x1 + log_x2, data = g2_data)

  result <- metafrontier(models = list(G1 = fit_g1, G2 = fit_g2))

  expect_s3_class(result, "metafrontier")
  expect_length(result$groups, 2)
  expect_true(all(result$te_group > 0))
})
