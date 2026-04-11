# Test | formula interface for inefficiency determinants

test_that("formula with | works for hnormal", {
  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 80, seed = 50)
  # Add a Z variable correlated with inefficiency
  sim$data$z1 <- sim$data$true_u + rnorm(nrow(sim$data), 0, 0.1)

  fit <- metafrontier(log_y ~ log_x1 + log_x2 | z1,
                      data = sim$data, group = "group",
                      method = "sfa", dist = "hnormal")

  expect_s3_class(fit, "metafrontier")
  expect_length(fit$te_group, nrow(sim$data))
  expect_true(all(fit$te_group > 0 & fit$te_group <= 1))

  # Check delta coefficients are stored in group models
  for (g in fit$groups) {
    gm <- fit$group_models[[g]]
    expect_true(!is.null(gm$delta))
    expect_true(!is.null(gm$Z))
  }
})


test_that("formula with | works for tnormal", {
  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 80, seed = 51)
  sim$data$z1 <- sim$data$true_u + rnorm(nrow(sim$data), 0, 0.1)

  fit <- metafrontier(log_y ~ log_x1 + log_x2 | z1,
                      data = sim$data, group = "group",
                      method = "sfa", dist = "tnormal")

  expect_s3_class(fit, "metafrontier")
  expect_length(fit$te_group, nrow(sim$data))
})


test_that("formula with | works for exponential", {
  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 80, seed = 52)
  sim$data$z1 <- sim$data$true_u + rnorm(nrow(sim$data), 0, 0.1)

  fit <- metafrontier(log_y ~ log_x1 + log_x2 | z1,
                      data = sim$data, group = "group",
                      method = "sfa", dist = "exponential")

  expect_s3_class(fit, "metafrontier")
  expect_length(fit$te_group, nrow(sim$data))
})


test_that("formula without | still works (backward compatibility)", {
  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 60, seed = 53)

  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = sim$data, group = "group",
                      method = "sfa")

  expect_s3_class(fit, "metafrontier")
  for (g in fit$groups) {
    expect_null(fit$group_models[[g]]$delta)
    expect_null(fit$group_models[[g]]$Z)
  }
})


test_that("multiple Z variables work", {
  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 80, seed = 54)
  sim$data$z1 <- rnorm(nrow(sim$data))
  sim$data$z2 <- rnorm(nrow(sim$data))

  fit <- metafrontier(log_y ~ log_x1 + log_x2 | z1 + z2,
                      data = sim$data, group = "group",
                      method = "sfa")

  expect_s3_class(fit, "metafrontier")
  # Z should have 3 columns (intercept + z1 + z2)
  for (g in fit$groups) {
    expect_equal(ncol(fit$group_models[[g]]$Z), 3L)
    expect_length(fit$group_models[[g]]$delta, 3L)
  }
})
