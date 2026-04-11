# Test confint() and predict() methods

test_that("confint works for stochastic metafrontier", {
  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 200, seed = 70)

  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = sim$data, group = "group",
                      method = "sfa", meta_type = "stochastic")

  ci <- confint(fit)
  expect_true(is.matrix(ci))
  expect_equal(ncol(ci), 2)
  expect_equal(nrow(ci), length(fit$meta_coef))

  # Lower bound <= upper bound
  expect_true(all(ci[, 1] <= ci[, 2]))
})


test_that("confint respects level argument", {
  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 200, seed = 71)

  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = sim$data, group = "group",
                      method = "sfa", meta_type = "stochastic")

  ci_95 <- confint(fit, level = 0.95)
  ci_99 <- confint(fit, level = 0.99)

  # 99% CI should be wider than 95% CI
  widths_95 <- ci_95[, 2] - ci_95[, 1]
  widths_99 <- ci_99[, 2] - ci_99[, 1]
  expect_true(all(widths_99 > widths_95))
})


test_that("confint errors for deterministic metafrontier", {
  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 100, seed = 72)

  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = sim$data, group = "group",
                      method = "sfa", meta_type = "deterministic")

  expect_error(confint(fit), "stochastic")
})


test_that("confint parm subsetting works", {
  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 200, seed = 73)

  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = sim$data, group = "group",
                      method = "sfa", meta_type = "stochastic")

  ci_sub <- confint(fit, parm = 1:2)
  expect_equal(nrow(ci_sub), 2)
})


test_that("predict returns training predictions without newdata", {
  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 100, seed = 74)

  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = sim$data, group = "group",
                      method = "sfa")

  pred_meta <- predict(fit, type = "meta")
  pred_group <- predict(fit, type = "group")

  expect_length(pred_meta, nrow(sim$data))
  expect_length(pred_group, nrow(sim$data))
  expect_equal(pred_meta, fit$meta_frontier)
  expect_equal(pred_group, fit$group_frontier)
})


test_that("predict with newdata works for SFA", {
  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 100, seed = 75)

  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = sim$data, group = "group",
                      method = "sfa")

  newdata <- data.frame(log_x1 = c(1, 2, 3), log_x2 = c(2, 3, 4))
  pred <- predict(fit, newdata = newdata, type = "meta")

  expect_length(pred, 3)
  expect_true(all(is.finite(pred)))
})


test_that("predict errors for DEA with newdata", {
  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 50, seed = 76)

  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = sim$data, group = "group",
                      method = "dea")

  newdata <- data.frame(log_x1 = 1, log_x2 = 2)
  expect_error(predict(fit, newdata = newdata), "SFA")
})
