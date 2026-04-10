test_that("very small groups raise error", {
  tiny_sim <- simulate_metafrontier(n_groups = 2, n_per_group = 2, seed = 99)
  expect_error(
    metafrontier(log_y ~ log_x1 + log_x2,
                 data = tiny_sim$data, group = "group")
  )
})

test_that("group as vector works", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data,
                      group = test_data$group)
  expect_s3_class(fit, "metafrontier")
})

test_that("half-normal distribution works", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data,
                      group = "group",
                      dist = "hnormal")
  expect_s3_class(fit, "metafrontier")
  expect_true(all(fit$te_group > 0 & fit$te_group <= 1))
})

test_that("truncated normal distribution works", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data,
                      group = "group",
                      dist = "tnormal")
  expect_s3_class(fit, "metafrontier")
  expect_true(all(fit$te_group > 0 & fit$te_group <= 1))
})

test_that("exponential distribution works", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data,
                      group = "group",
                      dist = "exponential")
  expect_s3_class(fit, "metafrontier")
  expect_true(all(fit$te_group > 0 & fit$te_group <= 1))
})
