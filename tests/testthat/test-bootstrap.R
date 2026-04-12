# Tests for bootstrap TGR confidence intervals (P2)

test_that("boot_tgr returns correct class and structure", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data, group = "group",
                      meta_type = "stochastic")

  boot <- boot_tgr(fit, R = 10, seed = 42, progress = FALSE)

  expect_s3_class(boot, "boot_tgr")
  expect_true(is.matrix(boot$tgr_boot))
  expect_equal(ncol(boot$tgr_boot), length(fit$tgr))
  expect_true(boot$R_effective > 0)
  expect_true(is.matrix(boot$ci))
  expect_equal(ncol(boot$ci), 2)
  expect_s3_class(boot$ci_group, "data.frame")
})

test_that("seed ensures reproducibility", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data, group = "group",
                      meta_type = "stochastic")

  b1 <- boot_tgr(fit, R = 5, seed = 123, progress = FALSE)
  b2 <- boot_tgr(fit, R = 5, seed = 123, progress = FALSE)

  expect_equal(b1$tgr_boot, b2$tgr_boot)
})

test_that("nonparametric bootstrap works", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data, group = "group",
                      meta_type = "stochastic")

  boot <- boot_tgr(fit, R = 5, type = "nonparametric",
                    seed = 42, progress = FALSE)

  expect_s3_class(boot, "boot_tgr")
  expect_equal(boot$type, "nonparametric")
  expect_true(boot$R_effective > 0)
})

test_that("nonparametric bootstrap works for DEA", {
  fit_dea <- metafrontier(log_y ~ log_x1 + log_x2,
                          data = test_data, group = "group",
                          method = "dea")

  boot <- boot_tgr(fit_dea, R = 5, type = "nonparametric",
                    seed = 42, progress = FALSE)

  expect_s3_class(boot, "boot_tgr")
  expect_true(boot$R_effective > 0)
})

test_that("parametric bootstrap errors for DEA", {
  fit_dea <- metafrontier(log_y ~ log_x1 + log_x2,
                          data = test_data, group = "group",
                          method = "dea")

  expect_error(
    boot_tgr(fit_dea, R = 5, type = "parametric", progress = FALSE),
    "Parametric bootstrap is not available for DEA"
  )
})

test_that("confint.boot_tgr returns correct structure", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data, group = "group",
                      meta_type = "stochastic")

  boot <- boot_tgr(fit, R = 10, seed = 42, progress = FALSE)
  ci <- confint(boot)

  expect_true(is.matrix(ci))
  expect_equal(nrow(ci), length(fit$tgr))
  expect_equal(ncol(ci), 2)
})

test_that("print.boot_tgr produces output", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data, group = "group",
                      meta_type = "stochastic")

  boot <- boot_tgr(fit, R = 5, seed = 42, progress = FALSE)

  expect_output(print(boot), "Bootstrap TGR")
  expect_output(print(boot), "Replications")
})

test_that("deterministic metafrontier bootstrap works", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data, group = "group",
                      meta_type = "deterministic")

  boot <- boot_tgr(fit, R = 5, seed = 42, progress = FALSE)

  expect_s3_class(boot, "boot_tgr")
  expect_true(boot$R_effective > 0)
})
