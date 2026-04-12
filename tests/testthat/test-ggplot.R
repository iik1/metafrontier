# Tests for ggplot2 autoplot methods (P9)

skip_if_not_installed("ggplot2")
library(ggplot2)

test_that("autoplot.metafrontier produces ggplot for tgr", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data, group = "group",
                      meta_type = "stochastic")

  p <- autoplot(fit, which = "tgr")
  expect_s3_class(p, "gg")
})

test_that("autoplot.metafrontier produces ggplot for efficiency", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data, group = "group")

  p <- autoplot(fit, which = "efficiency")
  expect_s3_class(p, "gg")
})

test_that("autoplot.metafrontier produces ggplot for decomposition", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data, group = "group")

  p <- autoplot(fit, which = "decomposition")
  expect_s3_class(p, "gg")
})

test_that("autoplot.metafrontier produces ggplot for frontier", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data, group = "group",
                      meta_type = "stochastic")

  p <- autoplot(fit, which = "frontier")
  expect_s3_class(p, "gg")
})

test_that("autoplot.malmquist_meta works", {
  set.seed(42)
  panels <- lapply(1:3, function(t) {
    sim <- simulate_metafrontier(
      n_groups = 2, n_per_group = 30,
      tech_gap = c(0, 0.3 + 0.05 * t),
      sigma_u = c(0.2, 0.3),
      seed = 42 + t
    )
    sim$data$time <- t
    sim$data$id <- seq_len(nrow(sim$data))
    sim$data
  })
  panel_data <- do.call(rbind, panels)

  malm <- malmquist_meta(
    log_y ~ log_x1 + log_x2,
    data = panel_data,
    group = "group",
    time = "time"
  )

  p1 <- autoplot(malm, which = "decomposition")
  expect_s3_class(p1, "gg")

  p2 <- autoplot(malm, which = "tgr_evolution")
  expect_s3_class(p2, "gg")

  p3 <- autoplot(malm, which = "mpi_trend")
  expect_s3_class(p3, "gg")
})

test_that("autoplot.boot_tgr works", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data, group = "group",
                      meta_type = "stochastic")

  boot <- boot_tgr(fit, R = 10, seed = 42, progress = FALSE)

  p1 <- autoplot(boot, which = "distribution")
  expect_s3_class(p1, "gg")

  p2 <- autoplot(boot, which = "ci")
  expect_s3_class(p2, "gg")
})
