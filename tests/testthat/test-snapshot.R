# Snapshot tests for print/summary output stability

test_that("metafrontier print output is stable", {
  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 50, seed = 90)
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = sim$data, group = "group",
                      method = "sfa")

  expect_snapshot(print(fit))
})


test_that("metafrontier summary output is stable", {
  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 50, seed = 90)
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = sim$data, group = "group",
                      method = "sfa")

  expect_snapshot(summary(fit))
})


test_that("malmquist_meta print output is stable", {
  panels <- lapply(1:2, function(t) {
    sim <- simulate_metafrontier(n_groups = 2, n_per_group = 20,
                                  tech_gap = c(0, 0.3 + 0.05 * t),
                                  seed = 90 + t)
    sim$data$time <- t
    sim$data$id <- seq_len(nrow(sim$data))
    sim$data
  })
  pd <- do.call(rbind, panels)

  malm <- malmquist_meta(log_y ~ log_x1 + log_x2,
                         data = pd, group = "group", time = "time")

  expect_snapshot(print(malm))
})


test_that("malmquist_meta summary output is stable", {
  panels <- lapply(1:2, function(t) {
    sim <- simulate_metafrontier(n_groups = 2, n_per_group = 20,
                                  tech_gap = c(0, 0.3 + 0.05 * t),
                                  seed = 90 + t)
    sim$data$time <- t
    sim$data$id <- seq_len(nrow(sim$data))
    sim$data
  })
  pd <- do.call(rbind, panels)

  malm <- malmquist_meta(log_y ~ log_x1 + log_x2,
                         data = pd, group = "group", time = "time")

  expect_snapshot(summary(malm))
})


test_that("DEA metafrontier print output is stable", {
  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 30, seed = 90)
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = sim$data, group = "group",
                      method = "dea")

  expect_snapshot(print(fit))
})
