# Test print/summary output structure (platform-independent)

test_that("metafrontier print output has correct structure", {
  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 50, seed = 90)
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = sim$data, group = "group",
                      method = "sfa")

  out <- capture.output(print(fit))
  expect_true(any(grepl("Metafrontier Model", out)))
  expect_true(any(grepl("Method:", out)))
  expect_true(any(grepl("Groups:", out)))
  expect_true(any(grepl("Total obs:", out)))
  expect_true(any(grepl("Mean TGR by group:", out)))
  expect_true(any(grepl("G1", out)))
  expect_true(any(grepl("G2", out)))
})


test_that("metafrontier summary output has correct structure", {
  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 50, seed = 90)
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = sim$data, group = "group",
                      method = "sfa")

  out <- capture.output(summary(fit))
  expect_true(any(grepl("Metafrontier Model Summary", out)))
  expect_true(any(grepl("Call:", out)))
  expect_true(any(grepl("Efficiency Decomposition", out)))
  expect_true(any(grepl("Technology Gap Ratio Summary", out)))
  expect_true(any(grepl("Mean_TE", out)))
  expect_true(any(grepl("Mean_TGR", out)))
})


test_that("malmquist_meta print output has correct structure", {
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

  out <- capture.output(print(malm))
  expect_true(any(grepl("Metafrontier Malmquist", out)))
  expect_true(any(grepl("Orientation:", out)))
  expect_true(any(grepl("RTS:", out)))
  expect_true(any(grepl("MPI", out)))
  expect_true(any(grepl("TEC", out)))
  expect_true(any(grepl("TGC", out)))
})


test_that("malmquist_meta summary output has correct structure", {
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

  out <- capture.output(summary(malm))
  expect_true(any(grepl("Three-Way Decomposition", out)))
  expect_true(any(grepl("Technology Gap Ratios", out)))
  expect_true(any(grepl("Group:", out)))
})


test_that("DEA metafrontier print output has correct structure", {
  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 30, seed = 90)
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = sim$data, group = "group",
                      method = "dea")

  out <- capture.output(print(fit))
  expect_true(any(grepl("Metafrontier Model", out)))
  expect_true(any(grepl("dea", out)))
  expect_true(any(grepl("Mean TGR", out)))
})
