test_that("simulate_metafrontier returns correct structure", {
  sim <- simulate_metafrontier(n_groups = 3, n_per_group = 50, seed = 1)

  expect_type(sim, "list")
  expect_named(sim, c("data", "params"))

  # Data frame structure
  expect_s3_class(sim$data, "data.frame")
  expect_equal(nrow(sim$data), 150)
  expect_true(all(c("log_x1", "log_x2", "log_y", "group",
                     "true_te", "true_tgr") %in% names(sim$data)))

  # Groups
  expect_equal(nlevels(sim$data$group), 3)
  expect_equal(as.numeric(table(sim$data$group)), rep(50, 3))
})

test_that("simulate_metafrontier respects seed", {
  sim1 <- simulate_metafrontier(seed = 42)
  sim2 <- simulate_metafrontier(seed = 42)
  expect_identical(sim1$data, sim2$data)
})

test_that("simulate_metafrontier produces valid efficiency scores", {
  sim <- simulate_metafrontier(seed = 1)
  expect_true(all(sim$data$true_te > 0 & sim$data$true_te <= 1))
  expect_true(all(sim$data$true_tgr > 0 & sim$data$true_tgr <= 1))
  expect_true(all(sim$data$true_te_star > 0 & sim$data$true_te_star <= 1))
})

test_that("simulate_metafrontier decomposition identity holds", {
  sim <- simulate_metafrontier(seed = 1)
  expect_equal(sim$data$true_te_star,
               sim$data$true_te * sim$data$true_tgr,
               tolerance = 1e-10)
})

test_that("simulate_metafrontier handles per-group sample sizes", {
  sim <- simulate_metafrontier(n_groups = 3,
                               n_per_group = c(50, 100, 200),
                               seed = 1)
  expect_equal(nrow(sim$data), 350)
  expect_equal(as.numeric(table(sim$data$group)), c(50, 100, 200))
})

test_that("simulate_metafrontier validates inputs", {
  expect_error(simulate_metafrontier(beta_meta = c(1, 2)),
               "n_inputs \\+ 1")
  expect_error(simulate_metafrontier(tech_gap = c(0.1)),
               "n_groups")
  expect_error(simulate_metafrontier(sigma_u = c(0.3)),
               "n_groups")
})
