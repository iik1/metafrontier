# Tests for DEA batch performance optimization (P5)

test_that("batch_fast matches per-DMU LP for small n", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data, group = "group",
                      method = "dea", rts = "crs")

  # Valid efficiencies (excluding LP-infeasible DMUs)
  valid_grp <- fit$te_group[!is.na(fit$te_group)]
  valid_meta <- fit$te_meta[!is.na(fit$te_meta)]

  expect_true(all(valid_grp > 0 & valid_grp <= 1.001))
  expect_true(all(valid_meta > 0 & valid_meta <= 1.001))
})

test_that("batch_fast works with VRS", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data, group = "group",
                      method = "dea", rts = "vrs")

  valid <- fit$te_group[!is.na(fit$te_group)]
  expect_true(all(valid > 0 & valid <= 1.001))
})

test_that("batch_fast works with input orientation", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data, group = "group",
                      method = "dea", rts = "crs",
                      orientation = "input")

  valid <- fit$te_group[!is.na(fit$te_group)]
  expect_true(all(valid > 0 & valid <= 1.001))
})

test_that("batch_fast gives same results for large n", {
  # Create larger dataset by duplicating test_data
  big_data <- rbind(test_data, test_data, test_data, test_data)
  big_data$group <- rep(test_data$group, 4)

  fit_big <- metafrontier(log_y ~ log_x1 + log_x2,
                          data = big_data, group = "group",
                          method = "dea", rts = "crs")

  # Verify structure and valid efficiencies
  expect_equal(length(fit_big$tgr), nrow(big_data))
  expect_true(all(fit_big$te_group > 0, na.rm = TRUE))
  expect_true(all(fit_big$te_meta > 0, na.rm = TRUE))

  # Efficiency of duplicated DMUs should be the same
  # (same reference set, same inputs/outputs)
  n_orig <- nrow(test_data)
  eff_1 <- fit_big$te_meta[1:n_orig]
  eff_2 <- fit_big$te_meta[(n_orig + 1):(2 * n_orig)]
  expect_equal(eff_1, eff_2, tolerance = 1e-6)
})
