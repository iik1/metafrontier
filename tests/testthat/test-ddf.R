# Tests for Directional Distance Functions (P8)

test_that("DDF metafrontier produces valid output", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data, group = "group",
                      method = "dea", type = "directional",
                      direction = "proportional")

  expect_s3_class(fit, "metafrontier")
  expect_true(length(fit$tgr) == nrow(test_data))
  valid <- fit$te_group[!is.na(fit$te_group)]
  expect_true(all(valid > 0 & valid <= 1.001))
})

test_that("DDF with output direction works", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data, group = "group",
                      method = "dea", type = "directional",
                      direction = "output")

  expect_s3_class(fit, "metafrontier")
  expect_true(length(fit$tgr) == nrow(test_data))
})

test_that("DDF beta values are non-negative", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data, group = "group",
                      method = "dea", type = "directional",
                      direction = "proportional")

  valid_grp <- fit$beta_group[!is.na(fit$beta_group)]
  expect_true(all(valid_grp >= -0.001))

  valid_meta <- fit$beta_meta[!is.na(fit$beta_meta)]
  expect_true(all(valid_meta >= -0.001))
})

test_that("DDF TGR additive decomposition", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data, group = "group",
                      method = "dea", type = "directional",
                      direction = "proportional")

  # Additive TGR: beta_meta - beta_group
  valid <- !is.na(fit$ddf_tgr) & !is.na(fit$beta_meta) &
           !is.na(fit$beta_group)
  if (sum(valid) > 0) {
    expect_equal(fit$ddf_tgr[valid],
                 (fit$beta_meta - fit$beta_group)[valid],
                 tolerance = 1e-10)
  }
})

test_that("radial DEA still works (backward compat)", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data, group = "group",
                      method = "dea", type = "radial")

  expect_s3_class(fit, "metafrontier")
  expect_true(length(fit$tgr) == nrow(test_data))
})
