test_that("metafrontier SFA estimation returns correct class", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data,
                      group = "group",
                      method = "sfa",
                      meta_type = "deterministic")

  expect_s3_class(fit, "metafrontier")
  expect_s3_class(fit, "metafrontier_sfa")
})

test_that("metafrontier SFA returns expected components", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data,
                      group = "group")

  expect_true(!is.null(fit$meta_coef))
  expect_true(!is.null(fit$group_coef))
  expect_true(!is.null(fit$tgr))
  expect_true(!is.null(fit$te_group))
  expect_true(!is.null(fit$te_meta))
  expect_true(!is.null(fit$groups))
  expect_equal(length(fit$groups), 2)
  expect_equal(as.numeric(fit$nobs["total"]), nrow(test_data))
})

test_that("deterministic metafrontier envelops group frontiers", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data,
                      group = "group",
                      meta_type = "deterministic")

  # TGR should be <= 1 (metafrontier weakly dominates)
  expect_true(all(fit$tgr <= 1 + 1e-6))

  # Metafrontier values should be >= group frontier values
  expect_true(all(fit$meta_frontier >= fit$group_frontier - 1e-6))
})

test_that("efficiency decomposition identity holds: TE* = TE x TGR", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data,
                      group = "group")

  expect_equal(fit$te_meta, fit$te_group * fit$tgr,
               tolerance = 1e-10)
})

test_that("efficiency scores are bounded in (0, 1]", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data,
                      group = "group")

  expect_true(all(fit$te_group > 0 & fit$te_group <= 1))
  expect_true(all(fit$te_meta > 0 & fit$te_meta <= 1))
  expect_true(all(fit$tgr > 0 & fit$tgr <= 1 + 1e-6))
})

test_that("stochastic metafrontier produces standard errors", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data,
                      group = "group",
                      meta_type = "stochastic")

  expect_true(!is.null(fit$meta_vcov))
  expect_true(!is.null(fit$meta_logLik))
  expect_true(is.finite(fit$meta_logLik))
})

test_that("convergence is achieved", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data,
                      group = "group")

  # Check group model convergence
  for (g in fit$groups) {
    expect_equal(fit$group_models[[g]]$convergence, 0L)
  }

  # Meta convergence
  expect_equal(fit$meta_convergence, 0L)
})

test_that("input validation works", {
  expect_error(metafrontier(), "formula.*data.*group")
  expect_error(
    metafrontier(log_y ~ log_x1, data = test_data, group = "nonexistent"),
    "not found"
  )
})

test_that("single group raises error", {
  d <- test_data[test_data$group == "G1", ]
  d$group <- factor(d$group)  # drop unused level
  expect_error(
    metafrontier(log_y ~ log_x1 + log_x2, data = d, group = "group"),
    "At least 2 groups"
  )
})
