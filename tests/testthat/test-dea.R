# Test DEA metafrontier estimation

test_that("DEA metafrontier returns correct structure", {
  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 30,
                                seed = 99)
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = sim$data, group = "group",
                      method = "dea", rts = "crs")

  expect_s3_class(fit, "metafrontier")
  expect_s3_class(fit, "metafrontier_dea")
  expect_length(fit$te_group, nrow(sim$data))
  expect_length(fit$te_meta, nrow(sim$data))
  expect_length(fit$tgr, nrow(sim$data))
  expect_null(fit$meta_coef)
  expect_null(fit$group_coef)
})


test_that("DEA efficiencies are bounded in (0, 1]", {
  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 40,
                                seed = 77)
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = sim$data, group = "group",
                      method = "dea", rts = "crs")

  expect_true(all(fit$te_group > 0, na.rm = TRUE))
  expect_true(all(fit$te_group <= 1 + 1e-10, na.rm = TRUE))
  expect_true(all(fit$te_meta > 0, na.rm = TRUE))
  expect_true(all(fit$te_meta <= 1 + 1e-10, na.rm = TRUE))
})


test_that("DEA TGR satisfies TE* = TE x TGR", {
  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 40,
                                seed = 88)
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = sim$data, group = "group",
                      method = "dea", rts = "crs")

  expect_equal(fit$te_meta, fit$te_group * fit$tgr,
               tolerance = 1e-10)
})


test_that("DEA TGR <= 1 (pooled cannot be worse than group)", {
  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 40,
                                seed = 55)
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = sim$data, group = "group",
                      method = "dea", rts = "crs")

  expect_true(all(fit$tgr <= 1 + 1e-10, na.rm = TRUE))
})


test_that("DEA with VRS returns different scores than CRS", {
  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 40,
                                seed = 33)
  fit_crs <- metafrontier(log_y ~ log_x1 + log_x2,
                          data = sim$data, group = "group",
                          method = "dea", rts = "crs")
  fit_vrs <- metafrontier(log_y ~ log_x1 + log_x2,
                          data = sim$data, group = "group",
                          method = "dea", rts = "vrs")

  # VRS efficiency >= CRS efficiency (VRS frontier is at least
  # as enveloping)
  expect_true(all(fit_vrs$te_group >= fit_crs$te_group - 1e-6))
})


test_that("DEA with DRS and IRS work without error", {
  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 20,
                                seed = 44)

  expect_no_error(
    metafrontier(log_y ~ log_x1 + log_x2,
                 data = sim$data, group = "group",
                 method = "dea", rts = "drs")
  )

  expect_no_error(
    metafrontier(log_y ~ log_x1 + log_x2,
                 data = sim$data, group = "group",
                 method = "dea", rts = "irs")
  )
})


test_that("DEA with 3 groups works", {
  sim <- simulate_metafrontier(n_groups = 3, n_per_group = 25,
                                seed = 66)
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = sim$data, group = "group",
                      method = "dea", rts = "crs")

  expect_length(fit$groups, 3)
  expect_equal(unname(fit$nobs["total"]), 75L)
})


test_that("DEA S3 methods work", {
  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 30,
                                seed = 11)
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = sim$data, group = "group",
                      method = "dea")

  # print
  expect_output(print(fit), "Metafrontier Model")

  # summary
  s <- summary(fit)
  expect_s3_class(s, "summary.metafrontier")

  # efficiencies
  te <- efficiencies(fit, type = "group")
  te_star <- efficiencies(fit, type = "meta")
  tgr <- efficiencies(fit, type = "tgr")
  expect_length(te, nrow(sim$data))
  expect_equal(te_star, te * tgr, tolerance = 1e-10)

  # nobs
  expect_equal(nobs(fit), nrow(sim$data))

  # tgr_summary
  ts <- tgr_summary(fit)
  expect_equal(nrow(ts), 2)
  expect_true("Mean" %in% names(ts))

  # plot (should not error)
  expect_no_error(plot(fit, which = "tgr"))
  expect_no_error(plot(fit, which = "efficiency"))
  expect_no_error(plot(fit, which = "decomposition"))
})


test_that("DEA single-input model works", {
  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 30,
                                n_inputs = 1, seed = 22)
  fit <- metafrontier(log_y ~ log_x1,
                      data = sim$data, group = "group",
                      method = "dea", rts = "crs")

  expect_s3_class(fit, "metafrontier_dea")
  expect_true(all(fit$te_group > 0, na.rm = TRUE))
})


test_that("DEA input orientation works", {
  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 30,
                                seed = 77)

  fit_out <- metafrontier(log_y ~ log_x1 + log_x2,
                          data = sim$data, group = "group",
                          method = "dea", orientation = "output")
  fit_in <- metafrontier(log_y ~ log_x1 + log_x2,
                         data = sim$data, group = "group",
                         method = "dea", orientation = "input")

  # Both should produce valid efficiency scores
  expect_true(all(fit_out$te_group > 0, na.rm = TRUE))
  expect_true(all(fit_in$te_group > 0, na.rm = TRUE))
  expect_true(all(fit_in$te_group <= 1 + 1e-10, na.rm = TRUE))
})
