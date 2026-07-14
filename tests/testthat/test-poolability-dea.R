# Test the DEA permutation poolability test and the data.name fix

test_that("permutation test rejects poolability when groups differ", {
  skip_on_cran()

  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 60,
                               tech_gap = c(0, 0.6), seed = 1)
  fit_dea <- metafrontier(log_y ~ log_x1 + log_x2,
                          data = sim$data, group = "group",
                          method = "dea", rts = "vrs")

  pt <- poolability_test(fit_dea, B = 99, seed = 123)

  expect_lt(pt$p.value, 0.05)
})


test_that("permutation test does not reject under a pooled technology", {
  skip_on_cran()

  # Sanity check under H0 with a fixed seed, not a power study
  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 60,
                               tech_gap = c(0, 0), seed = 2)
  fit_dea <- metafrontier(log_y ~ log_x1 + log_x2,
                          data = sim$data, group = "group",
                          method = "dea", rts = "vrs")

  pt <- poolability_test(fit_dea, B = 99, seed = 123)

  expect_gt(pt$p.value, 0.05)
})


test_that("permutation test returns a well-formed htest object", {
  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 25,
                               tech_gap = c(0, 0.4), seed = 7)
  fit_dea <- metafrontier(log_y ~ log_x1 + log_x2,
                          data = sim$data, group = "group",
                          method = "dea", rts = "vrs")

  pt <- poolability_test(fit_dea, B = 49, seed = 9)

  expect_s3_class(pt, "htest")
  expect_identical(
    pt$method,
    "Permutation test for poolability of group frontiers (DEA)"
  )
  expect_named(pt$parameter, "B")
  expect_lte(pt$parameter[["B"]], 49)
  expect_named(pt$statistic, "mean technology gap")
  expect_equal(pt$statistic[["mean technology gap"]],
               mean(1 - fit_dea$tgr, na.rm = TRUE))
  expect_true(pt$p.value > 0 && pt$p.value <= 1)
})


test_that("data.name is the symbol passed for DEA objects", {
  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 25,
                               tech_gap = c(0, 0.4), seed = 7)
  fit_dea <- metafrontier(log_y ~ log_x1 + log_x2,
                          data = sim$data, group = "group",
                          method = "dea", rts = "vrs")

  pt <- poolability_test(fit_dea, B = 9, seed = 1)

  expect_identical(pt$data.name, "fit_dea")
})


test_that("SFA objects still take the LR branch with fixed data.name", {
  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 60,
                               tech_gap = c(0, 0.4), seed = 5)
  fit_sfa <- metafrontier(log_y ~ log_x1 + log_x2,
                          data = sim$data, group = "group",
                          method = "sfa")

  pt <- poolability_test(fit_sfa)

  expect_s3_class(pt, "htest")
  expect_identical(
    pt$method,
    "Likelihood Ratio Test for Poolability of Group Frontiers"
  )
  expect_named(pt$statistic, "LR")
  expect_identical(pt$data.name, "fit_sfa")
})


test_that("permutation test is reproducible with a fixed seed", {
  skip_on_cran()

  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 25,
                               tech_gap = c(0, 0.4), seed = 7)
  fit_dea <- metafrontier(log_y ~ log_x1 + log_x2,
                          data = sim$data, group = "group",
                          method = "dea", rts = "vrs")

  pt1 <- poolability_test(fit_dea, B = 49, seed = 42)
  pt2 <- poolability_test(fit_dea, B = 49, seed = 42)

  expect_identical(pt1$p.value, pt2$p.value)
})


test_that("non-metafrontier objects give an informative error", {
  expect_error(
    poolability_test(list()),
    "SFA- or DEA-based metafrontier"
  )
})
