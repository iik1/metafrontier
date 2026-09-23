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

# Two groups with very different TGR levels in unsorted row order.
# Group B's technology is half of group A's, so TGR is about 1 in A and
# about 0.5 in B.
unsorted_two_level_data <- function() {
  set.seed(2024)
  grp <- rep(c("A", "B"), c(25, 35))
  n <- length(grp)
  x1 <- runif(n, 1, 10)
  x2 <- runif(n, 1, 10)
  y <- ifelse(grp == "A", 1, 0.5) * sqrt(x1 * x2) *
    exp(-abs(rnorm(n, 0, 0.3)))
  d <- data.frame(y = y, x1 = x1, x2 = x2, group = grp)
  d[sample(n), ]
}

test_that("nonparametric group intervals use position blocks", {
  d <- unsorted_two_level_data()
  d_sorted <- d[order(d$group), ]
  fit <- metafrontier(y ~ x1 + x2, data = d, group = "group",
                      method = "dea")
  fit_sorted <- metafrontier(y ~ x1 + x2, data = d_sorted,
                             group = "group", method = "dea")

  # Same seed and same within-group row order give the same resampled
  # data, so the bootstrap draws are identical
  b <- boot_tgr(fit, R = 49, type = "nonparametric", seed = 7,
                progress = FALSE)
  b_sorted <- boot_tgr(fit_sorted, R = 49, type = "nonparametric",
                       seed = 7, progress = FALSE)
  expect_equal(b$tgr_boot, b_sorted$tgr_boot)
  expect_equal(b$boot_group, rep(c("A", "B"), c(25, 35)))

  # Group intervals must not depend on the row order of the data
  expect_equal(b$ci_group, b_sorted$ci_group)
  expect_equal(b$ci_group_median, b_sorted$ci_group_median)

  # Each interval covers its group estimate, and the intervals separate
  # the two technology levels
  for (tab in list(b$ci_group, b$ci_group_median)) {
    expect_true(all(tab[[3]] <= tab[[2]] + 1e-8))
    expect_true(all(tab[[2]] <= tab[[4]] + 1e-8))
    expect_lt(tab[tab$Group == "B", 4], tab[tab$Group == "A", 3])
  }
})

test_that("nonparametric bootstrap reports no observation-level CIs", {
  fit <- metafrontier(y ~ x1 + x2, data = unsorted_two_level_data(),
                      group = "group", method = "dea")
  boot <- boot_tgr(fit, R = 5, type = "nonparametric", seed = 1,
                   progress = FALSE)

  expect_null(boot$ci)
  expect_error(confint(boot), "not available for the nonparametric")
  expect_named(boot$ci_group_median,
               c("Group", "Median_TGR", "2.5%", "97.5%"))
  expect_output(print(boot), "median TGR")
})

test_that("nonparametric bootstrap skips rows dropped for missing values", {
  d <- test_data
  d$log_x1[c(5, 200)] <- NA
  fit <- metafrontier(log_y ~ log_x1 + log_x2, data = d, group = "group",
                      meta_type = "deterministic")

  boot <- boot_tgr(fit, R = 5, type = "nonparametric", seed = 1,
                   progress = FALSE)

  expect_true(boot$R_effective > 0)
  expect_equal(ncol(boot$tgr_boot), length(fit$tgr))
  expect_equal(length(boot$boot_group), length(fit$tgr))
  expect_equal(length(boot$group_vec), length(fit$tgr))
})
