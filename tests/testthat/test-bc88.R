# Tests for the BC88 (Battese and Coelli, 1988) efficiency estimator

.bc88_test_data <- function(n = 200, seed = 1) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  u <- abs(rnorm(n, 0, 0.3))
  v <- rnorm(n, 0, 0.2)
  data.frame(y = 1 + 0.5 * x1 + 0.3 * x2 + v - u, x1 = x1, x2 = x2)
}

test_that("bc88 and jlms are both stored, bounded, and highly correlated", {
  d <- .bc88_test_data()
  f <- Formula::Formula(y ~ x1 + x2)

  gm <- metafrontier:::.fit_sfa_group(f, d, "hnormal", list())

  expect_true(all(gm$efficiency_bc88 > 0 & gm$efficiency_bc88 < 1))
  expect_true(all(gm$efficiency_jlms > 0 & gm$efficiency_jlms < 1))
  expect_gt(cor(gm$efficiency_bc88, gm$efficiency_jlms), 0.99)
})

test_that("bc88 matches the textbook formula recomputed by hand", {
  d <- .bc88_test_data()
  f <- Formula::Formula(y ~ x1 + x2)

  gm <- metafrontier:::.fit_sfa_group(f, d, "hnormal", list())

  # Half-normal conditional posterior of u given eps
  sv <- gm$sigma_v
  su <- gm$sigma_u
  eps <- gm$residuals
  s2 <- sv^2 + su^2
  mu_star <- -eps * su^2 / s2
  sigma_star <- sv * su / sqrt(s2)

  te_hand <- exp(-mu_star + 0.5 * sigma_star^2) *
    pnorm(mu_star / sigma_star - sigma_star) / pnorm(mu_star / sigma_star)

  expect_equal(gm$efficiency_bc88, as.numeric(te_hand), tolerance = 1e-10)
})

test_that("bc88 and jlms agree with sfaR::sfacross", {
  skip_if_not_installed("sfaR")

  d <- .bc88_test_data()
  f <- Formula::Formula(y ~ x1 + x2)

  gm <- metafrontier:::.fit_sfa_group(f, d, "hnormal", list())
  sm <- sfaR::sfacross(y ~ x1 + x2, data = d, udist = "hnormal")
  effs <- sfaR::efficiencies(sm)

  expect_gt(cor(gm$efficiency_bc88, effs$teBC), 0.999)
  expect_lt(max(abs(gm$efficiency_bc88 - effs$teBC)), 1e-3)
  expect_gt(cor(gm$efficiency_jlms, effs$teJLMS), 0.999)
  expect_lt(max(abs(gm$efficiency_jlms - effs$teJLMS)), 1e-3)
})

test_that("estimator argument selects the efficiency vector, default bc88", {
  d <- .bc88_test_data()
  f <- Formula::Formula(y ~ x1 + x2)

  gm_def <- metafrontier:::.fit_sfa_group(f, d, "hnormal", list())
  gm_jlms <- metafrontier:::.fit_sfa_group(f, d, "hnormal", list(),
                                           estimator = "jlms")

  expect_identical(gm_def$estimator, "bc88")
  expect_identical(gm_def$efficiency, gm_def$efficiency_bc88)
  expect_identical(gm_jlms$estimator, "jlms")
  expect_identical(gm_jlms$efficiency, gm_jlms$efficiency_jlms)
})

test_that("bc88 is computed in all distribution branches", {
  d <- .bc88_test_data()
  set.seed(2)
  d$z1 <- rnorm(nrow(d))
  f <- Formula::Formula(y ~ x1 + x2)
  fz <- Formula::Formula(y ~ x1 + x2 | z1)

  for (dist in c("hnormal", "tnormal", "exponential")) {
    gm <- metafrontier:::.fit_sfa_group(f, d, dist, list())
    expect_true(all(gm$efficiency_bc88 > 0 & gm$efficiency_bc88 <= 1),
                label = paste("homoscedastic", dist))
    gmz <- metafrontier:::.fit_sfa_group(fz, d, dist, list())
    expect_true(all(gmz$efficiency_bc88 > 0 & gmz$efficiency_bc88 <= 1),
                label = paste("heteroscedastic", dist))
  }
})
