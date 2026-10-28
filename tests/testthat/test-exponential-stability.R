# Numerical stability of the normal/exponential SFA likelihood and the
# flagging of numerically failed fits (code 99).
#
# Regression: on test_data group G2 the direct log-density lost all
# precision for large sigma_v/sigma_u, and BFGS climbed the rounding
# error to sigma_v ~ 4e91, sigma_u ~ 5e34 and logLik ~ 2e100, reported
# with convergence code 0.

fit_exp <- metafrontier(log_y ~ log_x1 + log_x2, data = test_data,
                        group = "group", meta_type = "deterministic",
                        dist = "exponential")


test_that("exponential log-density matches numerical integration across the switch", {
  # f(eps) = int_0^Inf phi((eps + u)/sigma_v)/sigma_v * dexp(u, 1/sigma_u) du;
  # with sigma_v = 1, sigma_u = 1/40 these eps give a from 35 to 45
  eps <- seq(-5, 5, by = 0.5)
  ref <- vapply(eps, function(e) {
    log(integrate(function(u) dnorm(e + u) * dexp(u, rate = 40),
                  0, 2, rel.tol = 1e-12)$value)
  }, numeric(1))

  ll <- metafrontier:::.ll_exponential(eps, 1, 1 / 40)
  expect_lt(max(abs(ll - ref)), 1e-9)
})


test_that("exponential log-density tends to the normal as sigma_u/sigma_v -> 0", {
  # u ~ Exp(mean sigma_u) with sigma_u negligible: log f(eps) is the
  # N(-sigma_u, sigma_v^2) log-density to within O((sigma_u/sigma_v)^2).
  # The direct form was off by ~0.4 per observation at ratio 1e-8.
  eps <- seq(-3, 3, by = 0.5)
  expect_equal(metafrontier:::.ll_exponential(eps, 1, 1e-8),
               dnorm(eps + 1e-8, log = TRUE), tolerance = 1e-12)

  # The parameter values G2 diverged to (true value ~ -212, not ~ +7e97)
  sv <- 4.441459e+91
  su <- 4.613645e+34
  expect_equal(metafrontier:::.ll_exponential(c(-0.3, 0, 0.3), sv, su),
               dnorm(c(-0.3, 0, 0.3) + su, sd = sv, log = TRUE),
               tolerance = 1e-12)

  # Vector sigma_u (heteroscedastic model) mixing both branches
  ll_vec <- metafrontier:::.ll_exponential(c(0.2, 0.2), 1, c(1e-8, 0.5))
  expect_equal(ll_vec,
               c(metafrontier:::.ll_exponential(0.2, 1, 1e-8),
                 metafrontier:::.ll_exponential(0.2, 1, 0.5)))
})


test_that("internal exponential SFA matches sfaR::sfacross on test_data", {
  skip_if_not_installed("sfaR")

  for (g in c("G1", "G2")) {
    gm <- fit_exp$group_models[[g]]
    sm <- sfaR::sfacross(log_y ~ log_x1 + log_x2, udist = "exponential",
                         data = test_data[test_data$group == g, ], S = 1)
    cf <- coef(sm)

    expect_equal(gm$logLik, sm$mlLoglik, tolerance = TOLERANCE_LL)
    expect_equal(unname(gm$coefficients), unname(cf[1:3]),
                 tolerance = TOLERANCE_COEF)
    expect_equal(unname(gm$sigma_v), exp(cf[["Zv_(Intercept)"]] / 2),
                 tolerance = TOLERANCE_COEF)
    expect_equal(unname(gm$sigma_u), exp(cf[["Zu_(Intercept)"]] / 2),
                 tolerance = TOLERANCE_COEF)
  }
})


test_that("internal exponential SFA on test_data gives finite, plausible estimates", {
  for (g in fit_exp$groups) {
    gm <- fit_exp$group_models[[g]]
    n <- length(gm$y)

    expect_equal(gm$convergence, 0L)
    expect_true(is.finite(gm$sigma_v) && is.finite(gm$sigma_u))
    expect_gt(gm$sigma_v, 0)
    expect_lt(gm$sigma_v, sd(gm$y))
    expect_gt(gm$sigma_u, 0)
    expect_lt(gm$sigma_u, sd(gm$y))

    # The normal regression is the sigma_u -> 0 limit, so the maximum is
    # at least the OLS log-likelihood; no model with N(0, sigma_v^2)
    # noise can exceed n * (-log(sigma_v) - log(2 * pi) / 2).
    ols_ll <- as.numeric(logLik(lm(gm$y ~ gm$X - 1)))
    expect_gte(gm$logLik, ols_ll)
    expect_lt(gm$logLik, n * (-log(gm$sigma_v) - 0.5 * log(2 * pi)))
  }
  expect_true(all(check_convergence(fit_exp)$converged))
})


test_that("the estimator flags an impossible or non-finite logLik with code 99", {
  d2 <- test_data[test_data$group == "G2", ]
  f <- Formula::Formula(log_y ~ log_x1 + log_x2)

  # A log-density of 10 per observation exceeds the bound at any
  # plausible sigma_v; optim() still reports code 0 on the flat surface
  local_mocked_bindings(.ll_exponential = function(eps, sigma_v, sigma_u) {
    rep(10, length(eps))
  })
  expect_warning(
    gm <- metafrontier:::.fit_sfa_group(f, d2, "exponential", list()),
    "failed numerically \\(code 99: logLik .* exceeds bound"
  )
  expect_equal(gm$convergence, 99L)

  local_mocked_bindings(.ll_exponential = function(eps, sigma_v, sigma_u) {
    rep(NaN, length(eps))
  })
  expect_warning(
    gm <- metafrontier:::.fit_sfa_group(f, d2, "exponential", list()),
    "failed numerically \\(code 99: logLik not finite"
  )
  expect_equal(gm$convergence, 99L)
})


test_that("check_convergence flags a stored impossible or non-finite logLik", {
  # A fit as returned before the fix: code 0 at an impossible logLik
  bad <- fit_exp
  bad$group_models$G2$sigma_v <- 4.441459e+91
  bad$group_models$G2$logLik <- 2.050548e+100
  conv <- check_convergence(bad)
  expect_equal(conv$code, c(0L, 99L, 0L))
  expect_equal(conv$converged, c(TRUE, FALSE, TRUE))
  expect_match(conv$note[2], "exceeds bound")
  expect_output(print(conv), "did not converge")

  bad$group_models$G2$logLik <- NaN
  conv <- check_convergence(bad)
  expect_equal(conv$code[2], 99L)
  expect_match(conv$note[2], "not finite")

  # Stochastic metafrontier stage
  fit_st <- metafrontier(log_y ~ log_x1 + log_x2, data = test_data,
                         group = "group", meta_type = "stochastic",
                         dist = "exponential")
  expect_true(all(check_convergence(fit_st)$converged))
  fit_st$meta_logLik <- 1e100
  conv <- check_convergence(fit_st)
  expect_equal(conv$code[3], 99L)
  expect_false(conv$converged[3])
  expect_match(conv$note[3], "exceeds bound")
})
