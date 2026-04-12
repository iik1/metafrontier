# Tests for Murphy-Topel variance correction (P1)

test_that("Murphy-Topel corrected SE >= uncorrected SE", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data, group = "group",
                      meta_type = "stochastic")

  v_uncorr <- vcov(fit, correction = "none")
  v_corr   <- vcov(fit, correction = "murphy-topel")

  expect_true(is.matrix(v_corr))
  expect_equal(nrow(v_corr), nrow(v_uncorr))
  expect_equal(ncol(v_corr), ncol(v_uncorr))

  se_uncorr <- sqrt(diag(v_uncorr))
  se_corr   <- sqrt(diag(v_corr))

  # Corrected SEs should be at least as large as uncorrected
 for (i in seq_along(se_corr)) {
    expect_gte(se_corr[i], se_uncorr[i] * 0.99,
               label = paste("SE for", names(se_corr)[i]))
  }
})

test_that("Murphy-Topel vcov is symmetric and PSD", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data, group = "group",
                      meta_type = "stochastic")

  v_corr <- vcov(fit, correction = "murphy-topel")
  expect_true(is.matrix(v_corr))

  # Symmetric
  expect_equal(v_corr, t(v_corr), tolerance = 1e-10)

  # Positive semi-definite (all eigenvalues >= 0)
  eig <- eigen(v_corr, symmetric = TRUE)$values
  expect_true(all(eig >= -1e-8),
              info = paste("Min eigenvalue:", min(eig)))
})

test_that("confint passes correction through", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data, group = "group",
                      meta_type = "stochastic")

  ci_uncorr <- confint(fit)
  ci_corr   <- confint(fit, correction = "murphy-topel")

  expect_equal(nrow(ci_uncorr), nrow(ci_corr))

  # Corrected intervals should be at least as wide
  width_uncorr <- ci_uncorr[, 2] - ci_uncorr[, 1]
  width_corr   <- ci_corr[, 2] - ci_corr[, 1]

  for (i in seq_along(width_corr)) {
    expect_gte(width_corr[i], width_uncorr[i] * 0.99,
               label = paste("CI width for", rownames(ci_corr)[i]))
  }
})

test_that("Murphy-Topel errors on deterministic metafrontier", {
  fit_det <- metafrontier(log_y ~ log_x1 + log_x2,
                          data = test_data, group = "group",
                          meta_type = "deterministic")

  expect_error(
    vcov(fit_det, correction = "murphy-topel"),
    "Murphy-Topel correction requires a stochastic metafrontier"
  )
})

test_that("Murphy-Topel errors on DEA metafrontier", {
  fit_dea <- metafrontier(log_y ~ log_x1 + log_x2,
                          data = test_data, group = "group",
                          method = "dea")

  expect_error(
    vcov(fit_dea, correction = "murphy-topel"),
    "not available for DEA"
  )
})

test_that("vcov correction='none' matches default", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data, group = "group",
                      meta_type = "stochastic")

  v_default <- vcov(fit)
  v_none    <- vcov(fit, correction = "none")

  expect_equal(v_default, v_none)
})

test_that("observation-level LL sums to total LL", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data, group = "group",
                      meta_type = "stochastic")

  for (g in fit$groups) {
    gm <- fit$group_models[[g]]
    params <- gm$all_params
    y <- gm$y
    X <- gm$X

    # Total LL from the group model
    ll_total <- gm$logLik

    # Sum of obs-level LL
    ll_obs <- metafrontier:::.loglik_hnormal_obs(params, y, X)
    expect_equal(sum(ll_obs), ll_total, tolerance = 1e-6,
                 label = paste("LL sum for group", g))
  }
})

test_that("score vector has correct dimensions", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data, group = "group",
                      meta_type = "stochastic")

  for (g in fit$groups) {
    gm <- fit$group_models[[g]]
    S <- metafrontier:::.score_vector_sfa(gm)

    expect_equal(nrow(S), length(gm$y),
                 label = paste("Score rows for group", g))
    expect_equal(ncol(S), length(gm$all_params),
                 label = paste("Score cols for group", g))

    # At the MLE, score columns should sum to ~0
    col_sums <- colSums(S)
    for (j in seq_along(col_sums)) {
      expect_lt(abs(col_sums[j]), 1.0,
                label = paste("Score sum col", j, "group", g))
    }
  }
})
