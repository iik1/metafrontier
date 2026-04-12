#!/usr/bin/env Rscript
# Reviewer 6 Stress Tests for metafrontier v0.2.0
# Petrochemical refinery simulation scenarios
# Tests: collinearity, extreme sigma_u, unbalanced groups, high-dim X,
#        convergence, Murphy-Topel PD, bootstrap failure handling

suppressWarnings(suppressMessages({
  library(metafrontier)
  library(numDeriv)
}))

cat("=== REVIEWER 6 STRESS TESTS ===\n")
cat("Package version:", as.character(packageVersion("metafrontier")), "\n\n")

errors_found <- character(0)
warnings_found <- character(0)

# Helper: generate petrochemical refinery data
gen_petrochem <- function(n1, n2, n_inputs = 2, sigma_u = c(0.3, 0.3),
                          sigma_v = 0.2, collinear = FALSE,
                          collinear_r = 0.99, seed = 123) {
  set.seed(seed)
  n <- n1 + n2
  beta_meta <- c(1.0, seq(0.6, 0.2, length.out = n_inputs))

  make_group <- function(n_g, gap, su) {
    if (collinear && n_inputs >= 2) {
      x1 <- runif(n_g, 0.5, 5)
      noise <- rnorm(n_g, 0, sqrt(1 - collinear_r^2) * sd(x1))
      x2 <- collinear_r * x1 + noise
      if (n_inputs > 2) {
        x_rest <- matrix(runif(n_g * (n_inputs - 2), 0.5, 5),
                         nrow = n_g)
      }
      X_inp <- if (n_inputs > 2) cbind(x1, x2, x_rest) else cbind(x1, x2)
    } else {
      X_inp <- matrix(runif(n_g * n_inputs, 0.5, 5), nrow = n_g)
    }
    colnames(X_inp) <- paste0("log_x", seq_len(n_inputs))
    X <- cbind(1, X_inp)
    beta_g <- beta_meta
    beta_g[1] <- beta_g[1] - gap
    frontier_y <- X %*% beta_g
    v <- rnorm(n_g, 0, sigma_v)
    u <- abs(rnorm(n_g, 0, su))
    log_y <- as.numeric(frontier_y + v - u)
    df <- as.data.frame(X_inp)
    df$log_y <- log_y
    df
  }

  g1 <- make_group(n1, 0.0, sigma_u[1])
  g1$group <- "Refinery_A"
  g2 <- make_group(n2, 0.4, sigma_u[2])
  g2$group <- "Refinery_B"
  rbind(g1, g2)
}

# ============================================================
# TEST 1: Highly collinear inputs (r > 0.95)
# ============================================================
cat("--- TEST 1: Highly collinear inputs (r > 0.95) ---\n")
tryCatch({
  d1 <- gen_petrochem(200, 200, n_inputs = 2, collinear = TRUE,
                      collinear_r = 0.99)
  actual_cor <- cor(d1$log_x1, d1$log_x2)
  cat("  Actual correlation:", round(actual_cor, 4), "\n")

  fit1 <- metafrontier(log_y ~ log_x1 + log_x2, data = d1,
                       group = "group", method = "sfa",
                       meta_type = "stochastic")
  cat("  Convergence (meta):", fit1$meta_convergence, "\n")

  # Check if Hessian is invertible for each group
  for (g in fit1$groups) {
    H <- fit1$group_models[[g]]$hessian
    eig <- eigen(-H, symmetric = TRUE)$values
    cat("  Group", g, "- Hessian eigenvalues:", round(range(eig), 6), "\n")
    if (any(eig <= 0)) {
      msg <- paste("Group", g, ": negative Hessian eigenvalue under collinearity")
      warnings_found <- c(warnings_found, msg)
      cat("  WARNING:", msg, "\n")
    }
  }

  # Check condition number of X'X
  for (g in fit1$groups) {
    X <- fit1$group_models[[g]]$X
    cn <- kappa(crossprod(X))
    cat("  Group", g, "- Condition number X'X:", round(cn, 1), "\n")
    if (cn > 1e6) {
      msg <- paste("Group", g, ": extreme condition number", round(cn))
      warnings_found <- c(warnings_found, msg)
      cat("  WARNING:", msg, "\n")
    }
  }

  cat("  TGR range:", round(range(fit1$tgr), 4), "\n")
  cat("  TEST 1: PASSED (no crash)\n\n")
}, error = function(e) {
  msg <- paste("TEST 1 FAILED:", conditionMessage(e))
  errors_found <<- c(errors_found, msg)
  cat(" ", msg, "\n\n")
})


# ============================================================
# TEST 2: Very small sigma_u (near-zero inefficiency)
# ============================================================
cat("--- TEST 2: Very small sigma_u (near-zero inefficiency) ---\n")
tryCatch({
  d2 <- gen_petrochem(200, 200, sigma_u = c(0.001, 0.001), sigma_v = 0.2)

  fit2 <- metafrontier(log_y ~ log_x1 + log_x2, data = d2,
                       group = "group", method = "sfa",
                       meta_type = "stochastic")
  cat("  Convergence (meta):", fit2$meta_convergence, "\n")

  for (g in fit2$groups) {
    su <- fit2$group_models[[g]]$sigma_u
    sv <- fit2$group_models[[g]]$sigma_v
    lam <- su / sv
    cat("  Group", g, "- sigma_u:", round(su, 6),
        "sigma_v:", round(sv, 6),
        "lambda:", round(lam, 6), "\n")
    if (lam < 0.01) {
      msg <- paste("Group", g, ": lambda near zero, SFA may be unidentified")
      warnings_found <- c(warnings_found, msg)
    }
  }

  # Check efficiency estimates -- should all be near 1
  te <- fit2$te_group
  cat("  TE range:", round(range(te), 6), "\n")
  cat("  TE mean:", round(mean(te), 6), "\n")

  # Check if TGR is sensible
  cat("  TGR range:", round(range(fit2$tgr), 4), "\n")
  cat("  TEST 2: PASSED (no crash)\n\n")
}, error = function(e) {
  msg <- paste("TEST 2 FAILED:", conditionMessage(e))
  errors_found <<- c(errors_found, msg)
  cat(" ", msg, "\n\n")
})


# ============================================================
# TEST 3: Very large sigma_u (extreme inefficiency)
# ============================================================
cat("--- TEST 3: Very large sigma_u (extreme inefficiency) ---\n")
tryCatch({
  d3 <- gen_petrochem(200, 200, sigma_u = c(3.0, 3.0), sigma_v = 0.2)

  fit3 <- metafrontier(log_y ~ log_x1 + log_x2, data = d3,
                       group = "group", method = "sfa",
                       meta_type = "stochastic")
  cat("  Convergence (meta):", fit3$meta_convergence, "\n")

  for (g in fit3$groups) {
    su <- fit3$group_models[[g]]$sigma_u
    sv <- fit3$group_models[[g]]$sigma_v
    cat("  Group", g, "- sigma_u:", round(su, 4),
        "sigma_v:", round(sv, 4), "\n")
  }

  te <- fit3$te_group
  cat("  TE range:", round(range(te), 6), "\n")
  cat("  TE mean:", round(mean(te), 6), "\n")

  # With large sigma_u, some TE should be very small
  cat("  Fraction TE < 0.1:", round(mean(te < 0.1), 4), "\n")

  # Check for any NaN or Inf in TGR
  if (any(!is.finite(fit3$tgr))) {
    msg <- "TEST 3: Non-finite TGR values with large sigma_u"
    errors_found <- c(errors_found, msg)
    cat("  ERROR:", msg, "\n")
  }
  cat("  TGR range:", round(range(fit3$tgr), 4), "\n")
  cat("  TEST 3: PASSED (no crash)\n\n")
}, error = function(e) {
  msg <- paste("TEST 3 FAILED:", conditionMessage(e))
  errors_found <<- c(errors_found, msg)
  cat(" ", msg, "\n\n")
})


# ============================================================
# TEST 4: Highly unbalanced groups (30 vs 500)
# ============================================================
cat("--- TEST 4: Unbalanced groups (30 vs 500) ---\n")
tryCatch({
  d4 <- gen_petrochem(30, 500)

  fit4 <- metafrontier(log_y ~ log_x1 + log_x2, data = d4,
                       group = "group", method = "sfa",
                       meta_type = "stochastic")
  cat("  Convergence (meta):", fit4$meta_convergence, "\n")

  for (g in fit4$groups) {
    n_g <- fit4$group_models[[g]]$nobs
    su <- fit4$group_models[[g]]$sigma_u
    sv <- fit4$group_models[[g]]$sigma_v
    cat("  Group", g, "(n=", n_g, ") - sigma_u:", round(su, 4),
        "sigma_v:", round(sv, 4), "\n")
  }

  cat("  TGR range:", round(range(fit4$tgr), 4), "\n")

  # Also check deterministic metafrontier with unbalanced groups
  fit4d <- metafrontier(log_y ~ log_x1 + log_x2, data = d4,
                        group = "group", method = "sfa",
                        meta_type = "deterministic")
  cat("  Deterministic - convergence:", fit4d$meta_convergence, "\n")
  cat("  Deterministic TGR range:", round(range(fit4d$tgr), 4), "\n")

  cat("  TEST 4: PASSED (no crash)\n\n")
}, error = function(e) {
  msg <- paste("TEST 4 FAILED:", conditionMessage(e))
  errors_found <<- c(errors_found, msg)
  cat(" ", msg, "\n\n")
})


# ============================================================
# TEST 5: High-dimensional X (6 inputs)
# ============================================================
cat("--- TEST 5: High-dimensional X (6 inputs) ---\n")
tryCatch({
  d5 <- gen_petrochem(200, 200, n_inputs = 6)

  f5 <- log_y ~ log_x1 + log_x2 + log_x3 + log_x4 + log_x5 + log_x6
  fit5 <- metafrontier(f5, data = d5,
                       group = "group", method = "sfa",
                       meta_type = "stochastic")
  cat("  Convergence (meta):", fit5$meta_convergence, "\n")

  for (g in fit5$groups) {
    conv <- fit5$group_models[[g]]$convergence
    cat("  Group", g, "- convergence:", conv, "\n")
    if (conv != 0) {
      msg <- paste("Group", g, ": did not converge with 6 inputs")
      warnings_found <- c(warnings_found, msg)
    }
  }

  cat("  Meta coef:", round(fit5$meta_coef, 4), "\n")
  cat("  TGR range:", round(range(fit5$tgr), 4), "\n")
  cat("  TEST 5: PASSED (no crash)\n\n")
}, error = function(e) {
  msg <- paste("TEST 5 FAILED:", conditionMessage(e))
  errors_found <<- c(errors_found, msg)
  cat(" ", msg, "\n\n")
})


# ============================================================
# TEST 6: Murphy-Topel correction PD check
# ============================================================
cat("--- TEST 6: Murphy-Topel correction positive-definiteness ---\n")
tryCatch({
  d6 <- gen_petrochem(200, 200)
  fit6 <- metafrontier(log_y ~ log_x1 + log_x2, data = d6,
                       group = "group", method = "sfa",
                       meta_type = "stochastic")

  V_mt <- vcov(fit6, correction = "murphy-topel")

  if (is.null(V_mt)) {
    msg <- "Murphy-Topel returned NULL"
    errors_found <- c(errors_found, msg)
    cat("  ERROR:", msg, "\n")
  } else {
    eig_mt <- eigen(V_mt, symmetric = TRUE)$values
    cat("  MT eigenvalues:", round(eig_mt, 8), "\n")
    if (any(eig_mt <= 0)) {
      msg <- paste("Murphy-Topel matrix NOT positive definite.",
                   "Min eigenvalue:", min(eig_mt))
      errors_found <- c(errors_found, msg)
      cat("  ERROR:", msg, "\n")
    } else {
      cat("  Murphy-Topel matrix is positive definite.\n")
    }

    # Compare uncorrected vs corrected SEs
    V_unc <- vcov(fit6, correction = "none")
    se_unc <- sqrt(diag(V_unc))
    se_mt <- sqrt(diag(V_mt))
    cat("  Uncorrected SEs:", round(se_unc, 6), "\n")
    cat("  Murphy-Topel SEs:", round(se_mt, 6), "\n")

    # MT SEs should generally be >= uncorrected
    if (any(se_mt < se_unc * 0.9)) {
      msg <- "Murphy-Topel SEs are SMALLER than uncorrected for some params"
      warnings_found <- c(warnings_found, msg)
      cat("  WARNING:", msg, "\n")
    }
  }
  cat("  TEST 6: PASSED (no crash)\n\n")
}, error = function(e) {
  msg <- paste("TEST 6 FAILED:", conditionMessage(e))
  errors_found <<- c(errors_found, msg)
  cat(" ", msg, "\n\n")
})


# ============================================================
# TEST 7: Murphy-Topel with collinear data
# ============================================================
cat("--- TEST 7: Murphy-Topel under collinearity ---\n")
tryCatch({
  d7 <- gen_petrochem(200, 200, collinear = TRUE, collinear_r = 0.98)
  fit7 <- metafrontier(log_y ~ log_x1 + log_x2, data = d7,
                       group = "group", method = "sfa",
                       meta_type = "stochastic")

  V_mt7 <- vcov(fit7, correction = "murphy-topel")
  if (is.null(V_mt7)) {
    msg <- "Murphy-Topel returned NULL under collinearity"
    errors_found <- c(errors_found, msg)
    cat("  ERROR:", msg, "\n")
  } else {
    eig_mt7 <- eigen(V_mt7, symmetric = TRUE)$values
    cat("  MT eigenvalues:", round(eig_mt7, 8), "\n")
    if (any(eig_mt7 <= 0)) {
      msg <- paste("Murphy-Topel NOT PD under collinearity.",
                   "Min eigenvalue:", min(eig_mt7))
      errors_found <- c(errors_found, msg)
      cat("  ERROR:", msg, "\n")
    }
  }
  cat("  TEST 7: PASSED (no crash)\n\n")
}, error = function(e) {
  msg <- paste("TEST 7 FAILED:", conditionMessage(e))
  errors_found <<- c(errors_found, msg)
  cat(" ", msg, "\n\n")
})


# ============================================================
# TEST 8: Bootstrap with small R, check failure handling
# ============================================================
cat("--- TEST 8: Bootstrap failure handling ---\n")
tryCatch({
  d8 <- gen_petrochem(50, 50, sigma_u = c(0.01, 0.01))
  fit8 <- metafrontier(log_y ~ log_x1 + log_x2, data = d8,
                       group = "group", method = "sfa",
                       meta_type = "stochastic")

  # Run bootstrap with very few reps
  boot8 <- boot_tgr(fit8, R = 20, type = "nonparametric",
                    seed = 42, progress = FALSE)
  cat("  Bootstrap R_effective:", boot8$R_effective, "/", boot8$R, "\n")
  cat("  Any NA in CI:", any(is.na(boot8$ci)), "\n")

  # Check CI bounds
  if (!is.null(boot8$ci)) {
    cat("  CI range [1,]:", round(boot8$ci[1,], 4), "\n")
    if (any(boot8$ci[,1] > boot8$ci[,2], na.rm = TRUE)) {
      msg <- "Bootstrap CI lower > upper for some observations"
      errors_found <- c(errors_found, msg)
      cat("  ERROR:", msg, "\n")
    }
  }
  cat("  TEST 8: PASSED (no crash)\n\n")
}, error = function(e) {
  msg <- paste("TEST 8 FAILED:", conditionMessage(e))
  errors_found <<- c(errors_found, msg)
  cat(" ", msg, "\n\n")
})


# ============================================================
# TEST 9: Parametric bootstrap with extreme sigma_u
# ============================================================
cat("--- TEST 9: Parametric bootstrap with extreme sigma_u ---\n")
tryCatch({
  d9 <- gen_petrochem(100, 100, sigma_u = c(2.0, 2.0))
  fit9 <- metafrontier(log_y ~ log_x1 + log_x2, data = d9,
                       group = "group", method = "sfa",
                       meta_type = "stochastic")

  boot9 <- boot_tgr(fit9, R = 15, type = "parametric",
                    seed = 99, progress = FALSE)
  cat("  Bootstrap R_effective:", boot9$R_effective, "/", boot9$R, "\n")
  fail_rate <- 1 - boot9$R_effective / boot9$R
  cat("  Failure rate:", round(fail_rate * 100, 1), "%\n")

  if (fail_rate > 0.5) {
    msg <- paste("Parametric bootstrap >50% failure rate with large sigma_u:",
                 round(fail_rate * 100, 1), "%")
    errors_found <- c(errors_found, msg)
    cat("  ERROR:", msg, "\n")
  }
  cat("  TEST 9: PASSED (no crash)\n\n")
}, error = function(e) {
  msg <- paste("TEST 9 FAILED:", conditionMessage(e))
  errors_found <<- c(errors_found, msg)
  cat(" ", msg, "\n\n")
})


# ============================================================
# TEST 10: Input validation - missing checks
# ============================================================
cat("--- TEST 10: Input validation ---\n")

# 10a: Single group should error
tryCatch({
  d10a <- gen_petrochem(100, 0)
  d10a <- d10a[d10a$group == "Refinery_A", ]
  fit10a <- metafrontier(log_y ~ log_x1 + log_x2, data = d10a,
                         group = "group", method = "sfa")
  msg <- "Single-group data did NOT raise an error"
  errors_found <- c(errors_found, msg)
  cat("  10a ERROR:", msg, "\n")
}, error = function(e) {
  cat("  10a: Correctly errors on single group:", conditionMessage(e), "\n")
})

# 10b: NA in response
tryCatch({
  d10b <- gen_petrochem(100, 100)
  d10b$log_y[1:5] <- NA
  fit10b <- metafrontier(log_y ~ log_x1 + log_x2, data = d10b,
                         group = "group", method = "sfa",
                         meta_type = "stochastic")
  # Check if NAs were handled
  cat("  10b: NAs in response - nobs:", fit10b$nobs["total"], "\n")
  if (fit10b$nobs["total"] == 200) {
    cat("  10b: Total obs = 200 despite 5 NAs -- NAs silently dropped in groups but nobs counts all rows\n")
    msg <- "NA handling: nobs[total] counts original rows, not fitted rows"
    warnings_found <- c(warnings_found, msg)
  }
}, error = function(e) {
  cat("  10b: NA handling error:", conditionMessage(e), "\n")
  msg <- paste("NA in response crashes:", conditionMessage(e))
  errors_found <- c(errors_found, msg)
})

# 10c: Negative sigma_v via control override
tryCatch({
  d10c <- gen_petrochem(100, 100)
  fit10c <- metafrontier(log_y ~ log_x1 + log_x2, data = d10c,
                         group = "group", method = "sfa",
                         meta_type = "stochastic",
                         control = list(maxit = 1))
  cat("  10c: maxit=1 convergence:", fit10c$meta_convergence, "\n")
  for (g in fit10c$groups) {
    cat("    Group", g, "convergence:", fit10c$group_models[[g]]$convergence, "\n")
  }
}, error = function(e) {
  cat("  10c ERROR:", conditionMessage(e), "\n")
})

cat("  TEST 10: Completed\n\n")


# ============================================================
# TEST 11: TGR > 1 check for stochastic metafrontier
# ============================================================
cat("--- TEST 11: TGR bounds check ---\n")
tryCatch({
  d11 <- gen_petrochem(200, 200, sigma_u = c(0.5, 0.5))
  fit11_det <- metafrontier(log_y ~ log_x1 + log_x2, data = d11,
                            group = "group", method = "sfa",
                            meta_type = "deterministic")
  fit11_sto <- metafrontier(log_y ~ log_x1 + log_x2, data = d11,
                            group = "group", method = "sfa",
                            meta_type = "stochastic")

  cat("  Deterministic TGR range:", round(range(fit11_det$tgr), 4), "\n")
  cat("  Stochastic TGR range:", round(range(fit11_sto$tgr), 4), "\n")

  # Deterministic should be bounded [0,1]
  if (any(fit11_det$tgr > 1.001)) {
    msg <- paste("Deterministic TGR > 1:", max(fit11_det$tgr))
    errors_found <- c(errors_found, msg)
    cat("  ERROR:", msg, "\n")
  }

  # Stochastic TGR can exceed 1 -- is this documented/expected?
  if (any(fit11_sto$tgr > 1.0)) {
    n_above <- sum(fit11_sto$tgr > 1.0)
    max_above <- max(fit11_sto$tgr)
    msg <- paste("Stochastic TGR > 1 for", n_above, "obs (max:", round(max_above, 4), ")")
    warnings_found <- c(warnings_found, msg)
    cat("  NOTE:", msg, "\n")
  }

  # Check TE* = TE x TGR identity
  te_star_manual <- fit11_sto$te_group * fit11_sto$tgr
  max_diff <- max(abs(te_star_manual - fit11_sto$te_meta))
  cat("  TE* decomposition max error:", max_diff, "\n")
  if (max_diff > 1e-10) {
    msg <- paste("TE* = TE x TGR decomposition fails, max error:", max_diff)
    errors_found <- c(errors_found, msg)
    cat("  ERROR:", msg, "\n")
  }

  cat("  TEST 11: PASSED\n\n")
}, error = function(e) {
  msg <- paste("TEST 11 FAILED:", conditionMessage(e))
  errors_found <<- c(errors_found, msg)
  cat(" ", msg, "\n\n")
})


# ============================================================
# TEST 12: Convergence with 5 groups
# ============================================================
cat("--- TEST 12: 5 groups convergence ---\n")
tryCatch({
  sim12 <- simulate_metafrontier(n_groups = 5, n_per_group = 80,
                                 sigma_u = c(0.2, 0.3, 0.4, 0.5, 0.6),
                                 seed = 77)
  fit12 <- metafrontier(log_y ~ log_x1 + log_x2,
                        data = sim12$data, group = "group",
                        method = "sfa", meta_type = "stochastic")
  cat("  Convergence (meta):", fit12$meta_convergence, "\n")
  for (g in fit12$groups) {
    cat("  Group", g, "- conv:", fit12$group_models[[g]]$convergence,
        "sigma_u:", round(fit12$group_models[[g]]$sigma_u, 4), "\n")
  }
  cat("  TGR range:", round(range(fit12$tgr), 4), "\n")
  cat("  TEST 12: PASSED\n\n")
}, error = function(e) {
  msg <- paste("TEST 12 FAILED:", conditionMessage(e))
  errors_found <<- c(errors_found, msg)
  cat(" ", msg, "\n\n")
})


# ============================================================
# TEST 13: Check OLS starting values with extreme data
# ============================================================
cat("--- TEST 13: OLS starting values stability ---\n")
tryCatch({
  set.seed(456)
  n <- 200
  d13 <- data.frame(
    log_x1 = runif(n, 0, 10),
    log_x2 = runif(n, 0, 10),
    group = rep(c("A", "B"), each = n/2)
  )
  # Create very noisy data with outliers
  d13$log_y <- 1 + 0.5 * d13$log_x1 + 0.3 * d13$log_x2 +
    rnorm(n, 0, 5) - abs(rnorm(n, 0, 0.3))
  # Add extreme outliers
  d13$log_y[1:3] <- c(100, -50, 200)

  fit13 <- metafrontier(log_y ~ log_x1 + log_x2, data = d13,
                        group = "group", method = "sfa",
                        meta_type = "stochastic")
  cat("  Convergence:", fit13$meta_convergence, "\n")
  cat("  Meta coef:", round(fit13$meta_coef, 4), "\n")
  cat("  TEST 13: PASSED (no crash)\n\n")
}, error = function(e) {
  msg <- paste("TEST 13 FAILED:", conditionMessage(e))
  errors_found <<- c(errors_found, msg)
  cat(" ", msg, "\n\n")
})


# ============================================================
# TEST 14: Poolability test edge cases
# ============================================================
cat("--- TEST 14: Poolability test ---\n")
tryCatch({
  # Groups with identical DGP (should NOT reject H0)
  sim14 <- simulate_metafrontier(n_groups = 2, n_per_group = 200,
                                 tech_gap = c(0, 0), sigma_u = c(0.3, 0.3),
                                 seed = 88)
  fit14 <- metafrontier(log_y ~ log_x1 + log_x2,
                        data = sim14$data, group = "group",
                        method = "sfa", meta_type = "deterministic")
  pt14 <- poolability_test(fit14)
  cat("  Identical DGP - LR stat:", round(pt14$statistic, 4),
      "p-value:", round(pt14$p.value, 4), "\n")

  # If LR statistic is negative, that's a problem
  if (pt14$statistic < 0) {
    msg <- paste("Poolability test LR statistic is negative:", pt14$statistic)
    errors_found <- c(errors_found, msg)
    cat("  ERROR:", msg, "\n")
  }

  cat("  TEST 14: PASSED\n\n")
}, error = function(e) {
  msg <- paste("TEST 14 FAILED:", conditionMessage(e))
  errors_found <<- c(errors_found, msg)
  cat(" ", msg, "\n\n")
})


# ============================================================
# SUMMARY
# ============================================================
cat("\n========= SUMMARY =========\n")
cat("Errors found:", length(errors_found), "\n")
for (e in errors_found) cat("  [ERROR] ", e, "\n")
cat("\nWarnings found:", length(warnings_found), "\n")
for (w in warnings_found) cat("  [WARN]  ", w, "\n")
cat("============================\n")
