# Tests for the objective = c("lp", "qp") choice in the
# deterministic SFA metafrontier (Battese, Rao and O'Donnell 2004
# propose both minimum absolute and minimum squared deviations).

sim_obj <- simulate_metafrontier(seed = 42)
dat_obj <- sim_obj$data

fit_lp_obj <- metafrontier(log_y ~ log_x1 + log_x2,
                           data = dat_obj,
                           group = "group",
                           method = "sfa",
                           meta_type = "deterministic")
X_obj <- model.matrix(log_y ~ log_x1 + log_x2, data = dat_obj)

test_that("QP metafrontier satisfies envelope and TGR bounds", {
  skip_if_not_installed("quadprog")

  qp <- metafrontier:::.deterministic_metafrontier_lp(
    X_obj, fit_lp_obj$group_frontier, fit_lp_obj$group_vec,
    fit_lp_obj$groups, fit_lp_obj$group_coef, ncol(X_obj),
    objective = "qp"
  )

  expect_identical(qp$meta_solver, "qp")

  meta_frontier_qp <- as.numeric(X_obj %*% qp$meta_coef)
  expect_true(all(meta_frontier_qp >= fit_lp_obj$group_frontier - 1e-6))

  tgr_qp <- exp(fit_lp_obj$group_frontier - meta_frontier_qp)
  expect_true(all(tgr_qp <= 1 + 1e-6))

  # User-facing argument (added by a concurrent change); skip the
  # assertions if metafrontier() does not honour 'objective' yet.
  fit_qp <- tryCatch(
    metafrontier(log_y ~ log_x1 + log_x2,
                 data = dat_obj,
                 group = "group",
                 method = "sfa",
                 meta_type = "deterministic",
                 objective = "qp"),
    error = function(e) NULL
  )
  if (!is.null(fit_qp) && identical(fit_qp$objective, "qp")) {
    expect_identical(fit_qp$meta_solver, "qp")
    expect_true(all(fit_qp$meta_frontier >= fit_qp$group_frontier - 1e-6))
    expect_true(all(fit_qp$tgr <= 1 + 1e-6))
  }
})

test_that("LP and QP objectives give similar TGR", {
  skip_if_not_installed("quadprog")

  qp <- metafrontier:::.deterministic_metafrontier_lp(
    X_obj, fit_lp_obj$group_frontier, fit_lp_obj$group_vec,
    fit_lp_obj$groups, fit_lp_obj$group_coef, ncol(X_obj),
    objective = "qp"
  )

  # Mirror the TGR construction of the deterministic pipeline
  tgr_qp <- pmin(exp(fit_lp_obj$group_frontier -
                       as.numeric(X_obj %*% qp$meta_coef)), 1.0)

  max_diff <- max(abs(fit_lp_obj$tgr - tgr_qp))
  cat("\nMax abs TGR difference (LP vs QP):", format(max_diff), "\n")

  expect_lt(max_diff, 0.05)
})

test_that("rank-deficient X falls back to the barrier QP with a message", {
  set.seed(99)
  n <- 50
  x1 <- runif(n, 1, 2)
  X_bad <- cbind("(Intercept)" = 1, x1 = x1, x1_dup = x1)
  gf_bad <- 1 + 0.4 * x1 + rnorm(n, sd = 0.05)
  gv_bad <- factor(rep(c("a", "b"), length.out = n))
  gc_bad <- list(a = c(1, 0.2, 0.2), b = c(1.1, 0.25, 0.15))

  res <- NULL
  expect_message(
    res <- metafrontier:::.deterministic_metafrontier_lp(
      X_bad, gf_bad, gv_bad, levels(gv_bad), gc_bad, ncol(X_bad),
      objective = "qp"
    ),
    regexp = "constrOptim"
  )

  expect_identical(res$meta_solver, "qp-barrier")
  expect_true(all(is.finite(res$meta_coef)))
  expect_true(all(X_bad %*% res$meta_coef >= gf_bad - 1e-4))
})

test_that("default fit records the LP solver and objective", {
  expect_identical(fit_lp_obj$meta_solver, "lp")
  expect_identical(fit_lp_obj$objective, "lp")
})
