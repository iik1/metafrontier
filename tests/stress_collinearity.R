#!/usr/bin/env Rscript
# Stress test: collinear and ill-conditioned inputs for metafrontier
# Run: Rscript tests/stress_collinearity.R

suppressPackageStartupMessages(library(metafrontier))

cat("===== METAFRONTIER COLLINEARITY STRESS TESTS =====\n\n")

results <- list()
pass <- 0L
fail <- 0L
warn_count <- 0L

run_test <- function(name, expr) {
  cat(sprintf("--- TEST: %s ---\n", name))
  res <- tryCatch(
    withCallingHandlers(
      {
        val <- eval(expr)
        list(status = "PASS", value = val, warnings = character(0))
      },
      warning = function(w) {
        cat(sprintf("  WARNING: %s\n", conditionMessage(w)))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      list(status = "ERROR", message = conditionMessage(e))
    }
  )
  cat(sprintf("  Result: %s\n", res$status))
  if (res$status == "ERROR") cat(sprintf("  Message: %s\n", res$message))
  cat("\n")
  res
}

# Helper: build data with 2 groups, custom X columns
make_data <- function(n_per_group, make_x_fn, n_groups = 2L) {
  frames <- list()
  for (g in seq_len(n_groups)) {
    X <- make_x_fn(n_per_group)
    k <- ncol(X)
    beta <- c(1.0, seq(0.5, 0.2, length.out = k))
    log_y <- cbind(1, X) %*% beta + rnorm(n_per_group, 0, 0.2) -
      abs(rnorm(n_per_group, 0, 0.3))
    df <- as.data.frame(X)
    names(df) <- paste0("log_x", seq_len(k))
    df$log_y <- as.numeric(log_y)
    df$group <- paste0("G", g)
    frames[[g]] <- df
  }
  do.call(rbind, frames)
}

# ===================================================================
# TEST 1: Perfect collinearity (log_x2 = 2 * log_x1)
# ===================================================================
set.seed(42)
run_test("Perfect collinearity (x2 = 2*x1)", quote({
  dat <- make_data(100, function(n) {
    x1 <- runif(n, 0, 5)
    cbind(x1, 2 * x1)  # perfect collinearity
  })
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = dat, group = "group",
                      method = "sfa", meta_type = "deterministic")
  cat(sprintf("  Meta coefs: %s\n",
              paste(round(fit$meta_coef, 4), collapse = ", ")))
  cat(sprintf("  Mean TGR: %.4f\n", mean(fit$tgr)))
  cat(sprintf("  TE range: [%.4f, %.4f]\n",
              min(fit$te_meta), max(fit$te_meta)))
  "completed"
}))

# ===================================================================
# TEST 2: Near-perfect collinearity (x2 = x1 + tiny noise)
# ===================================================================
set.seed(42)
run_test("Near-perfect collinearity (x2 = x1 + eps)", quote({
  dat <- make_data(100, function(n) {
    x1 <- runif(n, 0, 5)
    cbind(x1, x1 + rnorm(n, 0, 0.001))
  })
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = dat, group = "group",
                      method = "sfa", meta_type = "deterministic")
  cat(sprintf("  Meta coefs: %s\n",
              paste(round(fit$meta_coef, 4), collapse = ", ")))

  # Check SEs from stochastic metafrontier
  fit2 <- metafrontier(log_y ~ log_x1 + log_x2,
                       data = dat, group = "group",
                       method = "sfa", meta_type = "stochastic")
  if (!is.null(fit2$meta_vcov)) {
    ses <- sqrt(diag(fit2$meta_vcov))
    cat(sprintf("  Stochastic SEs: %s\n",
                paste(round(ses, 4), collapse = ", ")))
    cat(sprintf("  Any SE > 100? %s\n", any(ses > 100)))
  } else {
    cat("  Stochastic vcov: NULL (Hessian singular)\n")
  }
  "completed"
}))

# ===================================================================
# TEST 3: Many redundant inputs (10 cols, 8 are linear combos)
# ===================================================================
set.seed(42)
run_test("8 redundant inputs (10 total, 8 linear combos of 2)", quote({
  dat <- make_data(100, function(n) {
    x1 <- runif(n, 0, 5)
    x2 <- runif(n, 0, 5)
    X <- cbind(x1, x2,
               x1 + x2, x1 - x2, 2*x1 + x2, x1 + 2*x2,
               3*x1, 3*x2, x1 + 3*x2, 2*x1 + 3*x2)
    X
  })
  fml <- as.formula(paste("log_y ~",
                          paste0("log_x", 1:10, collapse = " + ")))
  fit <- metafrontier(fml, data = dat, group = "group",
                      method = "sfa", meta_type = "deterministic")
  cat(sprintf("  Meta coefs: %s\n",
              paste(round(fit$meta_coef, 4), collapse = ", ")))
  cat(sprintf("  Any Inf/NaN coef? %s\n",
              any(!is.finite(fit$meta_coef))))
  "completed"
}))

# ===================================================================
# TEST 4: Constant column (log_x3 = 1.0 for all)
# ===================================================================
set.seed(42)
run_test("Constant column (log_x3 = 1.0)", quote({
  dat <- make_data(100, function(n) {
    cbind(runif(n, 0, 5), runif(n, 0, 5), rep(1.0, n))
  })
  fit <- metafrontier(log_y ~ log_x1 + log_x2 + log_x3,
                      data = dat, group = "group",
                      method = "sfa", meta_type = "deterministic")
  cat(sprintf("  Meta coefs: %s\n",
              paste(round(fit$meta_coef, 4), collapse = ", ")))
  "completed"
}))

# ===================================================================
# TEST 5: Single input (log_y ~ log_x1 only)
# ===================================================================
set.seed(42)
run_test("Single input (log_y ~ log_x1)", quote({
  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 100,
                               n_inputs = 1, seed = 42)
  fit <- metafrontier(log_y ~ log_x1,
                      data = sim$data, group = "group",
                      method = "sfa", meta_type = "deterministic")
  cat(sprintf("  Meta coefs: %s\n",
              paste(round(fit$meta_coef, 4), collapse = ", ")))
  cat(sprintf("  Mean TGR: %.4f, Mean TE*: %.4f\n",
              mean(fit$tgr), mean(fit$te_meta)))
  "completed"
}))

# ===================================================================
# TEST 6: Wide data (p > n within group: 20 inputs, 10 obs/group)
# ===================================================================
set.seed(42)
run_test("Wide data (20 inputs, 10 obs per group)", quote({
  dat <- make_data(10, function(n) {
    matrix(runif(n * 20, 0, 5), nrow = n, ncol = 20)
  })
  fml <- as.formula(paste("log_y ~",
                          paste0("log_x", 1:20, collapse = " + ")))
  fit <- metafrontier(fml, data = dat, group = "group",
                      method = "sfa", meta_type = "deterministic")
  cat(sprintf("  Any Inf/NaN coef? %s\n",
              any(!is.finite(fit$meta_coef))))
  cat(sprintf("  Convergence: %d\n", fit$convergence))
  "completed"
}))

# ===================================================================
# TEST 7: High condition number (~1e12)
# ===================================================================
set.seed(42)
run_test("High condition number X (~1e12)", quote({
  dat <- make_data(100, function(n) {
    x1 <- runif(n, 0, 5)
    x2 <- x1 * 1e6 + rnorm(n, 0, 1e-6)  # condition number ~1e12
    cbind(x1, x2)
  })
  cond_num <- kappa(cbind(1, dat$log_x1, dat$log_x2))
  cat(sprintf("  Condition number: %.2e\n", cond_num))
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = dat, group = "group",
                      method = "sfa", meta_type = "stochastic")
  cat(sprintf("  Meta coefs: %s\n",
              paste(round(fit$meta_coef, 4), collapse = ", ")))
  if (!is.null(fit$meta_vcov)) {
    ses <- sqrt(diag(fit$meta_vcov))
    cat(sprintf("  SEs: %s\n", paste(round(ses, 6), collapse = ", ")))
    cat(sprintf("  Any SE > 1000? %s\n", any(ses > 1000)))
  } else {
    cat("  Vcov: NULL (singular Hessian)\n")
  }
  "completed"
}))

# ===================================================================
# TEST 8: DEA with perfect collinearity
# ===================================================================
set.seed(42)
run_test("DEA: Perfect collinearity (x2 = 2*x1)", quote({
  dat <- make_data(50, function(n) {
    x1 <- runif(n, 0, 5)
    cbind(x1, 2 * x1)
  })
  # DEA needs positive levels, not logs with possible negatives
  dat$log_x1 <- exp(dat$log_x1)
  dat$log_x2 <- exp(dat$log_x2)
  dat$log_y <- exp(dat$log_y)
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = dat, group = "group",
                      method = "dea", meta_type = "deterministic",
                      rts = "vrs")
  cat(sprintf("  Mean TGR: %.4f\n", mean(fit$tgr)))
  cat(sprintf("  TE range: [%.4f, %.4f]\n",
              min(fit$te_meta), max(fit$te_meta)))
  "completed"
}))

# ===================================================================
# TEST 9: DEA with wide data (p > n)
# ===================================================================
set.seed(42)
run_test("DEA: Wide data (20 inputs, 10 obs per group)", quote({
  dat <- make_data(10, function(n) {
    matrix(runif(n * 20, 1, 5), nrow = n, ncol = 20)
  })
  fml <- as.formula(paste("log_y ~",
                          paste0("log_x", 1:20, collapse = " + ")))
  dat$log_y <- exp(dat$log_y)
  for (j in 1:20) dat[[paste0("log_x", j)]] <- exp(dat[[paste0("log_x", j)]])
  fit <- metafrontier(fml, data = dat, group = "group",
                      method = "dea", meta_type = "deterministic",
                      rts = "vrs")
  cat(sprintf("  Mean TGR: %.4f\n", mean(fit$tgr)))
  cat(sprintf("  All TE = 1? %s\n",
              all(abs(fit$te_group - 1) < 1e-6)))
  "completed"
}))

# ===================================================================
# TEST 10: Stochastic metafrontier with perfect collinearity
# ===================================================================
set.seed(42)
run_test("Stochastic metafrontier: Perfect collinearity", quote({
  dat <- make_data(100, function(n) {
    x1 <- runif(n, 0, 5)
    cbind(x1, 2 * x1)
  })
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = dat, group = "group",
                      method = "sfa", meta_type = "stochastic")
  cat(sprintf("  Meta coefs: %s\n",
              paste(round(fit$meta_coef, 4), collapse = ", ")))
  cat(sprintf("  Convergence: %d\n", fit$convergence))
  if (!is.null(fit$meta_vcov)) {
    ses <- sqrt(diag(fit$meta_vcov))
    cat(sprintf("  SEs finite? %s\n", all(is.finite(ses))))
  } else {
    cat("  Vcov: NULL\n")
  }
  "completed"
}))

cat("\n===== ALL TESTS COMPLETE =====\n")
