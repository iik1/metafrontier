#' Internal: Estimate SFA-based metafrontier
#'
#' Implements the deterministic metafrontier of Battese, Rao, and
#' O'Donnell (2004) via linear programming and the stochastic
#' metafrontier of Huang, Huang, and Liu (2014) via second-stage SFA.
#'
#' @keywords internal
#' @noRd
.estimate_sfa_metafrontier <- function(formula, data, group_vec,
                                       group_levels, group_models,
                                       meta_type, dist, control) {

  n <- nrow(data)
  n_groups <- length(group_levels)

  # Collect group frontier predictions at each observation's input mix
  if (inherits(formula, "Formula")) {
    f_base <- formula(formula, rhs = 1)
  } else {
    f_base <- formula
  }
  mf <- model.frame(f_base, data = data, na.action = na.omit)
  X <- model.matrix(f_base, data = mf)
  k <- ncol(X)

  # Group frontier values: ln f(x_i; beta_j) for each obs i in group j
  group_frontier <- numeric(n)
  te_group <- numeric(n)

  for (g in group_levels) {
    idx <- which(group_vec == g)
    beta_g <- group_models[[g]]$coefficients
    group_frontier[idx] <- X[idx, , drop = FALSE] %*% beta_g
    te_group[idx] <- group_models[[g]]$efficiency
  }

  # Collect group coefficients
  group_coef <- lapply(group_models, function(m) m$coefficients)

  if (meta_type == "deterministic") {
    meta_result <- .deterministic_metafrontier_lp(
      X, group_frontier, group_vec, group_levels, group_coef, k
    )
  } else {
    meta_result <- .stochastic_metafrontier(
      X, group_frontier, group_vec, group_levels, dist, control
    )
  }

  # Compute TGR and metafrontier efficiency
  meta_frontier <- as.numeric(X %*% meta_result$meta_coef)
  tgr <- exp(group_frontier - meta_frontier)

  # Bound TGR to (0, 1] for deterministic metafrontier
  if (meta_type == "deterministic") {
    tgr <- pmin(tgr, 1.0)
  }

  te_meta <- te_group * tgr

  list(
    meta_coef = meta_result$meta_coef,
    meta_vcov = meta_result$meta_vcov,
    group_coef = group_coef,
    tgr = tgr,
    te_group = te_group,
    te_meta = te_meta,
    group_frontier = group_frontier,
    meta_frontier = meta_frontier,
    logLik_groups = sapply(group_models, function(m) m$logLik),
    meta_logLik = meta_result$meta_logLik,
    meta_convergence = meta_result$convergence
  )
}


#' Deterministic metafrontier via LP (Battese, Rao, O'Donnell 2004)
#'
#' Minimises sum of squared deviations of the metafrontier from
#' group frontiers, subject to the envelopment constraint that
#' the metafrontier weakly dominates all group frontiers at all
#' observed input mixes.
#'
#' @keywords internal
#' @noRd
.deterministic_metafrontier_lp <- function(X, group_frontier,
                                           group_vec, group_levels,
                                           group_coef, k) {

  n <- length(group_frontier)

  # The LP approach: for log-linear models, the metafrontier

  # coefficients beta* must satisfy:
  #   x_i' beta* >= x_i' beta_j  for all i in group j
  #
  # We minimise: sum_i (x_i' beta* - x_i' beta_j)^2
  #
  # This is a quadratic programming problem. We solve it using

  # an iterative approach or reformulate as LP.
  #
  # Standard LP reformulation:
  #   min sum_i d_i
  #   s.t. x_i' beta* - x_i' beta_j_hat <= d_i  for all i
  #        x_i' beta* - x_i' beta_j_hat >= 0     for all i
  #        d_i >= 0

  # For simplicity and numerical stability, use QP via optim with
  # linear inequality constraints via constrOptim

  # Objective: minimise sum_i (X_i' beta* - gf_i)^2
  # Subject to: X_i' beta* >= gf_i for all i (envelopment)

  # Use constrOptim (built-in R, no extra dependency needed)
  # Constraint: -X %*% beta* <= -gf  (i.e., X %*% beta* >= gf)

  obj_fn <- function(beta_star) {
    fitted_meta <- X %*% beta_star
    sum((fitted_meta - group_frontier)^2)
  }

  grad_fn <- function(beta_star) {
    fitted_meta <- X %*% beta_star
    2 * crossprod(X, fitted_meta - group_frontier)
  }

  # Starting value: OLS on group frontier values (unconstrained)
  beta_start <- lm.fit(X, group_frontier)$coefficients

  # Ensure starting point is feasible
  violations <- group_frontier - X %*% beta_start
  if (any(violations > 0)) {
    # Shift intercept up to make feasible
    max_violation <- max(violations)
    beta_start[1] <- beta_start[1] + max_violation + 0.01
  }

  # Constraints: -X %*% beta <= -gf
  ui <- X        # X %*% beta >= gf
  ci <- group_frontier  # lower bounds

  result <- constrOptim(
    theta = beta_start,
    f = obj_fn,
    grad = grad_fn,
    ui = ui,
    ci = ci - 1e-6,  # small slack for numerical stability
    method = "BFGS",
    control = list(maxit = 10000, reltol = 1e-12)
  )

  meta_coef <- result$par
  names(meta_coef) <- colnames(X)

  list(
    meta_coef = meta_coef,
    meta_vcov = NULL,   # no variance for deterministic metafrontier
    meta_logLik = NULL,
    convergence = result$convergence
  )
}


#' Stochastic metafrontier (Huang, Huang, Liu 2014)
#'
#' Estimates the metafrontier as a second-stage SFA where the
#' estimated group frontiers are treated as dependent variables.
#'
#' @keywords internal
#' @noRd
.stochastic_metafrontier <- function(X, group_frontier,
                                     group_vec, group_levels,
                                     dist, control) {

  n <- length(group_frontier)
  k <- ncol(X)

  # The stochastic metafrontier:
  #   ln f_hat(x_i; beta_j) = x_i' beta* + v*_i - u*_i
  #
  # where v* is noise and u* >= 0 is the technology gap
  # TGR_i = exp(-u*_i)

  # Fit SFA directly to the group frontier values
  # (X and group_frontier are already extracted, so we use the
  # log-likelihood functions directly rather than .fit_sfa_group)

  ols <- lm.fit(X, group_frontier)
  ols_resid <- ols$residuals
  sigma_ols <- sqrt(sum(ols_resid^2) / (n - k))

  start_params <- c(ols$coefficients,
                    log_sigma_v = log(sigma_ols * 0.5),
                    log_sigma_u = log(sigma_ols * 0.5))

  ctrl <- list(fnscale = -1, maxit = 5000, reltol = 1e-10)
  ctrl[names(control)] <- control

  opt <- optim(
    par = start_params,
    fn = .loglik_hnormal,
    y = group_frontier,
    X = X,
    method = "BFGS",
    control = ctrl,
    hessian = TRUE
  )

  meta_coef <- opt$par[seq_len(k)]
  names(meta_coef) <- colnames(X)
  sigma_v <- exp(opt$par["log_sigma_v"])
  sigma_u <- exp(opt$par["log_sigma_u"])

  # Variance-covariance from inverse of negative Hessian
  meta_vcov <- tryCatch(
    solve(-opt$hessian),
    error = function(e) {
      warning("Hessian is singular; variance-covariance matrix unavailable.",
              call. = FALSE)
      NULL
    }
  )

  list(
    meta_coef = meta_coef,
    meta_vcov = meta_vcov,
    meta_sigma_v = sigma_v,
    meta_sigma_u = sigma_u,
    meta_logLik = opt$value,
    convergence = opt$convergence
  )
}
