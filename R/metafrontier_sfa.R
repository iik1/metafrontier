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

  # Determine valid (non-NA) rows after na.omit
  na_action <- attr(mf, "na.action")
  if (!is.null(na_action)) {
    valid_rows <- seq_len(n)[-na_action]
  } else {
    valid_rows <- seq_len(n)
  }
  n_valid <- length(valid_rows)

  # Filter group_vec to only valid rows
  group_vec_valid <- group_vec[valid_rows]

  # Group frontier values: ln f(x_i; beta_j) for each valid obs i in group j
  group_frontier <- numeric(n_valid)
  te_group <- numeric(n_valid)

  for (g in group_levels) {
    idx <- which(group_vec_valid == g)
    beta_g <- group_models[[g]]$coefficients
    group_frontier[idx] <- X[idx, , drop = FALSE] %*% beta_g
    te_group[idx] <- group_models[[g]]$efficiency
  }

  # Collect group coefficients
  group_coef <- lapply(group_models, function(m) m$coefficients)

  if (meta_type == "deterministic") {
    meta_result <- .deterministic_metafrontier_lp(
      X, group_frontier, group_vec_valid, group_levels, group_coef, k
    )
  } else {
    meta_result <- .stochastic_metafrontier(
      X, group_frontier, group_vec_valid, group_levels, dist, control,
      group_models
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
    meta_convergence = meta_result$convergence,
    meta_opt = meta_result$meta_opt,
    meta_dist = meta_result$meta_dist,
    valid_rows = valid_rows,
    n_valid = n_valid
  )
}


#' Deterministic metafrontier via LP (Battese, Rao, O'Donnell 2004)
#'
#' Minimises the sum of deviations of the metafrontier from
#' group frontiers via linear programming, subject to the
#' envelopment constraint that the metafrontier weakly dominates
#' all group frontiers at all observed input mixes.
#'
#' LP formulation (BRO 2004):
#'   min  sum_i (x_i' beta* - gf_i)   [= colSums(X)' beta* - sum(gf)]
#'   s.t. x_i' beta* >= gf_i          for all i (envelopment)
#'
#' Falls back to QP via constrOptim if the LP solver fails.
#'
#' @keywords internal
#' @noRd
.deterministic_metafrontier_lp <- function(X, group_frontier,
                                           group_vec, group_levels,
                                           group_coef, k) {

  n <- length(group_frontier)

  # ----- Primary: LP via lpSolveAPI (BRO 2004 formulation) -----
  # Variables: beta* (k free variables)
  # Objective: min colSums(X)' beta* (constant -sum(gf) dropped)
  # Constraints: x_i' beta* >= gf_i for all i

  lp_result <- tryCatch({
    lp <- lpSolveAPI::make.lp(nrow = n, ncol = k)

    # Set columns (variables beta*)
    for (j in seq_len(k)) {
      lpSolveAPI::set.column(lp, j, X[, j])
    }

    # Objective: min sum_i x_i' beta* = colSums(X)' beta*
    lpSolveAPI::set.objfn(lp, colSums(X))
    lpSolveAPI::lp.control(lp, sense = "min")

    # Constraints: x_i' beta* >= gf_i
    for (i in seq_len(n)) {
      lpSolveAPI::set.constr.type(lp, type = ">=", constraints = i)
      lpSolveAPI::set.rhs(lp, b = group_frontier[i], constraints = i)
    }

    # beta* variables are free (unbounded)
    for (j in seq_len(k)) {
      lpSolveAPI::set.bounds(lp, lower = -1e30, upper = 1e30, columns = j)
    }

    status <- lpSolveAPI::solve.lpExtPtr(lp)
    if (status != 0L) stop("lpSolveAPI returned status ", status)

    meta_coef <- lpSolveAPI::get.variables(lp)
    names(meta_coef) <- colnames(X)

    list(
      meta_coef = meta_coef,
      meta_vcov = NULL,
      meta_logLik = NULL,
      convergence = 0L
    )
  }, error = function(e) {
    warning("LP solver failed (", conditionMessage(e),
            "); falling back to QP via constrOptim.", call. = FALSE)
    NULL
  })

  if (!is.null(lp_result)) return(lp_result)

  # ----- Fallback: QP via constrOptim -----
  obj_fn <- function(beta_star) {
    fitted_meta <- X %*% beta_star
    sum((fitted_meta - group_frontier)^2)
  }

  grad_fn <- function(beta_star) {
    fitted_meta <- X %*% beta_star
    2 * crossprod(X, fitted_meta - group_frontier)
  }

  beta_start <- lm.fit(X, group_frontier)$coefficients

  if (any(!is.finite(beta_start))) {
    avg_coef <- rowMeans(do.call(cbind, group_coef))
    if (all(is.finite(avg_coef))) {
      beta_start <- avg_coef
    } else {
      beta_start <- rep(0, k)
      beta_start[1] <- max(group_frontier) + 0.1
    }
  }

  violations <- group_frontier - X %*% beta_start
  if (any(violations > 0)) {
    beta_start[1] <- beta_start[1] + max(violations) + 0.01
  }

  if (!is.finite(obj_fn(beta_start))) {
    beta_start <- rep(0, k)
    beta_start[1] <- max(group_frontier) + 0.1
    violations <- group_frontier - X %*% beta_start
    if (any(violations > 0)) {
      beta_start[1] <- beta_start[1] + max(violations) + 0.01
    }
  }

  ui <- X
  ci <- group_frontier

  result <- constrOptim(
    theta = beta_start,
    f = obj_fn,
    grad = grad_fn,
    ui = ui,
    ci = ci - 1e-6,
    method = "BFGS",
    control = list(maxit = 10000, reltol = 1e-12)
  )

  meta_coef <- result$par
  names(meta_coef) <- colnames(X)

  list(
    meta_coef = meta_coef,
    meta_vcov = NULL,
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
                                     dist, control,
                                     group_models = NULL) {

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

  # Select log-likelihood function based on dist
  loglik_fn <- switch(dist,
    hnormal     = .loglik_hnormal,
    tnormal     = .loglik_tnormal,
    exponential = .loglik_exponential
  )

  ols <- lm.fit(X, group_frontier)
  ols_resid <- ols$residuals
  sigma_ols <- sqrt(sum(ols_resid^2) / (n - k))

  # Starting values with lower bound guard
  log_sv <- log(max(sigma_ols * 0.5, 0.05))
  log_su <- log(max(sigma_ols * 0.5, 0.05))

  if (dist == "tnormal") {
    start_params <- c(ols$coefficients,
                      log_sigma_v = log_sv,
                      mu = 0,
                      log_sigma_u = log_su)
  } else {
    start_params <- c(ols$coefficients,
                      log_sigma_v = log_sv,
                      log_sigma_u = log_su)
  }

  ctrl <- list(fnscale = -1, maxit = 5000, reltol = 1e-10)
  ctrl[names(control)] <- control

  # Optimise with BFGS, falling back to Nelder-Mead on failure
  opt <- tryCatch(
    optim(par = start_params, fn = loglik_fn, y = group_frontier,
          X = X, method = "BFGS", control = ctrl, hessian = TRUE),
    error = function(e) {
      tryCatch(
        optim(par = start_params, fn = loglik_fn, y = group_frontier,
              X = X, method = "Nelder-Mead",
              control = list(fnscale = -1, maxit = 10000),
              hessian = TRUE),
        error = function(e2) {
          stop("Stage 2 MLE optimisation failed. Original error: ",
               conditionMessage(e), call. = FALSE)
        }
      )
    }
  )

  # Check convergence
  if (opt$convergence != 0L) {
    warning("Stage 2 SFA optimisation did not converge (code ",
            opt$convergence, "). Results may be unreliable.",
            call. = FALSE)
  }

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
    convergence = opt$convergence,
    meta_opt = opt,
    meta_dist = dist,
    stage1_models = group_models
  )
}


#' Murphy-Topel variance correction for two-stage estimation
#'
#' Adjusts the Stage 2 variance-covariance matrix to account for
#' Stage 1 estimation uncertainty (the generated-regressor problem).
#'
#' @param object A fitted metafrontier object (stochastic).
#' @return A corrected variance-covariance matrix.
#' @references Murphy, K.M. and Topel, R.H. (1985). Estimation and
#'   inference in two-step econometric models. \emph{Journal of
#'   Business & Economic Statistics}, 3(4), 370--379.
#' @keywords internal
#' @noRd
.murphy_topel_correction <- function(object) {

  group_models <- object$group_models
  X <- model.matrix(
    formula(object$formula, rhs = 1),
    data = model.frame(formula(object$formula, rhs = 1),
                       data = object$data, na.action = na.omit)
  )
  group_vec <- object$group_vec
  group_levels <- object$groups
  group_frontier <- object$group_frontier
  meta_opt <- object$meta_opt

  n <- nrow(X)
  k <- ncol(X)

  # Select log-likelihood functions based on Stage 2 dist
  loglik_fn <- switch(object$meta_dist,
    hnormal     = .loglik_hnormal,
    tnormal     = .loglik_tnormal,
    exponential = .loglik_exponential
  )
  obs_ll_fn <- switch(object$meta_dist,
    hnormal     = .loglik_hnormal_obs,
    tnormal     = .loglik_tnormal_obs,
    exponential = .loglik_exponential_obs
  )

  # Stage 2 vcov (uncorrected)
  V2 <- tryCatch(solve(-meta_opt$hessian), error = function(e) NULL)
  if (is.null(V2)) return(NULL)

  # Stage 2 log-likelihood as a function of group_frontier
  stage2_ll <- function(gf) {
    loglik_fn(meta_opt$par, y = gf, X = X)
  }

  # For each group j, compute:
  #   R_j = dL2/d(beta_j) via chain rule through group_frontier
  #   V1_j = solve(-H1_j)
  #   S1_j = score matrix (n_j x p_j)

  # Accumulate the correction terms
  # C = sum_j R_j' V1_j R_j  (simplified Murphy-Topel)
  # Plus cross-term: sum_j R_j' V1_j S1_j' S2

  # Stage 2 score vector (n x p2)
  S2_obs <- obs_ll_fn(meta_opt$par, y = group_frontier, X = X)
  S2 <- numDeriv::jacobian(
    function(p) obs_ll_fn(p, y = group_frontier, X = X),
    meta_opt$par,
    method.args = list(eps = 1e-4)
  )
  # S2 is n x p2

  p2 <- length(meta_opt$par)
  C_mat <- matrix(0, p2, p2)
  A_mat <- matrix(0, p2, p2)

  for (g in group_levels) {
    idx <- which(group_vec == g)
    gm <- group_models[[g]]
    n_g <- length(idx)
    p1_g <- length(gm$all_params)
    k_g <- length(gm$coefficients)

    # Stage 1 vcov for this group
    V1_g <- tryCatch(solve(-gm$hessian), error = function(e) NULL)
    if (is.null(V1_g)) next

    # Stage 1 score matrix (n_g x p1_g)
    S1_g <- .score_vector_sfa(gm)

    # R_j: derivative of Stage 2 LL w.r.t. Stage 1 beta_j
    # group_frontier[idx] = X[idx,] %*% beta_j
    # So d(group_frontier)/d(beta_j) = X[idx,] (for the beta portion)
    # and 0 for sigma_v, sigma_u params of Stage 1
    # Use numerical Jacobian for robustness:
    R_g <- numDeriv::jacobian(
      function(theta1) {
        gf_mod <- group_frontier
        beta_g_new <- theta1[seq_len(k_g)]
        gf_mod[idx] <- X[idx, , drop = FALSE] %*% beta_g_new
        loglik_fn(meta_opt$par, y = gf_mod, X = X)
      },
      gm$all_params,
      method.args = list(eps = 1e-4)
    )
    # R_g is 1 x p1_g -> need gradient vector per Stage 2 param
    # Actually we need the p2 x p1_g cross-derivative matrix.
    # Use Jacobian of the Stage 2 score w.r.t. theta1:
    R_g_full <- numDeriv::jacobian(
      function(theta1) {
        gf_mod <- group_frontier
        beta_g_new <- theta1[seq_len(k_g)]
        gf_mod[idx] <- X[idx, , drop = FALSE] %*% beta_g_new
        numDeriv::grad(
          function(p) loglik_fn(p, y = gf_mod, X = X),
          meta_opt$par,
          method.args = list(eps = 1e-4)
        )
      },
      gm$all_params,
      method.args = list(eps = 1e-4)
    )
    # R_g_full is p2 x p1_g

    # C term: R' V1 R
    C_mat <- C_mat + R_g_full %*% V1_g %*% t(R_g_full)

    # A term: R' (S1' S2_g) where S2_g is the Stage 2 scores
    # for observations in group g
    S2_g <- S2[idx, , drop = FALSE]
    A_mat <- A_mat + R_g_full %*% V1_g %*% (t(S1_g) %*% S2_g)
  }

  # Murphy-Topel corrected variance:
  # V_MT = V2 + V2 (C - A - A') V2
  correction <- C_mat - A_mat - t(A_mat)
  V_MT <- V2 + V2 %*% correction %*% V2

  # Ensure symmetry
  V_MT <- (V_MT + t(V_MT)) / 2

  # PSD enforcement: project to nearest PSD matrix if needed
  eig <- eigen(V_MT, symmetric = TRUE)
  if (any(eig$values < 0)) {
    eig$values <- pmax(eig$values, 0)
    V_MT <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
    V_MT <- (V_MT + t(V_MT)) / 2
    warning("Murphy-Topel corrected matrix was not positive semi-definite. ",
            "Projected to nearest PSD matrix.", call. = FALSE)
  }

  V_MT
}
