#' Internal: Panel SFA estimation for a single group
#'
#' Fits a panel SFA model to a single group using MLE.
#' Implements BC92 (Battese and Coelli, 1992) with time-varying
#' inefficiency: u_it = u_i * exp(-eta*(t - T)) where T is the global
#' final period and u_i ~ |N(0, sigma_u^2)|.
#' Also implements BC95 (Battese and Coelli, 1995) with observation-specific
#' mean: u_it ~ N+(z_it'delta, sigma_u^2).
#'
#' @param formula a Formula object.
#' @param data data frame for this group, must contain panel id/time cols.
#' @param dist distribution of u (currently "hnormal" for BC92).
#' @param panel_dist panel model type: "bc92" or "bc95".
#' @param panel_info list with id and time column names.
#' @param control list of control parameters.
#' @param ... additional arguments.
#' @param estimator character. Technical efficiency estimator:
#'   \code{"bc88"} (Battese and Coelli, 1988; \code{E[exp(-u)|eps]})
#'   or \code{"jlms"} (Jondrow et al., 1982; \code{exp(-E[u|eps])}).
#'
#' @return A list compatible with .fit_sfa_group output plus panel fields.
#' @keywords internal
#' @noRd
# Safe Mills ratio to avoid division by zero
.safe_mills <- function(z) {
  p <- pnorm(z)
  p <- pmax(p, .Machine$double.eps)
  dnorm(z) / p
}

.fit_sfa_panel_group <- function(formula, data, dist, panel_dist,
                                 panel_info, control, ...,
                                 estimator = c("bc88", "jlms")) {

  estimator <- match.arg(estimator)

  id_col <- panel_info$id
  time_col <- panel_info$time

  if (!id_col %in% names(data) || !time_col %in% names(data)) {
    stop("Panel columns '", id_col, "' and/or '", time_col,
         "' not found in data.", call. = FALSE)
  }

  # Build model matrices (input row order is preserved throughout so
  # that returned vectors align with the caller's data)
  if (inherits(formula, "Formula")) {
    f <- formula
    has_z <- length(f)[2] >= 2L

    if (has_z) {
      mf <- model.frame(f, data = data, na.action = na.omit)
      y <- model.response(mf)
      X <- model.matrix(f, data = mf, rhs = 1)
      Z <- model.matrix(f, data = mf, rhs = 2)
      if (ncol(Z) == 0L) { Z <- NULL; has_z <- FALSE }
    } else {
      f_base <- formula(f, rhs = 1)
      mf <- model.frame(f_base, data = data, na.action = na.omit)
      y <- model.response(mf)
      X <- model.matrix(f_base, data = mf)
      Z <- NULL
    }
  } else {
    mf <- model.frame(formula, data = data, na.action = na.omit)
    y <- model.response(mf)
    X <- model.matrix(formula, data = mf)
    Z <- NULL
    has_z <- FALSE
  }

  n <- length(y)
  k <- ncol(X)

  # Panel structure; drop any rows removed from the model frame by
  # na.omit so that firms/times stay aligned with y and X
  firms <- data[[id_col]]
  times <- data[[time_col]]
  na_act <- attr(mf, "na.action")
  if (!is.null(na_act)) {
    firms <- firms[-na_act]
    times <- times[-na_act]
  }
  firm_ids <- unique(firms)
  n_firms <- length(firm_ids)

  # Build firm-level index: list of obs indices per firm
  firm_idx <- lapply(firm_ids, function(f) which(firms == f))
  T_i <- sapply(firm_idx, length)  # periods per firm
  T_max <- max(times)

  # Time index relative to the global final period T for BC92:
  # u_it = u_i * exp(-eta * (t - T)), matching the unbalanced-panel
  # formulation of Battese and Coelli (1992)
  t_rel <- times - T_max

  # OLS starting values
  ols <- lm.fit(X, y)
  ols_resid <- ols$residuals
  sigma_ols <- sqrt(sum(ols_resid^2) / (n - k))

  sigma_v_start <- max(sigma_ols * 0.5, 0.1)
  sigma_u_start <- max(sigma_ols * 0.5, 0.1)

  if (panel_dist == "bc92") {
    # BC92: params = [beta, log_sigma_v, log_sigma_u, eta]
    start_params <- c(ols$coefficients,
                      log_sigma_v = log(sigma_v_start),
                      log_sigma_u = log(sigma_u_start),
                      eta = 0.0)

    loglik_fn <- function(params) {
      .loglik_bc92(params, y, X, k, firms, firm_idx, T_i, t_rel)
    }

  } else if (panel_dist == "bc95") {
    # BC95: requires Z variables for mean function
    if (!has_z) {
      stop("BC95 model requires inefficiency determinants (Z variables) ",
           "in the formula (e.g., y ~ x1 + x2 | z1 + z2).", call. = FALSE)
    }
    p <- ncol(Z)
    delta_start <- rep(0, p)
    delta_start[1] <- 0

    start_params <- c(ols$coefficients,
                      log_sigma_v = log(sigma_v_start),
                      setNames(delta_start, paste0("d_", seq_len(p))),
                      log_sigma_u = log(sigma_u_start))

    loglik_fn <- function(params) {
      .loglik_bc95(params, y, X, Z, k)
    }
  } else {
    stop("Unknown panel_dist: ", panel_dist, call. = FALSE)
  }

  ctrl <- list(fnscale = -1, maxit = 5000, reltol = 1e-10)
  ctrl[names(control)] <- control

  opt <- tryCatch(
    optim(par = start_params, fn = loglik_fn,
          method = "BFGS", control = ctrl, hessian = TRUE),
    error = function(e) {
      tryCatch(
        optim(par = start_params, fn = loglik_fn,
              method = "Nelder-Mead",
              control = list(fnscale = -1, maxit = 10000),
              hessian = TRUE),
        error = function(e2) {
          stop("MLE optimisation failed. The data may have too few ",
               "observations or extreme values. Original error: ",
               conditionMessage(e), call. = FALSE)
        }
      )
    }
  )

  if (opt$convergence != 0) {
    warning("Panel SFA (", panel_dist, ") did not converge (code ",
            opt$convergence, ").", call. = FALSE)
  }

  # Extract parameters
  beta <- opt$par[seq_len(k)]
  names(beta) <- colnames(X)
  sigma_v <- exp(opt$par["log_sigma_v"])
  sigma_u <- exp(opt$par["log_sigma_u"])

  fitted_vals <- as.numeric(X %*% beta)

  if (panel_dist == "bc92") {
    eta <- opt$par["eta"]

    # Conditional efficiency for BC92, both JLMS-style and the
    # Battese-Coelli (1992) closed form
    eps <- y - fitted_vals
    te_list <- .bc92_efficiency(eps, sigma_v, sigma_u, eta,
                                firms, firm_idx, T_i, t_rel)
    te_jlms <- te_list$jlms
    te_bc88 <- te_list$bc88
    te <- if (estimator == "bc88") te_bc88 else te_jlms

    result <- list(
      coefficients = beta,
      sigma_v = as.numeric(sigma_v),
      sigma_u = as.numeric(sigma_u),
      eta = as.numeric(eta),
      logLik = opt$value,
      efficiency = te,
      efficiency_jlms = te_jlms,
      efficiency_bc88 = te_bc88,
      estimator = estimator,
      fitted = fitted_vals,
      residuals = eps,
      hessian = opt$hessian,
      convergence = opt$convergence,
      all_params = opt$par,
      y = y,
      X = X,
      dist = "hnormal",
      panel_dist = panel_dist,
      panel_info = panel_info,
      firms = firms,
      n_firms = n_firms,
      T_i = T_i
    )

  } else {
    # BC95: cross-sectional truncated-normal with Z
    # Efficiency via JLMS conditional mean
    p <- ncol(Z)
    delta <- opt$par[(k + 1 + 1):(k + 1 + p)]

    eps <- y - fitted_vals
    mu_i <- as.numeric(Z %*% delta)
    sigma_sq <- sigma_v^2 + sigma_u^2
    lambda <- sigma_u / sigma_v
    mu_star <- (mu_i * sigma_v^2 - eps * sigma_u^2) / sigma_sq
    sigma_star <- sigma_v * sigma_u / sqrt(sigma_sq)

    u_hat <- mu_star + sigma_star * .safe_mills(mu_star / sigma_star)
    te_jlms <- as.numeric(exp(-u_hat))

    # BC88 (Battese and Coelli, 1988): E[exp(-u)|eps]. Where
    # Phi(mu*/sigma*) underflows to zero the ratio is indeterminate;
    # those observations fall back to the JLMS value.
    bc_ratio <- mu_star / sigma_star
    bc_denom <- pnorm(bc_ratio)
    te_bc88 <- as.numeric(
      exp(-mu_star + 0.5 * sigma_star^2) *
        pnorm(bc_ratio - sigma_star) / bc_denom
    )
    bc_bad <- !is.finite(te_bc88) | bc_denom == 0
    te_bc88[bc_bad] <- te_jlms[bc_bad]

    te <- if (estimator == "bc88") te_bc88 else te_jlms

    result <- list(
      coefficients = beta,
      sigma_v = as.numeric(sigma_v),
      sigma_u = as.numeric(sigma_u),
      delta = delta,
      logLik = opt$value,
      efficiency = te,
      efficiency_jlms = te_jlms,
      efficiency_bc88 = te_bc88,
      estimator = estimator,
      fitted = fitted_vals,
      residuals = eps,
      hessian = opt$hessian,
      convergence = opt$convergence,
      all_params = opt$par,
      y = y,
      X = X,
      Z = Z,
      dist = "tnormal",
      panel_dist = panel_dist,
      panel_info = panel_info,
      firms = firms,
      n_firms = n_firms,
      T_i = T_i
    )
  }

  result
}


# ---------- BC92 log-likelihood (half-normal, time-varying) ----------

#' BC92 log-likelihood
#'
#' u_it = u_i * exp(-eta*(t - T)) with T the global final period,
#' u_i ~ |N(0, sigma_u^2)|. Integrated over u_i (closed-form for
#' half-normal).
#'
#' @keywords internal
#' @noRd
# Per-firm marginal log-density after integrating out u_i
# (Battese and Coelli, 1992, Eq. 8):
# ll_i = -T_i/2 log(2*pi) - T_i log sigma_v - sum_t eps_it^2/(2 sigma_v^2)
#        + log 2 - 0.5 log(1 + sigma_u^2 sum_t d_t^2 / sigma_v^2)
#        + 0.5 (mu_i*/sigma_i*)^2 + log Phi(mu_i*/sigma_i*).
# Parameter layout: params = c(beta[1:k], log_sigma_v, log_sigma_u, eta).
.loglik_bc92 <- function(params, y, X, k, firms, firm_idx, T_i, t_rel) {

  beta <- params[seq_len(k)]
  sigma_v <- exp(params["log_sigma_v"])
  sigma_u <- exp(params["log_sigma_u"])
  eta <- params["eta"]

  if (sigma_v < 1e-10 || sigma_u < 1e-10) return(-1e20)

  eps <- as.numeric(y - X %*% beta)
  n_firms <- length(firm_idx)
  ll <- 0

  for (ff in seq_along(firm_idx)) {
    idx <- firm_idx[[ff]]
    T_f <- length(idx)
    eps_f <- eps[idx]
    t_rel_f <- t_rel[idx]

    # d_t = exp(-eta * (t - T))
    d_t <- exp(-eta * t_rel_f)

    # Integrated LL for firm f:
    # sum_t [ -0.5*log(2*pi*sigma_v^2) - eps_ft^2/(2*sigma_v^2) ]
    # + log(2) + log(sigma_u) - log(sigma_star_f)
    # + 0.5*(mu_star_f/sigma_star_f)^2
    # + log(Phi(mu_star_f/sigma_star_f))
    # - log(sigma_u) - 0.5*log(1 + sigma_u^2 * sum(d_t^2) / sigma_v^2)

    sum_d2 <- sum(d_t^2)
    sum_eps_d <- sum(eps_f * d_t)

    sigma_star2 <- 1 / (1 / sigma_u^2 + sum_d2 / sigma_v^2)
    sigma_star <- sqrt(sigma_star2)
    mu_star <- -sigma_star2 * sum_eps_d / sigma_v^2

    ll_f <- -T_f * 0.5 * log(2 * pi) -
      T_f * log(sigma_v) -
      sum(eps_f^2) / (2 * sigma_v^2) +
      log(2) -
      0.5 * log(1 + sigma_u^2 * sum_d2 / sigma_v^2) +
      0.5 * (mu_star / sigma_star)^2 +
      pnorm(mu_star / sigma_star, log.p = TRUE)

    ll <- ll + ll_f
  }

  if (!is.finite(ll)) return(-1e20)
  ll
}


# ---------- BC92 efficiency estimation ----------

# Returns both the JLMS-style estimator exp(-E[u_i|eps] * d_t) and the
# Battese-Coelli (1992, Eq. 10) closed form
# TE_it = {Phi(mu_i*/sigma_i* - d_t*sigma_i*) / Phi(mu_i*/sigma_i*)}
#         * exp(-d_t*mu_i* + 0.5*d_t^2*sigma_i*^2),
# where d_t = exp(-eta*(t - T)) and (mu_i*, sigma_i*) are the per-firm
# conditional posterior parameters of u_i.
.bc92_efficiency <- function(eps, sigma_v, sigma_u, eta,
                             firms, firm_idx, T_i, t_rel) {
  n <- length(eps)
  te_jlms <- numeric(n)
  te_bc88 <- numeric(n)

  for (ff in seq_along(firm_idx)) {
    idx <- firm_idx[[ff]]
    eps_f <- eps[idx]
    t_rel_f <- t_rel[idx]

    d_t <- exp(-eta * t_rel_f)
    sum_d2 <- sum(d_t^2)
    sum_eps_d <- sum(eps_f * d_t)

    sigma_star2 <- 1 / (1 / sigma_u^2 + sum_d2 / sigma_v^2)
    sigma_star <- sqrt(sigma_star2)
    mu_star <- -sigma_star2 * sum_eps_d / sigma_v^2

    # E[u_i | eps] (conditional mean of firm effect)
    ratio <- mu_star / sigma_star
    E_ui <- mu_star + sigma_star * .safe_mills(ratio)

    # u_it = u_i * exp(-eta*(t-T))
    u_it <- E_ui * d_t
    te_jlms[idx] <- exp(-u_it)

    # BC92 closed form; where Phi(mu*/sigma*) underflows to zero the
    # ratio is indeterminate, so fall back to the JLMS value
    denom <- pnorm(ratio)
    te_f <- pnorm(ratio - d_t * sigma_star) / denom *
      exp(-d_t * mu_star + 0.5 * d_t^2 * sigma_star2)
    bad <- !is.finite(te_f) | denom == 0
    te_f[bad] <- te_jlms[idx][bad]
    te_bc88[idx] <- te_f
  }

  list(jlms = te_jlms, bc88 = te_bc88)
}


# ---------- BC95 log-likelihood (truncated-normal with Z) ----------

#' BC95 log-likelihood
#'
#' Same as cross-sectional truncated-normal with Z-variables,
#' since BC95 treats each observation individually with
#' u_it ~ N+(z_it'delta, sigma_u^2).
#'
#' @keywords internal
#' @noRd
# Observation-level log-density (Battese and Coelli, 1995):
# ll_it = log phi((eps_it + mu_it)/sigma) - log sigma
#         + log Phi(mu_it*/sigma*) - log Phi(mu_it/sigma_u),
# with mu_it = z_it'delta. Parameter layout:
# params = c(beta[1:k], log_sigma_v, delta[1:p], log_sigma_u).
.loglik_bc95 <- function(params, y, X, Z, k) {
  n <- length(y)
  p <- ncol(Z)

  beta <- params[seq_len(k)]
  sigma_v <- exp(params[k + 1])  # log_sigma_v
  delta <- params[(k + 2):(k + 1 + p)]
  sigma_u <- exp(params[k + 1 + p + 1])  # log_sigma_u

  if (sigma_v < 1e-10 || sigma_u < 1e-10) return(-1e20)

  sigma_sq <- sigma_v^2 + sigma_u^2
  sigma <- sqrt(sigma_sq)

  eps <- as.numeric(y - X %*% beta)
  mu_i <- as.numeric(Z %*% delta)

  mu_star <- (mu_i * sigma_v^2 - eps * sigma_u^2) / sigma_sq
  sigma_star <- sigma_v * sigma_u / sigma

  ll <- -0.5 * log(2 * pi) - log(sigma) -
    0.5 * ((eps + mu_i)^2 / sigma_sq) +
    pnorm(mu_star / sigma_star, log.p = TRUE) -
    pnorm(mu_i / sigma_u, log.p = TRUE)

  sum(ll)
}
