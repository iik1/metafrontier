#' Internal: Panel SFA estimation for a single group
#'
#' Fits a panel SFA model to a single group using MLE.
#' Implements BC92 (Battese and Coelli, 1992) with time-varying
#' inefficiency: u_it = u_i * exp(-eta*(t-T_i)) where u_i ~ |N(0, sigma_u^2)|.
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
                                 panel_info, control, ...) {

  id_col <- panel_info$id
  time_col <- panel_info$time

  if (!id_col %in% names(data) || !time_col %in% names(data)) {
    stop("Panel columns '", id_col, "' and/or '", time_col,
         "' not found in data.", call. = FALSE)
  }

  # Sort by (firm, time) for clean panel structure
  data <- data[order(data[[id_col]], data[[time_col]]), ]

  # Build model matrices
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

  # Panel structure
  firms <- data[[id_col]]
  times <- data[[time_col]]
  firm_ids <- unique(firms)
  n_firms <- length(firm_ids)

  # Build firm-level index: list of obs indices per firm
  firm_idx <- lapply(firm_ids, function(f) which(firms == f))
  T_i <- sapply(firm_idx, length)  # periods per firm
  T_max <- max(times)

  # Time index relative to T_i for BC92
  t_rel <- numeric(n)
  for (ff in seq_along(firm_ids)) {
    idx <- firm_idx[[ff]]
    T_f <- max(times[idx])
    t_rel[idx] <- times[idx] - T_f  # t - T_i (negative or zero)
  }

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

  opt <- optim(
    par = start_params,
    fn = loglik_fn,
    method = "BFGS",
    control = ctrl,
    hessian = TRUE
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

    # JLMS-style conditional mean for BC92
    # E[u_i | eps_i1, ..., eps_iT]
    eps <- y - fitted_vals
    te <- .bc92_efficiency(eps, sigma_v, sigma_u, eta,
                           firms, firm_idx, T_i, t_rel)

    result <- list(
      coefficients = beta,
      sigma_v = as.numeric(sigma_v),
      sigma_u = as.numeric(sigma_u),
      eta = as.numeric(eta),
      logLik = opt$value,
      efficiency = te,
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
    te <- exp(-u_hat)

    result <- list(
      coefficients = beta,
      sigma_v = as.numeric(sigma_v),
      sigma_u = as.numeric(sigma_u),
      delta = delta,
      logLik = opt$value,
      efficiency = as.numeric(te),
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
#' u_it = u_i * exp(-eta*(t-T_i)), u_i ~ |N(0, sigma_u^2)|
#' Integrated over u_i (closed-form for half-normal).
#'
#' @keywords internal
#' @noRd
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

.bc92_efficiency <- function(eps, sigma_v, sigma_u, eta,
                             firms, firm_idx, T_i, t_rel) {
  n <- length(eps)
  te <- numeric(n)

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
    te[idx] <- exp(-u_it)
  }

  te
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
