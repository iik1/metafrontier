#' Internal SFA estimation for a single group
#'
#' Fits a stochastic frontier model to a single group using MLE.
#' Implements the normal/half-normal, normal/truncated-normal, and
#' normal/exponential models. Supports heteroscedastic inefficiency
#' when the formula contains a second RHS part (separated by \code{|}).
#'
#' @param formula a \code{Formula} object. If it has two RHS parts
#'   (e.g., \code{y ~ x1 + x2 | z1 + z2}), the second part specifies
#'   inefficiency determinants.
#' @param data data frame for this group.
#' @param dist distribution of the inefficiency term.
#' @param control list of control parameters.
#' @param ... additional arguments.
#'
#' @return A list with components: coefficients, sigma_v, sigma_u,
#'   logLik, efficiency, fitted, residuals, hessian, convergence.
#'   When Z variables are present, also includes delta, Z, and
#'   optionally sigma_u_vec or mu_vec.
#'
#' @keywords internal
#' @noRd
.fit_sfa_group <- function(formula, data, dist, control, ...) {

  # Build model frame from the full Formula
  if (inherits(formula, "Formula")) {
    f <- formula
    has_z <- length(f)[2] >= 2L

    if (has_z) {
      # Full model frame including all RHS parts
      mf <- model.frame(f, data = data, na.action = na.omit)
      y <- model.response(mf)
      X <- model.matrix(f, data = mf, rhs = 1)
      Z <- model.matrix(f, data = mf, rhs = 2)
      if (ncol(Z) == 0L) {
        Z <- NULL
        has_z <- FALSE
      }
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

  # OLS starting values
  ols <- lm.fit(X, y)
  ols_resid <- ols$residuals
  sigma_ols <- sqrt(sum(ols_resid^2) / (n - k))

  # Starting values for sigma_v and sigma_u
  sigma_u_start <- max(sigma_ols * 0.5, 0.1)
  sigma_v_start <- max(sigma_ols * 0.5, 0.1)

  # ---- Set up parameters and log-likelihood ----
  if (has_z) {
    p <- ncol(Z)
    delta_start <- rep(0, p)

    if (dist == "hnormal") {
      # sigma_u_i = exp(Z_i' delta)
      delta_start[1] <- log(sigma_u_start)
      start_params <- c(ols$coefficients,
                        log_sigma_v = log(sigma_v_start),
                        setNames(delta_start, paste0("d_", seq_len(p))))
      loglik_fn <- function(params, y, X)
        .loglik_hnormal_z(params, y, X, Z, k)
    } else if (dist == "exponential") {
      # sigma_u_i = exp(Z_i' delta)
      delta_start[1] <- log(sigma_u_start)
      start_params <- c(ols$coefficients,
                        log_sigma_v = log(sigma_v_start),
                        setNames(delta_start, paste0("d_", seq_len(p))))
      loglik_fn <- function(params, y, X)
        .loglik_exponential_z(params, y, X, Z, k)
    } else if (dist == "tnormal") {
      # mu_i = Z_i' delta; sigma_u remains scalar
      start_params <- c(ols$coefficients,
                        log_sigma_v = log(sigma_v_start),
                        setNames(delta_start, paste0("d_", seq_len(p))),
                        log_sigma_u = log(sigma_u_start))
      loglik_fn <- function(params, y, X)
        .loglik_tnormal_z(params, y, X, Z, k)
    }
  } else {
    if (dist == "hnormal") {
      start_params <- c(ols$coefficients,
                        log_sigma_v = log(sigma_v_start),
                        log_sigma_u = log(sigma_u_start))
      loglik_fn <- .loglik_hnormal
    } else if (dist == "exponential") {
      start_params <- c(ols$coefficients,
                        log_sigma_v = log(sigma_v_start),
                        log_sigma_u = log(sigma_u_start))
      loglik_fn <- .loglik_exponential
    } else if (dist == "tnormal") {
      start_params <- c(ols$coefficients,
                        log_sigma_v = log(sigma_v_start),
                        mu = 0,
                        log_sigma_u = log(sigma_u_start))
      loglik_fn <- .loglik_tnormal
    }
  }

  # MLE via optim (maximise log-likelihood)
  ctrl <- list(fnscale = -1, maxit = 5000, reltol = 1e-10)
  ctrl[names(control)] <- control

  opt <- tryCatch(
    optim(par = start_params, fn = loglik_fn, y = y, X = X,
          method = "BFGS", control = ctrl, hessian = TRUE),
    error = function(e) {
      tryCatch(
        optim(par = start_params, fn = loglik_fn, y = y, X = X,
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

  # Check convergence
  if (opt$convergence != 0L) {
    warning("SFA optimisation did not converge (code ", opt$convergence,
            "). Results may be unreliable. Consider increasing maxit ",
            "via control or checking the data for outliers.", call. = FALSE)
  }

  # Extract results
  beta_hat <- opt$par[seq_len(k)]
  names(beta_hat) <- colnames(X)

  sigma_v <- exp(opt$par["log_sigma_v"])
  eps <- y - X %*% beta_hat

  # Mills ratio helper
  .mills <- function(z) {
    pp <- pnorm(z)
    pp <- pmax(pp, .Machine$double.eps)
    dnorm(z) / pp
  }

  # ---- Extract sigma_u and compute JLMS efficiency ----
  if (has_z) {
    p <- ncol(Z)
    delta <- opt$par[(k + 2):(k + 1 + p)]
    names(delta) <- colnames(Z)

    if (dist == "hnormal") {
      sigma_u_vec <- exp(as.numeric(Z %*% delta))
      sigma_u <- mean(sigma_u_vec)
      sigma_sq <- sigma_v^2 + sigma_u_vec^2
      mu_star <- -as.numeric(eps) * sigma_u_vec^2 / sigma_sq
      sigma_star <- sigma_v * sigma_u_vec / sqrt(sigma_sq)
      u_hat <- mu_star + sigma_star * .mills(mu_star / sigma_star)
    } else if (dist == "exponential") {
      sigma_u_vec <- exp(as.numeric(Z %*% delta))
      sigma_u <- mean(sigma_u_vec)
      mu_star <- -as.numeric(eps) - sigma_v^2 / sigma_u_vec
      u_hat <- mu_star + sigma_v * .mills(mu_star / sigma_v)
    } else if (dist == "tnormal") {
      mu_vec <- as.numeric(Z %*% delta)
      sigma_u <- exp(opt$par["log_sigma_u"])
      sigma_u_vec <- NULL
      sigma_sq <- sigma_v^2 + sigma_u^2
      mu_star <- (mu_vec * sigma_v^2 -
                    as.numeric(eps) * sigma_u^2) / sigma_sq
      sigma_star <- sigma_v * sigma_u / sqrt(sigma_sq)
      u_hat <- mu_star + sigma_star * .mills(mu_star / sigma_star)
    }
  } else {
    sigma_u <- exp(opt$par["log_sigma_u"])
    sigma_u_vec <- NULL
    delta <- NULL

    if (dist == "hnormal") {
      sigma_sq <- sigma_v^2 + sigma_u^2
      mu_star <- -as.numeric(eps) * sigma_u^2 / sigma_sq
      sigma_star <- sigma_v * sigma_u / sqrt(sigma_sq)
      u_hat <- mu_star + sigma_star * .mills(mu_star / sigma_star)
    } else if (dist == "exponential") {
      mu_star <- -as.numeric(eps) - sigma_v^2 / sigma_u
      u_hat <- mu_star + sigma_v * .mills(mu_star / sigma_v)
    } else if (dist == "tnormal") {
      mu_val <- opt$par["mu"]
      sigma_sq <- sigma_v^2 + sigma_u^2
      mu_star <- (mu_val * sigma_v^2 -
                    as.numeric(eps) * sigma_u^2) / sigma_sq
      sigma_star <- sigma_v * sigma_u / sqrt(sigma_sq)
      u_hat <- mu_star + sigma_star * .mills(mu_star / sigma_star)
    }
  }

  # Ensure non-negative inefficiency
  u_hat <- pmax(as.numeric(u_hat), 0)

  # JLMS (Jondrow et al., 1982) point efficiency: exp(-E[u|eps])
  te <- as.numeric(exp(-u_hat))

  # Frontier values
  fitted_vals <- as.numeric(X %*% beta_hat)

  lambda <- sigma_u / sigma_v

  result <- list(
    coefficients = beta_hat,
    sigma_v = sigma_v,
    sigma_u = sigma_u,
    lambda = lambda,
    logLik = opt$value,
    efficiency = te,
    inefficiency = as.numeric(u_hat),
    fitted = fitted_vals,
    residuals = as.numeric(eps),
    hessian = opt$hessian,
    convergence = opt$convergence,
    dist = dist,
    nobs = n,
    X = X,
    y = y,
    all_params = opt$par
  )

  if (has_z) {
    result$Z <- Z
    result$delta <- delta
    if (!is.null(sigma_u_vec)) result$sigma_u_vec <- sigma_u_vec
    if (dist == "tnormal") result$mu_vec <- mu_vec
  }

  result
}


# ==== Homoscedastic log-likelihoods (no Z) ====

#' @noRd
.loglik_hnormal <- function(params, y, X) {
  k <- ncol(X)
  beta <- params[seq_len(k)]
  sigma_v <- exp(params["log_sigma_v"])
  sigma_u <- exp(params["log_sigma_u"])

  if (sigma_v < 1e-10 || sigma_u < 1e-10) return(-1e20)

  sigma_sq <- sigma_v^2 + sigma_u^2
  sigma <- sqrt(sigma_sq)
  lambda <- sigma_u / sigma_v

  eps <- as.numeric(y - X %*% beta)

  ll <- -0.5 * log(2 * pi) + log(2) - log(sigma) -
    0.5 * (eps / sigma)^2 +
    pnorm(-eps * lambda / sigma, log.p = TRUE)

  result <- sum(ll)
  if (!is.finite(result)) return(-1e20)
  result
}

#' @noRd
.loglik_tnormal <- function(params, y, X) {
  k <- ncol(X)
  beta <- params[seq_len(k)]
  sigma_v <- exp(params["log_sigma_v"])
  mu <- params["mu"]
  sigma_u <- exp(params["log_sigma_u"])

  if (sigma_v < 1e-10 || sigma_u < 1e-10) return(-1e20)

  sigma_sq <- sigma_v^2 + sigma_u^2
  sigma <- sqrt(sigma_sq)

  eps <- as.numeric(y - X %*% beta)
  mu_star <- (mu * sigma_v^2 - eps * sigma_u^2) / sigma_sq
  sigma_star <- sigma_v * sigma_u / sqrt(sigma_sq)

  ll <- -0.5 * log(2 * pi) - log(sigma) -
    0.5 * ((eps + mu)^2 / sigma_sq) +
    pnorm(mu_star / sigma_star, log.p = TRUE) -
    pnorm(mu / sigma_u, log.p = TRUE)

  result <- sum(ll)
  if (!is.finite(result)) return(-1e20)
  result
}

#' @noRd
.loglik_exponential <- function(params, y, X) {
  k <- ncol(X)
  beta <- params[seq_len(k)]
  sigma_v <- exp(params["log_sigma_v"])
  sigma_u <- exp(params["log_sigma_u"])

  if (sigma_v < 1e-10 || sigma_u < 1e-10) return(-1e20)

  lambda_u <- 1 / sigma_u
  eps <- as.numeric(y - X %*% beta)

  ll <- log(lambda_u) + lambda_u * eps + 0.5 * lambda_u^2 * sigma_v^2 +
    pnorm(-(eps + lambda_u * sigma_v^2) / sigma_v, log.p = TRUE)

  result <- sum(ll)
  if (!is.finite(result)) return(-1e20)
  result
}


# ==== Heteroscedastic log-likelihoods (with Z) ====

#' Half-normal with observation-specific sigma_u_i = exp(Z_i' delta)
#' @noRd
.loglik_hnormal_z <- function(params, y, X, Z, k) {
  beta <- params[seq_len(k)]
  sigma_v <- exp(params[k + 1L])
  p <- ncol(Z)
  delta <- params[(k + 2L):(k + 1L + p)]

  if (sigma_v < 1e-10) return(-1e20)

  sigma_u <- exp(as.numeric(Z %*% delta))
  if (any(sigma_u < 1e-10)) return(-1e20)

  sigma_sq <- sigma_v^2 + sigma_u^2
  sigma <- sqrt(sigma_sq)
  lambda <- sigma_u / sigma_v

  eps <- as.numeric(y - X %*% beta)

  ll <- -0.5 * log(2 * pi) + log(2) - log(sigma) -
    0.5 * (eps / sigma)^2 +
    pnorm(-eps * lambda / sigma, log.p = TRUE)

  result <- sum(ll)
  if (!is.finite(result)) return(-1e20)
  result
}

#' Truncated-normal with observation-specific mu_i = Z_i' delta
#' @noRd
.loglik_tnormal_z <- function(params, y, X, Z, k) {
  beta <- params[seq_len(k)]
  sigma_v <- exp(params[k + 1L])
  p <- ncol(Z)
  delta <- params[(k + 2L):(k + 1L + p)]
  sigma_u <- exp(params[k + 1L + p + 1L])

  if (sigma_v < 1e-10 || sigma_u < 1e-10) return(-1e20)

  mu <- as.numeric(Z %*% delta)
  sigma_sq <- sigma_v^2 + sigma_u^2
  sigma <- sqrt(sigma_sq)

  eps <- as.numeric(y - X %*% beta)
  mu_star <- (mu * sigma_v^2 - eps * sigma_u^2) / sigma_sq
  sigma_star <- sigma_v * sigma_u / sqrt(sigma_sq)

  ll <- -0.5 * log(2 * pi) - log(sigma) -
    0.5 * ((eps + mu)^2 / sigma_sq) +
    pnorm(mu_star / sigma_star, log.p = TRUE) -
    pnorm(mu / sigma_u, log.p = TRUE)

  result <- sum(ll)
  if (!is.finite(result)) return(-1e20)
  result
}

#' Exponential with observation-specific sigma_u_i = exp(Z_i' delta)
#' @noRd
.loglik_exponential_z <- function(params, y, X, Z, k) {
  beta <- params[seq_len(k)]
  sigma_v <- exp(params[k + 1L])
  p <- ncol(Z)
  delta <- params[(k + 2L):(k + 1L + p)]

  if (sigma_v < 1e-10) return(-1e20)

  sigma_u <- exp(as.numeric(Z %*% delta))
  if (any(sigma_u < 1e-10)) return(-1e20)

  lambda_u <- 1 / sigma_u
  eps <- as.numeric(y - X %*% beta)

  ll <- log(lambda_u) + lambda_u * eps + 0.5 * lambda_u^2 * sigma_v^2 +
    pnorm(-(eps + lambda_u * sigma_v^2) / sigma_v, log.p = TRUE)

  result <- sum(ll)
  if (!is.finite(result)) return(-1e20)
  result
}


# ==== Observation-level log-likelihoods (return vector, not sum) ====
# Used by Murphy-Topel correction and latent class EM

#' @noRd
.loglik_hnormal_obs <- function(params, y, X) {
  k <- ncol(X)
  beta <- params[seq_len(k)]
  sigma_v <- exp(params["log_sigma_v"])
  sigma_u <- exp(params["log_sigma_u"])

  if (sigma_v < 1e-10 || sigma_u < 1e-10) return(rep(-1e20, length(y)))

  sigma_sq <- sigma_v^2 + sigma_u^2
  sigma <- sqrt(sigma_sq)
  lambda <- sigma_u / sigma_v
  eps <- as.numeric(y - X %*% beta)

  -0.5 * log(2 * pi) + log(2) - log(sigma) -
    0.5 * (eps / sigma)^2 +
    pnorm(-eps * lambda / sigma, log.p = TRUE)
}

#' @noRd
.loglik_tnormal_obs <- function(params, y, X) {
  k <- ncol(X)
  beta <- params[seq_len(k)]
  sigma_v <- exp(params["log_sigma_v"])
  mu <- params["mu"]
  sigma_u <- exp(params["log_sigma_u"])

  if (sigma_v < 1e-10 || sigma_u < 1e-10) return(rep(-1e20, length(y)))

  sigma_sq <- sigma_v^2 + sigma_u^2
  sigma <- sqrt(sigma_sq)
  eps <- as.numeric(y - X %*% beta)
  mu_star <- (mu * sigma_v^2 - eps * sigma_u^2) / sigma_sq
  sigma_star <- sigma_v * sigma_u / sqrt(sigma_sq)

  -0.5 * log(2 * pi) - log(sigma) -
    0.5 * ((eps + mu)^2 / sigma_sq) +
    pnorm(mu_star / sigma_star, log.p = TRUE) -
    pnorm(mu / sigma_u, log.p = TRUE)
}

#' @noRd
.loglik_exponential_obs <- function(params, y, X) {
  k <- ncol(X)
  beta <- params[seq_len(k)]
  sigma_v <- exp(params["log_sigma_v"])
  sigma_u <- exp(params["log_sigma_u"])

  if (sigma_v < 1e-10 || sigma_u < 1e-10) return(rep(-1e20, length(y)))

  lambda_u <- 1 / sigma_u
  eps <- as.numeric(y - X %*% beta)

  log(lambda_u) + lambda_u * eps + 0.5 * lambda_u^2 * sigma_v^2 +
    pnorm(-(eps + lambda_u * sigma_v^2) / sigma_v, log.p = TRUE)
}

#' @noRd
.loglik_hnormal_z_obs <- function(params, y, X, Z, k) {
  beta <- params[seq_len(k)]
  sigma_v <- exp(params[k + 1L])
  p <- ncol(Z)
  delta <- params[(k + 2L):(k + 1L + p)]

  if (sigma_v < 1e-10) return(rep(-1e20, length(y)))

  sigma_u <- exp(as.numeric(Z %*% delta))
  if (any(sigma_u < 1e-10)) return(rep(-1e20, length(y)))

  sigma_sq <- sigma_v^2 + sigma_u^2
  sigma <- sqrt(sigma_sq)
  lambda <- sigma_u / sigma_v
  eps <- as.numeric(y - X %*% beta)

  -0.5 * log(2 * pi) + log(2) - log(sigma) -
    0.5 * (eps / sigma)^2 +
    pnorm(-eps * lambda / sigma, log.p = TRUE)
}

#' @noRd
.loglik_tnormal_z_obs <- function(params, y, X, Z, k) {
  beta <- params[seq_len(k)]
  sigma_v <- exp(params[k + 1L])
  p <- ncol(Z)
  delta <- params[(k + 2L):(k + 1L + p)]
  sigma_u <- exp(params[k + 1L + p + 1L])

  if (sigma_v < 1e-10 || sigma_u < 1e-10) return(rep(-1e20, length(y)))

  mu <- as.numeric(Z %*% delta)
  sigma_sq <- sigma_v^2 + sigma_u^2
  sigma <- sqrt(sigma_sq)
  eps <- as.numeric(y - X %*% beta)
  mu_star <- (mu * sigma_v^2 - eps * sigma_u^2) / sigma_sq
  sigma_star <- sigma_v * sigma_u / sqrt(sigma_sq)

  -0.5 * log(2 * pi) - log(sigma) -
    0.5 * ((eps + mu)^2 / sigma_sq) +
    pnorm(mu_star / sigma_star, log.p = TRUE) -
    pnorm(mu / sigma_u, log.p = TRUE)
}

#' @noRd
.loglik_exponential_z_obs <- function(params, y, X, Z, k) {
  beta <- params[seq_len(k)]
  sigma_v <- exp(params[k + 1L])
  p <- ncol(Z)
  delta <- params[(k + 2L):(k + 1L + p)]

  if (sigma_v < 1e-10) return(rep(-1e20, length(y)))

  sigma_u <- exp(as.numeric(Z %*% delta))
  if (any(sigma_u < 1e-10)) return(rep(-1e20, length(y)))

  lambda_u <- 1 / sigma_u
  eps <- as.numeric(y - X %*% beta)

  log(lambda_u) + lambda_u * eps + 0.5 * lambda_u^2 * sigma_v^2 +
    pnorm(-(eps + lambda_u * sigma_v^2) / sigma_v, log.p = TRUE)
}


# ==== Score vector computation ====

#' Compute observation-level score matrix for a fitted group model
#'
#' Returns an n x p matrix where each row is the gradient of that
#' observation's log-likelihood contribution w.r.t. all parameters.
#'
#' @param group_model A fitted group model from .fit_sfa_group().
#' @return An n x p numeric matrix.
#' @keywords internal
#' @noRd
.score_vector_sfa <- function(group_model) {
  params <- group_model$all_params
  y <- group_model$y
  X <- group_model$X
  has_z <- !is.null(group_model$Z)
  k <- ncol(X)

  if (has_z) {
    Z <- group_model$Z
    obs_fn <- switch(group_model$dist,
      hnormal     = function(p) .loglik_hnormal_z_obs(p, y, X, Z, k),
      tnormal     = function(p) .loglik_tnormal_z_obs(p, y, X, Z, k),
      exponential = function(p) .loglik_exponential_z_obs(p, y, X, Z, k)
    )
  } else {
    obs_fn <- switch(group_model$dist,
      hnormal     = function(p) .loglik_hnormal_obs(p, y, X),
      tnormal     = function(p) .loglik_tnormal_obs(p, y, X),
      exponential = function(p) .loglik_exponential_obs(p, y, X)
    )
  }

  numDeriv::jacobian(obs_fn, params,
                     method.args = list(eps = 1e-4))
}
