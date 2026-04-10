#' Internal SFA estimation for a single group
#'
#' Fits a stochastic frontier model to a single group using MLE.
#' Implements the normal/half-normal, normal/truncated-normal, and
#' normal/exponential models.
#'
#' @param formula a \code{Formula} object.
#' @param data data frame for this group.
#' @param dist distribution of the inefficiency term.
#' @param control list of control parameters.
#' @param ... additional arguments.
#'
#' @return A list with components: coefficients, sigma_v, sigma_u,
#'   logLik, efficiency, fitted, residuals, hessian, convergence.
#'
#' @keywords internal
#' @noRd
.fit_sfa_group <- function(formula, data, dist, control, ...) {

  # Extract model matrices using base R (Formula model.matrix can
  # cause issues, so we convert to a standard formula for extraction)
  if (inherits(formula, "Formula")) {
    # Extract just the first RHS part as a base formula
    f_base <- formula(formula, rhs = 1)
  } else {
    f_base <- formula
  }

  mf <- model.frame(f_base, data = data, na.action = na.omit)
  y <- model.response(mf)
  X <- model.matrix(f_base, data = mf)
  n <- length(y)
  k <- ncol(X)

  # OLS starting values
  ols <- lm.fit(X, y)
  ols_resid <- ols$residuals
  sigma_ols <- sqrt(sum(ols_resid^2) / (n - k))

  # Starting values for sigma_v and sigma_u
  # Use conservative starting values to avoid numerical issues
  sigma_u_start <- max(sigma_ols * 0.5, 0.1)
  sigma_v_start <- max(sigma_ols * 0.5, 0.1)

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

  # MLE via optim (maximise log-likelihood)
  # First try BFGS; fall back to Nelder-Mead if non-finite values
  ctrl <- list(fnscale = -1, maxit = 5000, reltol = 1e-10)
  ctrl[names(control)] <- control

  opt <- tryCatch(
    optim(par = start_params, fn = loglik_fn, y = y, X = X,
          method = "BFGS", control = ctrl, hessian = TRUE),
    error = function(e) {
      tryCatch(
        # Fall back to Nelder-Mead (more robust, no gradient needed)
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

  # Extract results
  beta_hat <- opt$par[seq_len(k)]
  names(beta_hat) <- colnames(X)

  sigma_v <- exp(opt$par["log_sigma_v"])
  sigma_u <- exp(opt$par["log_sigma_u"])
  sigma_sq <- sigma_v^2 + sigma_u^2
  lambda <- sigma_u / sigma_v

  # Composed residuals
  eps <- y - X %*% beta_hat

  # Conditional inefficiency (JLMS estimator)
  # Helper: Mills ratio dnorm(z)/pnorm(z) with numerical guard
  .mills <- function(z) {
    p <- pnorm(z)
    p <- pmax(p, .Machine$double.eps)
    dnorm(z) / p
  }

  if (dist == "hnormal") {
    mu_star <- -eps * sigma_u^2 / sigma_sq
    sigma_star <- sigma_v * sigma_u / sqrt(sigma_sq)
    u_hat <- mu_star + sigma_star * .mills(mu_star / sigma_star)
  } else if (dist == "exponential") {
    mu_star <- -eps - sigma_v^2 / sigma_u
    u_hat <- mu_star + sigma_v * .mills(mu_star / sigma_v)
  } else if (dist == "tnormal") {
    mu_val <- opt$par["mu"]
    mu_star <- (mu_val * sigma_v^2 - eps * sigma_u^2) / sigma_sq
    sigma_star <- sigma_v * sigma_u / sqrt(sigma_sq)
    u_hat <- mu_star + sigma_star * .mills(mu_star / sigma_star)
  }

  # Ensure non-negative inefficiency
  u_hat <- pmax(as.numeric(u_hat), 0)

  # Battese-Coelli (1988) efficiency estimator
  te <- as.numeric(exp(-u_hat))

  # Frontier values
  fitted_vals <- as.numeric(X %*% beta_hat)

  list(
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
}


# ---- Log-likelihood: Normal/Half-Normal ----
#' @noRd
.loglik_hnormal <- function(params, y, X) {
  k <- ncol(X)
  beta <- params[seq_len(k)]
  sigma_v <- exp(params["log_sigma_v"])
  sigma_u <- exp(params["log_sigma_u"])

  # Guard against degenerate values
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

# ---- Log-likelihood: Normal/Truncated-Normal ----
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

# ---- Log-likelihood: Normal/Exponential ----
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
