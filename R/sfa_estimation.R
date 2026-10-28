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
#' @param estimator character. Technical efficiency estimator:
#'   \code{"bc88"} (Battese and Coelli, 1988; \code{E[exp(-u)|eps]})
#'   or \code{"jlms"} (Jondrow et al., 1982; \code{exp(-E[u|eps])}).
#'
#' @return A list with components: coefficients, sigma_v, sigma_u,
#'   logLik, efficiency, efficiency_jlms, efficiency_bc88, estimator,
#'   fitted, residuals, hessian, convergence.
#'   When Z variables are present, also includes delta, Z, and
#'   optionally sigma_u_vec or mu_vec.
#'
#' @keywords internal
#' @noRd
.fit_sfa_group <- function(formula, data, dist, control, ...,
                           estimator = c("bc88", "jlms"), sign = 1) {

  estimator <- match.arg(estimator)

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
  # Fit the production-equivalent model on sign * y (see .unflip_frontier)
  y <- sign * y

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

  # Check convergence. optim() can report success (code 0) at a point
  # where the log-likelihood is non-finite or impossibly large, i.e. at
  # a numerical failure rather than a maximum; flag those with code 99.
  ll_problem <- .sfa_loglik_problem(opt$value, exp(opt$par["log_sigma_v"]),
                                    n)
  if (!is.null(ll_problem)) {
    if (opt$convergence == 0L) opt$convergence <- 99L
    warning("SFA optimisation failed numerically (code ", opt$convergence,
            ": ", ll_problem, "; see ?check_convergence). The estimates ",
            "are not a maximum of the likelihood and should not be used.",
            call. = FALSE)
  } else if (opt$convergence != 0L) {
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
      sigma_star <- sigma_v
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
      sigma_star <- sigma_v
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
  te_jlms <- as.numeric(exp(-u_hat))

  # BC88 (Battese and Coelli, 1988) point efficiency: E[exp(-u)|eps]
  #   = exp(-mu* + sigma*^2/2) * Phi(mu*/sigma* - sigma*) / Phi(mu*/sigma*)
  # For the exponential distribution sigma* is sigma_v (set in the
  # branches above). Evaluated by .bc88_te(), which is stable in the
  # tails on platforms without long doubles.
  te_bc88 <- .bc88_te(as.numeric(mu_star), sigma_star, te_jlms)

  te <- if (estimator == "bc88") te_bc88 else te_jlms

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
    efficiency_jlms = te_jlms,
    efficiency_bc88 = te_bc88,
    estimator = estimator,
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

  .unflip_frontier(result, sign)
}


#' Map a production-equivalent fit back to the frontier's own scale
#'
#' A cost frontier y = x'beta + v + u is the production frontier
#' -y = x'(-beta) + v' - u with v' = -v, so the group fitters estimate
#' on sign * y (sign = 1 production, -1 cost) and this flips the
#' beta-type quantities back. Efficiencies, inefficiencies and the
#' variance and inefficiency-determinant parameters are unchanged.
#'
#' @keywords internal
#' @noRd
.unflip_frontier <- function(res, sign) {
  res$sign <- sign
  if (sign == 1) return(res)
  k <- length(res$coefficients)
  res$coefficients <- -res$coefficients
  res$fitted <- -res$fitted
  res$residuals <- -res$residuals
  res$y <- -res$y
  res$all_params <- .flip_k(res$all_params, k, sign)
  res$hessian <- .flip_k_mat(res$hessian, k, sign)
  res
}

# Flip the sign of the first k parameters (sign = -1), or of the
# matching rows and columns of a Hessian or covariance matrix
.flip_k <- function(x, k, sign) {
  if (sign == 1 || is.null(x)) return(x)
  x[seq_len(k)] <- -x[seq_len(k)]
  x
}
.flip_k_mat <- function(H, k, sign) {
  if (sign == 1 || is.null(H)) return(H)
  J <- rep(1, nrow(H))
  J[seq_len(k)] <- -1
  H * outer(J, J)
}


#' Battese-Coelli conditional efficiency E[exp(-c u) | eps]
#'
#' Computes E[exp(-c u) | eps] where u | eps is N(mu_star, sigma_star^2)
#' truncated at zero: the BC88 estimator (Battese and Coelli, 1988,
#' Eq. 6) for c = 1, and the BC92 time-varying form (Battese and
#' Coelli, 1992, Eq. 10) for c = d_t.
#'
#' The Phi ratio is evaluated on the log scale via
#' \code{pnorm(log.p = TRUE)} so it cannot underflow. On platforms
#' without long doubles the natural-scale ratio of two subnormal Phi
#' values can return 0, Inf, or values above one; the log scale avoids
#' this. In the far left tail (mu*/sigma* below -1e4) the two log Phi
#' terms cancel catastrophically, but there u | eps is asymptotically
#' exponential with rate -mu*/sigma*^2, giving the closed form
#' 1 / (1 + c sigma*^2 / (-mu*)). Any remaining non-finite values fall
#' back to the JLMS estimate.
#'
#' @param mu_star,sigma_star conditional posterior parameters of u.
#' @param te_jlms JLMS efficiencies, used as a last-resort fallback.
#' @param c_comp scalar or vector multiplier of u (1 for BC88, the
#'   temporal decay d_t for BC92).
#' @return numeric vector of efficiencies in (0, 1].
#' @keywords internal
#' @noRd
.bc88_te <- function(mu_star, sigma_star, te_jlms, c_comp = 1) {
  n <- max(length(mu_star), length(sigma_star), length(c_comp))
  mu_star <- rep_len(as.numeric(mu_star), n)
  sigma_star <- rep_len(as.numeric(sigma_star), n)
  c_comp <- rep_len(as.numeric(c_comp), n)

  ratio <- mu_star / sigma_star
  log_te <- -c_comp * mu_star + 0.5 * c_comp^2 * sigma_star^2 +
    pnorm(ratio - c_comp * sigma_star, log.p = TRUE) -
    pnorm(ratio, log.p = TRUE)
  te <- exp(pmin(log_te, 0))

  far <- which(ratio < -1e4)
  if (length(far)) {
    te[far] <- 1 / (1 + c_comp[far] * sigma_star[far]^2 / (-mu_star[far]))
  }

  bad <- !is.finite(te)
  te[bad] <- te_jlms[bad]
  te
}


#' Detect a numerically failed SFA log-likelihood
#'
#' In every SFA model of the package the noise v ~ N(0, sigma_v^2) is
#' added to an independent inefficiency term, so the density of each
#' observation is an average of normal densities with standard
#' deviation sigma_v and cannot exceed 1 / (sigma_v sqrt(2 pi)). A
#' log-likelihood above n (-log(sigma_v) - log(2 pi) / 2) is therefore
#' not a genuine value: it arises when rounding error in evaluating the
#' likelihood swamps the true value and the optimiser follows it. A
#' non-finite value, or the -1e20 the likelihood functions return in
#' place of one, likewise means no maximum was found.
#'
#' @param logLik log-likelihood at the optimiser's solution.
#' @param sigma_v estimated noise standard deviation (NA skips the
#'   bound check).
#' @param n number of observations.
#' @return \code{NULL} if the log-likelihood is admissible, otherwise a
#'   short description of the problem.
#' @keywords internal
#' @noRd
.sfa_loglik_problem <- function(logLik, sigma_v, n) {
  logLik <- as.numeric(logLik)[1]
  sigma_v <- as.numeric(sigma_v)[1]
  if (!is.finite(logLik) || logLik <= -1e20) {
    return("logLik not finite")
  }
  if (is.na(sigma_v) || n < 1L) return(NULL)
  bound <- n * (-log(sigma_v) - 0.5 * log(2 * pi))
  tol <- if (is.finite(bound)) 1e-8 * max(1, abs(bound)) else 0
  if (logLik > bound + tol) {
    return(sprintf("logLik %.3g exceeds bound %.3g", logLik, bound))
  }
  NULL
}


# ==== Homoscedastic log-likelihoods (no Z) ====

# Parameter-vector layout convention (all likelihoods below):
# params = c(beta[1:k], log_sigma_v, <distribution-specific tail>),
# where the tail is log_sigma_u (half-normal, exponential),
# mu then log_sigma_u (truncated-normal), or delta[1:p] followed by
# log_sigma_u where Z variables enter.

# Normal/half-normal composed error eps = v - u, u ~ |N(0, sigma_u^2)|:
# log f(eps) = log 2 - log sigma + log phi(eps/sigma)
#              + log Phi(-eps*lambda/sigma), sigma^2 = sigma_v^2 + sigma_u^2,
# lambda = sigma_u/sigma_v. Aigner, Lovell and Schmidt (1977, Eq. 8).
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

# Normal/truncated-normal, u ~ N+(mu, sigma_u^2):
# log f(eps) = log phi((eps + mu)/sigma) - log sigma
#              + log Phi(mu*/sigma*) - log Phi(mu/sigma_u),
# with mu* = (mu*sigma_v^2 - eps*sigma_u^2)/sigma^2 and
# sigma* = sigma_v*sigma_u/sigma. Stevenson (1980, Eq. 15).
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

# Normal/exponential, u ~ Exp(rate = 1/sigma_u):
# log f(eps) = -log sigma_u + eps/sigma_u + sigma_v^2/(2*sigma_u^2)
#              + log Phi(-a),   a = (eps + sigma_v^2/sigma_u)/sigma_v.
# Meeusen and van den Broeck (1977).
#
# For large a, log Phi(-a) is about -a^2/2 and cancels the positive
# terms before it: at sigma_v/sigma_u = 1e57 both sides are ~1e113, and
# the rounding error of their sum (~1e97, of either sign) replaced the
# true value (~ -200), so BFGS climbed it to log-likelihoods of 1e100.
# For a > 40 the density is evaluated in the equivalent form
#   log f(eps) = -log(2*pi)/2 - eps^2/(2*sigma_v^2) - log(sigma_u*a)
#                + log(a*M(a)),
# with M(a) = Phi(-a)/phi(a) the Mills ratio and the asymptotic series
# a*M(a) = 1 - 1/a^2 + 3/a^4 - 15/a^6 + 105/a^8, accurate to 1e-13 for
# a > 40, the size of the rounding error of the direct form at a = 40.
# Below the switch the direct form is used unchanged.
# Vectorised over eps and sigma_u; returns per-observation values.
#' @noRd
.ll_exponential <- function(eps, sigma_v, sigma_u) {
  lambda_u <- 1 / sigma_u
  a <- (eps + lambda_u * sigma_v^2) / sigma_v

  ll <- log(lambda_u) + lambda_u * eps + 0.5 * lambda_u^2 * sigma_v^2 +
    pnorm(-a, log.p = TRUE)

  big <- which(a > 40)
  if (length(big)) {
    ab <- a[big]
    r <- 1 / ab^2
    su <- rep_len(sigma_u, length(eps))[big]
    ll[big] <- -0.5 * log(2 * pi) - 0.5 * (eps[big] / sigma_v)^2 -
      log(su * ab) + log1p(r * (-1 + r * (3 + r * (-15 + 105 * r))))
  }
  ll
}

#' @noRd
.loglik_exponential <- function(params, y, X) {
  k <- ncol(X)
  beta <- params[seq_len(k)]
  sigma_v <- exp(params["log_sigma_v"])
  sigma_u <- exp(params["log_sigma_u"])

  if (sigma_v < 1e-10 || sigma_u < 1e-10) return(-1e20)

  eps <- as.numeric(y - X %*% beta)

  result <- sum(.ll_exponential(eps, sigma_v, sigma_u))
  if (!is.finite(result)) return(-1e20)
  result
}


# ==== Heteroscedastic log-likelihoods (with Z) ====
# Same densities as above with the distribution parameter made
# observation-specific: sigma_u_i = exp(Z_i' delta) for half-normal and
# exponential, mu_i = Z_i' delta for truncated-normal.

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

  eps <- as.numeric(y - X %*% beta)

  result <- sum(.ll_exponential(eps, sigma_v, sigma_u))
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

  eps <- as.numeric(y - X %*% beta)

  .ll_exponential(eps, sigma_v, sigma_u)
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

  eps <- as.numeric(y - X %*% beta)

  .ll_exponential(eps, sigma_v, sigma_u)
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
  X <- group_model$X
  has_z <- !is.null(group_model$Z)
  k <- ncol(X)
  # Evaluate on the production-equivalent scale of a cost fit (see
  # .unflip_frontier); the beta columns are mapped back below
  sign <- if (is.null(group_model$sign)) 1 else group_model$sign
  params <- .flip_k(group_model$all_params, k, sign)
  y <- sign * group_model$y

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

  S <- numDeriv::jacobian(obs_fn, params, method.args = list(eps = 1e-4))
  # Chain rule: d theta_internal / d theta = -1 for the betas of a cost fit
  S[, seq_len(k)] <- sign * S[, seq_len(k)]
  S
}
