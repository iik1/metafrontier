#' Spatial-autoregressive (SAR) SFA for a single group
#'
#' y = rho W y + x'beta + v - u (Glass, Kenjegalieva and Sickles, 2016),
#' estimated by maximum likelihood. For a given rho the model is an
#' ordinary stochastic frontier for y - rho W y, so rho is profiled: each
#' evaluation fits .fit_sfa_group() to y - rho W y and adds the
#' log-Jacobian log|I - rho W| = sum log(1 - rho lambda_i), with lambda
#' the eigenvalues of W computed once. The reported Hessian is the
#' numerical Hessian of the full log-likelihood in (beta, log_sigma_v,
#' log_sigma_u[, mu], rho).
#'
#' Direct, indirect and total effects of each input (LeSage and Pace,
#' 2009) are the average diagonal, average row sum minus diagonal, and
#' average row sum of (I - rho W)^-1 times beta_k, with delta-method
#' standard errors.
#'
#' @keywords internal
#' @noRd
.fit_sfa_spatial_group <- function(formula, data, W, dist, control,
                                   estimator = "bc88", sign = 1) {
  f_base <- formula(Formula::Formula(formula), rhs = 1)
  y <- model.response(model.frame(f_base, data = data))
  X <- model.matrix(f_base, data = data)
  k <- ncol(X)
  n <- length(y)
  Wy <- as.numeric(W %*% y)
  ev <- eigen(W, only.values = TRUE)$values
  logdet <- function(rho) Re(sum(log(1 - rho * ev)))

  # The frontier for a given rho, on the response y - rho W y
  data$.y_sar <- 0
  f_star <- Formula::Formula(stats::update(f_base, .y_sar ~ .))
  fit_at <- function(rho) {
    data$.y_sar <- y - rho * Wy
    .fit_sfa_group(f_star, data, dist, control, estimator = estimator,
                   sign = sign)
  }

  re <- Re(ev)
  if (all(abs(ev) < 1e-12)) {
    rho <- 0   # no neighbours: rho is not identified, nest the SFA
  } else {
    lo <- if (min(re) < 0) 1 / min(re) else -1
    hi <- 1 / max(re)
    profile <- function(r) suppressWarnings(logdet(r) + fit_at(r)$logLik)
    rho <- stats::optimize(profile, c(lo, hi) * (1 - 1e-6),
                           maximum = TRUE, tol = 1e-7)$maximum
  }
  fit <- fit_at(rho)

  if (rho != 0) {
    llf <- switch(dist, hnormal = .loglik_hnormal,
                  tnormal = .loglik_tnormal,
                  exponential = .loglik_exponential)
    full_ll <- function(p) {
      r <- p[length(p)]
      core <- .flip_k(p[-length(p)], k, sign)
      logdet(r) + llf(core, sign * (y - r * Wy), X)
    }
    par <- c(fit$all_params, rho = rho)
    fit$hessian <- numDeriv::hessian(full_ll, par)
    dimnames(fit$hessian) <- list(names(par), names(par))
    fit$all_params <- par
    fit$logLik <- fit$logLik + logdet(rho)
  }

  fit$rho <- rho
  fit$y <- as.numeric(y)
  fit$W <- W
  fit$spatial_lag <- rho * Wy
  fit$effects <- .sar_effects(fit, W, k)
  fit
}


# Direct, indirect and total effects of each non-intercept regressor
.sar_effects <- function(fit, W, k) {
  beta <- fit$coefficients
  inputs <- setdiff(names(beta), "(Intercept)")
  eff <- function(b, rho) {
    S <- solve(diag(nrow(W)) - rho * W)
    direct <- mean(diag(S)) * b[inputs]
    total <- mean(rowSums(S)) * b[inputs]
    c(direct, total - direct, total)
  }
  est <- eff(beta, fit$rho)
  se <- rep(NA_real_, length(est))
  V <- tryCatch(solve(-fit$hessian), error = function(e) NULL)
  if (!is.null(V) && fit$rho != 0) {
    idx <- c(seq_len(k), nrow(V))   # beta and rho
    J <- numDeriv::jacobian(function(p) eff(p[seq_len(k)], p[k + 1L]),
                            c(beta, fit$rho))
    se <- sqrt(pmax(diag(J %*% V[idx, idx] %*% t(J)), 0))
  }
  m <- length(inputs)
  data.frame(input = inputs,
             direct = est[seq_len(m)], indirect = est[m + seq_len(m)],
             total = est[2 * m + seq_len(m)],
             se_direct = se[seq_len(m)], se_indirect = se[m + seq_len(m)],
             se_total = se[2 * m + seq_len(m)],
             row.names = NULL)
}


# Validate and convert a spatial weights matrix for n observations
.spatial_weights <- function(W, n) {
  if (inherits(W, "listw")) {
    if (!requireNamespace("spdep", quietly = TRUE)) {
      stop("Package 'spdep' is required for 'listw' weights.", call. = FALSE)
    }
    W <- spdep::listw2mat(W)
  }
  W <- as.matrix(W)
  if (!is.numeric(W) || nrow(W) != n || ncol(W) != n) {
    stop("'spatial$W' must be a numeric n x n matrix with n = nrow(data) ",
         "= ", n, ".", call. = FALSE)
  }
  if (anyNA(W)) stop("'spatial$W' contains missing values.", call. = FALSE)
  if (any(diag(W) != 0)) {
    stop("'spatial$W' must have a zero diagonal.", call. = FALSE)
  }
  unname(W)
}
