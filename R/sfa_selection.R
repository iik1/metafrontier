#' Sample-selection SFA for a single group (Greene, 2010)
#'
#' Two-step estimator. Step 1 is a probit for the selection indicator d
#' on Z. Step 2 maximises, given the probit index, the likelihood of the
#' selected observations
#'   L_i = E_u[ phi(v_i / sigma_v) / sigma_v *
#'              Phi((z_i'gamma + rho v_i / sigma_v) / sqrt(1 - rho^2)) ],
#' with v_i = y_i - x_i'beta + u, corr(w, v) = rho and u independent
#' half-normal or exponential. The expectation over u is evaluated by
#' Gauss-Hermite quadrature, as in sfaR::sfaselectioncross(). Second-step
#' standard errors are not corrected for the first step.
#'
#' Rows of `data` with d = 0 enter only the probit; the returned
#' efficiencies, fitted values and residuals are for the selected rows,
#' in data order.
#'
#' @keywords internal
#' @noRd
.fit_sfa_selection_group <- function(formula, selection, data, dist,
                                     control, estimator = "bc88",
                                     sign = 1) {
  if (!dist %in% c("hnormal", "exponential")) {
    stop("The sample-selection model supports dist = \"hnormal\" or ",
         "\"exponential\".", call. = FALSE)
  }

  # Step 1: probit on all rows of the group
  smf <- model.frame(selection, data = data)
  d <- model.response(smf)
  Zs <- model.matrix(selection, data = smf)
  probit <- suppressWarnings(stats::glm.fit(
    Zs, d, family = stats::binomial(link = "probit")
  ))
  gamma <- setNames(probit$coefficients, colnames(Zs))
  selected <- d == 1

  # Step 2: frontier on the selected rows
  f_base <- if (inherits(formula, "Formula")) formula(formula, rhs = 1) else formula
  dsel <- data[selected, , drop = FALSE]
  mf <- model.frame(f_base, data = dsel)
  y <- sign * model.response(mf)
  X <- model.matrix(f_base, data = mf)
  zg <- as.numeric(Zs[selected, , drop = FALSE] %*% gamma)
  k <- ncol(X)

  gh <- .gauss_hermite(100L)   # as sfaR::sfaselectioncross
  z <- sqrt(2) * gh$nodes
  wts <- gh$weights / sqrt(pi)
  # u at the quadrature nodes for u = sigma_u * q, q the standardised
  # inefficiency: |W| (half-normal) or -log Phi(-W) (exponential)
  q <- if (dist == "hnormal") abs(z) else -pnorm(-z, log.p = TRUE)

  node_terms <- function(params) {
    beta <- params[seq_len(k)]
    sv <- exp(params[k + 1L])
    su <- exp(params[k + 2L])
    rho <- tanh(params[k + 3L])
    eps <- as.numeric(y - X %*% beta)
    U <- outer(rep(1, length(y)), su * q)
    V <- eps + U
    log_f <- dnorm(V / sv, log = TRUE) - log(sv) +
      pnorm((zg + rho * V / sv) / sqrt(1 - rho^2), log.p = TRUE)
    list(U = U, f = sweep(exp(log_f), 2, wts, "*"))
  }
  loglik <- function(params) {
    ll <- sum(log(rowSums(node_terms(params)$f)))
    if (is.finite(ll)) ll else -1e20
  }

  ols <- lm.fit(X, y)
  s0 <- log(max(sqrt(sum(ols$residuals^2) / (length(y) - k)) * 0.5, 0.05))
  start <- c(ols$coefficients, log_sigma_v = s0, log_sigma_u = s0,
             atanh_rho = 0)
  ctrl <- list(fnscale = -1, maxit = 5000, reltol = 1e-10)
  ctrl[names(control)] <- control
  opt <- tryCatch(
    optim(start, loglik, method = "BFGS", control = ctrl, hessian = TRUE),
    error = function(e) {
      optim(start, loglik, method = "Nelder-Mead",
            control = list(fnscale = -1, maxit = 10000), hessian = TRUE)
    }
  )
  if (opt$convergence != 0L) {
    warning("Sample-selection SFA did not converge (code ",
            opt$convergence, ").", call. = FALSE)
  }

  # Posterior expectations over the quadrature nodes
  nt <- node_terms(opt$par)
  L <- rowSums(nt$f)
  te_bc88 <- rowSums(exp(-nt$U) * nt$f) / L
  te_jlms <- exp(-rowSums(nt$U * nt$f) / L)

  beta <- setNames(opt$par[seq_len(k)], colnames(X))
  res <- list(
    coefficients = beta,
    sigma_v = exp(opt$par[["log_sigma_v"]]),
    sigma_u = exp(opt$par[["log_sigma_u"]]),
    rho = tanh(opt$par[["atanh_rho"]]),
    gamma = gamma,
    logLik = opt$value,
    efficiency = if (estimator == "bc88") te_bc88 else te_jlms,
    efficiency_bc88 = te_bc88,
    efficiency_jlms = te_jlms,
    estimator = estimator,
    fitted = as.numeric(X %*% beta),
    residuals = as.numeric(y - X %*% beta),
    hessian = opt$hessian,
    convergence = opt$convergence,
    probit_converged = isTRUE(probit$converged),
    dist = dist,
    nobs = length(y),
    n_total = length(d),
    X = X,
    y = y,
    all_params = opt$par,
    selection = TRUE
  )
  res <- .unflip_frontier(res, sign)
  # corr(w, v) changes sign with v on the production-equivalent scale
  if (sign == -1) {
    res$rho <- -res$rho
    j <- k + 3L
    res$all_params[j] <- -res$all_params[j]
    res$hessian[j, ] <- -res$hessian[j, ]
    res$hessian[, j] <- -res$hessian[, j]
  }
  res
}


# Gauss-Hermite nodes and weights for weight function exp(-x^2)
# (Golub-Welsch)
.gauss_hermite <- function(n) {
  i <- seq_len(n - 1L)
  J <- matrix(0, n, n)
  J[cbind(i, i + 1L)] <- J[cbind(i + 1L, i)] <- sqrt(i / 2)
  e <- eigen(J, symmetric = TRUE)
  list(nodes = e$values, weights = sqrt(pi) * e$vectors[1, ]^2)
}
