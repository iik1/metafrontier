#' Latent Class Metafrontier
#'
#' Estimates a metafrontier model where group membership is unobserved,
#' using an EM algorithm to jointly estimate class membership
#' probabilities, class-specific frontier parameters, and the
#' metafrontier.
#'
#' @param formula a \code{Formula} object (y ~ x1 + x2).
#' @param data a data frame.
#' @param n_classes integer. Number of latent classes (default 2).
#' @param dist distribution of the inefficiency term.
#' @param meta_type metafrontier type for Stage 2.
#' @param n_starts integer. Number of random initializations
#'   (default 10). The best is selected by log-likelihood.
#' @param max_iter integer. Maximum EM iterations (default 200).
#' @param tol numeric. Convergence tolerance on marginal LL (default 1e-6).
#' @param seed optional integer seed.
#' @param control list of control parameters for the optimiser.
#' @param ... additional arguments.
#'
#' @return An object of class \code{"lc_metafrontier"} containing:
#'   \describe{
#'     \item{class_assignment}{MAP class assignment per observation}
#'     \item{posterior}{n x C matrix of posterior probabilities}
#'     \item{pi}{class mixing proportions}
#'     \item{class_params}{list of class-specific parameter vectors}
#'     \item{class_models}{list of class-specific model summaries}
#'     \item{metafrontier}{the fitted metafrontier object on MAP classes}
#'     \item{marginal_ll}{marginal log-likelihood at convergence}
#'     \item{BIC}{Bayesian Information Criterion}
#'     \item{n_classes}{number of classes}
#'     \item{n_iter}{number of EM iterations used}
#'   }
#'
#' @details
#' The EM algorithm iterates between:
#' \itemize{
#'   \item \strong{E-step}: compute posterior class membership probabilities
#'     for each observation using Bayes' rule
#'   \item \strong{M-step}: update class-specific frontier parameters via
#'     weighted MLE, and update class proportions
#' }
#' Multiple random starts (\code{n_starts}) are used to avoid local optima.
#' The run with the highest marginal log-likelihood is selected.
#' After convergence, observations are assigned to classes via MAP
#' (maximum a posteriori), and a standard metafrontier is fitted on the
#' MAP classes.
#'
#' @examples
#' \donttest{
#' sim <- simulate_metafrontier(n_groups = 2, n_per_group = 80, seed = 42)
#' lc <- latent_class_metafrontier(
#'   log_y ~ log_x1 + log_x2,
#'   data = sim$data, n_classes = 2, n_starts = 3, seed = 123
#' )
#' print(lc)
#' summary(lc)
#' coef(lc, which = "meta")
#' efficiencies(lc, type = "tgr")
#' }
#'
#' @export
latent_class_metafrontier <- function(formula, data,
                                      n_classes = 2,
                                      dist = c("hnormal", "tnormal",
                                                "exponential"),
                                      meta_type = c("deterministic",
                                                     "stochastic"),
                                      n_starts = 10,
                                      max_iter = 200,
                                      tol = 1e-6,
                                      seed = NULL,
                                      control = list(),
                                      ...) {

  dist <- match.arg(dist)
  meta_type <- match.arg(meta_type)

  if (!is.null(seed)) set.seed(seed)

  f <- Formula::Formula(formula)
  f_base <- formula(f, rhs = 1)
  mf <- model.frame(f_base, data = data, na.action = na.omit)
  y <- model.response(mf)
  X <- model.matrix(f_base, data = mf)
  n <- length(y)
  k <- ncol(X)

  # Obs-level LL function selector
  obs_ll_fn <- switch(dist,
    hnormal     = .loglik_hnormal_obs,
    tnormal     = .loglik_tnormal_obs,
    exponential = .loglik_exponential_obs
  )

  # Number of params per class: beta + log_sigma_v + log_sigma_u (+ mu for tnormal)
  n_params_class <- if (dist == "tnormal") k + 3L else k + 2L

  # --- Multiple starts ---
  best_ll <- -Inf
  best_result <- NULL

  for (start in seq_len(n_starts)) {
    result <- tryCatch(
      .em_one_start(y, X, n_classes, k, n_params_class, dist,
                    obs_ll_fn, max_iter, tol, control, start),
      error = function(e) NULL
    )
    if (!is.null(result) && result$marginal_ll > best_ll) {
      best_ll <- result$marginal_ll
      best_result <- result
    }
  }

  if (is.null(best_result)) {
    stop("All EM starts failed to converge.", call. = FALSE)
  }

  # --- MAP assignment and metafrontier ---
  class_assign <- apply(best_result$posterior, 1, which.max)

  # Sort classes by mean efficiency (avoid label switching)
  class_order <- order(sapply(seq_len(n_classes), function(c) {
    params_c <- best_result$class_params[[c]]
    mean(exp(-.loglik_to_u_hat(params_c, y, X, k, dist)))
  }), decreasing = TRUE)

  class_assign_sorted <- match(class_assign, class_order)
  class_names <- paste0("LC", seq_len(n_classes))
  data$lc_group <- class_names[class_assign_sorted]

  # Fit metafrontier on MAP classes
  meta_fit <- metafrontier(
    formula = formula,
    data = data,
    group = "lc_group",
    method = "sfa",
    meta_type = meta_type,
    dist = dist,
    control = control,
    ...
  )

  # BIC
  n_total_params <- n_classes * n_params_class + (n_classes - 1)
  bic <- -2 * best_ll + n_total_params * log(n)

  out <- list(
    class_assignment = class_names[class_assign_sorted],
    posterior = best_result$posterior[, class_order],
    pi = best_result$pi[class_order],
    class_params = best_result$class_params[class_order],
    metafrontier = meta_fit,
    marginal_ll = best_ll,
    BIC = bic,
    n_classes = n_classes,
    n_iter = best_result$n_iter,
    n = n
  )
  class(out) <- "lc_metafrontier"
  out
}


#' Select Number of Latent Classes via BIC
#'
#' @param formula formula.
#' @param data data frame.
#' @param n_classes_range integer vector of class counts to try.
#' @param ... additional arguments passed to
#'   \code{\link{latent_class_metafrontier}}.
#'
#' @return A data frame with columns \code{n_classes}, \code{BIC},
#'   and \code{marginal_ll}.
#'
#' @details
#' Fits latent class metafrontier models for each value in
#' \code{n_classes_range} and returns BIC values. The optimal
#' number of classes minimises BIC.
#'
#' @examples
#' \donttest{
#' sim <- simulate_metafrontier(n_groups = 2, n_per_group = 80, seed = 42)
#' bic_table <- select_n_classes(
#'   log_y ~ log_x1 + log_x2,
#'   data = sim$data, n_classes_range = 2:3,
#'   n_starts = 3, seed = 42
#' )
#' print(bic_table)
#' }
#'
#' @export
select_n_classes <- function(formula, data,
                             n_classes_range = 2:5,
                             ...) {
  results <- data.frame(
    n_classes = integer(0),
    BIC = numeric(0),
    marginal_ll = numeric(0)
  )

  for (nc in n_classes_range) {
    fit <- tryCatch(
      latent_class_metafrontier(formula, data, n_classes = nc, ...),
      error = function(e) NULL
    )
    if (!is.null(fit)) {
      results <- rbind(results, data.frame(
        n_classes = nc,
        BIC = fit$BIC,
        marginal_ll = fit$marginal_ll
      ))
    }
  }

  results
}


# ---------- Internal: Single EM start ----------

.em_one_start <- function(y, X, n_classes, k, n_params_class, dist,
                          obs_ll_fn, max_iter, tol, control, start_id) {
  n <- length(y)

  # Initialize via k-means + perturbation
  km_features <- cbind(y, X[, -1, drop = FALSE])
  km <- kmeans(km_features, centers = n_classes, nstart = 1)
  tau <- matrix(0, n, n_classes)
  for (i in seq_len(n)) {
    tau[i, km$cluster[i]] <- 0.8
    tau[i, -km$cluster[i]] <- 0.2 / (n_classes - 1)
  }
  # Add randomness for different starts
  tau <- tau + matrix(runif(n * n_classes, 0, 0.1), n, n_classes)
  tau <- tau / rowSums(tau)

  pi_c <- colMeans(tau)

  # Initialize class parameters via weighted OLS
  class_params <- vector("list", n_classes)
  for (c in seq_len(n_classes)) {
    w <- tau[, c]
    beta_init <- tryCatch({
      wfit <- lm.wfit(X, y, w = w)
      wfit$coefficients
    }, error = function(e) lm.fit(X, y)$coefficients)
    resid <- y - X %*% beta_init
    sigma_ols <- sqrt(max(sum(w * resid^2) / sum(w), 0.01))

    if (dist == "tnormal") {
      class_params[[c]] <- c(
        as.numeric(beta_init),
        log_sigma_v = log(max(sigma_ols * 0.5, 0.05)),
        mu = 0,
        log_sigma_u = log(max(sigma_ols * 0.5, 0.05))
      )
    } else {
      class_params[[c]] <- c(
        as.numeric(beta_init),
        log_sigma_v = log(max(sigma_ols * 0.5, 0.05)),
        log_sigma_u = log(max(sigma_ols * 0.5, 0.05))
      )
    }
  }

  prev_ll <- -Inf
  ctrl <- list(fnscale = -1, maxit = 50, reltol = 1e-8)
  ctrl[names(control)] <- control

  for (iter in seq_len(max_iter)) {
    # --- E-step ---
    log_L <- matrix(NA_real_, n, n_classes)
    for (c in seq_len(n_classes)) {
      ll_obs <- obs_ll_fn(class_params[[c]], y, X)
      log_L[, c] <- ll_obs
    }

    # Posterior: tau_ic = pi_c * L_c / sum_c'
    log_tau <- sweep(log_L, 2, log(pi_c), "+")
    log_tau_max <- apply(log_tau, 1, max)
    log_tau <- log_tau - log_tau_max
    tau <- exp(log_tau)
    tau <- tau / rowSums(tau)

    # Marginal LL
    marginal_ll <- sum(log_tau_max + log(rowSums(exp(log_tau))))
    # Correction: we already subtracted max, so add it back
    # Actually: marginal_ll = sum(log(sum_c pi_c * L_c(y_i)))
    log_lik_mix <- log_tau_max + log(rowSums(exp(log_tau)))
    marginal_ll <- sum(log_lik_mix)

    # Convergence check
    if (iter > 1 && abs(marginal_ll - prev_ll) < tol * abs(prev_ll)) {
      break
    }
    prev_ll <- marginal_ll

    # --- M-step ---
    pi_c <- pmax(colMeans(tau), 0.01)
    pi_c <- pi_c / sum(pi_c)

    for (c in seq_len(n_classes)) {
      w <- tau[, c]
      if (sum(w) < 1) next

      # Weighted MLE
      weighted_ll <- function(params) {
        ll_obs <- obs_ll_fn(params, y, X)
        sum(w * ll_obs)
      }

      opt <- tryCatch(
        optim(class_params[[c]], weighted_ll,
              method = "BFGS", control = ctrl),
        error = function(e) NULL
      )
      if (!is.null(opt) && is.finite(opt$value)) {
        class_params[[c]] <- opt$par
      }
    }
  }

  list(
    posterior = tau,
    pi = pi_c,
    class_params = class_params,
    marginal_ll = marginal_ll,
    n_iter = iter
  )
}


# ---------- Helper: safe Mills ratio for latent class ----------

.safe_mills_lc <- function(z) {
  p <- pnorm(z)
  p <- pmax(p, .Machine$double.eps)
  dnorm(z) / p
}


# ---------- Helper: extract u_hat from params ----------

.loglik_to_u_hat <- function(params, y, X, k, dist) {
  beta <- params[seq_len(k)]
  sigma_v <- exp(params[k + 1])
  eps <- as.numeric(y - X %*% beta)

  if (dist == "exponential") {
    # Exponential JLMS conditional mean (Jondrow et al. 1982)
    sigma_u <- exp(params[k + 2])
    mu_star <- -eps - sigma_v^2 / sigma_u
    u_hat <- mu_star + sigma_v * .safe_mills_lc(mu_star / sigma_v)
  } else if (dist == "tnormal") {
    # Truncated-normal JLMS conditional mean
    mu <- params[k + 2]  # mu sits between log_sigma_v and log_sigma_u
    sigma_u <- exp(params[k + 3])
    sigma_sq <- sigma_v^2 + sigma_u^2
    mu_star <- (mu * sigma_v^2 - eps * sigma_u^2) / sigma_sq
    sigma_star <- sigma_v * sigma_u / sqrt(sigma_sq)
    u_hat <- mu_star + sigma_star * .safe_mills_lc(mu_star / sigma_star)
  } else {
    # Half-normal JLMS conditional mean (default)
    sigma_u <- exp(params[k + 2])
    sigma_sq <- sigma_v^2 + sigma_u^2
    mu_star <- -eps * sigma_u^2 / sigma_sq
    sigma_star <- sigma_v * sigma_u / sqrt(sigma_sq)
    u_hat <- mu_star + sigma_star * .safe_mills_lc(mu_star / sigma_star)
  }
  pmax(u_hat, 0)
}


# ---------- S3 Methods ----------

#' @export
print.lc_metafrontier <- function(x, digits = 4, ...) {
  cat("\nLatent Class Metafrontier\n")
  cat("========================\n")
  cat("Classes:       ", x$n_classes, "\n")
  cat("EM iterations: ", x$n_iter, "\n")
  cat("Marginal LL:   ", format(x$marginal_ll, digits = 6), "\n")
  cat("BIC:           ", format(x$BIC, digits = 6), "\n\n")

  cat("Class proportions:\n")
  for (c in seq_len(x$n_classes)) {
    cat("  LC", c, ": ", format(x$pi[c], digits = digits),
        " (n = ", sum(x$class_assignment == paste0("LC", c)), ")\n",
        sep = "")
  }
  cat("\n")
  invisible(x)
}


#' @export
summary.lc_metafrontier <- function(object, ...) {
  out <- list(
    n_classes = object$n_classes,
    n_iter = object$n_iter,
    marginal_ll = object$marginal_ll,
    BIC = object$BIC,
    pi = object$pi,
    class_sizes = table(object$class_assignment),
    meta_summary = summary(object$metafrontier)
  )
  class(out) <- "summary.lc_metafrontier"
  out
}


#' @export
print.summary.lc_metafrontier <- function(x, digits = 4, ...) {
  cat("\nLatent Class Metafrontier\n")
  cat("========================\n")
  cat("Classes:       ", x$n_classes, "\n")
  cat("EM iterations: ", x$n_iter, "\n")
  cat("Marginal LL:   ", format(x$marginal_ll, digits = 6), "\n")
  cat("BIC:           ", format(x$BIC, digits = 6), "\n\n")

  cat("Class proportions:\n")
  for (c in seq_len(x$n_classes)) {
    cat("  LC", c, ": ", format(x$pi[c], digits = digits),
        " (n = ", x$class_sizes[paste0("LC", c)], ")\n",
        sep = "")
  }
  cat("\n--- Metafrontier on MAP classes ---\n\n")
  print(x$meta_summary)
  invisible(x)
}


#' @export
coef.lc_metafrontier <- function(object,
                                 which = c("meta", "class"),
                                 ...) {
  which <- match.arg(which)
  if (which == "meta") {
    coef(object$metafrontier, which = "meta")
  } else {
    object$class_params
  }
}


#' @export
efficiencies.lc_metafrontier <- function(object,
                                         type = c("meta", "group",
                                                   "tgr"),
                                         ...) {
  type <- match.arg(type)
  efficiencies(object$metafrontier, type = type)
}
