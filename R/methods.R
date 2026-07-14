#' Print a Metafrontier Object
#'
#' Prints a compact overview of a fitted metafrontier model: the
#' estimation method and metafrontier type, the efficiency estimator
#' and identification objective (where applicable), the groups and
#' their sample sizes, group log-likelihoods, mean technology gap
#' ratio by group, and a one-line convergence status.
#'
#' @param x a \code{"metafrontier"} object.
#' @param ... additional arguments (currently unused).
#' @return Invisibly returns \code{x}.
#' @examples
#' sim <- simulate_metafrontier(n_groups = 2, n_per_group = 50, seed = 42)
#' fit <- metafrontier(log_y ~ log_x1 + log_x2, data = sim$data, group = "group")
#' print(fit)
#' @export
print.metafrontier <- function(x, ...) {
  cat("\nMetafrontier Model\n")
  cat("------------------\n")
  cat("Method:          ", x$method, "\n")
  cat("Metafrontier:    ", x$meta_type, "\n")
  if (!is.null(x$estimator)) {
    cat("Estimator:       ", x$estimator, "\n")
  }
  if (!is.null(x$objective)) {
    cat("Objective:       ", x$objective, "\n")
  }
  cat("Groups:          ", paste(x$groups, collapse = ", "), "\n")
  cat("Total obs:       ", x$nobs["total"], "\n")

  for (g in x$groups) {
    cat("  ", g, ": ", x$nobs[g], " obs\n", sep = "")
  }

  if (!is.null(x$logLik_groups)) {
    cat("\nGroup log-likelihoods:\n")
    for (g in x$groups) {
      cat("  ", g, ": ", format(x$logLik_groups[g], digits = 5), "\n",
          sep = "")
    }
  }

  cat("\nMean TGR by group:\n")
  for (g in x$groups) {
    idx <- x$group_vec == g
    cat("  ", g, ": ", format(mean(x$tgr[idx]), digits = 4), "\n", sep = "")
  }

  conv <- tryCatch(.convergence_table(x), error = function(e) NULL)
  if (!is.null(conv)) {
    bad <- conv$stage[which(!conv$converged)]
    if (length(bad) == 0L) {
      cat("\nConvergence: OK\n")
    } else {
      cat("\nConvergence: WARNING (", paste(bad, collapse = ", "), ")\n",
          sep = "")
    }
  }

  invisible(x)
}


#' Summary of a Metafrontier Model
#'
#' Computes group-level summaries of technical efficiency (TE),
#' technology gap ratio (TGR), and metafrontier efficiency (TE*),
#' full coefficient tables for each group frontier (including
#' variance parameters and, for BC92 panels, \code{eta}, all with
#' standard errors where a Hessian is available), the metafrontier
#' coefficient table (with Murphy-Topel corrected standard errors
#' where applicable), and a per-stage convergence table.
#'
#' @param object a \code{"metafrontier"} object.
#' @param ... additional arguments (currently unused).
#' @return An object of class \code{"summary.metafrontier"}: a list
#'   with components
#'   \describe{
#'     \item{call}{the matched call of the original fit}
#'     \item{method}{estimation method (\code{"sfa"} or \code{"dea"})}
#'     \item{meta_type}{metafrontier type (\code{"deterministic"} or
#'       \code{"stochastic"})}
#'     \item{groups}{character vector of group labels}
#'     \item{nobs}{named vector of observation counts (total and per
#'       group)}
#'     \item{group_tables}{named list of coefficient matrices, one per
#'       group, with columns \code{Estimate}, \code{Std. Error},
#'       \code{z value}, and \code{Pr(>|z|)} where standard errors are
#'       available (empty list for DEA fits)}
#'     \item{meta_table}{metafrontier coefficient matrix in the same
#'       format, or \code{NULL} for DEA fits}
#'     \item{tgr_summary}{data frame of TGR statistics by group, as
#'       returned by \code{\link{tgr_summary}}}
#'     \item{efficiency_summary}{data frame with mean TE, mean TGR,
#'       and mean TE* by group}
#'     \item{logLik_groups}{named vector of group log-likelihoods, or
#'       \code{NULL}}
#'     \item{meta_logLik}{Stage 2 log-likelihood of the stochastic
#'       metafrontier, or \code{NULL}}
#'     \item{convergence}{data frame with columns \code{stage},
#'       \code{code}, and \code{converged} recording the optimiser
#'       status of each estimation stage, or \code{NULL} if
#'       unavailable; see \code{\link{check_convergence}}}
#'   }
#' @examples
#' sim <- simulate_metafrontier(n_groups = 2, n_per_group = 50, seed = 42)
#' fit <- metafrontier(log_y ~ log_x1 + log_x2, data = sim$data,
#'                     group = "group", meta_type = "stochastic")
#' s <- summary(fit)
#' print(s)
#' @export
summary.metafrontier <- function(object, ...) {

  # Group coefficient tables with SEs (SFA only). Internally fitted
  # groups report the full MLE parameter vector (raw parameterisation:
  # log_sigma_v, log_sigma_u, mu, eta, ...) so that auxiliary
  # parameters such as eta (BC92) are visible with standard errors.
  group_tables <- list()
  if (object$method != "dea") {
    for (g in object$groups) {
      gm <- object$group_models[[g]]
      if (!is.null(gm$all_params) && !is.null(gm$hessian)) {
        est <- gm$all_params
        vcov_g <- tryCatch(
          solve(-gm$hessian),
          error = function(e) NULL
        )
        if (!is.null(vcov_g) && nrow(vcov_g) == length(est)) {
          se <- sqrt(pmax(diag(vcov_g), 0))
        } else {
          se <- rep(NA_real_, length(est))
        }
        zval <- est / se
        pval <- 2 * pnorm(-abs(zval))
        group_tables[[g]] <- cbind(
          Estimate = est,
          `Std. Error` = se,
          `z value` = zval,
          `Pr(>|z|)` = pval
        )
      } else if (!is.null(gm$hessian)) {
        vcov_g <- tryCatch(
          solve(-gm$hessian),
          error = function(e) NULL
        )
        if (!is.null(vcov_g)) {
          k <- length(gm$coefficients)
          se_all <- sqrt(pmax(diag(vcov_g), 0))
          se <- se_all[seq_len(k)]
          beta <- gm$coefficients
          zval <- beta / se
          pval <- 2 * pnorm(-abs(zval))
          group_tables[[g]] <- cbind(
            Estimate = beta,
            `Std. Error` = se,
            `z value` = zval,
            `Pr(>|z|)` = pval
          )
        } else {
          group_tables[[g]] <- cbind(Estimate = gm$coefficients)
        }
      } else {
        group_tables[[g]] <- cbind(Estimate = gm$coefficients)
      }
    }
  }

  # Metafrontier coefficient table
  if (!is.null(object$meta_coef)) {
    if (!is.null(object$meta_vcov)) {
      k <- length(object$meta_coef)
      se <- sqrt(pmax(diag(object$meta_vcov)[seq_len(k)], 0))
      zval <- object$meta_coef / se
      pval <- 2 * pnorm(-abs(zval))
      meta_table <- cbind(
        Estimate = object$meta_coef,
        `Std. Error` = se,
        `z value` = zval,
        `Pr(>|z|)` = pval
      )
    } else {
      meta_table <- cbind(Estimate = object$meta_coef)
    }
  } else {
    meta_table <- NULL
  }

  # TGR summary
  tgr_tab <- tgr_summary(object)

  # Efficiency summary
  eff_tab <- do.call(rbind, lapply(object$groups, function(g) {
    idx <- object$group_vec == g
    data.frame(
      Group = g,
      Mean_TE = mean(object$te_group[idx]),
      Mean_TGR = mean(object$tgr[idx]),
      Mean_TE_star = mean(object$te_meta[idx]),
      stringsAsFactors = FALSE
    )
  }))

  # Per-stage convergence status
  conv_tab <- tryCatch(.convergence_table(object), error = function(e) NULL)
  if (!is.null(conv_tab)) {
    conv_tab <- conv_tab[, c("stage", "code", "converged")]
  }

  out <- list(
    call = object$call,
    method = object$method,
    meta_type = object$meta_type,
    groups = object$groups,
    nobs = object$nobs,
    group_tables = group_tables,
    meta_table = meta_table,
    tgr_summary = tgr_tab,
    efficiency_summary = eff_tab,
    logLik_groups = object$logLik_groups,
    meta_logLik = object$meta_logLik,
    convergence = conv_tab
  )
  class(out) <- "summary.metafrontier"
  out
}


#' @export
print.summary.metafrontier <- function(x, digits = 4, ...) {
  cat("\nMetafrontier Model Summary\n")
  cat("==========================\n\n")
  cat("Call:\n")
  print(x$call)
  cat("\nMethod:       ", x$method, "\n")
  cat("Metafrontier: ", x$meta_type, "\n\n")

  # Group-specific frontiers
  if (length(x$group_tables) > 0) {
    for (g in x$groups) {
      cat("--- Group: ", g, " (n = ", x$nobs[g], ") ---\n", sep = "")
      if (ncol(x$group_tables[[g]]) >= 4) {
        printCoefmat(x$group_tables[[g]], digits = digits,
                     P.values = TRUE, has.Pvalue = TRUE)
      } else {
        print(round(x$group_tables[[g]], digits))
      }
      if (!is.null(x$logLik_groups)) {
        cat("Log-likelihood:", format(x$logLik_groups[g], digits = 5), "\n")
      }
      cat("\n")
    }
  } else {
    for (g in x$groups) {
      cat("--- Group: ", g, " (n = ", x$nobs[g], ") ---\n", sep = "")
      cat("  DEA efficiency (nonparametric)\n\n")
    }
  }

  # Metafrontier
  if (!is.null(x$meta_table)) {
    cat("--- Metafrontier ---\n")
    if (ncol(x$meta_table) >= 4) {
      printCoefmat(x$meta_table, digits = digits,
                   P.values = TRUE, has.Pvalue = TRUE)
    } else {
      print(round(x$meta_table, digits))
    }
    if (!is.null(x$meta_logLik)) {
      cat("Log-likelihood:", format(x$meta_logLik, digits = 5), "\n")
    }
    cat("\n")
  }

  # Efficiency decomposition
  cat("--- Efficiency Decomposition ---\n")
  eff_print <- x$efficiency_summary
  num_cols <- sapply(eff_print, is.numeric)
  eff_print[num_cols] <- lapply(eff_print[num_cols], round, digits = digits)
  print(eff_print, row.names = FALSE)

  # TGR summary
  cat("\n--- Technology Gap Ratio Summary ---\n")
  tgr_print <- x$tgr_summary
  tgr_print[, -1] <- round(tgr_print[, -1], digits)
  print(tgr_print, row.names = FALSE)

  # Convergence status
  if (!is.null(x$convergence)) {
    cat("\n--- Convergence ---\n")
    bad <- which(!x$convergence$converged)
    if (length(bad) == 0L) {
      cat("All estimation stages converged.\n")
    } else {
      print(x$convergence[bad, , drop = FALSE], row.names = FALSE)
      cat("See ?check_convergence.\n")
    }
  }

  cat("\n")
  invisible(x)
}


#' Extract Coefficients from a Metafrontier Model
#'
#' @param object a \code{"metafrontier"} object.
#' @param which character. \code{"meta"} (default) returns the
#'   metafrontier coefficients; \code{"group"} returns a named list
#'   of group-specific coefficient vectors.
#' @param extraPar logical. If \code{TRUE}, auxiliary parameters are
#'   included alongside the frontier coefficients. For
#'   \code{which = "group"} the variance parameters are
#'   back-transformed to their natural scale (\code{sigmaV},
#'   \code{sigmaU}), \code{mu} and \code{eta} are kept as estimated,
#'   and heteroscedastic Z coefficients are labelled with their
#'   column names. For \code{which = "meta"} the Stage 2 variance
#'   parameters are appended for stochastic metafrontiers.
#' @param ... additional arguments (currently unused).
#'
#' @return A named numeric vector (\code{which = "meta"}) or a named
#'   list of numeric vectors (\code{which = "group"}).
#'
#' @examples
#' sim <- simulate_metafrontier(n_groups = 2, n_per_group = 50, seed = 42)
#' fit <- metafrontier(log_y ~ log_x1 + log_x2, data = sim$data,
#'                     group = "group", meta_type = "stochastic")
#' coef(fit)
#' coef(fit, extraPar = TRUE)
#' coef(fit, which = "group", extraPar = TRUE)
#'
#' @export
coef.metafrontier <- function(object,
                              which = c("meta", "group"),
                              extraPar = FALSE,
                              ...) {
  which <- match.arg(which)
  if (object$method == "dea") {
    stop("coef() is not available for DEA-based metafrontiers ",
         "(nonparametric model).", call. = FALSE)
  }
  if (which == "meta") {
    if (!extraPar) return(object$meta_coef)
    if (identical(object$meta_type, "stochastic") &&
        !is.null(object$meta_opt)) {
      pars <- object$meta_opt$par
      k <- length(object$meta_coef)
      extra <- .relabel_extra_params(pars[-seq_len(k)])
      return(c(object$meta_coef, extra))
    }
    message("No auxiliary parameters exist for this metafrontier; ",
            "returning the frontier coefficients only.")
    return(object$meta_coef)
  }

  # which == "group"
  if (!extraPar) return(object$group_coef)

  missing_extra <- vapply(object$group_models,
                          function(gm) is.null(gm$all_params),
                          logical(1))
  if (any(missing_extra)) {
    warning("Extra parameters are unavailable for externally fitted ",
            "groups (",
            paste(object$groups[missing_extra], collapse = ", "),
            "); returning their frontier coefficients only.",
            call. = FALSE)
  }

  out <- lapply(object$groups, function(g) {
    gm <- object$group_models[[g]]
    if (is.null(gm$all_params)) return(gm$coefficients)
    params <- gm$all_params
    if (!is.null(gm$Z)) {
      idx_z <- grep("^d_[0-9]+$", names(params))
      if (length(idx_z) == ncol(gm$Z)) {
        names(params)[idx_z] <- colnames(gm$Z)
      }
    }
    .relabel_extra_params(params)
  })
  names(out) <- object$groups
  out
}


#' Back-transform and rename raw MLE variance parameters
#'
#' Replaces \code{log_sigma_v}/\code{log_sigma_u} with
#' \code{sigmaV}/\code{sigmaU} on the natural scale; all other
#' entries (frontier coefficients, mu, eta, Z coefficients) pass
#' through unchanged.
#'
#' @keywords internal
#' @noRd
.relabel_extra_params <- function(params) {
  nm <- names(params)
  sv <- nm == "log_sigma_v"
  su <- nm == "log_sigma_u"
  params[sv] <- exp(params[sv])
  params[su] <- exp(params[su])
  nm[sv] <- "sigmaV"
  nm[su] <- "sigmaU"
  names(params) <- nm
  params
}


#' Variance-Covariance Matrix for Metafrontier Coefficients
#'
#' Returns the variance-covariance matrix of the Stage 2
#' (metafrontier) coefficients, or, with \code{which = "group"}, the
#' per-group matrices from the Stage 1 maximum likelihood fits.
#' \code{NULL} is returned when no Stage 2 Hessian exists: the
#' deterministic metafrontier is fitted by LP/QP optimisation and has
#' no sampling variance in this framework, so \code{which = "meta"}
#' returns \code{NULL} with a warning; with \code{which = "group"},
#' list entries are \code{NULL} for groups without a stored Hessian
#' (e.g. externally fitted models). DEA-based metafrontiers are
#' nonparametric and \code{vcov()} signals an error; use
#' \code{\link{boot_tgr}} for inference instead.
#'
#' @param object a \code{"metafrontier"} object.
#' @param correction character. \code{"none"} (default) returns the
#'   Stage 2 variance-covariance matrix. \code{"murphy-topel"} applies
#'   the Murphy and Topel (1985) correction for first-stage estimation
#'   uncertainty (the generated-regressor problem). Only available for
#'   stochastic metafrontiers.
#' @param which character. \code{"meta"} (default) returns the
#'   metafrontier (Stage 2) variance-covariance matrix;
#'   \code{"group"} returns a named list with one full
#'   variance-covariance matrix per group (from the inverse negative
#'   Hessian of the group MLE), with \code{NULL} entries for groups
#'   without a stored Hessian.
#' @param extraPar logical. If \code{TRUE} and \code{which = "meta"},
#'   the full Stage 2 matrix is returned, including the rows and
#'   columns for the auxiliary parameters (raw MLE parameterisation,
#'   e.g. \code{log_sigma_v}); the default returns only the block for
#'   the frontier coefficients.
#' @param ... additional arguments (currently unused).
#'
#' @return A variance-covariance matrix (\code{which = "meta"}), a
#'   named list of matrices (\code{which = "group"}), or \code{NULL}
#'   if unavailable.
#'
#' @references Murphy, K.M. and Topel, R.H. (1985). Estimation and
#'   inference in two-step econometric models. \emph{Journal of
#'   Business & Economic Statistics}, 3(4), 370--379.
#'
#' @examples
#' sim <- simulate_metafrontier(n_groups = 2, n_per_group = 50, seed = 42)
#' fit <- metafrontier(log_y ~ log_x1 + log_x2, data = sim$data,
#'                     group = "group", meta_type = "stochastic")
#' vcov(fit)
#' \donttest{
#' vcov(fit, correction = "murphy-topel")
#' }
#'
#' @export
vcov.metafrontier <- function(object,
                              correction = c("none", "murphy-topel"),
                              which = c("meta", "group"),
                              extraPar = FALSE,
                              ...) {
  correction <- match.arg(correction)
  which <- match.arg(which)

  if (object$method == "dea") {
    stop("vcov() is not available for DEA-based metafrontiers ",
         "(nonparametric model).", call. = FALSE)
  }

  if (which == "group") {
    out <- lapply(object$group_models, function(gm) {
      if (is.null(gm$hessian)) return(NULL)
      v <- tryCatch(solve(-gm$hessian), error = function(e) NULL)
      if (!is.null(v) && !is.null(gm$all_params) &&
          nrow(v) == length(gm$all_params)) {
        rownames(v) <- colnames(v) <- names(gm$all_params)
      }
      v
    })
    names(out) <- object$groups
    return(out)
  }

  k <- length(object$meta_coef)

  if (correction == "murphy-topel") {
    if (object$meta_type != "stochastic") {
      stop("Murphy-Topel correction requires a stochastic metafrontier ",
           "(meta_type = 'stochastic').", call. = FALSE)
    }
    if (!requireNamespace("numDeriv", quietly = TRUE)) {
      stop("Package 'numDeriv' is required for the Murphy-Topel ",
           "correction. Install it with install.packages('numDeriv').",
           call. = FALSE)
    }
    v <- .murphy_topel_correction(object)
    if (is.null(v)) {
      warning("Murphy-Topel correction failed; returning uncorrected ",
              "variance-covariance matrix.", call. = FALSE)
      v <- object$meta_vcov
    }
  } else {
    if (is.null(object$meta_vcov)) {
      warning("Variance-covariance matrix not available for ",
              "deterministic metafrontier. Use meta_type = 'stochastic'.",
              call. = FALSE)
      return(NULL)
    }
    v <- object$meta_vcov
  }

  if (extraPar) {
    if (!is.null(object$meta_opt) &&
        length(object$meta_opt$par) == nrow(v)) {
      rownames(v) <- colnames(v) <- names(object$meta_opt$par)
    }
    return(v)
  }

  v <- v[seq_len(k), seq_len(k), drop = FALSE]
  rownames(v) <- colnames(v) <- names(object$meta_coef)
  v
}


#' @export
logLik.metafrontier <- function(object, ...) {
  if (object$method == "dea") {
    stop("logLik() is not available for DEA-based metafrontiers ",
         "(nonparametric model).", call. = FALSE)
  }
  if (!is.null(object$meta_logLik)) {
    val <- object$meta_logLik
  } else {
    val <- sum(object$logLik_groups)
  }
  attr(val, "df") <- length(object$meta_coef)
  attr(val, "nobs") <- object$nobs["total"]
  class(val) <- "logLik"
  val
}


#' @export
nobs.metafrontier <- function(object, ...) {
  unname(object$nobs["total"])
}


#' @export
fitted.metafrontier <- function(object, ...) {
  if (is.null(object$meta_frontier)) {
    stop("fitted() is not available for DEA-based metafrontiers. ",
         "Use efficiencies() to extract DEA efficiency scores.",
         call. = FALSE)
  }

  object$meta_frontier
}


#' @export
residuals.metafrontier <- function(object, ...) {
  if (object$method == "dea") {
    stop("residuals() is not available for DEA-based metafrontiers.",
         call. = FALSE)
  }
  if (is.null(object$meta_frontier)) return(NULL)
  y <- model.response(model.frame(object$formula, data = object$data))
  y - object$meta_frontier
}


#' Confidence Intervals for Metafrontier Coefficients
#'
#' Computes Wald-type confidence intervals for the metafrontier
#' coefficients using the variance-covariance matrix from the
#' stochastic metafrontier.
#'
#' @param object a \code{"metafrontier"} object.
#' @param parm a character or integer vector of parameter names or
#'   indices. If missing, all frontier coefficients are used.
#' @param level the confidence level (default 0.95).
#' @param correction character. Passed to \code{\link{vcov.metafrontier}}.
#'   \code{"none"} (default) uses the Stage 2 variance-covariance
#'   matrix; \code{"murphy-topel"} applies the Murphy and Topel (1985)
#'   correction for first-stage estimation uncertainty.
#' @param ... additional arguments (currently unused).
#'
#' @return A matrix with columns for the lower and upper bounds.
#'
#' @examples
#' sim <- simulate_metafrontier(n_groups = 2, n_per_group = 50, seed = 42)
#' fit <- metafrontier(log_y ~ log_x1 + log_x2, data = sim$data,
#'                     group = "group", meta_type = "stochastic")
#' confint(fit)
#' # With Murphy-Topel correction (requires numDeriv):
#' \donttest{
#' confint(fit, correction = "murphy-topel")
#' }
#'
#' @export
confint.metafrontier <- function(object, parm, level = 0.95,
                                 correction = c("none", "murphy-topel"),
                                 ...) {
  correction <- match.arg(correction)
  v <- vcov(object, correction = correction)
  if (is.null(v)) {
    stop("Confidence intervals require a stochastic metafrontier ",
         "(meta_type = 'stochastic').", call. = FALSE)
  }

  cf <- object$meta_coef
  k <- length(cf)
  se <- sqrt(pmax(diag(v), 0))

  a <- (1 - level) / 2
  z <- qnorm(1 - a)

  ci <- cbind(cf - z * se, cf + z * se)
  rownames(ci) <- names(cf)
  pct <- paste(format(100 * c(a, 1 - a), trim = TRUE, digits = 3), "%")
  colnames(ci) <- pct

  if (!missing(parm)) {
    ci <- ci[parm, , drop = FALSE]
  }
  ci
}


#' Predict Frontier Values from a Metafrontier Model
#'
#' Computes predicted frontier values at given input levels using
#' either the metafrontier or a group-specific frontier.
#'
#' @param object a \code{"metafrontier"} object.
#' @param newdata an optional data frame of new inputs. If omitted,
#'   the training data predictions are returned.
#' @param type character. \code{"meta"} (default) for metafrontier
#'   predictions, or \code{"group"} for group-specific frontier
#'   predictions (requires a \code{group} column in \code{newdata}).
#' @param ... additional arguments (currently unused).
#'
#' @return A numeric vector of predicted frontier values.
#'
#' @examples
#' sim <- simulate_metafrontier(n_groups = 2, n_per_group = 50, seed = 42)
#' fit <- metafrontier(log_y ~ log_x1 + log_x2, data = sim$data,
#'                     group = "group", meta_type = "stochastic")
#' pred <- predict(fit)
#' # Out-of-sample prediction:
#' newdata <- data.frame(log_x1 = c(1, 2), log_x2 = c(1.5, 2.5))
#' predict(fit, newdata = newdata)
#'
#' @export
predict.metafrontier <- function(object, newdata = NULL,
                                 type = c("meta", "group"), ...) {
  type <- match.arg(type)

  if (is.null(newdata)) {
    if (type == "meta") {
      return(object$meta_frontier)
    } else {
      return(object$group_frontier)
    }
  }

  if (object$method == "dea") {
    stop("predict() with newdata is only supported for SFA metafrontiers.",
         call. = FALSE)
  }

  # Build design matrix from newdata (no response needed)
  f <- object$formula
  if (inherits(f, "Formula")) {
    f_base <- formula(f, rhs = 1)
  } else {
    f_base <- f
  }
  tt <- delete.response(terms(f_base))
  mf_new <- model.frame(tt, data = newdata, na.action = na.omit)
  X_new <- model.matrix(tt, data = mf_new)

  if (type == "meta") {
    if (is.null(object$meta_coef)) {
      stop("Metafrontier coefficients not available (DEA model?).",
           call. = FALSE)
    }
    return(as.numeric(X_new %*% object$meta_coef))
  }

  # Group predictions
  if (!"group" %in% names(newdata) && is.null(object$group_vec)) {
    stop("'newdata' must contain a 'group' column for type = 'group'.",
         call. = FALSE)
  }

  group_col <- if ("group" %in% names(newdata)) newdata$group
               else object$group_vec
  preds <- numeric(nrow(X_new))
  for (g in object$groups) {
    idx <- which(group_col == g)
    if (length(idx) == 0) next
    beta_g <- object$group_coef[[g]]
    preds[idx] <- as.numeric(X_new[idx, , drop = FALSE] %*% beta_g)
  }
  preds
}
