#' Print a Metafrontier Object
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

  invisible(x)
}


#' Summary of a Metafrontier Model
#'
#' @param object a \code{"metafrontier"} object.
#' @param ... additional arguments (currently unused).
#' @return An object of class \code{"summary.metafrontier"}.
#' @examples
#' sim <- simulate_metafrontier(n_groups = 2, n_per_group = 50, seed = 42)
#' fit <- metafrontier(log_y ~ log_x1 + log_x2, data = sim$data,
#'                     group = "group", meta_type = "stochastic")
#' s <- summary(fit)
#' print(s)
#' @export
summary.metafrontier <- function(object, ...) {

  # Group coefficient tables with SEs (SFA only)
  group_tables <- list()
  if (object$method != "dea") {
    for (g in object$groups) {
      gm <- object$group_models[[g]]
      if (!is.null(gm$hessian)) {
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
    meta_logLik = object$meta_logLik
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

  cat("\n")
  invisible(x)
}


#' @export
coef.metafrontier <- function(object,
                              which = c("meta", "group"),
                              ...) {
  which <- match.arg(which)
  if (object$method == "dea") {
    stop("coef() is not available for DEA-based metafrontiers ",
         "(nonparametric model).", call. = FALSE)
  }
  if (which == "meta") {
    object$meta_coef
  } else {
    object$group_coef
  }
}


#' Variance-Covariance Matrix for Metafrontier Coefficients
#'
#' @param object a \code{"metafrontier"} object.
#' @param correction character. \code{"none"} (default) returns the
#'   Stage 2 variance-covariance matrix. \code{"murphy-topel"} applies
#'   the Murphy and Topel (1985) correction for first-stage estimation
#'   uncertainty (the generated-regressor problem). Only available for
#'   stochastic metafrontiers.
#' @param ... additional arguments (currently unused).
#'
#' @return A variance-covariance matrix, or \code{NULL} if unavailable.
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
                              ...) {
  correction <- match.arg(correction)

  if (object$method == "dea") {
    stop("vcov() is not available for DEA-based metafrontiers ",
         "(nonparametric model).", call. = FALSE)
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
      v <- object$meta_vcov[seq_len(k), seq_len(k)]
    } else {
      v <- v[seq_len(k), seq_len(k)]
    }
  } else {
    if (is.null(object$meta_vcov)) {
      warning("Variance-covariance matrix not available for ",
              "deterministic metafrontier. Use meta_type = 'stochastic'.",
              call. = FALSE)
      return(NULL)
    }
    v <- object$meta_vcov[seq_len(k), seq_len(k)]
  }

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
