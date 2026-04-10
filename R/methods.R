#' Print a Metafrontier Object
#'
#' @param x a \code{"metafrontier"} object.
#' @param ... additional arguments (currently unused).
#' @return Invisibly returns \code{x}.
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
#' @export
summary.metafrontier <- function(object, ...) {

  # Group coefficient tables with SEs
  group_tables <- list()
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
  if (which == "meta") {
    object$meta_coef
  } else {
    object$group_coef
  }
}


#' @export
vcov.metafrontier <- function(object, ...) {
  if (is.null(object$meta_vcov)) {
    warning("Variance-covariance matrix not available for ",
            "deterministic metafrontier. Use meta_type = 'stochastic'.",
            call. = FALSE)
    return(NULL)
  }
  k <- length(object$meta_coef)
  v <- object$meta_vcov[seq_len(k), seq_len(k)]
  rownames(v) <- colnames(v) <- names(object$meta_coef)
  v
}


#' @export
logLik.metafrontier <- function(object, ...) {
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
  object$meta_frontier
}


#' @export
residuals.metafrontier <- function(object, ...) {
  if (is.null(object$meta_frontier)) return(NULL)
  y <- model.response(model.frame(object$formula, data = object$data))
  y - object$meta_frontier
}
