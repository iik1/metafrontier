# Convergence diagnostics for metafrontier fits.


#' Internal: per-stage convergence table for a metafrontier object
#'
#' Builds one row per group-level stage plus one row for the
#' metafrontier stage. Shared by \code{check_convergence()},
#' \code{summary.metafrontier()}, and \code{print.metafrontier()}.
#'
#' @param object a \code{"metafrontier"} object.
#' @return A data frame with columns \code{stage}, \code{method},
#'   \code{code}, \code{converged}, and \code{note}.
#' @keywords internal
#' @noRd
.convergence_table <- function(object) {

  groups <- object$groups
  is_dea <- identical(object$method, "dea")
  n_stages <- length(groups) + 1L

  stage <- character(n_stages)
  method <- character(n_stages)
  code <- integer(n_stages)
  converged <- logical(n_stages)
  note <- character(n_stages)

  .dea_note <- function(n_na) {
    if (n_na > 0L) {
      paste0(n_na, " infeasible DEA programs (NA efficiency)")
    } else {
      ""
    }
  }

  for (i in seq_along(groups)) {
    g <- groups[i]
    gm <- object$group_models[[g]]
    stage[i] <- paste0("group: ", g)

    if (is_dea) {
      n_na <- sum(is.na(object$te_group[object$group_vec == g]))
      method[i] <- "DEA"
      code[i] <- NA_integer_
      converged[i] <- n_na == 0L
      note[i] <- .dea_note(n_na)
    } else if (is.null(gm$convergence)) {
      # Models-path import: the external fitter's convergence status
      # is not carried through the extractors.
      method[i] <- "external"
      code[i] <- NA_integer_
      converged[i] <- NA
      note[i] <- "fitted externally - convergence not tracked"
    } else {
      method[i] <- "MLE"
      code[i] <- as.integer(gm$convergence)
      converged[i] <- code[i] == 0L
      note[i] <- ""
    }
  }

  j <- n_stages
  stage[j] <- "metafrontier"

  if (is_dea) {
    n_na <- sum(is.na(object$te_meta))
    method[j] <- "DEA"
    code[j] <- NA_integer_
    converged[j] <- n_na == 0L
    note[j] <- .dea_note(n_na)
  } else {
    if (identical(object$meta_type, "stochastic")) {
      method[j] <- "MLE"
    } else {
      # Deterministic metafrontier: meta_solver records which solver
      # produced the estimate; NULL on objects from older versions,
      # where the LP was the only solver.
      solver <- object$meta_solver
      method[j] <- if (is.null(solver)) {
        "LP"
      } else {
        switch(solver,
               "lp" = "LP",
               "qp" = "QP",
               "qp-barrier" = "QP (barrier)",
               solver)
      }
    }
    code[j] <- if (is.null(object$meta_convergence)) {
      NA_integer_
    } else {
      as.integer(object$meta_convergence)
    }
    converged[j] <- if (is.na(code[j])) NA else code[j] == 0L
    note[j] <- ""
  }

  data.frame(
    stage = stage,
    method = method,
    code = code,
    converged = converged,
    note = note,
    stringsAsFactors = FALSE
  )
}


#' Check Convergence of All Estimation Stages
#'
#' Reports the convergence status of every estimation stage of a
#' fitted metafrontier model: each group-level frontier and the
#' metafrontier itself. This makes it easy to verify that all
#' optimisers (MLE) and mathematical programmes (LP/QP) finished
#' successfully before interpreting technology gap ratios,
#' confidence intervals, or efficiency decompositions.
#'
#' @param object a fitted model object.
#' @param ... additional arguments passed to methods.
#'
#' @return A data frame of class \code{"metafrontier_convergence"}
#'   with one row per estimation stage and columns:
#'   \describe{
#'     \item{stage}{stage label, e.g. \code{"group: G1"} or
#'       \code{"metafrontier"}}
#'     \item{method}{how the stage was estimated: \code{"MLE"},
#'       \code{"LP"}, \code{"QP"}, \code{"QP (barrier)"},
#'       \code{"DEA"}, or \code{"external"}}
#'     \item{code}{the integer convergence code returned by the
#'       optimiser (0 indicates success); \code{NA} for DEA stages
#'       and externally fitted groups}
#'     \item{converged}{logical convergence indicator; for DEA
#'       stages \code{TRUE} unless any efficiency score is
#'       \code{NA} (infeasible programme); \code{NA} for externally
#'       fitted groups}
#'     \item{note}{additional detail, e.g. the number of infeasible
#'       DEA programmes}
#'   }
#'
#' @details
#' For groups supplied via the \code{models} argument of
#' \code{\link{metafrontier}} the convergence status of the external
#' fitter is not tracked, and the corresponding rows carry
#' \code{NA} with an explanatory note.
#'
#' @examples
#' sim <- simulate_metafrontier(n_groups = 2, n_per_group = 50, seed = 42)
#' fit <- metafrontier(log_y ~ log_x1 + log_x2, data = sim$data,
#'                     group = "group")
#' check_convergence(fit)
#'
#' @export
check_convergence <- function(object, ...) {
  UseMethod("check_convergence")
}


#' @rdname check_convergence
#' @export
check_convergence.metafrontier <- function(object, ...) {
  tab <- .convergence_table(object)
  class(tab) <- c("metafrontier_convergence", "data.frame")
  tab
}


#' @rdname check_convergence
#' @export
check_convergence.lc_metafrontier <- function(object, ...) {
  # The latent class fit embeds a full metafrontier object estimated
  # on the MAP class allocation; report its stages (labelled by the
  # latent classes LC1, LC2, ...).
  check_convergence(object$metafrontier)
}


#' @rdname check_convergence
#' @export
check_convergence.malmquist_meta <- function(object, ...) {
  stop("check_convergence() is not available for 'malmquist_meta' ",
       "objects: the period-by-group frontier fits are not retained ",
       "on the returned object.", call. = FALSE)
}


#' @export
check_convergence.default <- function(object, ...) {
  stop("check_convergence() is not implemented for objects of class '",
       paste(class(object), collapse = "', '"),
       "'. It supports 'metafrontier' and 'lc_metafrontier' objects.",
       call. = FALSE)
}


#' @export
print.metafrontier_convergence <- function(x, ...) {
  cat("\nConvergence of estimation stages\n")
  cat("--------------------------------\n")
  print.data.frame(x, row.names = FALSE)
  if (any(!x$converged, na.rm = TRUE)) {
    cat("\nOne or more stages did not converge; interpret TGRs, ",
        "confidence intervals,\nand decompositions with caution.\n",
        sep = "")
  }
  invisible(x)
}
