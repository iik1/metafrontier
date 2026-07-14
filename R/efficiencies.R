#' Extract Efficiency Scores from a Metafrontier Model
#'
#' Extracts technical efficiency scores from a fitted metafrontier
#' model. Returns either group-specific efficiency (\eqn{TE}),
#' metafrontier efficiency (\eqn{TE^* = TE \times TGR}), or the
#' technology gap ratio (\eqn{TGR}).
#'
#' @param object a fitted \code{"metafrontier"} object.
#' @param type character. The type of efficiency to return:
#'   \code{"group"} for efficiency relative to the group frontier,
#'   \code{"meta"} (default) for efficiency relative to the
#'   metafrontier, or \code{"tgr"} for the technology gap ratio.
#' @param estimator optional character. Override the efficiency
#'   estimator used at fit time: \code{"bc88"} for the Battese-Coelli
#'   (1988) conditional expectation \eqn{E[\exp(-u)|\varepsilon]} or
#'   \code{"jlms"} for \eqn{\exp(-E[u|\varepsilon])} (Jondrow et al.,
#'   1982). Both are stored on SFA fits, so no refitting is needed;
#'   \code{type = "meta"} is recomputed as \eqn{TE \times TGR}. The
#'   TGR itself does not depend on the estimator. Ignored (with a
#'   warning) for DEA fits and externally fitted group models that do
#'   not carry both estimators. Default \code{NULL} returns the
#'   scores selected at fit time.
#' @param ... additional arguments (currently unused).
#'
#' @return A numeric vector of efficiency scores of length
#'   \code{nobs(object)}.
#'
#' @details
#' The fundamental metafrontier decomposition is:
#' \deqn{TE^*_i = TE_i \times TGR_i}
#' where \eqn{TE_i} is efficiency relative to the group frontier
#' (returned by \code{type = "group"}) and \eqn{TGR_i} is the
#' technology gap ratio (returned by \code{type = "tgr"}).
#'
#' @examples
#' set.seed(42)
#' sim <- simulate_metafrontier(n_groups = 2, n_per_group = 100)
#' fit <- metafrontier(log_y ~ log_x1 + log_x2,
#'                     data = sim$data, group = "group")
#'
#' # Group-level efficiency
#' te <- efficiencies(fit, type = "group")
#'
#' # Metafrontier efficiency
#' te_star <- efficiencies(fit, type = "meta")
#'
#' # Verify decomposition: TE* = TE x TGR
#' tgr <- efficiencies(fit, type = "tgr")
#' all.equal(te_star, te * tgr)
#'
#' @seealso \code{\link{technology_gap_ratio}}, \code{\link{metafrontier}}
#'
#' @export
efficiencies <- function(object, ...) {
  UseMethod("efficiencies")
}

#' @rdname efficiencies
#' @export
efficiencies.metafrontier <- function(object,
                                      type = c("meta", "group", "tgr"),
                                      estimator = NULL,
                                      ...) {
  type <- match.arg(type)

  if (is.null(estimator)) {
    return(switch(type,
      meta  = object$te_meta,
      group = object$te_group,
      tgr   = object$tgr
    ))
  }

  estimator <- match.arg(estimator, c("bc88", "jlms"))
  field <- paste0("efficiency_", estimator)

  # Rebuild group-level TE from the stored per-estimator vectors
  te_group <- object$te_group
  available <- TRUE
  for (g in object$groups) {
    gm <- object$group_models[[g]]
    if (is.null(gm[[field]])) {
      available <- FALSE
      break
    }
    te_group[object$group_vec == g] <- gm[[field]]
  }

  if (!available) {
    warning("Estimator-specific efficiencies are not stored for all ",
            "groups (DEA or externally fitted models); returning the ",
            "scores selected at fit time.", call. = FALSE)
    te_group <- object$te_group
  }

  switch(type,
    group = te_group,
    tgr   = object$tgr,
    meta  = te_group * object$tgr
  )
}
