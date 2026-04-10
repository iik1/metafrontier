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
                                      ...) {
  type <- match.arg(type)

  switch(type,
    meta  = object$te_meta,
    group = object$te_group,
    tgr   = object$tgr
  )
}
