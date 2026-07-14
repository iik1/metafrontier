#' Extract Technology Gap Ratios
#'
#' Extracts the technology gap ratios (TGR) from a fitted metafrontier
#' model. The TGR measures how close a group's production frontier is
#' to the global metafrontier at each observation's input mix.
#'
#' @param object a fitted \code{"metafrontier"} object.
#' @param by_group logical. If \code{TRUE} (default), returns a named
#'   list of TGR vectors, one per group. If \code{FALSE}, returns a
#'   single numeric vector.
#' @param ... additional arguments (currently unused).
#'
#' @return If \code{by_group = TRUE}, a named list of numeric vectors.
#'   If \code{by_group = FALSE}, a numeric vector of length
#'   \code{nobs(object)}.
#'
#' @details
#' The technology gap ratio is defined as:
#' \deqn{TGR_i = \frac{f(x_i; \hat\beta_j)}{f(x_i; \hat\beta^*)}}
#' for SFA-based metafrontiers, and
#' \deqn{TGR_i = \frac{TE^*_i}{TE^{group}_i}}
#' for DEA-based metafrontiers.
#'
#' Under the deterministic metafrontier and DEA, TGR lies in (0, 1]
#' by construction. Under the stochastic metafrontier of Huang, Huang,
#' and Liu (2014), TGR can exceed 1 for some observations because the
#' metafrontier need not envelop the group frontiers at every point.
#'
#' A TGR of 1 means the group frontier coincides with the
#' metafrontier at that input mix. Values less than 1 indicate a
#' technology gap. Since each group technology is a restricted subset
#' of the common metatechnology (Battese, Rao and O'Donnell, 2004),
#' the gap reflects the restrictions a group faces (regulation,
#' environment, endowments) rather than a fundamentally different
#' technology.
#'
#' @examples
#' set.seed(42)
#' sim <- simulate_metafrontier(n_groups = 2, n_per_group = 100)
#' fit <- metafrontier(log_y ~ log_x1 + log_x2,
#'                     data = sim$data, group = "group")
#' tgr <- technology_gap_ratio(fit)
#' lapply(tgr, summary)
#'
#' @seealso \code{\link{metafrontier}}, \code{\link{efficiencies.metafrontier}}
#'
#' @export
technology_gap_ratio <- function(object, by_group = TRUE, ...) {
  if (!inherits(object, "metafrontier")) {
    stop("'object' must be of class 'metafrontier'.", call. = FALSE)
  }

  tgr <- object$tgr

  if (by_group) {
    group_vec <- object$group_vec
    group_levels <- object$groups
    result <- lapply(group_levels, function(g) {
      tgr[group_vec == g]
    })
    names(result) <- group_levels
    return(result)
  }

  tgr
}


#' Summary of Technology Gap Ratios
#'
#' Prints a summary table of TGR statistics by group. The underlying
#' observation-level TGR values are the same as those returned by
#' \code{efficiencies(object, type = "tgr")}.
#'
#' @param object a fitted \code{"metafrontier"} object.
#' @param ... additional arguments (currently unused).
#'
#' @return A data frame with columns: Group, N, Mean, SD, Min, Q1,
#'   Median, Q3, Max.
#'
#' @seealso \code{\link{efficiencies.metafrontier}},
#'   \code{\link{technology_gap_ratio}}, \code{\link{boot_tgr}}
#'
#' @examples
#' sim <- simulate_metafrontier(n_groups = 2, n_per_group = 50,
#'                              seed = 42)
#' fit <- metafrontier(log_y ~ log_x1 + log_x2, data = sim$data,
#'                     group = "group", method = "sfa",
#'                     meta_type = "deterministic")
#' tgr_summary(fit)
#'
#' @export
tgr_summary <- function(object, ...) {
  if (!inherits(object, "metafrontier")) {
    stop("'object' must be of class 'metafrontier'.", call. = FALSE)
  }

  tgr_list <- technology_gap_ratio(object, by_group = TRUE)

  do.call(rbind, lapply(names(tgr_list), function(g) {
    x <- tgr_list[[g]]
    data.frame(
      Group = g,
      N = length(x),
      Mean = mean(x),
      SD = sd(x),
      Min = min(x),
      Q1 = quantile(x, 0.25),
      Median = median(x),
      Q3 = quantile(x, 0.75),
      Max = max(x),
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  }))
}
