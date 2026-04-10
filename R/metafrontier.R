#' Estimate a Metafrontier Production Function
#'
#' Estimates group-specific frontiers and a metafrontier that envelops
#' all group technologies. Supports both SFA-based (parametric) and
#' DEA-based (nonparametric) approaches, with deterministic (Battese,
#' Rao, and O'Donnell, 2004) or stochastic (Huang, Huang, and Liu, 2014)
#' metafrontier estimation.
#'
#' @param formula an object of class \code{\link[Formula]{Formula}}.
#'   The left-hand side specifies the (log) output variable. The first
#'   right-hand side part specifies inputs for the frontier. An optional
#'   second part (separated by \code{|}) specifies inefficiency
#'   determinants. Example: \code{log_y ~ log_x1 + log_x2 | z1 + z2}.
#'   Ignored if \code{models} is provided.
#' @param data a data frame containing all variables in the formula
#'   and the grouping variable. Ignored if \code{models} is provided.
#' @param group a character string naming the column in \code{data}
#'   that identifies technology groups, or a vector of group indicators
#'   of length \code{nrow(data)}. Ignored if \code{models} is provided.
#' @param method character. The frontier estimation method for
#'   group-specific models: \code{"sfa"} (default) for stochastic
#'   frontier analysis or \code{"dea"} for data envelopment analysis.
#' @param meta_type character. The method for estimating the
#'   metafrontier: \code{"deterministic"} (default) uses the linear
#'   programming approach of Battese, Rao, and O'Donnell (2004);
#'   \code{"stochastic"} uses the second-stage SFA approach of Huang,
#'   Huang, and Liu (2014).
#' @param dist character. Distribution of the one-sided inefficiency
#'   term in SFA models. One of \code{"hnormal"} (half-normal, default),
#'   \code{"tnormal"} (truncated normal), or \code{"exponential"}.
#'   Ignored when \code{method = "dea"}.
#' @param orientation character. For DEA: \code{"output"} (default) or
#'   \code{"input"} orientation. Ignored when \code{method = "sfa"}.
#' @param rts character. Returns to scale for DEA: \code{"crs"}
#'   (constant, default), \code{"vrs"} (variable), \code{"drs"}
#'   (decreasing), or \code{"irs"} (increasing). Ignored when
#'   \code{method = "sfa"}.
#' @param models an optional named list of pre-fitted group-specific
#'   frontier models (objects from \pkg{sfaR}, \pkg{frontier}, or
#'   \pkg{Benchmarking}). If provided, \code{formula}, \code{data},
#'   and \code{group} are ignored.
#' @param control a list of control parameters for the optimiser.
#'   See Details.
#' @param ... additional arguments passed to the group-level
#'   estimation functions.
#'
#' @return An object of class \code{"metafrontier"} (with subclass
#'   \code{"metafrontier_sfa"} or \code{"metafrontier_dea"}),
#'   containing:
#'   \describe{
#'     \item{call}{the matched function call}
#'     \item{group_models}{list of fitted group-specific models}
#'     \item{meta_coef}{estimated metafrontier parameters}
#'     \item{group_coef}{list of group-specific coefficient vectors}
#'     \item{tgr}{technology gap ratios for each observation}
#'     \item{te_group}{group-specific technical efficiency}
#'     \item{te_meta}{metafrontier technical efficiency (TE* = TE x TGR)}
#'     \item{logLik_groups}{log-likelihoods of group models}
#'     \item{nobs}{number of observations per group and total}
#'     \item{groups}{group labels}
#'     \item{method}{estimation method used}
#'     \item{meta_type}{metafrontier type used}
#'     \item{convergence}{convergence status}
#'   }
#'
#' @details
#' The metafrontier framework decomposes efficiency relative to a
#' global technology into two components:
#'
#' \deqn{TE^*_i = TE_i \times TGR_i}
#'
#' where \eqn{TE_i} is efficiency relative to the group frontier and
#' \eqn{TGR_i} is the technology gap ratio measuring how close the
#' group frontier is to the metafrontier.
#'
#' The deterministic metafrontier (Battese, Rao, and O'Donnell, 2004)
#' is estimated by minimising the sum of squared deviations from group
#' frontiers subject to the constraint that the metafrontier envelops
#' all group frontiers. The stochastic metafrontier (Huang, Huang, and
#' Liu, 2014) replaces this with a second-stage SFA, providing a
#' distributional framework for inference on the TGR.
#'
#' @references
#' Battese, G.E., Rao, D.S.P. and O'Donnell, C.J. (2004). A
#' metafrontier production function for estimation of technical
#' efficiencies and technology gaps for firms operating under different
#' technologies. \emph{Journal of Productivity Analysis}, 21(1),
#' 91--103. \doi{10.1023/B:PROD.0000016869.68578.b7}
#'
#' Huang, C.J., Huang, T.-H. and Liu, N.-H. (2014). A new approach
#' to estimating the metafrontier production function based on a
#' stochastic frontier framework. \emph{Journal of Productivity
#' Analysis}, 42(3), 241--254. \doi{10.1007/s11123-014-0397-1}
#'
#' O'Donnell, C.J., Rao, D.S.P. and Battese, G.E. (2008).
#' Metafrontier frameworks for the study of firm-level efficiencies
#' and technology ratios. \emph{Empirical Economics}, 34(2), 231--255.
#' \doi{10.1007/s00181-006-0095-0}
#'
#' @examples
#' # Simulate metafrontier data
#' set.seed(42)
#' sim <- simulate_metafrontier(n_groups = 2, n_per_group = 100)
#'
#' # Estimate deterministic SFA metafrontier (BRO 2004)
#' fit <- metafrontier(log_y ~ log_x1 + log_x2,
#'                     data = sim$data,
#'                     group = "group",
#'                     method = "sfa",
#'                     meta_type = "deterministic")
#' summary(fit)
#'
#' # Technology gap ratios
#' tgr <- technology_gap_ratio(fit)
#' summary(tgr)
#'
#' @export
metafrontier <- function(formula = NULL,
                         data = NULL,
                         group = NULL,
                         method = c("sfa", "dea"),
                         meta_type = c("deterministic", "stochastic"),
                         dist = c("hnormal", "tnormal", "exponential"),
                         orientation = c("output", "input"),
                         rts = c("crs", "vrs", "drs", "irs"),
                         models = NULL,
                         control = list(),
                         ...) {

  cl <- match.call()
  method <- match.arg(method)
  meta_type <- match.arg(meta_type)
  dist <- match.arg(dist)
  orientation <- match.arg(orientation)
  rts <- match.arg(rts)

  # ------ Input validation ------
  if (is.null(models)) {
    if (is.null(formula) || is.null(data) || is.null(group)) {
      stop("Either provide 'models' or all of 'formula', 'data', and 'group'.",
           call. = FALSE)
    }
    result <- .estimate_from_data(formula, data, group, method,
                                  meta_type, dist, orientation,
                                  rts, control, ...)
  } else {
    result <- .estimate_from_models(models, meta_type, control, ...)
  }

  result$call <- cl
  result$method <- method
  result$meta_type <- meta_type

  class(result) <- c(
    paste0("metafrontier_", method),
    "metafrontier"
  )

  result
}


# ---------- Internal: estimate from raw data ----------
.estimate_from_data <- function(formula, data, group, method,
                                meta_type, dist, orientation,
                                rts, control, ...) {

  # Parse the group variable
  if (is.character(group) && length(group) == 1L) {
    if (!group %in% names(data)) {
      stop("Column '", group, "' not found in data.", call. = FALSE)
    }
    group_vec <- data[[group]]
  } else {
    if (length(group) != nrow(data)) {
      stop("'group' must be a column name or a vector of length nrow(data).",
           call. = FALSE)
    }
    group_vec <- group
  }

  group_vec <- as.factor(group_vec)
  group_levels <- levels(group_vec)
  n_groups <- length(group_levels)

  if (n_groups < 2L) {
    stop("At least 2 groups are required for metafrontier analysis.",
         call. = FALSE)
  }

  # Parse formula
  f <- Formula::Formula(formula)

  # Estimate group-specific models
  group_models <- vector("list", n_groups)
  names(group_models) <- group_levels

  for (g in group_levels) {
    idx <- which(group_vec == g)
    n_g <- length(idx)
    if (n_g < 3L) {
      warning("Group '", g, "' has fewer than 3 observations. ",
              "Results may be unreliable.", call. = FALSE)
    }
    data_g <- data[idx, , drop = FALSE]

    if (method == "sfa") {
      group_models[[g]] <- .fit_sfa_group(f, data_g, dist, control, ...)
    } else {
      group_models[[g]] <- .fit_dea_group(f, data_g, orientation, rts, ...)
    }
  }

  # Estimate the metafrontier
  if (method == "sfa") {
    meta_result <- .estimate_sfa_metafrontier(
      f, data, group_vec, group_levels, group_models,
      meta_type, dist, control
    )
  } else {
    meta_result <- .estimate_dea_metafrontier(
      f, data, group_vec, group_levels, group_models,
      orientation, rts
    )
  }

  meta_result$group_models <- group_models
  meta_result$groups <- group_levels
  meta_result$group_vec <- group_vec
  meta_result$formula <- f
  meta_result$data <- data
  meta_result$nobs <- c(
    total = nrow(data),
    setNames(table(group_vec), group_levels)
  )

  meta_result
}


# ---------- Internal: estimate from pre-fitted models ----------
.estimate_from_models <- function(models, meta_type, control, ...) {
  if (!is.list(models) || is.null(names(models))) {
    stop("'models' must be a named list of fitted frontier objects.",
         call. = FALSE)
  }

  # TODO: Extract coefficients and data from sfaR/frontier/Benchmarking
  # objects and dispatch to the appropriate metafrontier estimator
  stop("Estimation from pre-fitted models is not yet implemented. ",
       "Please use the formula interface.", call. = FALSE)
}
