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
#'   (decreasing), \code{"irs"} (increasing), or \code{"fdh"} (free
#'   disposable hull, i.e. no convexity). Ignored when
#'   \code{method = "sfa"}.
#' @param models an optional named list of pre-fitted group-specific
#'   frontier models (objects from \pkg{sfaR} or \pkg{frontier}, or
#'   hand-built lists). Fitted model objects are converted
#'   automatically via \code{\link{as_metafrontier_model}}, so no
#'   manual conversion is required (pre-converting is harmless, the
#'   conversion is idempotent). Farrell objects from
#'   \pkg{Benchmarking} store neither coefficients nor data and
#'   cannot be used here; use the formula interface with
#'   \code{method = "dea"} instead. If \code{models} is provided,
#'   \code{formula}, \code{data}, and \code{group} are ignored.
#' @param panel an optional list with components \code{id} and
#'   \code{time} naming the panel identifier and time columns in
#'   \code{data}. When non-NULL, panel SFA models (BC92/BC95) are
#'   used at the group level.
#' @param panel_dist character. Panel SFA model: \code{"bc92"}
#'   (Battese and Coelli 1992, time-varying inefficiency, default)
#'   or \code{"bc95"} (Battese and Coelli 1995, observation-specific
#'   mean). Only used when \code{panel} is non-NULL.
#' @param type character. For DEA: \code{"radial"} (default) for
#'   standard radial DEA, \code{"directional"} for directional
#'   distance functions, or \code{"hyperbolic"} for hyperbolic
#'   (graph) efficiency, which contracts inputs and expands outputs
#'   simultaneously.
#' @param direction direction vector for DDF. Either a character
#'   preset (\code{"proportional"} (default), \code{"output"}, or
#'   \code{"input"}), a numeric vector of length m + s giving a
#'   common direction (first m elements for inputs, last s for
#'   outputs), or a numeric n x (m + s) matrix of firm-specific
#'   directions. With numeric directions the ratio-based TGR is not
#'   defined; the additive gap (\code{ddf_gap}) is reported instead.
#'   Only used when \code{type = "directional"}.
#' @param control a named list of control parameters passed to
#'   \code{\link[stats]{optim}}. Common options include
#'   \code{maxit} (maximum iterations, default 5000),
#'   \code{reltol} (relative convergence tolerance, default 1e-10),
#'   and \code{fnscale} (set to -1 internally for maximisation).
#' @param estimator character. Technical efficiency estimator for
#'   SFA models: \code{"bc88"} (default) for the conditional
#'   expectation \eqn{E[\exp(-u)|\varepsilon]} of Battese and Coelli
#'   (1988), which is the consistent estimator of technical
#'   efficiency, or \code{"jlms"} for
#'   \eqn{\exp(-E[u|\varepsilon])} following Jondrow et al. (1982).
#'   Both are stored on the fitted object; see
#'   \code{\link{efficiencies}}. Ignored when \code{method = "dea"}.
#' @param objective character. Identification criterion for the
#'   deterministic metafrontier: \code{"lp"} (default) minimises the
#'   sum of absolute deviations (a linear programme), \code{"qp"}
#'   minimises the sum of squared deviations (a quadratic programme,
#'   solved exactly via \pkg{quadprog} when available). Both criteria
#'   are proposed in Battese, Rao, and O'Donnell (2004). Only used
#'   when \code{method = "sfa"} and \code{meta_type =
#'   "deterministic"}.
#' @param engine character. Estimation backend for the group
#'   frontiers: \code{"internal"} (default) uses the package's own
#'   estimators; \code{"sfaR"} or \code{"frontier"} delegate the SFA
#'   group frontiers to \code{\link[sfaR]{sfacross}} or
#'   \code{\link[frontier]{sfa}} (cross-sectional, single-part
#'   formulas only); \code{"Benchmarking"} delegates the DEA group
#'   frontiers and the pooled metafrontier to
#'   \code{\link[Benchmarking]{dea}} (radial only), using its
#'   \code{XREF}/\code{YREF} external-reference facility for the
#'   metafrontier stage. The metafrontier stage for SFA methods is
#'   always estimated internally (the Murphy-Topel correction
#'   requires the internal likelihood).
#' @param slack logical. For radial DEA, compute second-stage input
#'   and output slacks (with the radial score held fixed) against
#'   both the group and the pooled reference sets. Default
#'   \code{FALSE}.
#' @param frontier character. \code{"production"} (default), with
#'   composed error \eqn{\varepsilon = v - u}, or \code{"cost"}, with
#'   \eqn{\varepsilon = v + u} and a cost metafrontier that lies on or
#'   below every group cost frontier. Efficiencies and technology gap
#'   ratios lie in (0, 1] in both cases. SFA methods only.
#' @param selection an optional two-sided formula \code{d ~ z1 + z2}
#'   naming a 0/1 selection indicator and its regressors. Each group
#'   frontier is then estimated as the sample-selection stochastic
#'   frontier of Greene (2010): a probit selection equation, and a
#'   frontier observed when \code{d = 1} whose noise is correlated
#'   (\code{rho}) with the selection error. Rows with \code{d = 0} enter
#'   only the probit and may have a missing outcome; efficiencies and
#'   technology gap ratios are reported for the selected rows.
#'   Half-normal or exponential inefficiency; cross-sectional SFA with
#'   the internal engine only.
#' @param spatial an optional list \code{list(W = , row_standardise =
#'   FALSE)}. \code{W} is an n x n spatial weights matrix
#'   for the rows of \code{data} (a matrix, a \pkg{Matrix} matrix or an
#'   \pkg{spdep} \code{listw} object) with a zero diagonal. Each group
#'   frontier is then the spatial-autoregressive stochastic frontier
#'   \eqn{y = \rho_g W_g y + x'\beta_g + v - u} (Glass, Kenjegalieva and
#'   Sickles, 2016), where \eqn{W_g} is the block of \code{W} for the
#'   group's own units (re-row-standardised when \code{row_standardise =
#'   TRUE}), so spillovers between groups are excluded. Each group model
#'   reports \code{rho} and the direct, indirect and total input effects
#'   (\code{effects}). The metafrontier envelops the long-run group
#'   frontiers \eqn{m_g x'\beta_g}, with \eqn{m_g} the mean row sum of
#'   \eqn{(I - \rho_g W_g)^{-1}} (\eqn{1 / (1 - \rho_g)} for
#'   row-standardised weights): with row-standardised weights the
#'   intercept and \eqn{\rho_g} are nearly collinear, so the direct
#'   technology \eqn{x'\beta_g} is poorly determined while the long-run
#'   frontier is not. Cross-sectional SFA without
#'   inefficiency determinants, internal engine only.
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
#'     \item{meta_convergence}{integer convergence code for the
#'       metafrontier stage (0 = success; \code{\link[stats]{optim}}
#'       codes for the stochastic metafrontier and the QP barrier
#'       fallback; 0 for a successful LP or DEA solution). Each SFA
#'       group model in \code{group_models} additionally carries its
#'       own \code{convergence} code. Use
#'       \code{\link{check_convergence}} to inspect all stages.}
#'     \item{estimator, objective, engine, meta_solver}{the
#'       estimation choices used for the fit}
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
#' The deterministic metafrontier is identified by one of the two
#' criteria proposed by Battese, Rao, and O'Donnell (2004), subject in
#' both cases to the constraint that the metafrontier envelops all
#' group frontiers: minimising the sum of absolute deviations, which
#' reduces to a linear programme because the envelope constraints
#' force every deviation to be non-negative (O'Donnell, Rao, and
#' Battese, 2008, Eqs. 23-25), or minimising the sum of squared
#' deviations, a convex quadratic programme. The LP
#' (\code{objective = "lp"}, the default) is solved via
#' \pkg{lpSolveAPI}; the QP (\code{objective = "qp"}) is solved
#' exactly via \pkg{quadprog} when available, with an adaptive-barrier
#' fallback via \code{constrOptim()}.
#' The stochastic metafrontier (Huang, Huang, and Liu, 2014) replaces
#' this with a second-stage SFA, providing a distributional framework
#' for inference on the TGR.
#'
#' \strong{Convergence and failure handling:} estimation stops with an
#' error only when no usable estimate exists (for example, when both
#' the BFGS and Nelder-Mead optimisers fail for a group frontier).
#' When an optimiser stops at a non-zero convergence code, the fitted
#' object is returned with a warning and the code is recorded; use
#' \code{\link{check_convergence}} or \code{summary()} to verify all
#' estimation stages before interpreting technology gap ratios,
#' confidence intervals, or productivity decompositions. Infeasible
#' DEA programmes yield \code{NA} efficiency scores, accompanied by a
#' warning and counted by \code{\link{check_convergence}}.
#'
#' \strong{Note on standard errors (stochastic metafrontier):}
#' The stochastic metafrontier is a two-stage estimator. Stage 2 treats
#' the fitted group frontier values as data, so the reported standard
#' errors, confidence intervals, and variance-covariance matrix do not
#' account for estimation uncertainty from Stage 1 (the
#' generated-regressor problem; see Murphy and Topel, 1985). Use
#' \code{vcov(fit, correction = "murphy-topel")} or bootstrap-based
#' confidence intervals via \code{\link{boot_tgr}} for corrected
#' inference.
#'
#' \strong{Note on frontier orientation (SFA path):}
#' The SFA path estimates production frontiers
#' (\eqn{\varepsilon = v - u}) by default and cost frontiers
#' (\eqn{\varepsilon = v + u}) with \code{frontier = "cost"}. For cost
#' frontiers the deterministic metafrontier envelops the group cost
#' frontiers from below, and TGR is the ratio of metatechnology to
#' group minimum cost. The DEA path supports both
#' \code{orientation = "output"} and \code{orientation = "input"}.
#'
#' @references
#' Battese, G.E., Rao, D.S.P. and O'Donnell, C.J. (2004). A
#' metafrontier production function for estimation of technical
#' efficiencies and technology gaps for firms operating under different
#' technologies. \emph{Journal of Productivity Analysis}, 21(1),
#' 91--103. \doi{10.1023/B:PROD.0000012454.06094.29}
#'
#' Huang, C.J., Huang, T.-H. and Liu, N.-H. (2014). A new approach
#' to estimating the metafrontier production function based on a
#' stochastic frontier framework. \emph{Journal of Productivity
#' Analysis}, 42(3), 241--254. \doi{10.1007/s11123-014-0402-2}
#'
#' O'Donnell, C.J., Rao, D.S.P. and Battese, G.E. (2008).
#' Metafrontier frameworks for the study of firm-level efficiencies
#' and technology ratios. \emph{Empirical Economics}, 34(2), 231--255.
#' \doi{10.1007/s00181-007-0119-4}
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
                         rts = c("crs", "vrs", "drs", "irs", "fdh"),
                         models = NULL,
                         panel = NULL,
                         panel_dist = c("bc92", "bc95"),
                         type = c("radial", "directional", "hyperbolic"),
                         direction = c("proportional", "output", "input"),
                         control = list(),
                         estimator = c("bc88", "jlms"),
                         objective = c("lp", "qp"),
                         engine = c("internal", "sfaR", "frontier", "Benchmarking"),
                         slack = FALSE,
                         frontier = c("production", "cost"),
                         selection = NULL,
                         spatial = NULL,
                         ...) {

  cl <- match.call()
  frontier_given <- !missing(frontier)
  frontier <- match.arg(frontier)
  method <- match.arg(method)
  meta_type <- match.arg(meta_type)
  dist <- match.arg(dist)
  orientation <- match.arg(orientation)
  rts <- match.arg(rts)
  panel_dist <- match.arg(panel_dist)
  type <- match.arg(type)
  if (is.character(direction)) {
    direction <- match.arg(direction)
  } else if (!is.numeric(direction)) {
    stop("'direction' must be a character preset or a numeric ",
         "vector/matrix.", call. = FALSE)
  }
  estimator <- match.arg(estimator)
  objective <- match.arg(objective)
  engine <- match.arg(engine)
  if (!is.logical(slack) || length(slack) != 1L || is.na(slack)) {
    stop("'slack' must be TRUE or FALSE.", call. = FALSE)
  }

  # ------ Engine / option compatibility ------
  if (engine == "Benchmarking") {
    if (method != "dea") {
      stop("engine = \"Benchmarking\" is only available for ",
           "method = \"dea\".", call. = FALSE)
    }
    if (type != "radial") {
      stop("engine = \"Benchmarking\" supports type = \"radial\" only; ",
           "use engine = \"internal\" for directional or hyperbolic ",
           "efficiency.", call. = FALSE)
    }
  }
  if (engine %in% c("sfaR", "frontier")) {
    if (method != "sfa") {
      stop("engine = \"", engine, "\" is only available for ",
           "method = \"sfa\".", call. = FALSE)
    }
    if (!is.null(panel)) {
      stop("External SFA engines support cross-sectional fits only; ",
           "use engine = \"internal\" for panel models.", call. = FALSE)
    }
    if (engine == "frontier" && dist != "hnormal") {
      stop("engine = \"frontier\" supports dist = \"hnormal\" only.",
           call. = FALSE)
    }
  }
  if (frontier == "cost" && method != "sfa") {
    stop("frontier = \"cost\" is available for SFA metafrontiers only; ",
         "cost DEA requires input prices.", call. = FALSE)
  }
  if (!is.null(selection) &&
      (method != "sfa" || !is.null(panel) || engine != "internal" ||
       !is.null(models))) {
    stop("'selection' is available for cross-sectional SFA with the ",
         "internal engine and the formula interface only.", call. = FALSE)
  }
  if (!is.null(selection) && length(Formula::Formula(formula))[2] > 1L) {
    stop("'selection' does not support inefficiency determinants.",
         call. = FALSE)
  }
  if (!is.null(spatial)) {
    if (method != "sfa" || !is.null(panel) || !is.null(selection) ||
        engine != "internal" || !is.null(models)) {
      stop("'spatial' is available for cross-sectional SFA with the ",
           "internal engine and the formula interface only, without ",
           "'selection'.", call. = FALSE)
    }
    if (length(Formula::Formula(formula))[2] > 1L) {
      stop("'spatial' does not support inefficiency determinants.",
           call. = FALSE)
    }
    spatial$W <- .spatial_weights(spatial$W, nrow(data))
  }
  if (objective == "qp" && (method != "sfa" ||
                            meta_type != "deterministic")) {
    warning("'objective' applies only to the deterministic SFA ",
            "metafrontier; ignored.", call. = FALSE)
  }
  if (slack && (method != "dea" || type != "radial")) {
    warning("'slack' is only used for radial DEA; ignored.",
            call. = FALSE)
    slack <- FALSE
  }

  # ------ Input validation ------
  if (is.null(models)) {
    if (is.null(formula) || is.null(data) || is.null(group)) {
      stop("Either provide 'models' or all of 'formula', 'data', and 'group'.",
           call. = FALSE)
    }
    result <- .estimate_from_data(formula, data, group, method,
                                  meta_type, dist, orientation,
                                  rts, control,
                                  panel = panel,
                                  panel_dist = panel_dist,
                                  type = type,
                                  direction = direction,
                                  estimator = estimator,
                                  objective = objective,
                                  engine = engine,
                                  slack = slack,
                                  frontier = frontier,
                                  selection = selection,
                                  spatial = spatial, ...)
  } else {
    result <- .estimate_from_models(models, meta_type, control,
                                    objective = objective, ...)
    # The frontier type of imported models is read from the models
    if (frontier_given && frontier != result$frontier) {
      stop("frontier = \"", frontier, "\" does not match the imported ",
           "models, which are ", result$frontier, " frontiers.",
           call. = FALSE)
    }
    frontier <- result$frontier
  }

  result$call <- cl
  result$method <- method
  result$meta_type <- meta_type
  result$orientation <- orientation
  result$rts <- rts
  result$type <- type
  result$direction <- direction
  result$slack <- slack
  result$engine <- engine
  result$frontier <- frontier
  result$selection <- selection
  if (!is.null(spatial)) result$spatial <- spatial["row_standardise"]
  if (method == "sfa") {
    result$estimator <- estimator
    if (meta_type == "deterministic") result$objective <- objective
  }
  if (is.character(group) && length(group) == 1L) {
    result$group_col <- group
  } else {
    result$group_col <- "group"
  }

  class(result) <- c(
    paste0("metafrontier_", method),
    "metafrontier"
  )

  result
}


# ---------- Internal: estimate from raw data ----------
.estimate_from_data <- function(formula, data, group, method,
                                meta_type, dist, orientation,
                                rts, control,
                                panel = NULL,
                                panel_dist = "bc92",
                                type = "radial",
                                direction = "proportional",
                                estimator = "bc88",
                                objective = "lp",
                                engine = "internal",
                                slack = FALSE,
                                frontier = "production",
                                selection = NULL,
                                spatial = NULL, ...) {

  sign <- if (frontier == "cost") -1 else 1

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

  # Parse formula
  f <- Formula::Formula(formula)

  # Drop rows with a missing group or model variable once, here, so the
  # group fits, the pooled stage and every per-row result line up
  keep <- complete.cases(model.frame(f, data = data, na.action = na.pass)) &
    !is.na(group_vec)
  if (!is.null(selection)) {
    # Unselected rows enter only the probit, so their frontier variables
    # (typically the outcome) may be missing
    smf <- model.frame(selection, data = data, na.action = na.pass)
    d <- model.response(smf)
    sel_ok <- complete.cases(smf)
    if (!all(d[sel_ok] %in% c(0, 1))) {
      stop("The selection indicator must be 0/1 or logical.", call. = FALSE)
    }
    keep <- sel_ok & !is.na(group_vec) & (keep | (sel_ok & d == 0))
  }
  if (!all(keep)) {
    data <- data[keep, , drop = FALSE]
    group_vec <- group_vec[keep]
    if (is.matrix(direction) && nrow(direction) == length(keep)) {
      direction <- direction[keep, , drop = FALSE]
    }
    if (!is.null(spatial)) spatial$W <- spatial$W[keep, keep, drop = FALSE]
  }

  group_vec <- droplevels(as.factor(group_vec))
  group_levels <- levels(group_vec)
  n_groups <- length(group_levels)

  if (n_groups < 2L) {
    stop("At least 2 groups are required for metafrontier analysis.",
         call. = FALSE)
  }

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
      if (engine %in% c("sfaR", "frontier")) {
        group_models[[g]] <- .fit_sfa_group_external(
          f, data_g, dist, engine, estimator, sign = sign
        )
      } else if (!is.null(panel)) {
        group_models[[g]] <- .fit_sfa_panel_group(
          f, data_g, dist, panel_dist, panel, control,
          estimator = estimator, sign = sign, ...
        )
      } else if (!is.null(selection)) {
        group_models[[g]] <- .fit_sfa_selection_group(
          f, selection, data_g, dist, control, estimator = estimator,
          sign = sign
        )
      } else if (!is.null(spatial)) {
        # Within-group weights only: no spillovers between technologies
        W_g <- spatial$W[idx, idx, drop = FALSE]
        if (isTRUE(spatial$row_standardise)) {
          rs <- rowSums(W_g)
          W_g[rs > 0, ] <- W_g[rs > 0, , drop = FALSE] / rs[rs > 0]
        }
        group_models[[g]] <- .fit_sfa_spatial_group(
          f, data_g, W_g, dist, control, estimator = estimator, sign = sign
        )
      } else {
        group_models[[g]] <- .fit_sfa_group(f, data_g, dist, control,
                                            estimator = estimator,
                                            sign = sign, ...)
      }
    } else {
      if (engine == "Benchmarking") {
        group_models[[g]] <- .fit_dea_group_benchmarking(
          f, data_g, orientation, rts
        )
      } else if (type == "directional") {
        # Firm-specific direction matrices are subset to the group rows
        dir_g <- if (is.matrix(direction)) {
          direction[idx, , drop = FALSE]
        } else {
          direction
        }
        group_models[[g]] <- .fit_ddf_group(f, data_g, rts, dir_g)
      } else if (type == "hyperbolic") {
        group_models[[g]] <- .fit_hyperbolic_group(f, data_g, rts)
      } else {
        group_models[[g]] <- .fit_dea_group(f, data_g, orientation, rts,
                                            slack = slack, ...)
      }
    }
  }

  # Under sample selection the metafrontier and every per-row result
  # refer to the selected rows
  if (!is.null(selection)) {
    selected <- model.response(model.frame(selection, data = data)) == 1
    data <- data[selected, , drop = FALSE]
    group_vec <- group_vec[selected]
  }

  # Estimate the metafrontier
  if (method == "sfa") {
    meta_result <- .estimate_sfa_metafrontier(
      f, data, group_vec, group_levels, group_models,
      meta_type, dist, control, objective = objective, sign = sign
    )
  } else {
    if (engine == "Benchmarking") {
      meta_result <- .estimate_dea_metafrontier_benchmarking(
        f, data, group_vec, group_levels, group_models,
        orientation, rts
      )
    } else if (type == "directional") {
      meta_result <- .estimate_ddf_metafrontier(
        f, data, group_vec, group_levels, group_models,
        rts, direction
      )
    } else if (type == "hyperbolic") {
      meta_result <- .estimate_hyperbolic_metafrontier(
        f, data, group_vec, group_levels, group_models, rts
      )
    } else {
      meta_result <- .estimate_dea_metafrontier(
        f, data, group_vec, group_levels, group_models,
        orientation, rts, slack = slack
      )
    }
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
.estimate_from_models <- function(models, meta_type, control,
                                  objective = "lp", ...) {
  if (!is.list(models) || is.null(names(models))) {
    stop("'models' must be a named list of fitted frontier objects.",
         call. = FALSE)
  }

  n_groups <- length(models)
  if (n_groups < 2L) {
    stop("At least 2 groups are required for metafrontier analysis.",
         call. = FALSE)
  }

  group_levels <- names(models)

  # Extract components from each model
  extracted <- lapply(models, .extract_model_components)

  # Validate that all models provide X, y, and beta

  for (g in names(extracted)) {
    ex <- extracted[[g]]
    if (is.null(ex$beta) || is.null(ex$X) || is.null(ex$y)) {
      stop("Model '", g, "' does not provide coefficients, X, or y. ",
           "Nonparametric (DEA) models from Benchmarking cannot be used ",
           "with the 'models' argument. Use the formula interface with ",
           "method = \"dea\" instead.",
           call. = FALSE)
    }
  }

  # A metafrontier envelops technologies of one kind: production or cost
  signs <- vapply(extracted, function(ex) {
    if (is.null(ex$sign)) 1 else ex$sign
  }, numeric(1))
  if (length(unique(signs)) > 1L) {
    stop("'models' mixes production and cost frontiers; a metafrontier ",
         "requires all group models to be of the same kind.", call. = FALSE)
  }
  sign <- signs[[1]]

  # Build internal group_models structure
  group_models <- vector("list", n_groups)
  names(group_models) <- group_levels

  all_X <- list()
  all_y <- list()
  group_vec_parts <- list()

  for (g in group_levels) {
    ex <- extracted[[g]]
    group_models[[g]] <- list(
      coefficients = ex$beta,
      efficiency = ex$te,
      sigma_v = ex$sigma_v,
      sigma_u = ex$sigma_u,
      logLik = ex$logLik,
      hessian = ex$hessian,
      nobs = ex$n,
      X = ex$X,
      y = ex$y,
      dist = ex$dist,
      sign = sign
    )
    all_X[[g]] <- ex$X
    all_y[[g]] <- ex$y
    group_vec_parts[[g]] <- rep(g, ex$n)
  }

  # Combine data across groups
  X_combined <- do.call(rbind, all_X)
  y_combined <- do.call(c, all_y)
  group_vec <- factor(unlist(group_vec_parts), levels = group_levels)
  n <- length(y_combined)
  k <- ncol(X_combined)

  # Compute group frontier and efficiency for all observations
  group_frontier <- numeric(n)
  te_group <- numeric(n)

  for (g in group_levels) {
    idx <- which(group_vec == g)
    beta_g <- group_models[[g]]$coefficients
    group_frontier[idx] <- X_combined[idx, , drop = FALSE] %*% beta_g
    te_group[idx] <- group_models[[g]]$efficiency
  }

  group_coef <- lapply(group_models, function(m) m$coefficients)

  # Estimate the metafrontier
  meta_result <- .meta_stage(X_combined, group_frontier, group_vec,
                             group_levels, group_coef, k, meta_type,
                             "hnormal", control, NULL, objective, sign)

  meta_frontier <- as.numeric(X_combined %*% meta_result$meta_coef)
  tgr <- if (meta_type == "deterministic") {
    pmin(exp(sign * (group_frontier - meta_frontier)), 1.0)
  } else {
    meta_result$tgr
  }
  te_meta <- te_group * tgr

  # Construct a minimal Formula for downstream methods
  input_names <- colnames(X_combined)
  input_names <- input_names[input_names != "(Intercept)"]
  f_str <- paste("y ~", paste(input_names, collapse = " + "))
  f <- Formula::Formula(as.formula(f_str))

  # Build combined data frame
  combined_data <- as.data.frame(X_combined)
  combined_data$y <- y_combined

  list(
    meta_coef = meta_result$meta_coef,
    meta_vcov = meta_result$meta_vcov,
    group_coef = group_coef,
    tgr = tgr,
    te_group = te_group,
    te_meta = te_meta,
    group_frontier = group_frontier,
    meta_frontier = meta_frontier,
    logLik_groups = sapply(group_models, function(m) m$logLik),
    meta_logLik = meta_result$meta_logLik,
    meta_convergence = meta_result$convergence,
    meta_solver = meta_result$meta_solver,
    group_models = group_models,
    groups = group_levels,
    group_vec = group_vec,
    formula = f,
    data = combined_data,
    nobs = c(
      total = n,
      setNames(table(group_vec), group_levels)
    ),
    frontier = if (sign == -1) "cost" else "production"
  )
}


#' Extract components from a pre-fitted frontier model
#' @noRd
.extract_model_components <- function(model) {
  as_metafrontier_model(model)
}


#' Extract from sfaR::sfacross object
#' @noRd
.extract_sfacross <- function(model, estimator = "bc88") {
  if (!requireNamespace("sfaR", quietly = TRUE)) {
    stop("Package 'sfaR' is required to extract from sfacross objects.",
         call. = FALSE)
  }

  # Frontier coefficients (beta only)
  all_coef <- model$mlParam
  if (is.null(all_coef)) all_coef <- coef(model)
  n_beta <- model$nXvar
  if (is.null(n_beta)) {
    n_beta <- length(attr(terms(model$formula), "term.labels")) + 1L
  }
  beta <- all_coef[seq_len(n_beta)]

  # Technical efficiency — efficiencies() may return a data.frame.
  # sfaR's teBC is the BC88 estimator E[exp(-u)|eps]; teJLMS is
  # exp(-E[u|eps]). Pick the column matching the requested estimator.
  te_obj <- sfaR::efficiencies(model)
  te_bc88 <- NULL
  te_jlms <- NULL
  if (is.data.frame(te_obj) || is.matrix(te_obj)) {
    te_obj <- as.data.frame(te_obj)
    if ("teBC" %in% names(te_obj)) te_bc88 <- as.numeric(te_obj[["teBC"]])
    if ("teJLMS" %in% names(te_obj)) te_jlms <- as.numeric(te_obj[["teJLMS"]])
    te <- if (identical(estimator, "jlms") && !is.null(te_jlms)) {
      te_jlms
    } else if (!is.null(te_bc88)) {
      te_bc88
    } else if (!is.null(te_jlms)) {
      te_jlms
    } else {
      as.numeric(te_obj[[1]])
    }
  } else {
    te <- as.numeric(te_obj)
  }

  # Extract X and y from the model's stored data via formula
  dat <- model$dataTable
  if (is.null(dat)) dat <- model$data
  f <- model$formula
  mf <- model.frame(f, data = dat)
  y <- model.response(mf)
  X <- model.matrix(f, data = dat)

  # Sigma parameters
  sigma_v <- NA_real_
  sigma_u <- NA_real_

  # sfaR (all releases, 0.0.91 to 1.0.1) reports log variances
  # ("Zv_(Intercept)", "Zu_(Intercept)")
  # and the truncation mean ("Zmu_(Intercept)"). Use them only when the
  # intercept is the sole term, i.e. without heteroscedasticity or
  # determinants of mu.
  only_term <- function(prefix) {
    terms <- grep(paste0("^", prefix, "_"), names(all_coef), value = TRUE)
    identical(terms, paste0(prefix, "_(Intercept)"))
  }
  if (only_term("Zv")) sigma_v <- exp(all_coef[["Zv_(Intercept)"]] / 2)
  if (only_term("Zu")) sigma_u <- exp(all_coef[["Zu_(Intercept)"]] / 2)
  mu <- if (only_term("Zmu")) all_coef[["Zmu_(Intercept)"]] else NULL

  udist <- model$udist
  if (is.null(udist)) udist <- "hnormal"

  list(
    beta = beta,
    te = te,
    te_bc88 = te_bc88,
    te_jlms = te_jlms,
    X = X,
    y = as.numeric(y),
    sigma_v = as.numeric(sigma_v),
    sigma_u = as.numeric(sigma_u),
    mu = mu,
    logLik = if (is.null(model$mlLoglik)) NA_real_ else model$mlLoglik,
    hessian = model$mlHessian,
    n = length(y),
    dist = switch(udist,
                  "hnormal" = "hnormal",
                  "tnormal" = "tnormal",
                  "exponential" = "exponential",
                  "hnormal"),
    # sfacross fits cost frontiers with S = -1
    sign = if (identical(model$S, -1)) -1 else 1
  )
}


# ---------- Internal: external engines ----------

#' Fit a group SFA frontier via an external engine (sfaR or frontier)
#' @noRd
.fit_sfa_group_external <- function(formula, data, dist, engine,
                                    estimator = "bc88", sign = 1) {
  if (inherits(formula, "Formula") && length(formula)[2] > 1L) {
    stop("External SFA engines support single-part formulas only ",
         "(no inefficiency determinants).", call. = FALSE)
  }
  f_base <- if (inherits(formula, "Formula")) {
    formula(formula, rhs = 1)
  } else {
    formula
  }

  if (engine == "sfaR") {
    if (!requireNamespace("sfaR", quietly = TRUE)) {
      stop("Package 'sfaR' is required for engine = \"sfaR\".",
           call. = FALSE)
    }
    fit <- sfaR::sfacross(formula = f_base, udist = dist, data = data,
                          S = sign)
    ex <- .extract_sfacross(fit, estimator = estimator)
  } else {
    if (!requireNamespace("frontier", quietly = TRUE)) {
      stop("Package 'frontier' is required for engine = \"frontier\".",
           call. = FALSE)
    }
    fit <- frontier::sfa(f_base, data = data, ineffDecrease = sign == 1)
    ex <- .extract_frontier(fit)
    ex$te_bc88 <- ex$te  # frontier::efficiencies() is the BC88-type estimator
    ex$te_jlms <- NULL
    if (identical(estimator, "jlms")) {
      warning("engine = \"frontier\" provides BC88-type efficiencies ",
              "only; 'estimator = \"jlms\"' is not available for this ",
              "engine.", call. = FALSE)
    }
  }

  list(
    coefficients = ex$beta,
    efficiency = ex$te,
    efficiency_bc88 = ex$te_bc88,
    efficiency_jlms = ex$te_jlms,
    estimator = estimator,
    sigma_v = ex$sigma_v,
    sigma_u = ex$sigma_u,
    mu = ex$mu,
    logLik = ex$logLik,
    hessian = ex$hessian,
    nobs = ex$n,
    X = ex$X,
    y = ex$y,
    # Group frontier at the observed inputs, as for internal fits (the
    # parametric bootstrap builds its responses from it)
    fitted = as.numeric(ex$X %*% ex$beta),
    dist = ex$dist,
    sign = sign,
    engine = engine
  )
}


#' Build DEA input/output matrices from a formula and data
#' @noRd
.dea_matrices_mf <- function(formula, data) {
  f_base <- if (inherits(formula, "Formula")) {
    formula(formula, rhs = 1)
  } else {
    formula
  }
  mf <- model.frame(f_base, data = data)
  y_raw <- model.response(mf)
  X_raw <- model.matrix(f_base, data = data)
  if (colnames(X_raw)[1] == "(Intercept)") {
    X_raw <- X_raw[, -1, drop = FALSE]
  }
  Y <- if (is.matrix(y_raw)) y_raw else matrix(y_raw, ncol = 1)
  list(X = X_raw, Y = Y)
}


#' Farrell efficiency via Benchmarking::dea with an external reference set
#' @noRd
.benchmarking_te <- function(X, Y, XREF, YREF, orientation, rts) {
  if (!requireNamespace("Benchmarking", quietly = TRUE)) {
    stop("Package 'Benchmarking' is required for ",
         "engine = \"Benchmarking\".", call. = FALSE)
  }
  ORIENT <- if (orientation == "input") "in" else "out"
  e <- Benchmarking::dea(X = X, Y = Y, RTS = rts, ORIENTATION = ORIENT,
                         XREF = XREF, YREF = YREF)
  scores <- as.numeric(Benchmarking::eff(e))
  # Benchmarking returns the Farrell output measure (>= 1) for "out";
  # convert to TE in (0, 1] to match the internal convention.
  if (ORIENT == "out") 1 / scores else scores
}


#' Fit a group DEA frontier via Benchmarking::dea
#' @noRd
.fit_dea_group_benchmarking <- function(formula, data, orientation, rts) {
  mats <- .dea_matrices_mf(formula, data)
  te <- .benchmarking_te(mats$X, mats$Y, mats$X, mats$Y, orientation, rts)
  list(
    efficiency = te,
    nobs = nrow(mats$X),
    orientation = orientation,
    rts = rts,
    X = mats$X,
    Y = mats$Y,
    engine = "Benchmarking"
  )
}


#' DEA metafrontier via Benchmarking::dea with XREF/YREF pooling
#' @noRd
.estimate_dea_metafrontier_benchmarking <- function(formula, data,
                                                    group_vec,
                                                    group_levels,
                                                    group_models,
                                                    orientation, rts) {
  mats <- .dea_matrices_mf(formula, data)
  n <- nrow(mats$X)

  te_group <- numeric(n)
  for (g in group_levels) {
    idx <- which(group_vec == g)
    te_group[idx] <- group_models[[g]]$efficiency
  }

  # Pooled metafrontier: every observation scored against the full
  # reference set via Benchmarking's XREF/YREF facility.
  te_meta <- .benchmarking_te(mats$X, mats$Y, mats$X, mats$Y,
                              orientation, rts)

  tgr <- te_meta / te_group

  list(
    meta_coef = NULL,
    meta_vcov = NULL,
    group_coef = NULL,
    tgr = tgr,
    te_group = te_group,
    te_meta = te_meta,
    group_frontier = NULL,
    meta_frontier = NULL,
    logLik_groups = NULL,
    meta_logLik = NULL,
    meta_convergence = 0L,
    meta_solver = "Benchmarking"
  )
}
