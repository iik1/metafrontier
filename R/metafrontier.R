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
#' @param panel an optional list with components \code{id} and
#'   \code{time} naming the panel identifier and time columns in
#'   \code{data}. When non-NULL, panel SFA models (BC92/BC95) are
#'   used at the group level.
#' @param panel_dist character. Panel SFA model: \code{"bc92"}
#'   (Battese and Coelli 1992, time-varying inefficiency, default)
#'   or \code{"bc95"} (Battese and Coelli 1995, observation-specific
#'   mean). Only used when \code{panel} is non-NULL.
#' @param type character. For DEA: \code{"radial"} (default) for
#'   standard radial DEA or \code{"directional"} for directional
#'   distance functions.
#' @param direction character. Direction vector for DDF:
#'   \code{"proportional"} (default), \code{"output"}, or
#'   \code{"input"}. Only used when \code{type = "directional"}.
#' @param control a named list of control parameters passed to
#'   \code{\link[stats]{optim}}. Common options include
#'   \code{maxit} (maximum iterations, default 5000),
#'   \code{reltol} (relative convergence tolerance, default 1e-10),
#'   and \code{fnscale} (set to -1 internally for maximisation).
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
#' is estimated by solving a linear program that minimises the total
#' envelope overshoot subject to the constraint that the metafrontier
#' envelops all group frontiers. BRO (2004) originally proposed a
#' constrained least-squares (QP) formulation; the LP yields the
#' tightest envelope and is solved via \pkg{lpSolveAPI}, with a QP
#' fallback via \code{constrOptim()} when the LP is infeasible.
#' The stochastic metafrontier (Huang, Huang, and Liu, 2014) replaces
#' this with a second-stage SFA, providing a distributional framework
#' for inference on the TGR.
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
#' The SFA estimation path assumes a production frontier
#' (\eqn{\varepsilon = v - u}). Cost frontiers (\eqn{\varepsilon = v + u})
#' are not currently supported via the SFA path. The DEA path supports
#' both \code{orientation = "output"} and \code{orientation = "input"}.
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
                         panel = NULL,
                         panel_dist = c("bc92", "bc95"),
                         type = c("radial", "directional"),
                         direction = c("proportional", "output", "input"),
                         control = list(),
                         ...) {

  cl <- match.call()
  method <- match.arg(method)
  meta_type <- match.arg(meta_type)
  dist <- match.arg(dist)
  orientation <- match.arg(orientation)
  rts <- match.arg(rts)
  panel_dist <- match.arg(panel_dist)
  type <- match.arg(type)
  direction <- match.arg(direction)

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
                                  direction = direction, ...)
  } else {
    result <- .estimate_from_models(models, meta_type, control, ...)
  }

  result$call <- cl
  result$method <- method
  result$meta_type <- meta_type
  result$orientation <- orientation
  result$rts <- rts

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
      dots <- list(...)
      if (!is.null(dots$panel)) {
        group_models[[g]] <- .fit_sfa_panel_group(
          f, data_g, dist, dots$panel_dist, dots$panel, control, ...
        )
      } else {
        group_models[[g]] <- .fit_sfa_group(f, data_g, dist, control, ...)
      }
    } else {
      dots <- if (!exists("dots", inherits = FALSE)) list(...) else dots
      if (!is.null(dots$type) && dots$type == "directional") {
        group_models[[g]] <- .fit_ddf_group(f, data_g, rts, dots$direction)
      } else {
        group_models[[g]] <- .fit_dea_group(f, data_g, orientation, rts, ...)
      }
    }
  }

  # Estimate the metafrontier
  if (method == "sfa") {
    meta_result <- .estimate_sfa_metafrontier(
      f, data, group_vec, group_levels, group_models,
      meta_type, dist, control
    )
  } else {
    dots <- list(...)
    if (!is.null(dots$type) && dots$type == "directional") {
      meta_result <- .estimate_ddf_metafrontier(
        f, data, group_vec, group_levels, group_models,
        rts, dots$direction
      )
    } else {
      meta_result <- .estimate_dea_metafrontier(
        f, data, group_vec, group_levels, group_models,
        orientation, rts
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
.estimate_from_models <- function(models, meta_type, control, ...) {
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
      dist = ex$dist
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
  if (meta_type == "deterministic") {
    meta_result <- .deterministic_metafrontier_lp(
      X_combined, group_frontier, group_vec, group_levels, group_coef, k
    )
  } else {
    meta_result <- .stochastic_metafrontier(
      X_combined, group_frontier, group_vec, group_levels, "hnormal", control
    )
  }

  meta_frontier <- as.numeric(X_combined %*% meta_result$meta_coef)
  tgr <- exp(group_frontier - meta_frontier)
  if (meta_type == "deterministic") tgr <- pmin(tgr, 1.0)
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
    group_models = group_models,
    groups = group_levels,
    group_vec = group_vec,
    formula = f,
    data = combined_data,
    nobs = c(
      total = n,
      setNames(table(group_vec), group_levels)
    )
  )
}


#' Extract components from a pre-fitted frontier model
#' @noRd
.extract_model_components <- function(model) {
  as_metafrontier_model(model)
}


#' Extract from sfaR::sfacross object
#' @noRd
.extract_sfacross <- function(model) {
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

  # Technical efficiency — efficiencies() may return a data.frame
  te_obj <- sfaR::efficiencies(model)
  if (is.data.frame(te_obj) || is.matrix(te_obj)) {
    if ("teBC" %in% names(te_obj)) {
      te <- as.numeric(te_obj[["teBC"]])
    } else if ("teJLMS" %in% names(te_obj)) {
      te <- as.numeric(te_obj[["teJLMS"]])
    } else {
      te <- as.numeric(te_obj[[1]])
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
  sv_name <- intersect(c("lnsigmaV", "lnSigmaV"), names(all_coef))
  su_name <- intersect(c("lnsigmaU", "lnSigmaU"), names(all_coef))
  sigma_v <- if (length(sv_name)) exp(all_coef[sv_name[1]]) else NA_real_
  sigma_u <- if (length(su_name)) exp(all_coef[su_name[1]]) else NA_real_

  udist <- model$udist
  if (is.null(udist)) udist <- "hnormal"

  list(
    beta = beta,
    te = te,
    X = X,
    y = as.numeric(y),
    sigma_v = as.numeric(sigma_v),
    sigma_u = as.numeric(sigma_u),
    logLik = if (is.null(model$mlLoglik)) NA_real_ else model$mlLoglik,
    hessian = model$mlHessian,
    n = length(y),
    dist = switch(udist,
                  "hnormal" = "hnormal",
                  "tnormal" = "tnormal",
                  "exponential" = "exponential",
                  "hnormal")
  )
}


#' Extract from frontier::sfa object
#' @noRd
.extract_frontier_sfa <- function(model) {
  if (!requireNamespace("frontier", quietly = TRUE)) {
    stop("Package 'frontier' is required to extract from sfa objects.",
         call. = FALSE)
  }

  # frontier::sfa stores MLE params in mleParam
  all_coef <- model$mleParam
  # Frontier coefficients are named "beta_*" or similar
  beta_idx <- grep("^beta", names(all_coef), ignore.case = TRUE)
  if (length(beta_idx) == 0) {
    # Fallback: OLS params have same structure
    beta_idx <- grep("^beta", names(model$olsParam), ignore.case = TRUE)
    n_beta <- length(beta_idx)
    beta <- all_coef[seq_len(n_beta)]
  } else {
    beta <- all_coef[beta_idx]
  }

  # Efficiency: frontier stores TE in $efficiencies
  te <- as.numeric(frontier::efficiencies(model))

  # Design matrix from the model
  X <- model.matrix(model$formula, data = model$dataTable)
  y <- model.response(model.frame(model$formula, data = model$dataTable))

  # Extract sigma parameters
  sigma_sq <- all_coef["sigmaSq"]
  gamma <- all_coef["gamma"]
  sigma_v <- sqrt(sigma_sq * (1 - gamma))
  sigma_u <- sqrt(sigma_sq * gamma)

  list(
    beta = beta,
    te = as.numeric(te),
    X = X,
    y = as.numeric(y),
    sigma_v = as.numeric(sigma_v),
    sigma_u = as.numeric(sigma_u),
    logLik = model$mleLogl,
    hessian = NULL,
    n = length(y),
    dist = "hnormal"
  )
}

