#' Metafrontier Malmquist Productivity Index
#'
#' Computes the metafrontier Malmquist total factor productivity
#' (TFP) index and its three-way decomposition into technical
#' efficiency change (TEC), technology gap change (TGC), and
#' metafrontier technical change (TC*) for panel data, following
#' O'Donnell, Rao, and Battese (2008).
#'
#' @param formula an object of class \code{\link[Formula]{Formula}}.
#'   Left-hand side specifies the output(s); right-hand side
#'   specifies the inputs. Example: \code{y ~ x1 + x2}.
#' @param data a data frame containing all variables, plus the
#'   grouping and time variables.
#' @param group a character string naming the column in \code{data}
#'   that identifies technology groups, or a vector of group
#'   indicators.
#' @param time a character string naming the column in \code{data}
#'   that identifies time periods, or a vector of time indicators.
#'   Periods must be consecutive integers or sortable.
#' @param id optional. A character string naming the column in
#'   \code{data} that identifies firms across periods, or a vector
#'   of firm identifiers. When supplied, firms are matched across
#'   consecutive periods by identifier within each group. When
#'   \code{NULL} (default), firms are matched by row position within
#'   each group, which is valid only for balanced panels sorted
#'   identically in every period (see Details).
#' @param method character. \code{"dea"} (default) for DEA-based
#'   distance functions or \code{"sfa"} for SFA-based parametric
#'   distance functions (an approximation; see Details).
#' @param dist character. Distribution of the inefficiency term
#'   when \code{method = "sfa"}: \code{"hnormal"} (default),
#'   \code{"tnormal"}, or \code{"exponential"}.
#' @param estimator character. Technical efficiency estimator used
#'   when \code{method = "sfa"}: \code{"bc88"} (default) for the
#'   Battese and Coelli (1988) estimator
#'   \eqn{E[\exp(-u)|\varepsilon]}, or \code{"jlms"} for the Jondrow
#'   et al. (1982) estimator \eqn{\exp(-E[u|\varepsilon])}. Passed to
#'   the group SFA fitter.
#' @param orientation character. \code{"output"} (default) or
#'   \code{"input"}.
#' @param rts character. Returns to scale assumption:
#'   \code{"crs"} (default), \code{"vrs"}, \code{"drs"},
#'   \code{"irs"}, or \code{"fdh"}.
#' @param control a list of control parameters for the SFA optimiser.
#' @param ... additional arguments (currently unused).
#'
#' @return An object of class \code{"malmquist_meta"}, a list
#'   with components:
#'   \describe{
#'     \item{malmquist}{data frame with columns: \code{id},
#'       \code{group}, \code{period_from}, \code{period_to},
#'       \code{MPI} (metafrontier Malmquist TFP index),
#'       \code{TEC} (technical efficiency change),
#'       \code{TGC} (technology gap change),
#'       \code{TC} (metafrontier technical change). The \code{id}
#'       column holds the supplied firm identifiers when \code{id}
#'       is given, and the within-group match position otherwise.}
#'     \item{group_malmquist}{data frame with the within-group
#'       Malmquist index decomposition: \code{MPI_group},
#'       \code{EC_group}, \code{TC_group}}
#'     \item{meta_malmquist}{data frame with the metafrontier
#'       Malmquist index: \code{MPI_meta}, \code{EC_meta},
#'       \code{TC_meta}}
#'     \item{tgr}{data frame with technology gap ratios at each
#'       period endpoint: \code{id}, \code{group},
#'       \code{period_from}, \code{period_to}, \code{TGR_from}
#'       (TGR at the start period), \code{TGR_to} (TGR at the
#'       end period), and \code{TGC} (technology gap change,
#'       \code{TGR_to / TGR_from})}
#'     \item{call}{the matched function call}
#'     \item{method}{the estimation method used (\code{"dea"} or
#'       \code{"sfa"})}
#'     \item{orientation}{the orientation used}
#'     \item{rts}{the returns to scale assumption}
#'     \item{groups}{group labels}
#'     \item{periods}{time periods}
#'     \item{n_infeasible}{total number of infeasible cross-period
#'       DEA programs (always \code{0} for \code{method = "sfa"})}
#'     \item{infeasible_by_period}{data frame with the number of
#'       infeasible cross-period DEA programs per period pair
#'       (\code{method = "dea"} only)}
#'   }
#'
#' @details
#' The metafrontier Malmquist TFP index decomposes productivity
#' change into three components:
#'
#' \deqn{M^* = TEC \times TGC \times TC^*}
#'
#' where:
#' \itemize{
#'   \item \eqn{TEC = TE^{group}_{t+1} / TE^{group}_t}: technical
#'     efficiency change relative to the group frontier
#'   \item \eqn{TGC = TGR_{t+1} / TGR_t}: technology gap change,
#'     capturing whether a group's frontier is catching up to or
#'     falling behind the metafrontier
#'   \item \eqn{TC^*}: metafrontier technical change, measuring
#'     the shift of the global production possibility frontier
#' }
#'
#' \strong{Firm matching:} when \code{id} is supplied, firms are
#' matched across consecutive periods by identifier within each
#' technology group. Duplicated (id, period) combinations within a
#' group are an error. Observations without a within-group match in
#' the adjacent period, either because the panel is unbalanced or
#' because a firm switches group between periods, are dropped, and a
#' single consolidated warning reports the number dropped per period
#' pair. When \code{id} is \code{NULL}, firms are matched by row
#' position within each group; this is valid only for balanced
#' panels sorted identically in every period, so a message is
#' emitted as a reminder, and a warning is issued when group sizes
#' differ across a period pair (the unmatched observations are
#' dropped). Supplying \code{id} is recommended.
#'
#' \strong{DEA-based computation (\code{method = "dea"}):} for each
#' consecutive pair of periods \eqn{(s, t)}, eight sets of LP
#' problems are solved: within-group and pooled efficiencies at each
#' period, plus cross-period evaluations for the geometric mean
#' formulation of technical change. Distances to the metafrontier
#' are exact distances to the pooled-data frontier, as in O'Donnell,
#' Rao and Battese (2008).
#'
#' \strong{SFA-based computation is an approximation
#' (\code{method = "sfa"}):} period-specific group SFA frontiers are
#' estimated, and each observation's metafrontier distance is
#' approximated by the pointwise maximum of the estimated group
#' frontier functions evaluated at its inputs; no enveloping
#' metafrontier is re-estimated. This coincides with the O'Donnell
#' et al. (2008) metafrontier wherever a single group frontier
#' dominates, but can understate the metafrontier where group
#' frontiers cross, which affects TGC and TC*. Prefer
#' \code{method = "dea"} when an exact decomposition is required.
#'
#' \strong{Infeasible cross-period programs:} under
#' \code{rts = "vrs"}, \code{"drs"}, \code{"irs"}, or \code{"fdh"},
#' cross-period LPs can be genuinely infeasible because the
#' reference technology cannot reach the evaluated observation. Such
#' cases yield \code{NA} (never \code{Inf}), are excluded from the
#' reported means, and are counted in a single consolidated warning;
#' the counts are stored in the \code{n_infeasible} and
#' \code{infeasible_by_period} components. \code{rts = "crs"} avoids
#' the issue, as does the hyperbolic orientation available in
#' \code{\link{metafrontier}}.
#'
#' Note that the standard Malmquist index is not a \sQuote{proper}
#' (multiplicatively complete and transitive) TFP index in the sense
#' of O'Donnell (2012), so chained comparisons of index levels
#' across more than two periods should be avoided.
#'
#' @references
#' O'Donnell, C.J., Rao, D.S.P. and Battese, G.E. (2008).
#' Metafrontier frameworks for the study of firm-level efficiencies
#' and technology ratios. \emph{Empirical Economics}, 34(2),
#' 231--255. \doi{10.1007/s00181-007-0119-4}
#'
#' O'Donnell, C.J. (2012). An aggregate quantity framework for
#' measuring and decomposing productivity change.
#' \emph{Journal of Productivity Analysis}, 38(3), 255--272.
#' \doi{10.1007/s11123-012-0275-1}
#'
#' @examples
#' # Simulate panel data for 2 groups, 3 time periods
#' set.seed(42)
#' panels <- lapply(1:3, function(t) {
#'   sim <- simulate_metafrontier(
#'     n_groups = 2, n_per_group = 30,
#'     tech_gap = c(0, 0.3 + 0.05 * t),
#'     sigma_u = c(0.2, 0.3),
#'     seed = 42 + t
#'   )
#'   sim$data$time <- t
#'   sim$data$id <- seq_len(nrow(sim$data))
#'   sim$data
#' })
#' panel_data <- do.call(rbind, panels)
#'
#' # Compute metafrontier Malmquist index, matching firms by id
#' malm <- malmquist_meta(
#'   log_y ~ log_x1 + log_x2,
#'   data = panel_data,
#'   group = "group",
#'   time = "time",
#'   id = "id"
#' )
#' summary(malm)
#'
#' @export
malmquist_meta <- function(formula = NULL,
                           data = NULL,
                           group = NULL,
                           time = NULL,
                           id = NULL,
                           method = c("dea", "sfa"),
                           dist = c("hnormal", "tnormal", "exponential"),
                           estimator = c("bc88", "jlms"),
                           orientation = c("output", "input"),
                           rts = c("crs", "vrs", "drs", "irs", "fdh"),
                           control = list(),
                           ...) {

  cl <- match.call()
  method <- match.arg(method)
  dist <- match.arg(dist)
  estimator <- match.arg(estimator)
  orientation <- match.arg(orientation)
  rts <- match.arg(rts)

  # --- Input validation ---
  if (is.null(formula) || is.null(data) || is.null(group) ||
      is.null(time)) {
    stop("All of 'formula', 'data', 'group', and 'time' are required.",
         call. = FALSE)
  }

  # Parse group variable
  if (is.character(group) && length(group) == 1L) {
    if (!group %in% names(data)) {
      stop("Column '", group, "' not found in data.", call. = FALSE)
    }
    group_vec <- data[[group]]
  } else {
    group_vec <- group
  }
  group_vec <- as.factor(group_vec)

  # Parse time variable
  if (is.character(time) && length(time) == 1L) {
    if (!time %in% names(data)) {
      stop("Column '", time, "' not found in data.", call. = FALSE)
    }
    time_vec <- data[[time]]
  } else {
    time_vec <- time
  }

  # Parse id variable
  if (is.null(id)) {
    id_vec <- NULL
  } else if (is.character(id) && length(id) == 1L) {
    if (!id %in% names(data)) {
      stop("Column '", id, "' not found in data.", call. = FALSE)
    }
    id_vec <- data[[id]]
  } else {
    id_vec <- id
  }
  if (!is.null(id_vec) && length(id_vec) != nrow(data)) {
    stop("'id' must supply one identifier per row of 'data'.",
         call. = FALSE)
  }

  group_levels <- levels(group_vec)
  time_levels <- sort(unique(time_vec))

  if (length(group_levels) < 2L) {
    stop("At least 2 groups are required.", call. = FALSE)
  }
  if (length(time_levels) < 2L) {
    stop("At least 2 time periods are required.", call. = FALSE)
  }

  if (!is.null(id_vec)) {
    for (g in group_levels) {
      for (tt in time_levels) {
        ids_gt <- id_vec[group_vec == g & time_vec == tt]
        if (anyDuplicated(ids_gt)) {
          dup <- unique(ids_gt[duplicated(ids_gt)])
          stop("Duplicated (id, period) combinations in group '", g,
               "' at period ", tt, ": ",
               paste(dup, collapse = ", "), ".", call. = FALSE)
        }
      }
    }
  } else {
    message("No 'id' supplied: firms are matched by row position ",
            "within each group, which is valid only for balanced ",
            "panels sorted identically in every period. Supply 'id' ",
            "to match firms explicitly.")
  }

  # Parse formula
  f <- Formula::Formula(formula)

  # --- SFA Malmquist path ---
  if (method == "sfa") {
    message("Note: method = 'sfa' approximates metafrontier distances ",
            "by the pointwise maximum over estimated group frontiers; ",
            "see ?malmquist_meta for details.")
    result <- .malmquist_sfa(f, data, group_vec, group_levels,
                             time_vec, time_levels, dist, control,
                             id_vec = id_vec, estimator = estimator)
    result$call <- cl
    result$orientation <- orientation
    result$rts <- rts
    result$groups <- group_levels
    result$periods <- time_levels
    result$method <- method
    class(result) <- "malmquist_meta"
    return(result)
  }

  # --- DEA Malmquist path ---
  mf <- model.frame(f, data = data)
  y_raw <- model.response(mf)
  X_raw <- model.matrix(f, data = data, rhs = 1)
  if (colnames(X_raw)[1] == "(Intercept)") {
    X_raw <- X_raw[, -1, drop = FALSE]
  }
  if (!is.matrix(y_raw)) {
    Y_raw <- matrix(y_raw, ncol = 1)
  } else {
    Y_raw <- y_raw
  }

  # --- Solve DEA LPs for all period pairs ---
  # Pre-allocate results list (upper bound: pairs x groups x max_obs)
  n_periods <- length(time_levels) - 1L
  n_groups  <- length(group_levels)
  max_obs_per_group <- max(table(group_vec)) # conservative upper bound per period
  results <- vector("list", n_periods * n_groups * max_obs_per_group)
  result_idx <- 0L

  # Per-period-pair accounting for consolidated warnings
  pair_from <- time_levels[seq_len(n_periods)]
  pair_to   <- time_levels[seq_len(n_periods) + 1L]
  infeas_by_pair  <- integer(n_periods)
  nprog_by_pair   <- integer(n_periods)
  dropped_by_pair <- integer(n_periods)
  n_skipped <- 0L

  for (tp in seq_len(length(time_levels) - 1L)) {
    t_s <- time_levels[tp]
    t_t <- time_levels[tp + 1L]

    idx_s <- which(time_vec == t_s)
    idx_t <- which(time_vec == t_t)

    X_s <- X_raw[idx_s, , drop = FALSE]
    Y_s <- Y_raw[idx_s, , drop = FALSE]
    X_t <- X_raw[idx_t, , drop = FALSE]
    Y_t <- Y_raw[idx_t, , drop = FALSE]

    grp_s <- group_vec[idx_s]
    grp_t <- group_vec[idx_t]

    ids_s <- if (!is.null(id_vec)) id_vec[idx_s] else NULL
    ids_t <- if (!is.null(id_vec)) id_vec[idx_t] else NULL

    n_s <- length(idx_s)
    n_t <- length(idx_t)

    # Storage for group-level DEA scores
    # D_j^s(s): group eff of period-s obs against period-s group tech
    d_grp_ss <- numeric(n_s)
    # D_j^t(t): group eff of period-t obs against period-t group tech
    d_grp_tt <- numeric(n_t)
    # D_j^s(t): period-t obs against period-s group tech
    d_grp_st <- numeric(n_t)
    # D_j^t(s): period-s obs against period-t group tech
    d_grp_ts <- numeric(n_s)

    for (g in group_levels) {
      idx_gs <- which(grp_s == g)
      idx_gt <- which(grp_t == g)

      if (length(idx_gs) == 0 || length(idx_gt) == 0) next

      X_gs <- X_s[idx_gs, , drop = FALSE]
      Y_gs <- Y_s[idx_gs, , drop = FALSE]
      X_gt <- X_t[idx_gt, , drop = FALSE]
      Y_gt <- Y_t[idx_gt, , drop = FALSE]

      # D_j^s(x_s, y_s) — same-period group efficiency
      d_grp_ss[idx_gs] <- .muffle_lp_warnings(.dea_batch_fast(
        X_gs, Y_gs, orientation, rts))

      # D_j^t(x_t, y_t) — same-period group efficiency
      d_grp_tt[idx_gt] <- .muffle_lp_warnings(.dea_batch_fast(
        X_gt, Y_gt, orientation, rts))

      # D_j^s(x_t, y_t) — cross-period: period-t obs against period-s group tech
      d_grp_st[idx_gt] <- .muffle_lp_warnings(.dea_batch_fast(
        X_gt, Y_gt, orientation, rts,
        X_ref = X_gs, Y_ref = Y_gs))

      # D_j^t(x_s, y_s) — cross-period: period-s obs against period-t group tech
      d_grp_ts[idx_gs] <- .muffle_lp_warnings(.dea_batch_fast(
        X_gs, Y_gs, orientation, rts,
        X_ref = X_gt, Y_ref = Y_gt))

      nprog_by_pair[tp] <- nprog_by_pair[tp] +
        length(idx_gt) + length(idx_gs)
      infeas_by_pair[tp] <- infeas_by_pair[tp] +
        sum(is.na(d_grp_st[idx_gt])) + sum(is.na(d_grp_ts[idx_gs]))
    }

    # Metafrontier (pooled) DEA scores. Infeasible cross-period LPs
    # yield NA; they are counted here and reported once per call.
    d_meta_ss <- .muffle_lp_warnings(
      .dea_batch(X_s, Y_s, X_s, Y_s, orientation, rts))
    d_meta_tt <- .muffle_lp_warnings(
      .dea_batch(X_t, Y_t, X_t, Y_t, orientation, rts))
    d_meta_st <- .muffle_lp_warnings(
      .dea_batch(X_t, Y_t, X_s, Y_s, orientation, rts))
    d_meta_ts <- .muffle_lp_warnings(
      .dea_batch(X_s, Y_s, X_t, Y_t, orientation, rts))

    nprog_by_pair[tp] <- nprog_by_pair[tp] + n_t + n_s
    infeas_by_pair[tp] <- infeas_by_pair[tp] +
      sum(is.na(d_meta_st)) + sum(is.na(d_meta_ts))

    # --- Identify matched firms across periods ---
    for (g in group_levels) {
      idx_gs <- which(grp_s == g)
      idx_gt <- which(grp_t == g)

      if (!is.null(id_vec)) {
        ids_gs <- ids_s[idx_gs]
        ids_gt <- ids_t[idx_gt]
        common <- intersect(ids_gs, ids_gt)
        dropped_by_pair[tp] <- dropped_by_pair[tp] +
          (length(ids_gs) - length(common)) +
          (length(ids_gt) - length(common))
        pos_s <- idx_gs[match(common, ids_gs)]
        pos_t <- idx_gt[match(common, ids_gt)]
        out_ids <- common
      } else {
        n_matched <- min(length(idx_gs), length(idx_gt))
        dropped_by_pair[tp] <- dropped_by_pair[tp] +
          abs(length(idx_gs) - length(idx_gt))
        pos_s <- idx_gs[seq_len(n_matched)]
        pos_t <- idx_gt[seq_len(n_matched)]
        out_ids <- seq_len(n_matched)
      }

      if (length(pos_s) == 0L) next

      for (i in seq_along(pos_s)) {
        is <- pos_s[i]
        it <- pos_t[i]

        # Guard against zero/NA
        if (is.na(d_grp_ss[is]) || is.na(d_grp_tt[it]) ||
            d_grp_ss[is] <= 0 || d_grp_tt[it] <= 0) {
          n_skipped <- n_skipped + 1L
          next
        }

        # -- Within-group Malmquist --
        ec_grp <- d_grp_tt[it] / d_grp_ss[is]

        tc_grp_1 <- if (d_grp_tt[it] > 0 && !is.na(d_grp_st[it]))
          d_grp_st[it] / d_grp_tt[it] else NA_real_
        tc_grp_2 <- if (d_grp_ts[is] > 0 && !is.na(d_grp_ts[is]))
          d_grp_ss[is] / d_grp_ts[is] else NA_real_

        tc_grp <- if (!is.na(tc_grp_1) && !is.na(tc_grp_2) &&
                      tc_grp_1 > 0 && tc_grp_2 > 0)
          sqrt(tc_grp_1 * tc_grp_2) else NA_real_

        mpi_grp <- ec_grp * tc_grp

        # -- Metafrontier Malmquist --
        ec_meta <- d_meta_tt[it] / d_meta_ss[is]

        tc_meta_1 <- if (d_meta_tt[it] > 0 && !is.na(d_meta_st[it]))
          d_meta_st[it] / d_meta_tt[it] else NA_real_
        tc_meta_2 <- if (d_meta_ts[is] > 0 && !is.na(d_meta_ts[is]))
          d_meta_ss[is] / d_meta_ts[is] else NA_real_

        tc_meta <- if (!is.na(tc_meta_1) && !is.na(tc_meta_2) &&
                       tc_meta_1 > 0 && tc_meta_2 > 0)
          sqrt(tc_meta_1 * tc_meta_2) else NA_real_

        mpi_meta <- ec_meta * tc_meta

        # -- Three-way decomposition (ORB 2008) --
        # TEC = within-group efficiency change
        tec <- ec_grp

        # TGR in [0,1]: Farrell TE ratio for DEA, distance function ratio for SFA
        # Both conventions yield TGR = TE_meta/TE_group <= 1 (BRO 2004 definition)
        # DEA path: d_meta and d_grp are Farrell efficiencies (<=1), so
        # TGR = TE_meta / TE_group <= 1 because TE_meta <= TE_group.
        tgr_s <- d_meta_ss[is] / d_grp_ss[is]
        tgr_t <- d_meta_tt[it] / d_grp_tt[it]

        # TGC = technology gap change
        tgc <- tgr_t / tgr_s

        # TC* = M* / (TEC * TGC) = metafrontier technical change
        tc_star <- if (!is.na(mpi_meta) && tec > 0 && !is.na(tgc) &&
                       tgc > 0)
          mpi_meta / (tec * tgc) else NA_real_

        result_idx <- result_idx + 1L
        results[[result_idx]] <- data.frame(
          id = out_ids[i],
          group = g,
          period_from = t_s,
          period_to = t_t,
          # Three-way metafrontier decomposition
          MPI = mpi_meta,
          TEC = tec,
          TGC = tgc,
          TC = tc_star,
          # Within-group decomposition
          MPI_group = mpi_grp,
          EC_group = ec_grp,
          TC_group = tc_grp,
          # Metafrontier decomposition
          MPI_meta = mpi_meta,
          EC_meta = ec_meta,
          TC_meta = tc_meta,
          # TGR at each period
          TGR_from = tgr_s,
          TGR_to = tgr_t,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (result_idx == 0L) {
    stop("No matched observations found across time periods.",
         call. = FALSE)
  }

  .warn_dropped(dropped_by_pair, pair_from, pair_to, is.null(id_vec))

  total_infeasible <- sum(infeas_by_pair)
  if (total_infeasible > 0L || n_skipped > 0L) {
    msg <- character(0)
    if (total_infeasible > 0L) {
      msg <- paste0(total_infeasible, " of ", sum(nprog_by_pair),
                    " cross-period DEA programs were infeasible (rts = \"",
                    rts, "\"); the affected TC and MPI values are NA and ",
                    "are excluded from reported means.")
    }
    if (n_skipped > 0L) {
      msg <- c(msg, paste0(n_skipped, " matched observation(s) were ",
                           "skipped because of missing or non-positive ",
                           "same-period distance scores."))
    }
    warning(paste(msg, collapse = " "), call. = FALSE)
  }

  malmquist_df <- do.call(rbind, results[seq_len(result_idx)])
  rownames(malmquist_df) <- NULL

  out <- list(
    malmquist = malmquist_df[, c("id", "group", "period_from",
                                  "period_to", "MPI", "TEC",
                                  "TGC", "TC")],
    group_malmquist = malmquist_df[, c("id", "group", "period_from",
                                        "period_to", "MPI_group",
                                        "EC_group", "TC_group")],
    meta_malmquist = malmquist_df[, c("id", "group", "period_from",
                                       "period_to", "MPI_meta",
                                       "EC_meta", "TC_meta")],
    tgr = malmquist_df[, c("id", "group", "period_from",
                            "period_to", "TGR_from", "TGR_to",
                            "TGC")],
    call = cl,
    method = method,
    orientation = orientation,
    rts = rts,
    groups = group_levels,
    periods = time_levels,
    n_infeasible = total_infeasible,
    infeasible_by_period = data.frame(
      period_from = pair_from,
      period_to = pair_to,
      n_infeasible = infeas_by_pair,
      n_programs = nprog_by_pair,
      stringsAsFactors = FALSE
    )
  )
  class(out) <- "malmquist_meta"
  out
}


#' Batch DEA solver
#'
#' Evaluates all rows of X_eval/Y_eval against the reference
#' technology X_ref/Y_ref.  Delegates to \code{.dea_batch_fast()}
#' which reuses a single LP object for performance.
#'
#' @keywords internal
#' @noRd
.dea_batch <- function(X_eval, Y_eval, X_ref, Y_ref,
                       orientation, rts) {
  .dea_batch_fast(X_eval, Y_eval, orientation, rts,
                  X_ref = X_ref, Y_ref = Y_ref)
}


#' Muffle per-DMU LP infeasibility warnings
#'
#' Infeasible programs return NA scores; they are counted by the
#' caller and reported in one consolidated warning per call instead
#' of one warning per DMU. Other warnings pass through untouched.
#'
#' @keywords internal
#' @noRd
.muffle_lp_warnings <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("LP infeasible|No dominating FDH reference point",
              conditionMessage(w))) {
      invokeRestart("muffleWarning")
    }
  })
}


#' Consolidated warning for observations dropped during matching
#'
#' @keywords internal
#' @noRd
.warn_dropped <- function(dropped_by_pair, pair_from, pair_to,
                          positional) {
  total_dropped <- sum(dropped_by_pair)
  if (total_dropped == 0L) return(invisible(NULL))
  nz <- dropped_by_pair > 0L
  detail <- paste0(dropped_by_pair[nz], " in ", pair_from[nz],
                   " -> ", pair_to[nz], collapse = "; ")
  if (positional) {
    warning("Group sizes differ across periods: positional matching ",
            "dropped ", total_dropped, " observation(s) (", detail,
            "). Supply 'id' to match firms explicitly.", call. = FALSE)
  } else {
    warning(total_dropped, " observation(s) could not be matched ",
            "across periods within their group and were dropped (",
            detail, ").", call. = FALSE)
  }
  invisible(NULL)
}


#' @export
print.malmquist_meta <- function(x, ...) {
  cat("\nMetafrontier Malmquist TFP Index\n")
  cat("================================\n")
  if (!is.null(x$method)) {
    cat("Method:      ", x$method, "\n")
  }
  cat("Orientation: ", x$orientation, "\n")
  cat("RTS:         ", x$rts, "\n")
  cat("Groups:      ", paste(x$groups, collapse = ", "), "\n")
  cat("Periods:     ", paste(x$periods, collapse = " -> "), "\n")
  cat("Observations:", nrow(x$malmquist), "\n\n")

  cat("Mean decomposition (M* = TEC x TGC x TC*):\n")
  means <- colMeans(x$malmquist[, c("MPI", "TEC", "TGC", "TC")],
                    na.rm = TRUE)
  cat("  MPI  =", format(means["MPI"], digits = 4), "\n")
  cat("  TEC  =", format(means["TEC"], digits = 4), "\n")
  cat("  TGC  =", format(means["TGC"], digits = 4), "\n")
  cat("  TC*  =", format(means["TC"], digits = 4), "\n")

  if (!is.null(x$n_infeasible) && x$n_infeasible > 0) {
    cat("\nInfeasible cross-period DEA programs:", x$n_infeasible,
        "(affected values are NA and excluded from means)\n")
  }
  if (!is.null(x$method) && x$method == "sfa") {
    cat("\nNote: SFA metafrontier distances are pointwise-maximum",
        "approximations\nover estimated group frontiers (see",
        "?malmquist_meta).\n")
  }

  invisible(x)
}


#' @export
summary.malmquist_meta <- function(object, ...) {
  decomp_cols <- c("MPI", "TEC", "TGC", "TC")

  # Compute per-group means
  by_group <- lapply(object$groups, function(g) {
    idx <- object$malmquist$group == g
    n_g <- sum(idx)
    means <- colMeans(object$malmquist[idx, decomp_cols, drop = FALSE],
                      na.rm = TRUE)
    list(n = n_g, means = means)
  })
  names(by_group) <- object$groups

  # Compute per-period means
  periods_from <- unique(object$malmquist$period_from)
  by_period <- lapply(periods_from, function(p) {
    idx <- object$malmquist$period_from == p
    period_to <- unique(object$malmquist$period_to[idx])
    means <- colMeans(object$malmquist[idx, decomp_cols, drop = FALSE],
                      na.rm = TRUE)
    list(period_from = p, period_to = period_to, means = means)
  })
  names(by_period) <- as.character(periods_from)

  # Compute per-group TGR summary
  tgr_summary <- lapply(object$groups, function(g) {
    idx <- object$tgr$group == g
    sub <- object$tgr[idx, ]
    c(TGR_from = mean(sub$TGR_from, na.rm = TRUE),
      TGR_to   = mean(sub$TGR_to,   na.rm = TRUE),
      TGC      = mean(sub$TGC,      na.rm = TRUE))
  })
  names(tgr_summary) <- object$groups

  out <- list(
    call        = object$call,
    method      = object$method,
    orientation = object$orientation,
    rts         = object$rts,
    groups      = object$groups,
    periods     = object$periods,
    n_obs       = nrow(object$malmquist),
    overall     = colMeans(object$malmquist[, decomp_cols, drop = FALSE],
                           na.rm = TRUE),
    by_group    = by_group,
    by_period   = by_period,
    tgr_summary = tgr_summary,
    n_infeasible = object$n_infeasible,
    infeasible_by_period = object$infeasible_by_period
  )
  class(out) <- "summary.malmquist_meta"
  out
}


#' @export
print.summary.malmquist_meta <- function(x, ...) {
  cat("\nMetafrontier Malmquist TFP Index Summary\n")
  cat("=========================================\n\n")
  cat("Call:\n")
  print(x$call)
  if (!is.null(x$method)) {
    cat("\nMethod:     ", x$method, "\n")
    cat("Orientation:", x$orientation, "\n")
  } else {
    cat("\nOrientation:", x$orientation, "\n")
  }
  cat("RTS:        ", x$rts, "\n")
  cat("Groups:     ", paste(x$groups, collapse = ", "), "\n")
  cat("Periods:    ", paste(x$periods, collapse = " -> "), "\n")
  cat("Observations:", x$n_obs, "\n")

  if (!is.null(x$n_infeasible) && x$n_infeasible > 0) {
    cat("Infeasible cross-period DEA programs:", x$n_infeasible,
        "(values NA, excluded from means)\n")
    if (!is.null(x$infeasible_by_period)) {
      tab <- x$infeasible_by_period
      tab <- tab[tab$n_infeasible > 0, , drop = FALSE]
      for (k in seq_len(nrow(tab))) {
        cat("  ", tab$period_from[k], "->", tab$period_to[k], ":",
            tab$n_infeasible[k], "of", tab$n_programs[k], "programs\n")
      }
    }
  }
  cat("\n")

  # Overall means
  cat("Overall means:\n")
  cat("  MPI  =", format(x$overall["MPI"], digits = 4), "\n")
  cat("  TEC  =", format(x$overall["TEC"], digits = 4), "\n")
  cat("  TGC  =", format(x$overall["TGC"], digits = 4), "\n")
  cat("  TC*  =", format(x$overall["TC"],  digits = 4), "\n\n")

  # By-group means
  cat("--- Three-Way Decomposition by Group ---\n")
  cat("M* = TEC x TGC x TC*\n\n")

  for (g in x$groups) {
    bg <- x$by_group[[g]]
    cat("Group:", g, "(n =", bg$n, ")\n")
    print(round(bg$means, 4))
    cat("\n")
  }

  # By-period means
  if (length(x$by_period) > 1L) {
    cat("--- By Period ---\n\n")
    for (bp in x$by_period) {
      cat("Period", bp$period_from, "->", bp$period_to, "\n")
      print(round(bp$means, 4))
      cat("\n")
    }
  }

  # Overall TGR change
  cat("--- Technology Gap Ratios ---\n\n")
  for (g in x$groups) {
    tgr <- x$tgr_summary[[g]]
    cat("Group:", g, "\n")
    cat("  Mean TGR (from):", format(tgr["TGR_from"], digits = 4), "\n")
    cat("  Mean TGR (to):  ", format(tgr["TGR_to"],   digits = 4), "\n")
    cat("  Mean TGC:       ", format(tgr["TGC"],       digits = 4), "\n\n")
  }

  invisible(x)
}


# ---------- SFA Malmquist ----------

#' SFA-based Malmquist decomposition
#'
#' Fits period-specific group SFA models and computes the three-way
#' decomposition using SFA distance functions. The metafrontier is
#' approximated by the pointwise maximum over the estimated group
#' frontiers (see the Details section of \code{malmquist_meta}).
#'
#' @keywords internal
#' @noRd
.malmquist_sfa <- function(formula, data, group_vec, group_levels,
                           time_vec, time_levels, dist, control,
                           id_vec = NULL, estimator = "bc88") {

  # Fit period x group SFA models
  models <- list()
  for (tt in time_levels) {
    models[[as.character(tt)]] <- list()
    for (g in group_levels) {
      idx <- which(group_vec == g & time_vec == tt)
      if (length(idx) < 5) {
        warning("Group '", g, "' at time ", tt, " has ", length(idx),
                " obs. Skipping.", call. = FALSE)
        next
      }
      data_gt <- data[idx, , drop = FALSE]
      models[[as.character(tt)]][[g]] <- tryCatch(
        .fit_sfa_group(formula, data_gt, dist, control,
                       estimator = estimator),
        error = function(e) NULL
      )
    }
  }

  # Build design matrices
  mf <- model.frame(formula(formula, rhs = 1), data = data,
                     na.action = na.omit)
  X_all <- model.matrix(formula(formula, rhs = 1), data = mf)
  y_all <- model.response(mf)

  # Pre-allocate results list (upper bound: pairs x groups x max_obs)
  sfa_n_periods <- length(time_levels) - 1L
  sfa_n_groups  <- length(group_levels)
  sfa_max_obs   <- max(table(group_vec))
  results <- vector("list", sfa_n_periods * sfa_n_groups * sfa_max_obs)
  result_idx <- 0L

  # Per-period-pair accounting for consolidated warnings
  pair_from <- time_levels[seq_len(sfa_n_periods)]
  pair_to   <- time_levels[seq_len(sfa_n_periods) + 1L]
  dropped_by_pair <- integer(sfa_n_periods)
  n_skipped <- 0L

  for (tp in seq_len(length(time_levels) - 1L)) {
    t_s <- time_levels[tp]
    t_t <- time_levels[tp + 1L]
    ts_chr <- as.character(t_s)
    tt_chr <- as.character(t_t)

    for (g in group_levels) {
      mod_s <- models[[ts_chr]][[g]]
      mod_t <- models[[tt_chr]][[g]]
      if (is.null(mod_s) || is.null(mod_t)) next

      beta_s <- mod_s$coefficients
      beta_t <- mod_t$coefficients

      # Obs indices for each period
      idx_s <- which(group_vec == g & time_vec == t_s)
      idx_t <- which(group_vec == g & time_vec == t_t)

      if (!is.null(id_vec)) {
        ids_gs <- id_vec[idx_s]
        ids_gt <- id_vec[idx_t]
        common <- intersect(ids_gs, ids_gt)
        dropped_by_pair[tp] <- dropped_by_pair[tp] +
          (length(ids_gs) - length(common)) +
          (length(ids_gt) - length(common))
        pos_s <- idx_s[match(common, ids_gs)]
        pos_t <- idx_t[match(common, ids_gt)]
        out_ids <- common
      } else {
        n_matched <- min(length(idx_s), length(idx_t))
        dropped_by_pair[tp] <- dropped_by_pair[tp] +
          abs(length(idx_s) - length(idx_t))
        pos_s <- idx_s[seq_len(n_matched)]
        pos_t <- idx_t[seq_len(n_matched)]
        out_ids <- seq_len(n_matched)
      }

      if (length(pos_s) == 0L) next

      for (i in seq_along(pos_s)) {
        is <- pos_s[i]
        it <- pos_t[i]

        x_s <- X_all[is, ]
        x_t <- X_all[it, ]
        y_s <- y_all[is]
        y_t <- y_all[it]

        # SFA distance: D(x,y|beta) = exp(x'beta - ln_y)
        # = exp(x'beta) / exp(ln_y) = frontier(x) / y
        # For log output, ln D = x'beta - y (since y is already log)
        d_grp_ss <- exp(sum(x_s * beta_s) - y_s)  # same-period
        d_grp_tt <- exp(sum(x_t * beta_t) - y_t)

        # Cross-period
        d_grp_st <- exp(sum(x_t * beta_s) - y_t)  # t obs, s tech
        d_grp_ts <- exp(sum(x_s * beta_t) - y_s)  # s obs, t tech

        # Group TE
        te_s <- mod_s$efficiency[match(is, idx_s)]
        te_t <- mod_t$efficiency[match(it, idx_t)]

        if (is.na(te_s) || is.na(te_t) || te_s <= 0 || te_t <= 0) {
          n_skipped <- n_skipped + 1L
          next
        }

        # TEC = TE_t / TE_s
        tec <- te_t / te_s

        # Group TC (geometric mean)
        tc_grp_1 <- d_grp_st / d_grp_tt
        tc_grp_2 <- d_grp_ss / d_grp_ts
        tc_grp <- sqrt(tc_grp_1 * tc_grp_2)

        ec_grp <- tec
        mpi_grp <- ec_grp * tc_grp

        # Pooled (meta) distances — use average coefficients across groups
        # as a simple pooled approximation
        all_betas_s <- lapply(group_levels, function(gg) {
          m <- models[[ts_chr]][[gg]]
          if (!is.null(m)) m$coefficients else NULL
        })
        all_betas_t <- lapply(group_levels, function(gg) {
          m <- models[[tt_chr]][[gg]]
          if (!is.null(m)) m$coefficients else NULL
        })
        all_betas_s <- all_betas_s[!sapply(all_betas_s, is.null)]
        all_betas_t <- all_betas_t[!sapply(all_betas_t, is.null)]

        # Metafrontier approximation: envelope of group-specific frontier predictions.
        # This is the pointwise maximum over group distance functions, which
        # approximates the O'Donnell-Rao-Battese (2008) meta-technology as the
        # dominant group technology at each input mix. For non-parallel frontiers
        # with 3+ groups, this may differ from the convex hull definition.
        meta_ss <- max(sapply(all_betas_s, function(b) exp(sum(x_s * b) - y_s)))
        meta_tt <- max(sapply(all_betas_t, function(b) exp(sum(x_t * b) - y_t)))
        meta_st <- max(sapply(all_betas_s, function(b) exp(sum(x_t * b) - y_t)))
        meta_ts <- max(sapply(all_betas_t, function(b) exp(sum(x_s * b) - y_s)))

        # TGR in [0,1]: Farrell TE ratio for DEA, distance function ratio for SFA
        # Both conventions yield TGR = TE_meta/TE_group <= 1 (BRO 2004 definition)
        # SFA path: d_grp and meta are Shephard distance function values (>=1), so
        # TGR = D_group / D_meta <= 1 because D_meta >= D_group.
        tgr_s <- d_grp_ss / meta_ss
        tgr_t <- d_grp_tt / meta_tt
        tgc <- tgr_t / tgr_s

        ec_meta <- meta_tt * te_t / (meta_ss * te_s)
        tc_meta <- sqrt((meta_st / meta_tt) * (meta_ss / meta_ts))
        mpi_meta <- ec_meta * tc_meta

        tc_star <- if (tec > 0 && tgc > 0 && is.finite(mpi_meta))
          mpi_meta / (tec * tgc) else NA_real_

        result_idx <- result_idx + 1L
        results[[result_idx]] <- data.frame(
          id = out_ids[i],
          group = g,
          period_from = t_s,
          period_to = t_t,
          MPI = mpi_meta,
          TEC = tec,
          TGC = tgc,
          TC = tc_star,
          MPI_group = mpi_grp,
          EC_group = ec_grp,
          TC_group = tc_grp,
          MPI_meta = mpi_meta,
          EC_meta = ec_meta,
          TC_meta = tc_meta,
          TGR_from = tgr_s,
          TGR_to = tgr_t,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (result_idx == 0L) {
    stop("No matched observations found for SFA Malmquist.",
         call. = FALSE)
  }

  .warn_dropped(dropped_by_pair, pair_from, pair_to, is.null(id_vec))

  if (n_skipped > 0L) {
    warning(n_skipped, " matched observation(s) were skipped because ",
            "of missing or non-positive efficiency scores.",
            call. = FALSE)
  }

  malmquist_df <- do.call(rbind, results[seq_len(result_idx)])
  rownames(malmquist_df) <- NULL

  list(
    malmquist = malmquist_df[, c("id", "group", "period_from",
                                  "period_to", "MPI", "TEC",
                                  "TGC", "TC")],
    group_malmquist = malmquist_df[, c("id", "group", "period_from",
                                        "period_to", "MPI_group",
                                        "EC_group", "TC_group")],
    meta_malmquist = malmquist_df[, c("id", "group", "period_from",
                                       "period_to", "MPI_meta",
                                       "EC_meta", "TC_meta")],
    tgr = malmquist_df[, c("id", "group", "period_from",
                            "period_to", "TGR_from", "TGR_to",
                            "TGC")],
    n_infeasible = 0L
  )
}
