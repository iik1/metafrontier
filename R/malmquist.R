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
#' @param orientation character. \code{"output"} (default) or
#'   \code{"input"}.
#' @param rts character. Returns to scale assumption:
#'   \code{"crs"} (default), \code{"vrs"}, \code{"drs"}, or
#'   \code{"irs"}.
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
#'       \code{TC} (metafrontier technical change)}
#'     \item{group_malmquist}{data frame with the within-group
#'       Malmquist index decomposition: \code{MPI_group},
#'       \code{EC_group}, \code{TC_group}}
#'     \item{meta_malmquist}{data frame with the metafrontier
#'       Malmquist index: \code{MPI_meta}, \code{EC_meta},
#'       \code{TC_meta}}
#'     \item{call}{the matched function call}
#'     \item{orientation}{the orientation used}
#'     \item{rts}{the returns to scale assumption}
#'     \item{groups}{group labels}
#'     \item{periods}{time periods}
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
#' Computation uses DEA-based distance functions. For each
#' consecutive pair of periods \eqn{(s, t)}, eight sets of LP
#' problems are solved: within-group and pooled efficiencies at
#' each period, plus cross-period evaluations for the geometric
#' mean formulation of technical change.
#'
#' @references
#' O'Donnell, C.J., Rao, D.S.P. and Battese, G.E. (2008).
#' Metafrontier frameworks for the study of firm-level efficiencies
#' and technology ratios. \emph{Empirical Economics}, 34(2),
#' 231--255. \doi{10.1007/s00181-006-0095-0}
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
#' # Compute metafrontier Malmquist index
#' malm <- malmquist_meta(
#'   log_y ~ log_x1 + log_x2,
#'   data = panel_data,
#'   group = "group",
#'   time = "time"
#' )
#' summary(malm)
#'
#' @export
malmquist_meta <- function(formula = NULL,
                           data = NULL,
                           group = NULL,
                           time = NULL,
                           orientation = c("output", "input"),
                           rts = c("crs", "vrs", "drs", "irs"),
                           ...) {

  cl <- match.call()
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

  group_levels <- levels(group_vec)
  time_levels <- sort(unique(time_vec))

  if (length(group_levels) < 2L) {
    stop("At least 2 groups are required.", call. = FALSE)
  }
  if (length(time_levels) < 2L) {
    stop("At least 2 time periods are required.", call. = FALSE)
  }

  # Parse formula and extract data matrices
  f <- Formula::Formula(formula)
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
  results <- list()

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
      for (i in seq_along(idx_gs)) {
        d_grp_ss[idx_gs[i]] <- .dea_solve_lp(
          X_gs[i, ], Y_gs[i, ], X_gs, Y_gs, orientation, rts
        )
      }

      # D_j^t(x_t, y_t) — same-period group efficiency
      for (i in seq_along(idx_gt)) {
        d_grp_tt[idx_gt[i]] <- .dea_solve_lp(
          X_gt[i, ], Y_gt[i, ], X_gt, Y_gt, orientation, rts
        )
      }

      # D_j^s(x_t, y_t) — period-t obs against period-s group tech
      for (i in seq_along(idx_gt)) {
        d_grp_st[idx_gt[i]] <- .dea_solve_lp(
          X_gt[i, ], Y_gt[i, ], X_gs, Y_gs, orientation, rts
        )
      }

      # D_j^t(x_s, y_s) — period-s obs against period-t group tech
      for (i in seq_along(idx_gs)) {
        d_grp_ts[idx_gs[i]] <- .dea_solve_lp(
          X_gs[i, ], Y_gs[i, ], X_gt, Y_gt, orientation, rts
        )
      }
    }

    # Metafrontier (pooled) DEA scores
    # D*^s(s): period-s obs against period-s pooled tech
    d_meta_ss <- .dea_batch(X_s, Y_s, X_s, Y_s, orientation, rts)
    # D*^t(t): period-t obs against period-t pooled tech
    d_meta_tt <- .dea_batch(X_t, Y_t, X_t, Y_t, orientation, rts)
    # D*^s(t): period-t obs against period-s pooled tech
    d_meta_st <- .dea_batch(X_t, Y_t, X_s, Y_s, orientation, rts)
    # D*^t(s): period-s obs against period-t pooled tech
    d_meta_ts <- .dea_batch(X_s, Y_s, X_t, Y_t, orientation, rts)

    # --- Identify matched firms across periods ---
    # We need firms present in both periods. Match by group.
    for (g in group_levels) {
      idx_gs <- which(grp_s == g)
      idx_gt <- which(grp_t == g)
      n_pairs <- min(length(idx_gs), length(idx_gt))

      if (n_pairs == 0) next

      for (i in seq_len(n_pairs)) {
        is <- idx_gs[i]
        it <- idx_gt[i]

        # Guard against zero/NA
        if (is.na(d_grp_ss[is]) || is.na(d_grp_tt[it]) ||
            d_grp_ss[is] <= 0 || d_grp_tt[it] <= 0) next

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

        # TGR at each period
        tgr_s <- d_meta_ss[is] / d_grp_ss[is]
        tgr_t <- d_meta_tt[it] / d_grp_tt[it]

        # TGC = technology gap change
        tgc <- tgr_t / tgr_s

        # TC* = M* / (TEC * TGC) = metafrontier technical change
        tc_star <- if (!is.na(mpi_meta) && tec > 0 && !is.na(tgc) &&
                       tgc > 0)
          mpi_meta / (tec * tgc) else NA_real_

        results[[length(results) + 1L]] <- data.frame(
          id = i,
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

  if (length(results) == 0) {
    stop("No matched observations found across time periods.",
         call. = FALSE)
  }

  malmquist_df <- do.call(rbind, results)
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
    orientation = orientation,
    rts = rts,
    groups = group_levels,
    periods = time_levels
  )
  class(out) <- "malmquist_meta"
  out
}


#' Batch DEA solver
#'
#' Evaluates all rows of X_eval/Y_eval against the reference
#' technology X_ref/Y_ref.
#'
#' @keywords internal
#' @noRd
.dea_batch <- function(X_eval, Y_eval, X_ref, Y_ref,
                       orientation, rts) {
  n <- nrow(X_eval)
  eff <- numeric(n)
  for (i in seq_len(n)) {
    eff[i] <- .dea_solve_lp(
      x_i = X_eval[i, ],
      y_i = Y_eval[i, ],
      X = X_ref,
      Y = Y_ref,
      orientation = orientation,
      rts = rts
    )
  }
  eff
}


#' @export
print.malmquist_meta <- function(x, ...) {
  cat("\nMetafrontier Malmquist TFP Index\n")
  cat("================================\n")
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

  invisible(x)
}


#' @export
summary.malmquist_meta <- function(object, ...) {
  cat("\nMetafrontier Malmquist TFP Index Summary\n")
  cat("=========================================\n\n")
  cat("Call:\n")
  print(object$call)
  cat("\nOrientation:", object$orientation, "\n")
  cat("RTS:        ", object$rts, "\n")
  cat("Groups:     ", paste(object$groups, collapse = ", "), "\n")
  cat("Periods:    ", paste(object$periods, collapse = " -> "), "\n\n")

  # By-group means
  cat("--- Three-Way Decomposition by Group ---\n")
  cat("M* = TEC x TGC x TC*\n\n")

  for (g in object$groups) {
    idx <- object$malmquist$group == g
    sub <- object$malmquist[idx, c("MPI", "TEC", "TGC", "TC")]
    cat("Group:", g, "(n =", sum(idx), ")\n")
    print(round(colMeans(sub, na.rm = TRUE), 4))
    cat("\n")
  }

  # By-period means
  periods_from <- unique(object$malmquist$period_from)
  if (length(periods_from) > 1L) {
    cat("--- By Period ---\n\n")
    for (p in periods_from) {
      idx <- object$malmquist$period_from == p
      sub <- object$malmquist[idx, c("MPI", "TEC", "TGC", "TC")]
      cat("Period", p, "->",
          unique(object$malmquist$period_to[idx]), "\n")
      print(round(colMeans(sub, na.rm = TRUE), 4))
      cat("\n")
    }
  }

  # Overall TGR change
  cat("--- Technology Gap Ratios ---\n\n")
  for (g in object$groups) {
    idx <- object$tgr$group == g
    sub <- object$tgr[idx, ]
    cat("Group:", g, "\n")
    cat("  Mean TGR (from):", format(mean(sub$TGR_from, na.rm = TRUE),
                                     digits = 4), "\n")
    cat("  Mean TGR (to):  ", format(mean(sub$TGR_to, na.rm = TRUE),
                                     digits = 4), "\n")
    cat("  Mean TGC:       ", format(mean(sub$TGC, na.rm = TRUE),
                                     digits = 4), "\n\n")
  }

  invisible(object)
}
