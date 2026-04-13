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
#' @param method character. \code{"dea"} (default) for DEA-based
#'   distance functions or \code{"sfa"} for SFA-based parametric
#'   distance functions.
#' @param dist character. Distribution of the inefficiency term
#'   when \code{method = "sfa"}: \code{"hnormal"} (default),
#'   \code{"tnormal"}, or \code{"exponential"}.
#' @param orientation character. \code{"output"} (default) or
#'   \code{"input"}.
#' @param rts character. Returns to scale assumption:
#'   \code{"crs"} (default), \code{"vrs"}, \code{"drs"}, or
#'   \code{"irs"}.
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
#'       \code{TC} (metafrontier technical change)}
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
#' \strong{Balanced panel assumption:} Firms are matched across
#' periods by position within each group. The data should contain
#' a balanced panel (the same firms observed in every period) with
#' consistent ordering. If group sizes differ across periods,
#' only the first \code{min(n_s, n_t)} firms per group are paired
#' and unmatched observations are silently dropped.
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
                           method = c("dea", "sfa"),
                           dist = c("hnormal", "tnormal", "exponential"),
                           orientation = c("output", "input"),
                           rts = c("crs", "vrs", "drs", "irs"),
                           control = list(),
                           ...) {

  cl <- match.call()
  method <- match.arg(method)
  dist <- match.arg(dist)
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

  # Parse formula
  f <- Formula::Formula(formula)

  # --- SFA Malmquist path ---
  if (method == "sfa") {
    result <- .malmquist_sfa(f, data, group_vec, group_levels,
                             time_vec, time_levels, dist, control)
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
      d_grp_ss[idx_gs] <- .dea_batch_fast(
        X_gs, Y_gs, orientation, rts)

      # D_j^t(x_t, y_t) — same-period group efficiency
      d_grp_tt[idx_gt] <- .dea_batch_fast(
        X_gt, Y_gt, orientation, rts)

      # D_j^s(x_t, y_t) — cross-period: period-t obs against period-s group tech
      d_grp_st[idx_gt] <- suppressWarnings(.dea_batch_fast(
        X_gt, Y_gt, orientation, rts,
        X_ref = X_gs, Y_ref = Y_gs))

      # D_j^t(x_s, y_s) — cross-period: period-s obs against period-t group tech
      d_grp_ts[idx_gs] <- suppressWarnings(.dea_batch_fast(
        X_gs, Y_gs, orientation, rts,
        X_ref = X_gt, Y_ref = Y_gt))
    }

    # Metafrontier (pooled) DEA scores
    # Suppress LP infeasibility warnings from cross-period evaluations
    # (expected when reference technology cannot envelop all eval DMUs)
    d_meta_ss <- suppressWarnings(
      .dea_batch(X_s, Y_s, X_s, Y_s, orientation, rts))
    d_meta_tt <- suppressWarnings(
      .dea_batch(X_t, Y_t, X_t, Y_t, orientation, rts))
    d_meta_st <- suppressWarnings(
      .dea_batch(X_t, Y_t, X_s, Y_s, orientation, rts))
    d_meta_ts <- suppressWarnings(
      .dea_batch(X_s, Y_s, X_t, Y_t, orientation, rts))

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

  if (result_idx == 0L) {
    stop("No matched observations found across time periods.",
         call. = FALSE)
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
    orientation = object$orientation,
    rts         = object$rts,
    groups      = object$groups,
    periods     = object$periods,
    n_obs       = nrow(object$malmquist),
    overall     = colMeans(object$malmquist[, decomp_cols, drop = FALSE],
                           na.rm = TRUE),
    by_group    = by_group,
    by_period   = by_period,
    tgr_summary = tgr_summary
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
  cat("\nOrientation:", x$orientation, "\n")
  cat("RTS:        ", x$rts, "\n")
  cat("Groups:     ", paste(x$groups, collapse = ", "), "\n")
  cat("Periods:    ", paste(x$periods, collapse = " -> "), "\n")
  cat("Observations:", x$n_obs, "\n\n")

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
#' decomposition using SFA distance functions.
#'
#' @keywords internal
#' @noRd
.malmquist_sfa <- function(formula, data, group_vec, group_levels,
                           time_vec, time_levels, dist, control) {

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
        .fit_sfa_group(formula, data_gt, dist, control),
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
      n_pairs <- min(length(idx_s), length(idx_t))
      if (n_pairs == 0) next

      for (i in seq_len(n_pairs)) {
        is <- idx_s[i]
        it <- idx_t[i]

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

        if (is.na(te_s) || is.na(te_t) || te_s <= 0 || te_t <= 0) next

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
          id = i,
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
                            "TGC")]
  )
}
