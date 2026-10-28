#' Fare-Primont and Lowe TFP Indexes with a Metafrontier Decomposition
#'
#' Computes multiplicatively complete and transitive total factor
#' productivity (TFP) indexes, the Fare-Primont and the Lowe index, and
#' decomposes TFP efficiency following O'Donnell (2012), both against
#' each group's technology and against the metatechnology (all groups
#' pooled).
#'
#' @param formula a formula with the outputs on the left and the inputs
#'   on the right, e.g. \code{cbind(y1, y2) ~ x1 + x2}. Quantities must
#'   be non-negative (levels, not logs).
#' @param data a data frame.
#' @param group,time,id column names in \code{data} (or vectors of
#'   length \code{nrow(data)}) identifying the technology group, the
#'   period and the unit.
#' @param index \code{"fare-primont"} (default) or \code{"lowe"}.
#' @param orientation \code{"output"} (default) or \code{"input"}:
#'   which efficiency components the TFP efficiency is decomposed into.
#' @param rts returns to scale of the technologies: \code{"vrs"}
#'   (default) or \code{"crs"}.
#' @param tech \code{"contemporaneous"} (default; the technology of
#'   period t is formed from the observations of period t) or
#'   \code{"sequential"} (from all observations up to period t, so no
#'   technical regress).
#' @param output_prices,input_prices column names of the output and
#'   input prices, in the order of the outputs and inputs in
#'   \code{formula}. Required for the Lowe index, which aggregates with
#'   their sample means; ignored for the Fare-Primont index.
#'
#' @details
#' Output and input aggregates are linear, \eqn{Q(y) = p'y} and
#' \eqn{X(x) = w'x}, with the same weights for every observation, so TFP
#' = Q/X is a level and the index between any two observations, the
#' ratio of their TFP levels, is transitive. The Fare-Primont weights
#' are the normalised shadow prices of the output and input distance
#' functions at the sample-mean observation, evaluated against the
#' pooled data (O'Donnell, 2012); the Lowe weights are the sample means
#' of the observed prices.
#'
#' TFP efficiency is TFPE = TFP / TFP*, where TFP* is the maximum TFP
#' attainable in the reference technology. With output orientation it is
#' decomposed as TFPE = OTE x OME x ROSE = OTE x OSE x RME (technical,
#' mix, residual scale, scale and residual mix efficiency), with scale-mix
#' efficiency OSME = OME x ROSE; with input orientation as TFPE = ITE x
#' IME x RISE = ITE x ISE x RME, with ISME = IME x RISE.
#'
#' Each component is computed against the technology of the unit's group
#' (suffix \code{_group}) and against the metatechnology formed from all
#' groups (suffix \code{_meta}). The TFP metatechnology ratio is MTR =
#' TFP*_group / TFP*_meta, in (0, 1], so that TFPE_meta = TFPE_group x
#' MTR.
#'
#' Programmes that are infeasible give \code{NA} and are counted in a
#' single warning.
#'
#' @return An object of class \code{"tfp_meta"}: a list with
#'   \describe{
#'     \item{levels}{one row per observation: id, group, time, the
#'       aggregates \code{AO} and \code{AI}, \code{TFP}, and for each
#'       reference (\code{_group}, \code{_meta}) \code{TFP_star},
#'       \code{TFPE} and the efficiency components, plus \code{MTR}}
#'     \item{changes}{one row per unit and pair of consecutive periods
#'       of the sample in which it is observed in both: the ratio
#'       (period_to over period_from)
#'       of every level column, prefixed with \code{d}}
#'     \item{weights}{the output and input aggregation weights}
#'     \item{settings, groups, periods, n_infeasible, call}{settings and
#'       metadata}
#'   }
#'
#' @references
#' O'Donnell, C.J. (2012). Nonparametric estimates of the components of
#' productivity and profitability change in U.S. agriculture.
#' \emph{American Journal of Agricultural Economics}, 94(4), 873--890.
#' \doi{10.1093/ajae/aas023}
#'
#' O'Donnell, C.J. (2018). \emph{Productivity and Efficiency Analysis}.
#' Springer. \doi{10.1007/978-981-13-2984-5}
#'
#' @examples
#' set.seed(1)
#' d <- expand.grid(firm = 1:6, year = 1:3)
#' d$region <- ifelse(d$firm <= 3, "A", "B")
#' d$x1 <- runif(nrow(d), 1, 10)
#' d$x2 <- runif(nrow(d), 1, 10)
#' d$y <- d$x1^0.5 * d$x2^0.3 * runif(nrow(d), 0.6, 1)
#' tfp <- tfp_index(y ~ x1 + x2, data = d, group = "region",
#'                  time = "year", id = "firm")
#' tfp
#' head(tfp$levels[, c("firm", "year", "TFP", "TFPE_meta", "MTR")])
#'
#' @seealso \code{\link{malmquist_meta}}
#' @export
tfp_index <- function(formula, data, group, time, id,
                      index = c("fare-primont", "lowe"),
                      orientation = c("output", "input"),
                      rts = c("vrs", "crs"),
                      tech = c("contemporaneous", "sequential"),
                      output_prices = NULL, input_prices = NULL) {

  cl <- match.call()
  index <- match.arg(index)
  orientation <- match.arg(orientation)
  rts <- match.arg(rts)
  tech <- match.arg(tech)

  column <- function(v, arg) {
    if (is.character(v) && length(v) == 1L) {
      if (!v %in% names(data)) {
        stop("Column '", v, "' not found in data.", call. = FALSE)
      }
      return(data[[v]])
    }
    if (length(v) != nrow(data)) {
      stop("'", arg, "' must be a column name or a vector of length ",
           "nrow(data).", call. = FALSE)
    }
    v
  }
  group_vec <- column(group, "group")
  time_vec <- column(time, "time")
  id_vec <- column(id, "id")

  price_cols <- c(output_prices, input_prices)
  if (index == "lowe") {
    if (is.null(output_prices) || is.null(input_prices)) {
      stop("The Lowe index requires 'output_prices' and 'input_prices'.",
           call. = FALSE)
    }
    missing_cols <- setdiff(price_cols, names(data))
    if (length(missing_cols)) {
      stop("Price column(s) not found in data: ",
           paste(missing_cols, collapse = ", "), ".", call. = FALSE)
    }
  }

  # Drop rows with a missing quantity, price, group, period or id once,
  # here, so every index lines up with its row
  f <- Formula::Formula(formula)
  keep <- complete.cases(model.frame(f, data = data, na.action = na.pass)) &
    !is.na(group_vec) & !is.na(time_vec) & !is.na(id_vec)
  if (index == "lowe") keep <- keep & complete.cases(data[price_cols])
  data <- data[keep, , drop = FALSE]
  group_vec <- as.character(group_vec[keep])
  time_vec <- time_vec[keep]
  id_vec <- id_vec[keep]

  mats <- .dea_matrices_mf(f, data)
  X <- mats$X
  Y <- mats$Y
  if (any(X < 0) || any(Y < 0)) {
    stop("Quantities must be non-negative (levels, not logs).",
         call. = FALSE)
  }
  if (anyDuplicated(data.frame(id_vec, time_vec))) {
    stop("Duplicated (id, time) combinations.", call. = FALSE)
  }

  # Aggregation weights, common to all observations
  if (index == "fare-primont") {
    weights <- .fp_weights(X, Y, rts)
  } else {
    if (length(output_prices) != ncol(Y) ||
        length(input_prices) != ncol(X)) {
      stop("Supply one price column per output (", ncol(Y), ") and per ",
           "input (", ncol(X), ").", call. = FALSE)
    }
    weights <- list(output = colMeans(data[output_prices]),
                    input = colMeans(data[input_prices]))
    if (any(data[price_cols] <= 0)) {
      stop("Prices must be positive.", call. = FALSE)
    }
  }
  AO <- as.numeric(Y %*% weights$output)
  AI <- as.numeric(X %*% weights$input)

  periods <- sort(unique(time_vec))
  groups <- sort(unique(group_vec))
  comp_names <- if (orientation == "output") {
    c("OTE", "OSE", "OME", "ROSE", "OSME", "RME")
  } else {
    c("ITE", "ISE", "IME", "RISE", "ISME", "RME")
  }
  decomp <- list(
    group = matrix(NA_real_, nrow(X), 2L + length(comp_names)),
    meta = matrix(NA_real_, nrow(X), 2L + length(comp_names))
  )

  for (i in seq_len(nrow(X))) {
    in_period <- if (tech == "contemporaneous") {
      time_vec == time_vec[i]
    } else {
      time_vec <= time_vec[i]
    }
    refs <- list(group = in_period & group_vec == group_vec[i],
                 meta = in_period)
    for (r in names(refs)) {
      decomp[[r]][i, ] <- .tfp_decompose(
        X[i, ], Y[i, ], AO[i], AI[i],
        X[refs[[r]], , drop = FALSE], Y[refs[[r]], , drop = FALSE],
        weights, orientation, rts
      )
    }
  }

  levels <- data.frame(id_vec, group = group_vec, time_vec,
                       AO = AO, AI = AI, TFP = AO / AI,
                       stringsAsFactors = FALSE)
  names(levels)[c(1, 3)] <- c(
    if (is.character(id)) id else "id",
    if (is.character(time)) time else "time"
  )
  for (r in c("group", "meta")) {
    colnames(decomp[[r]]) <- paste0(c("TFP_star", "TFPE", comp_names),
                                    "_", r)
    levels <- cbind(levels, decomp[[r]])
  }
  levels$MTR <- levels$TFP_star_group / levels$TFP_star_meta

  n_infeasible <- sum(is.na(levels$TFP_star_group)) +
    sum(is.na(levels$TFP_star_meta))
  if (n_infeasible > 0L) {
    warning(n_infeasible, " TFP-maximisation programme(s) were infeasible; ",
            "the affected values are NA.", call. = FALSE)
  }

  out <- list(
    levels = levels,
    changes = .tfp_changes(levels, id_vec, time_vec, periods),
    weights = weights,
    settings = list(index = index, orientation = orientation, rts = rts,
                    tech = tech),
    groups = groups,
    periods = periods,
    n_infeasible = n_infeasible,
    call = cl
  )
  class(out) <- "tfp_meta"
  out
}


# Solve a small LP given as dense constraint rows; NA when not optimal
.lp_dense <- function(obj, A, dir, rhs, sense, free = integer(0)) {
  lp <- lpSolveAPI::make.lp(0, length(obj))
  lpSolveAPI::lp.control(lp, sense = sense, verbose = "neutral")
  lpSolveAPI::set.objfn(lp, obj)
  for (r in seq_len(nrow(A))) {
    lpSolveAPI::add.constraint(lp, A[r, ], dir[r], rhs[r])
  }
  if (length(free)) {
    lpSolveAPI::set.bounds(lp, lower = rep(-Inf, length(free)),
                           columns = free)
  }
  if (lpSolveAPI::solve.lpExtPtr(lp) != 0L) {
    return(list(objective = NA_real_, solution = rep(NA_real_, length(obj))))
  }
  list(objective = lpSolveAPI::get.objective(lp),
       solution = lpSolveAPI::get.variables(lp))
}


#' Fare-Primont aggregation weights
#'
#' Normalised shadow prices of the output and input distance functions
#' at the sample-mean observation against the pooled data (O'Donnell,
#' 2012), as in productivity::fareprim().
#'
#' @keywords internal
#' @noRd
.fp_weights <- function(X, Y, rts) {
  n <- nrow(X)
  m <- ncol(X)
  s <- ncol(Y)
  x0 <- colMeans(X)
  y0 <- colMeans(Y)
  vrs <- rts == "vrs"
  w0 <- if (vrs) 1   # free intercept column under VRS

  # Output distance: min v'x0 (+ w) s.t. u'y0 = 1, v'x_j - u'y_j (+ w) >= 0
  A <- rbind(cbind(-Y, X, if (vrs) 1), c(y0, rep(0, m), if (vrs) 0))
  o <- .lp_dense(c(rep(0, s), x0, w0), A, c(rep(">=", n), "="),
                 c(rep(0, n), 1), "min", if (vrs) s + m + 1L)
  u <- o$solution[seq_len(s)]
  v <- o$solution[s + seq_len(m)]
  p <- u / (sum(v * x0) + if (vrs) o$solution[s + m + 1L] else 0)

  # Input distance: max u'y0 (+ w) s.t. v'x0 = 1, u'y_j - v'x_j (+ w) <= 0
  A <- rbind(cbind(Y, -X, if (vrs) 1), c(rep(0, s), x0, if (vrs) 0))
  i <- .lp_dense(c(y0, rep(0, m), w0), A, c(rep("<=", n), "="),
                 c(rep(0, n), 1), "max", if (vrs) s + m + 1L)
  u <- i$solution[seq_len(s)]
  v <- i$solution[s + seq_len(m)]
  w <- v / (sum(u * y0) + if (vrs) i$solution[s + m + 1L] else 0)

  list(output = setNames(p, colnames(Y)), input = setNames(w, colnames(X)))
}


#' TFP* and the O'Donnell (2012) decomposition of TFP efficiency
#'
#' Returns c(TFP_star, TFPE, components) for one observation against
#' the reference set (Xr, Yr).
#'
#' @keywords internal
#' @noRd
.tfp_decompose <- function(x, y, AO, AI, Xr, Yr, weights, orientation,
                           rts) {
  n <- nrow(Xr)
  m <- ncol(Xr)
  s <- ncol(Yr)
  p <- weights$output
  w <- weights$input
  vrs_row <- function(lead, lag = 0) c(rep(0, lead), rep(1, n), rep(0, lag))

  # TFP*: max p'y s.t. w'x = 1, x >= Xr'l, y <= Yr'l (VRS: sum l = t)
  A <- rbind(cbind(matrix(0, m, s), -diag(m), t(Xr), 0),
             cbind(-diag(s), matrix(0, s, m), t(Yr), 0),
             c(rep(0, s), w, rep(0, n + 1)))
  dir <- c(rep("<=", m), rep(">=", s), "=")
  rhs <- c(rep(0, m + s), 1)
  if (rts == "vrs") {
    A <- rbind(A, c(rep(0, s + m), rep(1, n), -1))
    dir <- c(dir, "=")
    rhs <- c(rhs, 0)
  }
  tfp_star <- .lp_dense(c(p, rep(0, m + n + 1)), A, dir, rhs,
                        "max")$objective
  tfpe <- (AO / AI) / tfp_star

  if (orientation == "output") {
    # Radial output efficiency, CRS then VRS: max phi s.t. Xr'l <= x,
    # phi y <= Yr'l
    A <- rbind(cbind(0, t(Xr)), cbind(-y, t(Yr)))
    dir <- c(rep("<=", m), rep(">=", s))
    rhs <- c(x, rep(0, s))
    ote_crs <- 1 / .lp_dense(c(1, rep(0, n)), A, dir, rhs, "max",
                             1L)$objective
    ote <- if (rts == "vrs") {
      1 / .lp_dense(c(1, rep(0, n)), rbind(A, vrs_row(1)), c(dir, "="),
                    c(rhs, 1), "max", 1L)$objective
    } else ote_crs
    ose <- ote_crs / ote
    # Mix: max p'y*/p'y s.t. Xr'l <= x, y* <= Yr'l
    A <- rbind(cbind(matrix(0, m, s), t(Xr)), cbind(-diag(s), t(Yr)))
    if (rts == "vrs") {
      A <- rbind(A, vrs_row(s))
      dir <- c(dir, "=")
      rhs <- c(rhs, 1)
    }
    mix <- .lp_dense(c(p / sum(p * y), rep(0, n)), A, dir, rhs,
                     "max")$objective
    ome <- 1 / (mix * ote)
    rose <- (AO / (ote * ome) / AI) / tfp_star
    c(tfp_star, tfpe, ote, ose, ome, rose, ome * rose, tfpe / ote / ose)
  } else {
    # Radial input efficiency, CRS then VRS: min theta s.t.
    # Xr'l <= theta x, Yr'l >= y
    A <- rbind(cbind(-x, t(Xr)), cbind(0, t(Yr)))
    dir <- c(rep("<=", m), rep(">=", s))
    rhs <- c(rep(0, m), y)
    ite_crs <- .lp_dense(c(1, rep(0, n)), A, dir, rhs, "min", 1L)$objective
    ite <- if (rts == "vrs") {
      .lp_dense(c(1, rep(0, n)), rbind(A, vrs_row(1)), c(dir, "="),
                c(rhs, 1), "min", 1L)$objective
    } else ite_crs
    ise <- ite_crs / ite
    # Mix: min w'x*/w'x s.t. Xr'l <= x*, Yr'l >= y
    A <- rbind(cbind(-diag(m), t(Xr)), cbind(matrix(0, s, m), t(Yr)))
    if (rts == "vrs") {
      A <- rbind(A, vrs_row(m))
      dir <- c(dir, "=")
      rhs <- c(rhs, 1)
    }
    mix <- .lp_dense(c(w / sum(w * x), rep(0, n)), A, dir, rhs,
                     "min")$objective
    ime <- mix / ite
    rise <- (AO / (AI * ime * ite)) / tfp_star
    c(tfp_star, tfpe, ite, ise, ime, rise, ime * rise, tfpe / ite / ise)
  }
}


# Ratios of every level column between consecutive observed periods of
# each unit
.tfp_changes <- function(levels, id_vec, time_vec, periods) {
  value_cols <- setdiff(names(levels), names(levels)[1:3])
  rows <- list()
  for (k in seq_len(length(periods) - 1L)) {
    from <- which(time_vec == periods[k])
    to <- which(time_vec == periods[k + 1L])
    to <- to[match(id_vec[from], id_vec[to])]
    ok <- !is.na(to)
    from <- from[ok]
    to <- to[ok]
    if (!length(from)) next
    ratio <- levels[to, value_cols, drop = FALSE] /
      levels[from, value_cols, drop = FALSE]
    names(ratio) <- paste0("d", value_cols)
    rows[[k]] <- data.frame(levels[from, 1:2, drop = FALSE],
                            period_from = periods[k],
                            period_to = periods[k + 1L], ratio,
                            row.names = NULL)
  }
  do.call(rbind, rows)
}


#' @export
print.tfp_meta <- function(x, digits = 4, ...) {
  idx_name <- if (x$settings$index == "lowe") "Lowe" else "Fare-Primont"
  cat("\n", idx_name, " TFP index with metafrontier decomposition\n",
      sep = "")
  cat("Orientation: ", x$settings$orientation, "  RTS: ", x$settings$rts,
      "  Technology: ", x$settings$tech, "\n", sep = "")
  cat("Observations:", nrow(x$levels), " Groups:",
      paste(x$groups, collapse = ", "), " Periods:",
      paste(x$periods, collapse = ", "), "\n\n")
  means <- aggregate(x$levels[c("TFP", "TFPE_group", "MTR", "TFPE_meta")],
                     by = list(group = x$levels$group), FUN = mean,
                     na.rm = TRUE)
  cat("Mean levels by group:\n")
  print(means, digits = digits, row.names = FALSE)
  invisible(x)
}


#' @export
summary.tfp_meta <- function(object, ...) {
  lv <- object$levels
  by <- list(group = lv$group, time = lv[[3]])
  cols <- c("TFP", "TFP_star_group", "TFPE_group", "TFP_star_meta",
            "TFPE_meta", "MTR")
  out <- list(
    settings = object$settings,
    levels = aggregate(lv[cols], by = by, FUN = mean, na.rm = TRUE),
    changes = if (!is.null(object$changes)) {
      ch <- object$changes
      aggregate(ch[paste0("d", cols)],
                by = list(group = ch$group, period_from = ch$period_from,
                          period_to = ch$period_to),
                FUN = mean, na.rm = TRUE)
    },
    n_infeasible = object$n_infeasible
  )
  class(out) <- "summary.tfp_meta"
  out
}


#' @export
print.summary.tfp_meta <- function(x, digits = 4, ...) {
  cat("\nTFP index summary (", x$settings$index, ", ",
      x$settings$orientation, " orientation, ", x$settings$rts, ", ",
      x$settings$tech, " technology)\n\n", sep = "")
  cat("Mean levels by group and period:\n")
  print(x$levels, digits = digits, row.names = FALSE)
  if (!is.null(x$changes)) {
    cat("\nMean changes by group and period pair:\n")
    print(x$changes, digits = digits, row.names = FALSE)
  }
  if (x$n_infeasible > 0L) {
    cat("\nInfeasible programmes:", x$n_infeasible, "\n")
  }
  invisible(x)
}
