#' Internal: DEA estimation for a single group
#'
#' Solves the DEA linear programme for a single technology group.
#'
#' @keywords internal
#' @noRd
.fit_dea_group <- function(formula, data, orientation, rts, ...) {

  mf <- model.frame(formula, data = data)
  y_raw <- model.response(mf)
  X_raw <- model.matrix(formula, data = data, rhs = 1)

  # Remove intercept for DEA (inputs are raw, not in a regression)
  if (colnames(X_raw)[1] == "(Intercept)") {
    X_raw <- X_raw[, -1, drop = FALSE]
  }

  # For DEA, we expect raw (non-logged) inputs and outputs
  # If the user passed log values, that is their choice
  n <- nrow(X_raw)
  m <- ncol(X_raw)  # number of inputs

  # Handle multiple outputs
  if (is.matrix(y_raw)) {
    Y <- y_raw
  } else {
    Y <- matrix(y_raw, ncol = 1)
  }
  s <- ncol(Y)  # number of outputs

  # Solve DEA LP for each DMU
  efficiency <- numeric(n)

  for (i in seq_len(n)) {
    efficiency[i] <- .dea_solve_lp(
      x_i = X_raw[i, ],
      y_i = Y[i, ],
      X = X_raw,
      Y = Y,
      orientation = orientation,
      rts = rts
    )
  }

  list(
    efficiency = efficiency,
    nobs = n,
    orientation = orientation,
    rts = rts,
    X = X_raw,
    Y = Y
  )
}


#' Solve a single DEA LP
#' @keywords internal
#' @noRd
.dea_solve_lp <- function(x_i, y_i, X, Y, orientation, rts) {
  n <- nrow(X)
  m <- ncol(X)
  s <- ncol(Y)

  if (orientation == "input") {
    # Input-oriented: min theta
    # s.t. X' lambda <= theta * x_i
    #      Y' lambda >= y_i
    #      RTS constraints on lambda
    #      lambda >= 0

    lp <- lpSolveAPI::make.lp(0, n + 1)  # n lambdas + theta
    lpSolveAPI::set.objfn(lp, c(rep(0, n), 1))  # min theta

    # Input constraints: sum_j X_mj lambda_j <= theta * x_im
    for (mm in seq_len(m)) {
      lpSolveAPI::add.constraint(lp, c(X[, mm], -x_i[mm]), "<=", 0)
    }

    # Output constraints: sum_j Y_sj lambda_j >= y_is
    for (ss in seq_len(s)) {
      lpSolveAPI::add.constraint(lp, c(Y[, ss], 0), ">=", y_i[ss])
    }

    # RTS constraints
    if (rts == "vrs") {
      lpSolveAPI::add.constraint(lp, c(rep(1, n), 0), "=", 1)
    } else if (rts == "drs") {
      lpSolveAPI::add.constraint(lp, c(rep(1, n), 0), "<=", 1)
    } else if (rts == "irs") {
      lpSolveAPI::add.constraint(lp, c(rep(1, n), 0), ">=", 1)
    }

    # Bounds: lambda >= 0, theta free (but typically in [0,1])
    lpSolveAPI::set.bounds(lp, lower = c(rep(0, n), -Inf))

    lpSolveAPI::lp.control(lp, sense = "min", verbose = "neutral")
    status <- lpSolveAPI::solve.lpExtPtr(lp)

    if (status == 0) {
      return(lpSolveAPI::get.objective(lp))
    } else {
      warning("LP infeasible for a DMU.", call. = FALSE)
      return(NA_real_)
    }

  } else {
    # Output-oriented: max phi
    # s.t. X' lambda <= x_i
    #      Y' lambda >= phi * y_i
    #      RTS constraints
    #      lambda >= 0

    lp <- lpSolveAPI::make.lp(0, n + 1)
    lpSolveAPI::set.objfn(lp, c(rep(0, n), 1))

    # Input constraints
    for (mm in seq_len(m)) {
      lpSolveAPI::add.constraint(lp, c(X[, mm], 0), "<=", x_i[mm])
    }

    # Output constraints: sum_j Y_sj lambda_j >= phi * y_is
    for (ss in seq_len(s)) {
      lpSolveAPI::add.constraint(lp, c(Y[, ss], -y_i[ss]), ">=", 0)
    }

    # RTS
    if (rts == "vrs") {
      lpSolveAPI::add.constraint(lp, c(rep(1, n), 0), "=", 1)
    } else if (rts == "drs") {
      lpSolveAPI::add.constraint(lp, c(rep(1, n), 0), "<=", 1)
    } else if (rts == "irs") {
      lpSolveAPI::add.constraint(lp, c(rep(1, n), 0), ">=", 1)
    }

    lpSolveAPI::set.bounds(lp, lower = c(rep(0, n), 1))

    lpSolveAPI::lp.control(lp, sense = "max", verbose = "neutral")
    status <- lpSolveAPI::solve.lpExtPtr(lp)

    if (status == 0) {
      return(1 / lpSolveAPI::get.objective(lp))  # Farrell efficiency
    } else {
      warning("LP infeasible for a DMU.", call. = FALSE)
      return(NA_real_)
    }
  }
}


#' Internal: Estimate DEA-based metafrontier
#'
#' Computes group DEA and pooled (meta) DEA, then derives TGR.
#'
#' @keywords internal
#' @noRd
.estimate_dea_metafrontier <- function(formula, data, group_vec,
                                       group_levels, group_models,
                                       orientation, rts) {

  mf <- model.frame(formula, data = data)
  y_raw <- model.response(mf)
  X_raw <- model.matrix(formula, data = data, rhs = 1)
  if (colnames(X_raw)[1] == "(Intercept)") {
    X_raw <- X_raw[, -1, drop = FALSE]
  }

  if (!is.matrix(y_raw)) {
    Y <- matrix(y_raw, ncol = 1)
  } else {
    Y <- y_raw
  }

  n <- nrow(X_raw)

  # Group-level efficiency
  te_group <- numeric(n)
  for (g in group_levels) {
    idx <- which(group_vec == g)
    te_group[idx] <- group_models[[g]]$efficiency
  }

  # Pooled DEA (metafrontier): use ALL observations
  te_meta <- numeric(n)
  for (i in seq_len(n)) {
    te_meta[i] <- .dea_solve_lp(
      x_i = X_raw[i, ],
      y_i = Y[i, ],
      X = X_raw,
      Y = Y,
      orientation = orientation,
      rts = rts
    )
  }

  # TGR = metafrontier efficiency / group efficiency
  tgr <- te_meta / te_group

  list(
    meta_coef = NULL,  # DEA is nonparametric
    meta_vcov = NULL,
    group_coef = NULL,
    tgr = tgr,
    te_group = te_group,
    te_meta = te_meta,
    group_frontier = NULL,
    meta_frontier = NULL,
    logLik_groups = NULL,
    meta_logLik = NULL,
    meta_convergence = 0L
  )
}
