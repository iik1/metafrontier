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

  # Solve DEA LP for each DMU (batch for performance)
  efficiency <- .dea_batch_fast(X_raw, Y, orientation, rts)

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


#' Batch DEA with LP reuse for performance
#'
#' Creates a single LP object and updates constraint coefficients
#' per DMU instead of rebuilding from scratch. Falls back to
#' per-DMU LP creation for small n where overhead is negligible.
#'
#' @keywords internal
#' @noRd
.dea_batch_fast <- function(X, Y, orientation, rts,
                            X_ref = NULL, Y_ref = NULL) {
  # If X_ref/Y_ref not supplied, evaluate X/Y against itself
  if (is.null(X_ref)) X_ref <- X
  if (is.null(Y_ref)) Y_ref <- Y

  n_eval <- nrow(X)
  n_ref  <- nrow(X_ref)
  m <- ncol(X)
  s <- ncol(Y)

  # For small n, per-DMU creation is fine
  if (n_eval <= 50) {
    eff <- numeric(n_eval)
    for (i in seq_len(n_eval)) {
      eff[i] <- .dea_solve_lp(X[i, ], Y[i, ], X_ref, Y_ref, orientation, rts)
    }
    return(eff)
  }

  # Build template LP once (lambdas correspond to reference set)
  n_vars <- n_ref + 1  # n_ref lambdas + theta/phi
  n_cons_base <- m + s
  has_rts <- rts != "crs"

  if (orientation == "input") {
    lp <- lpSolveAPI::make.lp(0, n_vars)
    lpSolveAPI::set.objfn(lp, c(rep(0, n_ref), 1))
    lpSolveAPI::lp.control(lp, sense = "min", verbose = "neutral")

    # Input constraints (m): X_ref[,mm] lambda - x_i[mm]*theta <= 0
    for (mm in seq_len(m)) {
      lpSolveAPI::add.constraint(lp, c(X_ref[, mm], 0), "<=", 0)
    }
    # Output constraints (s): Y_ref[,ss] lambda >= y_i[ss]
    for (ss in seq_len(s)) {
      lpSolveAPI::add.constraint(lp, c(Y_ref[, ss], 0), ">=", 0)
    }
    # RTS
    if (rts == "vrs") {
      lpSolveAPI::add.constraint(lp, c(rep(1, n_ref), 0), "=", 1)
    } else if (rts == "drs") {
      lpSolveAPI::add.constraint(lp, c(rep(1, n_ref), 0), "<=", 1)
    } else if (rts == "irs") {
      lpSolveAPI::add.constraint(lp, c(rep(1, n_ref), 0), ">=", 1)
    }
    lpSolveAPI::set.bounds(lp, lower = c(rep(0, n_ref), -Inf))

    eff <- numeric(n_eval)
    for (i in seq_len(n_eval)) {
      # Update input constraint RHS: the theta coefficient col
      for (mm in seq_len(m)) {
        lpSolveAPI::set.mat(lp, mm, n_vars, -X[i, mm])
        lpSolveAPI::set.rhs(lp, 0, mm)
      }
      # Update output constraint RHS
      for (ss in seq_len(s)) {
        lpSolveAPI::set.rhs(lp, Y[i, ss], m + ss)
      }

      status <- lpSolveAPI::solve.lpExtPtr(lp)
      if (status == 0) {
        eff[i] <- lpSolveAPI::get.objective(lp)
      } else {
        warning("LP infeasible for a DMU.", call. = FALSE)
        eff[i] <- NA_real_
      }
    }

  } else {
    # Output-oriented
    lp <- lpSolveAPI::make.lp(0, n_vars)
    lpSolveAPI::set.objfn(lp, c(rep(0, n_ref), 1))
    lpSolveAPI::lp.control(lp, sense = "max", verbose = "neutral")

    # Input constraints: X_ref[,mm] lambda <= x_i[mm]
    for (mm in seq_len(m)) {
      lpSolveAPI::add.constraint(lp, c(X_ref[, mm], 0), "<=", 0)
    }
    # Output constraints: Y_ref[,ss] lambda - y_i[ss]*phi >= 0
    for (ss in seq_len(s)) {
      lpSolveAPI::add.constraint(lp, c(Y_ref[, ss], 0), ">=", 0)
    }
    if (rts == "vrs") {
      lpSolveAPI::add.constraint(lp, c(rep(1, n_ref), 0), "=", 1)
    } else if (rts == "drs") {
      lpSolveAPI::add.constraint(lp, c(rep(1, n_ref), 0), "<=", 1)
    } else if (rts == "irs") {
      lpSolveAPI::add.constraint(lp, c(rep(1, n_ref), 0), ">=", 1)
    }
    lpSolveAPI::set.bounds(lp, lower = c(rep(0, n_ref), 1))

    eff <- numeric(n_eval)
    for (i in seq_len(n_eval)) {
      # Update input RHS
      for (mm in seq_len(m)) {
        lpSolveAPI::set.rhs(lp, X[i, mm], mm)
      }
      # Update output constraint phi coefficient
      for (ss in seq_len(s)) {
        lpSolveAPI::set.mat(lp, m + ss, n_vars, -Y[i, ss])
        lpSolveAPI::set.rhs(lp, 0, m + ss)
      }

      status <- lpSolveAPI::solve.lpExtPtr(lp)
      if (status == 0) {
        eff[i] <- 1 / lpSolveAPI::get.objective(lp)
      } else {
        warning("LP infeasible for a DMU.", call. = FALSE)
        eff[i] <- NA_real_
      }
    }
  }

  eff
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
  te_meta <- .dea_batch_fast(X_raw, Y, orientation, rts)

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


#' Solve a single DDF LP
#'
#' Computes the directional distance function:
#' max beta s.t. Y'lambda >= y_i + beta*g_y,
#'               X'lambda <= x_i - beta*g_x,
#'               RTS constraints, lambda >= 0, beta >= 0.
#'
#' @param x_i input vector for DMU i.
#' @param y_i output vector for DMU i.
#' @param X input matrix (reference set).
#' @param Y output matrix (reference set).
#' @param g_x input direction vector.
#' @param g_y output direction vector.
#' @param rts returns to scale assumption.
#'
#' @return The DDF value beta (inefficiency measure).
#' @keywords internal
#' @noRd
.dea_solve_ddf <- function(x_i, y_i, X, Y, g_x, g_y, rts) {
  n <- nrow(X)
  m <- ncol(X)
  s <- ncol(Y)

  # Variables: lambda_1..lambda_n, beta
  n_vars <- n + 1
  lp <- lpSolveAPI::make.lp(0, n_vars)
  lpSolveAPI::set.objfn(lp, c(rep(0, n), 1))  # max beta
  lpSolveAPI::lp.control(lp, sense = "max", verbose = "neutral")

  # Input constraints: X'lambda + beta*g_x <= x_i
  for (mm in seq_len(m)) {
    lpSolveAPI::add.constraint(lp, c(X[, mm], g_x[mm]), "<=", x_i[mm])
  }

  # Output constraints: Y'lambda - beta*g_y >= y_i
  for (ss in seq_len(s)) {
    lpSolveAPI::add.constraint(lp, c(Y[, ss], -g_y[ss]), ">=", y_i[ss])
  }

  # RTS constraints
  if (rts == "vrs") {
    lpSolveAPI::add.constraint(lp, c(rep(1, n), 0), "=", 1)
  } else if (rts == "drs") {
    lpSolveAPI::add.constraint(lp, c(rep(1, n), 0), "<=", 1)
  } else if (rts == "irs") {
    lpSolveAPI::add.constraint(lp, c(rep(1, n), 0), ">=", 1)
  }

  # Bounds
  lpSolveAPI::set.bounds(lp, lower = c(rep(0, n), 0))

  status <- lpSolveAPI::solve.lpExtPtr(lp)

  if (status == 0) {
    lpSolveAPI::get.objective(lp)
  } else {
    warning("DDF LP infeasible for a DMU.", call. = FALSE)
    NA_real_
  }
}


#' Batch DDF solver
#' @keywords internal
#' @noRd
.ddf_batch <- function(X, Y, g_x_mat, g_y_mat, rts) {
  n <- nrow(X)
  eff <- numeric(n)
  for (i in seq_len(n)) {
    eff[i] <- .dea_solve_ddf(
      X[i, ], Y[i, ], X, Y,
      g_x = g_x_mat[i, ], g_y = g_y_mat[i, ],
      rts = rts
    )
  }
  eff
}


#' DDF-based group efficiency
#' @keywords internal
#' @noRd
.fit_ddf_group <- function(formula, data, rts, direction, ...) {

  mf <- model.frame(formula, data = data)
  y_raw <- model.response(mf)
  X_raw <- model.matrix(formula, data = data, rhs = 1)

  if (colnames(X_raw)[1] == "(Intercept)") {
    X_raw <- X_raw[, -1, drop = FALSE]
  }

  n <- nrow(X_raw)
  m <- ncol(X_raw)

  if (is.matrix(y_raw)) {
    Y <- y_raw
  } else {
    Y <- matrix(y_raw, ncol = 1)
  }
  s <- ncol(Y)

  # Direction vectors
  if (direction == "proportional") {
    g_x_mat <- X_raw
    g_y_mat <- Y
  } else if (direction == "output") {
    g_x_mat <- matrix(0, n, m)
    g_y_mat <- Y
  } else if (direction == "input") {
    g_x_mat <- X_raw
    g_y_mat <- matrix(0, n, s)
  } else {
    stop("Unknown direction: ", direction, call. = FALSE)
  }

  efficiency <- .ddf_batch(X_raw, Y, g_x_mat, g_y_mat, rts)

  list(
    efficiency = efficiency,
    nobs = n,
    rts = rts,
    direction = direction,
    X = X_raw,
    Y = Y
  )
}


#' DDF-based metafrontier estimation
#' @keywords internal
#' @noRd
.estimate_ddf_metafrontier <- function(formula, data, group_vec,
                                       group_levels, group_models,
                                       rts, direction) {

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
  m <- ncol(X_raw)
  s <- ncol(Y)

  # Direction vectors
  if (direction == "proportional") {
    g_x_mat <- X_raw
    g_y_mat <- Y
  } else if (direction == "output") {
    g_x_mat <- matrix(0, n, m)
    g_y_mat <- Y
  } else {
    g_x_mat <- X_raw
    g_y_mat <- matrix(0, n, s)
  }

  # Group-level DDF
  beta_group <- numeric(n)
  for (g in group_levels) {
    idx <- which(group_vec == g)
    beta_group[idx] <- group_models[[g]]$efficiency
  }

  # Pooled DDF (metafrontier)
  beta_meta <- .ddf_batch(X_raw, Y, g_x_mat, g_y_mat, rts)

  # DDF TGR: additive decomposition
  # beta_meta = beta_group + TGR_DDF
  # TGR = beta_meta - beta_group (additive gap)
  # For compatibility, also compute ratio-based TE
  te_group <- 1 / (1 + beta_group)
  te_meta <- 1 / (1 + beta_meta)
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
    beta_group = beta_group,
    beta_meta = beta_meta,
    ddf_tgr = beta_meta - beta_group
  )
}
