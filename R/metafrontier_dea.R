#' Internal: DEA estimation for a single group
#'
#' Solves the DEA linear programme for a single technology group.
#' Supports \code{rts = "fdh"} (free disposal hull, solved by exact
#' enumeration). When \code{slack = TRUE}, second-stage slacks against
#' the group's own reference set are stored as \code{slack_x} (n x m)
#' and \code{slack_y} (n x s) matrices.
#'
#' @keywords internal
#' @noRd
.fit_dea_group <- function(formula, data, orientation, rts,
                           slack = FALSE, ...) {

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

  out <- list(
    efficiency = efficiency,
    nobs = n,
    orientation = orientation,
    rts = rts,
    X = X_raw,
    Y = Y
  )

  if (isTRUE(slack)) {
    sl <- .dea_slacks(X_raw, Y, efficiency, orientation, rts, X_raw, Y)
    out$slack_x <- sl$slack_x
    out$slack_y <- sl$slack_y
  }

  out
}


#' Solve a single DEA LP
#'
#' Under \code{rts = "fdh"} the problem is solved by exact enumeration
#' of dominating reference points rather than by linear programming.
#'
#' @keywords internal
#' @noRd
.dea_solve_lp <- function(x_i, y_i, X, Y, orientation, rts) {
  n <- nrow(X)
  m <- ncol(X)
  s <- ncol(Y)

  if (rts == "fdh") {
    return(.fdh_radial(matrix(x_i, nrow = 1), matrix(y_i, nrow = 1),
                       orientation, X, Y)$efficiency)
  }

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

    # phi is bounded below by 0, not 1: cross-period evaluation of
    # super-efficient DMUs can have phi* < 1, and a lower bound of 1
    # would make those LPs infeasible. Same-period scores are
    # unaffected since phi* >= 1 whenever the DMU is in its own
    # reference set.
    lpSolveAPI::set.bounds(lp, lower = c(rep(0, n), 0))

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

  if (rts == "fdh") {
    return(.fdh_radial(X, Y, orientation, X_ref, Y_ref)$efficiency)
  }

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
    # Lower bound 0 on phi (not 1): see .dea_solve_lp
    lpSolveAPI::set.bounds(lp, lower = c(rep(0, n_ref), 0))

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


#' Radial FDH efficiency by exact enumeration
#'
#' Computes Farrell radial efficiency under free disposability without
#' convexity. Input orientation: the minimum over reference points that
#' weakly dominate on all outputs of the maximum input ratio. Output
#' orientation: the maximum over reference points that use no more of
#' any input of the minimum output ratio, reported as 1/phi (Farrell
#' convention, matching the LP path). If no reference point dominates
#' (possible in cross-period evaluation), the score is \code{NA} with a
#' warning, mirroring the infeasible-LP behaviour.
#'
#' @return A list with \code{efficiency} (numeric vector) and
#'   \code{peer} (integer row indices into the reference set of the
#'   single best dominating peer, \code{NA} where none exists).
#' @keywords internal
#' @noRd
.fdh_radial <- function(X, Y, orientation, X_ref, Y_ref) {
  n_eval <- nrow(X)
  s <- ncol(Y_ref)
  m <- ncol(X_ref)

  eff <- rep(NA_real_, n_eval)
  peer <- rep(NA_integer_, n_eval)

  for (i in seq_len(n_eval)) {
    if (orientation == "input") {
      # Reference points that weakly dominate on all outputs
      dom <- which(colSums(t(Y_ref) >= Y[i, ]) == s)
      if (length(dom) == 0L) {
        warning("No dominating FDH reference point for a DMU.",
                call. = FALSE)
        next
      }
      scores <- apply(X_ref[dom, , drop = FALSE], 1,
                      function(r) max(r / X[i, ]))
      k <- which.min(scores)
      eff[i] <- scores[k]
      peer[i] <- dom[k]
    } else {
      # Reference points that use no more of any input
      dom <- which(colSums(t(X_ref) <= X[i, ]) == m)
      if (length(dom) == 0L) {
        warning("No dominating FDH reference point for a DMU.",
                call. = FALSE)
        next
      }
      phis <- apply(Y_ref[dom, , drop = FALSE], 1,
                    function(r) min(r / Y[i, ]))
      k <- which.max(phis)
      eff[i] <- 1 / phis[k]
      peer[i] <- dom[k]
    }
  }

  list(efficiency = eff, peer = peer)
}


#' Internal: Estimate DEA-based metafrontier
#'
#' Computes group DEA and pooled (meta) DEA, then derives TGR.
#' Supports \code{rts = "fdh"}. When \code{slack = TRUE}, second-stage
#' slacks against the pooled reference set are stored as
#' \code{slack_x_meta} and \code{slack_y_meta}.
#'
#' @keywords internal
#' @noRd
.estimate_dea_metafrontier <- function(formula, data, group_vec,
                                       group_levels, group_models,
                                       orientation, rts, slack = FALSE) {

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

  out <- list(
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

  if (isTRUE(slack)) {
    sl <- .dea_slacks(X_raw, Y, te_meta, orientation, rts, X_raw, Y)
    out$slack_x_meta <- sl$slack_x
    out$slack_y_meta <- sl$slack_y
  }

  out
}


#' Second-stage slack maximisation
#'
#' Standard two-stage slack computation with the radial score held
#' fixed. Input orientation: max sum(s_x) + sum(s_y) subject to
#' X_ref' lambda + s_x = theta_i * x_i and Y_ref' lambda - s_y = y_i.
#' Output orientation: X_ref' lambda + s_x = x_i and
#' Y_ref' lambda - s_y = phi_i * y_i. The equalities are implemented
#' as pairs of inequalities relaxed by a small tolerance because the
#' radially scaled targets are floating-point products. Under
#' \code{rts = "fdh"} slacks are measured against the single
#' dominating peer identified by the enumeration.
#'
#' @param scores Farrell efficiency scores as returned by the radial
#'   solvers (input orientation: theta; output orientation: 1/phi).
#' @return A list with matrices \code{slack_x} (n x m) and
#'   \code{slack_y} (n x s); rows are \code{NA} where the radial score
#'   is \code{NA}.
#' @keywords internal
#' @noRd
.dea_slacks <- function(X, Y, scores, orientation, rts, X_ref, Y_ref) {
  n <- nrow(X)
  m <- ncol(X)
  s <- ncol(Y)
  n_ref <- nrow(X_ref)

  slack_x <- matrix(NA_real_, n, m)
  slack_y <- matrix(NA_real_, n, s)

  if (rts == "fdh") {
    fdh <- .fdh_radial(X, Y, orientation, X_ref, Y_ref)
    for (i in seq_len(n)) {
      j <- fdh$peer[i]
      if (is.na(j)) next
      if (orientation == "input") {
        slack_x[i, ] <- pmax(scores[i] * X[i, ] - X_ref[j, ], 0)
        slack_y[i, ] <- pmax(Y_ref[j, ] - Y[i, ], 0)
      } else {
        slack_x[i, ] <- pmax(X[i, ] - X_ref[j, ], 0)
        slack_y[i, ] <- pmax(Y_ref[j, ] - Y[i, ] / scores[i], 0)
      }
    }
    return(list(slack_x = slack_x, slack_y = slack_y))
  }

  n_vars <- n_ref + m + s  # lambdas, then s_x, then s_y
  tol <- 1e-9

  for (i in seq_len(n)) {
    if (is.na(scores[i])) next

    if (orientation == "input") {
      x_target <- scores[i] * X[i, ]
      y_target <- Y[i, ]
    } else {
      x_target <- X[i, ]
      y_target <- Y[i, ] / scores[i]
    }

    lp <- lpSolveAPI::make.lp(0, n_vars)
    lpSolveAPI::set.objfn(lp, c(rep(0, n_ref), rep(1, m + s)))
    lpSolveAPI::lp.control(lp, sense = "max", verbose = "neutral")

    # X_ref' lambda + s_x = x_target (relaxed equality)
    for (mm in seq_len(m)) {
      coef <- c(X_ref[, mm], as.numeric(seq_len(m) == mm), rep(0, s))
      eps <- tol * (1 + abs(x_target[mm]))
      lpSolveAPI::add.constraint(lp, coef, "<=", x_target[mm] + eps)
      lpSolveAPI::add.constraint(lp, coef, ">=", x_target[mm] - eps)
    }
    # Y_ref' lambda - s_y = y_target (relaxed equality)
    for (ss in seq_len(s)) {
      coef <- c(Y_ref[, ss], rep(0, m), -as.numeric(seq_len(s) == ss))
      eps <- tol * (1 + abs(y_target[ss]))
      lpSolveAPI::add.constraint(lp, coef, "<=", y_target[ss] + eps)
      lpSolveAPI::add.constraint(lp, coef, ">=", y_target[ss] - eps)
    }

    if (rts == "vrs") {
      lpSolveAPI::add.constraint(lp, c(rep(1, n_ref), rep(0, m + s)), "=", 1)
    } else if (rts == "drs") {
      lpSolveAPI::add.constraint(lp, c(rep(1, n_ref), rep(0, m + s)), "<=", 1)
    } else if (rts == "irs") {
      lpSolveAPI::add.constraint(lp, c(rep(1, n_ref), rep(0, m + s)), ">=", 1)
    }

    status <- lpSolveAPI::solve.lpExtPtr(lp)
    if (status == 0) {
      sol <- lpSolveAPI::get.variables(lp)
      slack_x[i, ] <- pmax(sol[n_ref + seq_len(m)], 0)
      slack_y[i, ] <- pmax(sol[n_ref + m + seq_len(s)], 0)
    } else {
      warning("Slack LP infeasible for a DMU.", call. = FALSE)
    }
  }

  list(slack_x = slack_x, slack_y = slack_y)
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
#' @param rts returns to scale assumption; \code{"fdh"} is solved as a
#'   mixed-integer programme with binary intensity variables.
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
  } else if (rts == "fdh") {
    # Free disposal hull: exactly one reference point, so the lambdas
    # are binary and sum to one (mixed-integer programme)
    lpSolveAPI::add.constraint(lp, c(rep(1, n), 0), "=", 1)
    lpSolveAPI::set.type(lp, columns = seq_len(n), type = "binary")
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


#' Build DDF direction matrices
#'
#' Accepts a character preset (\code{"proportional"} g = (x_i, y_i),
#' \code{"output"} g = (0, y_i), \code{"input"} g = (x_i, 0)), a
#' numeric vector of length m + s giving a common direction (first m
#' elements g_x, last s elements g_y), or a numeric n x (m + s) matrix
#' whose row i is firm i's direction.
#'
#' @return A list with \code{g_x} (n x m), \code{g_y} (n x s) and
#'   \code{numeric} (\code{TRUE} for user-supplied numeric directions).
#' @keywords internal
#' @noRd
.ddf_direction_mats <- function(direction, X, Y) {
  n <- nrow(X)
  m <- ncol(X)
  s <- ncol(Y)

  if (is.character(direction) && length(direction) == 1L) {
    if (direction == "proportional") {
      return(list(g_x = X, g_y = Y, numeric = FALSE))
    } else if (direction == "output") {
      return(list(g_x = matrix(0, n, m), g_y = Y, numeric = FALSE))
    } else if (direction == "input") {
      return(list(g_x = X, g_y = matrix(0, n, s), numeric = FALSE))
    }
    stop("Unknown direction: ", direction, call. = FALSE)
  }

  if (is.numeric(direction)) {
    if (is.matrix(direction)) {
      if (nrow(direction) != n || ncol(direction) != m + s) {
        stop("A direction matrix must have dimensions n x (m + s), here ",
             n, " x ", m + s, ".", call. = FALSE)
      }
      g <- direction
    } else {
      if (length(direction) != m + s) {
        stop("A direction vector must have length m + s = ", m + s, ".",
             call. = FALSE)
      }
      g <- matrix(direction, n, m + s, byrow = TRUE)
    }
    if (any(!is.finite(g)) || any(g < 0)) {
      stop("Numeric directions must be finite and non-negative.",
           call. = FALSE)
    }
    if (any(rowSums(g) == 0)) {
      stop("Each direction must have at least one positive element.",
           call. = FALSE)
    }
    return(list(g_x = g[, seq_len(m), drop = FALSE],
                g_y = g[, m + seq_len(s), drop = FALSE],
                numeric = TRUE))
  }

  stop("'direction' must be a character preset, a numeric vector of ",
       "length m + s, or a numeric n x (m + s) matrix.", call. = FALSE)
}


#' DDF-based group efficiency
#'
#' \code{direction} may be a character preset, a numeric vector of
#' length m + s, or a numeric matrix with one row per observation of
#' \code{data} (see \code{.ddf_direction_mats}).
#'
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

  if (is.matrix(y_raw)) {
    Y <- y_raw
  } else {
    Y <- matrix(y_raw, ncol = 1)
  }

  # Direction vectors
  dirs <- .ddf_direction_mats(direction, X_raw, Y)

  efficiency <- .ddf_batch(X_raw, Y, dirs$g_x, dirs$g_y, rts)

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
#'
#' \code{direction} may be a character preset, a numeric vector of
#' length m + s, or a numeric n x (m + s) matrix (see
#' \code{.ddf_direction_mats}). Convention: the multiplicative
#' conversion te = 1/(1 + beta) is only meaningful for the character
#' presets, where the direction scales with the observation. For
#' user-supplied numeric directions \code{te_group}, \code{te_meta}
#' and \code{tgr} are set to \code{NA} and the additive fields
#' \code{beta_group}, \code{beta_meta} and
#' \code{ddf_gap = beta_meta - beta_group} (the additive technology
#' gap, non-negative up to solver tolerance) carry the results. The
#' additive fields are populated for the presets too, for
#' comparability.
#'
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

  # Direction vectors
  dirs <- .ddf_direction_mats(direction, X_raw, Y)

  # Group-level DDF
  beta_group <- numeric(n)
  for (g in group_levels) {
    idx <- which(group_vec == g)
    beta_group[idx] <- group_models[[g]]$efficiency
  }

  # Pooled DDF (metafrontier)
  beta_meta <- .ddf_batch(X_raw, Y, dirs$g_x, dirs$g_y, rts)

  if (dirs$numeric) {
    # te = 1/(1 + beta) is not meaningful for arbitrary directions;
    # only the additive decomposition below applies
    te_group <- rep(NA_real_, n)
    te_meta <- rep(NA_real_, n)
    tgr <- rep(NA_real_, n)
  } else {
    # DDF TGR: additive decomposition
    # beta_meta = beta_group + TGR_DDF
    # TGR = beta_meta - beta_group (additive gap)
    # For compatibility, also compute ratio-based TE
    te_group <- 1 / (1 + beta_group)
    te_meta <- 1 / (1 + beta_meta)
    tgr <- te_meta / te_group
  }

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
    ddf_tgr = beta_meta - beta_group,
    ddf_gap = beta_meta - beta_group
  )
}


#' Batch hyperbolic graph efficiency
#'
#' Computes hyperbolic efficiency gamma in (0, 1]: the smallest gamma
#' such that (gamma * x_i, y_i / gamma) remains in the technology.
#' Under CRS the exact closed form gamma = sqrt(theta) is used, where
#' theta is the input-oriented radial score against the same reference
#' set. Under FDH the minimum over reference points of the smallest
#' feasible gamma is exact. Under vrs/drs/irs gamma is found by
#' bisection on feasibility LPs, reusing a single LP object across
#' DMUs and bisection steps.
#'
#' @keywords internal
#' @noRd
.hyperbolic_batch <- function(X, Y, rts, X_ref = NULL, Y_ref = NULL) {
  if (is.null(X_ref)) X_ref <- X
  if (is.null(Y_ref)) Y_ref <- Y

  n_eval <- nrow(X)
  m <- ncol(X)
  s <- ncol(Y)
  n_ref <- nrow(X_ref)

  if (rts == "crs") {
    theta <- .dea_batch_fast(X, Y, "input", "crs", X_ref, Y_ref)
    return(sqrt(theta))
  }

  if (rts == "fdh") {
    # Reference point j admits gamma iff gamma >= X_jm / x_im for all
    # inputs and gamma >= y_is / Y_js for all outputs
    gamma <- rep(NA_real_, n_eval)
    for (i in seq_len(n_eval)) {
      req_x <- apply(X_ref, 1, function(r) max(r / X[i, ]))
      req_y <- apply(Y_ref, 1, function(r) max(Y[i, ] / r))
      gamma[i] <- min(pmax(req_x, req_y))
    }
    return(gamma)
  }

  # vrs / drs / irs: bisection on feasibility LPs
  lp <- lpSolveAPI::make.lp(0, n_ref)
  lpSolveAPI::set.objfn(lp, rep(0, n_ref))
  lpSolveAPI::lp.control(lp, sense = "min", verbose = "neutral")
  for (mm in seq_len(m)) {
    lpSolveAPI::add.constraint(lp, X_ref[, mm], "<=", 0)
  }
  for (ss in seq_len(s)) {
    lpSolveAPI::add.constraint(lp, Y_ref[, ss], ">=", 0)
  }
  if (rts == "vrs") {
    lpSolveAPI::add.constraint(lp, rep(1, n_ref), "=", 1)
  } else if (rts == "drs") {
    lpSolveAPI::add.constraint(lp, rep(1, n_ref), "<=", 1)
  } else if (rts == "irs") {
    lpSolveAPI::add.constraint(lp, rep(1, n_ref), ">=", 1)
  }

  feasible <- function(gamma, x_i, y_i) {
    for (mm in seq_len(m)) {
      lpSolveAPI::set.rhs(lp, gamma * x_i[mm], mm)
    }
    for (ss in seq_len(s)) {
      lpSolveAPI::set.rhs(lp, y_i[ss] / gamma, m + ss)
    }
    lpSolveAPI::solve.lpExtPtr(lp) == 0
  }

  gamma <- rep(NA_real_, n_eval)
  for (i in seq_len(n_eval)) {
    if (!feasible(1, X[i, ], Y[i, ])) {
      warning("Hyperbolic feasibility fails at gamma = 1 for a DMU.",
              call. = FALSE)
      next
    }
    lo <- 0
    hi <- 1
    for (iter in seq_len(40L)) {
      if (hi - lo < 1e-8) break
      mid <- (lo + hi) / 2
      if (feasible(mid, X[i, ], Y[i, ])) {
        hi <- mid
      } else {
        lo <- mid
      }
    }
    gamma[i] <- hi
  }
  gamma
}


#' Hyperbolic group efficiency
#'
#' Mirrors \code{.fit_dea_group} with hyperbolic (graph) efficiency:
#' the firm is projected to (gamma * x, y / gamma) with te = gamma.
#'
#' @keywords internal
#' @noRd
.fit_hyperbolic_group <- function(formula, data, rts, ...) {

  mf <- model.frame(formula, data = data)
  y_raw <- model.response(mf)
  X_raw <- model.matrix(formula, data = data, rhs = 1)

  if (colnames(X_raw)[1] == "(Intercept)") {
    X_raw <- X_raw[, -1, drop = FALSE]
  }

  if (is.matrix(y_raw)) {
    Y <- y_raw
  } else {
    Y <- matrix(y_raw, ncol = 1)
  }

  efficiency <- .hyperbolic_batch(X_raw, Y, rts)

  list(
    efficiency = efficiency,
    nobs = nrow(X_raw),
    rts = rts,
    X = X_raw,
    Y = Y
  )
}


#' Hyperbolic metafrontier estimation
#'
#' Mirrors \code{.estimate_dea_metafrontier}: group gammas are read
#' from the fitted group models, the meta gamma is computed against
#' the pooled reference set, and tgr = gamma_meta / gamma_group lies
#' in (0, 1] since the pooled reference set is a superset.
#'
#' @keywords internal
#' @noRd
.estimate_hyperbolic_metafrontier <- function(formula, data, group_vec,
                                              group_levels, group_models,
                                              rts) {

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

  # Pooled hyperbolic efficiency (metafrontier)
  te_meta <- .hyperbolic_batch(X_raw, Y, rts)

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
