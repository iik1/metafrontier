# Tests for v0.3.0 DEA features: phi-bound fix, FDH, user-defined
# DDF directions, hyperbolic efficiency, and second-stage slacks

# Small positive-valued production dataset shared across tests
.make_dea_feature_data <- function(n = 40, seed = 123) {
  set.seed(seed)
  x1 <- runif(n, 1, 10)
  x2 <- runif(n, 1, 10)
  grp <- rep(c("A", "B"), each = n / 2)
  tech <- ifelse(grp == "A", 1, 0.8)
  y <- tech * (x1^0.4 * x2^0.4) * runif(n, 0.6, 1)
  data.frame(y = y, x1 = x1, x2 = x2, group = grp)
}


# ---- (i) phi lower bound fix ----

test_that("output-oriented LP solves for super-efficient DMUs (phi < 1)", {
  # Reference technology: y = x under CRS. The evaluated point (1, 2)
  # lies strictly above it, so phi* = 0.5 < 1. With the old lower
  # bound phi >= 1 this LP was infeasible.
  X_ref <- matrix(c(1, 2, 3, 4, 5), ncol = 1)
  Y_ref <- matrix(c(1, 2, 3, 4, 5), ncol = 1)

  te <- metafrontier:::.dea_solve_lp(1, 2, X_ref, Y_ref, "output", "crs")
  expect_false(is.na(te))
  expect_equal(te, 2, tolerance = 1e-8)  # Farrell TE = 1/phi

  # The batch solver takes the same cross-period path
  te_b <- metafrontier:::.dea_batch_fast(matrix(1, 1, 1), matrix(2, 1, 1),
                                         "output", "crs",
                                         X_ref = X_ref, Y_ref = Y_ref)
  expect_equal(te_b, 2, tolerance = 1e-8)
})


test_that("same-period output-oriented scores remain in (0, 1]", {
  dat <- .make_dea_feature_data()
  X <- cbind(dat$x1, dat$x2)
  Y <- matrix(dat$y, ncol = 1)

  te <- metafrontier:::.dea_batch_fast(X, Y, "output", "crs")
  expect_true(all(te > 0))
  expect_true(all(te <= 1 + 1e-8))
})


# ---- (ii) FDH ----

test_that("FDH radial efficiency dominates VRS in both orientations", {
  skip_on_cran()
  dat <- .make_dea_feature_data()
  X <- cbind(dat$x1, dat$x2)
  Y <- matrix(dat$y, ncol = 1)

  for (orient in c("input", "output")) {
    e_fdh <- metafrontier:::.dea_batch_fast(X, Y, orient, "fdh")
    e_vrs <- metafrontier:::.dea_batch_fast(X, Y, orient, "vrs")
    expect_true(all(e_fdh >= e_vrs - 1e-8))
    expect_true(all(e_fdh > 0))
    expect_true(all(e_fdh <= 1 + 1e-8))
  }
})


test_that("FDH metafrontier TGR lies in (0, 1]", {
  skip_on_cran()
  dat <- .make_dea_feature_data()
  f <- Formula::Formula(y ~ x1 + x2)
  gvec <- factor(dat$group)
  glev <- levels(gvec)

  gm <- lapply(glev, function(g) {
    metafrontier:::.fit_dea_group(f, dat[gvec == g, ], "input", "fdh")
  })
  names(gm) <- glev

  res <- metafrontier:::.estimate_dea_metafrontier(f, dat, gvec, glev, gm,
                                                   "input", "fdh")
  expect_true(all(res$tgr > 0))
  expect_true(all(res$tgr <= 1 + 1e-8))
})


test_that("FDH returns NA with a warning when no reference point dominates", {
  # Cross-period evaluation: no reference point can produce y = 10
  expect_warning(
    v <- metafrontier:::.dea_batch_fast(matrix(0.5, 1, 1), matrix(10, 1, 1),
                                        "input", "fdh",
                                        X_ref = matrix(1:3, ncol = 1),
                                        Y_ref = matrix(1:3, ncol = 1)),
    "No dominating FDH reference point"
  )
  expect_true(is.na(v))
})


test_that("DDF with FDH (binary MIP) is no larger than the VRS beta", {
  skip_on_cran()
  dat <- .make_dea_feature_data(n = 20)
  X <- cbind(dat$x1, dat$x2)
  Y <- matrix(dat$y, ncol = 1)

  b_fdh <- metafrontier:::.dea_solve_ddf(X[1, ], Y[1, ], X, Y,
                                         g_x = X[1, ], g_y = Y[1, ],
                                         rts = "fdh")
  b_vrs <- metafrontier:::.dea_solve_ddf(X[1, ], Y[1, ], X, Y,
                                         g_x = X[1, ], g_y = Y[1, ],
                                         rts = "vrs")
  expect_true(b_fdh >= -1e-8)
  expect_true(b_fdh <= b_vrs + 1e-8)
})


# ---- (iii) user-defined DDF directions ----

test_that("numeric direction vector runs and reports additive gaps", {
  skip_on_cran()
  dat <- .make_dea_feature_data()
  f <- Formula::Formula(y ~ x1 + x2)
  gvec <- factor(dat$group)
  glev <- levels(gvec)
  gdir <- c(mean(dat$x1), mean(dat$x2), mean(dat$y))

  gm <- lapply(glev, function(g) {
    metafrontier:::.fit_ddf_group(f, dat[gvec == g, ], "crs", gdir)
  })
  names(gm) <- glev

  res <- metafrontier:::.estimate_ddf_metafrontier(f, dat, gvec, glev, gm,
                                                   "crs", gdir)

  # te = 1/(1 + beta) is undefined for arbitrary numeric directions
  expect_true(all(is.na(res$te_group)))
  expect_true(all(is.na(res$te_meta)))
  expect_true(all(is.na(res$tgr)))

  # Additive fields carry the results; the gap is non-negative because
  # the pooled reference set is a superset
  expect_length(res$beta_group, nrow(dat))
  expect_length(res$beta_meta, nrow(dat))
  expect_true(all(res$beta_meta >= res$beta_group - 1e-8))
  expect_equal(res$ddf_gap, res$beta_meta - res$beta_group)
})


test_that("firm-specific direction matrix runs and matches common vector", {
  skip_on_cran()
  dat <- .make_dea_feature_data()
  f <- Formula::Formula(y ~ x1 + x2)
  gvec <- factor(dat$group)
  glev <- levels(gvec)
  n <- nrow(dat)
  gdir <- c(mean(dat$x1), mean(dat$x2), mean(dat$y))
  dmat <- matrix(rep(gdir, n), n, 3, byrow = TRUE)

  gm <- lapply(glev, function(g) {
    metafrontier:::.fit_ddf_group(f, dat[gvec == g, ], "crs",
                                  dmat[gvec == g, , drop = FALSE])
  })
  names(gm) <- glev

  res_mat <- metafrontier:::.estimate_ddf_metafrontier(f, dat, gvec, glev,
                                                       gm, "crs", dmat)

  gm_vec <- lapply(glev, function(g) {
    metafrontier:::.fit_ddf_group(f, dat[gvec == g, ], "crs", gdir)
  })
  names(gm_vec) <- glev
  res_vec <- metafrontier:::.estimate_ddf_metafrontier(f, dat, gvec, glev,
                                                       gm_vec, "crs", gdir)

  expect_equal(res_mat$beta_meta, res_vec$beta_meta, tolerance = 1e-10)
  expect_equal(res_mat$beta_group, res_vec$beta_group, tolerance = 1e-10)
})


test_that("character presets still populate te and the additive fields", {
  dat <- .make_dea_feature_data(n = 20)
  f <- Formula::Formula(y ~ x1 + x2)
  gvec <- factor(dat$group)
  glev <- levels(gvec)

  gm <- lapply(glev, function(g) {
    metafrontier:::.fit_ddf_group(f, dat[gvec == g, ], "crs", "proportional")
  })
  names(gm) <- glev

  res <- metafrontier:::.estimate_ddf_metafrontier(f, dat, gvec, glev, gm,
                                                   "crs", "proportional")
  expect_false(any(is.na(res$te_meta)))
  expect_false(any(is.na(res$tgr)))
  expect_equal(res$ddf_gap, res$beta_meta - res$beta_group)
  expect_true(all(res$ddf_gap >= -1e-8))
})


test_that("malformed numeric directions are rejected", {
  X <- matrix(1:6, 3, 2)
  Y <- matrix(1:3, 3, 1)
  expect_error(metafrontier:::.ddf_direction_mats(c(1, 2), X, Y),
               "length m \\+ s")
  expect_error(metafrontier:::.ddf_direction_mats(matrix(1, 2, 3), X, Y),
               "n x \\(m \\+ s\\)")
  expect_error(metafrontier:::.ddf_direction_mats(c(-1, 1, 1), X, Y),
               "non-negative")
  expect_error(metafrontier:::.ddf_direction_mats(c(0, 0, 0), X, Y),
               "at least one positive")
})


# ---- (iv) hyperbolic efficiency ----

test_that("hyperbolic CRS closed form agrees with direct bisection", {
  skip_on_cran()
  dat <- .make_dea_feature_data(n = 20)
  X <- cbind(dat$x1, dat$x2)
  Y <- matrix(dat$y, ncol = 1)
  n <- nrow(X)

  g_closed <- metafrontier:::.hyperbolic_batch(X, Y, "crs")
  expect_true(all(g_closed > 0))
  expect_true(all(g_closed <= 1 + 1e-8))

  # Independent bisection on CRS feasibility LPs:
  # gamma feasible iff exists lambda >= 0 with
  # X' lambda <= gamma * x_i and Y' lambda >= y_i / gamma
  crs_feasible <- function(gamma, i) {
    lp <- lpSolveAPI::make.lp(0, n)
    lpSolveAPI::set.objfn(lp, rep(0, n))
    lpSolveAPI::lp.control(lp, sense = "min", verbose = "neutral")
    for (mm in 1:2) {
      lpSolveAPI::add.constraint(lp, X[, mm], "<=", gamma * X[i, mm])
    }
    lpSolveAPI::add.constraint(lp, Y[, 1], ">=", Y[i, 1] / gamma)
    lpSolveAPI::solve.lpExtPtr(lp) == 0
  }

  for (i in c(1, 5, 10, 15, 20)) {
    lo <- 0
    hi <- 1
    for (iter in 1:40) {
      if (hi - lo < 1e-9) break
      mid <- (lo + hi) / 2
      if (crs_feasible(mid, i)) hi <- mid else lo <- mid
    }
    expect_equal(g_closed[i], hi, tolerance = 1e-6)
  }
})


test_that("hyperbolic metafrontier works under VRS with TGR <= 1", {
  skip_on_cran()
  dat <- .make_dea_feature_data()
  f <- Formula::Formula(y ~ x1 + x2)
  gvec <- factor(dat$group)
  glev <- levels(gvec)

  gm <- lapply(glev, function(g) {
    metafrontier:::.fit_hyperbolic_group(f, dat[gvec == g, ], "vrs")
  })
  names(gm) <- glev

  res <- metafrontier:::.estimate_hyperbolic_metafrontier(f, dat, gvec,
                                                          glev, gm, "vrs")

  expect_true(all(res$te_group > 0))
  expect_true(all(res$te_group <= 1 + 1e-8))
  expect_true(all(res$te_meta > 0))
  expect_true(all(res$te_meta <= 1 + 1e-8))
  # Bisection tolerance (1e-8 per gamma) is amplified in the ratio
  expect_true(all(res$tgr > 0))
  expect_true(all(res$tgr <= 1 + 1e-6))
  expect_identical(res$meta_convergence, 0L)
  expect_null(res$meta_coef)
})


test_that("hyperbolic FDH and VRS envelop CRS", {
  skip_on_cran()
  dat <- .make_dea_feature_data(n = 20)
  X <- cbind(dat$x1, dat$x2)
  Y <- matrix(dat$y, ncol = 1)

  g_crs <- metafrontier:::.hyperbolic_batch(X, Y, "crs")
  g_vrs <- metafrontier:::.hyperbolic_batch(X, Y, "vrs")
  g_fdh <- metafrontier:::.hyperbolic_batch(X, Y, "fdh")

  expect_true(all(g_vrs >= g_crs - 1e-6))
  expect_true(all(g_fdh >= g_vrs - 1e-6))
  expect_true(all(g_fdh <= 1 + 1e-8))
})


# ---- (v) slacks ----

test_that("slack = TRUE returns non-negative slack matrices of correct dims", {
  skip_on_cran()
  dat <- .make_dea_feature_data()
  f <- Formula::Formula(y ~ x1 + x2)
  gvec <- factor(dat$group)
  glev <- levels(gvec)

  gm <- lapply(glev, function(g) {
    metafrontier:::.fit_dea_group(f, dat[gvec == g, ], "input", "vrs",
                                  slack = TRUE)
  })
  names(gm) <- glev

  for (g in glev) {
    n_g <- gm[[g]]$nobs
    expect_identical(dim(gm[[g]]$slack_x), c(n_g, 2L))
    expect_identical(dim(gm[[g]]$slack_y), c(n_g, 1L))
    expect_true(all(gm[[g]]$slack_x >= 0))
    expect_true(all(gm[[g]]$slack_y >= 0))
  }

  res <- metafrontier:::.estimate_dea_metafrontier(f, dat, gvec, glev, gm,
                                                   "input", "vrs",
                                                   slack = TRUE)
  expect_identical(dim(res$slack_x_meta), c(nrow(dat), 2L))
  expect_identical(dim(res$slack_y_meta), c(nrow(dat), 1L))
  expect_true(all(res$slack_x_meta >= 0))
  expect_true(all(res$slack_y_meta >= 0))
})


test_that("a radially efficient but dominated DMU shows positive slack", {
  # Classic example: A = (1, 4), B = (4, 1), C = (2, 2) efficient;
  # D = (1, 5) has theta* = 1 under VRS but slack of 1 in input 2
  X <- matrix(c(1, 4, 1, 2,
                4, 1, 5, 2), ncol = 2)
  Y <- matrix(1, 4, 1)

  theta <- metafrontier:::.dea_batch_fast(X, Y, "input", "vrs")
  expect_equal(theta[3], 1, tolerance = 1e-8)

  sl <- metafrontier:::.dea_slacks(X, Y, theta, "input", "vrs", X, Y)
  expect_true(all(sl$slack_x >= 0))
  expect_true(all(sl$slack_y >= 0))
  expect_equal(sl$slack_x[3, 2], 1, tolerance = 1e-6)
  expect_equal(sl$slack_x[3, 1], 0, tolerance = 1e-6)
})


test_that("FDH slacks are measured against the dominating peer", {
  X <- matrix(c(1, 4, 1, 2,
                4, 1, 5, 2), ncol = 2)
  Y <- matrix(1, 4, 1)

  theta <- metafrontier:::.dea_batch_fast(X, Y, "input", "fdh")
  sl <- metafrontier:::.dea_slacks(X, Y, theta, "input", "fdh", X, Y)

  expect_true(all(sl$slack_x >= 0))
  expect_true(all(sl$slack_y >= 0))
  # D = (1, 5) is dominated by peer A = (1, 4): slack of 1 in input 2
  expect_equal(sl$slack_x[3, 2], 1, tolerance = 1e-8)
})


test_that("output-oriented slacks are non-negative", {
  skip_on_cran()
  X <- matrix(c(1, 4, 1, 2,
                4, 1, 5, 2), ncol = 2)
  Y <- matrix(1, 4, 1)

  te <- metafrontier:::.dea_batch_fast(X, Y, "output", "vrs")
  sl <- metafrontier:::.dea_slacks(X, Y, te, "output", "vrs", X, Y)
  expect_true(all(sl$slack_x >= 0))
  expect_true(all(sl$slack_y >= 0))
})
