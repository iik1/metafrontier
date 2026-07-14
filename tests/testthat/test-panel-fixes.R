# Regression tests for panel SFA fixes: row order, NA handling,
# global-T decay anchoring, and unbalanced BC92 estimation

.panel_test_group <- function(seed = 42, eta = 0.1, n_firms = 30,
                              n_periods = 5) {
  sim <- simulate_panel_metafrontier(
    n_groups = 2, n_firms_per_group = n_firms,
    n_periods = n_periods, eta = eta, seed = seed
  )
  g1 <- sim$data[sim$data$group == "G1", ]
  # Firm IDs F1...F30: string sort order (F1, F10, F11, ...) differs
  # from numeric order, which is what exposed the released row-order bug
  g1$firm <- sub("^G1_", "", g1$firm)
  list(data = g1, params = sim$params)
}

.panel_formula <- function() Formula::Formula(log_y ~ log_x1 + log_x2)
.panel_info <- function() list(id = "firm", time = "year")

test_that("panel fit is invariant to row order", {
  g1 <- .panel_test_group()$data
  set.seed(7)
  g1_shuffled <- g1[sample(nrow(g1)), ]

  fit_sorted <- metafrontier:::.fit_sfa_panel_group(
    .panel_formula(), g1, "hnormal", "bc92", .panel_info(), list()
  )
  fit_shuffled <- metafrontier:::.fit_sfa_panel_group(
    .panel_formula(), g1_shuffled, "hnormal", "bc92", .panel_info(), list()
  )

  expect_equal(fit_sorted$coefficients, fit_shuffled$coefficients,
               tolerance = 1e-8)

  # Efficiency must agree per (firm, time) key, i.e. returned vectors
  # are in input-row order for both fits
  key_sorted <- paste(g1$firm, g1$year)
  key_shuffled <- paste(g1_shuffled$firm, g1_shuffled$year)
  expect_equal(fit_shuffled$efficiency,
               fit_sorted$efficiency[match(key_shuffled, key_sorted)],
               tolerance = 1e-8)
})

test_that("NA rows are dropped consistently across y, X and firm index", {
  g1 <- .panel_test_group()$data
  g1_na <- g1
  g1_na$log_x1[5] <- NA

  fit_na <- metafrontier:::.fit_sfa_panel_group(
    .panel_formula(), g1_na, "hnormal", "bc92", .panel_info(), list()
  )
  fit_dropped <- metafrontier:::.fit_sfa_panel_group(
    .panel_formula(), g1[-5, ], "hnormal", "bc92", .panel_info(), list()
  )

  expect_length(fit_na$efficiency, nrow(g1) - 1)
  expect_identical(as.character(fit_na$firms), as.character(g1$firm[-5]))
  expect_equal(fit_na$efficiency, fit_dropped$efficiency, tolerance = 1e-10)
  expect_equal(fit_na$coefficients, fit_dropped$coefficients,
               tolerance = 1e-10)
})

test_that("unbalanced BC92 converges and recovers eta", {
  gen <- .panel_test_group(seed = 456, eta = 0.1, n_firms = 50,
                           n_periods = 6)
  g1 <- gen$data
  set.seed(3)
  g1_unbal <- g1[-sample(nrow(g1), 30), ]

  fit <- metafrontier:::.fit_sfa_panel_group(
    .panel_formula(), g1_unbal, "hnormal", "bc92", .panel_info(), list()
  )

  expect_identical(fit$convergence, 0L)
  expect_true(all(fit$efficiency > 0 & fit$efficiency <= 1))
  # Loose tolerance: eta is weakly identified in short panels
  expect_lt(abs(fit$eta - gen$params$eta), 0.15)
})

test_that("BC92 decay is anchored at the global final period", {
  g1 <- .panel_test_group(seed = 42, eta = 0.15)$data
  # Remove firm F1's final period so its own last period < global T
  g1_unbal <- g1[!(g1$firm == "F1" & g1$year == max(g1$year)), ]

  fit <- metafrontier:::.fit_sfa_panel_group(
    .panel_formula(), g1_unbal, "hnormal", "bc92", .panel_info(), list()
  )

  # Recompute BC92 efficiency by hand with d_t = exp(-eta * (t - T)),
  # T the GLOBAL final period; must match the stored values exactly
  T_max <- max(g1_unbal$year)
  sv <- fit$sigma_v
  su <- fit$sigma_u
  eta_hat <- fit$eta
  idx_f1 <- which(fit$firms == "F1")
  eps_f1 <- fit$residuals[idx_f1]
  d_t <- exp(-eta_hat * (g1_unbal$year[g1_unbal$firm == "F1"] - T_max))

  sigma_star2 <- 1 / (1 / su^2 + sum(d_t^2) / sv^2)
  sigma_star <- sqrt(sigma_star2)
  mu_star <- -sigma_star2 * sum(eps_f1 * d_t) / sv^2
  ratio <- mu_star / sigma_star

  te_bc88_hand <- pnorm(ratio - d_t * sigma_star) / pnorm(ratio) *
    exp(-d_t * mu_star + 0.5 * d_t^2 * sigma_star2)

  expect_equal(fit$efficiency_bc88[idx_f1], te_bc88_hand, tolerance = 1e-10)

  # With eta_hat > 0 inefficiency decays over time, so F1's efficiency
  # path must be monotone in the direction implied by the sign of eta
  eff_path <- fit$efficiency_bc88[idx_f1][order(g1_unbal$year[g1_unbal$firm == "F1"])]
  if (fit$eta > 0) {
    expect_true(all(diff(eff_path) > 0))
  } else if (fit$eta < 0) {
    expect_true(all(diff(eff_path) < 0))
  }
})

test_that("panel estimator argument works like the cross-sectional one", {
  g1 <- .panel_test_group()$data

  fit_def <- metafrontier:::.fit_sfa_panel_group(
    .panel_formula(), g1, "hnormal", "bc92", .panel_info(), list()
  )
  fit_jlms <- metafrontier:::.fit_sfa_panel_group(
    .panel_formula(), g1, "hnormal", "bc92", .panel_info(), list(),
    estimator = "jlms"
  )

  expect_identical(fit_def$estimator, "bc88")
  expect_identical(fit_def$efficiency, fit_def$efficiency_bc88)
  expect_identical(fit_jlms$estimator, "jlms")
  expect_identical(fit_jlms$efficiency, fit_jlms$efficiency_jlms)
  expect_true(all(fit_def$efficiency_bc88 > 0 & fit_def$efficiency_bc88 <= 1))
  expect_gt(cor(fit_def$efficiency_bc88, fit_def$efficiency_jlms), 0.99)
})
