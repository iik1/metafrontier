# Spatial (SAR) SFA metafrontier: y = rho W y + x'beta + v - u, with W
# restricted to each group's own units

sp_sim <- simulate_metafrontier(n_groups = 2, n_per_group = 120,
                                sigma_v = 0.1, seed = 17)
sp_data <- sp_sim$data

test_that("an empty weights matrix reproduces the non-spatial fit", {
  W0 <- matrix(0, nrow(sp_data), nrow(sp_data))
  fit_sp <- metafrontier(log_y ~ log_x1 + log_x2, data = sp_data,
                         group = "group", spatial = list(W = W0))
  fit <- metafrontier(log_y ~ log_x1 + log_x2, data = sp_data,
                      group = "group")
  expect_equal(fit_sp$group_coef, fit$group_coef)
  expect_equal(fit_sp$tgr, fit$tgr)
  expect_equal(fit_sp$te_group, fit$te_group)
  expect_equal(unname(fit_sp$group_models$G1$rho), 0)
})

test_that("SAR frontiers recover rho, beta and the spillover effects", {
  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 400,
                               sigma_v = 0.1, seed = 18,
                               spatial = list(rho = 0.4, k = 5))
  expect_equal(dim(sim$W), c(800, 800))
  # No spillovers between groups in the simulated weights
  expect_equal(sum(sim$W[1:400, 401:800]), 0)
  fit <- metafrontier(log_y ~ log_x1 + log_x2, data = sim$data,
                      group = "group", spatial = list(W = sim$W))
  for (g in c("G1", "G2")) {
    gm <- fit$group_models[[g]]
    expect_equal(gm$rho, 0.4, tolerance = 0.1)
    expect_equal(unname(gm$coefficients[-1]),
                 sim$params$beta_groups[[g]][-1], tolerance = 0.1)
    # Row-standardised W: total effect = beta / (1 - rho)
    expect_equal(gm$effects$total,
                 unname(gm$coefficients[-1]) / (1 - gm$rho),
                 tolerance = 1e-6)
    expect_true(all(gm$effects$se_total > 0))
    expect_true(is.finite(sqrt(diag(solve(-gm$hessian)))["rho"]))
  }
  expect_true(all(fit$tgr > 0 & fit$tgr <= 1))
})

sp2 <- simulate_metafrontier(n_groups = 2, n_per_group = 80, sigma_v = 0.1,
                             seed = 19, spatial = list(rho = 0.3, k = 4))

test_that("dense, sparse and listw weights give the same fit", {
  skip_if_not_installed("Matrix")
  skip_if_not_installed("spdep")
  f <- log_y ~ log_x1 + log_x2
  fit_d <- metafrontier(f, data = sp2$data, group = "group",
                        spatial = list(W = sp2$W))
  fit_m <- metafrontier(f, data = sp2$data, group = "group",
                        spatial = list(W = Matrix::Matrix(sp2$W,
                                                          sparse = TRUE)))
  fit_l <- metafrontier(f, data = sp2$data, group = "group",
                        spatial = list(W = spdep::mat2listw(sp2$W,
                                                            style = "M")))
  expect_equal(fit_m$tgr, fit_d$tgr)
  expect_equal(fit_l$tgr, fit_d$tgr)
  expect_equal(fit_l$group_models$G1$rho, fit_d$group_models$G1$rho)
})

test_that("invalid weights and unsupported combinations fail loudly", {
  f <- log_y ~ log_x1 + log_x2
  d <- sp2$data
  W <- sp2$W
  expect_error(metafrontier(f, data = d, group = "group",
                            spatial = list(W = W[-1, ])), "n x n")
  W_diag <- W
  W_diag[1, 1] <- 1
  expect_error(metafrontier(f, data = d, group = "group",
                            spatial = list(W = W_diag)), "zero diagonal")
  W_na <- W
  W_na[1, 2] <- NA
  expect_error(metafrontier(f, data = d, group = "group",
                            spatial = list(W = W_na)), "missing")
  expect_error(metafrontier(f, data = d, group = "group",
                            spatial = list(W = W), method = "dea"),
               "spatial")
  fit <- metafrontier(f, data = d, group = "group",
                      spatial = list(W = W), meta_type = "stochastic")
  expect_error(boot_tgr(fit, R = 5, progress = FALSE), "spatial")
  expect_error(vcov(fit, correction = "murphy-topel"), "spatial")
  expect_error(poolability_test(fit), "spatial")
})

test_that("rows dropped for missing values are dropped from W too", {
  d_na <- sp2$data
  d_na$log_x1[7] <- NA
  f <- log_y ~ log_x1 + log_x2
  fit_na <- metafrontier(f, data = d_na, group = "group",
                         spatial = list(W = sp2$W))
  fit_cl <- metafrontier(f, data = d_na[-7, ], group = "group",
                         spatial = list(W = sp2$W[-7, -7]))
  expect_equal(fit_na$tgr, fit_cl$tgr)
  expect_equal(fit_na$group_models$G1$rho, fit_cl$group_models$G1$rho)
})
