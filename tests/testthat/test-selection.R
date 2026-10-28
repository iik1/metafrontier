# Sample-selection metafrontier (Greene, 2010): probit selection
# d* = z'gamma + w, frontier y = x'beta + v - u observed when d = 1,
# corr(w, v) = rho

make_sel_data <- function(n = 1000, rho = 0.5, seed = 3) {
  set.seed(seed)
  do.call(rbind, lapply(c("A", "B"), function(g) {
    z1 <- rnorm(n)
    x1 <- runif(n, 0, 3)
    x2 <- runif(n, 0, 3)
    w <- rnorm(n)
    v <- 0.2 * (rho * w + sqrt(1 - rho^2) * rnorm(n))
    u <- abs(rnorm(n, 0, 0.3))
    d <- as.numeric(0.3 + 0.8 * z1 + w > 0)
    y <- 1 - 0.3 * (g == "B") + 0.5 * x1 + 0.3 * x2 + v - u
    y[d == 0] <- NA
    data.frame(group = g, d = d, z1 = z1, x1 = x1, x2 = x2, y = y)
  }))
}
sel_data <- make_sel_data()

test_that("selection group frontiers match sfaR::sfaselectioncross", {
  skip_if_not_installed("sfaR")
  fit <- metafrontier(y ~ x1 + x2, data = sel_data, group = "group",
                      selection = d ~ z1)
  for (g in c("A", "B")) {
    dg <- sel_data[sel_data$group == g, ]
    dg$y[dg$d == 0] <- 0   # sfaR needs a value for unselected rows
    ref <- sfaR::sfaselectioncross(selectionF = d ~ z1,
                                   frontierF = y ~ x1 + x2, data = dg,
                                   udist = "hnormal", printInfo = FALSE)
    cf <- coef(ref)
    gm <- fit$group_models[[g]]
    expect_equal(unname(gm$coefficients), unname(cf[1:3]), tolerance = 1e-2)
    expect_equal(unname(gm$sigma_u), exp(cf[["Zu_(Intercept)"]] / 2),
                 tolerance = 1e-2)
    expect_equal(unname(gm$sigma_v), exp(cf[["Zv_(Intercept)"]] / 2),
                 tolerance = 1e-2)
    expect_equal(unname(gm$rho), cf[["rho"]], tolerance = 2e-2)
    # TE = E[exp(-u) | eps, d = 1] by direct numerical integration
    # (sfaR's efficiencies() misaligns the probit index in 1.0.1)
    zg <- as.numeric(cbind(1, dg$z1[dg$d == 1]) %*% gm$gamma)
    te_int <- function(e, z) {
      f <- function(u) dnorm((e + u) / gm$sigma_v) *
        pnorm((z + gm$rho * (e + u) / gm$sigma_v) / sqrt(1 - gm$rho^2)) *
        dnorm(u / gm$sigma_u)
      integrate(function(u) exp(-u) * f(u), 0, Inf)$value /
        integrate(f, 0, Inf)$value
    }
    i <- 1:20
    expect_equal(gm$efficiency[i], mapply(te_int, gm$residuals[i], zg[i]),
                 tolerance = 5e-3)
  }
  n_sel <- sum(sel_data$d == 1)
  expect_equal(nobs(fit), n_sel)
  expect_length(fit$tgr, n_sel)
  expect_true(all(fit$te_group > 0 & fit$te_group <= 1))
})

test_that("with rho = 0 the selection fit is close to the uncorrected fit", {
  d0 <- make_sel_data(rho = 0, seed = 9)
  fit_s <- metafrontier(y ~ x1 + x2, data = d0, group = "group",
                        selection = d ~ z1)
  fit_u <- metafrontier(y ~ x1 + x2, data = d0[d0$d == 1, ],
                        group = "group")
  expect_equal(fit_s$group_coef$A, fit_u$group_coef$A, tolerance = 0.05)
  expect_equal(mean(fit_s$tgr[fit_s$group_vec == "B"]),
               mean(fit_u$tgr[fit_u$group_vec == "B"]), tolerance = 0.02)
  expect_lt(abs(fit_s$group_models$A$rho), 0.3)
})

test_that("selection keeps unselected rows with a missing outcome and drops incomplete ones", {
  d_na <- sel_data
  d_na$z1[3] <- NA            # missing selection regressor: dropped
  d_na$x1[which(d_na$d == 1)[1]] <- NA   # selected, missing input: dropped
  drop <- c(3, which(d_na$d == 1)[1])
  fit_na <- metafrontier(y ~ x1 + x2, data = d_na, group = "group",
                         selection = d ~ z1)
  fit_cl <- metafrontier(y ~ x1 + x2, data = d_na[-drop, ],
                         group = "group", selection = d ~ z1)
  expect_equal(fit_na$tgr, fit_cl$tgr)
  expect_equal(fit_na$group_models$A$n_total, sum(d_na$group[-drop] == "A"))
  expect_error(
    metafrontier(y ~ x1 + x2, data = transform(sel_data, d = d * 2),
                 group = "group", selection = d ~ z1),
    "0/1"
  )
})

test_that("cost selection frontiers mirror the production fit", {
  fit_p <- metafrontier(y ~ x1 + x2, data = sel_data, group = "group",
                        selection = d ~ z1)
  fit_c <- metafrontier(I(-y) ~ x1 + x2, data = sel_data, group = "group",
                        selection = d ~ z1, frontier = "cost")
  gp <- fit_p$group_models$A
  gc <- fit_c$group_models$A
  expect_equal(unname(gc$coefficients), -unname(gp$coefficients),
               tolerance = 1e-6)
  expect_equal(gc$rho, -gp$rho, tolerance = 1e-6)
  expect_equal(gc$efficiency, gp$efficiency, tolerance = 1e-6)
})

test_that("unsupported selection combinations fail loudly", {
  expect_error(metafrontier(y ~ x1 + x2, data = sel_data, group = "group",
                            selection = d ~ z1, method = "dea"), "selection")
  expect_error(metafrontier(y ~ x1 + x2, data = sel_data, group = "group",
                            selection = d ~ z1, engine = "sfaR"),
               "selection")
  expect_error(metafrontier(y ~ x1 + x2 | z1, data = sel_data,
                            group = "group", selection = d ~ z1),
               "inefficiency determinants")
  fit <- metafrontier(y ~ x1 + x2, data = sel_data, group = "group",
                      selection = d ~ z1, meta_type = "stochastic")
  expect_error(boot_tgr(fit, R = 5, progress = FALSE), "selection")
  expect_error(vcov(fit, correction = "murphy-topel"), "selection")
  expect_error(poolability_test(fit), "selection")
})

test_that("simulate_metafrontier(selection = ) generates recoverable selection data", {
  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 1500,
                               sigma_v = 0.2, seed = 12,
                               selection = list(gamma = c(0.3, 0.8),
                                                rho = 0.6))
  d <- sim$data
  expect_true(all(c("z1", "d") %in% names(d)))
  expect_true(all(is.na(d$log_y[d$d == 0])))
  expect_true(all(!is.na(d$log_y[d$d == 1])))
  expect_equal(sim$params$selection$rho, 0.6)

  fit <- metafrontier(log_y ~ log_x1 + log_x2, data = d, group = "group",
                      selection = d ~ z1)
  rho_hat <- vapply(fit$group_models, function(m) m$rho, numeric(1))
  expect_true(all(rho_hat > 0.2))
  expect_equal(unname(fit$group_coef$G1), sim$params$beta_groups$G1,
               tolerance = 0.1)
  expect_output(print(summary(fit)), "Selection")
  expect_true(all(check_convergence(fit)$converged[1:2]))
})
