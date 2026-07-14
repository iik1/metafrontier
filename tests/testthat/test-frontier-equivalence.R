# Equivalence check against frontier::sfa()
# The package's custom MLE should match frontier::sfa() slope and intercept
# estimates to within 1e-4 on the same simulated data.

skip_if_not_installed("frontier")

test_that("group SFA coefficients match frontier::sfa() to 1e-4", {
  sim <- simulate_metafrontier(n_groups = 3, n_per_group = 150,
                               tech_gap = c(0, 0.3, 0.5), seed = 42)
  dat <- sim$data

  fit_mf <- metafrontier(log_y ~ log_x1 + log_x2, data = dat,
                         group = "group", method = "sfa",
                         meta_type = "deterministic")

  for (g in c("G2", "G3")) {
    dat_g <- dat[dat$group == g, , drop = FALSE]
    f_g <- frontier::sfa(log_y ~ log_x1 + log_x2, data = dat_g)

    b_mf <- coef(fit_mf$group_models[[g]])
    b_fr <- coef(f_g)[names(b_mf)]

    expect_lt(max(abs(b_mf - b_fr)), 1e-4,
              label = paste("max |coef diff| for group", g))
  }
})
