# Rows with a missing model variable or group are dropped once, up
# front, so a fit on data with NA rows equals the fit on the same data
# with those rows removed beforehand

sim_na <- simulate_metafrontier(n_groups = 2, n_per_group = 60,
                                seed = 31)$data
sim_na$log_x1[c(5, 70)] <- NA          # missing input in each group
sim_na$group <- as.character(sim_na$group)
sim_na$group[c(20, 100)] <- NA         # missing group label
na_rows <- c(5, 70, 20, 100)
sim_clean <- sim_na[-na_rows, ]

test_that("SFA metafrontier with NA rows equals the pre-cleaned fit", {
  for (mt in c("deterministic", "stochastic")) {
    fit_na <- metafrontier(log_y ~ log_x1 + log_x2, data = sim_na,
                           group = "group", meta_type = mt)
    fit_cl <- metafrontier(log_y ~ log_x1 + log_x2, data = sim_clean,
                           group = "group", meta_type = mt)
    expect_equal(fit_na$tgr, fit_cl$tgr)
    expect_equal(fit_na$te_group, fit_cl$te_group)
    expect_equal(nobs(fit_na), nobs(fit_cl))
    expect_length(fit_na$group_vec, length(fit_na$tgr))
    expect_no_error(summary(fit_na))
  }
})

test_that("DEA metafrontier with NA rows equals the pre-cleaned fit", {
  f <- exp(log_y) ~ exp(log_x1) + exp(log_x2)
  fit_na <- metafrontier(f, data = sim_na, group = "group",
                         method = "dea")
  fit_cl <- metafrontier(f, data = sim_clean, group = "group",
                         method = "dea")
  expect_equal(fit_na$tgr, fit_cl$tgr)
  expect_equal(fit_na$te_group, fit_cl$te_group)
  expect_true(all(fit_na$tgr > 0 & fit_na$tgr <= 1 + 1e-9))
})

test_that("DEA Malmquist with NA rows equals the pre-cleaned fit", {
  panel <- rbind(transform(sim_clean, time = 1, id = seq_len(nrow(sim_clean))),
                 transform(sim_clean, time = 2, id = seq_len(nrow(sim_clean)),
                           log_y = log_y + 0.1))
  panel_na <- panel
  panel_na$log_x2[3] <- NA
  panel_na$time[200] <- NA
  f <- exp(log_y) ~ exp(log_x1) + exp(log_x2)
  m_na <- suppressWarnings(malmquist_meta(f, data = panel_na,
                                          group = "group", time = "time",
                                          id = "id"))
  m_cl <- suppressWarnings(malmquist_meta(f, data = panel[-c(3, 200), ],
                                          group = "group", time = "time",
                                          id = "id"))
  expect_equal(m_na$malmquist, m_cl$malmquist)
})

test_that("latent class metafrontier drops NA rows instead of failing", {
  lc_na <- suppressWarnings(latent_class_metafrontier(log_y ~ log_x1 + log_x2,
                                     data = sim_na, n_starts = 2,
                                     max_iter = 20, seed = 1))
  lc_cl <- suppressWarnings(latent_class_metafrontier(log_y ~ log_x1 + log_x2,
                                     data = sim_na[-c(5, 70), ],
                                     n_starts = 2, max_iter = 20,
                                     seed = 1))
  expect_equal(lc_na$metafrontier$tgr, lc_cl$metafrontier$tgr)
})
