# Tests for convergence diagnostics and extended coef/vcov/summary methods

sim_cs <- simulate_metafrontier(n_groups = 2, n_per_group = 40, seed = 42)
fit_det <- metafrontier(log_y ~ log_x1 + log_x2, data = sim_cs$data,
                        group = "group")
fit_stoch <- metafrontier(log_y ~ log_x1 + log_x2, data = sim_cs$data,
                          group = "group", meta_type = "stochastic")

sim_panel <- simulate_panel_metafrontier(n_groups = 2,
                                         n_firms_per_group = 20,
                                         n_periods = 4,
                                         eta = 0.05, seed = 123)
fit_panel <- metafrontier(log_y ~ log_x1 + log_x2, data = sim_panel$data,
                          group = "group",
                          panel = list(id = "firm", time = "year"),
                          panel_dist = "bc92")


test_that("check_convergence reports all stages for a deterministic SFA fit", {
  conv <- check_convergence(fit_det)

  expect_s3_class(conv, "metafrontier_convergence")
  expect_s3_class(conv, "data.frame")
  expect_named(conv, c("stage", "method", "code", "converged", "note"))
  expect_equal(nrow(conv), length(fit_det$groups) + 1L)
  expect_equal(conv$stage,
               c(paste0("group: ", fit_det$groups), "metafrontier"))
  expect_true(all(conv$converged))
  expect_true(all(conv$method[seq_along(fit_det$groups)] == "MLE"))
  expect_true(conv$method[nrow(conv)] %in% c("LP", "QP", "QP (barrier)"))
  expect_true(is.integer(conv$code))
})


test_that("check_convergence on a stochastic fit reports MLE meta stage", {
  conv <- check_convergence(fit_stoch)
  expect_equal(conv$method[nrow(conv)], "MLE")
  expect_true(all(conv$converged))
})


test_that("check_convergence on a DEA fit uses method DEA and NA codes", {
  fit_dea <- metafrontier(log_y ~ log_x1 + log_x2, data = sim_cs$data,
                          group = "group", method = "dea", rts = "crs")
  conv <- check_convergence(fit_dea)

  expect_true(all(conv$method == "DEA"))
  expect_true(all(is.na(conv$code)))
  expect_true(all(conv$converged))
  expect_true(all(conv$note == ""))
})


test_that("check_convergence flags externally fitted groups", {
  mods <- lapply(fit_det$group_models, function(gm) {
    list(coefficients = gm$coefficients, efficiency = gm$efficiency,
         X = gm$X, y = gm$y)
  })
  fit_ext <- metafrontier(models = mods)
  conv <- check_convergence(fit_ext)

  n_grp <- length(fit_ext$groups)
  expect_true(all(conv$method[seq_len(n_grp)] == "external"))
  expect_true(all(is.na(conv$code[seq_len(n_grp)])))
  expect_true(all(is.na(conv$converged[seq_len(n_grp)])))
  expect_match(conv$note[1], "fitted externally")

  expect_warning(
    coef(fit_ext, which = "group", extraPar = TRUE),
    "externally"
  )
})


test_that("check_convergence.default errors clearly", {
  expect_error(check_convergence(lm(dist ~ speed, data = cars)),
               "not implemented")
})


test_that("coef(..., extraPar = TRUE) exposes eta for a BC92 panel fit", {
  cf <- coef(fit_panel, which = "group", extraPar = TRUE)

  expect_true(is.list(cf))
  expect_named(cf, fit_panel$groups)
  for (g in fit_panel$groups) {
    expect_true(all(c("sigmaV", "sigmaU", "eta") %in% names(cf[[g]])))
    expect_gt(unname(cf[[g]]["sigmaV"]), 0)
    expect_gt(unname(cf[[g]]["sigmaU"]), 0)
  }

  # Backward compatibility: default call unchanged
  expect_identical(coef(fit_panel, which = "group"), fit_panel$group_coef)
})


test_that("coef(..., extraPar = TRUE) appends stage-2 variances (meta)", {
  cf_plain <- coef(fit_stoch)
  cf_extra <- coef(fit_stoch, extraPar = TRUE)

  expect_gt(length(cf_extra), length(cf_plain))
  expect_true(all(c("sigmaV", "sigmaU") %in% names(cf_extra)))
  expect_equal(cf_extra[seq_along(cf_plain)], cf_plain)

  # Deterministic metafrontier has no auxiliary parameters
  expect_message(cf_det <- coef(fit_det, extraPar = TRUE),
                 "No auxiliary parameters")
  expect_identical(cf_det, coef(fit_det))
})


test_that("vcov(..., extraPar = TRUE) returns the full stage-2 matrix", {
  v <- vcov(fit_stoch)
  v_full <- vcov(fit_stoch, extraPar = TRUE)

  expect_gt(nrow(v_full), nrow(v))
  expect_equal(unname(v_full[seq_len(nrow(v)), seq_len(ncol(v))]),
               unname(v))
})


test_that("vcov(..., which = 'group') returns per-group matrices", {
  v_g <- vcov(fit_stoch, which = "group")

  expect_true(is.list(v_g))
  expect_named(v_g, fit_stoch$groups)
  for (g in fit_stoch$groups) {
    gm <- fit_stoch$group_models[[g]]
    expect_true(is.matrix(v_g[[g]]))
    expect_equal(nrow(v_g[[g]]), length(gm$all_params))
    expect_true("log_sigma_v" %in% rownames(v_g[[g]]))
  }
})


test_that("summary of a BC92 fit shows eta with a standard error", {
  s <- summary(fit_panel)

  for (g in fit_panel$groups) {
    tab <- s$group_tables[[g]]
    expect_true("eta" %in% rownames(tab))
    expect_true(is.finite(tab["eta", "Std. Error"]))
  }

  expect_output(print(s), "eta")
})


test_that("print methods report convergence status", {
  expect_output(print(fit_det), "Convergence: OK")
  expect_output(print(summary(fit_det)),
                "All estimation stages converged")
  expect_output(print(check_convergence(fit_det)),
                "Convergence of estimation stages")

  # A non-converged stage triggers warnings in all three displays
  fit_bad <- fit_det
  fit_bad$group_models[[1]]$convergence <- 1L
  expect_output(print(fit_bad), "Convergence: WARNING")
  expect_output(print(summary(fit_bad)), "See \\?check_convergence")
  expect_output(print(check_convergence(fit_bad)),
                "did not converge")
})


test_that("print shows Estimator and Objective lines when present", {
  fit2 <- fit_det
  if (is.null(fit2$estimator)) fit2$estimator <- "bc88"
  if (is.null(fit2$objective)) fit2$objective <- "lp"
  expect_output(print(fit2), "Estimator:\\s+bc88")
  expect_output(print(fit2), "Objective:\\s+lp")
})
