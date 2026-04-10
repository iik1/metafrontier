test_that("print method produces output", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data, group = "group")
  expect_output(print(fit), "Metafrontier Model")
  expect_output(print(fit), "Method:")
  expect_output(print(fit), "Mean TGR")
})

test_that("summary returns correct class", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data, group = "group")
  s <- summary(fit)
  expect_s3_class(s, "summary.metafrontier")
  expect_output(print(s), "Metafrontier Model Summary")
  expect_output(print(s), "Efficiency Decomposition")
})

test_that("coef extracts meta and group coefficients", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data, group = "group")

  # Meta coefficients
  mc <- coef(fit, which = "meta")
  expect_length(mc, 3)
  expect_named(mc)
  expect_true(is.numeric(mc))

  # Group coefficients
  gc <- coef(fit, which = "group")
  expect_type(gc, "list")
  expect_length(gc, 2)
})

test_that("vcov returns matrix for stochastic metafrontier", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data, group = "group",
                      meta_type = "stochastic")

  v <- vcov(fit)
  expect_true(is.matrix(v))
  expect_equal(nrow(v), ncol(v))
  expect_equal(nrow(v), length(coef(fit)))

  # Symmetric
  expect_equal(v, t(v))
})

test_that("vcov warns for deterministic metafrontier", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data, group = "group",
                      meta_type = "deterministic")
  expect_warning(vcov(fit), "not available")
})

test_that("logLik returns proper class", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data, group = "group")

  ll <- logLik(fit)
  expect_s3_class(ll, "logLik")
  expect_true(is.finite(as.numeric(ll)))
  expect_true(!is.null(attr(ll, "df")))
  expect_true(!is.null(attr(ll, "nobs")))

  # AIC and BIC should work automatically
  expect_true(is.finite(AIC(fit)))
  expect_true(is.finite(BIC(fit)))
})

test_that("nobs returns total observations", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data, group = "group")
  expect_equal(nobs(fit), nrow(test_data))
})

test_that("efficiencies extracts correct type", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data, group = "group")

  te <- efficiencies(fit, type = "group")
  expect_length(te, nrow(test_data))

  te_star <- efficiencies(fit, type = "meta")
  expect_length(te_star, nrow(test_data))

  tgr <- efficiencies(fit, type = "tgr")
  expect_length(tgr, nrow(test_data))

  # Decomposition
  expect_equal(te_star, te * tgr, tolerance = 1e-10)
})

test_that("technology_gap_ratio works by group", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data, group = "group")

  tgr_list <- technology_gap_ratio(fit, by_group = TRUE)
  expect_type(tgr_list, "list")
  expect_length(tgr_list, 2)
  expect_named(tgr_list, fit$groups)

  tgr_flat <- technology_gap_ratio(fit, by_group = FALSE)
  expect_length(tgr_flat, nrow(test_data))
})

test_that("tgr_summary returns data frame", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data, group = "group")

  tab <- tgr_summary(fit)
  expect_s3_class(tab, "data.frame")
  expect_equal(nrow(tab), 2)
  expect_true("Mean" %in% names(tab))
  expect_true("SD" %in% names(tab))
})

test_that("plot does not error", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data, group = "group")

  expect_no_error(plot(fit, which = "tgr"))
  expect_no_error(plot(fit, which = "efficiency"))
  expect_no_error(plot(fit, which = "decomposition"))
})

test_that("fitted and residuals work", {
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = test_data, group = "group")

  f <- fitted(fit)
  expect_length(f, nrow(test_data))

  r <- residuals(fit)
  expect_length(r, nrow(test_data))
})
