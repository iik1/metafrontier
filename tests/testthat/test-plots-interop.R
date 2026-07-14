# Tests for revised plots (Figures 4 and 5) and model interoperability

# ---------- Plot revisions ----------

test_that("autoplot.boot_tgr distribution shows CI bounds as vlines", {
  skip_if_not_installed("ggplot2")

  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 40, seed = 42)
  fit <- metafrontier(log_y ~ log_x1 + log_x2, data = sim$data,
                      group = "group")
  boot <- boot_tgr(fit, R = 20, seed = 1, progress = FALSE)

  p <- ggplot2::autoplot(boot, which = "distribution")
  expect_s3_class(p, "gg")

  is_vline <- vapply(p$layers, function(l) inherits(l$geom, "GeomVline"),
                     logical(1))
  expect_true(any(is_vline))

  # One dashed line per group and per CI bound
  vline_data <- ggplot2::layer_data(p, which(is_vline)[1])
  expect_equal(nrow(vline_data), 2L * length(boot$groups))
})

test_that("autoplot.malmquist_meta mpi_trend plots at end periods", {
  skip_if_not_installed("ggplot2")

  psim <- simulate_panel_metafrontier(n_groups = 2, n_firms_per_group = 15,
                                      n_periods = 3, seed = 7)
  malm <- suppressWarnings(suppressMessages(
    malmquist_meta(log_y ~ log_x1 + log_x2, data = psim$data,
                   group = "group", time = "year", method = "dea")
  ))

  p <- ggplot2::autoplot(malm, which = "mpi_trend")
  expect_s3_class(p, "gg")

  built <- ggplot2::ggplot_build(p)
  xvals <- sort(unique(unlist(lapply(built$data[1:2],
                                     function(d) d$x))))
  expect_equal(as.numeric(xvals), 2:3)
  expect_match(p$labels$x, "change from previous period")
})

# ---------- as_metafrontier_model interoperability ----------

test_that("as_metafrontier_model is idempotent and feeds metafrontier()", {
  mk_group <- function(dat) {
    X <- model.matrix(~ log_x1 + log_x2, dat)
    ols <- lm.fit(X, dat$log_y)
    res <- ols$residuals
    list(coefficients = ols$coefficients,
         efficiency = exp(res - max(res)),
         X = X, y = dat$log_y)
  }
  grps <- levels(test_data$group)
  conv1 <- as_metafrontier_model(mk_group(test_data[test_data$group ==
                                                      grps[1], ]))
  conv2 <- as_metafrontier_model(mk_group(test_data[test_data$group ==
                                                      grps[2], ]))

  expect_s3_class(conv1, "metafrontier_model")
  # Converting twice is a no-op
  expect_identical(as_metafrontier_model(conv1), conv1)

  fit <- metafrontier(models = list(A = conv1, B = conv2))
  expect_s3_class(fit, "metafrontier")
  expect_equal(fit$groups, c("A", "B"))
})

test_that("frontier::sfa objects pass directly through metafrontier(models=)", {
  skip_if_not_installed("frontier")

  grps <- levels(test_data$group)
  fits <- lapply(grps, function(g) {
    suppressWarnings(
      frontier::sfa(log_y ~ log_x1 + log_x2,
                    data = test_data[test_data$group == g, ])
    )
  })
  names(fits) <- grps

  expect_s3_class(fits[[1]], "frontier")

  fit <- metafrontier(models = fits)
  expect_s3_class(fit, "metafrontier")
  expect_equal(unname(fit$nobs["total"]), nrow(test_data))
  expect_true(all(fit$tgr > 0 & fit$tgr <= 1))
})

test_that("Farrell objects warn on conversion and error in metafrontier()", {
  skip_if_not_installed("Benchmarking")

  X <- exp(as.matrix(test_data[, c("log_x1", "log_x2")]))
  y <- matrix(exp(test_data$log_y), ncol = 1)
  d <- Benchmarking::dea(X, y, RTS = "vrs")

  expect_warning(as_metafrontier_model(d), "cannot be used")

  expect_error(
    suppressWarnings(metafrontier(models = list(A = d, B = d))),
    "does not provide"
  )
})

test_that("unsupported classes get an informative error", {
  expect_error(
    as_metafrontier_model(structure(1, class = "no_such_model")),
    "sfacross.*frontier.*Farrell"
  )
})
