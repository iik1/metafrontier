# Tests for SFA-based Malmquist index (P4)

# Create panel data for Malmquist tests
set.seed(42)
malm_panels <- lapply(1:3, function(t) {
  sim <- simulate_metafrontier(
    n_groups = 2, n_per_group = 40,
    tech_gap = c(0, 0.3 + 0.05 * t),
    sigma_u = c(0.2, 0.3),
    seed = 42 + t
  )
  sim$data$time <- t
  sim$data$id <- seq_len(nrow(sim$data))
  sim$data
})
malm_data <- do.call(rbind, malm_panels)

test_that("SFA Malmquist produces valid output", {
  malm <- malmquist_meta(
    log_y ~ log_x1 + log_x2,
    data = malm_data,
    group = "group",
    time = "time",
    method = "sfa"
  )

  expect_s3_class(malm, "malmquist_meta")
  expect_true(nrow(malm$malmquist) > 0)
  expect_true(all(c("MPI", "TEC", "TGC", "TC") %in% names(malm$malmquist)))
})

test_that("SFA Malmquist TEC x TGC x TC* approximates MPI", {
  malm <- malmquist_meta(
    log_y ~ log_x1 + log_x2,
    data = malm_data,
    group = "group",
    time = "time",
    method = "sfa"
  )

  m <- malm$malmquist
  valid <- complete.cases(m[, c("MPI", "TEC", "TGC", "TC")])
  if (sum(valid) > 0) {
    product <- m$TEC[valid] * m$TGC[valid] * m$TC[valid]
    expect_equal(product, m$MPI[valid], tolerance = 0.1)
  }
})

test_that("DEA Malmquist still works with method arg", {
  malm_dea <- malmquist_meta(
    log_y ~ log_x1 + log_x2,
    data = malm_data,
    group = "group",
    time = "time",
    method = "dea"
  )

  expect_s3_class(malm_dea, "malmquist_meta")
  expect_true(nrow(malm_dea$malmquist) > 0)
})

test_that("SFA Malmquist TGR values are positive", {
  malm <- malmquist_meta(
    log_y ~ log_x1 + log_x2,
    data = malm_data,
    group = "group",
    time = "time",
    method = "sfa"
  )

  valid <- malm$tgr$TGR_from[is.finite(malm$tgr$TGR_from)]
  expect_true(all(valid > 0))
})
