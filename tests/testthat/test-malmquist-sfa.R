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

test_that("SFA Malmquist recovers the direction of a known frontier shift", {
  # Both group frontiers shift up by 0.2 in logs between periods 1 and 2,
  # so technical change is exp(0.2) ~ 1.22 and there is no efficiency trend
  set.seed(7)
  n <- 80
  shift_data <- do.call(rbind, lapply(c("A", "B"), function(g) {
    do.call(rbind, lapply(1:2, function(t) {
      lx1 <- log(runif(n, 1, 10))
      lx2 <- log(runif(n, 1, 10))
      log_y <- 1 - 0.3 * (g == "B") + 0.5 * lx1 + 0.3 * lx2 +
        0.2 * (t - 1) - abs(rnorm(n, 0, 0.2)) + rnorm(n, 0, 0.05)
      data.frame(id = paste0(g, seq_len(n)), group = g, time = t,
                 log_x1 = lx1, log_x2 = lx2, log_y = log_y)
    }))
  }))

  malm <- suppressMessages(malmquist_meta(
    log_y ~ log_x1 + log_x2, data = shift_data, group = "group",
    time = "time", id = "id", method = "sfa"
  ))

  expect_equal(mean(malm$group_malmquist$TC_group), exp(0.2),
               tolerance = 0.05)
  expect_equal(mean(malm$meta_malmquist$TC_meta), exp(0.2),
               tolerance = 0.05)
  # Meta efficiency change is TE*_t / TE*_s = TEC x TGC
  expect_equal(malm$meta_malmquist$EC_meta,
               malm$malmquist$TEC * malm$malmquist$TGC)
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
