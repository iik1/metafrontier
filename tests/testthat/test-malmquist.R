# Test Malmquist metafrontier index

# Helper: create small panel dataset
make_panel <- function(n_groups = 2, n_per_group = 20,
                       n_periods = 3, seed = 42) {
  panels <- lapply(seq_len(n_periods), function(t) {
    sim <- simulate_metafrontier(
      n_groups = n_groups,
      n_per_group = n_per_group,
      tech_gap = c(0, 0.3 + 0.05 * t),
      sigma_u = c(0.2, 0.3),
      seed = seed + t
    )
    sim$data$time <- t
    sim$data$id <- seq_len(nrow(sim$data))
    sim$data
  })
  do.call(rbind, panels)
}


test_that("malmquist_meta returns correct structure", {
  pd <- make_panel(seed = 1)
  malm <- malmquist_meta(log_y ~ log_x1 + log_x2,
                         data = pd, group = "group",
                         time = "time")

  expect_s3_class(malm, "malmquist_meta")
  expect_true(is.data.frame(malm$malmquist))
  expect_true(is.data.frame(malm$group_malmquist))
  expect_true(is.data.frame(malm$meta_malmquist))
  expect_true(is.data.frame(malm$tgr))

  # Check required columns
  expect_true(all(c("id", "group", "period_from", "period_to",
                     "MPI", "TEC", "TGC", "TC") %in%
                    names(malm$malmquist)))
  expect_true(all(c("MPI_group", "EC_group", "TC_group") %in%
                    names(malm$group_malmquist)))
  expect_true(all(c("MPI_meta", "EC_meta", "TC_meta") %in%
                    names(malm$meta_malmquist)))
  expect_true(all(c("TGR_from", "TGR_to", "TGC") %in%
                    names(malm$tgr)))
})


test_that("malmquist_meta covers all period pairs", {
  pd <- make_panel(n_periods = 4, seed = 2)
  malm <- malmquist_meta(log_y ~ log_x1 + log_x2,
                         data = pd, group = "group",
                         time = "time")

  # 3 consecutive period transitions
  expect_equal(sort(unique(malm$malmquist$period_from)),
               c(1, 2, 3))
  expect_equal(sort(unique(malm$malmquist$period_to)),
               c(2, 3, 4))
})


test_that("MPI = TEC x TGC x TC decomposition holds", {
  pd <- make_panel(seed = 3)
  malm <- malmquist_meta(log_y ~ log_x1 + log_x2,
                         data = pd, group = "group",
                         time = "time")

  m <- malm$malmquist
  complete <- complete.cases(m[, c("MPI", "TEC", "TGC", "TC")])
  sub <- m[complete, ]

  if (nrow(sub) > 0) {
    recomposed <- sub$TEC * sub$TGC * sub$TC
    expect_equal(sub$MPI, recomposed, tolerance = 1e-6)
  }
})


test_that("TGC = TGR_to / TGR_from", {
  pd <- make_panel(seed = 4)
  malm <- malmquist_meta(log_y ~ log_x1 + log_x2,
                         data = pd, group = "group",
                         time = "time")

  tgr_df <- malm$tgr
  complete <- complete.cases(tgr_df)
  sub <- tgr_df[complete, ]

  if (nrow(sub) > 0) {
    expect_equal(sub$TGC, sub$TGR_to / sub$TGR_from,
                 tolerance = 1e-6)
  }
})


test_that("TGR values are in (0, 1]", {
  pd <- make_panel(seed = 5)
  malm <- malmquist_meta(log_y ~ log_x1 + log_x2,
                         data = pd, group = "group",
                         time = "time")

  tgr_df <- malm$tgr
  expect_true(all(tgr_df$TGR_from > 0, na.rm = TRUE))
  expect_true(all(tgr_df$TGR_from <= 1 + 1e-6, na.rm = TRUE))
  expect_true(all(tgr_df$TGR_to > 0, na.rm = TRUE))
  expect_true(all(tgr_df$TGR_to <= 1 + 1e-6, na.rm = TRUE))
})


test_that("malmquist_meta works with VRS", {
  pd <- make_panel(seed = 6)
  expect_no_error(
    malmquist_meta(log_y ~ log_x1 + log_x2,
                   data = pd, group = "group",
                   time = "time", rts = "vrs")
  )
})


test_that("malmquist_meta print and summary work", {
  pd <- make_panel(seed = 7)
  malm <- malmquist_meta(log_y ~ log_x1 + log_x2,
                         data = pd, group = "group",
                         time = "time")

  expect_output(print(malm), "Metafrontier Malmquist")
  s <- summary(malm)
  expect_s3_class(s, "summary.malmquist_meta")
  expect_output(print(s), "Three-Way Decomposition")
})


test_that("malmquist_meta errors on missing inputs", {
  pd <- make_panel(seed = 8)

  expect_error(
    malmquist_meta(log_y ~ log_x1 + log_x2,
                   data = pd, group = "group"),
    "required"
  )
  expect_error(
    malmquist_meta(log_y ~ log_x1 + log_x2,
                   data = pd, time = "time"),
    "required"
  )
})


test_that("malmquist_meta errors with single group or period", {
  pd <- make_panel(seed = 9)
  pd_one_group <- droplevels(pd[pd$group == "G1", ])

  expect_error(
    malmquist_meta(log_y ~ log_x1 + log_x2,
                   data = pd_one_group, group = "group",
                   time = "time"),
    "2 groups"
  )

  pd_one_time <- pd[pd$time == 1, ]
  expect_error(
    malmquist_meta(log_y ~ log_x1 + log_x2,
                   data = pd_one_time, group = "group",
                   time = "time"),
    "2 time periods"
  )
})


test_that("malmquist_meta with 3 groups and single input", {
  panels <- lapply(1:2, function(t) {
    sim <- simulate_metafrontier(
      n_groups = 3, n_per_group = 15,
      n_inputs = 1,
      tech_gap = c(0, 0.2, 0.4),
      seed = 100 + t
    )
    sim$data$time <- t
    sim$data$id <- seq_len(nrow(sim$data))
    sim$data
  })
  pd <- do.call(rbind, panels)

  malm <- malmquist_meta(log_y ~ log_x1,
                         data = pd, group = "group",
                         time = "time")

  expect_s3_class(malm, "malmquist_meta")
  expect_length(malm$groups, 3)
  expect_true(nrow(malm$malmquist) > 0)
})
