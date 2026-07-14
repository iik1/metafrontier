# Tests for firm matching via `id`, infeasibility accounting, and
# method visibility in malmquist_meta().

# Deterministic positive-level panel: 2 groups x n_per_group firms x
# n_periods periods, suitable for DEA on levels.
make_id_panel <- function(n_per_group = 6, n_periods = 3) {
  df <- expand.grid(firm = seq_len(2 * n_per_group),
                    time = seq_len(n_periods))
  df$group <- ifelse(df$firm <= n_per_group, "A", "B")
  df$x1 <- 1 + 0.1 * (df$firm %% 5) + 0.05 * df$time
  df$x2 <- 1.5 + 0.08 * (df$firm %% 7)
  eff <- 0.6 + 0.4 * ((df$firm * 7) %% 10) / 10
  df$y <- eff * df$x1^0.4 * df$x2^0.3 * (1 + 0.1 * df$time)
  df
}

test_that("id gives order-invariant results on a balanced panel", {
  pd <- make_id_panel()

  base <- suppressMessages(
    malmquist_meta(y ~ x1 + x2, data = pd,
                   group = "group", time = "time")
  )
  set.seed(1)
  pd_scr <- pd[sample(nrow(pd)), ]
  scr <- malmquist_meta(y ~ x1 + x2, data = pd_scr,
                        group = "group", time = "time", id = "firm")

  # Positional ids in the baseline are within-group positions; map
  # them back to firm numbers (group A = 1..6, group B = 7..12).
  base_m <- base$malmquist
  base_m$firm <- base_m$id + ifelse(base_m$group == "B", 6L, 0L)
  scr_m <- scr$malmquist

  base_m <- base_m[order(base_m$period_from, base_m$group, base_m$firm), ]
  scr_m <- scr_m[order(scr_m$period_from, scr_m$group, scr_m$id), ]

  expect_equal(scr_m$id, base_m$firm)
  expect_equal(scr_m$MPI, base_m$MPI, tolerance = 1e-8)
  expect_equal(scr_m$TEC, base_m$TEC, tolerance = 1e-8)
  expect_equal(scr_m$TGC, base_m$TGC, tolerance = 1e-8)
  expect_equal(scr_m$TC,  base_m$TC,  tolerance = 1e-8)

  # Method visibility
  expect_identical(scr$method, "dea")
  expect_output(print(scr), "Method")
})

test_that("unbalanced panels drop unmatched firms with a consolidated warning", {
  pd <- make_id_panel()

  # Make firms 3 (group A) and 9 (group B) strictly dominated (same
  # inputs as firms 2 and 8, much lower output) so that dropping them
  # leaves every DEA reference technology unchanged.
  for (tt in unique(pd$time)) {
    r3 <- pd$firm == 3 & pd$time == tt
    r2 <- pd$firm == 2 & pd$time == tt
    pd[r3, c("x1", "x2")] <- pd[r2, c("x1", "x2")]
    pd$y[r3] <- 0.4 * pd$y[r2]
    r9 <- pd$firm == 9 & pd$time == tt
    r8 <- pd$firm == 8 & pd$time == tt
    pd[r9, c("x1", "x2")] <- pd[r8, c("x1", "x2")]
    pd$y[r9] <- 0.4 * pd$y[r8]
  }

  bal <- malmquist_meta(y ~ x1 + x2, data = pd,
                        group = "group", time = "time", id = "firm")

  pd_unb <- pd[!(pd$firm == 3 & pd$time == 2) &
               !(pd$firm == 9 & pd$time == 1), ]
  # 3 dropped in total: firm 3 + firm 9 in pair 1 -> 2, firm 3 in 2 -> 3
  expect_warning(
    unb <- malmquist_meta(y ~ x1 + x2, data = pd_unb,
                          group = "group", time = "time", id = "firm"),
    "3 observation"
  )

  # The right firms are dropped
  expect_false(any(unb$malmquist$id == 3))
  expect_false(any(unb$malmquist$id == 9 & unb$malmquist$period_from == 1))
  expect_true(any(unb$malmquist$id == 9 & unb$malmquist$period_from == 2))

  # Retained firms match the balanced results exactly
  key_b <- paste(bal$malmquist$id, bal$malmquist$period_from)
  key_u <- paste(unb$malmquist$id, unb$malmquist$period_from)
  expect_true(all(key_u %in% key_b))
  pos <- match(key_u, key_b)
  expect_equal(unb$malmquist$MPI, bal$malmquist$MPI[pos], tolerance = 1e-8)
  expect_equal(unb$malmquist$TC,  bal$malmquist$TC[pos],  tolerance = 1e-8)
})

test_that("duplicated (id, period) combinations within a group error", {
  pd <- make_id_panel()
  pd_dup <- rbind(pd, pd[pd$firm == 1 & pd$time == 1, ])
  expect_error(
    malmquist_meta(y ~ x1 + x2, data = pd_dup,
                   group = "group", time = "time", id = "firm"),
    "Duplicated"
  )
  expect_error(
    malmquist_meta(y ~ x1 + x2, data = pd_dup,
                   group = "group", time = "time", id = "firm"),
    "period 1"
  )
})

test_that("id = NULL messages, and unequal group sizes warn", {
  pd <- make_id_panel()

  expect_message(
    malmquist_meta(y ~ x1 + x2, data = pd,
                   group = "group", time = "time"),
    "row position"
  )

  pd_unb <- pd[!(pd$firm == 3 & pd$time == 2), ]
  expect_warning(
    suppressMessages(
      malmquist_meta(y ~ x1 + x2, data = pd_unb,
                     group = "group", time = "time")
    ),
    "positional matching"
  )
})

test_that("infeasible cross-period programs are counted and warned under vrs", {
  pd <- make_id_panel(n_per_group = 5, n_periods = 2)

  # Anchor firms (one per group) define the CRS frontier in both
  # periods, so all cross-period CRS programs stay feasible.
  pd[pd$firm == 1, c("x1", "x2", "y")] <- list(1, 1, 2)
  pd[pd$firm == 6, c("x1", "x2", "y")] <- list(1, 1, 2)
  # Tiny firm whose period-1 inputs undercut every period-2
  # observation: its period-1 cross-period VRS programs (group and
  # meta) are genuinely infeasible.
  pd[pd$firm == 2 & pd$time == 1, c("x1", "x2", "y")] <-
    list(0.01, 0.01, 0.001)
  pd[pd$firm == 2 & pd$time == 2, c("x1", "x2", "y")] <-
    list(0.05, 0.05, 0.002)

  expect_warning(
    malm_vrs <- malmquist_meta(y ~ x1 + x2, data = pd,
                               group = "group", time = "time",
                               id = "firm", rts = "vrs"),
    "cross-period DEA programs were infeasible"
  )
  expect_equal(malm_vrs$n_infeasible, 2L)
  expect_true(any(is.na(malm_vrs$malmquist$TC)))
  expect_true(is.data.frame(malm_vrs$infeasible_by_period))
  expect_equal(sum(malm_vrs$infeasible_by_period$n_infeasible),
               malm_vrs$n_infeasible)
  expect_output(print(malm_vrs), "Infeasible")
  expect_output(print(summary(malm_vrs)), "Infeasible")

  # CRS on the same data has no infeasible programs and no NAs
  expect_no_warning(
    malm_crs <- malmquist_meta(y ~ x1 + x2, data = pd,
                               group = "group", time = "time",
                               id = "firm", rts = "crs")
  )
  expect_identical(malm_crs$n_infeasible, 0L)
  expect_false(any(is.na(malm_crs$malmquist$TC)))
  expect_false(any(is.na(malm_crs$malmquist$MPI)))
})

test_that("sfa path notes the approximation and carries id through", {
  set.seed(7)
  n <- 12
  df <- expand.grid(firm = seq_len(2 * n), time = 1:2)
  df$group <- ifelse(df$firm <= n, "A", "B")
  df$log_x1 <- rnorm(nrow(df))
  df$log_x2 <- rnorm(nrow(df))
  df$log_y <- 1 + 0.4 * df$log_x1 + 0.3 * df$log_x2 +
    0.2 * (df$group == "B") +
    rnorm(nrow(df), sd = 0.1) - abs(rnorm(nrow(df), sd = 0.3))

  suppressWarnings(
    expect_message(
      malm <- malmquist_meta(log_y ~ log_x1 + log_x2, data = df,
                             group = "group", time = "time",
                             id = "firm", method = "sfa"),
      "pointwise maximum"
    )
  )
  expect_identical(malm$method, "sfa")
  expect_true(all(malm$malmquist$id %in% df$firm))
  expect_output(print(malm), "pointwise-maximum")
})
