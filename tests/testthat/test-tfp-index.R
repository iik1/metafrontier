# Fare-Primont and Lowe TFP indexes (O'Donnell 2012) with a metafrontier
# decomposition

make_tfp_data <- function(n_firms = 8, n_periods = 3, seed = 21) {
  set.seed(seed)
  d <- expand.grid(firm = seq_len(n_firms), year = seq_len(n_periods))
  d$region <- ifelse(d$firm <= n_firms / 2, "North", "South")
  d$x1 <- runif(nrow(d), 1, 10)
  d$x2 <- runif(nrow(d), 1, 10)
  shift <- ifelse(d$region == "North", 1, 0.8) * (1 + 0.05 * (d$year - 1))
  eff <- runif(nrow(d), 0.6, 1)
  d$y1 <- shift * eff * d$x1^0.4 * d$x2^0.4 * runif(nrow(d), 0.8, 1.2)
  d$y2 <- shift * eff * d$x1^0.3 * d$x2^0.5 * runif(nrow(d), 0.8, 1.2)
  d$p1 <- runif(nrow(d), 1, 3)
  d$p2 <- runif(nrow(d), 1, 3)
  d$w1 <- runif(nrow(d), 1, 2)
  d$w2 <- runif(nrow(d), 1, 2)
  d
}
tfp_data <- make_tfp_data()

# productivity orders rows by year, then id
align <- function(levels) levels[order(levels$year, levels$firm), ]

test_that("Fare-Primont levels against the metatechnology match productivity::fareprim", {
  skip_if_not_installed("productivity")
  for (orient in c("output", "input")) {
    tfp <- tfp_index(cbind(y1, y2) ~ x1 + x2, data = tfp_data,
                     group = "region", time = "year", id = "firm",
                     index = "fare-primont", orientation = orient,
                     rts = "vrs")
    ref <- suppressWarnings(productivity::fareprim(
      tfp_data, id.var = "firm", time.var = "year",
      x.vars = c("x1", "x2"), y.vars = c("y1", "y2"), rts = "vrs",
      orientation = if (orient == "output") "out" else "in",
      scaled = FALSE))$Levels
    lv <- align(tfp$levels)
    cols <- if (orient == "output") {
      c("OTE", "OSE", "OME", "ROSE", "OSME", "RME")
    } else {
      c("ITE", "ISE", "IME", "RISE", "ISME", "RME")
    }
    expect_equal(lv$TFP, ref$TFP, tolerance = 1e-6)
    expect_equal(lv$TFP_star_meta, ref$MP, tolerance = 1e-6)
    expect_equal(lv$TFPE_meta, ref$TFPE, tolerance = 1e-6)
    for (cl in cols) {
      expect_equal(lv[[paste0(cl, "_meta")]], ref[[cl]], tolerance = 1e-6,
                   label = paste(orient, cl))
    }
  }
})

test_that("Lowe levels and sequential technologies match productivity::lowe", {
  skip_if_not_installed("productivity")
  for (tch in c("contemporaneous", "sequential")) {
    tfp <- tfp_index(cbind(y1, y2) ~ x1 + x2, data = tfp_data,
                     group = "region", time = "year", id = "firm",
                     index = "lowe", output_prices = c("p1", "p2"),
                     input_prices = c("w1", "w2"), tech = tch)
    ref <- suppressWarnings(productivity::lowe(
      tfp_data, id.var = "firm", time.var = "year",
      x.vars = c("x1", "x2"), y.vars = c("y1", "y2"),
      w.vars = c("w1", "w2"), p.vars = c("p1", "p2"), rts = "vrs",
      orientation = "out", scaled = FALSE,
      tech.reg = tch == "contemporaneous"))$Levels
    lv <- align(tfp$levels)
    expect_equal(lv$TFP, ref$TFP, tolerance = 1e-6)
    for (cl in c("OTE", "OSE", "OME", "ROSE", "RME")) {
      expect_equal(lv[[paste0(cl, "_meta")]], ref[[cl]], tolerance = 1e-6,
                   label = paste(tch, cl))
    }
  }
})

test_that("indexes are transitive and the metafrontier decomposition is exact", {
  tfp <- tfp_index(cbind(y1, y2) ~ x1 + x2, data = tfp_data,
                   group = "region", time = "year", id = "firm")
  lv <- tfp$levels
  # Identity, transitivity and multiplicative completeness: an index
  # between observations is the ratio of their TFP levels
  a <- 1; b <- 10; c <- 20
  tfpi <- function(i, j) lv$TFP[j] / lv$TFP[i]
  expect_equal(tfpi(a, a), 1)
  expect_equal(tfpi(a, c), tfpi(a, b) * tfpi(b, c))
  # TFP = TFP* x TFPE, TFPE = OTE x OME x ROSE = OTE x OSE x RME, for both
  # references (O'Donnell 2012)
  for (r in c("group", "meta")) {
    comp <- function(n) lv[[paste0(n, "_", r)]]
    expect_equal(lv$TFP, comp("TFP_star") * comp("TFPE"))
    expect_equal(comp("TFPE"), comp("OTE") * comp("OME") * comp("ROSE"))
    expect_equal(comp("OSME"), comp("OME") * comp("ROSE"))
    expect_equal(comp("TFPE"), comp("OTE") * comp("OSE") * comp("RME"))
  }
  # Metatechnology ratio: TFPE_meta = TFPE_group x MTR, MTR in (0, 1]
  expect_true(all(lv$MTR > 0 & lv$MTR <= 1 + 1e-9))
  expect_equal(lv$TFPE_meta, lv$TFPE_group * lv$MTR)
  # Changes are ratios of consecutive levels for the same unit
  ch <- tfp$changes
  r1 <- ch[ch$firm == 3 & ch$period_from == 1, ]
  expect_equal(r1$dTFP, lv$TFP[lv$firm == 3 & lv$year == 2] /
                 lv$TFP[lv$firm == 3 & lv$year == 1])
  expect_equal(nrow(ch), 8 * 2)
})

test_that("group-reference radial efficiency matches productivity on the group alone", {
  skip_if_not_installed("productivity")
  tfp <- tfp_index(cbind(y1, y2) ~ x1 + x2, data = tfp_data,
                   group = "region", time = "year", id = "firm")
  north <- tfp_data[tfp_data$region == "North", ]
  ref <- suppressWarnings(productivity::fareprim(
    north, id.var = "firm", time.var = "year", x.vars = c("x1", "x2"),
    y.vars = c("y1", "y2"), rts = "vrs", orientation = "out",
    scaled = FALSE))$Levels
  lv <- align(tfp$levels[tfp$levels$group == "North", ])
  # OTE and OSE do not depend on the aggregation weights
  expect_equal(lv$OTE_group, ref$OTE, tolerance = 1e-6)
  expect_equal(lv$OSE_group, ref$OSE, tolerance = 1e-6)
})

test_that("tfp_index drops incomplete rows and validates prices", {
  d_na <- tfp_data
  d_na$x1[5] <- NA
  d_na$region[20] <- NA
  tfp_na <- tfp_index(cbind(y1, y2) ~ x1 + x2, data = d_na,
                      group = "region", time = "year", id = "firm")
  tfp_cl <- tfp_index(cbind(y1, y2) ~ x1 + x2, data = d_na[-c(5, 20), ],
                      group = "region", time = "year", id = "firm")
  expect_equal(tfp_na$levels, tfp_cl$levels)
  expect_error(tfp_index(cbind(y1, y2) ~ x1 + x2, data = tfp_data,
                         group = "region", time = "year", id = "firm",
                         index = "lowe"), "output_prices")
  expect_error(tfp_index(cbind(y1, y2) ~ x1 + x2, data = tfp_data,
                         group = "region", time = "year", id = "firm",
                         index = "lowe", output_prices = "p1",
                         input_prices = c("w1", "w2")),
               "one price column per output")
  expect_output(print(tfp_index(cbind(y1, y2) ~ x1 + x2, data = tfp_data,
                                group = "region", time = "year",
                                id = "firm")), "Fare-Primont")
})

test_that("autoplot.tfp_meta draws level paths and the metafrontier decomposition", {
  skip_if_not_installed("ggplot2")
  tfp <- tfp_index(cbind(y1, y2) ~ x1 + x2, data = tfp_data,
                   group = "region", time = "year", id = "firm")
  p1 <- ggplot2::autoplot(tfp, which = "levels")
  p2 <- ggplot2::autoplot(tfp, which = "decomposition")
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_true(all(c("TFPE_group", "MTR", "TFPE_meta") %in%
                    unique(p2$data$component)))
})

test_that("summary.tfp_meta tabulates levels and changes by group and period", {
  tfp <- tfp_index(cbind(y1, y2) ~ x1 + x2, data = tfp_data,
                   group = "region", time = "year", id = "firm")
  s <- summary(tfp)
  expect_equal(nrow(s$levels), 2 * 3)
  expect_equal(nrow(s$changes), 2 * 2)
  expect_output(print(s), "Mean changes")
})
