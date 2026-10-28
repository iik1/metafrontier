# Cost frontiers: ln C = x'beta + v + u, efficiency exp(-u) in (0, 1]

make_cost_data <- function(n = 200, seed = 11) {
  set.seed(seed)
  do.call(rbind, lapply(c("A", "B"), function(g) {
    log_y <- log(runif(n, 1, 10))
    log_w <- log(runif(n, 1, 5))
    shift <- if (g == "A") 0 else 0.3   # B has higher minimum cost
    log_c <- 0.5 + shift + 0.8 * log_y + 0.6 * log_w +
      rnorm(n, 0, 0.1) + abs(rnorm(n, 0, 0.25))
    data.frame(group = g, log_c = log_c, log_y = log_y, log_w = log_w)
  }))
}
cost_data <- make_cost_data()

test_that("cost group frontiers match sfaR::sfacross(S = -1)", {
  skip_if_not_installed("sfaR")
  fit <- metafrontier(log_c ~ log_y + log_w, data = cost_data,
                      group = "group", frontier = "cost")
  for (g in c("A", "B")) {
    ref <- sfaR::sfacross(log_c ~ log_y + log_w, udist = "hnormal",
                          data = cost_data[cost_data$group == g, ], S = -1)
    gm <- fit$group_models[[g]]
    expect_equal(unname(gm$coefficients), unname(coef(ref)[1:3]),
                 tolerance = 1e-3)
    expect_equal(unname(gm$sigma_v), exp(coef(ref)[["Zv_(Intercept)"]] / 2),
                 tolerance = 1e-3)
    expect_equal(unname(gm$sigma_u), exp(coef(ref)[["Zu_(Intercept)"]] / 2),
                 tolerance = 1e-3)
    expect_equal(gm$efficiency, sfaR::efficiencies(ref)$teBC,
                 tolerance = 1e-3)
  }
  expect_identical(fit$frontier, "cost")
  expect_true(all(fit$te_group > 0 & fit$te_group <= 1))
})

test_that("deterministic cost metafrontier envelops the group frontiers from below", {
  for (obj in c("lp", "qp")) {
    fit <- metafrontier(log_c ~ log_y + log_w, data = cost_data,
                        group = "group", frontier = "cost", objective = obj)
    expect_true(all(fit$meta_frontier <= fit$group_frontier + 1e-8))
    expect_true(all(fit$tgr > 0 & fit$tgr <= 1))
    # The envelope touches the lower-cost group A somewhere
    expect_equal(max(fit$tgr[cost_data$group == "A"]), 1, tolerance = 1e-6)
    # Group B's minimum cost is about exp(0.3) above A's
    expect_equal(mean(fit$tgr[cost_data$group == "B"]), exp(-0.3),
                 tolerance = 0.05)
    expect_equal(fit$te_meta, fit$te_group * fit$tgr)
  }
})

test_that("stochastic cost TGR is E[exp(-u*)|eps*] with eps* = v* + u*", {
  fit <- metafrontier(log_c ~ log_y + log_w, data = cost_data,
                      group = "group", frontier = "cost",
                      meta_type = "stochastic")
  p <- fit$meta_opt$par
  sv <- exp(p[["log_sigma_v"]])
  su <- exp(p[["log_sigma_u"]])
  eps <- fit$group_frontier - fit$meta_frontier
  mu_star <- eps * su^2 / (sv^2 + su^2)
  sigma_star <- sv * su / sqrt(sv^2 + su^2)
  expected <- exp(-mu_star + sigma_star^2 / 2) *
    pnorm(mu_star / sigma_star - sigma_star) / pnorm(mu_star / sigma_star)

  expect_equal(fit$tgr, expected, tolerance = 1e-8)
  expect_true(all(fit$tgr > 0 & fit$tgr <= 1))
  # Stage 2 coefficients and their covariance are on the cost scale
  expect_equal(unname(fit$meta_coef), unname(p[1:3]))
  expect_true(all(diag(vcov(fit)) > 0))
})

test_that("Murphy-Topel on a cost fit mirrors the production fit to -ln C", {
  fit_c <- metafrontier(log_c ~ log_y + log_w, data = cost_data,
                        group = "group", frontier = "cost",
                        meta_type = "stochastic")
  fit_p <- metafrontier(I(-log_c) ~ log_y + log_w, data = cost_data,
                        group = "group", meta_type = "stochastic")
  v_c <- vcov(fit_c, correction = "murphy-topel", extraPar = TRUE)
  v_p <- vcov(fit_p, correction = "murphy-topel", extraPar = TRUE)
  J <- c(-1, -1, -1, 1, 1)
  expect_equal(unname(v_c), unname(v_p * outer(J, J)), tolerance = 1e-6)
})

test_that("simulate_metafrontier(frontier = 'cost') generates recoverable cost data", {
  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 400,
                               tech_gap = c(0, 0.4), sigma_u = c(0.3, 0.3),
                               sigma_v = 0.1, frontier = "cost", seed = 5)
  d <- sim$data
  # Costs lie above the group cost frontier: ln C - x'beta_g = v + u
  expect_gt(mean(d$true_u), 0)
  expect_equal(unique(round(d$true_tgr, 10)), round(exp(-c(0, 0.4)), 10))
  expect_identical(sim$params$frontier, "cost")

  fit <- metafrontier(log_y ~ log_x1 + log_x2, data = d, group = "group",
                      frontier = "cost")
  expect_equal(unname(fit$group_coef$G2), sim$params$beta_groups$G2,
               tolerance = 0.1)
  tgr_g2 <- mean(fit$tgr[d$group == "G2"])
  expect_equal(tgr_g2, exp(-0.4), tolerance = 0.05)
  expect_gt(cor(fit$te_group, d$true_te), 0.7)
})

test_that("parametric boot_tgr on a cost fit draws cost data and refits cost", {
  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 200,
                               tech_gap = c(0, 0.4), sigma_v = 0.1,
                               frontier = "cost", seed = 6)
  fit <- metafrontier(log_y ~ log_x1 + log_x2, data = sim$data,
                      group = "group", frontier = "cost")
  b <- boot_tgr(fit, R = 30, seed = 1, progress = FALSE)
  expect_equal(b$R_effective, 30L)
  ci_g2 <- unlist(b$ci_group[b$ci_group$Group == "G2", 3:4])
  expect_true(ci_g2[1] <= mean(fit$tgr[sim$data$group == "G2"]) &&
                mean(fit$tgr[sim$data$group == "G2"]) <= ci_g2[2])
  expect_true(all(b$tgr_boot <= 1))
})

test_that("poolability_test compares cost frontiers against a pooled cost frontier", {
  same <- simulate_metafrontier(n_groups = 2, n_per_group = 300,
                                tech_gap = c(0, 0), sigma_v = 0.1,
                                frontier = "cost", seed = 8)$data
  fit_same <- metafrontier(log_y ~ log_x1 + log_x2, data = same,
                           group = "group", frontier = "cost")
  expect_gt(poolability_test(fit_same)$p.value, 0.01)

  gap <- simulate_metafrontier(n_groups = 2, n_per_group = 300,
                               tech_gap = c(0, 0.4), sigma_v = 0.1,
                               frontier = "cost", seed = 8)$data
  fit_gap <- metafrontier(log_y ~ log_x1 + log_x2, data = gap,
                          group = "group", frontier = "cost")
  expect_lt(poolability_test(fit_gap)$p.value, 0.001)
})

test_that("imported and external-engine cost models give a cost metafrontier", {
  skip_if_not_installed("sfaR")
  skip_if_not_installed("frontier")
  fit_sfar <- metafrontier(log_c ~ log_y + log_w, data = cost_data,
                           group = "group", frontier = "cost",
                           engine = "sfaR")
  mods <- lapply(c(A = "A", B = "B"), function(g) {
    sfaR::sfacross(log_c ~ log_y + log_w, udist = "hnormal",
                   data = cost_data[cost_data$group == g, ], S = -1)
  })
  fit_mods <- metafrontier(models = mods)
  expect_identical(fit_mods$frontier, "cost")
  expect_equal(fit_mods$tgr, fit_sfar$tgr, tolerance = 1e-8)

  fit_front <- metafrontier(log_c ~ log_y + log_w, data = cost_data,
                            group = "group", frontier = "cost",
                            engine = "frontier")
  expect_equal(unname(fit_front$group_coef$B),
               unname(fit_sfar$group_coef$B), tolerance = 1e-2)
  expect_true(all(fit_front$te_group > 0 & fit_front$te_group <= 1))

  mods$A <- sfaR::sfacross(log_c ~ log_y + log_w, udist = "hnormal",
                           data = cost_data[cost_data$group == "A", ],
                           S = 1)
  expect_error(metafrontier(models = mods), "production and cost")
})

test_that("BC92 panel cost frontiers match frontier::sfa(ineffDecrease = FALSE)", {
  skip_if_not_installed("frontier")
  skip_if_not_installed("plm")
  p <- simulate_panel_metafrontier(n_groups = 2, n_firms_per_group = 40,
                                   n_periods = 5, seed = 4)$data
  # ln C = -ln y is a cost frontier with slopes -beta and u raising cost
  p$log_c <- -p$log_y
  fit <- metafrontier(log_c ~ log_x1 + log_x2, data = p, group = "group",
                      panel = list(id = "firm", time = "year"),
                      frontier = "cost")
  for (g in fit$groups) {
    pd <- plm::pdata.frame(p[p$group == g, ], index = c("firm", "year"))
    ref <- frontier::sfa(log_c ~ log_x1 + log_x2, data = pd,
                         timeEffect = TRUE, ineffDecrease = FALSE)
    gm <- fit$group_models[[g]]
    expect_equal(unname(gm$coefficients), unname(coef(ref)[1:3]),
                 tolerance = 1e-3)
    expect_equal(gm$eta, unname(coef(ref)["time"]), tolerance = 1e-2)
    # pdata.frame sorts rows by firm and year
    pg <- p[p$group == g, ]
    expect_equal(gm$efficiency[order(pg$firm, pg$year)],
                 as.numeric(frontier::efficiencies(ref, asInData = TRUE)),
                 tolerance = 1e-3)
  }
})

test_that("unsupported cost combinations fail loudly and output names the frontier", {
  expect_error(
    metafrontier(exp(log_c) ~ exp(log_y) + exp(log_w), data = cost_data,
                 group = "group", method = "dea", frontier = "cost"),
    "SFA"
  )
  expect_error(
    latent_class_metafrontier(log_c ~ log_y + log_w, data = cost_data,
                              frontier = "cost"),
    "production"
  )
  expect_error(
    malmquist_meta(log_c ~ log_y + log_w,
                   data = transform(cost_data, t = rep(1:2, 200), id = 1:400),
                   group = "group", time = "t", id = "id",
                   frontier = "cost"),
    "production"
  )
  fit <- metafrontier(log_c ~ log_y + log_w, data = cost_data,
                      group = "group", frontier = "cost")
  expect_output(print(fit), "Frontier:[[:space:]]+cost")
  expect_output(print(summary(fit)), "Frontier:[[:space:]]+cost")
})
