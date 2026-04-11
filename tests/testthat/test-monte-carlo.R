# Monte Carlo parameter recovery test

test_that("SFA metafrontier recovers true parameters (Monte Carlo)", {
  skip_on_cran()

  n_rep <- 30
  n_per_group <- 200
  true_beta <- c(1.0, 0.5, 0.3)
  true_gap <- c(0, 0.4)
  true_sigma_u <- c(0.2, 0.3)
  true_sigma_v <- 0.15

  beta_estimates <- matrix(NA, n_rep, length(true_beta))
  tgr_mean_g2 <- numeric(n_rep)

  for (r in seq_len(n_rep)) {
    sim <- simulate_metafrontier(
      n_groups = 2,
      n_per_group = n_per_group,
      beta_meta = true_beta,
      tech_gap = true_gap,
      sigma_u = true_sigma_u,
      sigma_v = true_sigma_v,
      seed = 1000 + r
    )

    fit <- tryCatch(
      metafrontier(log_y ~ log_x1 + log_x2,
                   data = sim$data, group = "group",
                   method = "sfa", meta_type = "deterministic"),
      error = function(e) NULL
    )

    if (!is.null(fit) && !is.null(fit$meta_coef)) {
      beta_estimates[r, ] <- fit$meta_coef
      idx_g2 <- sim$data$group == "G2"
      tgr_mean_g2[r] <- mean(fit$tgr[idx_g2])
    }
  }

  # Remove failed runs
  valid <- complete.cases(beta_estimates)
  expect_true(sum(valid) >= n_rep * 0.8,
              info = paste("Only", sum(valid), "of", n_rep,
                           "replications converged"))

  beta_est <- beta_estimates[valid, ]
  tgr_est <- tgr_mean_g2[valid]

  # Check bias: mean estimate should be within 30% of true value
  mean_beta <- colMeans(beta_est)
  for (j in seq_along(true_beta)) {
    if (abs(true_beta[j]) > 0.01) {
      bias_pct <- abs(mean_beta[j] - true_beta[j]) / abs(true_beta[j])
      expect_true(bias_pct < 0.30,
                  info = paste("Beta", j, "bias:", round(bias_pct * 100, 1),
                               "%, mean:", round(mean_beta[j], 3),
                               "true:", true_beta[j]))
    }
  }

  # TGR for G2 should reflect the technology gap
  # True TGR for G2 = exp(-0.4) ~ 0.67
  true_tgr_g2 <- exp(-true_gap[2])
  mean_tgr_g2 <- mean(tgr_est)
  expect_true(abs(mean_tgr_g2 - true_tgr_g2) < 0.15,
              info = paste("Mean TGR G2:", round(mean_tgr_g2, 3),
                           "True:", round(true_tgr_g2, 3)))
})
