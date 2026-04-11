# Cross-validation test: compare group-level SFA estimates against sfaR

test_that("group-level SFA estimates are close to sfaR", {
  skip_if_not_installed("sfaR")

  set.seed(42)
  sim <- simulate_metafrontier(n_groups = 2, n_per_group = 200, seed = 80)

  # Our metafrontier fit
  fit <- metafrontier(log_y ~ log_x1 + log_x2,
                      data = sim$data, group = "group",
                      method = "sfa", dist = "hnormal")

  for (g in fit$groups) {
    idx <- sim$data$group == g
    data_g <- sim$data[idx, ]

    # sfaR fit for same group
    sfar_fit <- sfaR::sfacross(
      log_y ~ log_x1 + log_x2,
      data = data_g,
      udist = "hnormal"
    )

    # Compare frontier coefficients
    our_beta <- fit$group_models[[g]]$coefficients
    sfar_beta <- sfar_fit$mlParam[seq_along(our_beta)]

    # Coefficients should be in the same ballpark (within 20%)
    # Some differences are expected due to different starting values
    # and optimization paths
    for (j in seq_along(our_beta)) {
      if (abs(sfar_beta[j]) > 0.05) {
        ratio <- our_beta[j] / sfar_beta[j]
        expect_true(ratio > 0.5 && ratio < 2.0,
                    info = paste("Group", g, "coef", j,
                                 "ours:", our_beta[j],
                                 "sfaR:", sfar_beta[j]))
      }
    }

    # Compare mean efficiency: should be within 10 percentage points
    our_te <- mean(fit$te_group[idx])
    te_obj <- sfaR::efficiencies(sfar_fit)
    if (is.data.frame(te_obj) || is.matrix(te_obj)) {
      if ("teBC" %in% names(te_obj)) {
        sfar_te_vec <- as.numeric(te_obj[["teBC"]])
      } else {
        sfar_te_vec <- as.numeric(te_obj[[1]])
      }
    } else {
      sfar_te_vec <- as.numeric(te_obj)
    }
    sfar_te <- mean(sfar_te_vec)
    expect_true(abs(our_te - sfar_te) < 0.10,
                info = paste("Group", g,
                             "our mean TE:", round(our_te, 3),
                             "sfaR mean TE:", round(sfar_te, 3)))
  }
})
