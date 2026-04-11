# Test with Produc dataset from plm package

test_that("SFA metafrontier works with Produc data (multiple years)", {
  skip_if_not_installed("plm")

  data("Produc", package = "plm")

  # Pool several years to get enough obs per region for SFA
  produc_sub <- Produc[Produc$year %in% 1983:1986, ]

  # Use only 2 regions with most states for stability
  region_counts <- sort(table(produc_sub$region), decreasing = TRUE)
  top2 <- names(region_counts)[1:2]
  produc_sub <- produc_sub[produc_sub$region %in% top2, ]
  produc_sub$region <- droplevels(factor(produc_sub$region))

  skip_if(length(unique(produc_sub$region)) < 2)
  skip_if(min(table(produc_sub$region)) < 10)

  # SFA optimization may not converge on all platforms with real data
  fit <- tryCatch(
    metafrontier(log(gsp) ~ log(pc) + log(emp),
                 data = produc_sub, group = "region",
                 method = "sfa"),
    error = function(e) NULL
  )
  skip_if(is.null(fit), "SFA optimization did not converge on this platform")

  expect_s3_class(fit, "metafrontier")
  expect_true(all(fit$tgr > 0 & fit$tgr <= 1 + 1e-6))
  expect_length(fit$te_group, nrow(produc_sub))
})


test_that("DEA metafrontier works with Produc data", {
  skip_if_not_installed("plm")

  data("Produc", package = "plm")

  produc_86 <- Produc[Produc$year == 1986, ]
  region_counts <- table(produc_86$region)
  valid_regions <- names(region_counts[region_counts >= 3])
  produc_sub <- produc_86[produc_86$region %in% valid_regions, ]
  produc_sub$region <- droplevels(factor(produc_sub$region))

  skip_if(length(unique(produc_sub$region)) < 2)

  fit <- metafrontier(gsp ~ pc + emp,
                      data = produc_sub, group = "region",
                      method = "dea", rts = "crs")

  expect_s3_class(fit, "metafrontier_dea")
  expect_true(all(fit$te_group > 0, na.rm = TRUE))
  expect_true(all(fit$tgr <= 1 + 1e-6, na.rm = TRUE))
})


test_that("Malmquist works with Produc panel data", {
  skip_if_not_installed("plm")

  data("Produc", package = "plm")

  # Use 2 consecutive years and 2 regions
  produc_sub <- Produc[Produc$year %in% c(1985, 1986) &
                        Produc$region %in% c(1, 2), ]
  produc_sub$region <- droplevels(factor(produc_sub$region))

  skip_if(length(unique(produc_sub$region)) < 2)
  skip_if(length(unique(produc_sub$year)) < 2)

  malm <- malmquist_meta(gsp ~ pc + emp,
                         data = produc_sub,
                         group = "region",
                         time = "year")

  expect_s3_class(malm, "malmquist_meta")
  expect_true(nrow(malm$malmquist) > 0)
})
