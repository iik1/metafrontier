# Tests for Latent Class Metafrontier (P7)

test_that("latent_class_metafrontier runs and returns correct class", {
  lc <- latent_class_metafrontier(
    log_y ~ log_x1 + log_x2,
    data = test_data,
    n_classes = 2,
    n_starts = 3,
    max_iter = 50,
    seed = 42
  )

  expect_s3_class(lc, "lc_metafrontier")
  expect_equal(lc$n_classes, 2)
  expect_equal(length(lc$class_assignment), nrow(test_data))
  expect_true(is.matrix(lc$posterior))
  expect_equal(ncol(lc$posterior), 2)
  expect_true(lc$marginal_ll > -Inf)
  expect_true(is.finite(lc$BIC))
})

test_that("posterior probabilities sum to 1", {
  lc <- latent_class_metafrontier(
    log_y ~ log_x1 + log_x2,
    data = test_data,
    n_classes = 2,
    n_starts = 2,
    max_iter = 30,
    seed = 123
  )

  row_sums <- rowSums(lc$posterior)
  expect_equal(row_sums, rep(1, nrow(test_data)), tolerance = 1e-10)
})

test_that("class proportions sum to 1", {
  lc <- latent_class_metafrontier(
    log_y ~ log_x1 + log_x2,
    data = test_data,
    n_classes = 2,
    n_starts = 2,
    max_iter = 30,
    seed = 456
  )

  expect_equal(sum(lc$pi), 1.0, tolerance = 1e-6)
})

test_that("metafrontier object is accessible", {
  lc <- latent_class_metafrontier(
    log_y ~ log_x1 + log_x2,
    data = test_data,
    n_classes = 2,
    n_starts = 2,
    max_iter = 30,
    seed = 789
  )

  expect_s3_class(lc$metafrontier, "metafrontier")
  mc <- coef(lc, which = "meta")
  expect_true(is.numeric(mc))
  expect_length(mc, 3)
})

test_that("print produces output", {
  lc <- latent_class_metafrontier(
    log_y ~ log_x1 + log_x2,
    data = test_data,
    n_classes = 2,
    n_starts = 2,
    max_iter = 30,
    seed = 42
  )

  expect_output(print(lc), "Latent Class Metafrontier")
  expect_output(print(lc), "BIC")
})

test_that("3 classes works", {
  lc <- latent_class_metafrontier(
    log_y ~ log_x1 + log_x2,
    data = test_data,
    n_classes = 3,
    n_starts = 2,
    max_iter = 30,
    seed = 42
  )

  expect_equal(lc$n_classes, 3)
  expect_equal(ncol(lc$posterior), 3)
})
