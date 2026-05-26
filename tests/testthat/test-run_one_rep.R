test_that("run_one_rep returns one-row data frame with expected columns and finite estimates", {
  skip_if_not_installed("corHMM")

  set.seed(1)
  tree <- ape::rtree(10)
  tree$edge.length <- tree$edge.length / max(ape::branching.times(tree))  # scale to depth 1

  Q <- matrix(c(-0.3, 0.3, 0.2, -0.2), nrow = 2, byrow = TRUE)
  rownames(Q) <- colnames(Q) <- c("0", "1")

  result <- run_one_rep(tree = tree, Q = Q, lambda = 0, seed = 42)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1L)
  expect_named(result, c("q01_true", "q10_true", "q01_hat", "q10_hat"), ignore.order = TRUE)
  expect_true(is.finite(result$q01_hat))
  expect_true(is.finite(result$q10_hat))
  expect_equal(result$q01_true, 0.3)
  expect_equal(result$q10_true, 0.2)
})

test_that("run_one_rep works with lambda > 0", {
  skip_if_not_installed("corHMM")

  set.seed(2)
  tree <- ape::rtree(10)
  tree$edge.length <- tree$edge.length / max(ape::branching.times(tree))

  Q <- matrix(c(-0.3, 0.3, 0.2, -0.2), nrow = 2, byrow = TRUE)
  rownames(Q) <- colnames(Q) <- c("0", "1")

  result <- run_one_rep(tree = tree, Q = Q, lambda = 0.5, seed = 42)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1L)
  expect_named(result, c("q01_true", "q10_true", "q01_hat", "q10_hat"), ignore.order = TRUE)
  expect_true(is.finite(result$q01_hat))
  expect_true(is.finite(result$q10_hat))
})
