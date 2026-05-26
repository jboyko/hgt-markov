test_that("run_one_rep returns 4-row data frame with required columns", {
  skip_if_not_installed("corHMM")

  set.seed(1)
  tree <- ape::rtree(10)
  tree$edge.length <- tree$edge.length / max(ape::branching.times(tree))

  Q <- matrix(c(-0.3, 0.3, 0.2, -0.2), nrow = 2, byrow = TRUE)
  rownames(Q) <- colnames(Q) <- c("0", "1")

  result <- run_one_rep(tree = tree, Q = Q, lambda = 0, seed = 42)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 4L)
  expect_named(
    result,
    c("model", "aicc", "loglik", "converged", "conv_note", "aicc_winner", "q01_true", "q10_true"),
    ignore.order = TRUE
  )
})

test_that("run_one_rep fits all four models", {
  skip_if_not_installed("corHMM")

  set.seed(1)
  tree <- ape::rtree(10)
  tree$edge.length <- tree$edge.length / max(ape::branching.times(tree))

  Q <- matrix(c(-0.3, 0.3, 0.2, -0.2), nrow = 2, byrow = TRUE)
  rownames(Q) <- colnames(Q) <- c("0", "1")

  result <- run_one_rep(tree = tree, Q = Q, lambda = 0, seed = 42)

  expect_setequal(result$model, c("ER", "ARD", "HMM2", "HMM3"))
})

test_that("run_one_rep aicc_winner is one of the four models", {
  skip_if_not_installed("corHMM")

  set.seed(1)
  tree <- ape::rtree(10)
  tree$edge.length <- tree$edge.length / max(ape::branching.times(tree))

  Q <- matrix(c(-0.3, 0.3, 0.2, -0.2), nrow = 2, byrow = TRUE)
  rownames(Q) <- colnames(Q) <- c("0", "1")

  result <- run_one_rep(tree = tree, Q = Q, lambda = 0, seed = 42)

  expect_true(all(result$aicc_winner %in% c("ER", "ARD", "HMM2", "HMM3")))
  # all rows report the same winner
  expect_equal(length(unique(result$aicc_winner)), 1L)
})

test_that("run_one_rep aicc values are finite for normal data", {
  skip_if_not_installed("corHMM")

  set.seed(3)
  tree <- ape::rtree(15)
  tree$edge.length <- tree$edge.length / max(ape::branching.times(tree))

  Q <- matrix(c(-0.5, 0.5, 0.3, -0.3), nrow = 2, byrow = TRUE)
  rownames(Q) <- colnames(Q) <- c("0", "1")

  result <- run_one_rep(tree = tree, Q = Q, lambda = 0, seed = 7)

  converged_rows <- result[result$converged, ]
  expect_true(all(is.finite(converged_rows$aicc)))
  expect_true(all(is.finite(converged_rows$loglik)))
})

test_that("run_one_rep true params are correct in all rows", {
  skip_if_not_installed("corHMM")

  set.seed(1)
  tree <- ape::rtree(10)
  tree$edge.length <- tree$edge.length / max(ape::branching.times(tree))

  Q <- matrix(c(-0.3, 0.3, 0.2, -0.2), nrow = 2, byrow = TRUE)
  rownames(Q) <- colnames(Q) <- c("0", "1")

  result <- run_one_rep(tree = tree, Q = Q, lambda = 0, seed = 42)

  expect_true(all(result$q01_true == 0.3))
  expect_true(all(result$q10_true == 0.2))
})

test_that("run_one_rep works with lambda > 0 (mk_hgt generator)", {
  skip_if_not_installed("corHMM")

  set.seed(2)
  tree <- ape::rtree(10)
  tree$edge.length <- tree$edge.length / max(ape::branching.times(tree))

  Q <- matrix(c(-0.3, 0.3, 0.2, -0.2), nrow = 2, byrow = TRUE)
  rownames(Q) <- colnames(Q) <- c("0", "1")

  result <- run_one_rep(tree = tree, Q = Q, generator = "mk_hgt", lambda = 0.5, seed = 42)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 4L)
  expect_setequal(result$model, c("ER", "ARD", "HMM2", "HMM3"))
})

test_that("ARD wins AICc more often than HMM3 on Mk-null strong-signal data", {
  skip_if_not_installed("corHMM")
  skip_on_cran()

  # Strong signal: fast rates, large tree
  Q <- matrix(c(-1.0, 1.0, 0.8, -0.8), nrow = 2, byrow = TRUE)
  rownames(Q) <- colnames(Q) <- c("0", "1")

  set.seed(99)
  trees <- lapply(1:20, function(i) {
    tr <- ape::rtree(40)
    tr$edge.length <- tr$edge.length / max(ape::branching.times(tr))
    tr
  })

  results <- lapply(seq_along(trees), function(i) {
    run_one_rep(tree = trees[[i]], Q = Q, generator = "mk_null", lambda = 0, seed = i)
  })

  winners <- vapply(results, function(r) r$aicc_winner[1], character(1))
  ard_wins  <- sum(winners == "ARD")
  hmm3_wins <- sum(winners == "HMM3")

  expect_gt(ard_wins, hmm3_wins)
})
