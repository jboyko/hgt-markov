test_that("run_one_rep returns a list with $fits, $ard_q_hat, $node_probs_ard, $node_states_true", {
  skip_if_not_installed("corHMM")

  set.seed(1)
  tree <- ape::rtree(10)
  tree$edge.length <- tree$edge.length / max(ape::branching.times(tree))

  Q <- matrix(c(-0.3, 0.3, 0.2, -0.2), nrow = 2, byrow = TRUE)
  rownames(Q) <- colnames(Q) <- c("0", "1")

  result <- run_one_rep(tree = tree, Q = Q, lambda = 0, seed = 42)

  expect_type(result, "list")
  expect_named(result, c("fits", "ard_q_hat", "node_probs_ard", "node_states_true"),
               ignore.order = TRUE)
})

test_that("run_one_rep $fits is 4-row data frame with required columns", {
  skip_if_not_installed("corHMM")

  set.seed(1)
  tree <- ape::rtree(10)
  tree$edge.length <- tree$edge.length / max(ape::branching.times(tree))

  Q <- matrix(c(-0.3, 0.3, 0.2, -0.2), nrow = 2, byrow = TRUE)
  rownames(Q) <- colnames(Q) <- c("0", "1")

  result <- run_one_rep(tree = tree, Q = Q, lambda = 0, seed = 42)

  expect_s3_class(result$fits, "data.frame")
  expect_equal(nrow(result$fits), 4L)
  expect_named(
    result$fits,
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

  expect_setequal(result$fits$model, c("ER", "ARD", "HMM2", "HMM3"))
})

test_that("run_one_rep aicc_winner is one of the four models", {
  skip_if_not_installed("corHMM")

  set.seed(1)
  tree <- ape::rtree(10)
  tree$edge.length <- tree$edge.length / max(ape::branching.times(tree))

  Q <- matrix(c(-0.3, 0.3, 0.2, -0.2), nrow = 2, byrow = TRUE)
  rownames(Q) <- colnames(Q) <- c("0", "1")

  result <- run_one_rep(tree = tree, Q = Q, lambda = 0, seed = 42)

  expect_true(all(result$fits$aicc_winner %in% c("ER", "ARD", "HMM2", "HMM3")))
  expect_equal(length(unique(result$fits$aicc_winner)), 1L)
})

test_that("run_one_rep aicc values are finite for normal data", {
  skip_if_not_installed("corHMM")

  set.seed(3)
  tree <- ape::rtree(15)
  tree$edge.length <- tree$edge.length / max(ape::branching.times(tree))

  Q <- matrix(c(-0.5, 0.5, 0.3, -0.3), nrow = 2, byrow = TRUE)
  rownames(Q) <- colnames(Q) <- c("0", "1")

  result <- run_one_rep(tree = tree, Q = Q, lambda = 0, seed = 7)

  converged_rows <- result$fits[result$fits$converged, ]
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

  expect_true(all(result$fits$q01_true == 0.3))
  expect_true(all(result$fits$q10_true == 0.2))
})

test_that("run_one_rep ard_q_hat contains finite rates on convergence", {
  skip_if_not_installed("corHMM")

  set.seed(1)
  tree <- ape::rtree(15)
  tree$edge.length <- tree$edge.length / max(ape::branching.times(tree))

  Q <- matrix(c(-0.5, 0.5, 0.3, -0.3), nrow = 2, byrow = TRUE)
  rownames(Q) <- colnames(Q) <- c("0", "1")

  result <- run_one_rep(tree = tree, Q = Q, lambda = 0, seed = 7)

  expect_named(result$ard_q_hat, c("q01", "q10"))
  expect_true(is.finite(result$ard_q_hat["q01"]))
  expect_true(is.finite(result$ard_q_hat["q10"]))
})

test_that("run_one_rep node_probs_ard has n_internal rows and 2 columns", {
  skip_if_not_installed("corHMM")

  set.seed(1)
  tree <- ape::rtree(15)
  tree$edge.length <- tree$edge.length / max(ape::branching.times(tree))

  Q <- matrix(c(-0.5, 0.5, 0.3, -0.3), nrow = 2, byrow = TRUE)
  rownames(Q) <- colnames(Q) <- c("0", "1")

  result <- run_one_rep(tree = tree, Q = Q, lambda = 0, seed = 7)

  expect_true(is.matrix(result$node_probs_ard))
  expect_equal(nrow(result$node_probs_ard), tree$Nnode)
  expect_equal(ncol(result$node_probs_ard), 2L)
})

test_that("run_one_rep node_states_true has length n_internal and values 0 or 1", {
  skip_if_not_installed("corHMM")

  set.seed(1)
  tree <- ape::rtree(15)
  tree$edge.length <- tree$edge.length / max(ape::branching.times(tree))

  Q <- matrix(c(-0.5, 0.5, 0.3, -0.3), nrow = 2, byrow = TRUE)
  rownames(Q) <- colnames(Q) <- c("0", "1")

  result <- run_one_rep(tree = tree, Q = Q, lambda = 0, seed = 7)

  expect_equal(length(result$node_states_true), tree$Nnode)
  expect_true(all(result$node_states_true %in% c(0L, 1L)))
})

test_that("run_one_rep works with lambda > 0 (mk_hgt generator)", {
  skip_if_not_installed("corHMM")

  set.seed(2)
  tree <- ape::rtree(10)
  tree$edge.length <- tree$edge.length / max(ape::branching.times(tree))

  Q <- matrix(c(-0.3, 0.3, 0.2, -0.2), nrow = 2, byrow = TRUE)
  rownames(Q) <- colnames(Q) <- c("0", "1")

  result <- run_one_rep(tree = tree, Q = Q, generator = "mk_hgt", lambda = 0.5, seed = 42)

  expect_type(result, "list")
  expect_equal(nrow(result$fits), 4L)
  expect_setequal(result$fits$model, c("ER", "ARD", "HMM2", "HMM3"))
})

test_that("ARD wins AICc more often than HMM3 on Mk-null strong-signal data", {
  skip_if_not_installed("corHMM")
  skip_on_cran()

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

  winners <- vapply(results, function(r) r$fits$aicc_winner[1], character(1))
  ard_wins  <- sum(winners == "ARD")
  hmm3_wins <- sum(winners == "HMM3")

  expect_gt(ard_wins, hmm3_wins)
})
