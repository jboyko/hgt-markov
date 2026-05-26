test_that("simulate_tree returns a phylo with non-ultrametric branch lengths", {
  skip_if_not_installed("TreeSim")
  tree <- simulate_tree(n_target = 20, seed = 1)

  expect_s3_class(tree, "phylo")
  expect_false(ape::is.ultrametric(tree, tol = 1e-8),
    label = "BD-sampling tree should be non-ultrametric")
})

test_that("simulate_tree tips have varied sampling times", {
  skip_if_not_installed("TreeSim")
  tree <- simulate_tree(n_target = 20, seed = 1)

  depths <- ape::node.depth.edgelength(tree)
  tip_depths <- depths[seq_len(length(tree$tip.label))]
  # At least two distinct tip depths
  expect_gt(length(unique(round(tip_depths, 6))), 1L,
    label = "tips should have varied sampling depths")
})

test_that("simulate_tree is reproducible with the same seed", {
  skip_if_not_installed("TreeSim")
  t1 <- simulate_tree(n_target = 20, seed = 99)
  t2 <- simulate_tree(n_target = 20, seed = 99)

  expect_equal(t1$edge.length, t2$edge.length)
  expect_equal(t1$tip.label,   t2$tip.label)
})

test_that("tip branch contributes to HGT exposure only up to its sampling time", {
  # Manually construct: root -> t1 (depth 1.0), root -> t2 (depth 0.5)
  # Phase [0, 0.5]: k=2 → HGT possible
  # Phase [0.5, 1.0]: k=1 → no HGT possible
  # Expected total events = Poisson(lambda * 2 * 0.5) = Poisson(lambda)
  tree <- ape::read.tree(text = "(t1:1.0,t2:0.5);")
  Q <- matrix(c(-0.3, 0.3, 0.2, -0.2), nrow = 2, byrow = TRUE)
  rownames(Q) <- colnames(Q) <- c("0", "1")

  lambda <- 10.0
  expected_events <- lambda * 2 * 0.5  # k=2, dt=0.5 for the active phase

  n_reps <- 2000
  counts <- vapply(seq_len(n_reps), function(i)
    simulate_hgt(tree, Q, lambda = lambda, seed = i)$n_hgt_events,
    integer(1L))

  obs <- mean(counts)
  se  <- sqrt(expected_events / n_reps)
  expect_true(abs(obs - expected_events) < 3 * se,
    label = sprintf(
      "obs=%.3f expected=%.3f 3se=%.4f — t2 branch should not contribute after its sampling time",
      obs, expected_events, 3 * se))
})

test_that("run_one_rep runs end-to-end on a freshly simulated BD-sampling tree at lambda > 0", {
  skip_if_not_installed("corHMM")
  skip_if_not_installed("TreeSim")

  Q <- matrix(c(-0.3, 0.3, 0.2, -0.2), nrow = 2, byrow = TRUE)
  rownames(Q) <- colnames(Q) <- c("0", "1")

  result <- run_one_rep(Q = Q, lambda = 0.5, seed = 7)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1L)
  expect_named(result, c("q01_true", "q10_true", "q01_hat", "q10_hat"), ignore.order = TRUE)
  expect_true(is.finite(result$q01_hat))
  expect_true(is.finite(result$q10_hat))
})
