test_that("lambda=0 tip distribution matches expm(Q*t) within 3-sigma MC tolerance", {
  set.seed(42)

  # 2-tip tree: root -> tip1 (t=0.5), root -> tip2 (t=0.5)
  tree <- ape::read.tree(text = "((tip1:0.5,tip2:0.5));")
  # drop the root-to-root edge; use a simple 2-tip balanced tree
  tree <- ape::read.tree(text = "(tip1:0.5,tip2:0.5);")

  Q <- matrix(c(-0.4, 0.4, 0.2, -0.2), nrow = 2, byrow = TRUE)
  rownames(Q) <- colnames(Q) <- c("0", "1")

  n_reps <- 10000
  root_freq <- c(0.5, 0.5)

  # Expected tip1 state distribution: marginalise over root state
  P <- expm::expm(Q * 0.5)
  # P[i,j] = prob of going from state i to j along t=0.5
  # marginal tip1 prob of state 1 (index 2): sum_i root_freq[i] * P[i,2]
  expected_p1 <- root_freq %*% P[, 2]  # prob tip1 is in state "1"

  # Simulate n_reps replicates, record tip1 state each time
  tip1_states <- vapply(seq_len(n_reps), function(i) {
    res <- simulate_mk(tree, Q, root_freq = root_freq, seed = i)
    as.integer(res$tip_states["tip1"])
  }, integer(1L))

  obs_p1 <- mean(tip1_states)
  se <- sqrt(expected_p1 * (1 - expected_p1) / n_reps)

  expect_true(
    abs(obs_p1 - expected_p1) < 3 * se,
    label = sprintf("obs=%.4f expected=%.4f diff=%.4f 3se=%.4f",
                    obs_p1, expected_p1, abs(obs_p1 - expected_p1), 3 * se)
  )
})
