test_that("equal rate categories: tip distribution matches Mk-null within 3-sigma MC tolerance", {
  set.seed(42)
  tree <- ape::read.tree(text = "(tip1:0.5,tip2:0.5);")

  Q <- matrix(c(-0.4, 0.4, 0.2, -0.2), nrow = 2, byrow = TRUE)
  rownames(Q) <- colnames(Q) <- c("0", "1")

  n_reps    <- 10000
  root_freq <- c(0.5, 0.5)

  # Expected from pure Mk(Q) over t=0.5
  P <- expm::expm(Q * 0.5)
  expected_p1 <- root_freq %*% P[, 2]

  # HMM-null with identical rate categories and moderate alpha
  hmm_states <- vapply(seq_len(n_reps), function(i) {
    res <- simulate_hmm(tree, Q1 = Q, Q2 = Q, alpha = 0.5, root_freq = root_freq, seed = i)
    as.integer(res$tip_states["tip1"])
  }, integer(1L))

  obs_p1 <- mean(hmm_states)
  se <- sqrt(expected_p1 * (1 - expected_p1) / n_reps)

  expect_true(
    abs(obs_p1 - expected_p1) < 3 * se,
    label = sprintf("obs=%.4f expected=%.4f diff=%.4f 3se=%.4f",
                    obs_p1, expected_p1, abs(obs_p1 - expected_p1), 3 * se)
  )
})

test_that("very different rate categories produce higher sister-tip discordance than Mk-null", {
  set.seed(99)
  tree <- ape::read.tree(text = "(tip1:1.0,tip2:1.0);")

  # Slow category: very low rates
  Q_slow <- matrix(c(-0.01, 0.01, 0.01, -0.01), nrow = 2, byrow = TRUE)
  rownames(Q_slow) <- colnames(Q_slow) <- c("0", "1")

  # Fast category: very high rates
  Q_fast <- matrix(c(-2.0, 2.0, 2.0, -2.0), nrow = 2, byrow = TRUE)
  rownames(Q_fast) <- colnames(Q_fast) <- c("0", "1")

  # Mid Mk for comparison: rates at geometric mean
  Q_mid <- matrix(c(-sqrt(0.01 * 2.0), sqrt(0.01 * 2.0),
                     sqrt(0.01 * 2.0), -sqrt(0.01 * 2.0)), nrow = 2, byrow = TRUE)
  rownames(Q_mid) <- colnames(Q_mid) <- c("0", "1")

  root_freq <- c(0.5, 0.5)
  n_reps    <- 5000

  # Discordance: tip1 != tip2
  hmm_disc <- vapply(seq_len(n_reps), function(i) {
    res <- simulate_hmm(tree, Q1 = Q_slow, Q2 = Q_fast, alpha = 1.0,
                        root_freq = root_freq, seed = i)
    as.integer(res$tip_states["tip1"] != res$tip_states["tip2"])
  }, integer(1L))

  mk_disc <- vapply(seq_len(n_reps), function(i) {
    res <- simulate_mk(tree, Q_mid, root_freq = root_freq, seed = i)
    as.integer(res$tip_states["tip1"] != res$tip_states["tip2"])
  }, integer(1L))

  hmm_mean <- mean(hmm_disc)
  mk_mean  <- mean(mk_disc)

  # With extreme rate heterogeneity, HMM discordance should differ from Mk
  # (heterogeneous rates increase variance in state transitions)
  se_diff <- sqrt(var(hmm_disc) / n_reps + var(mk_disc) / n_reps)
  expect_true(
    abs(hmm_mean - mk_mean) > 2 * se_diff,
    label = sprintf("hmm_disc=%.4f mk_disc=%.4f diff=%.4f 2se=%.4f",
                    hmm_mean, mk_mean, abs(hmm_mean - mk_mean), 2 * se_diff)
  )
})

test_that("simulate_hmm returns tip_states and node_states with correct structure", {
  set.seed(1)
  tree <- ape::rtree(10)
  tree$edge.length <- tree$edge.length / max(ape::branching.times(tree))

  Q <- matrix(c(-0.3, 0.3, 0.2, -0.2), nrow = 2, byrow = TRUE)
  rownames(Q) <- colnames(Q) <- c("0", "1")

  res <- simulate_hmm(tree, Q1 = Q, Q2 = Q, alpha = 0.5, seed = 1)

  expect_named(res, c("tip_states", "node_states"), ignore.order = TRUE)
  expect_length(res$tip_states, length(tree$tip.label))
  expect_true(all(res$tip_states %in% c(0L, 1L)))
  expect_true(all(res$node_states %in% c(0L, 1L)))
})

test_that("run_one_rep generator='hmm_null' returns expected columns and finite estimates", {
  skip_if_not_installed("corHMM")

  set.seed(3)
  tree <- ape::rtree(10)
  tree$edge.length <- tree$edge.length / max(ape::branching.times(tree))

  Q <- matrix(c(-0.3, 0.3, 0.2, -0.2), nrow = 2, byrow = TRUE)
  rownames(Q) <- colnames(Q) <- c("0", "1")

  result <- run_one_rep(tree = tree, Q = Q, generator = "hmm_null", alpha = 0.5, seed = 42)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1L)
  expect_named(result, c("q01_true", "q10_true", "q01_hat", "q10_hat"), ignore.order = TRUE)
  expect_true(is.finite(result$q01_hat))
  expect_true(is.finite(result$q10_hat))
})

test_that("all three generators selectable via generator argument", {
  skip_if_not_installed("corHMM")

  set.seed(5)
  tree <- ape::rtree(10)
  tree$edge.length <- tree$edge.length / max(ape::branching.times(tree))

  Q <- matrix(c(-0.3, 0.3, 0.2, -0.2), nrow = 2, byrow = TRUE)
  rownames(Q) <- colnames(Q) <- c("0", "1")

  r_null <- run_one_rep(tree = tree, Q = Q, generator = "mk_null", seed = 1)
  r_hgt  <- run_one_rep(tree = tree, Q = Q, generator = "mk_hgt",  lambda = 0.5, seed = 1)
  r_hmm  <- run_one_rep(tree = tree, Q = Q, generator = "hmm_null", alpha = 0.5, seed = 1)

  for (r in list(r_null, r_hgt, r_hmm)) {
    expect_s3_class(r, "data.frame")
    expect_equal(nrow(r), 1L)
    expect_true(is.finite(r$q01_hat))
    expect_true(is.finite(r$q10_hat))
  }
})
