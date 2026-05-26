test_that("lambda=0 simulate_hgt marginal tip distribution matches simulate_mk within 3-sigma", {
  tree <- ape::read.tree(text = "(tip1:0.5,tip2:0.5);")
  Q <- matrix(c(-0.4, 0.4, 0.2, -0.2), nrow = 2, byrow = TRUE)
  rownames(Q) <- colnames(Q) <- c("0", "1")

  n_reps <- 5000
  mk_tip1 <- vapply(seq_len(n_reps), function(i)
    simulate_mk(tree, Q, root_freq = c(0.5, 0.5), seed = i)$tip_states["tip1"],
    integer(1L))
  hgt_tip1 <- vapply(seq_len(n_reps), function(i)
    simulate_hgt(tree, Q, lambda = 0, root_freq = c(0.5, 0.5), seed = i)$tip_states["tip1"],
    integer(1L))

  p_mk  <- mean(mk_tip1)
  p_hgt <- mean(hgt_tip1)
  se    <- sqrt(p_mk * (1 - p_mk) / n_reps)
  expect_true(abs(p_mk - p_hgt) < 3 * se,
    label = sprintf("mk=%.4f hgt=%.4f diff=%.4f 3se=%.4f",
                    p_mk, p_hgt, abs(p_mk - p_hgt), 3 * se))
})

test_that("HGT event counts are Poisson(lambda*k*dt) within 3-sigma MC tolerance", {
  # Star tree: 4 tips, single epoch, k=4, dt=1
  tree <- ape::read.tree(text = "(t1:1,t2:1,t3:1,t4:1);")
  Q <- matrix(c(-0.3, 0.3, 0.2, -0.2), nrow = 2, byrow = TRUE)
  rownames(Q) <- colnames(Q) <- c("0", "1")

  lambda <- 2.0; k <- 4L; dt <- 1.0
  expected <- lambda * k * dt  # = 8

  n_reps <- 2000
  counts <- vapply(seq_len(n_reps), function(i)
    simulate_hgt(tree, Q, lambda = lambda, seed = i)$n_hgt_events,
    integer(1L))

  obs <- mean(counts)
  se  <- sqrt(expected / n_reps)
  expect_true(abs(obs - expected) < 3 * se,
    label = sprintf("obs=%.3f expected=%.3f 3se=%.4f", obs, expected, 3 * se))
})

test_that("silent transfers are counted as events, not discarded", {
  # Q with zero rates: states never change via Mk, all HGT events are silent
  tree <- ape::read.tree(text = "(t1:1,t2:1,t3:1,t4:1);")
  Q_zero <- matrix(0, nrow = 2, ncol = 2)
  rownames(Q_zero) <- colnames(Q_zero) <- c("0", "1")

  lambda <- 3.0; k <- 4L; dt <- 1.0
  expected <- lambda * k * dt  # = 12

  n_reps <- 1000
  counts <- vapply(seq_len(n_reps), function(i)
    simulate_hgt(tree, Q_zero, lambda = lambda, root_freq = c(1, 0), seed = i)$n_hgt_events,
    integer(1L))

  obs <- mean(counts)
  se  <- sqrt(expected / n_reps)

  # Events fired even though all were silent
  expect_true(abs(obs - expected) < 3 * se,
    label = sprintf("obs=%.3f expected=%.3f 3se=%.4f", obs, expected, 3 * se))

  # And tip states are all 0 (silent = no state change)
  res <- simulate_hgt(tree, Q_zero, lambda = lambda, root_freq = c(1, 0), seed = 1)
  expect_true(all(res$tip_states == 0L))
})

test_that("high lambda drives tip states toward stationary distribution faster than lambda=0", {
  # Very short branches so Mk barely mixes; root all in state 0
  tips <- paste0("t", 1:10, ":0.01", collapse = ",")
  tree <- ape::read.tree(text = paste0("(", tips, ");"))
  Q <- matrix(c(-0.4, 0.4, 0.2, -0.2), nrow = 2, byrow = TRUE)
  rownames(Q) <- colnames(Q) <- c("0", "1")
  pi1 <- 0.4 / (0.4 + 0.2)  # stationary P(state=1) = 2/3

  n_reps <- 500
  mean_tip1 <- function(lam) {
    mean(vapply(seq_len(n_reps), function(i)
      mean(simulate_hgt(tree, Q, lambda = lam, root_freq = c(1, 0), seed = i)$tip_states),
      numeric(1L)))
  }
  m_low  <- mean_tip1(0)
  m_high <- mean_tip1(100)

  expect_true(abs(m_high - pi1) < abs(m_low - pi1),
    label = sprintf("low=%.4f high=%.4f stationary=%.4f", m_low, m_high, pi1))
})
