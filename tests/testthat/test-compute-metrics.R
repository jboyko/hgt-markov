make_rep <- function(
  fits          = NULL,
  ard_q_hat     = c(q01 = 0.32, q10 = 0.19),
  node_probs_ard = matrix(c(1, 0, 0, 1, 1, 0), nrow = 3L, ncol = 2L, byrow = TRUE),
  node_states_true = c(0L, 1L, 0L)
) {
  if (is.null(fits)) {
    fits <- data.frame(
      model       = c("ER", "ARD", "HMM2", "HMM3"),
      aicc        = c(50.0, 45.0, 48.0, 52.0),
      loglik      = c(-22.0, -20.0, -21.0, -23.0),
      converged   = c(TRUE, TRUE, TRUE, TRUE),
      conv_note   = rep(NA_character_, 4L),
      aicc_winner = "ARD",
      q01_true    = 0.3,
      q10_true    = 0.2,
      stringsAsFactors = FALSE
    )
  }
  list(
    fits             = fits,
    ard_q_hat        = ard_q_hat,
    node_probs_ard   = node_probs_ard,
    node_states_true = node_states_true
  )
}

test_that("compute_metrics returns a one-row data frame", {
  result <- compute_metrics(make_rep())
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1L)
})

test_that("compute_metrics output has all four metric families", {
  result <- compute_metrics(make_rep())
  expect_true("aicc_winner"    %in% names(result))   # model selection
  expect_true("delta_aicc_ER"  %in% names(result))   # model selection delta
  expect_true("q01_bias"       %in% names(result))   # parameter bias
  expect_true("q10_bias"       %in% names(result))   # parameter bias
  expect_true("brier_asr"      %in% names(result))   # ASR Brier
  expect_true("pi0_error"      %in% names(result))   # stationary pi error
})

test_that("Brier score is 0 when posterior places prob 1 on true state at every node", {
  node_probs <- matrix(c(1, 0,
                         0, 1,
                         1, 0,
                         0, 1), nrow = 4L, ncol = 2L, byrow = TRUE)
  node_states_true <- c(0L, 1L, 0L, 1L)
  result <- compute_metrics(make_rep(node_probs_ard = node_probs,
                                     node_states_true = node_states_true))
  expect_equal(result$brier_asr, 0)
})

test_that("Brier score matches closed-form value when posterior is uniform (k=2 states)", {
  # uniform posterior: Brier = (0.5-1)^2 + (0.5-0)^2 = 0.5 regardless of true state
  n_nodes <- 5L
  node_probs <- matrix(0.5, nrow = n_nodes, ncol = 2L)
  node_states_true <- c(0L, 1L, 0L, 1L, 0L)
  result <- compute_metrics(make_rep(node_probs_ard = node_probs,
                                     node_states_true = node_states_true))
  expect_equal(result$brier_asr, 0.5)
})

test_that("Brier score is a mean across nodes, not a sum", {
  # 4 nodes, half perfect half uniform — mean ≠ sum
  node_probs <- matrix(c(1, 0,
                         1, 0,
                         0.5, 0.5,
                         0.5, 0.5), nrow = 4L, ncol = 2L, byrow = TRUE)
  node_states_true <- c(0L, 0L, 0L, 0L)
  result <- compute_metrics(make_rep(node_probs_ard = node_probs,
                                     node_states_true = node_states_true))
  # perfect nodes contribute 0, uniform contribute 0.5 each → mean = 0.25
  expect_equal(result$brier_asr, 0.25)
})

test_that("parameter bias columns are correct", {
  result <- compute_metrics(make_rep(ard_q_hat = c(q01 = 0.35, q10 = 0.15)))
  expect_equal(result$q01_bias, 0.35 - 0.3, tolerance = 1e-10)
  expect_equal(result$q10_bias, 0.15 - 0.2, tolerance = 1e-10)
})

test_that("delta_aicc columns are differences relative to ARD", {
  result <- compute_metrics(make_rep())
  expect_equal(result$delta_aicc_ER,   50.0 - 45.0)
  expect_equal(result$delta_aicc_HMM2, 48.0 - 45.0)
  expect_equal(result$delta_aicc_HMM3, 52.0 - 45.0)
})

test_that("stationary pi is computed from Q_hat (not tip frequencies)", {
  # q01 = q10 = 0.4 → pi0 = 0.4 / 0.8 = 0.5 regardless of tip states
  result <- compute_metrics(make_rep(ard_q_hat = c(q01 = 0.4, q10 = 0.4)))
  expect_equal(result$pi0_hat, 0.5, tolerance = 1e-10)
})

test_that("stationary pi_true computed from true Q, not estimated Q", {
  # q01_true=0.3, q10_true=0.2 → pi0_true = 0.2/(0.2+0.3) = 0.4
  result <- compute_metrics(make_rep())
  expect_equal(result$pi0_true, 0.2 / (0.2 + 0.3), tolerance = 1e-10)
})

test_that("pi0_error is pi0_hat - pi0_true", {
  result <- compute_metrics(make_rep(ard_q_hat = c(q01 = 0.4, q10 = 0.4)))
  expect_equal(result$pi0_error, result$pi0_hat - result$pi0_true, tolerance = 1e-10)
})

test_that("convergence failures propagate as NA metrics, not crashes", {
  failed_fits <- data.frame(
    model       = c("ER", "ARD", "HMM2", "HMM3"),
    aicc        = rep(NA_real_, 4L),
    loglik      = rep(NA_real_, 4L),
    converged   = rep(FALSE, 4L),
    conv_note   = rep("error", 4L),
    aicc_winner = NA_character_,
    q01_true    = 0.3,
    q10_true    = 0.2,
    stringsAsFactors = FALSE
  )
  rep_output <- list(
    fits             = failed_fits,
    ard_q_hat        = c(q01 = NA_real_, q10 = NA_real_),
    node_probs_ard   = NULL,
    node_states_true = c(0L, 1L, 0L)
  )
  result <- compute_metrics(rep_output)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1L)
  expect_true(is.na(result$q01_bias))
  expect_true(is.na(result$q10_bias))
  expect_true(is.na(result$brier_asr))
  expect_true(is.na(result$pi0_hat))
  expect_true(is.na(result$pi0_error))
  expect_true(is.na(result$delta_aicc_ER))
})
