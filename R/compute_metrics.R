#' Compute per-replicate metrics from a run_one_rep() output
#'
#' @param rep_output list returned by run_one_rep(); must contain:
#'   $fits (4-row data.frame), $ard_q_hat (named numeric(2) with q01/q10),
#'   $node_probs_ard (n_internal x 2 matrix or NULL), $node_states_true (integer vector)
#' @return one-row data.frame covering four metric families: model selection
#'   (aicc_winner, delta_aicc_*), parameter bias (q01_bias, q10_bias, rate_ratio_*),
#'   ASR Brier score (brier_asr), stationary frequency error (pi0_*)
compute_metrics <- function(rep_output) {
  fits   <- rep_output$fits
  q_hat  <- rep_output$ard_q_hat
  nprobs <- rep_output$node_probs_ard
  ns_true <- rep_output$node_states_true

  q01_true <- fits$q01_true[1L]
  q10_true <- fits$q10_true[1L]

  # --- 1. Parameter bias (matched ARD fit) ---
  q01_hat <- unname(if (!is.null(q_hat) && !is.na(q_hat["q01"])) q_hat["q01"] else NA_real_)
  q10_hat <- unname(if (!is.null(q_hat) && !is.na(q_hat["q10"])) q_hat["q10"] else NA_real_)
  q01_bias <- q01_hat - q01_true
  q10_bias <- q10_hat - q10_true
  rate_ratio_hat  <- if (!is.na(q01_hat) && !is.na(q10_hat) && q10_hat > 0) {
    q01_hat / q10_hat
  } else {
    NA_real_
  }
  rate_ratio_true <- if (q10_true > 0) q01_true / q10_true else NA_real_

  # --- 2. Model-selection ---
  aicc_winner <- fits$aicc_winner[1L]
  aicc_ard  <- fits$aicc[fits$model == "ARD"]
  aicc_er   <- fits$aicc[fits$model == "ER"]
  aicc_hmm2 <- fits$aicc[fits$model == "HMM2"]
  aicc_hmm3 <- fits$aicc[fits$model == "HMM3"]
  delta_aicc <- function(a, ref) {
    if (length(a) && length(ref) && is.finite(a) && is.finite(ref)) a - ref else NA_real_
  }
  delta_aicc_ER   <- delta_aicc(aicc_er,   aicc_ard)
  delta_aicc_HMM2 <- delta_aicc(aicc_hmm2, aicc_ard)
  delta_aicc_HMM3 <- delta_aicc(aicc_hmm3, aicc_ard)

  # --- 3. ASR Brier score ---
  brier_asr <- if (!is.null(nprobs) && !is.null(ns_true) &&
                   is.matrix(nprobs) && nrow(nprobs) == length(ns_true)) {
    brier_internal(nprobs, ns_true)
  } else {
    NA_real_
  }

  # --- 4. Stationary frequency from Q ---
  pi0_hat <- if (!is.na(q01_hat) && !is.na(q10_hat) && (q01_hat + q10_hat) > 0) {
    q10_hat / (q01_hat + q10_hat)
  } else {
    NA_real_
  }
  pi0_true  <- q10_true / (q01_true + q10_true)
  pi0_error <- pi0_hat - pi0_true

  data.frame(
    aicc_winner     = aicc_winner,
    delta_aicc_ER   = delta_aicc_ER,
    delta_aicc_HMM2 = delta_aicc_HMM2,
    delta_aicc_HMM3 = delta_aicc_HMM3,
    q01_hat         = q01_hat,
    q10_hat         = q10_hat,
    q01_true        = q01_true,
    q10_true        = q10_true,
    q01_bias        = q01_bias,
    q10_bias        = q10_bias,
    rate_ratio_hat  = rate_ratio_hat,
    rate_ratio_true = rate_ratio_true,
    brier_asr       = brier_asr,
    pi0_hat         = pi0_hat,
    pi0_true        = pi0_true,
    pi0_error       = pi0_error,
    stringsAsFactors = FALSE,
    row.names        = NULL
  )
}

# Mean multi-class Brier score across internal nodes
brier_internal <- function(node_probs, node_states_true) {
  n  <- nrow(node_probs)
  k  <- ncol(node_probs)
  bs <- numeric(n)
  for (i in seq_len(n)) {
    indicators        <- numeric(k)
    indicators[node_states_true[i] + 1L] <- 1  # 0-indexed → 1-indexed
    bs[i] <- sum((node_probs[i, ] - indicators)^2)
  }
  mean(bs)
}
