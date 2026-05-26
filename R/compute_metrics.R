#' Compute per-replicate metrics from a run_one_rep() output
#'
#' @param rep_output list returned by run_one_rep()
#' @return one-row data.frame covering four metric families for all four models:
#'   model selection (aicc_winner, delta_aicc_*), parameter bias (*_q01_bias,
#'   *_q10_bias, *_rate_ratio_*), ASR Brier score (brier_*),
#'   stationary frequency error (*_pi0_*)
compute_metrics <- function(rep_output) {
  fits     <- rep_output$fits
  q_hats   <- rep_output$q_hat
  nprobs   <- rep_output$node_probs
  ns_true  <- rep_output$node_states_true

  q01_true <- fits$q01_true[1L]
  q10_true <- fits$q10_true[1L]
  rate_ratio_true <- if (q10_true > 0) q01_true / q10_true else NA_real_
  pi0_true        <- q10_true / (q01_true + q10_true)

  model_names <- c("ER", "ARD", "HMM2", "HMM3")

  per_model <- lapply(model_names, function(m) {
    qh  <- q_hats[[m]]
    np  <- nprobs[[m]]
    q01_hat <- unname(qh["q01"])
    q10_hat <- unname(qh["q10"])

    brier <- if (!is.null(np) && !is.null(ns_true) &&
                 is.matrix(np) && nrow(np) == length(ns_true)) {
      brier_internal(np, ns_true)
    } else {
      NA_real_
    }

    rate_ratio_hat <- if (!is.na(q01_hat) && !is.na(q10_hat) && q10_hat > 0) {
      q01_hat / q10_hat
    } else {
      NA_real_
    }

    pi0_hat <- if (!is.na(q01_hat) && !is.na(q10_hat) && (q01_hat + q10_hat) > 0) {
      q10_hat / (q01_hat + q10_hat)
    } else {
      NA_real_
    }

    list(
      q01_hat        = q01_hat,
      q10_hat        = q10_hat,
      q01_bias       = q01_hat - q01_true,
      q10_bias       = q10_hat - q10_true,
      rate_ratio_hat = rate_ratio_hat,
      pi0_hat        = pi0_hat,
      pi0_error      = pi0_hat - pi0_true,
      brier          = brier
    )
  })
  names(per_model) <- model_names

  # --- Model selection (shared across models) ---
  aicc_winner <- fits$aicc_winner[1L]
  aicc_ard    <- fits$aicc[fits$model == "ARD"]
  delta_aicc <- function(m) {
    a <- fits$aicc[fits$model == m]
    if (length(a) && is.finite(a) && is.finite(aicc_ard)) a - aicc_ard else NA_real_
  }

  flat <- c(
    list(
      aicc_winner      = aicc_winner,
      delta_aicc_ER    = delta_aicc("ER"),
      delta_aicc_HMM2  = delta_aicc("HMM2"),
      delta_aicc_HMM3  = delta_aicc("HMM3"),
      q01_true         = q01_true,
      q10_true         = q10_true,
      rate_ratio_true  = rate_ratio_true,
      pi0_true         = pi0_true
    ),
    unlist(lapply(model_names, function(m) {
      x <- per_model[[m]]
      stats::setNames(
        list(x$q01_hat, x$q10_hat, x$q01_bias, x$q10_bias,
             x$rate_ratio_hat, x$pi0_hat, x$pi0_error, x$brier),
        paste0(tolower(m), c("_q01_hat", "_q10_hat", "_q01_bias", "_q10_bias",
                             "_rate_ratio_hat", "_pi0_hat", "_pi0_error", "_brier"))
      )
    }), recursive = FALSE)
  )

  as.data.frame(flat, stringsAsFactors = FALSE, row.names = NULL)
}

# Mean multi-class Brier score across internal nodes
brier_internal <- function(node_probs, node_states_true) {
  n  <- nrow(node_probs)
  k  <- ncol(node_probs)
  bs <- numeric(n)
  for (i in seq_len(n)) {
    indicators        <- numeric(k)
    indicators[node_states_true[i] + 1L] <- 1  # 0-indexed -> 1-indexed
    bs[i] <- sum((node_probs[i, ] - indicators)^2)
  }
  mean(bs)
}
