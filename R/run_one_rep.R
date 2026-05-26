#' Run one simulation-and-inference replicate
#'
#' @param Q 2x2 rate matrix (rownames/colnames "0","1")
#' @param generator one of "mk_null", "mk_hgt", "hmm_null"
#' @param lambda HGT rate; used when generator="mk_hgt"
#' @param alpha hidden-category transition rate; used when generator="hmm_null"
#' @param Q2 2x2 rate matrix for hidden category 2; defaults to Q when generator="hmm_null"
#' @param seed integer random seed
#' @param tree phylo object; if NULL a BD-sampling tree is simulated via simulate_tree()
#' @param n_target tip count passed to simulate_tree() when tree is NULL
#' @return list with:
#'   $fits: data.frame with four rows (one per model: ER, ARD, HMM2, HMM3) and columns
#'     model, aicc, loglik, converged, conv_note, aicc_winner, q01_true, q10_true;
#'   $ard_q_hat: named numeric(2) c(q01, q10) from ARD fit, or NA if failed;
#'   $node_probs_ard: n_internal x 2 matrix of marginal posteriors from ARD fit, or NULL;
#'   $node_states_true: integer vector of true states at internal nodes (0-indexed)
run_one_rep <- function(Q, generator = "mk_null", lambda = 0, alpha = 0.5, Q2 = NULL,
                        seed = NULL, tree = NULL, n_target = 100L) {
  if (is.null(tree)) tree <- simulate_tree(n_target = n_target, seed = seed)
  sim <- switch(generator,
    mk_null  = simulate_mk(tree, Q, root_freq = c(0.5, 0.5), seed = seed),
    mk_hgt   = simulate_hgt(tree, Q, lambda = lambda, root_freq = c(0.5, 0.5), seed = seed),
    hmm_null = simulate_hmm(tree, Q1 = Q, Q2 = if (is.null(Q2)) Q else Q2,
                            alpha = alpha, root_freq = c(0.5, 0.5), seed = seed),
    stop("generator must be one of 'mk_null', 'mk_hgt', 'hmm_null'")
  )

  tip_data <- data.frame(
    taxon = names(sim$tip_states),
    state = factor(sim$tip_states, levels = c(0L, 1L)),
    stringsAsFactors = FALSE
  )

  model_specs <- list(
    ER   = list(rate.cat = 1L, model = "ER"),
    ARD  = list(rate.cat = 1L, model = "ARD"),
    HMM2 = list(rate.cat = 2L, model = "ARD"),
    HMM3 = list(rate.cat = 3L, model = "ARD")
  )

  fit_model <- function(spec) {
    tryCatch({
      suppressMessages(utils::capture.output(
        fit <- corHMM::corHMM(
          phy         = tree,
          data        = tip_data,
          rate.cat    = spec$rate.cat,
          model       = spec$model,
          node.states = "marginal",
          nstarts     = 1L,
          use_RTMB    = TRUE
        )
      ))
      aicc <- fit$AICc
      ll   <- fit$loglik
      ok   <- is.finite(aicc) && is.finite(ll)
      list(aicc      = if (ok) aicc else NA_real_,
           loglik    = if (ok) ll   else NA_real_,
           converged = ok,
           conv_note = if (ok) NA_character_ else "non-finite AICc/loglik",
           raw_fit   = if (ok) fit else NULL)
    }, error = function(e) {
      list(aicc = NA_real_, loglik = NA_real_,
           converged = FALSE, conv_note = conditionMessage(e), raw_fit = NULL)
    })
  }

  fits <- lapply(model_specs, fit_model)

  aiccs <- vapply(fits, `[[`, numeric(1), "aicc")
  finite_idx <- which(is.finite(aiccs))
  winner <- if (length(finite_idx) > 0L) {
    names(model_specs)[finite_idx[which.min(aiccs[finite_idx])]]
  } else {
    NA_character_
  }

  fits_df <- data.frame(
    model       = names(model_specs),
    aicc        = aiccs,
    loglik      = vapply(fits, `[[`, numeric(1), "loglik"),
    converged   = vapply(fits, `[[`, logical(1), "converged"),
    conv_note   = vapply(fits, `[[`, character(1), "conv_note"),
    aicc_winner = winner,
    q01_true    = Q[1, 2],
    q10_true    = Q[2, 1],
    stringsAsFactors = FALSE,
    row.names   = NULL
  )

  ard_raw <- fits[["ARD"]]$raw_fit
  ard_q_hat <- if (!is.null(ard_raw)) {
    c(q01 = ard_raw$solution["0", "1"], q10 = ard_raw$solution["1", "0"])
  } else {
    c(q01 = NA_real_, q10 = NA_real_)
  }
  node_probs_ard <- if (!is.null(ard_raw)) ard_raw$states else NULL

  list(
    fits             = fits_df,
    ard_q_hat        = ard_q_hat,
    node_probs_ard   = node_probs_ard,
    node_states_true = unname(sim$node_states)
  )
}
