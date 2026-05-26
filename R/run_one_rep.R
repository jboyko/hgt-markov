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
#'   $node_probs: named list of n_internal x 2 marginal posterior matrices, one per model
#'     (NULL for a model that failed); for HMM models the matrix is collapsed to 2 columns
#'     by summing over hidden-rate categories;
#'   $q_hat: named list of c(q01, q10) per model; HMM values are means across categories;
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

  node_probs <- lapply(fits, function(f) collapse_node_probs(f$raw_fit))
  q_hat      <- lapply(fits, function(f) extract_q_hat(f$raw_fit))

  list(
    fits             = fits_df,
    node_probs       = node_probs,
    q_hat            = q_hat,
    node_states_true = unname(sim$node_states)
  )
}

# Collapse HMM node-prob matrix (n x 2k) to (n x 2) by summing over hidden cats.
# Returns NULL if raw_fit is NULL.
collapse_node_probs <- function(raw_fit) {
  if (is.null(raw_fit)) return(NULL)
  mat <- raw_fit$states
  if (is.null(mat)) return(NULL)
  k <- ncol(mat) / 2L
  if (k == 1L) return(mat)
  # columns alternate: state0_cat1, state1_cat1, state0_cat2, state1_cat2, ...
  state0_cols <- seq(1L, ncol(mat), by = 2L)
  state1_cols <- seq(2L, ncol(mat), by = 2L)
  cbind(rowSums(mat[, state0_cols, drop = FALSE]),
        rowSums(mat[, state1_cols, drop = FALSE]))
}

# Extract (q01, q10) from any corHMM fit. For HMM models, averages across
# rate categories. Returns c(q01=NA, q10=NA) if raw_fit is NULL.
extract_q_hat <- function(raw_fit) {
  if (is.null(raw_fit)) return(c(q01 = NA_real_, q10 = NA_real_))
  sol <- raw_fit$solution
  rn  <- rownames(sol)
  zero_idx <- grep("^0", rn)
  one_idx  <- grep("^1", rn)
  q01_vals <- mapply(function(z, o) sol[z, o], zero_idx, one_idx)
  q10_vals <- mapply(function(z, o) sol[o, z], zero_idx, one_idx)
  c(q01 = mean(q01_vals, na.rm = TRUE),
    q10 = mean(q10_vals, na.rm = TRUE))
}
