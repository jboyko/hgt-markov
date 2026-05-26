#' Run one simulation-and-inference replicate
#'
#' @param Q 2x2 rate matrix (rownames/colnames "0","1")
#' @param lambda HGT rate (0 = pure Mk)
#' @param seed integer random seed
#' @param tree phylo object; if NULL a BD-sampling tree is simulated via simulate_tree()
#' @param n_target tip count passed to simulate_tree() when tree is NULL
#' @return one-row data.frame with q01_true, q10_true, q01_hat, q10_hat
run_one_rep <- function(Q, lambda = 0, seed = NULL, tree = NULL, n_target = 100L) {
  if (is.null(tree)) tree <- simulate_tree(n_target = n_target, seed = seed)
  sim <- if (lambda == 0) {
    simulate_mk(tree, Q, root_freq = c(0.5, 0.5), seed = seed)
  } else {
    simulate_hgt(tree, Q, lambda = lambda, root_freq = c(0.5, 0.5), seed = seed)
  }

  tip_data <- data.frame(
    taxon = names(sim$tip_states),
    state = factor(sim$tip_states, levels = c(0L, 1L)),
    stringsAsFactors = FALSE
  )

  fit <- suppressMessages(
    utils::capture.output(
      fit_out <- corHMM::corHMM(
        phy        = tree,
        data       = tip_data,
        rate.cat   = 1,
        model      = "ARD",
        node.states = "marginal",
        nstarts    = 1,
        use_RTMB   = TRUE
      )
    )
  )
  fit <- fit_out

  # corHMM ARD rate matrix: off-diagonals are [q01, q10]
  # rate.cat=1, 2 states: solution$solution is 2x2
  sol <- fit$solution
  q01_hat <- sol[1, 2]
  q10_hat <- sol[2, 1]

  data.frame(
    q01_true = Q[1, 2],
    q10_true = Q[2, 1],
    q01_hat  = q01_hat,
    q10_hat  = q10_hat
  )
}
