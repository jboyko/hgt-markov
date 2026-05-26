#' Run one simulation-and-inference replicate
#'
#' @param tree phylo object
#' @param Q 2x2 rate matrix (rownames/colnames "0","1")
#' @param lambda HGT rate (only 0 supported in this version)
#' @param seed integer random seed
#' @return one-row data.frame with q01_true, q10_true, q01_hat, q10_hat
run_one_rep <- function(tree, Q, lambda = 0, seed = NULL) {
  if (lambda != 0) stop("lambda != 0 not yet implemented")

  sim <- simulate_mk(tree, Q, root_freq = c(0.5, 0.5), seed = seed)

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
