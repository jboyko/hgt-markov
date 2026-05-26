#' Simulate a binary character on a phylogeny under the Mk model
#'
#' @param tree phylo object (ape)
#' @param Q 2x2 rate matrix with named rows/cols
#' @param root_freq numeric(2) root state frequencies
#' @param seed integer random seed
#' @return list with tip_states (named integer) and node_states (named integer)
simulate_mk <- function(tree, Q, root_freq = c(0.5, 0.5), seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  diag(Q) <- 0
  diag(Q) <- -rowSums(Q)

  n_states <- nrow(Q)
  tree <- ape::reorder.phylo(tree, "postorder")
  n_tips <- length(tree$tip.label)
  n_nodes <- tree$Nnode
  root_node <- n_tips + 1L

  states <- integer(n_tips + n_nodes)
  states[root_node] <- sample.int(n_states, 1L, prob = root_freq / sum(root_freq))

  anc <- tree$edge[, 1L]
  des <- tree$edge[, 2L]
  el  <- tree$edge.length
  n_edges <- nrow(tree$edge)

  for (i in n_edges:1L) {
    P <- expm::expm(Q * el[i])
    states[des[i]] <- sample.int(n_states, 1L, prob = P[states[anc[i]], ])
  }

  tip_states <- states[seq_len(n_tips)] - 1L  # 0-indexed states to match Q rownames "0","1"
  names(tip_states) <- tree$tip.label
  node_states <- states[root_node:(n_tips + n_nodes)] - 1L
  names(node_states) <- as.character(root_node:(n_tips + n_nodes))

  list(tip_states = tip_states, node_states = node_states)
}
