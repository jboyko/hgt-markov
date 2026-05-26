# Returns a closure that computes expm(Q*t), using eigen decomp with expm fallback.
make_expm_fn <- function(Q) {
  e     <- tryCatch(eigen(Q), error = function(err) NULL)
  V_inv <- if (!is.null(e)) tryCatch(solve(e$vectors), error = function(err) NULL) else NULL

  if (!is.null(e) && !is.null(V_inv) && !anyNA(e$values)) {
    V    <- e$vectors
    vals <- e$values
    function(t) {
      result <- tryCatch(
        Re(V %*% diag(exp(vals * t)) %*% V_inv),
        error = function(err) NULL
      )
      if (is.null(result) || anyNA(result)) expm::expm(Q * t)
      else pmax(result, 0)
    }
  } else {
    function(t) expm::expm(Q * t)
  }
}

#' Simulate binary character under Mk+HGT via forward-time epoch sweep
#'
#' @param tree phylo object (ape, time-calibrated)
#' @param Q 2x2 rate matrix with named rows/cols ("0","1")
#' @param lambda per-lineage HGT receiving rate (0 = pure Mk)
#' @param root_freq numeric(2) root state frequencies
#' @param seed integer random seed
#' @return list: tip_states, node_states (0-indexed integers), n_hgt_events
simulate_hgt <- function(tree, Q, lambda = 0, root_freq = c(0.5, 0.5), seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  diag(Q) <- 0
  diag(Q) <- -rowSums(Q)

  expm_t  <- make_expm_fn(Q)
  n_states <- nrow(Q)
  n_tips   <- length(tree$tip.label)
  n_nodes  <- tree$Nnode
  root_node <- n_tips + 1L

  # Node depths (distance from root); root = 0
  node_depths <- ape::node.depth.edgelength(tree)

  states         <- integer(n_tips + n_nodes)
  lineage_state  <- integer(nrow(tree$edge))  # current state per edge (1-indexed)
  states[root_node] <- sample.int(n_states, 1L, prob = root_freq / sum(root_freq))

  anc_nodes <- tree$edge[, 1L]
  des_nodes <- tree$edge[, 2L]
  t_breaks  <- sort(unique(node_depths))
  n_hgt_events <- 0L

  for (ep in seq_len(length(t_breaks) - 1L)) {
    t_start <- t_breaks[ep]
    t_end   <- t_breaks[ep + 1L]
    dt      <- t_end - t_start

    active <- which(node_depths[anc_nodes] <= t_start + 1e-10 &
                    node_depths[des_nodes] >= t_end   - 1e-10)
    k <- length(active)

    # Initialize lineages entering this epoch from their parent node state
    new_e <- active[abs(node_depths[anc_nodes[active]] - t_start) < 1e-10]
    lineage_state[new_e] <- states[anc_nodes[new_e]]

    # Forward HGT sweep within epoch
    remaining <- dt
    while (k >= 2L && lambda > 0 && remaining > 1e-14) {
      t_hgt <- rexp(1L, lambda * k)
      if (t_hgt >= remaining) break

      # Advance all lineages by t_hgt under Mk
      P <- expm_t(t_hgt)
      for (e in active) {
        lineage_state[e] <- sample.int(n_states, 1L, prob = P[lineage_state[e], ])
      }

      # Donor -> recipient state copy
      idx_r <- sample.int(k, 1L)
      idx_d <- sample.int(k - 1L, 1L)
      if (idx_d >= idx_r) idx_d <- idx_d + 1L
      lineage_state[active[idx_r]] <- lineage_state[active[idx_d]]

      n_hgt_events <- n_hgt_events + 1L
      remaining    <- remaining - t_hgt
    }

    # Advance remaining time under pure Mk
    if (remaining > 1e-14) {
      P_rem <- expm_t(remaining)
      for (e in active) {
        lineage_state[e] <- sample.int(n_states, 1L, prob = P_rem[lineage_state[e], ])
      }
    }

    # Record states for lineages terminating at this epoch
    ends_here <- active[abs(node_depths[des_nodes[active]] - t_end) < 1e-10]
    states[des_nodes[ends_here]] <- lineage_state[ends_here]
  }

  tip_states  <- states[seq_len(n_tips)] - 1L
  names(tip_states) <- tree$tip.label
  node_states <- states[root_node:(n_tips + n_nodes)] - 1L
  names(node_states) <- as.character(root_node:(n_tips + n_nodes))

  list(tip_states = tip_states, node_states = node_states, n_hgt_events = n_hgt_events)
}

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
