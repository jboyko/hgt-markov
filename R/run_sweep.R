format_elapsed <- function(sec) {
  sprintf("%dm%02ds", as.integer(sec) %/% 60L, as.integer(sec) %% 60L)
}

#' Run a parameter sweep over a grid of generators and lambda values
#'
#' @param cells data.frame with columns `generator` (character) and `lambda` (numeric)
#' @param reps_per_cell integer; number of independent replicates per cell
#' @param n_cores integer; number of parallel workers (uses future.apply if >1)
#' @param cache_dir character; directory for per-rep RDS cache files
#' @param Q 2x2 rate matrix (rownames/colnames "0","1"); default mildly asymmetric ARD
#' @param n_target integer; number of tips for BD-sampled trees
#' @return tidy data.frame: one row per (cell, rep), metrics from compute_metrics()
#'   plus columns generator, lambda, rep
run_sweep <- function(cells, reps_per_cell = 5L, n_cores = 1L,
                      cache_dir = "cache",
                      Q = NULL, n_target = 100L) {
  if (is.null(Q)) {
    Q <- matrix(c(-0.1, 0.1, 0.3, -0.3), nrow = 2, byrow = TRUE)
    rownames(Q) <- colnames(Q) <- c("0", "1")
  }

  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

  jobs <- expand.grid(
    cell_idx = seq_len(nrow(cells)),
    rep      = seq_len(reps_per_cell),
    stringsAsFactors = FALSE
  )

  cache_path <- function(generator, lambda, rep) {
    lam_str <- formatC(lambda, format = "f", digits = 6, flag = "0")
    fname   <- paste0(generator, "__lam", lam_str, "__rep", rep, ".rds")
    file.path(cache_dir, fname)
  }

  n_jobs   <- nrow(jobs)
  n_cells  <- nrow(cells)
  t_start  <- proc.time()[["elapsed"]]

  run_job <- function(job_row, job_i = NA_integer_) {
    cell      <- cells[job_row$cell_idx, ]
    generator <- cell$generator
    lambda    <- cell$lambda
    rep       <- job_row$rep
    path      <- cache_path(generator, lambda, rep)

    cached <- file.exists(path)
    if (!cached) {
      seed <- as.integer((job_row$cell_idx * 1e4 + rep) %% .Machine$integer.max)
      rep_out <- run_one_rep(
        Q         = Q,
        generator = generator,
        lambda    = lambda,
        seed      = seed,
        n_target  = n_target
      )
      metrics <- compute_metrics(rep_out)
      row <- cbind(
        data.frame(generator = generator, lambda = lambda, rep = rep,
                   stringsAsFactors = FALSE),
        metrics
      )
      saveRDS(row, path)
    } else {
      row <- readRDS(path)
    }

    if (!is.na(job_i)) {
      elapsed  <- proc.time()[["elapsed"]] - t_start
      per_job  <- if (job_i > 1L) elapsed / (job_i - 1L) else NA_real_
      eta_sec  <- if (!is.na(per_job)) per_job * (n_jobs - job_i) else NA_real_
      eta_str  <- if (!is.na(eta_sec)) {
        sprintf("%dm%02ds", as.integer(eta_sec) %/% 60L, as.integer(eta_sec) %% 60L)
      } else "?"
      status   <- if (cached) "cached" else "done"
      message(sprintf("[%d/%d] %s | gen=%-8s lam=%.4g rep=%d | %s | ETA %s",
                      job_i, n_jobs, status,
                      generator, lambda, rep, format_elapsed(elapsed), eta_str))
    }

    row
  }

  if (n_cores > 1L) {
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      warning("future.apply not available; falling back to sequential execution")
      rows <- lapply(seq_len(n_jobs), function(i) run_job(jobs[i, ], i))
    } else {
      future::plan(future::multisession, workers = n_cores)
      on.exit(future::plan(future::sequential), add = TRUE)
      job_list <- lapply(seq_len(n_jobs), function(i) jobs[i, ])
      # progress not reliable across workers; run silently in parallel
      rows <- future.apply::future_lapply(job_list, run_job, future.seed = TRUE)
    }
  } else {
    rows <- lapply(seq_len(n_jobs), function(i) run_job(jobs[i, ], i))
  }

  do.call(rbind, c(rows, list(make.row.names = FALSE)))
}
