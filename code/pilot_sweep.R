#!/usr/bin/env Rscript
# Pilot sweep for issue #9: ~20 reps x 3 generators x 10 lambda values
#
# Usage:
#   Rscript code/pilot_sweep.R [n_cores]          # run sweep (skips cached reps)
#   Rscript code/pilot_sweep.R --reload           # load cache/pilot_results.rds, replot only
#   Rscript code/pilot_sweep.R --bust-cache       # delete per-rep cache, rerun everything
#   Rscript code/pilot_sweep.R --bust-cache 4     # same, using 4 cores

args <- commandArgs(trailingOnly = TRUE)

reload     <- "--reload"     %in% args
bust_cache <- "--bust-cache" %in% args
core_arg   <- args[!args %in% c("--reload", "--bust-cache")]
n_cores    <- if (length(core_arg) >= 1L) as.integer(core_arg[1L]) else 1L

RESULTS_RDS <- "cache/pilot_results.rds"
CACHE_DIR   <- "cache/pilot"

# Source all pipeline functions from R/
r_files <- c("simulate.R", "run_one_rep.R", "compute_metrics.R",
             "run_sweep.R", "plots.R")
for (f in r_files) source(file.path("R", f))

if (reload) {
  if (!file.exists(RESULTS_RDS)) stop("No saved results at ", RESULTS_RDS,
                                      " — run without --reload first.")
  cat("Loading saved results from", RESULTS_RDS, "\n")
  tidy <- readRDS(RESULTS_RDS)
  cat(sprintf("Loaded %d rows\n", nrow(tidy)))
} else {
  if (bust_cache && dir.exists(CACHE_DIR)) {
    cat("Busting per-rep cache at", CACHE_DIR, "\n")
    unlink(list.files(CACHE_DIR, full.names = TRUE))
  }

  lambda_grid <- 10^seq(log10(0.001), log10(1), length.out = 10)

  cells <- expand.grid(
    generator = c("mk_null", "mk_hgt", "hmm_null"),
    lambda    = c(0, lambda_grid),
    stringsAsFactors = FALSE
  )
  row.names(cells) <- NULL

  cat(sprintf("Pilot sweep: %d cells, 20 reps each (%d total jobs)\n",
              nrow(cells), nrow(cells) * 20L))
  cat(sprintf("Running on %d core(s)\n", n_cores))

  tidy <- run_sweep(
    cells         = cells,
    reps_per_cell = 20L,
    n_cores       = n_cores,
    cache_dir     = CACHE_DIR,
    n_target      = 100L
  )

  dir.create(dirname(RESULTS_RDS), showWarnings = FALSE, recursive = TRUE)
  saveRDS(tidy, RESULTS_RDS)
  cat(sprintf("Saved %d rows to %s\n", nrow(tidy), RESULTS_RDS))
}

plots <- save_plots(tidy, out_dir = "figures/pilot")
cat("Plots written to figures/pilot/\n")
