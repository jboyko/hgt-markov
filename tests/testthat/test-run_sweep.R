Q_default <- function() {
  Q <- matrix(c(-0.1, 0.1, 0.3, -0.3), nrow = 2, byrow = TRUE)
  rownames(Q) <- colnames(Q) <- c("0", "1")
  Q
}

small_cells <- function() {
  data.frame(
    generator = c("mk_null", "mk_hgt", "hmm_null"),
    lambda    = c(0, 0.5, 0),
    stringsAsFactors = FALSE
  )
}

test_that("run_sweep returns a data frame with expected columns", {
  skip_if_not_installed("corHMM")
  td <- withr::local_tempdir()

  result <- run_sweep(
    cells        = small_cells(),
    reps_per_cell = 1L,
    n_cores      = 1L,
    cache_dir    = td,
    Q            = Q_default(),
    n_target     = 10L
  )

  expect_s3_class(result, "data.frame")
  expect_true(all(c("generator", "lambda", "rep", "aicc_winner") %in% names(result)))
})

test_that("run_sweep returns one row per (cell, rep)", {
  skip_if_not_installed("corHMM")
  td <- withr::local_tempdir()

  cells <- small_cells()
  reps  <- 2L

  result <- run_sweep(
    cells        = cells,
    reps_per_cell = reps,
    n_cores      = 1L,
    cache_dir    = td,
    Q            = Q_default(),
    n_target     = 10L
  )

  expect_equal(nrow(result), nrow(cells) * reps)
})

test_that("run_sweep writes one RDS cache file per (cell, rep)", {
  skip_if_not_installed("corHMM")
  td <- withr::local_tempdir()

  cells <- small_cells()
  reps  <- 2L

  run_sweep(
    cells        = cells,
    reps_per_cell = reps,
    n_cores      = 1L,
    cache_dir    = td,
    Q            = Q_default(),
    n_target     = 10L
  )

  rds_files <- list.files(td, pattern = "\\.rds$", recursive = TRUE)
  expect_equal(length(rds_files), nrow(cells) * reps)
})

test_that("second run_sweep call makes zero new fits (all cache hits)", {
  skip_if_not_installed("corHMM")
  td <- withr::local_tempdir()

  args <- list(
    cells        = small_cells(),
    reps_per_cell = 1L,
    n_cores      = 1L,
    cache_dir    = td,
    Q            = Q_default(),
    n_target     = 10L
  )

  do.call(run_sweep, args)
  files_after_first <- list.files(td, pattern = "\\.rds$")

  # Touch mtimes so we can detect any re-writes
  Sys.sleep(0.05)
  mtime_before <- file.info(file.path(td, files_after_first))$mtime

  do.call(run_sweep, args)
  mtime_after <- file.info(file.path(td, files_after_first))$mtime

  expect_equal(mtime_before, mtime_after)
})

test_that("deleting one cache file causes only that rep to re-run", {
  skip_if_not_installed("corHMM")
  td <- withr::local_tempdir()

  args <- list(
    cells        = small_cells(),
    reps_per_cell = 2L,
    n_cores      = 1L,
    cache_dir    = td,
    Q            = Q_default(),
    n_target     = 10L
  )

  do.call(run_sweep, args)
  all_files <- list.files(td, pattern = "\\.rds$")
  expect_equal(length(all_files), nrow(small_cells()) * 2L)

  Sys.sleep(0.05)
  target <- all_files[1]
  file.remove(file.path(td, target))

  mtime_before <- file.info(file.path(td, all_files[-1]))$mtime

  do.call(run_sweep, args)

  mtime_after <- file.info(file.path(td, all_files[-1]))$mtime
  expect_equal(mtime_before, mtime_after)

  expect_true(file.exists(file.path(td, target)))
})

test_that("smoke: 3 generators x 3 lambdas x 2 reps produces correct row count", {
  skip_if_not_installed("corHMM")
  skip_on_cran()
  td <- withr::local_tempdir()

  cells <- data.frame(
    generator = c("mk_null", "mk_hgt", "mk_hgt", "mk_hgt", "hmm_null"),
    lambda    = c(0, 0.1, 1.0, 5.0, 0),
    stringsAsFactors = FALSE
  )

  result <- run_sweep(
    cells        = cells,
    reps_per_cell = 2L,
    n_cores      = 1L,
    cache_dir    = td,
    Q            = Q_default(),
    n_target     = 15L
  )

  expect_equal(nrow(result), nrow(cells) * 2L)
  expect_true(all(c("generator", "lambda", "rep", "aicc_winner",
                    "q01_bias", "q10_bias", "brier_asr") %in% names(result)))
})
