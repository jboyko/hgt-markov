library(ggplot2)

# Helper: minimal tidy sweep data frame covering all three generators
make_sweep_tidy <- function() {
  generators <- c("mk_null", "mk_hgt", "hmm_null")
  lambdas    <- c(0, 0.01, 0.1, 1)

  all_rows <- list()
  for (gen in generators) {
    lams <- if (gen == "mk_hgt") lambdas else 0
    for (lam in lams) {
      for (r in 1:3) {
        all_rows[[length(all_rows) + 1L]] <- data.frame(
          generator       = gen,
          lambda          = lam,
          rep             = r,
          q01_bias        = rnorm(1, 0, 0.05),
          q10_bias        = rnorm(1, 0, 0.05),
          q01_true        = 0.1,
          q10_true        = 0.3,
          q01_hat         = 0.1 + rnorm(1, 0, 0.05),
          q10_hat         = 0.3 + rnorm(1, 0, 0.05),
          aicc_winner     = sample(c("ARD", "ER", "HMM2", "HMM3"), 1),
          delta_aicc_ER   = rnorm(1, 2, 1),
          delta_aicc_HMM2 = rnorm(1, -1, 1),
          delta_aicc_HMM3 = rnorm(1,  1, 1),
          brier_asr       = runif(1, 0, 0.5),
          rate_ratio_hat  = runif(1, 0.2, 2),
          rate_ratio_true = 0.1 / 0.3,
          pi0_hat         = runif(1, 0.3, 0.7),
          pi0_true        = 0.3 / (0.1 + 0.3),
          pi0_error       = rnorm(1, 0, 0.05),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, c(all_rows, list(make.row.names = FALSE)))
}

tidy <- make_sweep_tidy()

# --- plot_bias_vs_lambda ---

test_that("plot_bias_vs_lambda returns a ggplot object", {
  p <- plot_bias_vs_lambda(tidy)
  expect_s3_class(p, "ggplot")
})

test_that("plot_bias_vs_lambda has log x-axis", {
  p <- plot_bias_vs_lambda(tidy)
  x_scales <- p$scales$scales
  is_log <- any(vapply(x_scales, function(s) {
    inherits(s, "ScaleContinuousPosition") &&
      !is.null(s$trans) &&
      grepl("log", s$trans$name, ignore.case = TRUE)
  }, logical(1)))
  expect_true(is_log)
})

test_that("plot_bias_vs_lambda renders without error", {
  p <- plot_bias_vs_lambda(tidy)
  expect_no_error(ggplot_build(p))
})

# --- plot_model_selection ---

test_that("plot_model_selection returns a ggplot object", {
  p <- plot_model_selection(tidy)
  expect_s3_class(p, "ggplot")
})

test_that("plot_model_selection renders without error", {
  p <- plot_model_selection(tidy)
  expect_no_error(ggplot_build(p))
})

# --- plot_brier_vs_lambda ---

test_that("plot_brier_vs_lambda returns a ggplot object", {
  p <- plot_brier_vs_lambda(tidy)
  expect_s3_class(p, "ggplot")
})

test_that("plot_brier_vs_lambda has log x-axis", {
  p <- plot_brier_vs_lambda(tidy)
  x_scales <- p$scales$scales
  is_log <- any(vapply(x_scales, function(s) {
    inherits(s, "ScaleContinuousPosition") &&
      !is.null(s$trans) &&
      grepl("log", s$trans$name, ignore.case = TRUE)
  }, logical(1)))
  expect_true(is_log)
})

test_that("plot_brier_vs_lambda renders without error", {
  p <- plot_brier_vs_lambda(tidy)
  expect_no_error(ggplot_build(p))
})

# --- save_plots ---

test_that("save_plots writes PDF and PNG files", {
  tmp <- tempdir()
  save_plots(tidy, out_dir = tmp)
  pdfs <- list.files(tmp, pattern = "\\.pdf$")
  pngs <- list.files(tmp, pattern = "\\.png$")
  expect_gte(length(pdfs), 3L)
  expect_gte(length(pngs), 3L)
})
