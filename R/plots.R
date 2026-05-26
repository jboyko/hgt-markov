#' Bias of q01 and q10 estimates vs. lambda
#'
#' @param tidy tidy data frame from run_sweep()
#' @return ggplot object
plot_bias_vs_lambda <- function(tidy) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")

  df <- tidy[, c("generator", "lambda", "q01_bias", "q10_bias")]
  long <- rbind(
    data.frame(generator = df$generator, lambda = df$lambda,
               param = "q01", bias = df$q01_bias, stringsAsFactors = FALSE),
    data.frame(generator = df$generator, lambda = df$lambda,
               param = "q10", bias = df$q10_bias, stringsAsFactors = FALSE)
  )

  # add tiny offset so lambda=0 is visible on log scale
  long$lambda_plot <- ifelse(long$lambda == 0, 1e-4, long$lambda)

  agg <- aggregate(bias ~ generator + lambda_plot + param, data = long,
                   FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                                       p10  = quantile(x, 0.10, na.rm = TRUE),
                                       p90  = quantile(x, 0.90, na.rm = TRUE)))
  agg <- do.call(data.frame, agg)
  names(agg) <- c("generator", "lambda_plot", "param", "mean", "p10", "p90")

  ggplot2::ggplot(agg, ggplot2::aes(
    x = lambda_plot, y = mean,
    colour = generator, fill = generator,
    ymin = p10, ymax = p90
  )) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_ribbon(alpha = 0.2, colour = NA) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::scale_x_log10(
      name = expression(lambda ~ "(events / lineage / unit time)")
    ) +
    ggplot2::scale_y_continuous(name = "Signed bias") +
    ggplot2::facet_wrap(~ param, labeller = ggplot2::label_parsed) +
    ggplot2::labs(colour = "Generator", fill = "Generator") +
    ggplot2::theme_bw()
}

#' Stacked-area P(AICc winner = X) vs. lambda
#'
#' @param tidy tidy data frame from run_sweep()
#' @return ggplot object
plot_model_selection <- function(tidy) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")

  df <- tidy[, c("generator", "lambda", "aicc_winner")]
  df$lambda_plot <- ifelse(df$lambda == 0, 1e-4, df$lambda)

  # P(winner = model) per (generator, lambda_plot)
  key <- paste(df$generator, df$lambda_plot, df$aicc_winner, sep = "|")
  n_per <- table(key)
  key_df <- do.call(rbind, strsplit(names(n_per), "|", fixed = TRUE))
  counts <- data.frame(
    generator  = key_df[, 1],
    lambda_plot = as.numeric(key_df[, 2]),
    aicc_winner = key_df[, 3],
    n           = as.integer(n_per),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  totals <- aggregate(n ~ generator + lambda_plot, data = counts, FUN = sum)
  names(totals)[3] <- "total"
  counts <- merge(counts, totals, by = c("generator", "lambda_plot"))
  counts$prob <- counts$n / counts$total

  ggplot2::ggplot(counts, ggplot2::aes(
    x = lambda_plot, y = prob,
    fill = aicc_winner
  )) +
    ggplot2::geom_area(position = "stack") +
    ggplot2::scale_x_log10(
      name = expression(lambda ~ "(events / lineage / unit time)")
    ) +
    ggplot2::scale_y_continuous(name = "P(AICc winner)", limits = c(0, 1)) +
    ggplot2::facet_wrap(~ generator) +
    ggplot2::labs(fill = "AICc winner") +
    ggplot2::theme_bw()
}

#' Mean ASR Brier score vs. lambda
#'
#' @param tidy tidy data frame from run_sweep()
#' @return ggplot object
plot_brier_vs_lambda <- function(tidy) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")

  df <- tidy[, c("generator", "lambda", "brier_asr")]
  df$lambda_plot <- ifelse(df$lambda == 0, 1e-4, df$lambda)

  agg <- aggregate(brier_asr ~ generator + lambda_plot, data = df,
                   FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                                       p10  = quantile(x, 0.10, na.rm = TRUE),
                                       p90  = quantile(x, 0.90, na.rm = TRUE)))
  agg <- do.call(data.frame, agg)
  names(agg) <- c("generator", "lambda_plot", "mean", "p10", "p90")

  ggplot2::ggplot(agg, ggplot2::aes(
    x = lambda_plot, y = mean,
    colour = generator, fill = generator,
    ymin = p10, ymax = p90
  )) +
    ggplot2::geom_ribbon(alpha = 0.2, colour = NA) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::scale_x_log10(
      name = expression(lambda ~ "(events / lineage / unit time)")
    ) +
    ggplot2::scale_y_continuous(name = "Mean Brier score") +
    ggplot2::labs(colour = "Generator", fill = "Generator") +
    ggplot2::theme_bw()
}

#' Save all three headline plots to out_dir as PDF and PNG
#'
#' @param tidy tidy data frame from run_sweep()
#' @param out_dir directory to write files into (created if absent)
#' @return invisibly, named list of ggplot objects
save_plots <- function(tidy, out_dir = "figures") {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  plots <- list(
    bias_vs_lambda    = plot_bias_vs_lambda(tidy),
    model_selection   = plot_model_selection(tidy),
    brier_vs_lambda   = plot_brier_vs_lambda(tidy)
  )

  for (nm in names(plots)) {
    p <- plots[[nm]]
    ggplot2::ggsave(file.path(out_dir, paste0(nm, ".pdf")), p,
                    width = 8, height = 5)
    ggplot2::ggsave(file.path(out_dir, paste0(nm, ".png")), p,
                    width = 8, height = 5, dpi = 150)
  }

  invisible(plots)
}
