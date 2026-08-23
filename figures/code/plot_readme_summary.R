#!/usr/bin/env Rscript

required_packages <- c("ggplot2", "patchwork", "svglite", "ragg")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

command <- commandArgs(trailingOnly = FALSE)
script_argument <- grep("^--file=", command, value = TRUE)
if (length(script_argument) != 1L) {
  stop("Unable to resolve the plotting-script location.")
}
script_path <- normalizePath(sub("^--file=", "", script_argument))
repository_root <- normalizePath(file.path(dirname(script_path), "..", ".."))
source_path <- file.path(
  repository_root,
  "figures",
  "source_data",
  "readme_summary_metrics.csv"
)
output_stem <- file.path(repository_root, "figures", "readme_results_summary")

metrics <- read.csv(source_path, check.names = FALSE, stringsAsFactors = FALSE)
required_columns <- c("evidence", "category", "n_units", "AUROC", "AUPRC")
if (!identical(names(metrics), required_columns)) {
  stop("Unexpected source-data columns: ", paste(names(metrics), collapse = ", "))
}
if (any(!is.finite(metrics$AUROC)) || any(!is.finite(metrics$AUPRC))) {
  stop("All performance values must be finite.")
}
if (any(metrics$AUROC < 0 | metrics$AUROC > 1) ||
    any(metrics$AUPRC < 0 | metrics$AUPRC > 1)) {
  stop("Performance values must lie in [0, 1].")
}

to_long <- function(data) {
  rbind(
    transform(data, metric = "AUROC", value = AUROC),
    transform(data, metric = "AUPRC", value = AUPRC)
  )
}

method_order <- c("ConquerBN", "MeanProb", "IDR", "RRA", "RankProduct")
task_order <- c("TP vs TN", "TP vs FP", "FP vs TN")
metric_shapes <- c(AUROC = 16, AUPRC = 15)
blue <- "#0072B2"
grey <- "#9AA0A6"
dark_grey <- "#4A4F55"
task_colours <- c(
  "TP vs TN" = "#0072B2",
  "TP vs FP" = "#009E73",
  "FP vs TN" = "#D55E00"
)

theme_set(
  theme_classic(base_size = 7.5, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.text = element_text(size = 6.8, colour = "black"),
      axis.title = element_text(size = 7.2, colour = "black"),
      legend.title = element_blank(),
      legend.text = element_text(size = 6.5),
      plot.title = element_text(size = 8.0, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 6.5, colour = dark_grey, hjust = 0),
      plot.margin = margin(5, 12, 5, 5)
    )
)

make_benchmark_panel <- function(data, title, subtitle) {
  data$category <- factor(data$category, levels = rev(method_order))
  data$highlight <- ifelse(data$category == "ConquerBN", "ConquerBN", "Baseline")
  long <- to_long(data)

  ggplot(data, aes(y = category)) +
    geom_segment(
      aes(x = AUPRC, xend = AUROC, yend = category, colour = highlight),
      linewidth = 1.0,
      alpha = 0.65,
      show.legend = FALSE
    ) +
    geom_point(
      data = long,
      aes(x = value, shape = metric, colour = highlight),
      size = 2.2,
      stroke = 0.35,
      show.legend = c(shape = TRUE, colour = FALSE)
    ) +
    geom_text(
      data = subset(long, metric == "AUROC"),
      aes(x = value, label = sprintf("%.3f", value), colour = highlight),
      nudge_x = 0.016,
      nudge_y = -0.13,
      size = 2.0,
      hjust = 0,
      show.legend = FALSE
    ) +
    geom_text(
      data = subset(long, metric == "AUPRC"),
      aes(x = value, label = sprintf("%.3f", value), colour = highlight),
      nudge_x = 0.016,
      nudge_y = 0.13,
      size = 2.0,
      hjust = 0,
      show.legend = FALSE
    ) +
    scale_colour_manual(values = c(ConquerBN = blue, Baseline = grey)) +
    scale_shape_manual(values = metric_shapes) +
    scale_x_continuous(
      limits = c(0.25, 1.09),
      breaks = c(0.25, 0.50, 0.75, 1.00),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(title = title, subtitle = subtitle, x = "Macro performance", y = NULL) +
    theme(legend.position = "bottom")
}

make_rf_panel <- function(data) {
  data$category <- factor(data$category, levels = rev(task_order))
  long <- to_long(data)

  ggplot(data, aes(y = category, colour = category)) +
    geom_segment(
      aes(x = AUPRC, xend = AUROC, yend = category),
      linewidth = 1.0,
      alpha = 0.65,
      show.legend = FALSE
    ) +
    geom_point(
      data = long,
      aes(x = value, shape = metric),
      size = 2.2,
      stroke = 0.35,
      show.legend = c(shape = TRUE, colour = FALSE)
    ) +
    geom_text(
      data = subset(long, metric == "AUROC"),
      aes(x = value, label = sprintf("%.3f", value)),
      nudge_x = 0.016,
      nudge_y = -0.13,
      size = 2.0,
      hjust = 0,
      show.legend = FALSE
    ) +
    geom_text(
      data = subset(long, metric == "AUPRC"),
      aes(x = value, label = sprintf("%.3f", value)),
      nudge_x = 0.016,
      nudge_y = 0.13,
      size = 2.0,
      hjust = 0,
      show.legend = FALSE
    ) +
    scale_colour_manual(values = task_colours) +
    scale_shape_manual(values = metric_shapes) +
    scale_x_continuous(
      limits = c(0.25, 1.09),
      breaks = c(0.25, 0.50, 0.75, 1.00),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      title = "Direct-posterior RF screening",
      subtitle = "8 TP-vs-FP, 7 FP-vs-TN and 11 TP-vs-TN technologies",
      x = "Technology-macro performance",
      y = NULL
    ) +
    theme(legend.position = "bottom")
}

panel_a <- make_benchmark_panel(
  subset(metrics, evidence == "LOTO"),
  "Leave-one-technology-out",
  "14 held-out technologies; common five-method domain"
)
panel_b <- make_benchmark_panel(
  subset(metrics, evidence == "LOFO"),
  "Leave-one-family-out",
  "5 independent held-out technology families"
)
panel_c <- make_rf_panel(subset(metrics, evidence == "RF"))

figure <- ((panel_a | panel_b) / panel_c) +
  plot_layout(guides = "collect", heights = c(1.35, 0.9)) +
  plot_annotation(
    title = "ConquerBN reconstruction and posterior-derived error modelling",
    subtitle = "Total > 10; values are descriptive macro averages",
    tag_levels = "a",
    theme = theme(
      plot.title = element_text(size = 9.2, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 7.0, colour = dark_grey, hjust = 0),
      plot.tag = element_text(size = 9.0, face = "bold")
    )
  ) &
  theme(legend.position = "bottom")

width_mm <- 183
height_mm <- 135
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4

svglite::svglite(
  paste0(output_stem, ".svg"),
  width = width_in,
  height = height_in,
  bg = "white"
)
print(figure)
dev.off()

grDevices::pdf(
  paste0(output_stem, ".pdf"),
  width = width_in,
  height = height_in,
  family = "Helvetica",
  bg = "white",
  useDingbats = FALSE,
  compress = TRUE
)
print(figure)
dev.off()

grDevices::cairo_pdf(
  paste0(output_stem, "_cairo.pdf"),
  width = width_in,
  height = height_in,
  family = "sans",
  bg = "white"
)
print(figure)
dev.off()

ragg::agg_png(
  paste0(output_stem, ".png"),
  width = width_in,
  height = height_in,
  units = "in",
  res = 300,
  background = "white"
)
print(figure)
dev.off()

cat(
  "Saved summary figure to: ",
  output_stem,
  ".{svg,pdf,_cairo.pdf,png}\n",
  sep = ""
)
