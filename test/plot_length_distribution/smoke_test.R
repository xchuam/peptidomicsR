# Smoke test for plot_length_distribution
# Run with: Rscript test/plot_length_distribution/smoke_test.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggpubr)
  library(scales)
})

r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
invisible(lapply(r_files, source))

result <- processPeptides(
  peptides_file          = "data/Yogurtexample_QR188-205.csv",
  intensity_columns_file = "data/Intensity_columns.csv",
  protein_mapping_file   = "data/protein_mapping.csv"
)

p_bar <- plot_length_distribution(
  result,
  metric      = "count",
  facet_rows  = "Digest.stage",
  facet_cols  = "Yogurt"
)
stopifnot(inherits(p_bar, "ggplot"))

p_density <- plot_length_distribution(
  result,
  type      = "reps",
  metric    = "count",
  plot_mode = "density",
  facet_rows = "Digest.stage"
)
stopifnot(inherits(p_density, "ggplot"))
stopifnot(identical(p_density$labels$colour, "Yogurt × Replicate"))

cat("plot_length_distribution smoke test completed\n")
