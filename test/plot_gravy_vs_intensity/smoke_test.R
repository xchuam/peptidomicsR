# Smoke test for plot_gravy_vs_intensity
# Run with: Rscript test/plot_gravy_vs_intensity/smoke_test.R

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

p_scatter <- plot_gravy_vs_intensity(result)
stopifnot(inherits(p_scatter, "ggplot"))

p_density <- plot_gravy_vs_intensity(
  result,
  type          = "reps",
  plot_mode     = "density",
  facet_rows    = "Digest.stage",
  filter_params = list(Yogurt = "Y1")
)
stopifnot(inherits(p_density, "ggplot"))
stopifnot(identical(p_density$labels$colour, "Replicate"))

cat("plot_gravy_vs_intensity smoke test completed\n")
