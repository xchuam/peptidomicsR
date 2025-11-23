#' Plot peptide length distribution by intensity or count
#'
#' @description
#' Summarizes and visualizes the distribution of peptide lengths,
#' either by total intensity or by peptide counts, stacking by protein group or name,
#' or overlaying sample-level density curves, with optional filtering, faceting,
#' and scientific notation on the y-axis.
#'
#' @param result List. Output of \code{processPeptides()}, containing at least:
#'   \itemize{
#'     \item \code{dt.peptides.int}: mean intensities,
#'     \item \code{dt.peptides.int.reps}: replicate-level intensities,
#'     \item \code{dt.peptides.count}: mean counts,
#'     \item \code{dt.peptides.count.reps}: replicate-level counts,
#'     \item \code{grp_cols}: character vector of grouping column names.
#'   }
#' @param type Character. Which table to plot: \code{"reps"} (replicate‐level) or \code{"mean"} (group‐mean). Default: \code{"mean"}.
#' @param metric Character. Which metric to plot: \code{"intensity"} for summed intensities,
#'   or \code{"count"} for peptide counts. Default: \code{"intensity"}.
#' @param color_by Character. Fill aesthetic: \code{"Protein.group"}, \code{"Protein.name"}, or \code{"none"}.
#'   Default: \code{"Protein.group"}. Ignored when \code{plot_mode = "density"}.
#' @param filter_params Named list, or \code{NULL}.  Each element’s name is a grouping column,
#'   and its value is a vector of values to include.  Multiple names impose an AND filter.
#'   For example: \code{list(Lipid = c("N","S"), Digest.stage = "G")}
#'   Default: \code{NULL} (no filtering).
#' @param facet_rows Character(1) or \code{NULL}.  Name of grouping column(s) for row facets.
#'   You can combine multiple variables with a plus: e.g.
#'   \code{"Casein.ratio+Digest.stage"}.  Default: \code{NULL}.
#' @param facet_cols Character(1) or \code{NULL}.  Name of grouping column(s) for column facets.
#'   Defaults to \code{"Replicate"} if \code{type = "reps"} and not explicitly set.
#'   You can combine multiple variables with a plus, e.g.
#'   \code{"Lipid+Replicate"}.  Default: \code{NULL}.
#' @param scientific_10_y Logical.  If \code{TRUE}, use scientific notation for y-axis.
#'   Default: \code{TRUE}. Ignored when \code{plot_mode = "density"}.
#' @param plot_mode Character. Choose between stacked \code{"bar"} (default) or overlaid
#'   \code{"density"} plots. In density mode, samples are distinguished by color in a single panel,
#'   ignoring \code{color_by}, \code{facet_rows}, and \code{facet_cols}.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' result <- processPeptides(
#'   peptides_file          = "data/Yogurtexample_QR188-205.csv",
#'   intensity_columns_file = "data/Intensity_columns.csv",
#'   protein_mapping_file   = "data/protein_mapping.csv"
#' )
#' # Plot mean Intensity vs. length, filter Lipid == "N", facet by Casein.ratio and Digest.stage
#' p1 <- plot_length_distribution(
#'   result,
#'   metric           = "intensity",
#'   filter_params    = list(Lipid = "N"),
#'   facet_rows       = "Casein.ratio",
#'   facet_cols       = "Digest.stage"
#' )
#' # Plot replicate‐level count vs. length,
#' colored by protein name, not scientific y-axis,
#' facet row by Casein.ratio+Digest.stage, col by Lipid+Replicate
#' p2 <- plot_length_distribution(
#'   result,
#'   type             = "reps",
#'   metric           = "count",
#'   color_by         = "Protein.name",
#'   facet_rows       = "Casein.ratio+Digest.stage",
#'   facet_cols       = "Lipid+Replicate",
#'   scientific_10_y  = FALSE
#' )
#'
#' # Plot replicate-level density curves colored by sample metadata
#' p3 <- plot_length_distribution(
#'   result,
#'   type       = "reps",
#'   metric     = "count",
#'   plot_mode  = "density"
#' )
#' }
#'
#' @import ggplot2
#' @import data.table
#' @importFrom ggpubr theme_pubr
#' @importFrom scales scientific_format hue_pal
#' @export
plot_length_distribution <- function(result,
                                     type             = c("mean", "reps"),
                                     metric           = c("intensity", "count"),
                                     color_by         = "Protein.group",
                                     filter_params    = NULL,
                                     facet_rows       = NULL,
                                     facet_cols       = NULL,
                                     scientific_10_y  = TRUE,
                                     plot_mode        = c("bar", "density")) {
  # validate type and metric
  type   <- match.arg(type)
  metric <- match.arg(metric)
  plot_mode <- match.arg(plot_mode)

  # select raw data based on type and metric
  if (type == "mean") {
    if (metric == "intensity") {
      dt  <- copy(result$dt.peptides.int)
      y_var <- "Mean.Intensity"
    } else {
      dt  <- copy(result$dt.peptides.count)
      y_var <- "Mean.Count"
    }
  } else {
    if (metric == "intensity") {
      dt  <- copy(result$dt.peptides.int.reps)
      y_var <- "Intensity"
    } else {
      dt  <- copy(result$dt.peptides.count.reps)
      y_var <- "Count"
    }
  }
  grp_cols <- result$grp_cols

  # default color_by validation for bar mode
  if (plot_mode == "bar") {
    if (!color_by %in% c("Protein.group","Protein.name","none")) {
      stop("color_by must be 'Protein.group', 'Protein.name', or 'none'.")
    }
  }

  # apply filters
  if (!is.null(filter_params)) {
    stopifnot(is.list(filter_params))
    for (col in names(filter_params)) {
      dt <- dt[get(col) %in% filter_params[[col]]]
    }
  }

  if (plot_mode == "bar") {
    # set default facet_cols for replicate-level
    if (type == "reps" &&
        (is.null(facet_rows) || !grepl("Replicate", facet_rows)) &&
        is.null(facet_cols)) {
      facet_cols <- "Replicate"
    }

    # build plot
    if (color_by == "none") {
      p <- ggplot(dt, aes(x = Length,
                          y = .data[[y_var]])) +
        geom_bar(stat = "identity", position = "stack",fill = def_color)+
        theme_pubr()
    } else {
      p <- ggplot(dt, aes(x = Length,
                          y = .data[[y_var]],
                          fill = .data[[color_by]])) +
        geom_bar(stat = "identity", position = "stack") +
        scale_fill_manual(values = protein_color)+
        theme_pubr()
    }

    # optional log scale
    if (scientific_10_y) {
      p <- p + scale_y_continuous(
        labels = scientific_10
      )
    }

    # apply faceting
    if (!is.null(facet_rows) || !is.null(facet_cols)) {
      rows <- if (!is.null(facet_rows)) facet_rows else "."
      cols <- if (!is.null(facet_cols)) facet_cols else "."
      p <- p + facet_grid(as.formula(paste(rows, "~", cols)))
    }

    # warn if any multi‐level group‐vars aren’t being distinguished anywher
    # build the full set of group‐vars we’d ideally split on
    vars_distinct <- grp_cols
    if (type == "reps") vars_distinct <- c(vars_distinct, "Replicate")
    # drop any that no longer vary in the data (only one unique value)
    vars_vary <- vars_distinct[
      sapply(vars_distinct, function(v) length(unique(dt[[v]])) > 1)]
    # collect what the user is actually splitting by
    used_vars <- c()
    if (!is.null(facet_rows)) used_vars <- c(used_vars, strsplit(facet_rows, "\\+")[[1]])
    if (!is.null(facet_cols)) used_vars <- c(used_vars, strsplit(facet_cols, "\\+")[[1]])
    used_vars <- unique(used_vars)
    # find which truly‐varying group‐vars are missing
    missing_vars <- setdiff(vars_vary, used_vars)
    if (length(missing_vars)) {
      warning(sprintf(
        "Grouping variable(s) %s not used in facets; data may be aggregated across them.",
        paste(missing_vars, collapse = ", ")
      ))
    }
  } else {
    # density mode ignores faceting and protein color mapping
    if (!is.null(facet_rows) || !is.null(facet_cols)) {
      warning("Density mode plots all samples in one panel; faceting arguments are ignored.")
    }
    sample_vars <- grp_cols
    if (type == "reps") sample_vars <- c(sample_vars, "Replicate")
    if (length(sample_vars) == 0) {
      dt[, Sample := "Sample"]
    } else {
      dt[, Sample := interaction(.SD, sep = "_", drop = TRUE),
         .SDcols = sample_vars]
    }
    sample_levels <- sort(unique(dt$Sample))
    sample_colors <- setNames(
      hue_pal()(length(sample_levels)),
      sample_levels
    )
    p <- ggplot(dt, aes(x = Length,
                        color = Sample)) +
      geom_density(aes(weight = .data[[y_var]]), linewidth = 1) +
      scale_color_manual(values = sample_colors) +
      labs(color = "Sample", y = "Density") +
      theme_pubr()
  }

  p
}
