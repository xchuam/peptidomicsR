#' Plot GRAVY score against log10 intensity
#'
#' @description
#' Visualize the relationship between peptide GRAVY score and log10-transformed intensity,
#' either at the mean (across replicates) or individual replicate level.
#' Alternatively, display density curves of GRAVY score distributions weighted by intensity.
#'
#' @param result List. Output of \code{processPeptides()}, containing at least:
#'   \itemize{
#'     \item \code{dt.peptides.int}: mean intensities with GRAVY scores,
#'     \item \code{dt.peptides.int.reps}: replicate-level intensities with GRAVY scores,
#'     \item \code{grp_cols}: character vector of grouping column names.
#'   }
#' @param type Character. Which data to plot: \code{"mean"} for group means,
#'   or \code{"reps"} for replicate-level points.  Default: \code{"mean"}
#' @param color_by Character. Color aesthetic for scatter plots: \code{"Protein.group"}, \code{"Protein.name"}, or \code{"none"}.
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
#' @param alpha_value Numeric. Transparency level for scatter plot points (0 = fully transparent, 1 = fully opaque).
#'   Default: \code{0.8}. Ignored when \code{plot_mode = "density"}.
#' @param plot_mode Character. Choose between point \code{"scatter"} (default) or GRAVY score
#'   \code{"density"} plots. Density mode colors lines by the minimal set of grouping variables
#'   (excluding facets) that still vary after filtering, and faceting remains available.
#'
#' @examples
#' \dontrun{
#' result <- processPeptides(
#'   peptides_file          = "data/Yogurtexample_QR188-205.csv",
#'   intensity_columns_file = "data/Intensity_columns.csv",
#'   protein_mapping_file   = "data/protein_mapping.csv"
#' )
#' # Scatter plot of GRAVY score vs log10 mean intensity
#' plot_gravy_vs_intensity(result)
#'
#' # Density plot per replicate, faceted by Digest.stage and filtered to Yogurt Y1
#' plot_gravy_vs_intensity(
#'   result,
#'   type          = "reps",
#'   plot_mode     = "density",
#'   facet_rows    = "Digest.stage",
#'   filter_params = list(Yogurt = "Y1")
#' )
#' }
#'
#' @return A \code{ggplot} object.
#'
#' @import ggplot2
#' @import data.table
#' @importFrom ggpubr theme_pubr
#' @importFrom scales hue_pal
#' @export
plot_gravy_vs_intensity <- function(result,
                                    type          = c("mean", "reps"),
                                    color_by      = "Protein.group",
                                    filter_params = NULL,
                                    facet_rows    = NULL,
                                    facet_cols    = NULL,
                                    alpha_value   = 0.8,
                                    plot_mode     = c("scatter", "density")
) {
  # validate arguments
  type     <- match.arg(type)
  plot_mode <- match.arg(plot_mode)

  # select data
  if (type == "mean") {
    dt    <- copy(result$dt.peptides.int)
    y_var <- "log10.Mean.Intensity"
    weight_var <- "Mean.Intensity"
  } else {
    dt    <- copy(result$dt.peptides.int.reps)
    y_var <- "log10.Intensity"
    weight_var <- "Intensity"
  }
  grp_cols <- result$grp_cols

  # default color_by
  if (plot_mode == "scatter") {
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

  # set default facet_cols for replicate-level
  if (plot_mode == "scatter" &&
      type == "reps" &&
      (is.null(facet_rows) || !grepl("Replicate", facet_rows)) &&
      is.null(facet_cols)) {
    facet_cols <- "Replicate"
  }

  split_facets <- function(x) {
    if (is.null(x) || identical(x, ".")) return(character())
    vars <- unique(unlist(strsplit(x, "\\+")))
    vars <- trimws(vars)
    vars[nzchar(vars)]
  }
  facet_vars <- unique(c(split_facets(facet_rows), split_facets(facet_cols)))

  dt_plot <- dt
  grouping_colors <- character(0)

  if (plot_mode == "scatter") {
    if (color_by == "none") {
      p <- ggplot(dt, aes(x = GRAVY.score,
                          y = .data[[y_var]])) +
        geom_point(color = def_color, alpha = alpha_value) +
        theme_pubr()
    } else {
      p <- ggplot(dt, aes(x = GRAVY.score,
                          y = .data[[y_var]],
                          color = .data[[color_by]])) +
        geom_point(alpha = alpha_value) +
        scale_color_manual(values = protein_color) +
        theme_pubr()
    }
  } else {
    dt_plot <- dt[!is.na(GRAVY.score) & !is.na(get(weight_var))]
    sample_vars <- grp_cols
    if (type == "reps") sample_vars <- c(sample_vars, "Replicate")
    if (length(sample_vars)) {
      sample_vars <- sample_vars[
        vapply(sample_vars, function(v) length(unique(dt_plot[[v]])) > 1, logical(1))
      ]
    }
    sample_vars <- setdiff(sample_vars, facet_vars)
    if (length(sample_vars) == 0) {
      dt_plot[, Sample := "Sample"]
      legend_title <- "Sample"
    } else {
      dt_plot[, Sample := interaction(.SD, sep = "_", drop = TRUE),
              .SDcols = sample_vars]
      legend_title <- paste(sample_vars, collapse = " × ")
    }
    sample_levels <- sort(unique(dt_plot$Sample))
    if (length(sample_levels) == 0) {
      sample_colors <- setNames(def_color, "Sample")
    } else {
      palette_vals <- hue_pal()(max(1, length(sample_levels)))
      sample_colors <- setNames(palette_vals[seq_along(sample_levels)], sample_levels)
    }
    grouping_colors <- sample_vars
    p <- ggplot(dt_plot, aes(x = GRAVY.score,
                             color = Sample)) +
      geom_density(aes(weight = .data[[weight_var]]), linewidth = 1) +
      scale_color_manual(values = sample_colors) +
      labs(color = legend_title, y = "Density") +
      theme_pubr()
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
    sapply(vars_distinct, function(v) length(unique(dt_plot[[v]])) > 1)]
  # collect what the user is actually splitting/colouring by
  used_vars <- unique(c(facet_vars, grouping_colors))
  # find which truly‐varying group‐vars are missing
  missing_vars <- setdiff(vars_vary, used_vars)
  if (length(missing_vars)) {
    warning(sprintf(
      "Grouping variable(s) %s not used in facets; data may be aggregated across them.",
      paste(missing_vars, collapse = ", ")
    ))
  }

  return(p)
}
