#' Plot GRAVY score against log10 intensity
#'
#' @description
#' Visualize the relationship between peptide GRAVY score and log10-transformed intensity,
#' either at the mean (across replicates) or individual replicate level.
#'
#' @param result List. Output of \code{processPeptides()}, containing at least:
#'   \itemize{
#'     \item \code{dt.peptides.int}: mean intensities with GRAVY scores,
#'     \item \code{dt.peptides.int.reps}: replicate-level intensities with GRAVY scores,
#'     \item \code{grp_cols}: character vector of grouping column names.
#'   }
#' @param type Character. Which data to plot: \code{"mean"} for group means,
#'   or \code{"reps"} for replicate-level points.  Default: \code{"mean"}
#' @param color_by Character. Fill aesthetic: \code{"Protein.group"}, \code{"Protein.name"}, or \code{"none"}.
#'   Default: \code{"Protein.group"}.
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
#' @param alpha_value Numeric. Transparency level for points (0 = fully transparent, 1 = fully opaque).
#'   Default: \code{0.8}.
#'
#' @return A \code{ggplot} object.
#'
#' @import ggplot2
#' @import data.table
#' @importFrom ggpubr theme_pubr
#' @export
plot_gravy_vs_intensity <- function(result,
                                    type          = c("mean", "reps"),
                                    color_by      = "Protein.group",
                                    filter_params = NULL,
                                    facet_rows    = NULL,
                                    facet_cols    = NULL,
                                    alpha_var     = 0.8
) {
  # validate arguments
  type     <- match.arg(type)

  # select data
  if (type == "mean") {
    dt    <- copy(result$dt.peptides.int)
    y_var <- "log10.Mean.Intensity"
  } else {
    dt    <- copy(result$dt.peptides.int.reps)
    y_var <- "log10.Intensity"
  }
  grp_cols <- result$grp_cols

  # default color_by
  if (!color_by %in% c("Protein.group","Protein.name","none")) {
    stop("color_by must be 'Protein.group', 'Protein.name', or 'none'.")
  }

  # apply filters
  if (!is.null(filter_params)) {
    stopifnot(is.list(filter_params))
    for (col in names(filter_params)) {
      dt <- dt[get(col) %in% filter_params[[col]]]
    }
  }

  # set default facet_cols for replicate-level
  if (type == "reps" &&
      (is.null(facet_rows) || !grepl("Replicate", facet_rows)) &&
      is.null(facet_cols)) {
    facet_cols <- "Replicate"
  }

  # build base plot
  if (color_by == "none") {
    p <- ggplot(dt, aes(x = GRAVY.score,
                        y = .data[[y_var]])) +
      geom_point(color = def_color, alpha = alpha_var) +
      theme_pubr()
  } else {
    p <- ggplot(dt, aes(x = GRAVY.score,
                        y = .data[[y_var]],
                        color = .data[[color_by]])) +
      geom_point(alpha = alpha_var) +
      scale_color_manual(values = protein_color) +
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

  return(p)
}
