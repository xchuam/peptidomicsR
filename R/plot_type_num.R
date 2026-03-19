#' Plot peptide type numbers with optional scientific notation y-axis, flexible filtering and faceting
#'
#' @description
#' Generates stacked-bar charts of numbers of unique peptide types (per-replicate or group-mean)
#' with:
#' * User-defined inclusion filters on any grouping columns,
#' * Dynamic x-axis chosen from grouping variables,
#' * Flexible coloring by protein group or protein name,
#' * Optional scientific notation of the y-axis,
#' * Flexible faceting by one or more grouping variables.
#'
#' @param result List. Output of \code{processPeptides()}, containing at least:
#'   \itemize{
#'     \item \code{dt.peptides.typenum.reps}: replicate-level peptide type numbers,
#'     \item \code{dt.peptides.typenum}: group-mean peptide type numbers,
#'     \item \code{grp_cols}: character vector of grouping column names.
#'   }
#' @param type Character. Which table to plot: \code{"reps"} (replicate-level) or \code{"mean"} (group-mean). Default: \code{"mean"}.
#' @param x_var Character. Name of the grouping column to use on the x-axis.
#'   Must be one of \code{result$grp_cols} (and for \code{type = "reps"} may also be \code{"Replicate"}).
#'   Defaults to the first element of \code{result$grp_cols}.
#' @param color_by Character. Fill aesthetic: \code{"Protein.group"}, \code{"Protein.name"}, or \code{"none"}.
#'   Default: \code{"Protein.group"}.
#' @param filter_params Named list, or \code{NULL}. Each element's name is a grouping column,
#'   and its value is a vector of values to include. Multiple names impose an AND filter.
#'   For example: \code{list(Yogurt = c("Y1", "Y2"), Digest.stage = "G120")}
#'   Default: \code{NULL} (no filtering).
#' @param facet_rows Character(1) or \code{NULL}. Name of grouping column(s) for row facets.
#'   You can combine multiple variables with a plus: e.g. \code{"Yogurt+Digest.stage"}.
#'   Default: \code{NULL}.
#' @param facet_cols Character(1) or \code{NULL}. Name of grouping column(s) for column facets.
#'   Defaults to \code{"Replicate"} if \code{type = "reps"} and not explicitly set.
#'   You can combine multiple variables with a plus, e.g. \code{"Yogurt+Replicate"}.
#'   Default: \code{NULL}.
#' @param scientific_10_y Logical. If \code{TRUE}, use scientific notation for y-axis.
#'   Default: \code{FALSE}.
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
#'
#' # 1) Simple replicate-level, default x, default facet = "Replicate"
#' p1 <- plot_type_num(result, type = "reps")
#'
#' # 2) Simple mean-level, x = Yogurt, color by Protein.name
#' p2 <- plot_type_num(
#'   result,
#'   type     = "mean",
#'   x_var    = "Yogurt",
#'   color_by = "Protein.name"
#' )
#'
#' # 3) Filter Yogurt == "Y1" AND Digest.stage == "I120"
#' p3 <- plot_type_num(
#'   result,
#'   filter_params   = list(Yogurt = "Y1", Digest.stage = "I120"),
#'   scientific_10_y = TRUE
#' )
#' }
#'
#' @import ggplot2
#' @import data.table
#' @importFrom ggpubr theme_pubr
#' @importFrom scales scientific_format
#' @export
plot_type_num <- function(result,
                          type             = c("mean", "reps"),
                          x_var            = NULL,
                          color_by         = "Protein.group",
                          filter_params    = NULL,
                          facet_rows       = NULL,
                          facet_cols       = NULL,
                          scientific_10_y  = FALSE) {
  .plot_type_num_impl(
    result          = result,
    type            = type,
    x_var           = x_var,
    color_by        = color_by,
    filter_params   = filter_params,
    facet_rows      = facet_rows,
    facet_cols      = facet_cols,
    scientific_10_y = scientific_10_y
  )
}

.plot_type_num_impl <- function(result,
                                type             = c("mean", "reps"),
                                x_var            = NULL,
                                color_by         = "Protein.group",
                                filter_params    = NULL,
                                facet_rows       = NULL,
                                facet_cols       = NULL,
                                scientific_10_y  = FALSE) {
  type <- match.arg(type)

  dt <- switch(type,
               reps = copy(result$dt.peptides.typenum.reps),
               mean = copy(result$dt.peptides.typenum))
  grp_cols <- result$grp_cols

  if (is.null(x_var)) {
    x_var <- grp_cols[1]
  }

  if (type == "reps") {
    if (!x_var %in% c("Replicate", grp_cols)) {
      stop("For type = 'reps', x_var must be 'Replicate' or one of result$grp_cols")
    }
  } else {
    if (!x_var %in% grp_cols) {
      stop("For type = 'mean', x_var must be one of result$grp_cols")
    }
  }

  if (!color_by %in% c("Protein.group", "Protein.name", "none")) {
    stop("color_by must be 'Protein.group', 'Protein.name', or 'none'.")
  }

  if (!is.null(filter_params)) {
    stopifnot(is.list(filter_params))
    for (col in names(filter_params)) {
      dt <- dt[get(col) %in% filter_params[[col]]]
    }
  }

  if (type == "reps" &&
      x_var != "Replicate" &&
      (is.null(facet_rows) || !grepl("Replicate", facet_rows)) &&
      is.null(facet_cols)) {
    facet_cols <- "Replicate"
  }

  y_var <- if (type == "reps") "Peptides.type.number" else "Mean.Peptides.type.number"

  if (color_by == "none") {
    p <- ggplot(dt, aes(x = .data[[x_var]],
                        y = .data[[y_var]])) +
      geom_bar(stat = "identity", position = "stack", fill = def_color) +
      theme_pubr()
  } else {
    p <- ggplot(dt, aes(x = .data[[x_var]],
                        y = .data[[y_var]],
                        fill = .data[[color_by]])) +
      geom_bar(stat = "identity", position = "stack") +
      scale_fill_manual(values = protein_color) +
      theme_pubr()
  }

  if (scientific_10_y) {
    p <- p + scale_y_continuous(labels = scientific_10)
  }

  if (!is.null(facet_rows) || !is.null(facet_cols)) {
    rows <- if (!is.null(facet_rows)) facet_rows else "."
    cols <- if (!is.null(facet_cols)) facet_cols else "."
    p <- p + facet_grid(as.formula(paste(rows, "~", cols)))
  }

  vars_distinct <- grp_cols
  if (type == "reps") vars_distinct <- c(vars_distinct, "Replicate")
  vars_vary <- vars_distinct[
    sapply(vars_distinct, function(v) length(unique(dt[[v]])) > 1)
  ]
  used_vars <- x_var
  if (!is.null(facet_rows)) used_vars <- c(used_vars, strsplit(facet_rows, "\\+")[[1]])
  if (!is.null(facet_cols)) used_vars <- c(used_vars, strsplit(facet_cols, "\\+")[[1]])
  used_vars <- unique(used_vars)
  missing_vars <- setdiff(vars_vary, used_vars)
  if (length(missing_vars)) {
    warning(sprintf(
      "Grouping variable(s) %s not used in x_var/facets; data may be aggregated across them.",
      paste(missing_vars, collapse = ", ")
    ))
  }

  p
}
