#' Plot peptide intensities with optional scientific notation y-axis, flexible filtering and faceting
#'
#' @description
#' Generates stacked‐bar charts of peptide intensities (per‐replicate or group‐mean)
#' with:
#' * User‐defined inclusion filters on any grouping columns,
#' * Dynamic x‐axis chosen from grouping variables,
#' * Flexible coloring by protein group or protein name,
#' * Optional scientific notationn of the y‐axis,
#' * Flexible faceting by one or more grouping variables.
#'
#' @param result List. Output of \code{processPeptides()}, containing at least:
#'   \itemize{
#'     \item \code{dt.peptides.int.reps}: replicate intensities,
#'     \item \code{dt.peptides.int}: mean intensities,
#'     \item \code{grp_cols}: character vector of grouping column names.
#'   }
#' @param type Character. Which table to plot: \code{"reps"} (replicate‐level) or \code{"mean"} (group‐mean). Default: \code{"mean"}.
#' @param x_var Character. Name of the grouping column to use on the x‐axis.
#'   Must be one of \code{result$grp_cols} (and for \code{type="reps"} may also be \code{"Replicate"}).
#'   Defaults to the first element of \code{result$grp_cols}.
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
#' @param scientific_10_y Logical.  If \code{TRUE}, use scientific notation for y-axis.
#'   Default: \code{TRUE}.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' result <- processPeptides(
#'   peptides_file          = "../Data/peptides.txt",
#'   intensity_columns_file = "../Data/Intensity_columns.csv",
#'   protein_mapping_file   = "../Data/protein_mapping.csv"
#' )
#'
#' # 1) Simple replicate‐level, default x, default facet="Replicate"
#' p1 <- plot_int(result, type = "reps")
#'
#' # 2) Simple mean‐level, x = Lipid, color by Protein.name, y axis not scientific notation
#' p2 <- plot_int(result,
#'                type             = "mean",
#'                x_var            = "Lipid",
#'                color_by         = "Protein.name",
#'                scientific_10_y  = FALSE)
#'
#' # 3) Complex filter: Lipid in N or S AND Digest.stage == "G"
#' p3 <- plot_int(result,
#'                type          = "reps",
#'                filter_params = list(Lipid = c("N","S"),
#'                                     Digest.stage = "G"))
#'
#' # 4) Complex filter to exclude a level (NOT): all Lipids except "N"
#' p4 <- plot_int(result,
#'                type          = "reps",
#'                filter_params = list(
#'                Lipid         = setdiff(unique(result$dt.peptides.int.reps$Lipid), "N")
#'                ))
#'
#' # 5) Multi‐variable faceting: rows = Casein.ratio+Digest.stage, cols = Lipid+Replicate
#' p5 <- plot_int(result,
#'                type        = "reps",
#'                facet_rows  = "Casein.ratio+Digest.stage",
#'                facet_cols  = "Lipid+Replicate")
#' }
#'
#' @import ggplot2
#' @import data.table
#' @importFrom ggpubr theme_pubr
#' @importFrom scales scientific_format
#' @export
plot_int <- function(result,
                     type    = c("mean", "reps"),
                     x_var   = NULL,
                     color_by = "Protein.group",
                     filter_params = NULL,
                     facet_rows    = NULL,
                     facet_cols    = NULL,
                     scientific_10_y   = TRUE) {
  # validate type
  type <- match.arg(type)

  # choose data
  dt      <- switch(type,
                    reps = copy(result$dt.peptides.int.reps),
                    mean = copy(result$dt.peptides.int))

  grp_cols <- result$grp_cols

  # set default x_var
  if (is.null(x_var)) {
    x_var <- grp_cols[1]
  }
  # validate x_var
  if (type == "reps") {
    if (!x_var %in% c("Replicate", grp_cols)) {
      stop("For type = 'reps', x_var must be 'Replicate' or one of result$grp_cols")
    }
  } else {
    if (!x_var %in% grp_cols) {
      stop("For type = 'mean', x_var must be one of result$grp_cols")
    }
  }

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
      x_var != "Replicate" &&
      (is.null(facet_rows) || !grepl("Replicate", facet_rows)) &&
      is.null(facet_cols)) {
    facet_cols <- "Replicate"
  }

  # set y variable
  y_var <- if (type == "reps") "Intensity" else "Mean.Intensity"

  # build the plot
  if (color_by == "none") {
    p <- ggplot(dt, aes(x = .data[[x_var]],
                        y = .data[[y_var]])) +
      geom_bar(stat = "identity", position = "stack",fill = def_color)+
      theme_pubr()
  } else {
    p <- ggplot(dt, aes(x = .data[[x_var]],
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
  used_vars <- x_var
  if (!is.null(facet_rows)) used_vars <- c(used_vars, strsplit(facet_rows, "\\+")[[1]])
  if (!is.null(facet_cols)) used_vars <- c(used_vars, strsplit(facet_cols, "\\+")[[1]])
  used_vars <- unique(used_vars)
  # find which truly‐varying group‐vars are missing
  missing_vars <- setdiff(vars_vary, used_vars)
  if (length(missing_vars)) {
    warning(sprintf(
      "Grouping variable(s) %s not used in x_var/facets; data may be aggregated across them.",
      paste(missing_vars, collapse = ", ")
    ))
  }

  p
}
