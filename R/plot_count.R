#' Plot peptide type numbers (deprecated wrapper)
#'
#' @description
#' \code{plot_count()} is kept for backward compatibility, but it has been replaced by
#' \code{plot_type_num()}. It plots numbers of unique peptide types using the same
#' functionality as \code{plot_type_num()}.
#'
#' @inheritParams plot_type_num
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
#' plot_count(result, type = "reps")
#' }
#'
#' @export
plot_count <- function(result,
                       type             = c("mean", "reps"),
                       x_var            = NULL,
                       color_by         = "Protein.group",
                       filter_params    = NULL,
                       facet_rows       = NULL,
                       facet_cols       = NULL,
                       scientific_10_y  = FALSE) {
  warning(
    "plot_count() is deprecated and has been replaced by plot_type_num(). Please use plot_type_num() instead.",
    call. = FALSE
  )

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
