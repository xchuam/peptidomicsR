#' Convert axis breaks to expression-format scientific notation (internal)
#'
#' @description
#' Internal helper to format numeric axis break values into plotmath expressions
#' suitable for use with \\code{scale_y_continuous(labels = scientific_10)} when
#' using a log10 transformation. Values greater than 10 are displayed in 10^x format.
#'
#' @param x Numeric vector of axis break values.
#' @return A list of plotmath expressions for axis labels, or NA for non-numeric/NA.
#'
#' @details
#' This function is not exported and is intended for internal use by plotting functions
#' within this package. It leverages \\code{scales::scientific_format()} to generate
#' scientific notation strings, then converts them to plotmath expressions.
#'
#' @noRd
#' @keywords internal
#' @importFrom scales scientific_format
scientific_10 <- function(x) {
  sapply(x, function(xi) {
    if (is.na(xi) || !is.numeric(xi)) {
      NA
    } else if (xi > 10) {
      # Convert to scientific notation for xi > 10 and then to the desired format
      formatted <- scientific_format()(xi)
      # Replace 'e+' with ' %*% 10^' for plotmath
      expr_text <- gsub("e\\+", " %*% 10^", formatted)
      parse(text = expr_text)
    } else {
      # Directly use integer values without exponentiation
      as.expression(xi)
    }
  }, USE.NAMES = FALSE)
}
