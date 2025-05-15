#' Internal constant color mapping for protein names
#'
#' A named character vector specifying colors for each protein name in plots.
#'
#' @noRd
protein_name_color = c("Beta-casein" =  "#E41A1C",
                       "Kappa-casein" = "#377EB8",
                       "Alpha-S1-casein" = "#4DAF4A",
                       "Alpha-S2-casein"= "#984EA3",
                       "Alpha-lactalbumin" = "#FF7F00",
                       "Beta-lactoglobulin" = "#FCC88F",
                       "Others" = "gray40")

#' Internal constant color mapping for protein groups
#'
#' A named character vector specifying colors for each protein group.
#'
#' @noRd
protein_group_color <- c(
  "Casein" = "#E8A03F",
  "Whey"   = "#0073B4"
)

#' Internal combined color mapping for proteins
#'
#' A named character vector merging name- and group-level colors for use in plotting.
#'
#' @noRd
protein_color <- c(
  protein_name_color,
  protein_group_color
)

#' Internal color for no mapping
#'
#' A default colors for use in plotting without mapping.
#'
#' @noRd
def_color <- "#9AC9DB"
