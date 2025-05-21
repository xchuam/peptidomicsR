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

#' Kyte–Doolittle hydropathy indices for amino acids
#'
#' Named numeric vector of length 20 giving the hydropathy (hydrophobicity) index
#' for each of the standard amino acids, on the Kyte–Doolittle scale.
#'
#' @noRd
hydropathy <- c(
  A = 1.8, R = -4.5, N = -3.5, D = -3.5, C = 2.5,
  Q = -3.5, E = -3.5, G = -0.4, H = -3.2, I = 4.5,
  L = 3.8, K = -3.9, M = 1.9, F = 2.8, P = -1.6,
  S = -0.8, T = -0.7, W = -0.9, Y = -1.3, V = 4.2
)
