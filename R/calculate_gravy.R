#' Calculate the GRAVY (Grand Average of Hydropathy) score for a peptide
#'
#' @description
#' Computes the average hydropathy index of a peptide sequence using the
#' \code{hydropathy} vector (Kyte–Doolittle scale). A positive value indicates
#' a hydrophobic peptide; negative indicates hydrophilic.
#'
#' @param peptide Character(1). A single amino acid sequence.
#'
#' @return Numeric(1). The GRAVY score: the sum of the per-residue hydropathy
#'   values divided by the peptide length.
#'
#' @examples
#' \dontrun{
#' calculate_gravy("ACDEFGHIK")         # for a 9-mer
#' peptides <- c("ALWKTML", "GGGGGGG")
#' sapply(peptides, calculate_gravy)
#' }
#'
#' @noRd
#' @keywords internal
calculate_gravy <- function(peptide) {
  # Split the peptide sequence into individual amino acids
  amino_acids <- unlist(strsplit(peptide, split = ""))
  # Sum the hydropathy values of the amino acids in the peptide
  total_hydropathy <- sum(hydropathy[amino_acids])
  # Calculate the GRAVY score
  gravy_score <- total_hydropathy / length(amino_acids)
  return(gravy_score)
}
