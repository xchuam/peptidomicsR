#' Process MaxQuant peptide data and compute mean intensities across replicates
#'
#' @description
#' Read peptide information from MaxQuant output file, intensity columns, and protein mapping,
#' filter peptides, compute mean intensities across replicates, and return
#' both a wide and a long `data.table`.
#'
#' @param peptides_file Character. Path to the MaxQuant peptide data file.
#' @param intensity_columns_file Character. Path to the CSV file describing intensity columns and other metadata.
#' @param protein_mapping_file Character. Path to the CSV file containing protein mapping information.
#'
#' @return A list with two elements:
#' * `dt.peptides`: A `data.table` of filtered peptide data with intensity columns.
#' * `dt.peptides.long`: A `data.table` in long format with grouping columns and `mean.intensity`
#'
#' @import data.table
#' @export
#'
#' @examples
#' \dontrun{
#' result <- processPeptides(
#'   peptides_file          = "../Data/peptides.txt",
#'   intensity_columns_file = "../Data/Intensity_columns.csv",
#'   protein_mapping_file   = "../Data/protein_mapping.csv"
#' )
#' dt_wide <- result$dt.peptides
#' dt_long <- result$dt.peptides.long
#' }
processPeptides <- function(peptides_file,
                            intensity_columns_file,
                            protein_mapping_file) {
  # ensure data.table is available
  requireNamespace("data.table", quietly = TRUE)

  # 1. Read input files
  dt.peptides <- fread(peptides_file, integer64 = "double")
  dt.int_col  <- fread(intensity_columns_file)
  dt.p_map    <- fread(protein_mapping_file)

  # 2. Clean column names: replace spaces with dots
  setnames(dt.peptides,
           old = names(dt.peptides),
           new = gsub(" ", ".", names(dt.peptides)))
  dt.int_col[, Intensity.column := gsub(" ", ".", Intensity.column)]

  # 3. Define columns to select: basic peptide info + all intensity columns
  peptides_select_col.basic <- c(
    "Sequence",
    "Leading.razor.protein",
    "Length",
    "Start.position",
    "End.position",
    "Amino.acid.before",
    "First.amino.acid",
    "Last.amino.acid",
    "Amino.acid.after"
  )
  peptides_select_col <- c(
    peptides_select_col.basic,
    dt.int_col$Intensity.column
  )

  # 4. Subset to selected columns, then filter and map proteins
  dt.peptides <- dt.peptides[
    , ..peptides_select_col
  ][
    # remove contaminants and reverse hits
    !(`Leading.razor.protein` %like% "^CON_" |
        `Leading.razor.protein` %like% "^REV_")
  ][
    # map protein name and group
    dt.p_map, on = "Leading.razor.protein",
    `:=`(
      Protein.name  = i.Protein.name,
      Protein.group = i.Protein.group
    )
  ][
    # fill missing mappings
    is.na(Protein.name),  Protein.name  := "Others"
  ][
    is.na(Protein.group), Protein.group := "whey"
  ]

  # 5. Identify grouping columns (exclude metadata cols)
  grp_cols <- setdiff(names(dt.int_col), c("Intensity.column", "Replicate"))

  # 6. Build nested list of replicates at deepest level
  build_nested_list <- function(meta, grp_cols) {
    this_col <- grp_cols[1]
    vals     <- unique(meta[[this_col]])
    out      <- setNames(vector("list", length(vals)), as.character(vals))

    if (length(grp_cols) > 1) {
      for (v in vals) {
        sub_meta <- meta[get(this_col) == v]
        out[[as.character(v)]] <- build_nested_list(sub_meta, grp_cols[-1])
      }
    } else {
      for (v in vals) {
        sub_meta  <- meta[get(this_col) == v]
        cols      <- c(peptides_select_col.basic,
                       sub_meta[["Intensity.column"]])
        dt.rep    <- dt.peptides[, ..cols]
        n_reps    <- length(sub_meta[["Intensity.column"]])
        rep_names <- paste0("R", seq_len(n_reps))
        setnames(dt.rep,
                 old = cols,
                 new = c(peptides_select_col.basic, rep_names))

        # keep peptides with at least half replicates >0
        keep_rows <- rowSums(dt.rep[, ..rep_names] > 0) >= ceiling(n_reps/2)
        dt.rep    <- dt.rep[keep_rows]

        # mean intensity of non-zero replicates
        dt.rep[, mean.intensity := rowSums(.SD) / rowSums(.SD > 0),
               .SDcols = rep_names]

        out[[as.character(v)]] <- dt.rep
      }
    }
    out
  }

  nested.peptides <- build_nested_list(dt.int_col, grp_cols)

  # 7. Flatten nested list back to long data.table
  flatten_nested <- function(nlist, grp_cols) {
    out <- list()
    rec <- function(node, keys = list()) {
      if (inherits(node, "data.table")) {
        dt <- copy(node)
        for (nm in names(keys)) dt[, (nm) := keys[[nm]]]
        out[[length(out) + 1]] <<- dt
      } else {
        this_group <- grp_cols[length(keys) + 1]
        for (val in names(node)) {
          rec(node[[val]], c(keys, setNames(list(val), this_group)))
        }
      }
    }
    rec(nlist)
    rbindlist(out, fill = TRUE)
  }

  dt.peptides.long <- flatten_nested(nested.peptides, grp_cols)

  # 8. Return both wide and long tables
  list(
    dt.peptides      = dt.peptides,
    dt.peptides.long = dt.peptides.long
  )
}
