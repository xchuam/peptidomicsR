#' Process and summarize MaxQuant peptide data
#'
#' @description
#' Read peptide information from MaxQuant output file, intensity column meatadata, and protein mapping;
#' filters out contaminants; computes per-replicate and per-group mean intensities;
#' and counts peptides per replicate and per group .
#'
#' @param peptides_file
#'   A \code{character} specifying the path to the MaxQuant peptide data file,
#'   or a \code{data.frame}/\code{data.table} already loaded in R.
#' @param intensity_columns_file
#'   A \code{character} specifying the path to the CSV of intensity-column metadata,
#'   or a \code{data.frame}/\code{data.table} already loaded in R.
#' @param protein_mapping_file
#'   A \code{character} specifying the path to the CSV of protein mapping information,
#'   or a \code{data.frame}/\code{data.table} already loaded in R.
#'
#' @return A named list with the following elements:
#' * `dt.peptides`:  Filtered wide-format `data.table` of peptides with basic columns, intensity measurements, `Protein.name`, and `Protein.group`.
#' * `dt.peptides.int`: Data.table with `Mean.Intensity` per group across replicates.
#' * `dt.peptides.int.reps`: Data.table with `Intensity` per replicate per group.
#' * `dt.peptides.count`: Counts of peptides per group in `dt.peptides.int`.
#' * `dt.peptides.count.reps`: Counts of peptides per replicate per group in `dt.peptides.int.reps`.
#' * `dt.int_col`: The intensity column metadata read from `intensity_columns_file`.
#' * `grp_cols`: Character vector of grouping column names used.
#' * `peptides_select_col.basic`: Character vector of basic peptide columns selected initially.
#'
#' @import data.table
#' @export
#'
#' @examples
#' \dontrun{
#' # 1) Reading from files
#' result <- processPeptides(
#'   peptides_file          = "../Data/peptides.txt",
#'   intensity_columns_file = "../Data/Intensity_columns.csv",
#'   protein_mapping_file   = "../Data/protein_mapping.csv"
#' )
#'
#' # 2) Reading from in-memory files
#' library(data.table)
#' dt_pep <- fread("../Data/peptides.txt", integer64 = "double")
#' dt_lfq <- fread("../Data/Intensity_columns.csv")
#' dt_map <- fread("../Data/protein_mapping.csv")
#' result_2 <- processPeptides(
#'   peptides_file          = dt_pep,
#'   intensity_columns_file = dt_lfq,
#'   protein_mapping_file   = dt_map
#' )
#' }
processPeptides <- function(peptides_file,
                            intensity_columns_file,
                            protein_mapping_file) {
  # ensure data.table is available
  requireNamespace("data.table", quietly = TRUE)

  # 1. Read inputs: allow a path or an in-memory table
  # peptides
  if (is.character(peptides_file) && length(peptides_file) == 1) {
    if (!file.exists(peptides_file)) {
      stop("peptides_file path does not exist: ", peptides_file)
    }
    dt.peptides <- fread(peptides_file, integer64 = "double")
  } else if (data.table::is.data.table(peptides_file)) {
    dt.peptides <- data.table::copy(peptides_file)
  } else if (is.data.frame(peptides_file)) {
    dt.peptides <- data.table::as.data.table(peptides_file)
  } else {
    stop("`peptides_file` must be either a file path or a data.frame/data.table")
  }

  # intensity columns metadata
  if (is.character(intensity_columns_file) && length(intensity_columns_file) == 1) {
    if (!file.exists(intensity_columns_file)) {
      stop("intensity_columns_file path does not exist: ", intensity_columns_file)
    }
    dt.int_col <- fread(intensity_columns_file)
  } else if (data.table::is.data.table(intensity_columns_file)) {
    dt.int_col <- data.table::copy(intensity_columns_file)
  } else if (is.data.frame(intensity_columns_file)) {
    dt.int_col <- data.table::as.data.table(intensity_columns_file)
  } else {
    stop("`intensity_columns_file` must be either a file path or a data.frame/data.table")
  }

  # protein mapping
  if (is.character(protein_mapping_file) && length(protein_mapping_file) == 1) {
    if (!file.exists(protein_mapping_file)) {
      stop("protein_mapping_file path does not exist: ", protein_mapping_file)
    }
    dt.p_map <- fread(protein_mapping_file)
  } else if (data.table::is.data.table(protein_mapping_file)) {
    dt.p_map <- data.table::copy(protein_mapping_file)
  } else if (is.data.frame(protein_mapping_file)) {
    dt.p_map <- data.table::as.data.table(protein_mapping_file)
  } else {
    stop("`protein_mapping_file` must be either a file path or a data.frame/data.table")
  }

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
    is.na(Protein.group), Protein.group := "Whey"
  ][,
    `:=`(
      Protein.name  = factor(Protein.name,
                             levels = c(
                               unique(dt.p_map[["Protein.name"]]),
                               "Others"
                               )),
      Protein.group = factor(Protein.group,
                             levels = c(
                               unique(dt.p_map[["Protein.group"]])
                             ))
    )]

  # 5. Identify grouping columns (exclude Replicate cols)
  grp_cols <- setdiff(names(dt.int_col), c("Intensity.column", "Replicate"))

  # 6. Build nested lists
  # deepest level is list(mean = <mean DT>, reps = <each rep DT>)
  build_nested_list <- function(meta, grp_cols) {
    this_col <- grp_cols[1]
    vals     <- unique(meta[[this_col]])
    out      <- setNames(vector("list", length(vals)), as.character(vals))

    if (length(grp_cols) > 1) { #more levels need to go
      for (v in vals) {
        sub_meta <- meta[get(this_col) == v]
        out[[as.character(v)]] <- build_nested_list(sub_meta, grp_cols[-1])
      }
    } else { #at the deepest level
      for (v in vals) {
        # 1) prepare the table for mean intensity
        sub_meta  <- meta[get(this_col) == v]
        cols      <- c(peptides_select_col.basic,
                       "Protein.name",
                       "Protein.group",
                       sub_meta[["Intensity.column"]])
        dt.mean    <- dt.peptides[, ..cols]

        # 2) rename to R1…Rn
        n_reps    <- length(sub_meta[["Intensity.column"]])
        rep_names <- paste0("R", seq_len(n_reps))
        setnames(dt.mean,
                 old = cols,
                 new = c(peptides_select_col.basic,
                         "Protein.name",
                         "Protein.group",
                         rep_names))

        # 3) keep peptides with at least half replicates >0
        keep_rows <- rowSums(dt.mean[, ..rep_names] > 0) >= ceiling(n_reps/2)
        dt.mean    <- dt.mean[keep_rows]

        # 4) mean intensity of non-zero replicates
        dt.mean[, Mean.Intensity := rowSums(.SD) / rowSums(.SD > 0),
               .SDcols = rep_names]

        # 5) make the table of each replicate
        dt.rep <- melt(
          copy(dt.mean),
          id.vars       = c(peptides_select_col.basic,
                            "Protein.name","Protein.group"),
          measure.vars = rep_names,
          variable.name = "Replicate",
          value.name    = "Intensity"
        )
        # drop any NA or zero intensities
        dt.rep <- dt.rep[!is.na(Intensity) & Intensity > 0]

        out[[as.character(v)]] <- list(mean = dt.mean,
                                       reps = dt.rep)
      }
    }
    out
  }

  nested.peptides <- build_nested_list(dt.int_col, grp_cols)

  # 7. Flatten nested list back to long data.table
  flatten_nested <- function(nl, grp_cols, leaf_name) {
    out <- list()
    rec <- function(node, keys = list()) {
      # if at a leaf, it will be a list with an element called leaf_name
      if (is.data.table(node[[leaf_name]])) {
        dt <- copy(node[[leaf_name]])
        # attach the grouping columns
        for (nm in names(keys)) {
          dt[, (nm) := keys[[nm]]]
        }
        out[[length(out) + 1]] <<- dt
      } else {
        # still in the nested-list world: descend one level
        this_group <- grp_cols[length(keys) + 1]
        for (val in names(node)) {
          rec(node[[val]],
              c(keys, setNames(list(val), this_group)))
        }
      }
    }
    rec(nl)
    rbindlist(out, fill = TRUE)
  }

  dt.peptides.int <- flatten_nested(nested.peptides, grp_cols, "mean")
  dt.peptides.int.reps <- flatten_nested(nested.peptides, grp_cols, "reps")

  # 8. Convert grouping columns to ordered factors based on original file order
  for (col in grp_cols) {
    levels_vec <- unique(dt.int_col[[col]])
    dt.peptides.int[, (col) := factor(get(col), levels = levels_vec)]
    dt.peptides.int.reps[, (col) := factor(get(col), levels = levels_vec)]
  }

  # 9. compute count of peptides per Length, Protein.name, Protein.group, grp_cols, annd/or Replicate
  dt.peptides.count <- dt.peptides.int[
    , .(Mean.Count = .N)
    , by = c("Length", "Protein.name", "Protein.group", grp_cols)
  ]
  dt.peptides.count.reps <- dt.peptides.int.reps[
    , .(Count = .N)
    , by = c("Length", "Protein.name", "Protein.group", "Replicate", grp_cols)
  ]


  # 9. Return outputs
  list(
    dt.peptides               = dt.peptides,
    dt.peptides.int           = dt.peptides.int,
    dt.peptides.int.reps      = dt.peptides.int.reps,
    dt.peptides.count         = dt.peptides.count,
    dt.peptides.count.reps    = dt.peptides.count.reps,
    dt.int_col                = dt.int_col,
    grp_cols                  = grp_cols,
    peptides_select_col.basic = peptides_select_col.basic

  )
}
