#' Filter peptide results for downstream analysis
#'
#' @description
#' Subsets the full `processPeptides()` result to only those peptides you’re interested in,
#' either by exact sequence, by regex pattern, or by any combination of your grouping variables.
#' Returns the same list-of-tables structure so it will work with `plot_int()`, `plot_count()`, etc.
#'
#' @param result List. The output of \code{processPeptides()} (must contain at least
#'   \code{dt.peptides}, \code{dt.peptides.int}, \code{dt.peptides.int.reps},
#'   \code{dt.peptides.count}, \code{dt.peptides.count.reps}, \code{grp_cols},
#'   \code{dt.int_col}, and \code{peptides_select_col.basic}).
#' @param seqs Character vector, or \code{NULL}.  If provided, only exact-matching \code{Sequence}
#'   values are kept.
#' @param seq_pattern Character(1), or \code{NULL}.  If provided, a regular expression;
#'   any \code{Sequence} matching this pattern is kept.
#'   (You may supply both \code{seqs} and \code{seq_pattern}, in which case the union is kept.)
#' @param filter_params Named list, or \code{NULL}.
#' Names are grouping columns (elements of \code{result$grp_cols}), values are
#'   vectors of levels to include (multiple entries are ANDed).
#'
#' @return A list with the identical named elements as the original \code{result}, but
#'   each of the five data.tables (\code{dt.peptides}, \code{dt.peptides.int},
#'   \code{dt.peptides.int.reps}, \code{dt.peptides.count}, \code{dt.peptides.count.reps}, and corresponding metadata \code{dt.int_col})
#'   is filtered to only the selected sequences/groups.  The rest part
#'   (\code{grp_cols}, \code{peptides_select_col.basic}) is passed through untouched.
#'
#' @examples
#' \dontrun{
#' result <- processPeptides(
#'   peptides_file          = "../Data/peptides.txt",
#'   intensity_columns_file = "../Data/Intensity_columns.csv",
#'   protein_mapping_file   = "../Data/protein_mapping.csv"
#' )
#'
#' # 1) exact sequence list
#' sub_result1 <- filterPeptides(
#'   result = result,
#'   seqs   = c("AAEVIREA","AAEVIREALQGITDPLFK")
#' )
#'
#' # 2) regex pattern
#' sub_result2 <- filterPeptides(
#'   result      = result,
#'   seq_pattern = "^AEVIR"
#' )
#'
#' # 3) combine with grouping filters
#' sub_result3 <- filterPeptides(
#'   result        = result,
#'   seq_pattern   = "PLFK$",
#'   filter_params = list(Lipid = "N", Digest.stage = "G")
#' )
#'
#' # Then you can do e.g.
#' plot_int(sub3)
#' plot_count(sub3, type = "reps")
#' }
#'
#' @import data.table
#' @export
filterPeptides <- function(result,
                           seqs        = NULL,
                           seq_pattern = NULL,
                           filter_params = NULL) {
  requireNamespace("data.table", quietly=TRUE)
  # copy required tables
  dt_p       <- data.table::copy(result$dt.peptides)
  dt_mean    <- data.table::copy(result$dt.peptides.int)
  dt_mean_r  <- data.table::copy(result$dt.peptides.int.reps)
  dt_int_col <- data.table::copy(result$dt.int_col)
  grp_cols   <- data.table::copy(result$grp_cols)

  # build a mask of sequences to keep
  keep_seqs <- character(0) # empty
  all_seqs <- unique(dt_p$Sequence)

  if (!is.null(seqs)) {
    keep_seqs <- union(keep_seqs, seqs)
  }
  if (!is.null(seq_pattern)) {
    keep_seqs <- union(keep_seqs, grep(seq_pattern,
                                       all_seqs,
                                       value=TRUE))
  }
  # if neither specified, default to all
  if (length(keep_seqs)==0 && is.null(filter_params)) {
    keep_seqs <- all_seqs
  }

  # apply Sequence filter
  dt_p       <- dt_p[Sequence %in% keep_seqs]
  dt_mean    <- dt_mean[Sequence %in% keep_seqs]
  dt_mean_r  <- dt_mean_r[Sequence %in% keep_seqs]

  # compute count of peptides per Length, Protein.name, Protein.group, grp_cols, annd/or Replicate
  dt_count <- dt_mean[
    , .(Mean.Count = .N)
    , by = c("Length", "Protein.name", "Protein.group", grp_cols)
  ]
  dt_count_r <- dt_mean_r[
    , .(Count = .N)
    , by = c("Length", "Protein.name", "Protein.group", "Replicate", grp_cols)
  ]

  # apply grouping filters if requested
  if (!is.null(filter_params)) {
    stopifnot(is.list(filter_params))
    for (col in names(filter_params)) {
      vals <- filter_params[[col]]

      #wild table filter
      dt_mean    <- dt_mean[get(col) %in% vals]
      dt_mean_r  <- dt_mean_r[get(col) %in% vals]
      dt_count   <- dt_count[get(col) %in% vals]
      dt_count_r <- dt_count_r[get(col) %in% vals]

      #intensity‐column filter for long table
      dt_int_col <- dt_int_col[get(col) %in% vals]
    }
    #long table filter
    dt_p <- dt_p[, c(result$peptides_select_col.basic,
                     dt_int_col[["Intensity.column"]],
                     "Protein.name",
                     "Protein.group"), with = FALSE]
  }

  # return same structure
  list(
    dt.peptides              = dt_p,
    dt.peptides.int          = dt_mean,
    dt.peptides.int.reps     = dt_mean_r,
    dt.peptides.count        = dt_count,
    dt.peptides.count.reps   = dt_count_r,
    dt.int_col               = dt_int_col,
    grp_cols                 = grp_cols,
    peptides_select_col.basic= result$peptides_select_col.basic
  )
}
