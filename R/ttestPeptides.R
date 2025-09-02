#' Compare peptide intensities between user-defined groups (t-test)
#'
#' @description
#' Performs two-sample t-tests on peptide replicate intensities between
#' user-specified groups. Groups are defined by levels in the columns listed
#' in `result$grp_cols`, and replicate values are taken from columns `R1..Rn`
#' in `result$dt.peptides.int.ttest`. Users can specify selectors as
#' underscore-joined tokens (e.g., `"C1_B1_T1"`), in any order (e.g.,
#' `"B1_T1_C1"`), and may omit some grouping columns (e.g., `"C1_B1"`). If a
#' selector omits some grouping columns, replicates across the omitted
#' dimensions are pooled (a warning is issued).
#'
#' Zero intensities are replaced by a pseudocount provided by the user with default value 1.
#'
#' @param result List. Output of \code{processPeptides()}. Must contain:
#'   \itemize{
#'     \item \code{dt.peptides.int.ttest}: data.table with peptide intensities for t.test analysis.
#'     \item \code{grp_cols}: character vector of grouping columns to define sample groups.
#'   }
#' @param comparisons Either:
#'   * a list of character pairs, e.g. `list(list(c("C40_G_N", "C40_G_S")), list(c("C40_G", "C40_G")))`,
#'   or
#'   * a single character vector of length 2, e.g. `c("C40_G_N", "C40_G_S")`.
#'   Each selector is a character string describing one set of samples, using the levels from
#'   `result$grp_cols` joined by underscores. Order of tokens is flexible.
#'   Examples:
#'   * `"C40_G_N"` vs `"C40_G_S"`: full specification (Casein.ratio = C40,
#'     Digest.stage = G, Lipid = N vs S).
#'   * `"C40_N_G"` vs `"C40_S_G"`: same as above, token order does not matter.
#'   * `"C40_N"` vs `"C40_S"`: partial specification. Here, samples are pooled
#'     across unspecified grouping columns (warning issued).
#'   Disambiguation: If a token could belong to more than one grouping column,
#'   use the `col=value` form. For example, `"Digest.stage=G_Lipid=G"` explicitly
#'   selects `Digest.stage = G` and `Lipid = G`.
#' @param pseudocount Numeric. Value used to replace zero intensities before
#'   statistical testing. This helps avoid issues with log2 transformation and division by zero.
#'   Default: 1.
#' @param equal_var Logical. If `TRUE`, use Student's t-test with equal
#'   variances; otherwise Welch's test (default).
#' @param alternative Character. `"two.sided"` (default), `"less"`, or `"greater"`.
#' @param adjust Character. Method passed to [stats::p.adjust()], default `"BH"`.
#' @param min_reps_per_side Integer. Minimum non-NA observations required on
#'   each side to perform a test. Default: 2.
#'
#' @return A named list of `data.table`, one per comparison. Each table has one
#'   row per peptide and includes:
#'   `t`, `df`, `p.value`, `log2FC`, `p.adj`.
#'
#' @examples
#' \dontrun{
#' # Prepare data
#' result <- processPeptides(
#'   peptides_file          = "../Data/peptides.txt",
#'   intensity_columns_file = "../Data/Intensity_columns.csv",
#'   protein_mapping_file   = "../Data/protein_mapping.csv"
#' )
#'
#'
#' #  1) One t.test with specification of groups
#' ttest.1 <- ttestPeptides(
#'   result,
#'   comparisons = list(c("C40_G_N", "C40_I_N"))
#' )
#'
#' #  2) For one t.test, the comparisons can be a a single character vector
#' ttest.2 <- ttestPeptides(
#'   result,
#'   comparisons = c("C40_G_N", "C40_I_N")
#' )
#'
#' #  3) Token order does not matter
#' # Same comparison, different token order
#' ttest.3 <- ttestPeptides(
#'   result,
#'   comparisons = list(c("G_N_C40", "I_N_C40"))
#' )
#'
#' #  4) More than one t.test can be done at one time
#' ttest.4 <- ttestPeptides(
#'   result,
#'   comparisons = list(
#'     c("C40_G_N", "C40_I_N"),
#'     c("C40_G_N", "C80_I_N")
#'   )
#' )
#'
#' #  5) Partial specification (pools over omitted column)
#' ttest.5 <- ttestPeptides(
#'   result,
#'   comparisons = list(c("C40_G", "C80_G"))
#' )
#'
#' #  6) Explicit disambiguation using col=value syntax
#' # Equivalent to Example 1 but written with explicit column names
#' ttest.6 <- ttestPeptides(
#'   result,
#'   comparisons = list(c("Casein.ratio=C40_Digest.stage=G_Lipid=N",
#'                        "Casein.ratio=C40_Digest.stage=I_Lipid=N"))
#' )
#'
#'
#' # Access the first comparison result table
#' ttest.1[[1]][order(p.adj)][1:10]
#' }
#'
#' @import data.table
#' @importFrom stats t.test p.adjust
#' @export
ttestPeptides <- function(
    result,
    comparisons,
    pseudocount       = 1,
    equal_var         = FALSE,
    alternative       = c("two.sided","less","greater"),
    adjust            = "BH",
    min_reps_per_side = 2
) {
  alternative <- match.arg(alternative)

  requireNamespace("data.table", quietly = TRUE)

  dt <- data.table::copy(result$dt.peptides.int.ttest)
  grp_cols <- result$grp_cols
  id_cols <- c(result$peptides_select_col.basic,
               "Protein.name","Protein.group")

  # --- locate replicate columns R1..Rn ---
  rep_cols <- grep("^R\\d+$", names(dt), value = TRUE)
  if (length(rep_cols) == 0L)
    stop("No replicate columns 'R1..Rn' found in result$dt.peptides.int.")

  # --- make long table of replicate intensities (keep zeros here) ---
  long <- data.table::melt(
    dt,
    id.vars      = c(id_cols, grp_cols),
    measure.vars = rep_cols,
    variable.name = "Replicate",
    value.name    = "Intensity"
  )

  # zero replacement by pseudocount
  long[, Intensity_adj := ifelse(Intensity == 0, pseudocount, Intensity)]
  # choose the working value (optionally log2)
  long[, value := log2(Intensity_adj)]

  # helper: unique levels by column
  levels_by_col <- lapply(grp_cols, function(cn) unique(as.character(dt[[cn]])))
  names(levels_by_col) <- grp_cols

  # parse a selector like "C1_B1_T1" or "B1_T1_C1" or "time=T1_batch=B1"
  parse_selector <- function(sel) {
    parts <- strsplit(sel, "_", fixed = TRUE)[[1]]
    mapping <- list()
    used <- character(0)

    for (tok in parts) {
      if (grepl("[:=]", tok)) {
        sp <- strsplit(tok, "[:=]")[[1]]
        if (length(sp) != 2) stop("Malformed token: '", tok, "'. Use 'col=val'.")
        col <- sp[1]; val <- sp[2]
        if (!col %in% grp_cols) stop("Unknown grouping column '", col, "' in selector '", sel, "'. Check result$grp_cols.")
        mapping[[col]] <- val
        used <- c(used, col)
        next
      }
      # infer column by value
      hits <- names(Filter(function(vs) tok %in% vs, levels_by_col))
      if (length(hits) == 1L) {
        col <- hits
        mapping[[col]] <- tok
        used <- c(used, col)
      } else if (length(hits) == 0L) {
        stop(tok, "' not found in any grouping column as indicated by result$grp_cols.")
      } else {
        stop("Token '", tok, "' is ambiguous (matches columns: ",
             paste(hits, collapse = ", "),
             "). Disambiguate using 'col=val', e.g. 'time=", tok, "'.")
      }
    }
    if (any(duplicated(used))) stop("Selector '", sel, "' assigns the same column more than once.")
    mapping
  }

  #function to filter dt_long based on the condition translated by parse_selector()
  subset_by_selector <- function(dt_long, sel_map) {
    out <- dt_long
    for (nm in names(sel_map)) {
      out <- out[get(nm) == sel_map[[nm]]]
    }
    out
  }

  make_label <- function(sel_map_raw) {
    # pretty label in declared grp_cols order; only show specified pieces
    paste(paste(names(sel_map_raw), sel_map_raw, sep = "="), collapse = ", ")
  }

  # normalize comparisons input to list of pairs
  if (is.character(comparisons) && length(comparisons) == 2L) {
    comparisons <- list(comparisons)
  }
  if (!is.list(comparisons) || !all(vapply(comparisons, function(x) is.character(x) && length(x) == 2L, logical(1)))) {
    stop("`comparisons` must be a list of character pairs, or a character vector of length 2.")
  }

  res_list <- vector("list", length(comparisons))
  names(res_list) <- vapply(
    comparisons,
    function(ab) paste0(ab[1], "_vs_", ab[2]),
    character(1)
  )

  for (i in seq_along(comparisons)) {
    ab <- comparisons[[i]]
    selA_raw <- ab[1]; selB_raw <- ab[2]
    selA <- parse_selector(selA_raw)
    selB <- parse_selector(selB_raw)

    # warn if partial (pooling across omitted columns)
    if (length(selA) < length(grp_cols) || length(selB) < length(grp_cols)) {
      warning("One or both selectors omit grouping columns; replicates will be pooled across omitted dimensions.\n",
              "  A: ", make_label(selA), "\n",
              "  B: ", make_label(selB))
    }

    A <- subset_by_selector(long, selA)
    B <- subset_by_selector(long, selB)

    # summarize per peptide and run t-tests
    # prepare A
    A_sum <- A[, .(
      A_vals = list(value[is.finite(value)]),
      n_A    = sum(is.finite(value)),
      mean_A = mean(value, na.rm = TRUE)
    ), by = id_cols]
    # prepare B
    B_sum <- B[, .(
      B_vals = list(value[is.finite(value)]),
      n_B    = sum(is.finite(value)),
      mean_B = mean(value, na.rm = TRUE)
    ), by = id_cols]

    # merge; keep all peptides that appear in at least one side
    M <- merge(A_sum, B_sum, by = id_cols, all = TRUE)

    # remove rows with only pseudocount
    M <- M[!(mean_A == log2(pseudocount) & mean_B == log2(pseudocount))]

    # rowwise t-test
    M[, c("t","df","p.value") := {
      a <- unlist(A_vals); b <- unlist(B_vals)
      if (length(a) >= min_reps_per_side && length(b) >= min_reps_per_side) {
        tt <- try(stats::t.test(a, b, var.equal = equal_var, alternative = alternative), silent = TRUE)
        if (inherits(tt, "try-error")) list(NA_real_, NA_real_, NA_real_)
        else list(unname(tt$statistic), unname(tt$parameter), unname(tt$p.value))
      } else {
        list(NA_real_, NA_real_, NA_real_)
      }
    }, by = id_cols]

    # log2FC difference of means on log2 scale
    M[, log2FC := mean_A - mean_B]

    # p.adjust within this contrast
    M[, p.adj := stats::p.adjust(p.value, method = adjust)]

    #rename some columns
    setnames(M,
             c("A_vals", "n_A", "mean_A"),
             c(paste0("vals.", selA_raw),
               paste0("n.", selA_raw),
               paste0("mean.", selA_raw)))
    setnames(M,
             c("B_vals", "n_B", "mean_B"),
             c(paste0("vals.", selB_raw),
               paste0("n.", selB_raw),
               paste0("mean.", selB_raw)))

    # sort by p.adj
    data.table::setorder(M, p.adj)

    res_list[[i]] <- M[]
  }

  res_list
}
