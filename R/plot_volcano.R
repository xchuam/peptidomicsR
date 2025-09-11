#' Volcano plot(s) for peptide differential tests (t.test or limma TREAT)
#'
#' @description
#' Create one or multiple volcano plots from a \code{ttestPeptides()} result list.
#' Each plot shows \code{log2FC} on the x-axis and \eqn{-log10(adjusted~p)} on the y-axis,
#' optionally with a log-scaled y-axis. Points are colored
#' by significance using the \code{sig} column produced by \code{ttestPeptides()}.
#' Threshold lines are drawn depending on \code{test_method}:
#' \itemize{
#'   \item \code{"treat"}: a horizontal line at \eqn{-log10(\alpha)} (log2 fold-change threshold is baked into TREAT).
#'   \item \code{"plain"}: horizontal line at \eqn{-log10(\alpha)} plus vertical lines at
#'         \eqn{\pm} \code{lfc_thresh}.
#' }
#'
#' @details
#' \strong{Selecting comparisons.} \code{comparisons} supports two modes:
#' \enumerate{
#'   \item \emph{Index mode}: a single integer, numeric vector, or list of integers.
#'         e.g., \code{1}, \code{c(1,3)}, or \code{list(2,4)}.
#'         These index into \code{ttest_result}. The output list keeps the corresponding names
#'         from \code{ttest_result} for those indices.
#'   \item \emph{Name-pair mode}: a character vector of length 2 (e.g., \code{c("A","B")}) or a list of such vectors.
#'         For each pair \code{c("A","B")}, the function looks up the element named \code{"A_vs_B"} in \code{ttest_result}.
#' }
#'
#' \strong{Labels and highlights.} The function support to label and highlight the points of selected sequences.
#' \itemize{
#'   \item \code{label_seqs}: sequences to annotate with non-overlapping text.
#'   \item \code{highlight_seqs}: a \emph{named list} mapping border color \eqn{\rightarrow} character vector of sequences.
#'         Each color creates a ring around the selected points.
#'         e.g., \code{list("red" = c(Seq1, Seq2), "blue" = c(Seq3, Seq4))}
#' }
#'
#' @param ttest_result A named list of \code{data.table}s as returned by \code{ttestPeptides()},
#'   where each element is one comparison result table containing at least
#'   \code{Sequence}, \code{log2FC}, \code{p.adj}, and \code{sig}.
#' @param comparisons Comparison selector. Either:
#'   \enumerate{
#'     \item index mode: an integer, numeric vector, or list of integers indexing \code{ttest_result}, or
#'     \item name-pair mode: a character vector of length 2 (e.g., \code{c("A","B")}) or a list of such vectors,
#'           which will be resolved to names like \code{"A_vs_B"} in \code{ttest_result}.
#'   }
#' @param test_method Character; either \code{"treat"} (default) or \code{"plain"}.
#'   Controls which threshold lines are drawn, matching how the test was performed.
#'   This should match how the tests were performed in \code{ttestPeptides()}.
#' @param show_threshold Logical; draw threshold lines. Default \code{TRUE}.
#' @param lfc_thresh Numeric; absolute log2 fold-change threshold for the \code{"plain"} method
#'   (ignored for drawing when \code{test_method = "treat"}).
#'   This should match how the tests were performed in \code{ttestPeptides()}.
#'   Default \code{1}.
#' @param alpha Numeric in (0,1]; adjusted p-value cutoff for threshold line.
#'   This should match how the tests were performed in \code{ttestPeptides()}.
#'   Default \code{0.05}.
#' @param fill_values Named character vector of length 2 giving colors for significant and not sinificant points.
#'   \code{c(no = "...", yes = "...")}, used by \code{scale_fill_manual()} and \code{scale_color_manual()}.
#'   Default \code{c(no = "grey75", yes = "#FFC010")}.
#' @param point_size Numeric; point size for the base layer. Default \code{2}.
#' @param point_alpha Numeric in [0,1]; base point transparency. Default \code{0.85}.
#' @param label_seqs Character vector of peptide sequences to label by text (non-overlapping) and rings. Default \code{NULL}.
#' @param label_size Numeric; text size for labels. Default \code{3}.
#' @param label_col Character; label text and marker border color. Default \code{"black"}.
#' @param highlight_seqs Named list where each name is a border color and each element is a character
#'   vector of sequences to highlight by rings. Default \code{NULL}.
#' @param highlight_size Numeric; size of the highlight rings (defaults to \code{point_size + 0.6} if \code{NULL}).
#' @param highlight_stroke Numeric; thickness of the highlight rings. Default \code{1.2}.
#' @param y_log_scale Logical; if \code{TRUE}, use a log10 axis on \eqn{-log10(p.adj)} with tick labels
#'   shown in the original \eqn{-log10(p.adj)} scale. Default \code{FALSE}.
#'
#' @return A \emph{named list} of \code{ggplot} objects, one per requested comparison.
#'   Names are the same as selected \code{ttest_result} element names.
#'
#' @examples
#' \dontrun{
#' # Prepare data
#' result <- processPeptides(
#'   peptides_file          = "../Data/peptides.txt",
#'   intensity_columns_file = "../Data/Intensity_columns.csv",
#'   protein_mapping_file   = "../Data/protein_mapping.csv"
#' )
#' # suppose you ran ttestPeptides() like this with default test_method = "treat",
#'   lfc_thresh = 1 and alpha = 0.05:
#' ttest.1 <- ttestPeptides(
#'   result,
#'   comparisons = list(
#'     c("C40_G_N", "C40_I_N"),
#'     c("C40_G_N", "C80_I_N")
#'   )
#' )
#'
#' # 1) Plot by indices (index mode)
#' v1 <- plot_volcano(
#'   ttest_result = ttest.1,
#'   comparisons  = c(1, 2),                 # or list(1L, 2L)
#' )
#' v1$C40_G_N_vs_C80_I_N   # show one
#'
#' # 2) Plot by name pairs and try different colors
#' v2 <- plot_volcano(
#'  ttest_result = ttest.1,
#'  comparisons = list(
#'    c("C40_G_N", "C40_I_N"),
#'    c("C40_G_N", "C80_I_N")
#'  ),
#'  fill_values  = c(no = "grey80", yes = "#E41A1C")
#' )
#' v2$C40_G_N_vs_C40_I_N
#'
#' # 3) Label a few sequences and enable y-axis log scaling
#' some_peps <- c("AAEVIREA", "AAEVIREALQGITDPLFK")
#' v3 <- plot_volcano(
#'   ttest_result = ttest.1,
#'   comparisons  = 1,
#'   label_seqs   = some_peps,
#'   label_col    = "black",
#'   y_log_scale  = TRUE
#' )
#' v3$C40_G_N_vs_C40_I_N
#'
#' # 4) Highlight selected peptides with colored rings, set highlight size = 2 and stroke = 2
#' hl <- list(black = c("AAEVIREA"),
#'            "#377EB8" = c("AAEVIREALQGITDPLFK"))
#' v4 <- plot_volcano(
#' ttest_result   = ttest.1,
#'   comparisons    = 1,
#'   test_method    = "plain",
#'   highlight_seqs = hl,
#'   highlight_size = 2,
#'   highlight_stroke = 2
#' )
#' v4$C40_G_N_vs_C40_I_N
#' }
#'
#' # 5) suppose you ran ttestPeptides() not with the default settings:
#' my_lfc_thresh <- 2
#' my_alpha = 0.05
#' ttest.2 <- ttestPeptides(
#'   result,
#'   comparisons = list(
#'     c("C40_G_N", "C40_I_N"),
#'     c("C40_G_N", "C80_I_N")
#'   ),
#'   test_method  = "plain",
#'   lfc_thresh = my_lfc_thresh,
#'   alpha      = my_alpha
#' )
#' v5 <- plot_volcano(
#'   ttest_result   = ttest.2,
#'   comparisons    = c(1,2),
#'   test_method    = "plain",
#'   lfc_thresh     = my_lfc_thresh,
#'   alpha          = my_alpha
#' )
#' v5$C40_G_N_vs_C40_I_N
#'
#' @seealso \code{\link{ttestPeptides}}}
#'
#' @import ggplot2
#' @import data.table
#' @import ggpubr
#' @importFrom ggrepel geom_text_repel
#' @export

plot_volcano <- function(
    ttest_result,
    comparisons,
    test_method        = c("treat","plain"),
    show_threshold     = TRUE,
    lfc_thresh         = 1,
    alpha              = 0.05,
    fill_values        = c(no = "grey75", yes = "#FFC010"),
    point_size         = 2,
    point_alpha        = 0.85,
    label_seqs         = NULL,
    label_size         = 3,
    label_col          = "black",
    highlight_seqs     = NULL,
    highlight_size     = NULL,
    highlight_stroke   = 1.2,
    y_log_scale        = FALSE) {
  stopifnot(is.list(ttest_result), length(ttest_result) >= 1)
  test_method <- match.arg(test_method)

  # ------- helper: are we dealing with index-based selection? -------
  is_indices <- FALSE
  if (is.numeric(comparisons)) {
    is_indices <- TRUE
    idx <- as.integer(comparisons)
  } else if (is.list(comparisons) && length(comparisons) > 0 &&
             all(vapply(comparisons, function(x) is.numeric(x) && length(x) == 1, logical(1)))) {
    is_indices <- TRUE
    idx <- as.integer(unlist(comparisons, use.names = FALSE))

  } else if (is.character(comparisons) && length(comparisons) == 2L) { #named-pair selection
    pair_list <-list(comparisons)
  } else if (is.list(comparisons) && length(comparisons) > 0 &&
             all(vapply(comparisons, function(x) is.character(x) && length(x) == 2L, logical(1)))) {
    pair_list <- comparisons
  } else {
    stop("`comparisons` must be either numeric indices (single/int vector/list-of-ints) ",
         "or character pairs (c('A','B') or list of such vectors).")
  }

  # ------- case (a): index-based selection -------
  if (is_indices) {
    if (length(idx) == 0) stop("Empty numeric `comparisons` provided.")
    if (any(idx < 1L | idx > length(ttest_result))) {
      bad <- idx[idx < 1L | idx > length(ttest_result)]
      stop("Index/indices out of range in `comparisons`: ", paste(unique(bad), collapse = ", "))
    }

    out_names <- names(ttest_result)
    chosen_names <- out_names[idx]

    res_list <- vector("list", length(idx))
    names(res_list) <- chosen_names
  } else {
    # ------- case (b): named-pair selection -------
    # Normalize to list of length-2 character vectors
    comp_names <- vapply(pair_list, function(ab) paste0(ab[1], "_vs_", ab[2]), character(1))

    if (!all(comp_names %in% names(ttest_result))){
      bad <- comp_names[!comp_names %in% names(ttest_result)]
      stop("Comparison name(s) not found in `ttest_result`: ",
           paste(bad, collapse = ", "), ". Available: ",
           paste(names(ttest_result), collapse = ", "))
    }

    res_list   <- vector("list", length(pair_list))
    names(res_list) <- comp_names

  }

  # plot each comparison
  for (nm in names(res_list)) {
    dt <- ttest_result[[nm]]
    if (is.null(dt)) stop("`ttest_result[['", nm, "']]` is NULL.")

    # Build y variable
    padj <- dt[["p.adj"]]
    # Guard against 0 p-values (replace with 0.1 * smallest positive)
    if (any(padj == 0, na.rm = TRUE)) {
      min_pos <- suppressWarnings(min(padj[padj > 0], na.rm = TRUE))
      if (!is.finite(min_pos)) min_pos <- .Machine$double.xmin
      padj <- ifelse(padj == 0, min_pos * 0.1, padj)
    }
    dt[, yvar := -log10(padj)]
    y_lab <- expression(-log[10]("adjusted p-value"))

    ## Option to log-scale the axis of yvar (i.e., log10(-log10(p)))
    ## Set y_log_scale = TRUE/FALSE in your function args
    if (isTRUE(y_log_scale)) {
      # log10 axis can’t show zeros; push any non-positive y slightly above 0
      y_min_pos <- suppressWarnings(min(dt$yvar[dt$yvar > 0], na.rm = TRUE))
      if (!is.finite(y_min_pos)) y_min_pos <- 1e-3
      eps <- y_min_pos * 0.1
      dt[, yvar_plot := ifelse(yvar > 0, yvar, eps)]

      # decade breaks across observed range of yvar
      rng <- range(dt$yvar_plot, na.rm = TRUE)
      e_min <- floor(log10(rng[1]))
      e_max <- ceiling(log10(rng[2]))
      y_breaks <- 10^(e_min:e_max)

      y_scale <- ggplot2::scale_y_continuous(
        trans  = "log10",
        breaks = y_breaks,
        labels = function(b) formatC(b, format = "fg", digits = 1, flag = "#")
      )
    } else {
      dt[, yvar_plot := yvar]
      y_scale <- ggplot2::scale_y_continuous()
    }

    # factor significance
    dt[, sig := factor(sig, levels = c("yes","no"))]

    # Base plot: shape 21 to allow filled interior + separate border color
    p <- ggplot(dt, aes(x = .data[["log2FC"]], y = .data[["yvar_plot"]])) +
      geom_point(aes(fill = .data[["sig"]], colour = .data[["sig"]]),
                 shape = 21,
                 size  = point_size,
                 alpha = point_alpha)+
      scale_fill_manual(values = fill_values, name = "Significant") +
      scale_color_manual(values = fill_values, name = "Significant") +
      xlab("log2 fold change") + ylab(y_lab) +
      y_scale +
      theme_pubr()+
      theme(legend.position = "bottom")

    # Threshold lines
    if (isTRUE(show_threshold)) {
      # Horizontal p-value line
      y_thr <- -log10(alpha)
      p <- p + geom_hline(yintercept = y_thr, linetype = "dashed", linewidth = 0.5)
      # Vertical LFC lines for plain t.test only
      if (test_method == "plain") {
        p <- p + geom_vline(xintercept = c(-lfc_thresh, lfc_thresh),
                            linetype = "dashed", linewidth = 0.5)
      }
    }

    # Labels (avoid overlap)
    if (!is.null(label_seqs)) {
      lab_dt <- dt[Sequence %in% label_seqs]
      if (nrow(lab_dt)) {
        p <- p + ggrepel::geom_text_repel(
          data = lab_dt,
          aes(label = Sequence),
          size = label_size,
          col = label_col,
          max.overlaps = Inf,
          min.segment.length = 0
        )+ geom_point(
          data   = lab_dt,
          shape  = 21,
          size   = point_size,
          fill   = NA,
          color  = label_col
        )

      }
    }

    # Highlights: colored borders for selected sequences (rings)
    if (is.null(highlight_size)) highlight_size <- point_size + 0.6
    if (!is.null(highlight_seqs) && length(highlight_seqs)) {
      stopifnot(is.list(highlight_seqs), !is.null(names(highlight_seqs)))
      for (col_name in names(highlight_seqs)) {
        seqs <- highlight_seqs[[col_name]]
        hdt  <- dt[Sequence %in% seqs]
        if (nrow(hdt)) {
          p <- p + geom_point(
            data   = hdt,
            shape  = 21,
            size   = highlight_size,
            stroke = highlight_stroke,
            fill   = NA,
            color  = col_name
          )
        }
      }
    }

    res_list[[nm]] <- p
  }

  res_list
}

