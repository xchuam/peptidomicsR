#' Plot N/C-terminal peptide cleavage site sequqnce logos by group
#'
#' @description
#' Visualize amino-acid preferences at peptide cleavage sites using sequence logos,
#' aggregated from the \code{processPeptides()} result. You can plot the N- or
#' C-terminal residue distribution (or both), measure either total intensity or
#' peptide counts, average replicates or show each replicate separately, and
#' optionally drop grouping columns that are constant after filtering so the
#' x-axis only reflects variables that truly vary.
#'
#' @param result List. Output of \code{processPeptides()}. Must contain:
#'   \itemize{
#'     \item \code{dt.peptides.int.reps}: data.table with per-replicate peptide intensities
#'           (columns include grouping variables in \code{grp_cols}, \code{Replicate},
#'           \code{Intensity}, \code{First.amino.acid}, \code{Last.amino.acid}).
#'     \item \code{grp_cols}: character vector of grouping columns to define sample groups.
#'   }
#' @param terminal Character. Which terminal(s) to plot: one of \code{"both"}, \code{"N"}, or \code{"C"}. Default \code{"both"}.
#' @param measure Character. What to aggregate: \code{"intensity"} (sum of intensities per replicate) or \code{"count"}
#'   (number of peptides with non-zero intensity per replicate). Default \code{"intensity"}.
#' @param replicate_mode Character. How to handle replicates:
#'   \code{"mean"} (average replicate values within each group) or
#'   \code{"reps"} (treat each replicate as its own x-axis level).
#'   Default: \code{"mean"}.
#' @param filter_params Named list, or \code{NULL}.  Each element’s name is a grouping column,
#'   and its value is a vector of values to include.  Multiple names impose an AND filter.
#'   For example: \code{list(Lipid = c("N","S"), Digest.stage = "G")}
#'   Default: \code{NULL} (no filtering).
#' @param scientific_10_y Logical.  If \code{TRUE}, use scientific notation for y-axis.
#'   Default: \code{TRUE}.
#' @param drop_constant_groups Logical. If \code{TRUE}, any grouping columns from
#'   \code{result$grp_cols} (and \code{Replicate} when \code{replicate_mode = "reps"})
#'   that have only one remaining level after filtering are dropped from the x-axis
#'   grouping. Default: \code{TRUE}.
#'
#' @return A \code{ggplot} object when \code{terminal} is \code{"N"} or \code{"C"};
#'   otherwise a combined grob (two panels) when \code{terminal = "both"}.
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
#' # 1) N-terminal, mean over replicates, intensity
#' p1 <- plot_cleavage_site(
#'   result,
#'   terminal       = "N",
#'   measure        = "intensity",
#'   replicate_mode = "mean"
#' )
#'
#' # 2) C-terminal, show each replicate separately on x (labels include Replicate)
#' p2 <- plot_cleavage_site(
#'   result,
#'   terminal       = "C",
#'   measure        = "count",
#'   replicate_mode = "reps"
#' )
#'
#' # 3) Both terminals side by side; filter and drop constant grouping columns on x
#' p3 <- plot_cleavage_site(
#'   result,
#'   terminal            = "both",
#'   measure             = "intensity",
#'   replicate_mode      = "mean",
#'   filter_params       = list(Lipid = "N", Digest.stage = "G"),
#'   drop_constant_groups = TRUE
#' )
#'
#' # 4) Keep all grouping columns on x even if constant; single x label when all collapse
#' p4 <- plot_cleavage_site(
#'   result,
#'   terminal              = "N",
#'   measure               = "intensity",
#'   replicate_mode        = "mean",
#'   filter_params         = list(Lipid = "N", Digest.stage = "G"),
#'   drop_constant_groups  = FALSE
#' )
#' }
#'
#' @import data.table
#' @import ggplot2
#' @import gridExtra
#' @importFrom cowplot get_legend
#' @importFrom ggseqlogo ggseqlogo
#' @export
plot_cleavage_site <- function(result,
                               terminal = c("both","N","C"),
                               measure = c("intensity","count"),
                               replicate_mode = c("mean","reps"),
                               filter_params = NULL,
                               scientific_10_y   = TRUE,
                               drop_constant_groups = TRUE) {
  terminal <- match.arg(terminal)
  measure  <- match.arg(measure)
  replicate_mode <- match.arg(replicate_mode)

  # ---- fetch table from result ----
  grp_cols <- result$grp_cols
  dt <- copy(result$dt.peptides.int.reps)

  # apply filters
  if (!is.null(filter_params)) {
    stopifnot(is.list(filter_params))
    for (col in names(filter_params)) {
      dt <- dt[get(col) %in% filter_params[[col]]]
    }
  }

  # when plotting per-replicate, Replicate can also be a grouping column
  if (replicate_mode == "reps") grp_cols <- c(grp_cols, "Replicate")

  # choose which grouping columns still vary (post-filter)
  if (drop_constant_groups) {
    varies <- vapply(grp_cols, function(g) uniqueN(dt[[g]]) > 1L, logical(1))
    active_grp_cols <- grp_cols[varies]
  } else {
    active_grp_cols <- grp_cols
  }

  # helper to aggregate per terminal
  make_mat <- function(side = c("N","C")) {
    side <- match.arg(side)
    aa_col <- if (side == "N") "First.amino.acid" else "Last.amino.acid"

    # per-replicate aggregate by AA and group columns
    if (measure == "intensity") {
      dt_rep <- dt[, .(value = sum(Intensity, na.rm = TRUE)), by = c(aa_col, unique(active_grp_cols, "Replicate"))]
    } else { # count
      dt_rep <- dt[Intensity > 0 & !is.na(Intensity), .N, by = c(aa_col, unique(active_grp_cols, "Replicate"))]
      setnames(dt_rep, "N", "value")
    }

    # if replicate_mode == "mean", average across replicate (becasue "Replicate" in not in active_grp_cols)
    dt_group <- dt_rep[, .(value = mean(value, na.rm = TRUE)), by = c(aa_col, active_grp_cols)]

    # Cast to wide: rows = AA, cols = interaction of active_grp_cols
    rhs <- paste(active_grp_cols, collapse = "+")
    mat_dt <- dcast(dt_group,
                    formula = as.formula(paste(aa_col, "~", rhs)),
                    value.var = "value",
                    fill = 0)

    setnames(mat_dt, aa_col, "AA")
    setkey(mat_dt, "AA")

    rn <- mat_dt$AA
    mat <- as.matrix(mat_dt[, -1, with = FALSE])
    rownames(mat) <- rn

    mat
  }

  # Build requested plot(s)
  make_logo <- function(mat, side_label) {
    p <- suppressMessages(
      ggseqlogo(mat, method = 'custom', seq_type = 'aa') +
      scale_x_continuous(breaks = seq_len(ncol(mat)), labels = colnames(mat)) +
      labs(x = "Sample group", y = paste0(side_label, " terminal ",
                                          if (measure == "intensity") "Intensity" else "Count")) +
      theme_classic() +
      theme(axis.line = element_line(linewidth = 0.2, colour = "gray"),
            axis.text.x = element_text(angle = 45, hjust = 1))
    )
    # optional log scale
    if (scientific_10_y) {
      p <- p + scale_y_continuous(
        labels = scientific_10
      )
    }
    p
  }

  ## helper to build a combined legend from N & C
  make_combined_legend <- function(matN,
                                   matC) {
    # AA present in either plot (non-zero rows)
    aas <- union(rownames(matN),rownames(matC))

    # Dummy PWM with one column; each present AA gets a positive weight
    pwm <- matrix(1, nrow = length(aas), ncol = 1,
                  dimnames = list(aas, "legend"))

    p_leg_sr <- make_logo(pwm,"N")+theme(legend.position  = "bottom")

    get_legend(p_leg_sr)
  }

  if (terminal == "N") {
    matN <- make_mat("N")
    return(make_logo(matN, "N"))
  } else if (terminal == "C") {
    matC <- make_mat("C")
    return(make_logo(matC, "C"))
  } else {
    matN <- make_mat("N")
    matC <- make_mat("C")

    pN <- make_logo(matN, "N")+
      theme(legend.position = "none")
    pC <- make_logo(matC, "C")+
      theme(legend.position = "none")

    p_legend <- make_combined_legend(matN,
                                     matC)

    # arrange side-by-side
    grid.arrange(
      arrangeGrob(pN,pC, ncol=2),
      p_legend,
      heights = c(10,1)
    )

  }
}

