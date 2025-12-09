#' Plot peptide alignment along a protein sequence
#'
#' @description
#' Visualize peptides mapped to a protein sequence
#' bar on top and stacked peptide rectangles below. Peptides are stacked to avoid
#' overlaps, lettered, colored by intensity, and grouped by remaining grouping
#' variables (and \code{Replicate} for replicate-level plots). Protein sequence is
#' looked up in \code{peptidomics_protein_mapping_example} when available, or
#' can be supplied manually via \code{protein_seq}.
#'
#' @param result List. Output of \code{processPeptides()} or \code{filterPeptides()},
#'   containing at least \code{dt.peptides.int}, \code{dt.peptides.int.reps},
#'   and \code{grp_cols}.
#' @param protein_name Character(1). Leading razor protein identifier to plot.
#'   Must be present in the selected result table list.
#' @param protein_seq Character(1), or \code{NULL}. Full protein sequence. If
#'   \code{NULL}, the function attempts to find the sequence in
#'   \code{peptidomics_protein_mapping_example}; when not found, a warning is
#'   emitted and the x-axis range falls back to the observed peptide endpoints.
#' @param type Character. Which table to plot: \code{"mean"} for group means or
#'   \code{"reps"} for replicate-level peptides. Default: \code{"mean"}.
#' @param filter_params Named list, or \code{NULL}.  Each element’s name is a
#'   grouping column (from \code{result$grp_cols}, and optionally
#'   \code{"Replicate"} when \code{type = "reps"}), and its value is a vector of
#'   values to include. Multiple names impose an AND filter. Default: \code{NULL}.
#' @param x_interval Numeric. Spacing between x-axis breaks (protein positions).
#'   Default: \code{25}.
#' @param seq_col Character. Column containing peptide sequences. Default:
#'   \code{"Sequence"}.
#' @param start_col Character. Column containing peptide start positions.
#'   Default: \code{"Start.position"}.
#' @param end_col Character. Column containing peptide end positions.
#'   Default: \code{"End.position"}.
#' @param x_range Numeric vector of length 2, or \code{NULL}. Custom x-axis
#'   limits. Default: \code{NULL} (computed from the protein length or peptide
#'   endpoints).
#' @param y_range Numeric vector of length 2, or \code{NULL}. Custom y-axis
#'   limits for the peptide panel (before reversal). Default: \code{NULL}
#'   (automatically sized to the number of stacked rows).
#' @param intensity_col Character(1), or \code{NULL}. Column for intensity-based
#'   fill. Defaults to \code{"Mean.Intensity"} for mean plots and
#'   \code{"Intensity"} for replicate plots.
#' @param label_seq Named list, or \code{NULL}. Optional sequence highlighting:
#'   each element’s name is a label and its value are these
#'   sequences to outline. Defaults to \code{NULL} (no labels).
#' @param label_col Named list, or \code{NULL}. Optional colors for each label
#'   name in \code{label_seq}. Any labels not listed here use internal defaults.
#'
#' @return A combined \code{ggplot} object with protein and peptide panels.
#'
#' @examples
#' \dontrun{
#' result <- processPeptides(
#'   peptides_file          = "data/Yogurtexample_QR188-205.csv",
#'   intensity_columns_file = "data/Intensity_columns.csv",
#'   protein_mapping_file   = "data/protein_mapping.csv"
#' )
#'
#' # 1) Minimal call — sequence auto-looked-up if present in the packaged mapping
#' plot_pep_align(result, protein_name = "P02662")
#'
#' # 2) Provide the protein sequence manually
#' data(peptidomics_protein_mapping_example)
#' alpha_seq <- peptidomics_protein_mapping_example[
#'   peptidomics_protein_mapping_example[["Leading.razor.protein"]] == "P02662",
#'   "Protein.seq"
#' ]
#' plot_pep_align(
#'   result,
#'   protein_name = "P02662",
#'   protein_seq  = alpha_seq
#' )
#'
#' # 3) Mean (default) vs. replicate view, with filters on grouping variables
#' plot_pep_align(
#'   result,
#'   protein_name  = "P02662",
#'   filter_params = list(Digest.stage = "G120")
#' )
#' plot_pep_align(
#'   result,
#'   protein_name  = "P02662",
#'   type          = "reps",
#'   filter_params = list(Yogurt = "Y1")
#' )
#'
#' # 4) Plot detail options: tick spacing and axis ranges
#' plot_pep_align(
#'   result,
#'   protein_name = "P02662",
#'   x_interval   = 20,
#'   x_range      = c(0, 260),
#'   y_range      = c(0, 12)
#' )
#'
#' # 5) Custom column names (if your tables differ from MaxQuant defaults)
#' plot_pep_align(
#'   result,
#'   protein_name = "P02662",
#'   seq_col      = "Sequence",
#'   start_col    = "Start.position",
#'   end_col      = "End.position"
#' )
#'
#' # 6) Highlight specific sequences with optional custom colors
#' plot_pep_align(
#'   result,
#'   protein_name = "P02662",
#'   label_seq    = list(highlight_1 = c("DQAMEDIKQ", "DQAMEDIKQM"),
#'                       highlight_2 = c("AMEDIKQM")),
#'   label_col    = list(highlight_1 = "red",
#'                      highlight_2 = "blue")
#' )
#' }
#'
#' @import ggplot2
#' @import data.table
#' @importFrom ggpubr theme_pubr
#' @importFrom utils data
#' @importFrom patchwork wrap_plots plot_spacer
#' @export
plot_pep_align <- function(result,
                           protein_name,
                           protein_seq = NULL,
                           type = c("mean", "reps"),
                           filter_params = NULL,
                           x_interval = 25,
                           seq_col = "Sequence",
                           start_col = "Start.position",
                           end_col = "End.position",
                           x_range = NULL,
                           y_range = NULL,
                           intensity_col = NULL,
                           label_seq = NULL,
                           label_col = NULL) {
  type <- match.arg(type)

  dt <- switch(
    type,
    mean = data.table::copy(result$dt.peptides.int),
    reps = data.table::copy(result$dt.peptides.int.reps)
  )
  grp_cols <- result$grp_cols

  if (is.null(intensity_col)) {
    intensity_col <- if (type == "mean") "Mean.Intensity" else "Intensity"
  }

  required_cols <- c("Leading.razor.protein", seq_col, start_col, end_col, intensity_col)
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols)) {
    stop("Missing required column(s): ", paste(missing_cols, collapse = ", "))
  }

  # apply filters on grouping variables
  if (!is.null(filter_params)) {
    stopifnot(is.list(filter_params))
    allowed_filters <- grp_cols
    if (type == "reps") {
      allowed_filters <- c(allowed_filters, "Replicate")
    }
    invalid_filters <- setdiff(names(filter_params), allowed_filters)
    if (length(invalid_filters)) {
      stop("Filter column(s) not found in the result: ",
           paste(invalid_filters, collapse = ", "))
    }
    for (col in names(filter_params)) {
      dt <- dt[get(col) %in% filter_params[[col]]]
    }
  }

  dt <- dt[Leading.razor.protein == protein_name]
  dt <- dt[!is.na(get(start_col)) & !is.na(get(end_col))]
  if (nrow(dt) == 0) {
    stop("No peptides found for protein_name after filtering.")
  }

  # look up protein sequence if not provided
  if (is.null(protein_seq)) {
    env <- new.env(parent = emptyenv())
    data_loaded <- try(utils::data("peptidomics_protein_mapping_example",
                                   package = "peptidomicsR",
                                   envir = env),
                       silent = TRUE)
    if (inherits(data_loaded, "try-error") ||
        !exists("peptidomics_protein_mapping_example", envir = env)) {
      data_loaded <- try(utils::data("peptidomics_protein_mapping_example",
                                     envir = env),
                         silent = TRUE)
    }
    if ((inherits(data_loaded, "try-error") ||
         !exists("peptidomics_protein_mapping_example", envir = env)) &&
        file.exists("data/peptidomics_protein_mapping_example.rda")) {
      try(load("data/peptidomics_protein_mapping_example.rda", envir = env),
          silent = TRUE)
    }
    if (exists("peptidomics_protein_mapping_example", envir = env)) {
      mapping_dt <- data.table::as.data.table(
        get("peptidomics_protein_mapping_example", envir = env)
      )
      seq_hit <- mapping_dt[Leading.razor.protein == protein_name, Protein.seq]
      if (length(seq_hit) && !is.na(seq_hit[1])) {
        protein_seq <- seq_hit[1]
      }
    }
    if (is.null(protein_seq)) {
      warning("Protein ", protein_name,
              " not found in peptidomics_protein_mapping_example; ",
              "using peptide endpoints for axis scaling. ",
              "Provide protein_seq to silence this warning.")
    }
  }

  protein_length <- if (!is.null(protein_seq)) {
    nchar(protein_seq)
  } else {
    max(dt[[end_col]], na.rm = TRUE)
  }
  if (is.na(protein_length) || protein_length <= 0) {
    stop("Unable to determine protein length; please provide protein_seq.")
  }

  # derive panel labels from remaining grouping variables
  facet_vars <- grp_cols
  if (type == "reps") {
    facet_vars <- c(facet_vars, "Replicate")
  }
  if (!is.null(filter_params)) {
    facet_vars <- setdiff(facet_vars, names(filter_params))
  }
  facet_vars <- facet_vars[
    vapply(facet_vars,
           function(v) length(unique(dt[[v]])) > 1,
           logical(1))
  ]
  if (length(facet_vars) == 0) {
    dt[, sample_panel := "Sample"]
  } else {
    dt[, sample_panel := interaction(.SD, sep = " | ", drop = TRUE),
       .SDcols = facet_vars]
  }
  dt[, sample_panel := factor(sample_panel, levels = unique(sample_panel))]

  # ensure numeric for stacking inputs
  dt[[start_col]] <- as.numeric(dt[[start_col]])
  dt[[end_col]]   <- as.numeric(dt[[end_col]])
  dt[[intensity_col]] <- as.numeric(dt[[intensity_col]])

  # drop non-positive intensities for log scale safety
  dt <- dt[!is.na(get(intensity_col)) & get(intensity_col) > 0]

  # assign stacking rows per sample
  split_samples <- split(dt, dt$sample_panel)
  stacked_list <- lapply(split_samples, function(sdt) {
    sdt <- sdt[order(get(start_col), -get(end_col))]
    n <- nrow(sdt)
    if (n == 0) return(NULL)

    # stacking groups
    stack_group <- integer(n)
    stack_group[1] <- 1L
    current_group <- 1L
    current_end <- sdt[[end_col]][1]
    if (n > 1) {
      for (i in 2:n) {
        if (sdt[[start_col]][i] <= current_end + 1) {
          stack_group[i] <- current_group
          current_end <- max(current_end, sdt[[end_col]][i])
        } else {
          current_group <- current_group + 1L
          stack_group[i] <- current_group
          current_end <- sdt[[end_col]][i]
        }
      }
    }
    sdt$stack_group <- stack_group

    # stacking rows within each group
    sdt$stack_row <- NA_integer_
    for (g in sort(unique(stack_group))) {
      idx <- which(stack_group == g)
      gdt <- sdt[idx, ]
      gdt <- gdt[order(gdt[[start_col]])]
      row_assign <- integer(nrow(gdt))
      row_ends <- numeric(0)
      for (j in seq_len(nrow(gdt))) {
        st <- gdt[[start_col]][j]
        en <- gdt[[end_col]][j]
        eligible <- which(st > (row_ends + 1))
        if (length(eligible)) {
          chosen <- eligible[which.max(row_ends[eligible])]
          row_assign[j] <- chosen
          row_ends[chosen] <- en
        } else {
          row_assign[j] <- length(row_ends) + 1
          row_ends <- c(row_ends, en)
        }
      }
      sdt$stack_row[idx] <- row_assign
    }
    sdt
  })
  dt_stacked <- data.table::rbindlist(stacked_list, fill = TRUE)

  # sequence highlighting with optional custom colors and legend entries
  dt_stacked[, `:=`(label_color = NA_character_, label_name = NA_character_)]
  label_colors <- character(0)
  if (!is.null(label_seq)) {
    stopifnot(is.list(label_seq))
    default_cols <- c("black", "#4daf4a", "#984ea3",
                      "#a65628", "#f781bf")
    lab_names <- names(label_seq)
    if (is.null(lab_names) || any(lab_names == "")) {
      stop("label_seq elements must be named.")
    }
    for (i in seq_along(lab_names)) {
      lbl <- lab_names[i]
      seqs_lbl <- label_seq[[i]]
      col_lbl <- if (!is.null(label_col) && !is.null(label_col[[lbl]])) {
        label_col[[lbl]]
      } else {
        default_cols[(i - 1) %% length(default_cols) + 1]
      }
      label_colors[lbl] <- col_lbl
      dt_stacked[get(seq_col) %in% seqs_lbl,
                 `:=`(
                   label_color = col_lbl,
                   label_name = lbl
                 )]
    }
  }

  # rectangle boundaries for peptides
  dt_stacked[, rect_xmax := get(end_col) + 1]

  # letter data for peptide sequences
  letter_dt <- data.table::rbindlist(lapply(seq_len(nrow(dt_stacked)), function(i) {
    row <- dt_stacked[i]
    seq_txt <- as.character(row[[seq_col]])
    if (is.na(seq_txt) || nchar(seq_txt) == 0) return(NULL)
    L <- nchar(seq_txt)
    data.table::data.table(
      x = row[[start_col]] + (seq_len(L)) - 0.5,
      y = row$stack_row,
      letter = substring(seq_txt, seq_len(L), seq_len(L)),
      sample_panel = row$sample_panel
    )
  }))

  # protein sequence plot (top bar)
  aa_letters <- if (!is.null(protein_seq)) {
    strsplit(protein_seq, "")[[1]]
  } else {
    rep("X", protein_length)
  }
  protein_df <- data.table::data.table(
    start = seq_len(protein_length),
    end   = seq_len(protein_length) + 1,
    x     = seq_len(protein_length) + 0.5,
    letter = aa_letters
  )
  protein_plot <- ggplot(protein_df) +
    geom_rect(aes(xmin = start, xmax = end, ymin = 0.5, ymax = 1.5),
              fill = "gray95", color = "gray30", linewidth = 0.1) +
    geom_text(aes(x = x, y = 1, label = letter),
              size = 3, vjust = 0.5, color = "black") +
    scale_x_continuous(
      limits = c(1, protein_length + 1),
      breaks = c(1, seq(x_interval, protein_length + x_interval, by = x_interval)),
      expand = expansion(mult = c(0, 0.01)),
      position = "top"
    ) +
    coord_cartesian(xlim = if (is.null(x_range)) NULL else x_range) +
    labs(x = "Protein Position", y = NULL, tag = "P") +
    theme_pubr() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.tag.position = c(0, 0.2),
      plot.tag = element_text(angle = 90, vjust = 0, hjust = 0.5, size = 10),
      plot.margin = margin(b = 0.1)
    )

  # per-panel peptide plots with proportional heights
  xlim_use <- if (!is.null(x_range)) x_range else NULL
  panel_list <- split(dt_stacked, dt_stacked$sample_panel)
  intensity_range <- range(dt_stacked[[intensity_col]], na.rm = TRUE)
  panel_heights <- vapply(panel_list, function(sdt) {
    if (is.null(y_range)) max(sdt$stack_row, na.rm = TRUE) else abs(y_range[1] - y_range[2]) + 1
  }, numeric(1))
  panel_heights <- pmax(1, panel_heights)

  max_panel_h <- max(panel_heights)

  build_panel <- function(sdt, h) {
    panel_id <- unique(sdt$sample_panel)
    letters_sub <- letter_dt[sample_panel == panel_id]
    ymax_use <- if (is.null(y_range)) max(sdt$stack_row, na.rm = TRUE) + 0.5 else rev(y_range)[1]
    ymin_use <- if (is.null(y_range)) -0.5 else rev(y_range)[2]
    break_end <- if (!is.null(xlim_use)) {
      max(xlim_use, na.rm = TRUE)
    } else {
      max(protein_length + 1, max(sdt[[end_col]] + 1, na.rm = TRUE))
    }

    p <- ggplot(sdt,
                aes(xmin = .data[[start_col]],
                    xmax = rect_xmax,
                    ymin = stack_row - 0.45,
                    ymax = stack_row + 0.45,
                    fill = .data[[intensity_col]])) +
      geom_rect(color = NA) +
      geom_text(data = letters_sub,
                aes(x = x, y = y, label = letter, fill = NULL),
                inherit.aes = FALSE,
                size = 3,
                vjust = 0.5,
                color = "black") +
      scale_x_continuous(
        limits = if (is.null(xlim_use)) c(1, break_end) else xlim_use,
        breaks = c(1, seq(x_interval, break_end + x_interval, by = x_interval)),
        expand = expansion(mult = c(0, 0.01))
      ) +
      scale_y_reverse() +
      scale_fill_distiller(
        palette = "RdYlBu",
        direction = -1,
        trans = "log10",
        labels = scientific_10,
        limits = intensity_range,
        guide = "colourbar"
      ) +
      coord_cartesian(ylim = c(ymax_use, ymin_use)) +
      labs(x = NULL, y = NULL, fill = intensity_col, tag = panel_id) +
      theme_pubr() +
      theme(
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_blank(),
        plot.tag.position = c(0, 0.5),
        plot.tag = element_text(angle = 90, vjust = 0, hjust = 0.5, size = 10),
        legend.position = "bottom",
      legend.key.width = grid::unit(2, "line"),
      legend.key.height = grid::unit(0.5, "line"),
      plot.margin = margin(b = 0, t = 0)
      )
    if (length(label_colors) && any(!is.na(sdt$label_name))) {
      label_rects <- sdt[!is.na(label_name)]
      dummy_rects <- data.table::data.table(
        label_name = names(label_colors),
        xmin = 1, xmax = 1,
        ymin = ymin_use, ymax = ymin_use
      )
      p <- p +
        geom_rect(
          data = label_rects,
          aes(colour = label_name),
          fill = NA,
          linewidth = 2,
          show.legend = FALSE
        ) +
        geom_rect(
          data = dummy_rects,
          aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, colour = label_name),
          inherit.aes = FALSE,
          fill = NA,
          alpha = 0,
          linewidth = 0,
          show.legend = TRUE
        ) +
        scale_colour_manual(
          values = label_colors,
          drop = TRUE,
          na.translate = FALSE,
          guide = guide_legend(title = "Label", override.aes = list(fill = NA))
        )
    } else {
      p <- p + scale_colour_identity(guide = "none")
    }
    p
  }

  panel_plots <- mapply(function(sdt, h) {
    p <- build_panel(sdt, h)
    p
  }, panel_list, panel_heights, SIMPLIFY = FALSE)

  p_pep <- patchwork::wrap_plots(panel_plots, ncol = 1, heights = panel_heights)

  # combine panels with collected legend in its own row
  combined <- patchwork::wrap_plots(
    list(protein_plot, p_pep, patchwork::plot_spacer()),
    ncol = 1,
    heights = c(1, sum(panel_heights), 0.1)
  ) +
    patchwork::plot_layout(guides = "collect") &
    theme(
      legend.position = "bottom",
      legend.margin = margin(t = 0)
    )

  combined
}
