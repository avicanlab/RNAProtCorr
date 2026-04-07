# Essential genes analyses -----------------------------------------------------
#
# Description:
#   Functions for comparing expression distributions and mRNA-protein
#   correlations between essential/stimulon and non-essential genes.
#
# Requirements:
#   library(ggplot2)
#   library(patchwork)
#   library(dplyr)

# ── Constants ──────────────────────────────────────────────────────────────────

#' Mapping from essential-gene file abbreviations to full species names
SPECIES_ESSENTIAL_TO_TPM_MAP <- c(
  "Senterica"           = "Salmonella_enterica",
  "Saureus"             = "Staphylococcus_aureus",
  "Ypseudotuberculosis" = "Yersinia_pseudotuberculosis"
)

# ── Internal helpers ───────────────────────────────────────────────────────────

#' Annotate a Data Frame With Essential / Non-Essential Group Labels
#'
#' @param all_df       Tibble. All genes with at least `Species` and `Protein_id`.
#' @param essential_df Tibble. Essential gene IDs per species (`Species`,
#'   `Protein_id`).
#' @param value_col    Character. Column in `all_df` used for filtering finite
#'   values.
#'
#' @return Tibble with a `group` factor column (`"Essential"` /
#'   `"Non-Essential"`).
.annotate_essential_group <- function(all_df, essential_df, value_col) {
  all_df |>
    dplyr::left_join(
      dplyr::mutate(essential_df, group = "Essential"),
      by = c("Species", "Protein_id")
    ) |>
    dplyr::mutate(
      group = factor(
        tidyr::replace_na(group, "Non-Essential"),
        levels = c("Non-Essential", "Essential")
      )
    ) |>
    dplyr::filter(is.finite(.data[[value_col]]))
}

#' Build the x-axis label for essential-gene distribution plots
.essential_x_label <- function(x_lab, type = c("mean", "sd")) {
  type <- match.arg(type)
  if (type == "mean") {
    bquote(.(ifelse(x_lab == "TPM", "mRNA", "Protein")) ~
             "- Mean " ~ .(x_lab) ~ (Log[2]))
  } else {
    if (x_lab == "TPM") {
      bquote(atop("Standard deviation of mRNA -", log[2](mean ~ TPM)))
    } else {
      bquote(atop("Standard deviation of protein -", log[2](mean ~ .(x_lab))))
    }
  }
}

# ── Distribution plots ─────────────────────────────────────────────────────────

#' Plot Essential vs Non-Essential Expression Distribution per Treatment
#'
#' @description
#' Overlapping frequency histograms comparing log2 mean expression between
#' essential and non-essential genes, faceted by treatment condition within
#' each species. A two-sided Wilcoxon p-value is annotated on each panel.
#'
#' @param all_df       Tibble. All genes with a global log2 mean column.
#' @param essential_df Tibble. Essential gene IDs per species.
#' @param log2_col     Character. Name of the log2 mean column in `all_df`.
#' @param x_lab        Character. Short label for the data type (e.g. `"TPM"`,
#'   `"iBAQ"`).
#' @param output_path  Character. Directory for saving output files.
#'
#' @return A patchwork plot of all species combined (invisible).
plot_essential_distribution_per_treatment <- function(
  all_df, essential_df, log2_col, x_lab, output_path
) {
  condition_df <- .annotate_essential_group(all_df, essential_df, log2_col)
  x_label      <- .essential_x_label(x_lab, type = "mean")

  sp_plots <- condition_df |>
    dplyr::group_by(Species) |>
    dplyr::group_map(function(sp_df, sp_key) {
      species   <- sp_key$Species
      title_sp  <- format_species_title(species)
      x_lim     <- range(sp_df[[log2_col]], na.rm = TRUE)

      treatment_plots <- sp_df |>
        dplyr::group_by(Treatment) |>
        dplyr::group_map(function(treat_df, treat_key) {
          treatment <- treat_key$Treatment
          p_value   <- stats::wilcox.test(
            stats::as.formula(paste(log2_col, "~ group")),
            data = treat_df
          )$p.value

          ggplot2::ggplot(
            treat_df,
            ggplot2::aes(x = .data[[log2_col]], y = ggplot2::after_stat(density), fill = group)
          ) +
            ggplot2::geom_histogram(position = "identity", alpha = 0.5, bins = 40) +
            ggplot2::scale_x_continuous(limits = x_lim) +
            ggplot2::scale_fill_manual(
              values = c("Non-Essential" = "#F4A8A0", "Essential" = "#80CEC8")
            ) +
            ggplot2::annotate(
              "text",
              x = Inf, y = Inf,
              label  = format.pval(p_value, digits = 2, eps = .Machine$double.xmin),
              hjust  = 1.1, vjust = 1.5, size = 6, fontface = "italic"
            ) +
            ggplot2::labs(
              title = treatment,
              x     = x_label,
              y     = "Frequency",
              fill  = NULL
            ) +
            theme_publication() +
            ggplot2::theme(
              plot.title    = ggplot2::element_text(face = "italic", hjust = 0.5),
              legend.position = "none"
            )
        })

      sp_plot <- patchwork::wrap_plots(treatment_plots, ncol = 3) +
        patchwork::plot_annotation(
          title = title_sp,
          theme = ggplot2::theme(
            plot.title = ggplot2::element_text(
              face   = "bold.italic",
              hjust  = 0.5,
              margin = ggplot2::margin(b = 20)
            )
          )
        ) +
        patchwork::plot_layout(guides = "collect", axes = "collect") &
        ggplot2::theme(legend.position = "right")

      save_plot(
        sp_plot,
        file.path(output_path, species, paste(x_lab, "per_treatment_essential_distribution", sep = "_")),
        width  = 21,
        height = 24
      )
      sp_plot
    })

  patchwork::wrap_plots(sp_plots, ncol = 1) +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "right")
}

#' Plot Essential vs Non-Essential Expression Distribution (Global)
#'
#' @description
#' Overlapping frequency histograms (one panel per species) comparing log2
#' mean expression between essential and non-essential genes. A two-sided
#' Wilcoxon p-value is annotated on each panel.
#'
#' @param all_df       Tibble. All genes with a global log2 mean column.
#' @param essential_df Tibble. Essential gene IDs per species.
#' @param log2_col     Character. Name of the log2 mean column in `all_df`.
#' @param x_lab        Character. Short label for the data type.
#' @param output_path  Character. Directory for saving output files.
#'
#' @return A patchwork plot of all species combined (invisible).
plot_essential_distribution <- function(
  all_df, essential_df, log2_col, x_lab, output_path
) {
  condition_df <- .annotate_essential_group(all_df, essential_df, log2_col)
  x_label      <- .essential_x_label(x_lab, type = "mean")

  plots <- condition_df |>
    dplyr::group_by(Species) |>
    dplyr::group_map(function(sp_df, sp_key) {
      species  <- sp_key$Species
      p_value  <- stats::wilcox.test(
        stats::as.formula(paste(log2_col, "~ group")),
        data = sp_df
      )$p.value
      title_sp <- format_species_title(species)

      p <- ggplot2::ggplot(
        sp_df,
        ggplot2::aes(x = .data[[log2_col]], y = ggplot2::after_stat(density), fill = group)
      ) +
        ggplot2::geom_histogram(position = "identity", alpha = 0.5, bins = 40) +
        ggplot2::scale_fill_manual(
          values = c("Non-Essential" = "#F4A8A0", "Essential" = "#80CEC8")
        ) +
        ggplot2::annotate(
          "text",
          x = Inf, y = Inf,
          label  = format.pval(p_value, digits = 2, eps = .Machine$double.xmin),
          hjust  = 1.1, vjust = 1.5, size = 3.5, fontface = "italic"
        ) +
        ggplot2::labs(
          title = title_sp,
          x     = x_label,
          y     = "Frequency",
          fill  = NULL
        ) +
        theme_publication() +
        ggplot2::theme(legend.position = c(0.85, 0.85))

      save_plot(
        p,
        file.path(output_path, species, paste0(x_lab, "_essential_distribution")),
        width  = 12,
        height = 4
      )
      p
    })

  patchwork::wrap_plots(plots, nrow = 1) +
    patchwork::plot_layout(guides = "collect", axes = "collect_y") &
    ggplot2::theme(legend.position = "top")
}

#' Plot Essential vs Non-Essential SD Distribution
#'
#' @description
#' Frequency polygon comparing log2 expression standard deviation between
#' essential and non-essential genes, with a two-sided Wilcoxon p-value
#' annotated on each species panel.
#'
#' @param all_df       Tibble. All genes with an SD column (e.g. from
#'   `log2_mean_sd()`).
#' @param essential_df Tibble. Essential gene IDs per species.
#' @param sd_col       Character. Name of the SD column in `all_df`.
#' @param x_lab        Character. Short label for the data type.
#' @param output_path  Character. Directory for saving output files.
#'
#' @return A patchwork plot of all species combined (invisible).
plot_sd_distribution <- function(
  all_df, essential_df, sd_col, x_lab, output_path
) {
  condition_df <- .annotate_essential_group(all_df, essential_df, sd_col)
  x_label      <- .essential_x_label(x_lab, type = "sd")

  plots <- condition_df |>
    dplyr::group_by(Species) |>
    dplyr::group_map(function(sp_df, sp_key) {
      species  <- sp_key$Species
      p_value  <- stats::wilcox.test(
        stats::as.formula(paste(sd_col, "~ group")),
        data = sp_df
      )$p.value
      title_sp <- format_species_title(species)

      p <- ggplot2::ggplot(sp_df, ggplot2::aes(x = .data[[sd_col]], color = group)) +
        ggplot2::geom_step(
          ggplot2::aes(y = ggplot2::after_stat(density)),
          stat      = "bin",
          binwidth  = 0.05,
          linewidth = 0.8,
          direction = "mid"
        ) +
        ggplot2::scale_color_manual(
          values = c("Non-Essential" = "#F4A8A0", "Essential" = "#80CEC8")
        ) +
        ggplot2::annotate(
          "text",
          x = Inf, y = Inf,
          label    = sprintf("p-value: %.2e", p_value),
          hjust    = 1.1, vjust = 1.5, size = 8, fontface = "italic"
        ) +
        ggplot2::labs(
          title = title_sp,
          x     = x_label,
          y     = "Frequency",
          color = NULL
        ) +
        theme_publication(legend_position = c(0.75, 0.85)) +
        ggplot2::theme(
          plot.title        = ggplot2::element_text(hjust = 0.5),
          legend.background = ggplot2::element_rect(fill = "white", color = NA)
        )

      save_plot(
        p,
        file.path(output_path, species, paste0(x_lab, "_sd_distribution")),
        width  = 16,
        height = 8
      )
      p
    })

  patchwork::wrap_plots(plots, nrow = 1) +
    patchwork::plot_layout(guides = "collect", axes = "collect_y") &
    ggplot2::theme(legend.position = "top")
}

# ── Correlation comparison plots ───────────────────────────────────────────────

#' Plot Gene-Subset vs All-Gene mRNA-Protein Correlation
#'
#' @description
#' Scatter plot comparing per-treatment mRNA-protein Pearson R values
#' computed on all genes versus a specific gene subset (e.g. essential genes
#' or stimulon members). One panel per species. A diagonal reference line
#' (y = x) aids interpretation.
#'
#' This single generic function replaces the former separate
#' `plot_essential_vs_all_correlation()` and
#' `plot_stimulon_vs_all_correlation()` functions, which were identical
#' except for axis labels and file-name prefix.
#'
#' @param corr_all_df    Tibble. All-gene correlation results with columns
#'   `Species`, `Treatment`, `R`.
#' @param corr_subset_df Tibble. Subset-gene correlation results with the same
#'   column structure.
#' @param measurement    Character. Quantification label used in the output file
#'   name and y-axis (e.g. `"iBAQ"`).
#' @param subset_label   Character. Short label for the subset used in axis
#'   titles and the output file name (e.g. `"essential"`, `"stimulon"`).
#' @param output_path    Character. Directory for saving output files.
#'
#' @return List of ggplot objects (one per species), assembled with
#'   `patchwork::wrap_plots()`.
plot_subset_vs_all_correlation <- function(
  corr_all_df, corr_subset_df, measurement, subset_label, output_path
) {
  r_subset_col <- paste0("R_", subset_label)

  plot_df <- dplyr::inner_join(
    dplyr::rename(dplyr::distinct(corr_all_df,    Species, Treatment, R), R_all = R),
    dplyr::rename(dplyr::distinct(corr_subset_df, Species, Treatment, R), !!r_subset_col := R),
    by = c("Species", "Treatment")
  )

  plot_df |>
    dplyr::group_by(Species) |>
    dplyr::group_map(function(sp_df, sp_key) {
      species  <- sp_key$Species
      title_sp <- format_species_title(species)

      p <- ggplot2::ggplot(sp_df, ggplot2::aes(x = R_all, y = .data[[r_subset_col]], color = Treatment)) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                             color = "grey60", linewidth = 0.6) +
        ggplot2::geom_point(size = 4, alpha = 0.9) +
        ggplot2::scale_color_manual(
          values = setNames(scales::hue_pal()(nrow(sp_df)), sp_df$Treatment)
        ) +
        ggplot2::coord_fixed(ratio = 1, xlim = c(0.45, 0.9), ylim = c(0.45, 0.9)) +
        ggplot2::labs(
          title = title_sp,
          x     = "All genes'\nmRNA-protein level correlation",
          y     = paste0(stringr::str_to_title(subset_label), " genes'\nmRNA-protein level correlation"),
          color = "Treatment"
        ) +
        theme_publication(legend_position = "right") +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5),
          axis.title = ggplot2::element_text(lineheight = 0.5)
        )

      save_plot(
        p,
        file.path(output_path, species, paste0(measurement, "_", subset_label, "_vs_all_correlation")),
        width  = 5,
        height = 5
      )
      p
    })
}
