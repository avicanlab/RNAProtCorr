# CDS detection and TPM group distribution ------------------------------------
#
# Description:
#   Functions for reporting annotated CDS detection rates across RNA and
#   proteomics datasets, and for visualising TPM distributions split by
#   protein-detection status.
#
# Requirements:
#   library(dplyr)
#   library(tidyr)
#   library(ggplot2)
#   library(purrr)

# ── CDS detection summary ──────────────────────────────────────────────────────

#' Compute Detected CDS Counts and Percentages per Species
#'
#' @description
#' For each species in `annotation_data`, counts the number of detected CDS
#' (unique `Protein_id`) in RNA and proteomics datasets relative to the total
#' CDS in the annotation reference.
#'
#' @param annotation_data Tibble. Annotation reference with at least `Species`
#'   and `Protein_id` columns.
#' @param rna_data        Tibble. Detected genes with at least `Species` and
#'   `Protein_id` columns (typically long-format TPM data).
#' @param prot_data       Tibble. Detected proteins with at least `Species` and
#'   `Protein_id` columns (typically filtered proteomics data).
#'
#' @return Tibble with one row per species and columns:
#'   - `Species`       : species name
#'   - `total_cds`     : total annotated CDS in the reference
#'   - `rna_detected`  : CDS detected in RNA data
#'   - `prot_detected` : CDS detected in proteomics data
#'   - `pct_rna`       : percentage of CDS detected in RNA
#'   - `pct_prot`      : percentage of CDS detected in proteomics
get_detected_cds <- function(annotation_data, rna_data, prot_data) {
  purrr::map_dfr(unique(annotation_data$Species), function(sp) {
    total_cds    <- dplyr::filter(annotation_data, Species == sp) |>
      dplyr::distinct(Protein_id) |>
      nrow()
    rna_detected  <- dplyr::filter(rna_data,  Species == sp) |>
      dplyr::distinct(Protein_id) |>
      nrow()
    prot_detected <- dplyr::filter(prot_data, Species == sp) |>
      dplyr::distinct(Protein_id) |>
      nrow()

    tibble::tibble(
      Species       = sp,
      total_cds     = total_cds,
      rna_detected  = rna_detected,
      prot_detected = prot_detected,
      pct_rna       = 100 * rna_detected  / total_cds,
      pct_prot      = 100 * prot_detected / total_cds
    )
  })
}

#' Plot Percentage of Detected CDS per Species
#'
#' @description
#' Bar chart showing the percentage of annotated CDS detected in RNA and
#' proteomics data per species. Supports two positioning modes:
#'
#' - `"overlay"` (default): RNA bar is drawn first (wider, in the background);
#'   the protein bar overlaps it (narrower, in the foreground). This keeps both
#'   bars visible when protein coverage is substantially lower than RNA
#'   coverage.
#' - `"dodge"`: bars are placed side by side.
#'
#' @param detected_cds_data Tibble. Output of `get_detected_cds()`. Must
#'   contain `Species`, `pct_rna`, `pct_prot`.
#' @param position          Character. `"overlay"` or `"dodge"`. Default
#'   `"overlay"`.
#' @param show_percent      Logical. If TRUE, percentage values are annotated
#'   on the bars. Default FALSE.
#'
#' @return ggplot object.
plot_barplot_detected_cds <- function(
  detected_cds_data,
  position     = c("overlay", "dodge"),
  show_percent = FALSE
) {
  position <- match.arg(position)

  plot_df <- detected_cds_data |>
    tidyr::pivot_longer(
      cols      = c(pct_rna, pct_prot),
      names_to  = "type",
      values_to = "pct"
    ) |>
    dplyr::mutate(
      type    = factor(type,
                       levels = c("pct_rna", "pct_prot"),
                       labels = c("Transcriptome", "Proteome")),
      Species = factor(Species)
    )

  species_levels <- levels(plot_df$Species)
  x_labels <- setNames(
    lapply(species_levels, format_species_title, linebreak = TRUE),
    species_levels
  )

  fill_values <- c("Transcriptome" = "#C2E3F2", "Proteome" = "#FFD892")

  if (position == "overlay") {
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Species, y = pct, fill = type)) +
      ggplot2::geom_col(
        data     = ~ dplyr::filter(.x, type == "Transcriptome"),
        width    = 0.8, alpha = 1, position = "identity"
      ) +
      ggplot2::geom_col(
        data     = ~ dplyr::filter(.x, type == "Proteome"),
        width    = 0.8, alpha = 1, position = "identity"
      )

    if (show_percent) {
      p <- p +
        ggplot2::geom_text(
          data    = ~ dplyr::filter(.x, type == "Transcriptome"),
          ggplot2::aes(label = sprintf("%.1f%%", pct), y = pct),
          vjust   = -0.4, size = 3.5, color = "grey30"
        ) +
        ggplot2::geom_text(
          data    = ~ dplyr::filter(.x, type == "Proteome"),
          ggplot2::aes(label = sprintf("%.1f%%", pct), y = pct),
          vjust   = 1.4, size = 3.5, color = "grey30"
        )
    }
  } else {
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Species, y = pct, fill = type)) +
      ggplot2::geom_col(
        position = ggplot2::position_dodge(width = 0.75),
        width    = 0.65, alpha = 0.85
      )

    if (show_percent) {
      p <- p +
        ggplot2::geom_text(
          ggplot2::aes(label = sprintf("%.1f%%", pct)),
          position = ggplot2::position_dodge(width = 0.75),
          vjust    = -0.4, size = 3.5, color = "grey30"
        )
    }
  }

  p +
    ggplot2::scale_fill_manual(values = fill_values, name = NULL) +
    ggplot2::scale_x_discrete(labels = x_labels) +
    ggplot2::scale_y_continuous(
      limits = c(0, 100),
      expand = ggplot2::expansion(mult = c(0, 0.02)),
      breaks = seq(0, 100, by = 20)
    ) +
    ggplot2::labs(x = NULL, y = "CDS detected (%)") +
    theme_publication(legend_position = "top") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))
}

# ── TPM group classification ───────────────────────────────────────────────────

#' Classify Genes Into TPM Detection Groups for One Species
#'
#' @description
#' Splits genes into three groups:
#'   - `"mRNA"`        : all genes detected in TPM (full set, duplicated for
#'     overlay histogram)
#'   - `"Protein"`     : genes detected in both TPM and proteomics
#'   - `"Not protein"` : genes detected in TPM only (absent from proteomics)
#'
#' @param tpm_data  Tibble. Must contain `Species`, `Protein_id`, and
#'   `mean_log2_TPM`.
#' @param prot_data Tibble. Must contain `Species` and `Protein_id`.
#' @param sp        Character. Species name to filter on.
#'
#' @return Tibble with columns: `Protein_id`, `mean_log2_TPM`, `group`
#'   (factor with levels `"mRNA"`, `"Not protein"`, `"Protein"`).
classify_tpm_groups <- function(tpm_data, prot_data, sp) {
  tpm_sp  <- dplyr::filter(tpm_data,  Species == sp)
  prot_sp <- dplyr::filter(prot_data, Species == sp)

  common_ids <- extract_common_ids(tpm_sp, prot_sp)

  group_levels <- c("mRNA", "Not protein", "Protein")

  tpm_unique <- tpm_sp |>
    dplyr::filter(is.finite(mean_log2_TPM)) |>
    dplyr::group_by(Protein_id) |>
    dplyr::summarise(mean_log2_TPM = mean(mean_log2_TPM, na.rm = TRUE), .groups = "drop")

  # Overlay-style: duplicate the full set as "mRNA", then add specific groups
  dplyr::bind_rows(
    dplyr::mutate(
      tpm_unique,
      group = factor("mRNA", levels = group_levels)
    ),
    dplyr::mutate(
      tpm_unique,
      group = factor(
        ifelse(Protein_id %in% common_ids, "Protein", "Not protein"),
        levels = group_levels
      )
    )
  )
}

#' Plot TPM Distribution by Detection Group for One Species
#'
#' @description
#' Overlapping histogram with three layers (`mRNA`, `Not protein`, `Protein`)
#' drawn back-to-front so smaller distributions remain visible. The species
#' name is rendered as an italic panel title via `format_species_title()`.
#'
#' @param classified_df Tibble. Output of `classify_tpm_groups()` for one
#'   species.
#' @param species       Character. Species name for the italic panel title.
#' @param x_max         Numeric or NULL. Maximum x-axis value; inferred from
#'   data when NULL. Default NULL.
#'
#' @return ggplot object.
plot_tpm_detection_species <- function(classified_df, species, x_max = NULL) {
  if (is.null(x_max)) {
    x_max <- ceiling(max(classified_df$mean_log2_TPM, na.rm = TRUE))
  }

  title_expr <- format_species_title(species)

  layer_order  <- c("mRNA", "Not protein", "Protein")
  group_colours <- c(
    "mRNA"        = "#C2E3F2",
    "Not protein" = "#A3BECC",
    "Protein"     = "#FFD892"
  )

  classified_df |>
    dplyr::mutate(group = factor(group, levels = layer_order)) |>
    ggplot2::ggplot(ggplot2::aes(x = mean_log2_TPM, fill = group)) +
    ggplot2::geom_histogram(
      data     = ~ dplyr::filter(.x, group == "mRNA"),
      binwidth = 0.5, alpha = 0.85, colour = NA
    ) +
    ggplot2::geom_histogram(
      data     = ~ dplyr::filter(.x, group == "Not protein"),
      binwidth = 0.5, alpha = 0.85, colour = NA
    ) +
    ggplot2::geom_histogram(
      data     = ~ dplyr::filter(.x, group == "Protein"),
      binwidth = 0.5, alpha = 0.85, colour = NA
    ) +
    ggplot2::scale_fill_manual(values = group_colours, breaks = layer_order, name = NULL) +
    ggplot2::scale_x_continuous(
      limits = c(0, x_max),
      expand = ggplot2::expansion(mult = c(0, 0.02)),
      breaks = seq(0, x_max, by = 2.5)
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::labs(
      title = title_expr,
      x     = bquote("Mean TPM " ~ (Log[2])),
      y     = "Number of genes"
    ) +
    theme_publication() +
    ggplot2::theme(
      plot.title        = ggplot2::element_text(hjust = 0.5),
      panel.grid.major  = ggplot2::element_line(colour = "grey92"),
      legend.position   = c(0.72, 0.82),
      legend.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
}
