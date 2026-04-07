# Protein quantification from TMT MS data -------------------------------------
#
# Description:
#   Processes TMT (Tandem Mass Tag) proteomics data:
#     1. In silico trypsin digestion from FASTA files
#     2. Computes Intensity, RI (razor-peptide intensity), and iBAQ/riBAQ
#        quantifications per sample channel
#     3. Visualises protein-length and peptide-length distributions
#
# Requirements:
#   library(Biostrings)
#   library(cleaver)
#   library(ggdist)
#   library(dplyr)
#   library(tidyr)
#   library(purrr)
#   library(ggplot2)

library(Biostrings)
library(cleaver)
library(ggdist)

# ── In silico trypsin digest ───────────────────────────────────────────────────

#' Filter Peptides to a Configured Length Range
#'
#' @param peps            Tibble with at least a `peptide` column.
#' @param min_peptide_len Integer. Minimum peptide length (inclusive).
#'   Default 6.
#' @param max_peptide_len Integer. Maximum peptide length (inclusive).
#'   Default 30.
#'
#' @return Filtered tibble.
filter_peptides <- function(peps, min_peptide_len = 6, max_peptide_len = 30) {
  lengths <- nchar(peps$peptide)
  peps[lengths >= min_peptide_len & lengths <= max_peptide_len, ]
}

#' Run In Silico Trypsin Digestion on One or More FASTA Files
#'
#' @description
#' Digests proteins using a single enzyme configuration, filters peptides by
#' length, writes comparison and summary tables to `output_path`, and returns
#' both the raw digest tibble (all peptides) and observable peptide counts
#' (filtered peptides per protein).
#'
#' @param fasta_paths      Character vector. Paths to FASTA files.
#' @param output_path      Character. Directory for summary TSV outputs.
#' @param enzym            Character. One of `"trypsin"`, `"trypsin-high"`,
#'   or `"trypsin-low"`. Default `"trypsin"`.
#' @param missed_cleavages Integer or integer vector. Allowed missed cleavages.
#'   Default 0.
#' @param min_peptide_len  Integer. Minimum peptide length. Default 6.
#' @param max_peptide_len  Integer. Maximum peptide length. Default 30.
#'
#' @return Named list:
#'   - `all_raw_digests`     : tibble of all (unfiltered) peptides.
#'   - `observable_peptides` : tibble of filtered peptide counts per protein.
process_tryptic <- function(
  fasta_paths,
  output_path,
  enzym             = "trypsin",
  missed_cleavages  = 0,
  min_peptide_len   = 6,
  max_peptide_len   = 30
) {
  all_fasta_ids        <- character(0)
  all_raw_digests      <- tibble::tibble()
  all_filtered_digests <- tibble::tibble()

  for (fasta_path in fasta_paths) {
    message("\nProcessing: ", fasta_path)

    species   <- tools::file_path_sans_ext(
      stringr::str_extract(basename(fasta_path), "^[A-Za-z]+_[A-Za-z]+")
    )
    fasta     <- Biostrings::readAAStringSet(fasta_path)
    fasta_ids <- sub(" .*", "", names(fasta))

    duplicated_ids <- intersect(fasta_ids, all_fasta_ids)
    if (length(duplicated_ids) > 0) {
      warning(
        length(duplicated_ids), " protein ID(s) in ", basename(fasta_path),
        " already present in a previous FASTA file.\n",
        "Duplicated IDs: ", paste(utils::head(duplicated_ids, 5), collapse = ", "),
        if (length(duplicated_ids) > 5) " ..." else ""
      )
    }
    all_fasta_ids <- c(all_fasta_ids, fasta_ids)

    peptides_per_protein <- cleaver::cleave(
      fasta,
      enzym            = enzym,
      missedCleavages  = missed_cleavages
    )

    raw_digests <- purrr::map_dfr(names(peptides_per_protein), function(protein_name) {
      tibble::tibble(
        Species    = species,
        Protein_id = sub(" .*", "", protein_name),
        peptide    = as.character(peptides_per_protein[[protein_name]])
      )
    })

    filtered_digests     <- filter_peptides(raw_digests, min_peptide_len, max_peptide_len)
    all_raw_digests      <- dplyr::bind_rows(all_raw_digests,      raw_digests)
    all_filtered_digests <- dplyr::bind_rows(all_filtered_digests, filtered_digests)
  }

  # Peptide filtering comparison table — one row per species
  comparison_table <- all_raw_digests |>
    dplyr::count(Species, name = "n_before") |>
    dplyr::left_join(
      dplyr::count(all_filtered_digests, Species, name = "n_after"),
      by = "Species"
    ) |>
    dplyr::mutate(
      n_removed = n_before - n_after,
      pct_kept  = round(100 * n_after / n_before, 1)
    )

  message("\nPeptide filtering summary (length ", min_peptide_len, "-", max_peptide_len, " aa):")
  out_compare <- file.path(output_path, "peptide_filter_comparison.tsv")
  utils::write.table(comparison_table, out_compare, sep = "\t", row.names = FALSE, quote = FALSE)

  # Observable peptide counts per protein (after filtering)
  observable_peptides <- all_filtered_digests |>
    dplyr::count(Species, Protein_id, name = "n_peptides")

  summary_table <- observable_peptides |>
    dplyr::group_by(Species) |>
    dplyr::summarise(
      Min    = min(n_peptides),
      Q1     = stats::quantile(n_peptides, 0.25),
      Median = stats::median(n_peptides),
      Mean   = round(mean(n_peptides), 2),
      Q3     = stats::quantile(n_peptides, 0.75),
      Max    = max(n_peptides),
      .groups = "drop"
    )

  message("\nSummary of observable peptide counts after filtering:")
  out_summary <- file.path(output_path, "peptide_summary.tsv")
  utils::write.table(summary_table, out_summary, sep = "\t", row.names = FALSE, quote = FALSE)

  list(all_raw_digests = all_raw_digests, observable_peptides = observable_peptides)
}

# ── Visualisation: protein and peptide length distributions ───────────────────

#' Plot Protein Length Distribution (Raincloud)
#'
#' @description
#' Raincloud plot (half-violin + jitter + boxplot) of protein length
#' per species, optionally log10-transformed. Median values are annotated
#' inside the boxplot.
#'
#' @param fasta_lengths_df Tibble with columns: `Species`, `protein_length`.
#'   Build this before calling `process_tryptic()`:
#'   ```r
#'   fasta_lengths_df <- purrr::map_dfr(fasta_files, function(fp) {
#'     sp <- tools::file_path_sans_ext(
#'       stringr::str_extract(basename(fp), "^[A-Za-z]+_[A-Za-z]+"))
#'     fa <- Biostrings::readAAStringSet(fp)
#'     tibble::tibble(Species = sp, protein_length = Biostrings::width(fa))
#'   })
#'   ```
#' @param log_scale    Logical. If TRUE (default), uses a log10 y-axis.
#' @param point_alpha  Numeric. Jitter point transparency. Default 0.25.
#' @param point_size   Numeric. Jitter point size. Default 0.6.
#' @param jitter_width Numeric. Jitter horizontal spread. Default 0.08.
#' @param title        Character or NULL. Optional plot title. Default NULL.
#'
#' @return ggplot object.
plot_protein_length <- function(
  fasta_lengths_df,
  log_scale    = TRUE,
  point_alpha  = 0.25,
  point_size   = 0.6,
  jitter_width = 0.08,
  title        = NULL
) {
  species_lvls <- sort(unique(fasta_lengths_df$Species))
  n_sp  <- length(species_lvls)
  pal   <- setNames(species_colours[seq_len(n_sp)], species_lvls)
  x_labels <- setNames(
    lapply(species_lvls, format_species_title, linebreak = TRUE),
    species_lvls
  )

  df <- fasta_lengths_df |>
    dplyr::filter(!is.na(protein_length), protein_length > 10) |>
    dplyr::mutate(Species = factor(Species, levels = species_lvls))

  medians <- df |>
    dplyr::group_by(Species) |>
    dplyr::summarise(med = stats::median(protein_length), .groups = "drop")

  ggplot2::ggplot(df, ggplot2::aes(x = Species, y = protein_length, fill = Species, colour = Species)) +
    ggdist::stat_halfeye(
      ggplot2::aes(fill = Species),
      adjust        = 0.8,
      width         = 0.5,
      .width        = 0,
      point_colour  = NA,
      alpha         = 0.55,
      justification = -0.25,
      normalize     = "groups"
    ) +
    ggplot2::geom_jitter(
      width  = jitter_width,
      height = 0,
      size   = point_size,
      alpha  = point_alpha
    ) +
    ggplot2::geom_boxplot(
      ggplot2::aes(fill = Species),
      width         = 0.12,
      outlier.shape = NA,
      alpha         = 0.7,
      colour        = "white",
      linewidth     = 0.5
    ) +
    ggplot2::geom_text(
      data        = medians,
      ggplot2::aes(x = Species, y = med, label = scales::comma(round(med))),
      vjust       = 0.8,
      hjust       = 1.3,
      size        = 3.5,
      colour      = "white",
      inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::scale_colour_manual(values = pal) +
    ggplot2::scale_x_discrete(labels = x_labels) +
    ggplot2::scale_y_continuous(
      labels = scales::comma,
      trans  = if (log_scale) "log10" else "identity",
      name   = if (log_scale) "Protein length (aa, log10)" else "Protein length (aa)"
    ) +
    ggplot2::coord_flip() +
    ggplot2::labs(title = title, x = NULL) +
    theme_publication(base_size = 32) +
    ggplot2::theme(legend.position = "none")
}

#' Plot Peptide Length Distribution After In Silico Digestion
#'
#' @description
#' Histogram of all peptide lengths from `all_raw_digests`, with a shaded
#' band marking the retained length window and bars coloured by retained /
#' filtered status.
#'
#' @param all_raw_digests  Tibble. Raw digest output from `process_tryptic()`
#'   (the `all_raw_digests` element). Must have a `peptide` column.
#' @param min_peptide_len  Integer. Lower bound of the retention window.
#'   Default 6.
#' @param max_peptide_len  Integer. Upper bound of the retention window.
#'   Default 30.
#'
#' @return ggplot object.
plot_peptide_length <- function(
  all_raw_digests,
  min_peptide_len = 6,
  max_peptide_len = 30
) {
  df <- all_raw_digests |>
    dplyr::mutate(
      pep_len  = nchar(peptide),
      filtered = pep_len >= min_peptide_len & pep_len <= max_peptide_len
    )

  n_total  <- nrow(df)
  n_kept   <- sum(df$filtered)
  pct_kept <- round(100 * n_kept / n_total, 1)

  # Cap display at slightly beyond the filter window
  x_max_disp <- max_peptide_len + 15L

  counts <- df |>
    dplyr::filter(pep_len <= x_max_disp) |>
    dplyr::count(pep_len, filtered)

  window_df <- data.frame(
    xmin = min_peptide_len - 0.5,
    xmax = max_peptide_len + 0.5
  )

  ggplot2::ggplot(counts, ggplot2::aes(x = pep_len, y = n, fill = filtered)) +
    ggplot2::geom_rect(
      data        = window_df,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill        = "#E1F5EE",
      alpha       = 0.55
    ) +
    ggplot2::geom_col(width = 0.85, alpha = 0.85) +
    ggplot2::geom_vline(
      xintercept = min_peptide_len - 0.5,
      colour     = "#0F6E56", linewidth = 0.5, linetype = "dashed"
    ) +
    ggplot2::geom_vline(
      xintercept = max_peptide_len + 0.5,
      colour     = "#0F6E56", linewidth = 0.5, linetype = "dashed"
    ) +
    ggplot2::annotate(
      "text",
      x      = (min_peptide_len + max_peptide_len) / 2,
      y      = Inf,
      label  = sprintf(
        "Retained window\n%d\u2013%d aa  (%s%%)",
        min_peptide_len, max_peptide_len, pct_kept
      ),
      vjust  = 1.25, size = 3, colour = "#0F6E56", family = "sans"
    ) +
    ggplot2::scale_fill_manual(
      values = c("TRUE" = "#1D9E75", "FALSE" = "#B4B2A9"),
      labels = c("TRUE" = "Retained", "FALSE" = "Filtered out"),
      name   = NULL
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(0, x_max_disp, by = 5),
      limits = c(0, x_max_disp + 1),
      name   = "Peptide length (aa)"
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::comma,
      expand = ggplot2::expansion(mult = c(0, 0.12)),
      name   = "Number of peptides"
    ) +
    ggplot2::labs(
      subtitle = sprintf(
        "In silico digest  \u2022  total = %s peptides  \u2022  %s%% retained after length filter",
        scales::comma(n_total), pct_kept
      )
    ) +
    theme_publication(base_size = 32) +
    ggplot2::theme(
      legend.position    = c(0.82, 0.88),
      legend.key.size    = ggplot2::unit(0.55, "lines"),
      legend.background  = ggplot2::element_rect(fill = "white", colour = NA)
    )
}

# ── Quantification methods ─────────────────────────────────────────────────────

#' Compute Intensity Quantification per Protein per Sample Channel
#'
#' @description
#' Proportionally distributes total protein intensity across TMT channels,
#' then adds log2, log-ratio, normalised, and parts-per-billion columns.
#' Rows with zero or missing intensity are removed.
#'
#' @param tmt_data Tibble. Unified TMT data with `Species`, `Replicate`,
#'   `Treatment`, and `channel_*` columns. Requires `TMT_label` and
#'   `COL_TOTAL_INTENSITY` to be defined in scope.
#'
#' @return Tibble with additional columns: `Intensity`, `Intensity_log2`,
#'   `Intensity_meanlog2`, `Intensity_log2ratio`, `Intensity_norm`,
#'   `Intensity_log`, `Intensity_PpB`.
intensity_quantification <- function(tmt_data) {
  message("Computing intensity quantification for each sample...")

  channel_cols <- grep("^channel_", colnames(tmt_data), value = TRUE)

  tmt_data |>
    dplyr::mutate(
      channel_sum = rowSums(dplyr::across(dplyr::all_of(channel_cols)), na.rm = TRUE)
    ) |>
    tidyr::pivot_longer(
      cols      = dplyr::all_of(channel_cols),
      names_to  = "channel_col",
      values_to = "channel_intensity"
    ) |>
    dplyr::filter(channel_col == paste0("channel_", TMT_label)) |>
    dplyr::mutate(
      Intensity = .data[[COL_TOTAL_INTENSITY]] * channel_intensity / channel_sum
    ) |>
    dplyr::filter(dplyr::if_all(Intensity, ~ . > 0 & !is.na(.))) |>
    dplyr::group_by(Species, Protein_id, Treatment) |>
    dplyr::mutate(
      Intensity_log2      = log2(Intensity),
      Intensity_meanlog2  = mean(Intensity_log2, na.rm = TRUE),
      Intensity_log2ratio = Intensity_log2 - Intensity_meanlog2,
      Intensity_norm      = Intensity / sum(Intensity, na.rm = TRUE),
      Intensity_log       = 10 + log10(Intensity_norm),
      Intensity_PpB       = Intensity_norm * 1e9
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-channel_sum, -channel_intensity)
}

#' Compute Razor-Peptide Relative Intensity (RI) per Sample Channel
#'
#' @description
#' Divides total intensity by razor-peptide count (RI), then proportionally
#' distributes across TMT channels and adds log-transformed columns.
#' Rows with zero or missing RI are removed.
#'
#' @param tmt_data Tibble. See `intensity_quantification()` for requirements.
#'   Additionally requires `COL_RAZOR_PEPTIDES` to be defined in scope.
#'
#' @return Tibble with additional columns: `RI`, `RI_log2`, `RI_meanlog2`,
#'   `RI_log2ratio`, `RI_norm`, `RI_log`, `RI_PpB`.
ri_quantification <- function(tmt_data) {
  message("Computing relative intensity quantification for each sample...")

  channel_cols <- grep("^channel_", colnames(tmt_data), value = TRUE)

  tmt_data |>
    dplyr::mutate(
      channel_sum = rowSums(dplyr::across(dplyr::all_of(channel_cols)), na.rm = TRUE),
      RI_total    = .data[[COL_TOTAL_INTENSITY]] / .data[[COL_RAZOR_PEPTIDES]]
    ) |>
    tidyr::pivot_longer(
      cols      = dplyr::all_of(channel_cols),
      names_to  = "channel_col",
      values_to = "channel_intensity"
    ) |>
    dplyr::filter(channel_col == paste0("channel_", TMT_label)) |>
    dplyr::mutate(RI = RI_total * (channel_intensity / channel_sum)) |>
    dplyr::filter(dplyr::if_all(RI, ~ . > 0 & !is.na(.))) |>
    dplyr::group_by(Species, Protein_id, Treatment) |>
    dplyr::mutate(
      RI_log2      = log2(RI),
      RI_meanlog2  = mean(RI_log2, na.rm = TRUE),
      RI_log2ratio = RI_log2 - RI_meanlog2,
      RI_norm      = RI / sum(RI, na.rm = TRUE),
      RI_log       = 10 + log10(RI_norm),
      RI_PpB       = RI_norm * 1e9
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-channel_sum, -channel_intensity)
}

#' Compute iBAQ and riBAQ Intensities per Protein per Sample Channel
#'
#' @description
#' Joins observable (tryptic) peptide counts from `process_tryptic()` and
#' divides total intensity by the number of observable peptides (iBAQ).
#' The result is then normalised within each sample (riBAQ) and log-transformed.
#' Proteins without a FASTA match receive `iBAQ = NA` (and are optionally
#' written to disk).
#'
#' @param tmt_data            Tibble. Unified TMT data with `Species`,
#'   `Protein_id`, `channel_*`, and `total_intensity` columns.
#' @param observable_peptides Tibble. Output of `process_tryptic()` containing
#'   `Species`, `Protein_id`, and `n_peptides` (or the column named by
#'   `observable_col`).
#' @param observable_col      Character. Column in `observable_peptides` to use
#'   as denominator. Default `"n_tryptic"`.
#' @param output_path         Character or NULL. If provided, unmatched proteins
#'   are written to `iBAQ_unmatched_proteins.tsv` in this directory.
#'
#' @return Tibble with additional columns: `iBAQ`, `iBAQ_log2`,
#'   `iBAQ_meanlog2`, `iBAQ_log2ratio`, `riBAQ`, `iBAQ_log`, `iBAQ_PpB`.
ibaq_quantification <- function(
  tmt_data,
  observable_peptides,
  observable_col = "n_tryptic",
  output_path    = NULL
) {
  message("Computing iBAQ for each sample...")

  channel_cols <- grep("^channel_", colnames(tmt_data), value = TRUE)

  tmt_joined <- tmt_data |>
    dplyr::left_join(observable_peptides, by = c("Species", "Protein_id")) |>
    dplyr::mutate(unmatched = is.na(.data[[observable_col]]))

  unmatched_proteins <- dplyr::distinct(dplyr::filter(tmt_joined, unmatched), Species, Protein_id)
  n_unmatched        <- nrow(unmatched_proteins)

  if (n_unmatched > 0) {
    warning(n_unmatched, " proteins had no match in FASTA and will have iBAQ = NA")
    if (!is.null(output_path)) {
      unmatched_path <- file.path(output_path, "iBAQ_unmatched_proteins.tsv")
      utils::write.table(
        dplyr::arrange(unmatched_proteins, Species, Protein_id),
        unmatched_path, sep = "\t", row.names = FALSE, quote = FALSE
      )
      message("Unmatched proteins written to: ", unmatched_path)
    }
  }

  tmt_joined |>
    dplyr::mutate(
      iBAQ_total  = total_intensity / .data[[observable_col]],
      channel_sum = rowSums(dplyr::across(dplyr::all_of(channel_cols)), na.rm = TRUE)
    ) |>
    tidyr::pivot_longer(
      cols      = dplyr::all_of(channel_cols),
      names_to  = "channel_col",
      values_to = "channel_intensity"
    ) |>
    dplyr::filter(channel_col == paste0("channel_", TMT_label)) |>
    dplyr::mutate(iBAQ = iBAQ_total * (channel_intensity / channel_sum)) |>
    dplyr::group_by(Species, Protein_id, Treatment) |>
    dplyr::mutate(
      iBAQ_log2      = log2(iBAQ + 1),
      iBAQ_meanlog2  = mean(iBAQ_log2, na.rm = TRUE),
      iBAQ_log2ratio = log2(iBAQ) - iBAQ_meanlog2,
      riBAQ          = iBAQ / sum(iBAQ, na.rm = TRUE),
      iBAQ_log       = 10 + log10(riBAQ),
      iBAQ_PpB       = riBAQ * 1e9
    ) |>
    dplyr::ungroup()
}

# ── Replicate correlation plots ────────────────────────────────────────────────

#' Build One Pairwise Replicate Scatter Panel (Internal)
#'
#' @description
#' Scatter plot of log2 iBAQ for one replicate pair within one treatment,
#' with a Pearson R annotation and an identity diagonal.
#'
#' @param df_wide  Tibble. Wide-format data (one column per sample).
#' @param sample_x Character. Column name for the x-axis sample.
#' @param sample_y Character. Column name for the y-axis sample.
#' @param colour   Character. Hex fill/colour for points.
#'
#' @return ggplot object.
.one_replicate_panel <- function(df_wide, sample_x, sample_y, colour) {
  d <- df_wide |>
    dplyr::select(x = dplyr::all_of(sample_x), y = dplyr::all_of(sample_y)) |>
    dplyr::filter(is.finite(x), is.finite(y))

  if (nrow(d) < 3) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5,
                          label = "insufficient data", colour = "#888780") +
        ggplot2::theme_void()
    )
  }

  ct    <- stats::cor.test(d$x, d$y, method = "pearson")
  R_val <- round(ct$estimate, 2)
  p_val <- ct$p.value
  p_fmt <- if (p_val < 0.001) "0.00e+00" else sprintf("%.3f", p_val)

  label_txt <- sprintf("R: %.2f\np-val: %s", R_val, p_fmt)

  lo    <- floor(min(c(d$x, d$y)))
  hi    <- ceiling(max(c(d$x, d$y)))
  ann_x <- lo + 0.03 * (hi - lo)
  ann_y <- hi - 0.03 * (hi - lo)

  ggplot2::ggplot(d, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         colour = "#cccccc", linewidth = 0.35, linetype = "solid") +
    ggplot2::geom_point(colour = colour, alpha = 0.35, size = 0.55, shape = 16) +
    ggplot2::annotate(
      "text",
      x          = ann_x, y = ann_y,
      label      = label_txt,
      hjust      = 0, vjust = 1,
      size       = 10,
      colour     = "#2c2c2a",
      lineheight = 0.35,
      family     = "Arial"
    ) +
    ggplot2::scale_x_continuous(limits = c(lo, hi), breaks = scales::pretty_breaks(n = 4)) +
    ggplot2::scale_y_continuous(limits = c(lo, hi), breaks = scales::pretty_breaks(n = 4)) +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      title = sprintf("%s vs %s", sample_x, sample_y),
      x     = sample_x,
      y     = sample_y
    ) +
    theme_publication(base_size = 40)
}

#' Plot All Pairwise Replicate Correlations for One Species
#'
#' @description
#' Arranges pairwise scatter panels (log2 iBAQ) in a grid: treatments as rows,
#' replicate pairs as columns. Returns a patchwork figure.
#'
#' @param df           Tibble. Data for one species with columns: `Species`,
#'   `Protein_id`, `Treatment`, `Replicate`, `iBAQ_log2`.
#' @param species_name Character. Species label for the figure title.
#' @param sample_col   Character or NULL. If not NULL, this column is used as
#'   the sample label instead of constructing `Treatment_Replicate`.
#' @param colour_map   Named character vector or NULL. Maps treatment name to
#'   hex colour. If NULL, colours are auto-assigned from `treatment_colours`.
#'
#' @return patchwork ggplot object (or NULL invisibly if no panels generated).
plot_replicate_correlations <- function(
  df,
  species_name,
  sample_col   = NULL,
  colour_map   = NULL
) {
  if (!is.null(sample_col)) {
    df <- dplyr::rename(df, .sample = dplyr::all_of(sample_col))
  } else {
    df <- dplyr::mutate(df, .sample = paste0(Treatment, "_", Replicate))
  }

  treatments <- sort(unique(df$Treatment))
  n_tr       <- length(treatments)

  if (is.null(colour_map)) {
    pal        <- treatment_colours[((seq_len(n_tr) - 1L) %% length(treatment_colours)) + 1L]
    colour_map <- setNames(pal, treatments)
  }

  wide <- df |>
    dplyr::select(Protein_id, .sample, iBAQ_log2) |>
    tidyr::pivot_wider(
      names_from  = .sample,
      values_from = iBAQ_log2,
      values_fn   = mean    # collapse any duplicates
    )

  all_rows <- purrr::map(treatments, function(trt) {
    trt_df  <- dplyr::filter(df, Treatment == trt)
    samples <- sort(unique(trt_df$.sample))

    if (length(samples) < 2) {
      message("  Skipping ", trt, ": fewer than 2 replicates.")
      return(NULL)
    }

    pairs  <- utils::combn(samples, 2, simplify = FALSE)
    panels <- purrr::map(pairs, function(pr) {
      .one_replicate_panel(wide, pr[1], pr[2], colour_map[[trt]])
    })

    patchwork::wrap_plots(panels, ncol = 3)
  })

  all_rows <- purrr::compact(all_rows)
  if (length(all_rows) == 0) {
    warning("No panels generated for: ", species_name)
    return(invisible(NULL))
  }

  fig <- patchwork::wrap_plots(all_rows, ncol = 2) &
    ggplot2::theme(plot.margin = ggplot2::margin(2, 4, 2, 4))

  fig +
    patchwork::plot_annotation(
      title = species_name,
      theme = theme_publication(base_size = 32) +
        ggplot2::theme(
          plot.title       = ggplot2::element_text(face = "bold", colour = "#2c2c2a"),
          plot.caption     = ggplot2::element_text(colour = "#888780"),
          plot.background  = ggplot2::element_rect(fill = "white", colour = NA)
        )
    )
}
