# ============================================================================
# Protein quantification from TMT MS data
# ============================================================================
#
# Description:
#   This script processes TMT (Tandem Mass Tag) proteomics data. It:
#     1. Filters TMT results per species and replicate
#     2. Performs in silico trypsin digestion from FASTA files
#     3. Computes relative intensity, iBAQ, and riBAQ quantifications
#
# Requirements:
#   install.packages("tidyverse")
#   install.packages("readxl")
#   if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#   BiocManager::install("Biostrings")
#   BiocManager::install("cleaver")
# ============================================================================

library(Biostrings)
library(cleaver)
library(purrr)

# ============================================================================
# CONFIGURATION
# ============================================================================

# Column names in TMT results files
col_protein <- "Protein"
col_protein_id <- "Protein ID"
col_entry_name <- "Entry Name"
col_gene <- "Gene"
col_organism <- "Organism"
col_unique_peptides <- "Unique Peptides"
col_razor_peptides <- "Razor Peptides"
col_total_intensity <- "Total Intensity"

keep_cols <- c(
  col_protein,
  col_protein_id,
  col_entry_name,
  col_gene,
  col_organism,
  col_unique_peptides,
  col_razor_peptides,
  col_total_intensity
)

# ============================================================================
# IN SILICO TRYPSIN DIGEST FROM FASTA
# ============================================================================
#' Filter peptides to those within the configured length range
#'
#' @param peps A tibble with at least a `peptide` column
#' @param min_peptide_len Minimum peptide length (default 6)
#' @param max_peptide_len Maximum peptide length (default 30)
#' @return A filtered tibble
filter_peptides <- function(peps, min_peptide_len = 6, max_peptide_len = 30) {
    lengths <- nchar(peps$peptide)
    peps[lengths >= min_peptide_len & lengths <= max_peptide_len, ]
}

#' Run in silico trypsin digestion on one or more FASTA files
#'
#' Digests proteins using a single enzyme configuration, filters peptides by
#' length, and writes summary tables to output_path.
#'
#' @param fasta_paths      A character vector of paths to FASTA files
#' @param output_path      Path to the output directory (must end with "/")
#' @param enzym            Enzyme to use for digestion. One of:
#'   - "trypsin"      : cleaves after K/R, not before P (standard)
#'   - "trypsin-high" : cleaves after K/R including before P (high specificity)
#'   - "trypsin-low"  : cleaves after K/R, more permissive rules (low specificity)
#' @param missed_cleavages Integer or integer vector of allowed missed cleavages
#'   (e.g. 0 for none, 2 for up to 0:2). Default: 0.
#' @param min_peptide_len  Minimum peptide length to retain. Default: 6.
#' @param max_peptide_len  Maximum peptide length to retain. Default: 30.
#'
#' @return A tibble with Species, Protein_id and observable peptide count
process_tryptic <- function(fasta_paths, output_path,
                            enzym           = "trypsin",
                            missed_cleavages = 0,
                            min_peptide_len  = 6,
                            max_peptide_len  = 30) {

    all_fasta_ids      <- c()
    all_raw_digests    <- tibble()
    all_filtered_digests <- tibble()

    for (fasta_path in fasta_paths) {
        message("\nProcessing:", fasta_path, "\n")

        species   <- tools::file_path_sans_ext(
            str_extract(basename(fasta_path), "^[A-Za-z]+_[A-Za-z]+")
        )
        fasta     <- readAAStringSet(fasta_path)
        fasta_ids <- sub(" .*", "", names(fasta))

        # Warn on duplicate protein IDs across FASTA files (not expected)
        duplicated_ids <- intersect(fasta_ids, all_fasta_ids)
        if (length(duplicated_ids) > 0) {
            warning(
                length(duplicated_ids), " protein ID(s) in ", basename(fasta_path),
                " already present in a previous FASTA file.\n",
                "Duplicated IDs: ", paste(head(duplicated_ids, 5), collapse = ", "),
                if (length(duplicated_ids) > 5) " ..." else ""
            )
        }

        all_fasta_ids <- c(all_fasta_ids, fasta_ids)

        peptides_per_protein <- cleave(fasta, enzym = enzym, missedCleavages = missed_cleavages)
        raw_digests <- map_dfr(names(peptides_per_protein), function(protein_name) {
            tibble(
                Species    = species,
                Protein_id = sub(" .*", "", protein_name),
                peptide    = as.character(peptides_per_protein[[protein_name]])
            )
        })

        filtered_digests <- filter_peptides(raw_digests, min_peptide_len, max_peptide_len)

        all_raw_digests      <- bind_rows(all_raw_digests,      raw_digests)
        all_filtered_digests <- bind_rows(all_filtered_digests, filtered_digests)
    }

    # Peptide filtering comparison table — one row per species
    comparison_table <- all_raw_digests %>%
        count(Species, name = "n_before") %>%
        left_join(
            all_filtered_digests %>% count(Species, name = "n_after"),
            by = "Species"
        ) %>%
        mutate(
            n_removed = n_before - n_after,
            pct_kept  = round(100 * n_after / n_before, 1)
        )

    message("\nPeptide filtering summary (length ", min_peptide_len, "-", max_peptide_len, " aa):\n")
    write.table(
        comparison_table,
        paste0(output_path, "peptide_filter_comparison.tsv"),
        sep = "\t", row.names = FALSE, quote = FALSE
    )

    # Observable peptide counts per protein
    observable_peptides <- all_filtered_digests %>%
        count(Species, Protein_id, name = "n_peptides")

    # Per-species summary statistics
    summary_table <- observable_peptides %>%
        group_by(Species) %>%
        summarise(
            Min    = min(n_peptides),
            Q1     = quantile(n_peptides, 0.25),
            Median = median(n_peptides),
            Mean   = round(mean(n_peptides), 2),
            Q3     = quantile(n_peptides, 0.75),
            Max    = max(n_peptides),
            .groups = "drop"
        )

    message("\nSummary of observable peptide counts after filtering:\n")
    write.table(
        summary_table,
        paste0(output_path, "peptide_summary.tsv"),
        sep = "\t", row.names = FALSE, quote = FALSE
    )

    observable_peptides
}
# ============================================================================
# QUANTIFICATION
# ============================================================================

#' Compute relative intensity per protein per sample
#'
#' Divides total intensity by the number of razor peptides, then normalises
#' within each Species/Replicate group. Adds log2, log10, and parts-per-billion
#' columns.
#'
#' @param tmt_data A unified TMT tibble with Species and Replicate columns
#' @return The input tibble with additional TopN quantification columns
Intensity_quantification <- function(tmt_data) {
  message("Computing intensity quantification for each sample...\n")

  # Compute per-row channel sum (across all channel_ columns)
  channel_cols <- grep("^channel_", colnames(tmt_data), value = TRUE)

  tmt_data %>%
    mutate(
      channel_sum = rowSums(across(all_of(channel_cols)), na.rm = TRUE),
    ) %>%
    pivot_longer(
      cols = all_of(channel_cols),
      names_to = "channel_col",
      values_to = "channel_intensity"
    ) %>%
    # Keep only the channel matching this row's treatment
    filter(channel_col == paste0("channel_", TMT_label)) %>%
    mutate(
      Intensity = .data[[col_total_intensity]] * channel_intensity / channel_sum,
    ) %>%
    filter(if_all(Intensity, ~ . > 0 & !is.na(.))) %>%
    group_by(Species, Protein_id, Treatment) %>%
    mutate(
      Intensity_log2 = log2(Intensity),
      Intensity_meanlog2 = mean(Intensity_log2, na.rm = TRUE),
      Intensity_log2ratio = Intensity_log2 - Intensity_meanlog2,
      Intensity_norm = Intensity / sum(Intensity, na.rm = TRUE),
      Intensity_log = 10 + log10(Intensity_norm),
      Intensity_PpB = Intensity_norm * 1e9
    ) %>%
    ungroup() %>%
    select(-channel_sum, -channel_intensity)
}

#' Compute relative intensity per protein per sample
#'
#' Divides total intensity by the number of razor peptides, then normalises
#' within each Species/Replicate group. Adds log2, log10, and parts-per-billion
#' columns.
#'
#' @param tmt_data A unified TMT tibble with Species and Replicate columns
#' @return The input tibble with additional TopN quantification columns
RI_quantification <- function(tmt_data) {
  message("Computing relative intensity quantification for each sample...\n")

  # Compute per-row channel sum (across all channel_ columns)
  channel_cols <- grep("^channel_", colnames(tmt_data), value = TRUE)

  tmt_data %>%
    mutate(
      channel_sum = rowSums(across(all_of(channel_cols)), na.rm = TRUE),
      RI_total = (.data[[col_total_intensity]] / .data[[col_razor_peptides]])
    ) %>%
    pivot_longer(
      cols = all_of(channel_cols),
      names_to = "channel_col",
      values_to = "channel_intensity"
    ) %>%
    # Keep only the channel matching this row's treatment
    filter(channel_col == paste0("channel_", TMT_label)) %>%
    mutate(
      RI = RI_total * (channel_intensity / channel_sum)
    ) %>%
    filter(if_all(RI, ~ . > 0 & !is.na(.))) %>%
    group_by(Species, Protein_id, Treatment) %>%
    mutate(
      RI_log2 = log2(RI),
      RI_meanlog2 = mean(RI_log2, na.rm = TRUE),
      RI_log2ratio = RI_log2 - RI_meanlog2,
      RI_norm = RI / sum(RI, na.rm = TRUE),
      RI_log = 10 + log10(RI_norm),
      RI_PpB = RI_norm * 1e9
    ) %>%
    ungroup() %>%
    select(-channel_sum, -channel_intensity)
}

#' Compute iBAQ and riBAQ intensities per protein per sample
#'
#' Joins observable peptide counts and divides total intensity by the number
#' of observable peptides (iBAQ). Then normalises within each Species/Replicate
#' group (riBAQ) and adds log-transformed columns.
#'
#' @param tmt_data          A unified TMT tibble with Species and Replicate columns
#' @param observable_peptides A tibble from process_tryptic() with protein_id and n_* columns
#' @param observable_col    Name of the column to use as denominator (default: "n_tryptic")
#' @return The input tibble joined with iBAQ quantification columns
ibaq_quantification <- function(
  tmt_data,
  observable_peptides,
  observable_col = "n_tryptic",
  output_path    = NULL
  ) {
  message("Computing iBAQ for each sample...\n")

  channel_cols <- grep("^channel_", colnames(tmt_data), value = TRUE)

  tmt_data %>%
    left_join(observable_peptides, by = c("Species", "Protein_id")) %>%
    mutate(
      unmatched = is.na(.data[[observable_col]])
    ) %>%
    {
      unmatched_proteins <- distinct(filter(., unmatched), Species, Protein_id)
      n_unmatched <- nrow(unmatched_proteins)

      if (n_unmatched > 0) {
        warning(n_unmatched, " proteins had no match in FASTA and will have iBAQ = NA")

        if (!is.null(output_path)) {
          unmatched_path <- file.path(output_path, "iBAQ_unmatched_proteins.tsv")
          write.table(
            unmatched_proteins %>% arrange(Species, Protein_id),
            unmatched_path,
            sep = "\t", row.names = FALSE, quote = FALSE
          )
          message("Unmatched proteins written to: ", unmatched_path)
        }
      }
      .
    } %>%
    mutate(
      iBAQ_total  = total_intensity / .data[[observable_col]],
      channel_sum = rowSums(across(all_of(channel_cols)), na.rm = TRUE)
    ) %>%
    pivot_longer(
      cols      = all_of(channel_cols),
      names_to  = "channel_col",
      values_to = "channel_intensity"
    ) %>%
    filter(channel_col == paste0("channel_", TMT_label)) %>%
    mutate(iBAQ = iBAQ_total * (channel_intensity / channel_sum)) %>%
    group_by(Species, Protein_id, Treatment) %>%
    mutate(
      iBAQ_log2      = log2(iBAQ),
      iBAQ_meanlog2  = mean(iBAQ_log2, na.rm = TRUE),
      iBAQ_log2ratio = log2(iBAQ) - iBAQ_meanlog2,
      riBAQ          = iBAQ / sum(iBAQ, na.rm = TRUE),
      iBAQ_log       = 10 + log10(riBAQ),
      iBAQ_PpB       = riBAQ * 1e9
    ) %>%
    ungroup()
}