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

library(tidyverse)
library(readxl)
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

# Trypsin digest parameters
missed_cleavages <- 2
min_peptide_len <- 6
max_peptide_len <- 30

# ============================================================================
# TMT DATA PROCESSING
# ============================================================================

#' Read a TMT result TSV file
#'
#' @param path Path to the TSV file
#' @return A tibble with the relevant columns
read_tmt <- function(path) {
  readr::read_tsv(
    path,
    col_names = TRUE,
    show_col_types = FALSE,
    col_select = c(all_of(keep_cols), starts_with("channel_")),
    skip_empty_rows = TRUE
  )
}

#' Filter TMT data for a given species and replicate
#'
#' Applies sequential filters:
#'   1. Keep rows matching the species
#'   2. Remove contaminants
#'   3. Keep proteins with > 1 unique peptide
#'   4. Remove rows with missing Gene
#'   5. Deduplicate by keeping the row with highest unique peptides and intensity per gene
#'
#' @param tmt_df     A tibble from read_tmt()
#' @param curr_species   Species name (as it appears in the Organism column)
#' @param curr_replicate Replicate identifier (e.g. "R1")
#' @return A named list with:
#'   - `filtered_tmt`: filtered tibble
#'   - `filter_log`:   named vector of row counts at each step
filter_tmt <- function(tmt_df, curr_species, curr_replicate) {
  step0 <- tmt_df %>%
    mutate(
      total_intensity = as.numeric(.data[[col_total_intensity]]),
      total_intensity = na_if(total_intensity, 0),
      Protein_id      = .data[[col_protein_id]],
      Organism        = str_replace(.data[[col_organism]], " ", "_")
    )

  step1 <- step0 %>%
    filter(grepl(curr_species, .data[[col_organism]]))

  step2 <- step1 %>%
    filter(!grepl("contam_sp", .data[[col_protein]]))

  step3 <- step2 %>%
    filter(.data[[col_unique_peptides]] > 1)

  step4 <- step3 %>%
    filter(total_intensity > 0)

  step5 <- step4 %>%
    filter(!is.na(.data[[col_gene]]))

  step6 <- step5 %>%
    group_by(.data[[col_gene]]) %>%
    filter(
      .data[[col_unique_peptides]] == max(.data[[col_unique_peptides]]),
      .data[[col_total_intensity]] == max(.data[[col_total_intensity]])
    ) %>%
    ungroup()

  filter_log <- c(
    species = curr_species,
    replicate = curr_replicate,
    initial = nrow(step0),
    remove_non_sp = nrow(step1),
    remove_contam = nrow(step2),
    remove_unique_pep = nrow(step3),
    remove_intensity_0 = nrow(step4),
    remove_na_gene = nrow(step5),
    remove_dedup = nrow(step6)
  )

  list(filtered_tmt = step6, filter_log = filter_log)
}

#' Process all TMT TSV files in a directory
#'
#' Reads, filters, and combines all TSV files found recursively under tmt_path.
#' Expects filenames of the form: <Species>_R<n>.tsv
#' Writes a filter log TSV to output_path.
#'
#' @param tmt_path    Path to the directory containing TSV files
#' @param output_path Path to the output directory (must end with "/")
#' @return A unified tibble with Species and Replicate columns prepended
process_tmt <- function(tmt_path, output_path) {
  tsv_files <- list.files(
    tmt_path,
    pattern    = "\\.tsv$",
    recursive  = TRUE,
    full.names = TRUE
  )

  filter_log <- list()
  filtered_tmt <- list()

  for (tsv in tsv_files) {
    curr_species <- sub(".*/(\\w+)_R\\d+\\.tsv$", "\\1", tsv)
    curr_replicate <- sub(".*/\\w+_(R\\d+)\\.tsv$", "\\1", tsv)

    tmt_data <- read_tmt(tsv)
    filter_result <- filter_tmt(tmt_data, curr_species, curr_replicate)

    filtered_tmt[[curr_species]][[curr_replicate]] <- filter_result$filtered_tmt
    filter_log[[curr_species]][[curr_replicate]] <- filter_result$filter_log
  }

  # Combine all species/replicates into a single tibble
  unified_tmt <- map_dfr(names(filtered_tmt), function(species) {
    map_dfr(names(filtered_tmt[[species]]), function(rep) {
      filtered_tmt[[species]][[rep]] %>%
        mutate(Species = species, Replicate = rep)
    })
  }) %>%
    select(Species, Replicate, everything())

  # Build and write filter log
  filter_log_df <- map_dfr(filter_log, ~ bind_rows(.x)) %>%
    mutate(across(-c("species", "replicate"), as.integer))

  message("TMT filtering summary:\n")
  print(filter_log_df)

  write.table(
    filter_log_df,
    paste0(output_path, "filter_log.tsv"),
    sep       = "\t",
    row.names = FALSE,
    quote     = FALSE
  )

  unified_tmt
}

#' Remove proteins with incomplete data across replicates and treatments
#'
#' Keeps only proteins that have at least one value for every unique
#' Species-Replicate-Treatment combination in the dataset.
#'
#' @param tmt_data A tibble with Species, Replicate, and Treatment columns
#' @return Filtered tibble with complete proteins only
filter_complete_proteins <- function(tmt_data) {
  # Get the total number of unique Species-Replicate-Treatment combinations
  n_conditions <- tmt_data %>%
    distinct(Replicate, Treatment) %>%
    nrow()

  all_proteins <- tmt_data %>%
    distinct(Species, Protein_id)

  n_before <- nrow(all_proteins)

  # Keep only proteins that appear in all conditions
  filtered <- tmt_data %>%
    group_by(Species, Protein_id) %>%
    filter(n_distinct(paste(Replicate, Treatment, sep = "_")) == n_conditions) %>%
    ungroup()

  kept_proteins <- filtered %>%
    distinct(Species, Protein_id)

  n_after <- nrow(kept_proteins)
  n_removed <- n_before - n_after

  # Find removed proteins
  removed_proteins <- all_proteins %>%
    anti_join(kept_proteins, by = c("Species", "Protein_id")) %>%
    # Get additional info for removed proteins from original data
    left_join(
      tmt_data,
      by = c("Species", "Protein_id")
    )

  message("Protein completeness filtering:\n")
  message("  Total unique conditions: ", n_conditions, "\n")
  message("  Proteins before filtering: ", n_before, "\n")
  message("  Proteins after filtering: ", n_after, "\n")
  message("  Proteins removed: ", n_removed, "\n\n")

  list(filtered_data = filtered, removed_proteins = removed_proteins)
}

# ============================================================================
# IN SILICO TRYPSIN DIGEST FROM FASTA
# ============================================================================

#' Filter peptides to those within the configured length range
#'
#' @param peps A character vector of peptide sequences
#' @return A filtered character vector
filter_peptides <- function(peps) {
  lengths <- nchar(peps)
  peps[lengths >= min_peptide_len & lengths <= max_peptide_len]
}

#' Run in silico trypsin digestion on one or more FASTA files
#'
#' Digests proteins using 6 trypsin configurations (base, high, low, each
#' with and without missed cleavages). Filters peptides by length and writes
#' summary tables to output_path. Warns if any protein ID appears in multiple
#' FASTA files.
#'
#' @param fasta_paths A character vector of paths to FASTA files
#' @param output_path Path to the output directory (must end with "/")
#' @return A tibble with protein_id and observable peptide counts per config
process_tryptic <- function(fasta_paths, output_path) {
  # Digest configurations: name -> list(enzym, missedCleavages)
  configs <- list(
    tryptic         = list(enzym = "trypsin", missedCleavages = 0),
    tryptic_high    = list(enzym = "trypsin-high", missedCleavages = 0),
    tryptic_low     = list(enzym = "trypsin-low", missedCleavages = 0),
    tryptic_mc      = list(enzym = "trypsin", missedCleavages = 0:missed_cleavages),
    tryptic_mc_high = list(enzym = "trypsin-high", missedCleavages = 0:missed_cleavages),
    tryptic_mc_low  = list(enzym = "trypsin-low", missedCleavages = 0:missed_cleavages)
  )

  all_fasta_ids <- c()
  all_raw_digests <- list()
  all_filtered_digests <- list()

  for (fasta_path in fasta_paths) {
    message("\nProcessing:", fasta_path, "\n")

    fasta <- readAAStringSet(fasta_path)
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

    # Digest and filter for each config
    raw_digests <- map(configs, ~ as.list(cleave(fasta, enzym = .x$enzym, missedCleavages = .x$missedCleavages)))
    filtered_digests <- map(raw_digests, ~ map(.x, filter_peptides))

    # Accumulate results across FASTA files
    for (cfg in names(configs)) {
      all_raw_digests[[cfg]] <- c(all_raw_digests[[cfg]], raw_digests[[cfg]])
      all_filtered_digests[[cfg]] <- c(all_filtered_digests[[cfg]], filtered_digests[[cfg]])
    }
  }

  # Peptide filtering comparison table (before vs after length filter)
  comparison_table <- map_dfr(names(configs), function(cfg) {
    n_before <- sum(map_int(all_raw_digests[[cfg]], length))
    n_after <- sum(map_int(all_filtered_digests[[cfg]], length))
    tibble(
      config    = cfg,
      n_before  = n_before,
      n_after   = n_after,
      n_removed = n_before - n_after,
      pct_kept  = round(100 * n_after / n_before, 1)
    )
  })

  message("\nPeptide filtering summary (length", min_peptide_len, "-", max_peptide_len, "aa):\n")
  print(comparison_table)
  write.table(
    comparison_table,
    paste0(output_path, "peptide_filter_comparison.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )

  # Observable peptide counts per protein
  observable_peptides <- tibble(Protein_id = all_fasta_ids) %>%
    bind_cols(
      map_dfc(names(all_filtered_digests), ~ tibble(
        !!paste0("n_", .x) := map_int(all_filtered_digests[[.x]], length)
      ))
    )

  # Per-config summary statistics
  summary_table <- map_dfr(names(all_filtered_digests), function(cfg) {
    s <- summary(observable_peptides[[paste0("n_", cfg)]])
    as_tibble(t(as.matrix(s))) %>%
      mutate(config = cfg) %>%
      select("config", everything())
  })

  message("\nSummary of observable peptide counts after filtering:\n")
  print(summary_table)
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
RI_computation <- function(tmt_data) {
  message("Computing relative intensity for each sample...\n")

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
    group_by(Species, Replicate, Treatment) %>%
    mutate(
      RI_log2     = log2(RI),
      RI_log2mean = log2(RI) - mean(log2(RI), na.rm = TRUE),
      RI_norm     = RI / sum(RI, na.rm = TRUE),
      RI_log      = 10 + log10(RI_norm),
      RI_PpB      = RI_norm * 1e9
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
ibaq_computation <- function(tmt_data, observable_peptides, observable_col = "n_tryptic") {
  message("Computing iBAQ for each sample...\n")

  # Compute per-row channel sum (across all channel_ columns)
  channel_cols <- grep("^channel_", colnames(tmt_data), value = TRUE)
  tmt_data %>%
    left_join(observable_peptides, by = "Protein_id") %>%
    mutate(
      unmatched = is.na(.data[[observable_col]])
    ) %>%
    {
      n_unmatched <- sum(.$unmatched)
      if (n_unmatched > 0) {
        warning(n_unmatched, " proteins had no match in FASTA and will have iBAQ = NA")
      }
      .
    } %>%
    mutate(
      iBAQ_total = total_intensity / .data[[observable_col]],
      channel_sum = rowSums(across(all_of(channel_cols)), na.rm = TRUE)
    ) %>%
    pivot_longer(
      cols = all_of(channel_cols),
      names_to = "channel_col",
      values_to = "channel_intensity"
    ) %>%
    # Keep only the channel matching this row's treatment
    filter(channel_col == paste0("channel_", TMT_label)) %>%
    mutate(
      iBAQ = iBAQ_total * (channel_intensity / channel_sum)
    ) %>%
    group_by(Species, Replicate) %>%
    mutate(
      iBAQ_log2      = log2(iBAQ),
      iBAQ_log2_mean = log2(iBAQ) - mean(log2(iBAQ), na.rm = TRUE),
      riBAQ          = iBAQ / sum(iBAQ, na.rm = TRUE),
      iBAQ_log       = 10 + log10(riBAQ),
      iBAQ_PpB       = riBAQ * 1e9
    ) %>%
    ungroup()
}

# ============================================================================
# OUTPUT HELPERS
# ============================================================================

#' Write a tibble to a TSV file
#'
#' @param data        A tibble or data frame
#' @param output_path Output directory path (must end with "/")
#' @param filename    Output file name (e.g. "results.tsv")
write_tsv_table <- function(data, output_path, filename) {
  write.table(
    data,
    paste0(output_path, filename),
    sep       = "\t",
    row.names = FALSE,
    quote     = FALSE
  )
}

# ============================================================================
# MAIN ORCHESTRATION
# ============================================================================

#' Run the full TMT proteomics analysis pipeline
#'
#' @param tmt_path      Directory containing TMT TSV files (searched recursively)
#' @param fasta_path    Directory containing FASTA files (*.fa)
#' @param metadata_path Path to the sample mapping Excel file
#' @param output_path   Output directory (must end with "/")
main_analysis <- function(tmt_path, fasta_path, metadata_path, output_path) {
  # --- Metadata ---
  metadata <- read_excel(metadata_path, col_names = TRUE) %>%
    mutate(Replicate = paste0("R", Replicate))

  # --- TMT filtering ---
  filtered_tmt <- process_tmt(tmt_path, output_path)

  unified_tmt <- filtered_tmt %>%
    left_join(metadata, by = c("Species", "Replicate"), relationship = "many-to-many")

  # --- Filter for complete proteins ---
  completeness_result <- filter_complete_proteins(unified_tmt)
  unified_filtered_tmt <- completeness_result$filtered_data
  removed_proteins <- completeness_result$removed_proteins

  # Write removed proteins table for manual review
  if (nrow(removed_proteins) > 0) {
    write.table(
      removed_proteins %>% arrange(Species, Protein_id),
      paste0(output_path, "removed_proteins_incomplete.tsv"),
      sep = "\t", row.names = FALSE, quote = FALSE
    )
    message("Removed proteins written to: removed_proteins_incomplete.tsv\n")
  }

  # --- TopN quantification ---
  RI_data <- RI_computation(unified_filtered_tmt)

  # --- In silico digest ---
  fasta_files <- list.files(fasta_path, pattern = "\\.fa$", full.names = TRUE)
  protein2tryptic <- process_tryptic(fasta_files, output_path)

  # --- iBAQ quantification ---
  iBAQ_data <- ibaq_computation(unified_tmt, protein2tryptic, "n_tryptic")
  iBAQ_mc_data <- ibaq_computation(unified_tmt, protein2tryptic, "n_tryptic_mc")

  # --- Write outputs ---
  write_tsv_table(
    RI_data %>% select(
      Species, Replicate, Protein_id, Treatment, total_intensity,
      RI, RI_log2, RI_log2mean, RI_norm, RI_log, RI_PpB
    ),
    output_path, "RI_data.tsv"
  )

  ibaq_cols <- c(
    "Species", "Replicate", "Protein_id", "Treatment",
    "iBAQ", "iBAQ_log2", "iBAQ_log2_mean", "riBAQ", "iBAQ_log", "iBAQ_PpB"
  )

  write_tsv_table(iBAQ_data %>% select(all_of(ibaq_cols)), output_path, "iBAQ_data.tsv")
  write_tsv_table(iBAQ_mc_data %>% select(all_of(ibaq_cols)), output_path, "iBAQ_mc_data.tsv")

  message("\nAnalysis complete. Results written to:", output_path, "\n")
}

# ============================================================================
# ENTRY POINT
# ============================================================================

#' Run a test analysis with default paths
test_analysis <- function() {
  main_analysis(
    tmt_path      = "DATA/Proteome/",
    fasta_path    = "DATA/Fasta_files/",
    metadata_path = "DATA/Proteome/sample_mapping.xlsx",
    output_path   = "Analyses/"
  )
}

if (interactive()) {
  test_analysis()
}
