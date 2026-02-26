# ============================================================================
# Protein quantification from TMT MS data
# ============================================================================

# Requirements: install packages if needed
# install.packages("tidyverse")
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("Biostrings")
# BiocManager::install("cleaver")

library(tidyverse)
library(readxl)
library(Biostrings)
library(cleaver)
library(dplyr)
library(purrr)

# ============================================================================
# CONFIGURATION
# ============================================================================

# Column names in TMT results
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
# PROCESS TMT DATA
# ============================================================================

read_tmt <- function(path) {
  tmt_df <- readr::read_tsv(
    path,
    col_names = TRUE,
    show_col_types = FALSE,
    col_select = all_of(keep_cols),
    skip_empty_rows = TRUE
  )

  tmt_df
}

filter_tmt <- function(tmt_df, curr_species, curr_replicate) {
  # Replace 0s and empty strings with NA in intensity column
  step0 <- tmt_df %>%
    mutate(
      total_intensity = as.numeric(.data[[col_total_intensity]]),
      total_intensity = na_if(.data$total_intensity, 0),
      protein_id = .data[[col_protein_id]],
      Organism = .data[[col_organism]] %>% str_replace(" ", "_")
    )

  step1 <- step0 %>%
    filter(grepl(curr_species, .data[[col_organism]]))

  step2 <- step1 %>%
    filter(!grepl("contam_sp", col_protein))

  step3 <- step2 %>%
    filter(.data[[col_unique_peptides]] > 1)

  # From now filtering seems to have no effect, but we keep it for consistency with original script
  step4 <- step3 %>%
    filter(!is.na(.data[[col_gene]]))

  step5 <- step4 %>%
    group_by(.data[[col_gene]]) %>%
    filter(
      .data[[col_unique_peptides]] == max(.data[[col_unique_peptides]]),
      .data[[col_total_intensity]] == max(.data[[col_total_intensity]])
    ) %>%
    ungroup()

  # Store counts for this species
  filter_log <- c(
    species           = curr_species,
    replicate         = curr_replicate,
    initial           = nrow(step0),
    remove_non_sp     = nrow(step1),
    remove_contam     = nrow(step2),
    remove_unique_pep = nrow(step3),
    remove_na_gene    = nrow(step4),
    remove_dedup      = nrow(step5)
  )

  filtered_tmt <- step5

  list(filtered_tmt = filtered_tmt, filter_log = filter_log)
}

process_tmt <- function(tmt_path, output_path) {
  tsv_files <- list.files(
    tmt_path,
    pattern = "\\.tsv$",
    recursive = TRUE,
    full.names = TRUE
  )

  filter_log <- list()
  filtered_tmt <- list()

  for (tsv in tsv_files) {
    # cat("Processing file:", tsv, "\n")
    curr_species <- sub(".*/(\\w+)_R\\d+\\.tsv$", "\\1", tsv)
    curr_replicate <- sub(".*/\\w+_(R\\d+)\\.tsv$", "\\1", tsv)

    tmt_data <- read_tmt(tsv)

    filter_result <- filter_tmt(tmt_data, curr_species, curr_replicate)

    filtered_tmt[[curr_species]][[curr_replicate]] <- filter_result$filtered_tmt
    filter_log[[curr_species]][[curr_replicate]] <- filter_result$filter_log
  }

  unified_tmt <- map_dfr(names(filtered_tmt), function(species) {
    map_dfr(names(filtered_tmt[[species]]), function(rep) {
      filtered_tmt[[species]][[rep]] %>%
        mutate(Species = species, Replicate = rep)
    })
  }) %>%
    select(Species, Replicate, everything())

  filter_log_df <- map_dfr(filter_log, ~ bind_rows(.x)) %>%
    mutate(across(-c("species", "replicate"), as.integer))

  cat("TMT filtering summary:\n")
  print(filter_log_df)

  write.table(
    filter_log_df,
    paste0(output_path, "filter_log.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  unified_tmt
}

# ============================================================================
# IN SILICO TRYPSIN DIGEST FROM FASTA
# ============================================================================

# Filter peptides within length range
filter_peptides <- function(peps) {
  lengths <- nchar(peps)
  peps[lengths >= min_peptide_len & lengths <= max_peptide_len]
}

process_tryptic <- function(fasta_path, output_path) {
  fasta <- readAAStringSet(fasta_path)
  fasta_ids <- sub(" .*", "", names(fasta))

  # Define enzymes and missed cleavage configs
  configs <- list(
    tryptic         = list(enzym = "trypsin", missedCleavages = 0),
    tryptic_high    = list(enzym = "trypsin-high", missedCleavages = 0),
    tryptic_low     = list(enzym = "trypsin-low", missedCleavages = 0),
    tryptic_mc      = list(enzym = "trypsin", missedCleavages = 0:missed_cleavages),
    tryptic_mc_high = list(enzym = "trypsin-high", missedCleavages = 0:missed_cleavages),
    tryptic_mc_low  = list(enzym = "trypsin-low", missedCleavages = 0:missed_cleavages)
  )

  # Run all digestions
  raw_digests <- map(
    configs,
    ~ as.list(cleave(
      fasta,
      enzym = .x$enzym,
      missedCleavages = .x$missedCleavages
    ))
  )
  filtered_digests <- map(raw_digests, ~ map(.x, filter_peptides))

  # Build comparison table: before vs after filtering
  comparison_table <- map_dfr(names(configs), function(cfg) {
    n_before <- sum(map_int(raw_digests[[cfg]], length))
    n_after <- sum(map_int(filtered_digests[[cfg]], length))
    tibble(
      config = cfg,
      n_before = n_before,
      n_after = n_after,
      n_removed = n_before - n_after,
      pct_kept = round(100 * n_after / n_before, 1)
    )
  })

  cat("\nPeptide filtering summary (length", min_peptide_len, "-", max_peptide_len, "aa):\n")
  print(comparison_table)

  write.table(
    comparison_table,
    paste0(output_path, "peptide_filter_comparison.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  # Build final observable peptides table
  observable_peptides <- tibble(protein_id = fasta_ids) %>%
    bind_cols(
      map_dfc(names(filtered_digests), ~ tibble(
        !!paste0("n_", .x) := map_int(filtered_digests[[.x]], length)
      ))
    )

  cat("\nSummary of tryptic peptide counts after filtering:\n")

  summary_table <- map_dfr(names(filtered_digests), function(cfg) {
    s <- summary(observable_peptides[[paste0("n_", cfg)]])
    as_tibble(t(as.matrix(s))) %>%
      mutate(config = cfg) %>%
      select("config", everything())
  })

  print(summary_table)

  write.table(
    summary_table,
    paste0(output_path, "peptide_summary.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  observable_peptides
}

# ============================================================================
# Relative TopN computation (total intensity / number of razor peptides)
# ============================================================================
relative_TopN_computation <- function(tmt_data) {
  # Compute relative TopN for each sample by dividing the total intensity of each protein by the number of razor peptides
  print("Computing relative TopN for each sample...")
  tmt_data <- tmt_data %>%
    mutate(relative_TopN = .data[[col_total_intensity]] / .data[[col_razor_peptides]]) %>%
    group_by(Species, Replicate) %>%
    mutate(
      relative_TopN_log2      = log2(relative_TopN),
      relative_TopN_log2mean  = log2(relative_TopN) - mean(log2(relative_TopN), na.rm = TRUE),
      relative_TopN_norm      = relative_TopN / sum(relative_TopN, na.rm = TRUE),
      relative_TopN_log       = 10 + log10(relative_TopN_norm),
      relative_TopN_PpB       = relative_TopN_norm * 1e9
    ) %>%
    ungroup()
  tmt_data
}

# ============================================================================
# iBAQ / rIBAQ computation
# ============================================================================
ibaq_computation <- function(tmt_data, observable_peptides, observable_col = "n_tryptic") {
  iBAQ_data <- tmt_data %>%
    # Join observable peptide counts
    left_join(observable_peptides, by = "protein_id") %>%
    # Warn if any protein_id had no match in FASTA
    mutate(
      unmatched = is.na(.data[[observable_col]])
    ) %>%
    {
      n_unmatched <- sum(.$unmatched)
      if (n_unmatched > 0) warning(n_unmatched, " proteins had no match in FASTA and will have iBAQ = NA")
      .
    } %>%
    # Compute iBAQ: total intensity / number of observable peptides
    mutate(
      iBAQ = total_intensity / .data[[observable_col]]
    ) %>%
    # Compute riBAQ and normalized values per sample
    group_by(Species, Replicate) %>%
    mutate(
      iBAQ_log2       = log2(iBAQ),
      iBAQ_log2_mean  = log2(iBAQ) - mean(log2(iBAQ), na.rm = TRUE),
      riBAQ           = iBAQ / sum(iBAQ, na.rm = TRUE),
      iBAQ_log        = 10 + log10(riBAQ),
      iBAQ_PpB        = riBAQ * 1e9
    ) %>%
    ungroup()

  iBAQ_data
}

# =============================================================================
# 6. EXECUTION
# =============================================================================

# (Sections for iBAQ, transformations and exports are left commented
#  as in the original script.)

# 5. MAIN ORCHESTRATION FUNCTION
main_analysis <- function(tmt_path, fasta_path, metadata_path, output_path) {
  # Read metadata
  metadata <- read_excel(metadata_path, col_names = TRUE) %>%
    mutate(Replicate = paste0("R", Replicate))

  # Process the TMT data
  unified_filtered_tmt <- process_tmt(tmt_path, output_path)

  unified_metadata_filtered_tmt <- unified_filtered_tmt %>%
    left_join(metadata, by = c("Species", "Replicate"), relationship = "many-to-many")

  # Get the relative TopN matrix for each sample
  TopN_data <- relative_TopN_computation(unified_metadata_filtered_tmt)

  # Generate tryptic peptide counts from FASTA
  protein2tryptic <- process_tryptic(fasta_path, output_path)
  iBAQ_data <- ibaq_computation(unified_metadata_filtered_tmt, protein2tryptic)
  iBAQ_mc_data <- ibaq_computation(unified_metadata_filtered_tmt, protein2tryptic, "n_tryptic_mc")

  write.table(
    TopN_data %>%
      select(
        Species,
        Replicate,
        protein_id,
        Treatment,
        total_intensity,
        relative_TopN,
        relative_TopN_log2,
        relative_TopN_log2mean,
        relative_TopN_norm,
        relative_TopN_log,
        relative_TopN_PpB,
      ),
    paste0(output_path, "TOPN_data.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  write.table(
    iBAQ_data %>%
      select(
        Species,
        Replicate,
        protein_id,
        Treatment,
        iBAQ,
        iBAQ_log2,
        iBAQ_log2_mean,
        riBAQ,
        iBAQ_log,
        iBAQ_PpB
      ),
    paste0(output_path, "iBAQ_data.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  write.table(
    iBAQ_mc_data %>%
      select(
        Species,
        Replicate,
        protein_id,
        Treatment,
        iBAQ,
        iBAQ_log2,
        iBAQ_log2_mean,
        riBAQ,
        iBAQ_log,
        iBAQ_PpB
      ),
    paste0(output_path, "iBAQ_mc_data.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}

# 6. TEST FUNCTION
test_analysis <- function() {
  # Run main analysis for testing
  tmt_path <- "DATA/Proteome/"
  fasta_path <- "DATA/Fasta_files/Salmonella_enterica_subsp._enterica_serovar_Typhimurium_str._SL1344.fa"
  metadata_path <- "DATA/Proteome/sample_mapping.xlsx"
  output_path <- "Analyses/"
  main_analysis(
    tmt_path = tmt_path,
    fasta_path = fasta_path,
    metadata_path = metadata_path,
    output_path = output_path
  )
}

# Execute test when script runs interactively
if (interactive()) {
  test_analysis()
}
