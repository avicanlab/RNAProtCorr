# Function to read transcriptomic and proteomic data --------------------------
#
# Description:
#   Functions for loading and filtering RNA-seq (TPM) and proteomics (TMT)
#   data from Excel and TSV files into tidy long-format tibbles.
#
# Requirements:
#   library(readxl)
#   library(readr)
#   library(dplyr)
#   library(tidyr)
#   library(stringr)
#   library(purrr)

# ── TMT column name constants ──────────────────────────────────────────────────
# These match the exact header names used in the TMT result TSV files.

COL_PROTEIN          <- "Protein"
COL_PROTEIN_ID       <- "Protein ID"
COL_ENTRY_NAME       <- "Entry Name"
COL_GENE             <- "Gene"
COL_ORGANISM         <- "Organism"
COL_UNIQUE_PEPTIDES  <- "Unique Peptides"
COL_RAZOR_PEPTIDES   <- "Razor Peptides"
COL_TOTAL_INTENSITY  <- "Total Intensity"

TMT_KEEP_COLS <- c(
  COL_PROTEIN,
  COL_PROTEIN_ID,
  COL_ENTRY_NAME,
  COL_GENE,
  COL_ORGANISM,
  COL_UNIQUE_PEPTIDES,
  COL_RAZOR_PEPTIDES,
  COL_TOTAL_INTENSITY
)

# ── Transcriptomic data ────────────────────────────────────────────────────────

#' Read TPM Data from Excel
#'
#' @description
#' Loads Transcript Per Million (TPM) expression data from an Excel file,
#' selects relevant columns, and reshapes to long format with separate
#' `Treatment` and `Replicate` columns parsed from the column names.
#'
#' Column names are expected to match the pattern:
#'   `<ABBV>_<Treatment>_<Replicate> (GE) - TPM`
#'
#' @param dataset_path Character. Path to the Excel file containing TPM data.
#' @param sp_abbv      Character. Species abbreviation used in column names
#'   (e.g. `"SALMT"`).
#'
#' @return Tibble with columns: `Species`, `Protein_id`, `Treatment`,
#'   `Replicate`, `TPM`.
read_tpm <- function(dataset_path, sp_abbv) {
  select_columns <- c("Species", "New_locus_tag", "Old_locus_tag")

  # Regex: ABBV_Treatment_Replicate (GE) - TPM
  pattern <- paste0(
    "^", sp_abbv, "_(",
    paste(TREATMENTS_LIST_RNA, collapse = "|"),
    ")_(\\d+) \\(GE\\) - TPM$"
  )

  # Column-selection regex: keep any column that mentions a treatment and TPM
  tpm_col_pattern <- paste0(
    "(", paste(TREATMENTS_LIST_RNA, collapse = "|"), ").*TPM",
    "|TPM.*(", paste(TREATMENTS_LIST_RNA, collapse = "|"), ")"
  )

  readxl::read_excel(dataset_path, col_names = TRUE, na = "N/A") |>
    dplyr::select(
      dplyr::all_of(select_columns),
      dplyr::matches(tpm_col_pattern)
    ) |>
    dplyr::mutate(
      New_locus_tag = dplyr::if_else(
        is.na(New_locus_tag),
        Old_locus_tag,
        New_locus_tag
      )
    ) |>
    dplyr::select(-Old_locus_tag) |>
    dplyr::rename(Protein_id = New_locus_tag) |>
    tidyr::pivot_longer(
      cols      = -c(Species, Protein_id),
      names_to  = "col_name",
      values_to = "TPM"
    ) |>
    dplyr::mutate(
      Treatment = stringr::str_match(col_name, pattern)[, 2],
      Replicate = stringr::str_match(col_name, pattern)[, 3]
    ) |>
    dplyr::select(Species, Protein_id, Treatment, Replicate, TPM)
}

# ── Proteomic data ─────────────────────────────────────────────────────────────

#' Read a TMT Result TSV File
#'
#' @description
#' Reads a single TMT quantification file, retaining only the columns defined
#' in `TMT_KEEP_COLS` plus any `channel_*` intensity columns.
#'
#' @param path Character. Path to the TSV file.
#'
#' @return Tibble with protein metadata and channel intensity columns.
read_tmt <- function(path) {
  readr::read_tsv(
    path,
    col_names       = TRUE,
    show_col_types  = FALSE,
    col_select      = c(dplyr::all_of(TMT_KEEP_COLS), dplyr::starts_with("channel_")),
    skip_empty_rows = TRUE
  )
}

#' Filter a TMT Tibble for One Species and Replicate
#'
#' @description
#' Applies the following sequential filters:
#'   1. Keep only rows matching `curr_species` in the Organism column.
#'   2. Remove contaminant entries (`contam_sp` prefix).
#'   3. Keep proteins with more than one unique peptide.
#'   4. Remove rows where total intensity is zero or missing.
#'   5. Remove rows with a missing Gene identifier.
#'   6. Deduplicate by keeping, per gene, the row with the highest unique
#'      peptide count and highest total intensity.
#'
#' @param tmt_df        Tibble. Output of `read_tmt()`.
#' @param curr_species  Character. Species name as it appears in the Organism
#'   column (e.g. `"Salmonella_enterica"`).
#' @param curr_replicate Character. Replicate identifier (e.g. `"R1"`).
#'
#' @return Named list with two elements:
#'   - `filtered_tmt` : filtered tibble.
#'   - `filter_log`   : named character vector of row counts at each step.
filter_tmt <- function(tmt_df, curr_species, curr_replicate) {
  step0 <- tmt_df |>
    dplyr::mutate(
      total_intensity = as.numeric(.data[[COL_TOTAL_INTENSITY]]),
      total_intensity = dplyr::na_if(total_intensity, 0),
      Protein_id      = .data[[COL_PROTEIN_ID]],
      Organism        = stringr::str_replace(.data[[COL_ORGANISM]], " ", "_")
    )

  step1 <- step0 |>
    dplyr::filter(grepl(curr_species, .data[[COL_ORGANISM]]))

  step2 <- step1 |>
    dplyr::filter(!grepl("contam_sp", .data[[COL_PROTEIN]]))

  step3 <- step2 |>
    dplyr::filter(.data[[COL_UNIQUE_PEPTIDES]] > 1)

  step4 <- step3 |>
    dplyr::filter(total_intensity > 0)

  step5 <- step4 |>
    dplyr::filter(!is.na(.data[[COL_GENE]]))

  step6 <- step5 |>
    dplyr::group_by(.data[[COL_GENE]]) |>
    dplyr::filter(
      .data[[COL_UNIQUE_PEPTIDES]] == max(.data[[COL_UNIQUE_PEPTIDES]]),
      .data[[COL_TOTAL_INTENSITY]] == max(.data[[COL_TOTAL_INTENSITY]])
    ) |>
    dplyr::ungroup()

  filter_log <- c(
    species            = curr_species,
    replicate          = curr_replicate,
    initial            = nrow(step0),
    remove_non_sp      = nrow(step1),
    remove_contam      = nrow(step2),
    remove_unique_pep  = nrow(step3),
    remove_intensity_0 = nrow(step4),
    remove_na_gene     = nrow(step5),
    remove_dedup       = nrow(step6)
  )

  list(filtered_tmt = step6, filter_log = filter_log)
}

#' Remove Proteins With Incomplete Data Across Replicates and Treatments
#'
#' @description
#' Keeps only proteins that have at least one observation for every unique
#' Species × Replicate × Treatment combination in the dataset.
#'
#' @param tmt_data    Tibble with `Species`, `Replicate`, and `Treatment` columns.
#' @param protQ_value Character or NULL. Optional quantification column name; if
#'   provided, rows where this column is `NA` are additionally excluded before
#'   the completeness check. Default NULL.
#'
#' @return Named list:
#'   - `filtered_data`    : tibble of complete proteins only.
#'   - `removed_proteins` : tibble of removed proteins with their original rows.
filter_complete_proteins <- function(tmt_data, protQ_value = NULL) {
  n_conditions <- tmt_data |>
    dplyr::distinct(Replicate, Treatment) |>
    nrow()

  all_proteins <- tmt_data |> dplyr::distinct(Species, Protein_id)
  n_before     <- nrow(all_proteins)

  filtered <- if (!is.null(protQ_value)) {
    dplyr::filter(tmt_data, !is.na(.data[[protQ_value]]))
  } else {
    tmt_data
  }

  filtered <- filtered |>
    dplyr::group_by(Species, Protein_id) |>
    dplyr::filter(
      dplyr::n_distinct(paste(Replicate, Treatment, sep = "_")) == n_conditions
    ) |>
    dplyr::ungroup()

  kept_proteins <- filtered |> dplyr::distinct(Species, Protein_id)
  n_after       <- nrow(kept_proteins)
  n_removed     <- n_before - n_after

  removed_proteins <- all_proteins |>
    dplyr::anti_join(kept_proteins, by = c("Species", "Protein_id")) |>
    dplyr::left_join(tmt_data,       by = c("Species", "Protein_id"))

  message("Protein completeness filtering:")
  message("  Total unique conditions : ", n_conditions)
  message("  Proteins before         : ", n_before)
  message("  Proteins after          : ", n_after)
  message("  Proteins removed        : ", n_removed)

  list(filtered_data = filtered, removed_proteins = removed_proteins)
}

# ── Annotation / gene-list helpers ─────────────────────────────────────────────

#' Read Essential Genes from an Excel File
#'
#' @description
#' Loads an essential-genes Excel file and returns a flat character vector of
#' all gene IDs found across all columns.
#'
#' @param filepath Character. Path to the Excel file.
#'
#' @return Character vector of essential gene IDs.
read_essential_vec <- function(filepath) {
  readxl::read_excel(filepath) |>
    as.list() |>
    unlist() |>
    as.vector()
}

#' Read Stimulon Locus Tags from an Excel File
#'
#' @description
#' Loads a stimulon Excel file (skipping the first header row) and returns
#' a flat character vector of locus tags from the `"Locus Tag"` column.
#'
#' @param filepath Character. Path to the Excel file.
#'
#' @return Character vector of locus tag IDs.
read_stimulon_vec <- function(filepath) {
  readxl::read_excel(filepath, skip = 1) |>
    dplyr::select("Locus Tag") |>
    as.list() |>
    unlist() |>
    as.vector()
}

#' Read DEG Results From a Read-Counts Excel File
#'
#' @description
#' Loads a read-counts Excel file for one species and returns a single
#' long-format tibble with one row per gene × treatment combination.
#' Columns for logFC, FDR, and max group means are parsed from the expected
#' naming convention `"<Treatment> vs. Ctrl - <metric>"`.
#'
#' @param filepath Character. Path to the read-counts Excel file.
#'
#' @return Tibble with columns: `Protein_id`, `Treatment`, `log2FC`, `fdr`,
#'   `max_mean`.
read_deg_data <- function(filepath) {
  col_lfc  <- function(trt) paste0(trt, " vs. Ctrl - Log fold change")
  col_fdr  <- function(trt) paste0(trt, " vs. Ctrl - FDR p-value")
  col_mean <- function(trt) paste0(trt, " vs. Ctrl - Max group means")

  raw <- readxl::read_excel(filepath, col_names = TRUE) |>
    dplyr::mutate(
      Protein_id = dplyr::if_else(
        is.na(New_locus_tag) | New_locus_tag == "N/A",
        Old_locus_tag,
        New_locus_tag
      )
    ) |>
    dplyr::filter(!is.na(Protein_id))

  treatments_present <- STRESS_TREATMENTS_RNA[
    purrr::map_lgl(STRESS_TREATMENTS_RNA, ~ col_lfc(.x) %in% colnames(raw))
  ]

  purrr::map_dfr(treatments_present, function(trt) {
    raw |>
      dplyr::select(
        Protein_id,
        log2FC   = dplyr::all_of(col_lfc(trt)),
        fdr      = dplyr::all_of(col_fdr(trt)),
        max_mean = dplyr::all_of(col_mean(trt))
      ) |>
      dplyr::mutate(
        # Harmonise treatment label: Mig → Hyp
        Treatment = ifelse(trt == "Mig", "Hyp", trt),
        log2FC    = as.numeric(log2FC),
        fdr       = as.numeric(fdr),
        max_mean  = as.numeric(max_mean)
      )
  })
}
