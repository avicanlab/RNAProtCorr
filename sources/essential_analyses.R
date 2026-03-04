# ============================================================================
# Essential genes analyses
# ============================================================================
#
# Description:
# ============================================================================

library(tidyverse)
library(readxl)

# ============================================================================
# CONFIGURATION
# ============================================================================

#' Species Abbreviation Mapping
#' Maps full species names to abbreviations for data parsing
species_essential2tpm_map <- c(
    "Senterica"             = "Salmonella_enterica",
    "Saureus"               = "Staphylococcus_aureus",
    "Ypseudotuberculosis"   = "Yersinia_pseudotuberculosis"
)

#' Species Abbreviation Mapping
#' Maps full species names to abbreviations for data parsing
species_abbv_map <- c(
    "Salmonella_enterica"         = "SALMT",
    "Staphylococcus_aureus"       = "MSSA",
    "Yersinia_pseudotuberculosis" = "YPSTB"
)

#' Treatment List
#' Control and stress treatment conditions
treatment_list <- c(
    "Ctrl", "As", "Bs", "Mig", "Li", "Nd", "Ns", "Oss", "Oxs", "Sp", "Tm"
)

read_tpm <- function(dataset_path, sp_abbv) {
    select_columns <- c("Species", "New_locus_tag", "Old_locus_tag")

    # Build regex pattern: ABBV_Treatment_Replicate (GE) - TPM
    pattern <- paste0(
        "^", sp_abbv, "_(",
        paste(treatment_list, collapse = "|"),
        ")_(\\d+) \\(GE\\) - TPM$"
    )

    data <- read_excel(dataset_path, col_names = TRUE) %>%
        select(
            all_of(select_columns),
            matches(
                paste0(
                    "(",
                    paste(treatment_list, collapse = "|"),
                    ").*TPM|TPM.*(", paste(treatment_list, collapse = "|"), ")"
                )
            )
        ) %>%
        mutate(
            New_locus_tag = if_else(
                is.na(New_locus_tag) | New_locus_tag == "N/A",
                Old_locus_tag,
                New_locus_tag
            )
        )

    data %>%
        select(-Old_locus_tag) %>%
        rename(Protein_id = New_locus_tag) %>%
        pivot_longer(
            cols      = -c(Species, Protein_id),
            names_to  = "col_name",
            values_to = "TPM"
        ) %>%
        mutate(
            Treatment = str_match(col_name, pattern)[, 2],
            Replicate = str_match(col_name, pattern)[, 3]
        ) %>%
        select(Species, Protein_id, Treatment, Replicate, TPM)
}

# ============================================================================
# MAIN ORCHESTRATION
# ============================================================================

#' Process Single Species Analysis
#'
#' @description
#' Complete pipeline for one species: loads and transforms TPM data,
#' filters protein quantification tables, computes correlations, and
#' generates visualizations for multiple quantification measures.
#'
#' @param essential_file Character. Path to species essential genes Excel file
#' @param tpm_file Character. Path to species TPM Excel file
#' @param topN_data Tibble. TopN quantification data (all species)
#' @param iBAQ_data Tibble. iBAQ quantification data (all species)
#' @param iBAQ_mc_data Tibble. iBAQ MC quantification data (all species)
#' @param output_path Character. Directory for saving plots
#'
#' @return NULL (invisibly). Generates plots and console messages
process_species <- function(
  essential_vec, tpm_df, RI_data, iBAQ_data, iBAQ_mc_data, output_path
) {
    print(essential_vec)
    print(tpm_df)
    print(RI_data)
}

#' Run the full essential genes analyses
#'
#'
#' @param output_path  Root output directory.
#' @return NULL invisibly.
main_analysis <- function(
  essential_path, tpm_path, RI_path, iBAQ_path, iBAQ_mc_path, output_path
) {
    # Load protein quantification reference tables (contain all species)
    RI_data <- read_tsv(RI_path,
        col_names = TRUE,
        show_col_types = FALSE,
        skip_empty_rows = TRUE
    )
    iBAQ_data <- read_tsv(iBAQ_path,
        col_names = TRUE,
        show_col_types = FALSE,
        skip_empty_rows = TRUE
    )
    iBAQ_mc_data <- read_tsv(iBAQ_mc_path,
        col_names = TRUE,
        show_col_types = FALSE,
        skip_empty_rows = TRUE
    )

    extract_species_tpm <- function(filepath) {
        tools::file_path_sans_ext(
            str_extract(basename(filepath), "^[A-Za-z]+_[A-Za-z]+")
        )
    }
    extract_species_essential <- function(filepath) {
        species_name <- tools::file_path_sans_ext(
            str_extract(basename(filepath), "Essential_genes_([A-Za-z]+)", group = 1)
        )
        species_essential2tpm_map[species_name]
    }
    # Build a named list of file pairs, keyed by species
    tpm_files <- list.files(tpm_path, pattern = "\\.xlsx$", full.names = TRUE)
    essential_files <- list.files(essential_path, pattern = "\\.xlsx$", full.names = TRUE)
    species_files <- map(set_names(map_chr(tpm_files, extract_species_tpm)), function(species) {
        tpm_idx <- match(species, map_chr(tpm_files, extract_species_tpm))
        essential_idx <- match(species, map_chr(essential_files, extract_species_essential))
        list(
            tpm       = tpm_files[[tpm_idx]],
            essential = if (!is.na(essential_idx)) essential_files[[essential_idx]] else NULL
        )
    })


    # Process each species with error handling
    results <- lapply(names(species_files), function(species) {
        message("Processing species: ", species)
        files <- species_files[[species]]
        abbv <- species_abbv_map[species]
        if (is.na(abbv)) {
            warning(
                "No abbreviation found for species '", species,
                "', skipping file: ", basename(tpm_file)
            )
            return(invisible(NULL))
        }
        tpm_df <- read_tpm(files$tpm, abbv) %>%
            mutate(Species = species) %>%
            # Update treatment name for consistency ("Mig" to "Hyp")
            mutate(Treatment = ifelse(Treatment == "Mig", "Hyp", Treatment)) %>%
            # For Y.pseudotuberculosis: prefix locus tags with 'PI'
            mutate(Protein_id = sub("^pi", "PI", Protein_id, ignore.case = FALSE)) %>%
            filter(!is.na(Protein_id))

        essential_vec <- as.vector(
            read_excel(files$essential, col_names = FALSE)
        ) %>%
            unlist()

        process_species(
            essential_vec, tpm_df, RI_data, iBAQ_data, iBAQ_mc_data, output_path
        )
    })

    message("All species processed.")
}


# ============================================================================
# ENTRY POINT
# ============================================================================

#' Run analysis with default development paths
#' @return NULL invisibly.
test_analysis <- function() {
    main_analysis(
        essential_path  = "DATA/Essential genes",
        tpm_path        = "DATA/Read_counts/",
        RI_path         = "Analyses/RI_data.tsv",
        iBAQ_path       = "Analyses/iBAQ_data.tsv",
        iBAQ_mc_path    = "Analyses/iBAQ_mc_data.tsv",
        output_path     = "Analyses"
    )
}

if (interactive()) {
    test_analysis()
}
