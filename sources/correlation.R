# ============================================================================
# mRNA-Protein Correlation Analysis
# ============================================================================
#
# Description:
#   This script processes Transcript read counts data.
#
# Requirements:
#   install.packages("tidyverse")
#   install.packages("readxl")
# ============================================================================

library(tidyverse)
library(readxl)

# ============================================================================
# CONFIGURATION
# ============================================================================

species_abbv_map <- c(
    "Salmonella_enterica"           = "SALMT",
    "Staphylococcus_aureus"         = "MSSA",
    "Yersinia_pseudotuberculosis"   = "YPSTB"
)
treatment_list <- c(
    "Ctrl", "As", "Bs", "Hyp", "Li", "Nd", "Ns", "Oss", "Oxs", "Sp", "Tm"
)

# ============================================================================
# Process transcriptomic data
# ============================================================================

read_tpm <- function(dataset_path, sp_abbv) {
    select_columns <- c("Species", "Protein_id")
    # Build pattern: ABBV_Treatment_Replicate (GE) - TPM
    pattern <- paste0(
        "^", sp_abbv, "_(",
        paste(treatment_list, collapse = "|"),
        ")_(\\d+) \\(GE\\) - TPM$"
    )
    data <- read_excel(dataset_path, col_names = TRUE) %>%
        select(
            all_of(select_columns),
            matches(paste0("(", paste(treatment_list, collapse = "|"), ").*TPM|TPM.*(", paste(treatment_list, collapse = "|"), ")"))
        ) %>%
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

process_tpm <- function(tpm_files) {
    # Process TPM data
    # Return: processed TPM data frame
    map_dfr(tpm_files, function(path) {
        # Extract species name from filename
        curr_species <- str_extract(basename(path), "^[A-Za-z]+_[A-Za-z]+")

        # Get abbreviation from map
        abbv <- species_abbv_map[curr_species]
        if (is.na(abbv)) {
            warning("No abbreviation found for species '", curr_species, "', skipping file: ", basename(path))
            return(tibble())
        }

        cat("Processing:", basename(path), "| Species:", curr_species, "\n")

        read_tpm(path, abbv) %>%
            mutate(Species = curr_species)
    })
}

read_ibaq <- function() {
    # Load iBAQ data from file
    # Return: data frame with iBAQ values
}

# 2. DATA TRANSFORMATION FUNCTIONS
log2_transform <- function(data) {
    # Apply log2 transformation to data
    # Handle mean calculation if needed
    # Return: transformed data frame
}

# 3. CORRELATION ANALYSIS FUNCTIONS
calculate_correlation <- function(tpm_data, ibaq_data) {
    # Calculate correlation between TPM and iBAQ
    # Return: correlation results (coefficients, p-values, etc.)
}

# 4. VISUALIZATION FUNCTIONS
plot_correlation <- function(correlation_results) {
    # Create plots for 11 conditions
    # Return: visualization(s)
}

# 5. MAIN ORCHESTRATION FUNCTION
main_analysis <- function(TPM_PATH, TOPN_PATH, IBAQ_PATH, IBAQ_MC_PATH) {
    # Load data
    tpm_files <- list.files(TPM_PATH, pattern = "\\.xlsx$", full.names = TRUE)

    tpm_data <- process_tpm(tpm_files)

    topN_data <- read_tsv(
        TOPN_PATH,
        col_names = TRUE,
        show_col_types = FALSE,
        skip_empty_rows = TRUE
    )

    iBAQ_data <- read_tsv(
        IBAQ_PATH,
        col_names = TRUE,
        show_col_types = FALSE,
        skip_empty_rows = TRUE
    )

    iBAQ_mc_data(
        IBAQ_MC_PATH,
        col_names = TRUE,
        show_col_types = FALSE,
        skip_empty_rows = TRUE
    )

    # # Transform data
    # tpm_log2 <- log2_transform(tpm)
    # ibaq_log2 <- log2_transform(ibaq)

    # # Calculate and visualize
    # results <- calculate_correlation(tpm_log2, ibaq_log2)
    # plot_correlation(results)

    # return(results)
    print("Done")
}

# 6. TEST FUNCTION
test_analysis <- function() {
    # Run main analysis for testing
    main_analysis(
        "DATA/Read_counts/",
        "Analyses/TopN_data.tsv",
        "Analyses/iBAQ_data.tsv",
        "Analyses/iBAQ_mc_data.tsv"
    )
}

# Execute test when script runs independently
if (interactive()) {
    test_analysis()
}
