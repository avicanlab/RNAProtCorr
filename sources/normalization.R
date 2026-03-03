# ============================================================================
# Protein Normalization
# ============================================================================
#
# Description
# Comprehensive pipeline for analyzing correlations between transcript
# abundance (TPM) and protein quantification (iBAQ, TopN) across multiple
# bacterial species and experimental conditions. Generates summary statistics
# and visualization panels.

# ============================================================================
# LOAD DEPENDENCIES
# ============================================================================

library(readxl)
library(tidyverse)
library(NormalyzerDE)

# ============================================================================
# CONFIGURATION
# ============================================================================

norm_methods <- c(
    "Log2", "VSN", "CycLoess", "Quantile", "Median", "Mean", "GI", "RLR"
)
# ============================================================================
# NORMALIZATION
# ============================================================================
create_design_for_matrix <- function(data_path, tmp_dir = tempdir()) {
    mat <- read.table(data_path, sep = "\t", header = TRUE, check.names = FALSE)
    sample_cols <- colnames(mat)[-1] # exclude Protein_id

    design_df <- data.frame(
        sample = sample_cols,
        group = sub("^R[0-9]+_", "", sample_cols) # strip replicate -> treatment only
    )

    temp_file <- tempfile(pattern = "design_", tmpdir = tmp_dir, fileext = ".tsv")
    write.table(design_df, temp_file, sep = "\t", row.names = FALSE, quote = FALSE)
    return(temp_file)
}

transform_data <- function(data, protQ, tmp_dir = tempdir()) {
    # Placeholder for data transformation logic
    # This function would apply necessary transformations (e.g., log2)
    # to the raw iBAQ data before normalization.
    message("Transforming", protQ, "data", "\n")
    data_wide <- data %>%
        mutate(sample = paste0(Replicate, "_", Treatment)) %>%
        pivot_wider(
            id_cols = Protein_id,
            names_from = sample,
            values_from = protQ
        )
    data_wide <- data_wide %>%
        select(Protein_id, sort(colnames(data_wide)[-1]))
    temp_file <- tempfile(
        pattern = paste0("transformed_", protQ, "_"),
        tmpdir = tmp_dir,
        fileext = ".tsv"
    )
    write.table(
        data_wide, temp_file,
        sep = "\t", row.names = FALSE, quote = FALSE
    )
    return(temp_file) # Return transformed data (placeholder)
}

normalize_data <- function(data, design_path, protQ, output_path) {
    # Placeholder for normalization logic using NormalyzerDE
    # This function would take raw iBAQ data and apply the specified
    # normalization method, returning normalized data ready for analysis.
    message("Normalizing", protQ, "data", "\n")

    job_name <- paste0("Normalization_", protQ)
    exp_obj <- setupRawDataObject(data, design_path)
    norm_obj <- getVerifiedNormalyzerObject(job_name, exp_obj)
    norm_results <- normMethods(norm_obj)
    norm_results_with_eval <- analyzeNormalizations(norm_results)
    writeNormalizedDatasets(
        norm_results_with_eval,
        output_path,
        # includeCvCol = TRUE,
        # includeAnovaP = TRUE,
        normSuffix = paste0("_", protQ, ".tsv")
    )
    generatePlots(norm_results_with_eval, output_path)
    generatePlots(norm_results_with_eval, output_path, writeAsPngs = TRUE)
}

# ============================================================================
# VISUALIZATION
# ============================================================================

plot_pca <- function(norm_matrix, design_table, method_name) {
    # Protein_id is a column, set it as rownames and keep only numeric columns
    mat <- norm_matrix %>%
        select(-one_of(c("CV", "anovaP"))) %>%
        column_to_rownames("Protein_id") %>%
        select(where(is.numeric)) %>%
        as.matrix()
    # Remove rows with NA
    mat <- mat[complete.cases(mat), ]

    # PCA (samples as rows, so transpose)
    pca <- prcomp(t(mat), scale. = TRUE, center = TRUE)

    design_df <- read.table(design_table, header = TRUE, sep = "\t") %>%
        mutate(
            Treatment = group,
            Replicate = sub("_[a-zA-Z]+", "", sample) # strip treatment -> replicate only
        )

    # Build plot df
    pca_df <- as.data.frame(pca$x[, 1:2]) %>%
        rownames_to_column("sample") %>%
        left_join(design_df, by = "sample")


    # Variance explained
    var_exp <- round(summary(pca)$importance[2, 1:2] * 100, 1)

    ggplot(pca_df, aes(x = PC1, y = PC2, color = Treatment, shape = factor(Replicate), label = sample)) +
        geom_point(size = 3) +
        labs(
            title = paste("PCA —", method_name),
            x = paste0("PC1 (", var_exp[1], "%)"),
            y = paste0("PC2 (", var_exp[2], "%)"),
            shape = "Replicate"
        ) +
        theme_bw()
    return(list(df = pca_df, plot = last_plot()))
}

process_plot_pca <- function(norm_results, design_df, output_path, method_name) {
    norm_matrix <- read.table(
        norm_results,
        sep = "\t", header = TRUE, check.names = FALSE
    )
    pca_res <- plot_pca(norm_matrix, design_df, method_name)
    output_tsv <- paste0(output_path, "/", method_name, "_PCA_coordinates.tsv")
    write.table(
        pca_res$df,
        file = output_tsv,
        sep = "\t", row.names = FALSE, quote = FALSE
    )
    message("PCA coordinates saved to: ", output_tsv, "\n")

    output_png <- paste0(output_path, "/", method_name, ".png")
    ggsave(
        filename = output_png,
        plot = pca_res$plot,
        width = 6,
        height = 5
    )
    message("PCA plot saved to: ", output_png, "\n")
}

# ============================================================================
# MAIN ORCHESTRATION
# ============================================================================
#' Main Analysis Pipeline
#'
#' @description
#' Orchestrates complete mRNA-protein correlation analysis workflow.
#' Loads all data files, iterates through species, and generates
#' correlation statistics and visualization panels.
#'
#' @param TPM_PATH Character. Directory containing TPM Excel files
#' @param TOPN_PATH Character. Path to TopN TSV file
#' @param IBAQ_PATH Character. Path to iBAQ TSV file
#' @param IBAQ_MC_PATH Character. Path to iBAQ MC TSV file
#' @param output_path Character. Directory for saving outputs
#'
#' @return NULL (invisibly). Generates console output and plot files
main_analysis <- function(RI_PATH, IBAQ_PATH, IBAQ_MC_PATH, design_path, output_path) {
    # Load data files
    RI_data <- read_tsv(RI_PATH,
        col_names = TRUE,
        show_col_types = FALSE,
        skip_empty_rows = TRUE
    )
    iBAQ_data <- read_tsv(IBAQ_PATH,
        col_names = TRUE,
        show_col_types = FALSE,
        skip_empty_rows = TRUE
    )
    iBAQ_mc_data <- read_tsv(IBAQ_MC_PATH,
        col_names = TRUE,
        show_col_types = FALSE,
        skip_empty_rows = TRUE
    )

    # Extract unique species from RI data
    species_list <- unique(RI_data$Species)

    temp_dir <- tempdir(check = FALSE)
    output_dir <- paste(output_path, "NormTest", sep = "/")
    dir.create(output_dir, showWarnings = FALSE)

    # Iterate through each species and perform analysis
    for (species in species_list) {
        message("Processing species:", species, "\n")
        output_dir_sp <- paste0(output_dir, "/", species)
        dir.create(output_dir_sp, showWarnings = FALSE)
        # Subset data for current species
        RI_subset <- RI_data %>%
            filter(Species == species) %>%
            transform_data(protQ = "RI", tmp_dir = temp_dir)
        iBAQ_subset <- iBAQ_data %>%
            filter(Species == species) %>%
            transform_data(protQ = "iBAQ", tmp_dir = temp_dir)
        iBAQ_mc_subset <- iBAQ_mc_data %>%
            filter(Species == species) %>%
            transform_data(protQ = "iBAQ", tmp_dir = temp_dir)
        # Create design matched to each matrix
        design_RI <- create_design_for_matrix(RI_subset, tmp_dir = temp_dir)
        design_iBAQ <- create_design_for_matrix(iBAQ_subset, tmp_dir = temp_dir)
        design_iBAQ_MC <- create_design_for_matrix(iBAQ_mc_subset, tmp_dir = temp_dir)

        # Perform normalization on each subset
        output_dir_sp_RI <- paste0(output_dir_sp, "/RI")
        dir.create(output_dir_sp_RI, showWarnings = FALSE)
        normalize_data(RI_subset, design_RI, "RI", output_dir_sp_RI)

        output_dir_sp_iBAQ <- paste0(output_dir_sp, "/iBAQ")
        dir.create(output_dir_sp_iBAQ, showWarnings = FALSE)
        normalize_data(iBAQ_subset, design_iBAQ, "iBAQ", output_dir_sp_iBAQ)

        output_dir_sp_iBAQ_MC <- paste0(output_dir_sp, "/iBAQ_MC")
        dir.create(output_dir_sp_iBAQ_MC, showWarnings = FALSE)
        normalize_data(iBAQ_mc_subset, design_iBAQ_MC, "iBAQ_MC", output_dir_sp_iBAQ_MC)

        output_dir_sp_RI_PCA <- paste0(output_dir_sp_RI, "/PCA")
        dir.create(output_dir_sp_RI_PCA, showWarnings = FALSE)
        output_dir_sp_iBAQ_PCA <- paste0(output_dir_sp_iBAQ, "/PCA")
        dir.create(output_dir_sp_iBAQ_PCA, showWarnings = FALSE)
        output_dir_sp_iBAQ_MC_PCA <- paste0(output_dir_sp_iBAQ_MC, "/PCA")
        dir.create(output_dir_sp_iBAQ_MC_PCA, showWarnings = FALSE)

        for (norm_method in norm_methods) {
            message("Processing:", norm_method, "\n")
            norm_file_RI <- paste0(output_dir_sp_RI, "/", norm_method, "_RI.tsv")
            norm_file_iBAQ <- paste0(output_dir_sp_iBAQ, "/", norm_method, "_iBAQ.tsv")
            norm_file_iBAQ_MC <- paste0(output_dir_sp_iBAQ_MC, "/", norm_method, "_iBAQ_MC.tsv")

            if (file.exists(norm_file_RI)) {
                process_plot_pca(
                    norm_file_RI,
                    design_RI,
                    output_dir_sp_RI_PCA,
                    paste0("RI_", norm_method)
                )
            }
            if (file.exists(norm_file_iBAQ)) {
                process_plot_pca(
                    norm_file_iBAQ,
                    design_iBAQ,
                    output_dir_sp_iBAQ_PCA,
                    paste0("iBAQ_", norm_method)
                )
            }
            if (file.exists(norm_file_iBAQ_MC)) {
                process_plot_pca(
                    norm_file_iBAQ_MC,
                    design_iBAQ_MC,
                    output_dir_sp_iBAQ_MC_PCA,
                    paste0("iBAQ_MC_", norm_method)
                )
            }
        }
        message("- Correlation analysis for ", species, " completed.")
        message("- Plots for ", species, " saved to ", output_dir_sp)
    }

    invisible(NULL)
}
# ============================================================================
# TEST WRAPPER
# ============================================================================

#' Test Analysis with Default Parameters
#'
#' @description
#' Convenience function for running complete analysis with pre-configured
#' paths. Useful for development and testing.
#'
#' @return NULL (invisibly). Runs main_analysis() with default paths
test_analysis <- function() {
    main_analysis(
        "Analyses/RI_data.tsv",
        "Analyses/iBAQ_data.tsv",
        "Analyses/iBAQ_mc_data.tsv",
        "DATA/Proteome/sample_mapping.xlsx",
        "Analyses"
    )
}

# ============================================================================
# EXECUTION
# ============================================================================

# Run analysis when script is executed in interactive mode
if (interactive()) {
    test_analysis()
}
