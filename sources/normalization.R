# ============================================================================
# Protein Normalization Pipeline
# ============================================================================
#
# Description:
#   Pipeline for normalizing protein quantification data (RI, iBAQ, iBAQ_MC)
#   across multiple species and experimental conditions using NormalyzerDE.
#   Generates normalized datasets and PCA visualizations for method selection.
#
# Expected input columns: Species, Replicate, Treatment, Protein_id, <protQ>
# ============================================================================

library(tidyverse)
library(NormalyzerDE)

NORM_METHODS <- c("Log2", "VSN", "CycLoess", "Quantile", "Median", "Mean", "GI", "RLR")


# ============================================================================
# NORMALIZATION HELPERS
# ============================================================================

#' Create a NormalyzerDE design file matched to a data matrix
#'
#' Reads column names from the wide-format matrix and builds a design table
#' where group is the treatment (replicate prefix stripped from sample name).
#'
#' @param data_path  Path to wide-format TSV (Protein_id + Replicate_Treatment cols).
#' @param tmp_dir    Directory for temp output. Default tempdir().
#' @return Path to the written design TSV.
create_design_for_matrix <- function(data_path, tmp_dir = tempdir()) {
    sample_cols <- read.table(data_path,
        sep = "\t", header = TRUE,
        check.names = FALSE, nrows = 1
    ) %>%
        select(-Protein_id) %>%
        colnames()

    design_df <- tibble(
        sample = sample_cols,
        group  = str_remove(sample_cols, "^R\\d+_") # R1_Ctrl -> Ctrl
    )

    temp_file <- tempfile(pattern = "design_", tmpdir = tmp_dir, fileext = ".tsv")
    write.table(design_df, temp_file, sep = "\t", row.names = FALSE, quote = FALSE)
    temp_file
}


#' Pivot long-format protein data to wide format for NormalyzerDE
#'
#' Creates a Replicate_Treatment sample column, pivots to wide format, and
#' sorts columns alphabetically for consistent ordering across species.
#'
#' @param data    Long-format data frame: Protein_id, Replicate, Treatment, <protQ>.
#' @param protQ   Name of the quantification column (e.g. "RI", "iBAQ").
#' @param tmp_dir Directory for temp output. Default tempdir().
#' @return Path to the written wide-format TSV.
transform_data <- function(data, protQ, tmp_dir = tempdir()) {
    message("Transforming ", protQ, " data...")

    data_wide <- data %>%
        mutate(sample = paste0(Replicate, "_", Treatment)) %>%
        pivot_wider(
            id_cols     = Protein_id,
            names_from  = sample,
            values_from = all_of(protQ)
        ) %>%
        select(Protein_id, sort(colnames(.)[-1]))

    temp_file <- tempfile(
        pattern = paste0("transformed_", protQ, "_"),
        tmpdir  = tmp_dir,
        fileext = ".tsv"
    )
    write.table(data_wide, temp_file, sep = "\t", row.names = FALSE, quote = FALSE)
    temp_file
}


#' Normalize protein data using NormalyzerDE
#'
#' Runs all normalization methods, evaluates them, and writes normalized
#' datasets and diagnostic plots (PDF + PNG) to output_path.
#'
#' @param data_path   Path to wide-format TSV (from transform_data).
#' @param design_path Path to design TSV (from create_design_for_matrix).
#' @param protQ       Quantification label used in output filenames.
#' @param output_path Output directory (created if absent).
#' @return NULL invisibly.
normalize_data <- function(data_path, design_path, protQ, output_path) {
    message("Normalizing ", protQ, " data...")
    dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

    normalyzer(
        jobName = paste0(protQ),
        dataPath = data_path,
        designPath = design_path,
        outputDir = output_path,
        writeReportAsPngs = TRUE,
        zeroToNA = TRUE,
    )

    invisible(NULL)
}


# ============================================================================
# VISUALIZATION
# ============================================================================

#' Compute PCA on a normalized protein matrix
#'
#' Strips non-numeric columns, filters complete cases, and runs PCA on the
#' transposed matrix (samples as observations). Returns coordinates and plot.
#'
#' @param norm_matrix  Wide-format data frame (Protein_id + sample columns).
#' @param design_table Path to design TSV with columns sample and group.
#' @param method_name  Label used in the plot title.
#' @return Named list: df (PCA coordinates tibble), plot (ggplot2 object).
plot_pca <- function(norm_matrix, design_table, method_name) {
    mat <- norm_matrix %>%
        select(-any_of(c("CV", "anovaP"))) %>%
        column_to_rownames("Protein_id") %>%
        select(where(is.numeric)) %>%
        as.matrix()
    mat <- mat[complete.cases(mat), ]

    pca <- prcomp(t(mat), scale. = TRUE, center = TRUE)
    var_exp <- round(summary(pca)$importance[2, 1:2] * 100, 1)

    design_df <- read.table(design_table, header = TRUE, sep = "\t") %>%
        as_tibble() %>%
        mutate(
            Treatment = group,
            Replicate = str_remove(sample, "_[a-zA-Z]+$")
        )

    pca_df <- as.data.frame(pca$x[, 1:2]) %>%
        rownames_to_column("sample") %>%
        left_join(design_df, by = "sample")

    pca_plot <- ggplot(pca_df, aes(
        x = PC1, y = PC2,
        color = Treatment,
        shape = factor(Replicate),
        label = sample
    )) +
        geom_point(size = 3) +
        labs(
            title = paste("PCA \u2014", method_name),
            x     = paste0("PC1 (", var_exp[1], "%)"),
            y     = paste0("PC2 (", var_exp[2], "%)"),
            shape = "Replicate"
        ) +
        theme_bw()

    list(df = pca_df, plot = pca_plot)
}


#' Save PCA plot and coordinate table for a normalized dataset
#'
#' @param norm_path   Path to normalized TSV.
#' @param design_path Path to design TSV.
#' @param output_path Directory for saved outputs.
#' @param method_name Label used in filenames and plot title.
#' @return NULL invisibly.
process_plot_pca <- function(norm_path, design_path, output_path, method_name) {
    norm_matrix <- read.table(norm_path, sep = "\t", header = TRUE, check.names = FALSE)
    pca_res <- plot_pca(norm_matrix, design_path, method_name)

    tsv_out <- file.path(output_path, paste0(method_name, "_PCA_coordinates.tsv"))
    write.table(pca_res$df, file = tsv_out, sep = "\t", row.names = FALSE, quote = FALSE)
    message("PCA coordinates saved to: ", tsv_out)

    png_out <- file.path(output_path, paste0(method_name, ".png"))
    ggsave(filename = png_out, plot = pca_res$plot, width = 6, height = 5)
    message("PCA plot saved to: ", png_out)

    invisible(NULL)
}


# ============================================================================
# MAIN ORCHESTRATION
# ============================================================================

#' Run the full protein normalization and PCA pipeline
#'
#' For each species in the RI dataset:
#'   1. Transforms RI, iBAQ, iBAQ_MC data to wide format.
#'   2. Creates matched NormalyzerDE design files per quantification type.
#'   3. Normalizes with all methods via NormalyzerDE.
#'   4. Generates PCA plots per normalization method.
#'
#' @param RI_PATH      Path to RI data TSV (long format).
#' @param IBAQ_PATH    Path to iBAQ data TSV (long format).
#' @param IBAQ_MC_PATH Path to iBAQ_MC data TSV (long format).
#' @param output_path  Root output directory.
#' @return NULL invisibly.
main_analysis <- function(RI_PATH, IBAQ_PATH, IBAQ_MC_PATH, output_path) {
    # ---- Load data -----------------------------------------------------------
    read_data <- function(path) {
        read_tsv(path, col_names = TRUE, show_col_types = FALSE, skip_empty_rows = TRUE)
    }

    RI_data <- read_data(RI_PATH)
    iBAQ_data <- read_data(IBAQ_PATH)
    iBAQ_mc_data <- read_data(IBAQ_MC_PATH)

    species_list <- unique(RI_data$Species)
    tmp_dir <- tempdir(check = FALSE)
    output_dir <- file.path(output_path, "NormTest")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    # ---- Per-species loop ----------------------------------------------------
    for (species in species_list) {
        message("\n========== Processing species: ", species, " ==========")

        sp_dir <- file.path(output_dir, species)

        # -- Transform to wide format --------------------------------------------
        data_paths <- list(
            RI      = RI_data %>% filter(Species == species) %>% transform_data("RI", tmp_dir),
            iBAQ    = iBAQ_data %>% filter(Species == species) %>% transform_data("iBAQ", tmp_dir),
            iBAQ_MC = iBAQ_mc_data %>% filter(Species == species) %>% transform_data("iBAQ", tmp_dir)
        )

        # -- Create matched design files -----------------------------------------
        design_paths <- map(data_paths, ~ create_design_for_matrix(.x, tmp_dir))

        # -- Normalize -----------------------------------------------------------
        iwalk(data_paths, ~ normalize_data(.x, design_paths[[.y]], .y, sp_dir))


        # -- PCA per normalization method ----------------------------------------
        dirs <- list(
            RI      = file.path(sp_dir, "RI"),
            iBAQ    = file.path(sp_dir, "iBAQ"),
            iBAQ_MC = file.path(sp_dir, "iBAQ_MC")
        )
        walk(dirs, ~ dir.create(file.path(.x, "PCA"), showWarnings = FALSE))
        for (norm_method in NORM_METHODS) {
            message("  Generating PCA for: ", norm_method)

            iwalk(data_paths, function(path, quant_type) {
                norm_file <- file.path(dirs[[quant_type]], paste0(norm_method, "-normalized.txt"))
                if (!file.exists(norm_file)) {
                    warning("Normalized file not found for ", quant_type, " with method ", norm_method, ": ", norm_file)
                    return(invisible(NULL))
                }

                process_plot_pca(
                    norm_path   = norm_file,
                    design_path = design_paths[[quant_type]],
                    output_path = file.path(dirs[[quant_type]], "PCA"),
                    method_name = paste0(quant_type, "_", norm_method)
                )
            })
        }

        message("Species ", species, " complete. Outputs saved to: ", sp_dir)
    }

    invisible(NULL)
}


# ============================================================================
# ENTRY POINT
# ============================================================================

#' Run analysis with default development paths
#' @return NULL invisibly.
test_analysis <- function() {
    main_analysis(
        RI_PATH      = "Analyses/RI_data.tsv",
        IBAQ_PATH    = "Analyses/iBAQ_data.tsv",
        IBAQ_MC_PATH = "Analyses/iBAQ_mc_data.tsv",
        output_path  = "Analyses"
    )
}

if (interactive()) {
    test_analysis()
}
