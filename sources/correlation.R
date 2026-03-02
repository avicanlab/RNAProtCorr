# ============================================================================
# mRNA-Protein Correlation Analysis
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

library(MASS)
library(readxl)
library(tidyverse)

# ============================================================================
# CONFIGURATION
# ============================================================================

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
    "Ctrl", "As", "Bs", "Hyp", "Li", "Nd", "Ns", "Oss", "Oxs", "Sp", "Tm"
)

# ============================================================================
# TPM DATA LOADING AND PROCESSING
# ============================================================================

#' Read TPM Data from Excel
#'
#' @description
#' Loads Transcript Per Million (TPM) expression data from Excel file,
#' extracts relevant columns, and reshapes to long format with treatment
#' and replicate information.
#'
#' @param dataset_path Character. Path to Excel file containing TPM data
#' @param sp_abbv Character. Species abbreviation for column pattern matching
#'
#' @return Tibble with columns: Species, Protein_id, Treatment, Replicate, TPM
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

#' Process Multiple TPM Files
#'
#' @description
#' Combines TPM data from multiple species files, handling species mapping
#' and validation. Reports progress and skips files with unmapped species.
#'
#' @param tpm_files Character vector. Paths to TPM Excel files
#'
#' @return Tibble combining TPM data from all valid files
process_tpm <- function(tpm_files) {
    map_dfr(tpm_files, function(path) {
        # Extract species name from filename (format: Genus_species)
        curr_species <- str_extract(basename(path), "^[A-Za-z]+_[A-Za-z]+")

        # Look up species abbreviation
        abbv <- species_abbv_map[curr_species]
        if (is.na(abbv)) {
            warning(
                "No abbreviation found for species '", curr_species,
                "', skipping file: ", basename(path)
            )
            return(tibble())
        }

        cat("Processing:", basename(path), "| Species:", curr_species, "\n")

        read_tpm(path, abbv) %>%
            mutate(Species = curr_species)
    })
}

# ============================================================================
# DATA TRANSFORMATION
# ============================================================================

#' Log2 Transform by Group
#'
#' @description
#' Computes mean log2-transformed values grouped by Protein_id and Treatment.
#' Adds pseudocount of 1 to handle zero values.
#'
#' @param data Tibble with expression data
#' @param value_col Character. Name of column containing raw values
#'
#' @return Tibble with mean_log2 values per protein-treatment combination
log2_transform <- function(data, value_col) {
    data %>%
        group_by(Protein_id, Treatment) %>%
        summarise(
            mean_log2 = mean(log2(!!sym(value_col) + 1), na.rm = TRUE),
            .groups = "drop"
        )
}

# ============================================================================
# CORRELATION ANALYSIS
# ============================================================================

#' Calculate Pearson Correlations
#'
#' @description
#' Computes Pearson correlation coefficient and p-value between TPM and
#' protein quantification measures within each treatment condition.
#'
#' @param tpm_data Tibble. Log2-transformed TPM data with mean_Log2_TPM column
#' @param protQ_data Tibble. Protein quantification data with treatment groups
#' @param mean_protQ_col Character. Column name for mean protein quantification values
#'
#' @return Tibble with original data plus correlation stats and labels per treatment
calculate_correlation <- function(tpm_data, protQ_data, mean_protQ_col) {
    df <- inner_join(tpm_data, protQ_data, by = c("Protein_id", "Treatment"))

    stats_df <- df %>%
        group_by(Treatment) %>%
        summarise(
            R = cor(
                mean_Log2_TPM,
                .data[[mean_protQ_col]],
                method = "pearson",
                use = "complete.obs"
            ),
            p_value = cor.test(
                mean_Log2_TPM,
                .data[[mean_protQ_col]],
                method = "pearson"
            )$p.value,
            .groups = "drop"
        ) %>%
        mutate(
            label = sprintf("R: %.2f\np-value: %.2e", R, p_value)
        )

    left_join(df, stats_df, by = "Treatment")
}

# ============================================================================
# VISUALIZATION
# ============================================================================

compute_point_density <- function(x, y, n_bins = 50) {
    # Remove non-finite values
    keep <- is.finite(x) & is.finite(y)
    x_clean <- x[keep]
    y_clean <- y[keep]

    # Define grid breaks
    x_breaks <- seq(min(x_clean), max(x_clean), length.out = n_bins + 1)
    y_breaks <- seq(min(y_clean), max(y_clean), length.out = n_bins + 1)

    # Assign each point to a bin
    ix <- pmin(.bincode(x_clean, x_breaks, include.lowest = TRUE), n_bins)
    iy <- pmin(.bincode(y_clean, y_breaks, include.lowest = TRUE), n_bins)

    # Count number of points per bin
    bin_key <- paste(ix, iy, sep = "_")
    bin_counts <- table(bin_key)

    # Map count back to each point
    point_counts <- as.integer(bin_counts[bin_key])

    # Create output vector with NA for non-finite values
    density_out <- rep(NA_real_, length(x))
    density_out[keep] <- point_counts

    density_out
}

#' Plot Correlation Panels
#'
#' @description
#' Creates scatter plot panels for each treatment condition showing correlation
#' between TPM and protein abundance, with density-weighted colors, fitted line,
#' and correlation statistics.
#'
#' @param correlation_df Tibble. Data with columns: mean_Log2_TPM,
#'   Treatment, and protein quantification column specified by mean_protQ_col
#' @param mean_protQ_col Character. Column name for mean protein quantification
#' @param protQ_name Character. Display name for protein quantification measure
#' @param species Character. Species name for plot title
#' @param output_path Character. Directory path for saving plots
#'
#' @return NULL (invisibly). Saves PDF and PNG files to output_path
plot_correlation <- function(correlation_df, mean_protQ_col, protQ_name,
                             species, output_path) {
    # Compute point density for color scale
    correlation_df <- correlation_df %>%
        nest(data = -Treatment) %>%
        mutate(
            data = map(data, function(df) {
                df %>%
                    mutate(
                        density = compute_point_density(mean_Log2_TPM, !!sym(mean_protQ_col))
                    )
            })
        ) %>%
        unnest(data)

    # Extract correlation labels per treatment
    label_df <- correlation_df %>%
        group_by(Treatment) %>%
        summarise(
            x = min(mean_Log2_TPM, na.rm = TRUE),
            y = max(!!sym(mean_protQ_col), na.rm = TRUE),
            label = first(label),
            .groups = "drop"
        )

    # Order treatments: Control first, then alphabetical
    treatments <- unique(correlation_df$Treatment)
    ordered_treatments <- c("Ctrl", sort(treatments[treatments != "Ctrl"]))
    correlation_df <- correlation_df %>%
        mutate(Treatment = factor(Treatment, levels = ordered_treatments))
    label_df <- label_df %>%
        mutate(Treatment = factor(Treatment, levels = ordered_treatments))

    # Create plot
    p <- ggplot(
        correlation_df,
        aes(x = mean_Log2_TPM, y = !!sym(mean_protQ_col))
    ) +
        geom_point(aes(colour = density), alpha = 0.6, size = 1.2) +
        scale_colour_gradientn(
            colours = c(
                "#c8d8e8",
                "#a8d8ea",
                "#8B7BB5",
                "#6B3A6B",
                "#8B2020",
                "#C0392B",
                "#e8c8cc",
                "#c8d8e8"
            ),
            name = "counts"
        ) +
        geom_smooth(
            method = "lm",
            se = FALSE,
            colour = "red",
            linewidth = 0.8
        ) +
        geom_text(
            data = label_df,
            aes(x = x, y = y, label = label),
            hjust = 0, vjust = 1,
            size = 2.8,
            family = "mono"
        ) +
        facet_wrap(~Treatment, scales = "free", ncol = 6) +
        labs(
            x = "Mean TPM (log2)",
            y = paste("Mean", protQ_name, "(log2)", sep = " "),
            title = species
        ) +
        theme_bw(base_size = 10) +
        theme(
            plot.title = element_text(face = "italic"),
            strip.background = element_blank(),
            strip.text = element_text(hjust = 0, face = "bold"),
            panel.grid.minor = element_blank(),
            legend.position = "right"
        )

    # Construct output filename
    filename <- paste0(
        output_path, "/",
        paste(species, "TPM", protQ_name, "correlation_panels", sep = "_")
    )

    # Save plots
    ggsave(
        paste0(filename, ".pdf"),
        plot = p,
        width = 24, height = 6, units = "in", dpi = 300
    )
    message("Saved: ", paste0(filename, ".pdf"))

    ggsave(
        paste0(filename, ".png"),
        plot = p,
        width = 24, height = 6, units = "in", dpi = 300
    )
    message("Saved: ", paste0(filename, ".png"))
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
#' @param tpm_file Character. Path to species TPM Excel file
#' @param topN_data Tibble. TopN quantification data (all species)
#' @param iBAQ_data Tibble. iBAQ quantification data (all species)
#' @param iBAQ_mc_data Tibble. iBAQ MC quantification data (all species)
#' @param output_path Character. Directory for saving plots
#'
#' @return NULL (invisibly). Generates plots and console messages
process_species <- function(tpm_file, RI_data, iBAQ_data, iBAQ_mc_data, output_path) {
    # Extract species name from filename
    species_name <- tools::file_path_sans_ext(
        str_extract(basename(tpm_file), "^[A-Za-z]+_[A-Za-z]+")
    )
    message("Processing species: ", species_name)

    # Load and transform TPM data
    tpm_species <- process_tpm(list(tpm_file)) %>%
        filter(!is.na(Protein_id)) %>%
        log2_transform(value_col = "TPM") %>%
        rename(mean_Log2_TPM = mean_log2) %>%
        # For Y.pseudotuberculosis: prefix locus tags with 'PI'
        mutate(Protein_id = sub("^pi", "PI", Protein_id, ignore.case = FALSE))

    # Filter protein quantification data to current species
    RI_species <- RI_data %>%
        filter(Species == species_name)
    iBAQ_species <- iBAQ_data %>%
        filter(Species == species_name)
    iBAQ_mc_species <- iBAQ_mc_data %>%
        filter(Species == species_name)

    # Guard: skip if no matched data
    if (nrow(RI_species) == 0 &&
        nrow(iBAQ_species) == 0 &&
        nrow(iBAQ_mc_species) == 0) {
        warning("No data found for species: ", species_name, " — skipping.")
        return(invisible(NULL))
    }

    # ---- TOTAL INTENSITY ----
    if (nrow(RI_species) > 0) {
        intensity_species <- RI_species %>%
            log2_transform(value_col = "total_intensity") %>%
            rename(mean_Log2_intensity = mean_log2)

        corr_intensity <- calculate_correlation(tpm_species, intensity_species, "mean_Log2_intensity")
        plot_correlation(corr_intensity, "mean_Log2_intensity", "Intensity", species_name, output_path)
    }

    # ---- TOP-N QUANTIFICATION ----
    if (nrow(RI_species) > 0) {
        RI_transformed <- RI_species %>%
            log2_transform(value_col = "RI") %>%
            rename(mean_Log2_RI = mean_log2)

        corr_RI <- calculate_correlation(tpm_species, RI_transformed, "mean_Log2_RI")
        plot_correlation(corr_RI, "mean_Log2_RI", "RI", species_name, output_path)
    }

    # ---- iBAQ QUANTIFICATION ----
    if (nrow(iBAQ_species) > 0) {
        iBAQ_transformed <- iBAQ_species %>%
            log2_transform(value_col = "iBAQ") %>%
            rename(mean_Log2_iBAQ = mean_log2)

        corr_iBAQ <- calculate_correlation(tpm_species, iBAQ_transformed, "mean_Log2_iBAQ")
        plot_correlation(corr_iBAQ, "mean_Log2_iBAQ", "iBAQ", species_name, output_path)
    }

    # ---- iBAQ MC QUANTIFICATION ----
    if (nrow(iBAQ_mc_species) > 0) {
        iBAQ_mc_transformed <- iBAQ_mc_species %>%
            log2_transform(value_col = "iBAQ") %>%
            rename(mean_Log2_iBAQ = mean_log2)

        corr_iBAQ_mc <- calculate_correlation(tpm_species, iBAQ_mc_transformed, "mean_Log2_iBAQ")
        plot_correlation(corr_iBAQ_mc, "mean_Log2_iBAQ", "iBAQ_MC", species_name, output_path)
    }

    message("Done: ", species_name)
}

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
main_analysis <- function(TPM_PATH, RI_PATH, IBAQ_PATH, IBAQ_MC_PATH, output_path) {
    # Load protein quantification reference tables (contain all species)
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

    # Find all TPM files (one per species)
    tpm_files <- list.files(TPM_PATH, pattern = "\\.xlsx$", full.names = TRUE)

    # Process each species with error handling
    lapply(tpm_files, function(tpm_file) {
        tryCatch(
            process_species(tpm_file, RI_data, iBAQ_data, iBAQ_mc_data, output_path),
            error = function(e) {
                warning(
                    "Failed for file: ", basename(tpm_file),
                    "\n  Error: ", e$message
                )
            }
        )
    })

    message("All species processed.")
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
        "DATA/Read_counts/",
        "Analyses/RI_data.tsv",
        "Analyses/iBAQ_data.tsv",
        "Analyses/iBAQ_mc_data.tsv",
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
