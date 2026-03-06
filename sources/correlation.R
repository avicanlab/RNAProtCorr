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
library(ggstar)

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
    "Ctrl", "As", "Bs", "Mig", "Li", "Nd", "Ns", "Oss", "Oxs", "Sp", "Tm"
)

treatment_to_exposure_time <- c(
    "Ctrl" = 10,
    "As"   = 10,
    "Bs"   = 10,
    "Hyp"  = 240,
    "Li"   = 10,
    "Nd"   = 30,
    "Ns"   = 10,
    "Oss"  = 10,
    "Oxs"  = 10,
    "Sp"   = 241,
    "Tm"   = 20
)


#' Exposure Time Label Mapping
#' Maps numeric exposure times (minutes) to display labels
EXPOSURE_TIME_LABELS <- c(
    "10" = "10",
    "20" = "20 min",
    "30" = "30 min",
    "240" = "4 hours",
    "241" = "\u22654 hours"
)

#' Measurement Display Labels
MEASUREMENT_LABELS <- c(
    "corr_intensity" = "Intensity",
    "corr_RI"        = "RI",
    "corr_iBAQ"      = "iBAQ",
    "corr_iBAQ_mc"   = "iBAQ MC"
)

#' Treatment Shape Mapping (ggstar)
TREATMENT_SHAPES <- c(
    "As"   = 1,
    "Bs"   = 14,
    "Ctrl" = 15,
    "Hyp"  = 11,
    "Li"   = 7,
    "Nd"   = 8,
    "Ns"   = 4,
    "Oss"  = 22,
    "Oxs"  = 21,
    "Sp"   = 23,
    "Tm"   = 2
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
        filter(.data[[value_col]] > 0) %>%
        group_by(Protein_id, Treatment) %>%
        summarise(
            mean_log2 = mean(log2(!!sym(value_col)), na.rm = TRUE),
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
# CORRELATION VS EXPOSURE TIME ANALYSIS
# ============================================================================

# --- Data Preparation --------------------------------------------------------

#' Build Exposure Time Dataframe
#'
#' @description
#' Extracts per-treatment Pearson R values from all species and measurements,
#' then maps treatments to numeric exposure times and ordered display labels.
#'
#' @param results Named list. Output of main_analysis() loop. Each entry must
#'   contain: species_name (chr) and corr_* tibbles with Treatment and R columns.
#'
#' @return Tibble with columns:
#'   Treatment (chr), R (dbl), species (chr), measurement (chr),
#'   exposure_time (dbl), exposure_label (ordered factor)
build_exposure_time_df <- function(results) {
    map_dfr(results, function(sp) {
        map_dfr(names(MEASUREMENT_LABELS), function(meas_name) {
            sp[[meas_name]] %>%
                distinct(Treatment, R) %>%
                mutate(
                    species      = sp$species_name,
                    measurement  = meas_name
                )
        })
    }) %>%
        mutate(
            exposure_time = treatment_to_exposure_time[Treatment],
            exposure_label = factor(
                EXPOSURE_TIME_LABELS[as.character(exposure_time)],
                levels = EXPOSURE_TIME_LABELS
            )
        )
}

# --- Statistics --------------------------------------------------------------

#' Compute Meta-Correlation per Measurement
#'
#' @description
#' Computes a global Pearson correlation between numeric exposure time and
#' mRNA-protein R values across all species and treatments, per measurement.
#'
#' @param df Tibble. Output of build_exposure_time_df(), must contain
#'   columns: measurement (chr), exposure_time (dbl), R (dbl).
#'
#' @return Tibble with columns:
#'   measurement (chr), meta_R (dbl), meta_pval (dbl), label (chr)
compute_exposure_correlations <- function(df) {
    df %>%
        group_by(measurement) %>%
        summarise(
            meta_R = cor(exposure_time, R, method = "pearson", use = "complete.obs"),
            meta_pval = tryCatch(
                cor.test(exposure_time, R, method = "pearson")$p.value,
                error = function(e) NA_real_
            ),
            .groups = "drop"
        ) %>%
        mutate(label = sprintf("R: %.2f\np: %.2e", meta_R, meta_pval))
}

#' Build Label Position for Top-Right Annotation
#'
#' @description
#' Computes x/y coordinates for placing the correlation label at the
#' top-right corner of the plot, then joins the stat label text.
#'
#' @param df Tibble. Filtered to a single measurement. Must contain
#'   columns: exposure_label (factor), R (dbl).
#' @param stats Tibble. Single-row output of compute_exposure_correlations()
#'   filtered to the same measurement. Must contain: label (chr).
#'
#' @return Single-row tibble with columns: x (chr), y (dbl), label (chr)
build_label_position <- function(df, stats) {
    tibble(
        x = last(levels(df$exposure_label)),
        y = max(df$R, na.rm = TRUE)
    ) %>%
        bind_cols(stats %>% select(label))
}

# --- Plotting ----------------------------------------------------------------

#' Build Single Exposure Time Scatter Plot
#'
#' @description
#' Creates a ggplot scatter plot of Pearson R vs exposure time for one
#' measurement, colored by species and shaped by treatment, with a
#' top-right annotation showing the global meta-correlation.
#'
#' @param df Tibble. Data for one measurement. Output of build_exposure_time_df()
#'   filtered to a single measurement value.
#' @param label_pos Tibble. Single-row tibble from build_label_position().
#' @param title Chr. Plot title (measurement display name).
#'
#' @return ggplot object
build_exposure_plot <- function(df, label_pos, title) {
    ggplot(df, aes(x = exposure_label, y = R, fill = species)) +
        geom_star(
            aes(starshape = Treatment),
            size = 3,
            alpha = 0.8
        ) +
        scale_starshape_manual(values = TREATMENT_SHAPES) +
        geom_text(
            data = label_pos,
            aes(x = x, y = y, label = label),
            hjust = 0.5,
            vjust = 1,
            size = 3,
            family = "mono",
            color = "black",
            inherit.aes = FALSE
        ) +
        labs(
            title  = title,
            x      = "Exposure time",
            y      = "Pearson R",
            fill   = "Species",
            shape  = "Treatment"
        ) +
        theme_minimal() +
        theme(
            plot.title      = element_text(face = "bold"),
            strip.text      = element_text(face = "bold"),
            legend.position = "right"
        )
}

#' Save Plot to PDF and PNG
#'
#' @description
#' Convenience wrapper around ggsave() that saves both PDF and PNG
#' versions of a plot to the specified output path.
#'
#' @param plot ggplot object. The plot to save.
#' @param filepath Chr. Full path without extension.
#' @param width Num. Plot width in inches (default: 10).
#' @param height Num. Plot height in inches (default: 6).
#'
#' @return NULL (invisibly)
save_plot <- function(plot, filepath, width = 10, height = 6) {
    walk(c(".pdf", ".png"), function(ext) {
        path <- paste0(filepath, ext)
        ggsave(path, plot = plot, width = width, height = height, dpi = 300)
        message("Saved: ", path)
    })
}

# --- Orchestration -----------------------------------------------------------

#' Plot and Save All Exposure Time Correlation Figures
#'
#' @description
#' Iterates over all 4 measurements, builds one scatter plot each,
#' and saves PDF + PNG files named correlation_vs_exposure_time_{measurement}.
#'
#' @param df Tibble. Full output of build_exposure_time_df().
#' @param stats_df Tibble. Output of compute_exposure_correlations().
#' @param output_path Chr. Directory for saving output files.
#'
#' @return NULL (invisibly)
plot_all_exposure_figures <- function(df, stats_df, output_path) {
    walk(names(MEASUREMENT_LABELS), function(meas_name) {
        meas_df <- df %>% filter(measurement == meas_name)
        meas_stats <- stats_df %>% filter(measurement == meas_name)
        label_pos <- build_label_position(meas_df, meas_stats)

        p <- build_exposure_plot(
            df        = meas_df,
            label_pos = label_pos,
            title     = MEASUREMENT_LABELS[meas_name]
        )

        save_plot(
            plot     = p,
            filepath = file.path(output_path, paste0("correlation_vs_exposure_time_", meas_name))
        )
    })
}

#' Process Correlation vs Exposure Time for All Species
#'
#' @description
#' Top-level function that builds the exposure time dataframe, computes
#' global meta-correlations, and generates one plot per measurement.
#' Called at the end of main_analysis().
#'
#' @param results Named list. Output of main_analysis() loop. Each entry
#'   must contain: species_name (chr) and corr_* tibbles.
#' @param output_path Chr. Directory for saving output files.
#'
#' @return Tibble. Combined exposure time data for all species and
#'   measurements (invisibly returned for optional downstream use).
process_correlation_vs_exposure_time <- function(results, output_path) {
    df <- build_exposure_time_df(results)
    stats_df <- compute_exposure_correlations(df)

    plot_all_exposure_figures(df, stats_df, output_path)

    invisible(df)
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
process_species <- function(
  tpm_file, intensity_data, RI_data, iBAQ_data, iBAQ_mc_data, output_path
) {
    # Extract species name from filename
    species_name <- tools::file_path_sans_ext(
        str_extract(basename(tpm_file), "^[A-Za-z]+_[A-Za-z]+")
    )
    message("Processing species: ", species_name)
    # Look up species abbreviation
    abbv <- species_abbv_map[species_name]
    if (is.na(abbv)) {
        warning(
            "No abbreviation found for species '", species_name,
            "', skipping file: ", basename(tpm_file)
        )
        return(invisible(NULL))
    }
    # Load and transform TPM data
    tpm_species <- read_tpm(tpm_file, abbv) %>%
        mutate(Species = species_name) %>%
        # Update treatment name for consistency ("Mig" to "Hyp")
        mutate(Treatment = ifelse(Treatment == "Mig", "Hyp", Treatment)) %>%
        # For Y.pseudotuberculosis: prefix locus tags with 'PI'
        mutate(Protein_id = sub("^pi", "PI", Protein_id, ignore.case = FALSE)) %>%
        filter(!is.na(Protein_id)) %>%
        log2_transform(value_col = "TPM") %>%
        rename(mean_Log2_TPM = mean_log2)

    # Filter protein quantification data to current species
    intensity_species <- intensity_data %>%
        filter(Species == species_name)
    RI_species <- RI_data %>%
        filter(Species == species_name)
    iBAQ_species <- iBAQ_data %>%
        filter(Species == species_name)
    iBAQ_mc_species <- iBAQ_mc_data %>%
        filter(Species == species_name)

    # ---- TOTAL INTENSITY ----
    if (nrow(intensity_species) > 0) {
        corr_intensity <- calculate_correlation(tpm_species, intensity_species, "Intensity_meanlog2")
        plot_correlation(corr_intensity, "Intensity_meanlog2", "Intensity", species_name, output_path)
    }

    # ---- TOP-N QUANTIFICATION ----
    if (nrow(RI_species) > 0) {
        corr_RI <- calculate_correlation(tpm_species, RI_species, "RI_meanlog2")
        plot_correlation(corr_RI, "RI_meanlog2", "RI", species_name, output_path)
    }

    # ---- iBAQ QUANTIFICATION ----
    if (nrow(iBAQ_species) > 0) {
        corr_iBAQ <- calculate_correlation(tpm_species, iBAQ_species, "iBAQ_meanlog2")
        plot_correlation(corr_iBAQ, "iBAQ_meanlog2", "iBAQ", species_name, output_path)
    }

    # ---- iBAQ MC QUANTIFICATION ----
    if (nrow(iBAQ_mc_species) > 0) {
        corr_iBAQ_mc <- calculate_correlation(tpm_species, iBAQ_mc_species, "iBAQ_meanlog2")
        plot_correlation(corr_iBAQ_mc, "iBAQ_meanlog2", "iBAQ_MC", species_name, output_path)
    }

    message("Done: ", species_name)
    return(list(
        species_name = species_name,
        tpm_data = tpm_species,
        corr_intensity = corr_intensity,
        corr_RI = corr_RI,
        corr_iBAQ = corr_iBAQ,
        corr_iBAQ_mc = corr_iBAQ_mc
    ))
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
main_analysis <- function(
  TPM_PATH, Intensity_PATH, RI_PATH, IBAQ_PATH, IBAQ_MC_PATH, output_path
) {
    # Load protein quantification reference tables (contain all species)
    intensity_data <- read_tsv(Intensity_PATH,
        col_names = TRUE,
        show_col_types = FALSE,
        skip_empty_rows = TRUE
    )

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
    results <- list()
    lapply(tpm_files, function(tpm_file) {
        tryCatch(
            {
                result <- process_species(
                    tpm_file,
                    intensity_data,
                    RI_data,
                    iBAQ_data,
                    iBAQ_mc_data,
                    output_path
                )
                results[[result$species_name]] <<- result
            },
            error = function(e) {
                warning(
                    "Failed for file: ", basename(tpm_file),
                    "\n  Error: ", e$message
                )
            }
        )
    })

    message("All species processed.")

    process_correlation_vs_exposure_time(results, output_path)

    message("Exposure time correlation analysis complete.")
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
        "Analyses/Intensity_data.tsv",
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
