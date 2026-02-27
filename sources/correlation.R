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

library(MASS)
library(readxl)
library(tidyverse)

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
    select_columns <- c("Species", "New_locus_tag")
    # Build pattern: ABBV_Treatment_Replicate (GE) - TPM
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
        rename_at("New_locus_tag", ~"Protein_id") %>%
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

# 2. DATA TRANSFORMATION FUNCTIONS
log2_transform <- function(data, value_col) {
    # ------------------------------------------------------------
    # 2. COMPUTE MEANS LOG2 PER PROTEIN × CONDITION
    # ------------------------------------------------------------

    mean_df <- data %>%
        group_by(Protein_id, Treatment) %>%
        summarise(mean_log2 = mean(log2(!!sym(value_col) + 1), na.rm = TRUE), .groups = "drop")
}

# 3. CORRELATION ANALYSIS FUNCTIONS
calculate_correlation <- function(tpm_data, protQ_data, mean_protQ_col) {
    # ------------------------------------------------------------
    # 4. PEARSON R AND P-VALUE PER CONDITION
    # ------------------------------------------------------------
    df <- inner_join(tpm_data, protQ_data, by = c("Protein_id", "Treatment"))

    stats_df <- df %>%
        group_by(Treatment) %>%
        summarise(
            R       = cor(mean_Log2_TPM, .data[[mean_protQ_col]], method = "pearson", use = "complete.obs"),
            p_value = cor.test(mean_Log2_TPM, .data[[mean_protQ_col]], method = "pearson")$p.value,
            .groups = "drop"
        ) %>%
        mutate(
            label = sprintf("R: %.2f\np-value: %.2e", R, p_value)
        )

    df <- left_join(df, stats_df, by = "Treatment")
}

# 4. VISUALIZATION FUNCTIONS
plot_correlation <- function(correlation_df, mean_protQ_col, output_path, protQ_name = NULL) {
    # ============================================================
    # Scatter plots: Mean TPM (log2) vs Mean iBAQ (log2) per condition
    # Pearson R and p-value displayed on each panel
    # ============================================================

    # DENSITY (hex-bin counts) FOR COLOUR SCALE
    correlation_df$density <- NA_real_

    for (trt in unique(correlation_df$Treatment)) {
        idx <- correlation_df$Treatment == trt
        x <- correlation_df$mean_Log2_TPM[idx]
        y <- correlation_df[[mean_protQ_col]][idx]

        # Extra safety check
        keep <- is.finite(x) & is.finite(y)

        dens <- kde2d(x[keep], y[keep],
            n = 50,
            lims = c(range(x[keep]), range(y[keep]))
        )

        ix <- findInterval(x, dens$x)
        iy <- findInterval(y, dens$y)
        ix <- pmin(pmax(ix, 1), 50)
        iy <- pmin(pmax(iy, 1), 50)

        correlation_df$density[idx] <- mapply(function(i, j) dens$z[i, j], ix, iy)
    }

    # LABEL DATA FRAME FOR ANNOTATION
    label_df <- correlation_df %>%
        group_by(Treatment) %>%
        summarise(
            x = min(mean_Log2_TPM, na.rm = TRUE),
            y = max(.data[[mean_protQ_col]], na.rm = TRUE),
            label = unique(label),
            .groups = "drop"
        )

    # PLOT
    p <- ggplot(correlation_df, aes(x = mean_Log2_TPM, y = !!sym(mean_protQ_col))) +
        geom_point(aes(colour = density), alpha = 0.6, size = 1.2) +
        scale_colour_gradient2(low = "#c8d8e8", mid = "#1a006e", high = "red", name = "counts") +
        geom_smooth(method = "lm", se = FALSE, colour = "red", linewidth = 0.8) +
        geom_text(
            data = label_df,
            aes(x = x, y = y, label = label),
            hjust = 0, vjust = 1,
            size = 2.8,
            family = "mono"
        ) +
        facet_wrap(~Treatment, scales = "free", ncol = 6) +
        labs(
            x     = "Mean TPM (log2)",
            y     = "Mean iBAQ (log2)",
            title = "S. aureus"
        ) +
        theme_bw(base_size = 10) +
        theme(
            plot.title       = element_text(face = "italic"),
            strip.background = element_blank(),
            strip.text       = element_text(hjust = 0, face = "bold"),
            panel.grid.minor = element_blank(),
            legend.position  = "right"
        )

    if (!is.null(protQ_name)) {
        filename <- paste0(output_path, "/", paste("TPM", protQ_name, "correlation_panels", sep = "_"))
    } else {
        filename <- paste0(output_path, "/", paste("TPM", mean_protQ_col, "correlation_panels", sep = "_"))
    }

    ggsave(paste0(filename, ".pdf"),
        plot = p,
        width = 24, height = 6, units = "in", dpi = 300
    )

    message("Saved: ", paste0(filename, ".pdf"))

    ggsave(paste0(filename, ".png"),
        plot = p,
        width = 24, height = 6, units = "in", dpi = 300
    )

    message("Saved: ", paste0(filename, ".png"))
}

# 5. MAIN ORCHESTRATION FUNCTION
main_analysis <- function(TPM_PATH, TOPN_PATH, IBAQ_PATH, IBAQ_MC_PATH, output_path) {
    # Load data
    tpm_files <- list.files(TPM_PATH, pattern = "\\.xlsx$", full.names = TRUE)

    tpm_data <- process_tpm(tpm_files) %>%
        filter(!is.na(Protein_id)) %>%
        log2_transform(value_col = "TPM") %>%
        rename_at(vars(mean_log2), ~"mean_Log2_TPM")

    topN_data <- read_tsv(
        TOPN_PATH,
        col_names = TRUE,
        show_col_types = FALSE,
        skip_empty_rows = TRUE
    )

    total_intensity_data <- topN_data %>%
        log2_transform(value_col = "total_intensity") %>%
        rename_at(vars(mean_log2), ~"mean_Log2_intensity")

    corr_tpm_intensity_df <- calculate_correlation(tpm_data, total_intensity_data, "mean_Log2_intensity")
    plot_correlation(corr_tpm_intensity_df, "mean_Log2_intensity", output_path, protQ_name = "Intensity")

    topN_data <- topN_data %>%
        log2_transform(value_col = "relative_TopN") %>%
        rename_at(vars(mean_log2), ~"mean_Log2_TopN")

    corr_tpm_topN_df <- calculate_correlation(tpm_data, topN_data, "mean_Log2_TopN")
    plot_correlation(corr_tpm_topN_df, "mean_Log2_TopN", output_path, protQ_name = "TopN")

    iBAQ_data <- read_tsv(
        IBAQ_PATH,
        col_names = TRUE,
        show_col_types = FALSE,
        skip_empty_rows = TRUE
    ) %>%
        log2_transform(value_col = "iBAQ") %>%
        rename_at(vars(mean_log2), ~"mean_Log2_iBAQ")

    corr_tpm_iBAQ_df <- calculate_correlation(tpm_data, iBAQ_data, "mean_Log2_iBAQ")
    plot_correlation(corr_tpm_iBAQ_df, "mean_Log2_iBAQ", output_path, protQ_name = "iBAQ")

    iBAQ_mc_data <- read_tsv(
        IBAQ_MC_PATH,
        col_names = TRUE,
        show_col_types = FALSE,
        skip_empty_rows = TRUE
    ) %>%
        log2_transform(value_col = "iBAQ") %>%
        rename_at(vars(mean_log2), ~"mean_Log2_iBAQ")
    corr_tpm_iBAQ_mc_df <- calculate_correlation(tpm_data, iBAQ_mc_data, "mean_Log2_iBAQ")
    plot_correlation(corr_tpm_iBAQ_mc_df, "mean_Log2_iBAQ", output_path, protQ_name = "iBAQ_MC")


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
        "Analyses/iBAQ_mc_data.tsv",
        "Analyses"
    )
}

# Execute test when script runs independently
if (interactive()) {
    test_analysis()
}
