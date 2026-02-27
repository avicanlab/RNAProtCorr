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
    select_columns <- c("Species", "New_locus_tag", "Old_locus_tag")
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
        mutate(
            New_locus_tag = if_else(is.na(New_locus_tag) | New_locus_tag == "N/A", Old_locus_tag, New_locus_tag)
        )
    data %>%
        select(-Old_locus_tag) %>%
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
plot_correlation <- function(correlation_df, mean_protQ_col, protQ_name, species, output_path) {
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

    # Force Ctrl first, then alphabetical
    treatments <- unique(correlation_df$Treatment)
    ordered_treatments <- c("Ctrl", sort(treatments[treatments != "Ctrl"]))
    correlation_df$Treatment <- factor(correlation_df$Treatment, levels = ordered_treatments)
    label_df$Treatment <- factor(label_df$Treatment, levels = ordered_treatments)

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
            title = species
        ) +
        theme_bw(base_size = 10) +
        theme(
            plot.title       = element_text(face = "italic"),
            strip.background = element_blank(),
            strip.text       = element_text(hjust = 0, face = "bold"),
            panel.grid.minor = element_blank(),
            legend.position  = "right"
        )

    filename <- paste0(output_path, "/", paste(species, "TPM", protQ_name, "correlation_panels", sep = "_"))

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
# ============================================================
# Process one species: computes correlations and plots
# ============================================================

process_species <- function(tpm_file, topN_data, iBAQ_data, iBAQ_mc_data, output_path) {
    species_name <- tools::file_path_sans_ext(str_extract(basename(tpm_file), "^[A-Za-z]+_[A-Za-z]+"))
    message("Processing species: ", species_name)

    # --- TPM for this species ---
    tpm_species <- process_tpm(list(tpm_file)) %>%
        filter(!is.na(Protein_id)) %>%
        log2_transform(value_col = "TPM") %>%
        rename_at(vars(mean_log2), ~"mean_Log2_TPM") %>%
        mutate(Protein_id = sub("^pi", "PI", Protein_id, ignore.case = FALSE)) # For Y.pseudotuberculosis

    # --- Filter protQ tables to this species ---
    topN_species <- topN_data %>% filter(Species == species_name)
    iBAQ_species <- iBAQ_data %>% filter(Species == species_name)
    iBAQ_mc_species <- iBAQ_mc_data %>% filter(Species == species_name)

    # Guard: skip if no matching rows
    if (nrow(topN_species) == 0 & nrow(iBAQ_species) == 0 & nrow(iBAQ_mc_species) == 0) {
        warning("No data found for species: ", species_name, " — skipping.")
        return(invisible(NULL))
    }

    # --- Total Intensity ---
    intensity_species <- topN_species %>%
        log2_transform(value_col = "total_intensity") %>%
        rename_at(vars(mean_log2), ~"mean_Log2_intensity")

    corr_intensity <- calculate_correlation(tpm_species, intensity_species, "mean_Log2_intensity")
    plot_correlation(corr_intensity, "mean_Log2_intensity", "Intensity", species_name, output_path)

    # --- TopN ---
    topN_species <- topN_species %>%
        log2_transform(value_col = "relative_TopN") %>%
        rename_at(vars(mean_log2), ~"mean_Log2_TopN")

    corr_topN <- calculate_correlation(tpm_species, topN_species, "mean_Log2_TopN")
    plot_correlation(corr_topN, "mean_Log2_TopN", "TopN", species_name, output_path)

    # --- iBAQ ---
    iBAQ_species <- iBAQ_species %>%
        log2_transform(value_col = "iBAQ") %>%
        rename_at(vars(mean_log2), ~"mean_Log2_iBAQ")

    corr_iBAQ <- calculate_correlation(tpm_species, iBAQ_species, "mean_Log2_iBAQ")
    plot_correlation(corr_iBAQ, "mean_Log2_iBAQ", "iBAQ", species_name, output_path)

    # --- iBAQ MC ---
    iBAQ_mc_species <- iBAQ_mc_species %>%
        log2_transform(value_col = "iBAQ") %>%
        rename_at(vars(mean_log2), ~"mean_Log2_iBAQ")

    corr_iBAQ_mc <- calculate_correlation(tpm_species, iBAQ_mc_species, "mean_Log2_iBAQ")
    plot_correlation(corr_iBAQ_mc, "mean_Log2_iBAQ", "iBAQ_MC", species_name, output_path)

    message("Done: ", species_name)
}


# ============================================================
# Main: iterates over all TPM files and species
# ============================================================

main_analysis <- function(TPM_PATH, TOPN_PATH, IBAQ_PATH, IBAQ_MC_PATH, output_path) {
    # --- Load protQ tables once (they contain all species) ---
    topN_data <- read_tsv(TOPN_PATH,
        col_names = TRUE,
        show_col_types = FALSE, skip_empty_rows = TRUE
    )
    iBAQ_data <- read_tsv(IBAQ_PATH,
        col_names = TRUE,
        show_col_types = FALSE, skip_empty_rows = TRUE
    )
    iBAQ_mc_data <- read_tsv(IBAQ_MC_PATH,
        col_names = TRUE,
        show_col_types = FALSE, skip_empty_rows = TRUE
    )

    # --- List all TPM files (one per species) ---
    tpm_files <- list.files(TPM_PATH, pattern = "\\.xlsx$", full.names = TRUE)

    # --- Process each species ---
    lapply(tpm_files, function(tpm_file) {
        tryCatch(
            process_species(tpm_file, topN_data, iBAQ_data, iBAQ_mc_data, output_path),
            error = function(e) {
                warning("Failed for file: ", basename(tpm_file), "\n  Error: ", e$message)
            }
        )
    })

    message("All species processed.")
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
