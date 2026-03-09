# ============================================================================
# Differential Expression Analysis — DEG (TPM) and DEP (Proteomics)
# ============================================================================
#
# Description:
#   Calls differentially expressed genes (DEG) and proteins (DEP) per
#   stress treatment vs Control, then plots the percentage per condition.
#
#   DEG: Uses pre-computed log2FC and FDR p-value columns from the
#        read_counts Excel file. No limma needed.
#        Threshold: FDR < 0.05 AND |log2FC| > 0.5
#
#   DEP: Runs limma on pre-normalised proteomics data (wide format,
#        one column per Replicate_Treatment).
#        Threshold: BH adj. p-value < 0.05 AND |log2FC| > 0.5
#
# Inputs:
#   - Read counts Excel : DATA/Read_counts/<Species>.xlsx
#     Key columns (per treatment):
#       "{trt} vs. Ctrl - Log fold change"
#       "{trt} vs. Ctrl - FDR p-value"
#       "{trt} vs. Ctrl - Max group means"
#     Gene ID columns: New_locus_tag, Old_locus_tag
#     Raw count columns: R1_As, R1_Ctrl, ... (used by correlation.R only)
#
#   - Normalised proteomics TSV: Analyses/NormTest/<Species>/<Method>.tsv
#     Columns: Protein_id, R1_As, R1_Bs, R1_Ctrl, ..., R3_Tm
#
# Outputs (per species):
#   - <Species>_DEG_results.tsv
#   - <Species>_DEP_results.tsv
#   - <Species>_DEG_DEP_barplot.pdf / .png
#
# Requirements:
#   BiocManager::install("limma")
#   install.packages("tidyverse")
#   install.packages("readxl")
# ============================================================================

SOURCED_AS_MODULE <- TRUE
source("sources/correlation.R") # species_abbv_map, treatment_list, save_plot

library(limma)
library(readxl)
library(tidyverse)

# ============================================================================
# CONFIGURATION
# ============================================================================

DEG_PVAL_THRESHOLD <- 0.05 # FDR p-value cutoff for DEG
DEP_PVAL_THRESHOLD <- 0.05 # BH adj. p-value cutoff for DEP
DE_FC_THRESHOLD <- 0.5 # |log2FC| cutoff for both DEG and DEP

# Stress treatments to test (all vs Ctrl)
STRESS_TREATMENTS <- c("As", "Bs", "Li", "Nd", "Ns", "Oss", "Oxs", "Sp", "Tm")
STRESS_TREATMENTS_RNA <- c(STRESS_TREATMENTS, "Mig")
STRESS_TREATMENTS_PROT <- c(STRESS_TREATMENTS, "Hyp")

PROTQ_METHODS <- c("Intensity", "RI", "iBAQ", "iBAQ_MC")

# Y-axis order matching reference figure (reverse-alpha: Sp top, As bottom)
TREATMENT_ORDER <- rev(sort(STRESS_TREATMENTS_PROT))

# Plot colours matching reference figure
DE_COLOURS <- c("DEGs" = "#89CFF0", "DEPs" = "#F5A623")

# Column name templates in read_counts file
col_lfc <- function(trt) paste0(trt, " vs. Ctrl - Log fold change")
col_fdr <- function(trt) paste0(trt, " vs. Ctrl - FDR p-value")
col_mean <- function(trt) paste0(trt, " vs. Ctrl - Max group means")

# ============================================================================
# DATA LOADING
# ============================================================================

#' Read DEG Results from Read Counts Excel File
#'
#' @description
#' Loads the read counts Excel file, resolves the gene ID from
#' New_locus_tag (falling back to Old_locus_tag when missing/NA),
#' and extracts pre-computed log2FC and FDR p-value columns for each
#' stress treatment vs Ctrl.
#'
#' @param filepath Chr. Path to the read counts Excel file.
#'
#' @return Tibble: Protein_id, Treatment, log2FC, fdr, max_mean
read_deg_data <- function(filepath) {
    raw <- read_excel(filepath, col_names = TRUE)

    # Resolve gene ID: prefer New_locus_tag, fall back to Old_locus_tag
    raw <- raw %>%
        mutate(
            Protein_id = if_else(
                is.na(New_locus_tag) | New_locus_tag == "N/A",
                Old_locus_tag,
                New_locus_tag
            )
        ) %>%
        filter(!is.na(Protein_id))

    # Extract per-treatment columns for treatments present in file
    treatments_present <- intersect(
        STRESS_TREATMENTS_RNA,
        # A treatment is present if its log2FC column exists
        STRESS_TREATMENTS_RNA[map_lgl(STRESS_TREATMENTS_RNA, ~ col_lfc(.x) %in% colnames(raw))]
    )

    map_dfr(treatments_present, function(trt) {
        raw %>%
            dplyr::select(
                Protein_id,
                log2FC   = all_of(col_lfc(trt)),
                fdr      = all_of(col_fdr(trt)),
                max_mean = all_of(col_mean(trt))
            ) %>%
            mutate(
                Treatment = ifelse(trt == "Mig", "Hyp", trt),
                log2FC    = as.numeric(log2FC),
                fdr       = as.numeric(fdr),
                max_mean  = as.numeric(max_mean)
            )
    }) %>%
        dplyr::select(Protein_id, Treatment, log2FC, fdr, max_mean)
}

#' Read Normalised Proteomics TSV
#'
#' @description
#' Loads a wide-format normalised proteomics file and reshapes to long format.
#' Expects columns: Protein_id, R1_As, R1_Bs, ..., R3_Tm.
#'
#' @param filepath Chr. Path to the TSV file.
#'
#' @return Tibble: Protein_id, Replicate, Treatment, value
read_norm_proteomics <- function(filepath) {
    read_tsv(filepath, col_names = TRUE, show_col_types = FALSE) %>%
        pivot_longer(
            cols      = -Protein_id,
            names_to  = "sample",
            values_to = "value"
        ) %>%
        mutate(
            Replicate = str_extract(sample, "^R\\d+"),
            Treatment = str_extract(sample, "[A-Za-z]+$")
        ) %>%
        dplyr::select(Protein_id, Replicate, Treatment, value)
}

# ============================================================================
# DEG CALLING
# ============================================================================

#' Call DEGs from Pre-Computed Fold Change and FDR Columns
#'
#' @description
#' Applies FDR < DEG_PVAL_THRESHOLD AND |log2FC| > DE_FC_THRESHOLD
#' to flag differentially expressed genes. Summarises per treatment.
#'
#' @param deg_data       Tibble. Output of read_deg_data().
#' @param pval_threshold Num. FDR threshold. Default DEG_PVAL_THRESHOLD.
#' @param fc_threshold   Num. |log2FC| threshold. Default DE_FC_THRESHOLD.
#'
#' @return Tibble: Treatment, n_de, n_total, pct_de, de_ids (list-col)
call_deg <- function(
  deg_data,
  pval_threshold = DEG_PVAL_THRESHOLD,
  fc_threshold = DE_FC_THRESHOLD
) {
    deg_data %>%
        filter(!is.na(log2FC), !is.na(fdr)) %>%
        group_by(Treatment) %>%
        summarise(
            n_de = sum(fdr < pval_threshold & abs(log2FC) > fc_threshold,
                na.rm = TRUE
            ),
            n_total = n(),
            pct_de = 100 * n_de / n_total,
            de_ids = list(Protein_id[fdr < pval_threshold &
                abs(log2FC) > fc_threshold &
                !is.na(fdr)]),
            .groups = "drop"
        )
}

# ============================================================================
# LIMMA (DEP)
# ============================================================================

#' Build Wide Expression Matrix from Long-Format Data
#'
#' Pivots to wide matrix (features x samples). Drops rows with any NA
#' (complete cases only). Returns matrix and aligned sample metadata.
#'
#' @param data       Tibble: Protein_id, Replicate, Treatment, value.
#' @param treatments Chr vector. Treatments to include (NULL = all present).
#'
#' @return List:
#'   - mat      : numeric matrix, rownames = Protein_id
#'   - metadata : tibble with columns sample, Replicate, Treatment
build_expression_matrix <- function(data, treatments = NULL) {
    if (!is.null(treatments)) {
        data <- data %>% filter(Treatment %in% treatments)
    }

    wide <- data %>%
        mutate(sample = paste(Replicate, Treatment, sep = "_")) %>%
        dplyr::select(Protein_id, sample, value) %>%
        group_by(Protein_id, sample) %>%
        summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
        pivot_wider(names_from = sample, values_from = value) %>%
        column_to_rownames("Protein_id")

    wide <- wide[complete.cases(wide), ]
    mat <- as.matrix(wide)

    metadata <- tibble(sample = colnames(mat)) %>%
        mutate(
            Replicate = str_extract(sample, "^R\\d+"),
            Treatment = str_extract(sample, "[A-Za-z]+$")
        )

    list(mat = mat, metadata = metadata)
}

#' Run limma DE Analysis: All Stress Treatments vs Ctrl
#'
#' @description
#' Builds a design matrix with Replicate as a blocking factor.
#' Runs lmFit + eBayes on pre-normalised log2 proteomics data.
#' Applies BH correction per contrast.
#'
#' @param mat            Numeric matrix. Rows = proteins, cols = samples.
#' @param metadata       Tibble. Columns: sample, Replicate, Treatment.
#' @param pval_threshold Num. BH adj. p-value threshold.
#' @param fc_threshold   Num. |log2FC| threshold.
#'
#' @return Tibble: Treatment, n_de, n_total, pct_de, de_ids (list-col).
#'   Returns NULL if Ctrl absent or no stress treatments found.
run_limma_dep <- function(
  mat,
  metadata,
  pval_threshold = DEP_PVAL_THRESHOLD,
  fc_threshold = DE_FC_THRESHOLD
) {
    treatments_present <- unique(metadata$Treatment)

    if (!"Ctrl" %in% treatments_present) {
        warning("No Ctrl samples found — cannot run DEP analysis.")
        return(NULL)
    }

    stress_present <- intersect(STRESS_TREATMENTS_PROT, treatments_present)
    if (length(stress_present) == 0) {
        warning("No stress treatments found in proteomics data.")
        return(NULL)
    }

    metadata <- metadata %>%
        mutate(
            Treatment = factor(Treatment, levels = c("Ctrl", stress_present)),
            Replicate = factor(Replicate)
        )

    # ~ 0 + Treatment + Replicate: Replicate blocks inter-replicate variation
    design <- model.matrix(~ 0 + Treatment + Replicate, data = metadata)
    colnames(design) <- str_replace(colnames(design), "^Treatment", "")

    fit <- lmFit(mat, design)

    contrast_mat <- makeContrasts(
        contrasts = paste0(stress_present, " - Ctrl"),
        levels    = design
    )
    colnames(contrast_mat) <- stress_present

    fit2 <- eBayes(contrasts.fit(fit, contrast_mat))

    map_dfr(stress_present, function(trt) {
        tt <- topTable(fit2,
            coef = trt, number = Inf,
            sort.by = "none", adjust.method = "BH"
        )

        sig_mask <- !is.na(tt$adj.P.Val) &
            tt$adj.P.Val < pval_threshold &
            abs(tt$logFC) > fc_threshold

        tibble(
            Treatment = trt,
            n_de      = sum(sig_mask),
            n_total   = nrow(tt),
            pct_de    = 100 * sum(sig_mask) / nrow(tt),
            de_ids    = list(rownames(tt)[sig_mask])
        )
    })
}

# ============================================================================
# VISUALIZATION
# ============================================================================

#' Plot DEG vs DEP Percentage Horizontal Bar Chart
#'
#' @description
#' Reproduces reference figure C: horizontal grouped bars, DEGs (blue) and
#' DEPs (orange). Y-axis in reverse-alphabetical order (Sp top, As bottom).
#'
#' @param deg_df      Tibble. call_deg() output.
#' @param dep_df      Tibble. run_limma_dep() output.
#' @param species     Chr. Species name for italic title.
#' @param output_path Chr. Output directory.
#'
#' @return NULL invisibly.
plot_de_barplot <- function(deg_df, dep_df, species, output_path) {
    plot_df <- bind_rows(
        deg_df %>% mutate(type = "DEGs"),
        dep_df %>% mutate(type = "DEPs")
    ) %>%
        filter(Treatment %in% STRESS_TREATMENTS_PROT) %>%
        mutate(
            Treatment = factor(Treatment, levels = TREATMENT_ORDER),
            type      = factor(type, levels = c("DEGs", "DEPs"))
        )

    x_max <- ceiling(max(plot_df$pct_de, na.rm = TRUE) / 10) * 10
    x_max <- max(x_max, 20)

    p <- ggplot(plot_df, aes(x = pct_de, y = Treatment, fill = type)) +
        geom_col(
            position = position_dodge(width = 0.7, reverse = TRUE),
            width    = 0.6,
            colour   = NA
        ) +
        scale_fill_manual(values = DE_COLOURS, name = NULL) +
        scale_x_continuous(
            limits = c(0, x_max),
            expand = expansion(mult = c(0, 0.02)),
            breaks = seq(0, x_max, by = 20)
        ) +
        scale_y_discrete(drop = FALSE) +
        labs(
            title = bquote(italic(.(gsub("_", " ", species)))),
            x     = "Percentage",
            y     = "Stress Conditions"
        ) +
        theme_bw(base_size = 11) +
        theme(
            plot.title         = element_text(face = "italic", hjust = 0.5, size = 11),
            axis.title         = element_text(size = 10),
            axis.text          = element_text(size = 9),
            panel.grid.major.y = element_blank(),
            panel.grid.minor   = element_blank(),
            panel.border       = element_rect(colour = "grey70"),
            legend.position    = c(0.82, 0.92),
            legend.background  = element_rect(fill = "white", colour = NA),
            legend.key.size    = unit(0.4, "cm"),
            legend.text        = element_text(size = 9)
        )

    save_plot(
        plot     = p,
        filepath = file.path(output_path, paste0(species, "_DEG_DEP_barplot")),
        width    = 5,
        height   = 5
    )

    message("  Plot saved: ", species, "_DEG_DEP_barplot.pdf/.png")
    invisible(NULL)
}

# ============================================================================
# OUTPUT HELPERS
# ============================================================================

#' Write DE Summary Table to TSV
#'
#' @param de_df       Tibble. call_deg() or run_limma_dep() output.
#' @param output_path Chr. Output directory.
#' @param filename    Chr. Output filename.
write_de_table <- function(de_df, output_path, filename) {
    de_df %>%
        dplyr::select(Treatment, n_de, n_total, pct_de) %>%
        write.table(
            file.path(output_path, filename),
            sep = "\t", row.names = FALSE, quote = FALSE
        )
}

# ============================================================================
# SPECIES ORCHESTRATION
# ============================================================================

#' Process DEG and DEP Analysis for One Species
#'
#' @description
#' Loads read counts and normalised proteomics, calls DEGs from pre-computed
#' columns, runs limma for DEPs, writes summary TSVs, and plots the bar chart.
#'
#' @param counts_file    Chr. Path to read counts Excel file.
#' @param norm_path      Chr. Path to normalised proteomics TSV.
#' @param species        Chr. Species name.
#' @param output_path    Chr. Output directory.
#' @param pval_threshold Num. FDR/BH threshold. Default DEG_PVAL_THRESHOLD.
#' @param fc_threshold   Num. |log2FC| threshold. Default DE_FC_THRESHOLD.
#'
#' @return List: species, deg_results, dep_results. NULL on failure.
process_species_de <- function(
  deg_data,
  norm_path,
  species,
  output_path,
  pval_threshold = DEG_PVAL_THRESHOLD,
  fc_threshold = DE_FC_THRESHOLD,
  common_ids_only = TRUE
) {
    # ---- DEP from normalised proteomics via limma ---------------------------
    dep_results <- NULL

    if (!file.exists(norm_path)) {
        warning("  Normalised proteomics not found: ", norm_path)
    } else {
        message("  Loading normalised proteomics: ", basename(norm_path))
        prot_long <- read_norm_proteomics(norm_path) %>%
            filter(
                Treatment %in% c("Ctrl", STRESS_TREATMENTS_PROT),
                is.finite(value)
            )
    }

    if (common_ids_only) {
        common_ids <- Reduce(intersect, list(
            deg_data %>% pull(Protein_id) %>% unique(),
            prot_long %>% pull(Protein_id) %>% unique()
        ))
        message("  Common IDs across all datasets: ", length(common_ids))

        deg_data <- deg_data %>% filter(Protein_id %in% common_ids)
        prot_long <- prot_long %>% filter(Protein_id %in% common_ids)
    }

    deg_results <- call_deg(deg_data, pval_threshold, fc_threshold)

    write_de_table(deg_results, output_path, paste0(species, "_DEG_results.tsv"))
    message(
        "  DEG: ", sum(deg_results$n_de), " total DE genes across ",
        nrow(deg_results), " conditions"
    )

    prot_expr <- build_expression_matrix(prot_long)
    dep_results <- run_limma_dep(
        mat            = prot_expr$mat,
        metadata       = prot_expr$metadata,
        pval_threshold = DEP_PVAL_THRESHOLD,
        fc_threshold   = fc_threshold
    )

    if (!is.null(dep_results)) {
        write_de_table(
            dep_results, output_path,
            paste0(species, "_DEP_results.tsv")
        )
        message(
            "  DEP: ", sum(dep_results$n_de), " total DE proteins across ",
            nrow(dep_results), " conditions"
        )
    }

    # ---- Plot ---------------------------------------------------------------
    if (!is.null(dep_results)) {
        plot_de_barplot(deg_results, dep_results, species, output_path)
    } else {
        message("  Skipping plot — DEP results unavailable.")
    }

    list(species = species, deg_results = deg_results, dep_results = dep_results)
}

# ============================================================================
# MAIN
# ============================================================================

#' Run Full DEG / DEP Pipeline Across All Species
#'
#' @description
#' Discovers species from read counts files, pairs each with its normalised
#' proteomics file, and runs the full analysis.
#'
#' @param counts_path    Chr. Directory containing read counts Excel files.
#' @param norm_base      Chr. Base path for normalised proteomics.
#'   Expected file: <norm_base>/<Species>/<norm_method>.tsv
#' @param norm_method    Chr. Normalisation method name (no .tsv extension).
#' @param output_path    Chr. Output directory.
#' @param pval_threshold Num. FDR/BH threshold. Default 0.05.
#' @param fc_threshold   Num. |log2FC| threshold. Default 0.5.
#'
#' @return Named list of per-species results (invisibly).
main_analysis_de <- function(
  counts_path,
  norm_base,
  norm_method,
  output_path,
  pval_threshold = DEG_PVAL_THRESHOLD,
  fc_threshold = DE_FC_THRESHOLD
) {
    counts_files <- list.files(counts_path, pattern = "\\.xlsx$", full.names = TRUE)

    output_de <- paste(output_path, "DE", sep = "/")
    dir.create(output_de)

    if (length(counts_files) == 0) {
        stop("No read counts Excel files found in: ", counts_path)
    }

    results <- list()

    for (counts_file in counts_files) {
        species <- tools::file_path_sans_ext(
            str_extract(basename(counts_file), "^[A-Za-z]+_[A-Za-z]+")
        )
        # ---- DEG from pre-computed fold changes ---------------------------------
        message("  Loading read counts: ", basename(counts_file))
        deg_data <- read_deg_data(counts_file)

        for (protQ in PROTQ_METHODS) {
            message("Process ", protQ, " for ", species, "...")
            norm_path <- file.path(
                norm_base, species, protQ, paste0(norm_method, "-normalized.txt")
            )
            output_sp_protq <- paste(
                output_de, species, protQ, norm_method,
                sep = "/"
            )
            dir.create(output_sp_protq, recursive = TRUE)
            result <- tryCatch(
                process_species_de(
                    deg_data       = deg_data,
                    norm_path      = norm_path,
                    species        = species,
                    output_path    = output_sp_protq,
                    pval_threshold = pval_threshold,
                    fc_threshold   = fc_threshold
                ),
                error = function(e) {
                    warning("Failed for '", species, "': ", e$message)
                    NULL
                }
            )

            if (!is.null(result)) results[[species]] <- result
        }
    }

    message("\nDE analysis complete. Results written to: ", output_path)
    invisible(results)
}

# ============================================================================
# ENTRY POINT
# ============================================================================

#' Run test analysis with default paths
#'
#' @param norm_method    Chr. Normalisation method (folder/file name, no .tsv).
#' @param pval_threshold Num. FDR/BH threshold.
#' @param fc_threshold   Num. |log2FC| threshold.
test_analysis_de <- function(
  norm_method = "VSN",
  pval_threshold = DEG_PVAL_THRESHOLD,
  fc_threshold = DE_FC_THRESHOLD
) {
    main_analysis_de(
        counts_path    = "DATA/Read_counts/",
        norm_base      = "Analyses/NormTest",
        norm_method    = norm_method,
        output_path    = "Analyses",
        pval_threshold = pval_threshold,
        fc_threshold   = fc_threshold
    )
}

if (interactive()) {
    test_analysis_de()
}
