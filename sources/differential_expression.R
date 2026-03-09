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
library(tidyverse)
library(readxl)

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

#' Plot % DE vs mRNA-Protein Correlation Scatter
#'
#' @description
#' Scatter plot with one point per treatment showing the relationship between
#' the percentage of differentially expressed features (x-axis) and the
#' per-treatment Pearson R from mRNA-protein level correlation (y-axis).
#' Pearson R and Student's t-test p-value are computed across treatments and
#' displayed as a label inside the panel.
#' Call this function once for DEGs and once for DEPs.
#'
#' @param de_df       Tibble. call_deg() or run_limma_dep() output.
#'   Must contain: Treatment, pct_de.
#' @param corr_df     Tibble. calculate_correlation() output from correlation.R.
#'   Must contain: Treatment, R (per-treatment Pearson R mRNA vs protein).
#' @param de_type     Chr. Label for the DE type, e.g. "DEGs" or "DEPs".
#' @param protq_name  Chr. Proteomics quantification method, e.g. "iBAQ".
#'   Used in axis label and filename.
#' @param species     Chr. Species name for italic title.
#' @param output_path Chr. Output directory.
#'
#' @return NULL invisibly.
plot_de_vs_correlation <- function(
  de_df,
  corr_df,
  de_type,
  protq_name,
  species,
  output_path
) {
    # Join % DE with per-treatment mRNA-protein R
    plot_df <- de_df %>%
        dplyr::select(Treatment, pct_de) %>%
        inner_join(
            corr_df %>% distinct(Treatment, R),
            by = "Treatment"
        ) %>%
        filter(is.finite(pct_de), is.finite(R))

    if (nrow(plot_df) < 3) {
        warning(
            "  Too few points for DE vs correlation scatter (",
            de_type, ", ", protq_name, ") — skipping."
        )
        return(invisible(NULL))
    }

    # Pearson correlation + Student's t-test p-value across treatments
    pearson_test <- cor.test(plot_df$pct_de, plot_df$R, method = "pearson")
    r_val <- pearson_test$estimate
    p_val <- pearson_test$p.value
    label <- sprintf("R = %.2f\np = %.2e", r_val, p_val)

    # Label position: top-right inside panel
    label_x <- max(plot_df$pct_de) + diff(range(plot_df$pct_de)) * 0.05
    label_y <- max(plot_df$R)

    p <- ggplot(plot_df, aes(x = pct_de, y = R, label = Treatment)) +
        geom_smooth(
            method    = "lm",
            formula   = y ~ x,
            se        = TRUE,
            colour    = "grey50",
            fill      = "grey85",
            linewidth = 0.7
        ) +
        geom_point(
            size   = 3,
            colour = if (de_type == "DEGs") DE_COLOURS[["DEGs"]] else DE_COLOURS[["DEPs"]]
        ) +
        geom_text(
            vjust  = -0.8,
            size   = 3,
            colour = "grey30"
        ) +
        annotate(
            "label",
            x          = label_x,
            y          = label_y,
            label      = label,
            hjust      = 1,
            vjust      = 1,
            size       = 3,
            fill       = "white",
            label.size = 0.3,
            fontface   = "italic"
        ) +
        labs(
            title = bquote(italic(.(gsub("_", " ", species)))),
            x = paste0("% ", de_type),
            y = bquote("mRNA-protein correlation" ~ italic(R) ~
                "(" * .(protq_name) * ")")
        ) +
        theme_bw(base_size = 11) +
        theme(
            plot.title = element_text(face = "italic", hjust = 0.5, size = 11),
            axis.title = element_text(size = 10),
            axis.text = element_text(size = 9),
            panel.grid.minor = element_blank()
        )

    save_plot(
        plot = p,
        filepath = file.path(
            output_path,
            paste0(species, "_", de_type, "_vs_correlation_", protq_name)
        ),
        width = 5,
        height = 5
    )

    message(
        "  Plot saved: ", species, "_", de_type,
        "_vs_correlation_", protq_name, ".pdf/.png"
    )
    invisible(NULL)
}

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
#' Calls DEGs from pre-computed columns, runs limma for DEPs, writes summary
#' TSVs, plots the bar chart, and generates scatter plots of % DE vs
#' mRNA-protein correlation (once for DEGs, once for DEPs).
#'
#' @param deg_data       Tibble. Output of read_deg_data().
#' @param norm_path      Chr. Path to normalised proteomics TSV.
#' @param corr_df        Tibble. calculate_correlation() output for this species
#'   and proteomics method. Must contain Treatment and R columns.
#'   Pass NULL to skip the DE vs correlation scatter plots.
#' @param protq_name     Chr. Proteomics quantification method label (e.g. "iBAQ").
#'   Used in scatter plot axis label and filename.
#' @param species        Chr. Species name.
#' @param output_path    Chr. Output directory.
#' @param pval_threshold Num. FDR/BH threshold. Default DEG_PVAL_THRESHOLD.
#' @param fc_threshold   Num. |log2FC| threshold. Default DE_FC_THRESHOLD.
#' @param common_ids_only Logical. If TRUE restricts to IDs present in both
#'   DEG and DEP datasets. Default TRUE.
#'
#' @return List: species, deg_results, dep_results.
process_species_de <- function(
  deg_data,
  tpm_data,
  norm_path,
  protq_path,
  corr_df = NULL,
  protq_name = NULL,
  species,
  output_path,
  pval_threshold = DEG_PVAL_THRESHOLD,
  fc_threshold = DE_FC_THRESHOLD,
  common_ids_only = TRUE
) {
    message("\nProcessing DE for: ", species)

    # ---- DEP from normalised proteomics via limma ---------------------------
    dep_results <- NULL
    prot_long <- NULL

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

    # ---- Optionally restrict to common IDs ----------------------------------
    if (common_ids_only && !is.null(prot_long)) {
        common_ids <- intersect(
            deg_data %>% pull(Protein_id) %>% unique(),
            prot_long %>% pull(Protein_id) %>% unique()
        )
        message("  Common IDs (DEG ∩ DEP): ", length(common_ids))
        deg_data <- deg_data %>% filter(Protein_id %in% common_ids)
        prot_long <- prot_long %>% filter(Protein_id %in% common_ids)
    }

    # ---- DEG calling --------------------------------------------------------
    deg_results <- call_deg(deg_data, pval_threshold, fc_threshold)

    write_de_table(deg_results, output_path, paste0(species, "_DEG_results.tsv"))
    message(
        "  DEG: ", sum(deg_results$n_de), " total DE genes across ",
        nrow(deg_results), " conditions"
    )

    # ---- DEP via limma ------------------------------------------------------
    if (!is.null(prot_long)) {
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
    }

    # ---- Bar chart ----------------------------------------------------------
    if (!is.null(dep_results)) {
        plot_de_barplot(deg_results, dep_results, species, output_path)
    } else {
        message("  Skipping bar chart — DEP results unavailable.")
    }

    # ---- DE vs mRNA-protein correlation scatter plots -----------------------
    ## ---- Optionally restrict to common IDs ---------------------------------
    protq_data <- read_tsv(protq_path,
        col_names = TRUE,
        show_col_types = FALSE, skip_empty_rows = TRUE
    ) %>%
        filter(Species == species)

    if (common_ids_only && !is.null(protq_data)) {
        common_ids <- intersect(
            tpm_data %>% pull(Protein_id) %>% unique(),
            protq_data %>% pull(Protein_id) %>% unique()
        )
        message("  Common IDs (DEG ∩ DEP): ", length(common_ids))
        tpm_data <- tpm_data %>% filter(Protein_id %in% common_ids)
        protq_data <- protq_data %>% filter(Protein_id %in% common_ids)
    }

    tpm_data <- tpm_data %>%
        log2_transform(value_col = "TPM", pseudocount = TRUE) %>%
        rename(mean_Log2_TPM = mean_log2)

    corr_df <- calculate_correlation(
        tpm_data, protq_data,
        paste0(ifelse(protq_name == "iBAQ_MC", "iBAQ", protq_name), "_meanlog2")
    )

    plot_de_vs_correlation(
        de_df       = deg_results,
        corr_df     = corr_df,
        de_type     = "DEGs",
        protq_name  = protq_name,
        species     = species,
        output_path = output_path
    )
    plot_de_vs_correlation(
        de_df       = dep_results,
        corr_df     = corr_df,
        de_type     = "DEPs",
        protq_name  = protq_name,
        species     = species,
        output_path = output_path
    )


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
#' @param correlation_results List. Output of main_analysis() from correlation.R.
#'   Structured as: correlation_results[[species]][[corr_<method>]].
#'   Used to supply per-treatment R values for the DE vs correlation scatter.
#'   Pass NULL to skip scatter plots.
#'
#' @return Named list of per-species results (invisibly).
main_analysis_de <- function(
  counts_path,
  norm_base,
  norm_method,
  protQ_dir,
  output_path,
  pval_threshold = DEG_PVAL_THRESHOLD,
  fc_threshold = DE_FC_THRESHOLD
) {
    counts_files <- list.files(counts_path, pattern = "\\.xlsx$", full.names = TRUE)

    if (length(counts_files) == 0) {
        stop("No read counts Excel files found in: ", counts_path)
    }

    output_de <- file.path(output_path, "DE")
    dir.create(output_de, showWarnings = FALSE)

    results <- list()

    for (counts_file in counts_files) {
        species <- tools::file_path_sans_ext(
            str_extract(basename(counts_file), "^[A-Za-z]+_[A-Za-z]+")
        )
        message("Processing species: ", species)

        # Look up species abbreviation
        abbv <- species_abbv_map[species]
        if (is.na(abbv)) {
            warning(
                "No abbreviation found for species '", species,
                "', skipping file: ", basename(tpm_file)
            )
            return(invisible(NULL))
        }

        # Load raw TPM (before averaging — kept for downstream use)
        tpm_data <- read_tpm(counts_file, abbv) %>%
            mutate(
                Species    = species,
                Treatment  = ifelse(Treatment == "Mig", "Hyp", Treatment),
                Protein_id = sub("^pi", "PI", Protein_id, ignore.case = FALSE)
            ) %>%
            filter(!is.na(Protein_id), !is.na(TPM))

        message("  Loading read counts: ", basename(counts_file))
        deg_data <- read_deg_data(counts_file)

        for (protq in PROTQ_METHODS) {
            message("Processing ", protq, " for ", species, "...")

            norm_path <- file.path(
                norm_base, species, protq,
                paste0(norm_method, "-normalized.txt")
            )
            protq_path <- file.path(
                protQ_dir, paste0(protq, "_data.tsv")
            )
            output_sp_protq <- file.path(output_de, species, protq, norm_method)
            dir.create(output_sp_protq, recursive = TRUE, showWarnings = FALSE)

            result <- tryCatch(
                process_species_de(
                    deg_data        = deg_data,
                    tpm_data        = tpm_data,
                    norm_path       = norm_path,
                    protq_path      = protq_path,
                    corr_df         = corr_df,
                    protq_name      = protq,
                    species         = species,
                    output_path     = output_sp_protq,
                    pval_threshold  = pval_threshold,
                    fc_threshold    = fc_threshold
                ),
                error = function(e) {
                    warning("Failed for '", species, "' / '", protq, "': ", e$message)
                    NULL
                }
            )

            if (!is.null(result)) results[[species]][[protq]] <- result
        }
    }

    message("\nDE analysis complete. Results written to: ", output_de)
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
  correlation_results = NULL,
  pval_threshold = DEG_PVAL_THRESHOLD,
  fc_threshold = DE_FC_THRESHOLD
) {
    # correlation_results can be obtained by running:
    #   corr <- main_analysis(TPM_PATH, Intensity_PATH, RI_PATH,
    #                         IBAQ_PATH, IBAQ_MC_PATH, output_path = "Analyses")
    # then passing corr here to enable DE vs correlation scatter plots.
    main_analysis_de(
        counts_path         = "DATA/Read_counts/",
        norm_base           = "Analyses/NormTest",
        norm_method         = norm_method,
        protQ_dir           = "Analyses",
        output_path         = "Analyses",
        pval_threshold      = pval_threshold,
        fc_threshold        = fc_threshold
    )
}

if (interactive()) {
    test_analysis_de()
}
