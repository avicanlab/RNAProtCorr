# ============================================================================
# Differential Expression Analysis — DEG (TPM) and DEP (Proteomics)
# ============================================================================

library(limma)

# ============================================================================
# CONFIGURATION
# ============================================================================

TREATMENT_ORDER <- rev(sort(STRESS_TREATMENTS_PROT))
DE_COLOURS <- c("DEGs" = "#89CFF0", "DEPs" = "#F5A623")

# ============================================================================
# LIMMA (DEP)
# ============================================================================

#' Build Wide Expression Matrix from Long-Format Data
#'
#' @param data       Tibble: Protein_id, Replicate, Treatment, value.
#' @param treatments Chr vector. Treatments to include (NULL = all).
#' @return List: mat (matrix), metadata (tibble)
build_expression_matrix <- function(data, treatments = NULL) {
    if (!is.null(treatments)) data <- data %>% filter(Treatment %in% treatments)

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

#' Run limma DEP Analysis for One Species
#'
#' @param mat            Numeric matrix. Rows = proteins, cols = samples.
#' @param metadata       Tibble. Columns: sample, Replicate, Treatment.
#' @param species        Chr. Species name (added to output).
#' @param pval_threshold Num. BH adj. p-value threshold.
#' @param fc_threshold   Num. |log2FC| threshold.
#' @return Tibble: Species, Treatment, n_de, n_total, pct_de, de_ids
run_limma_dep <- function(
  mat,
  metadata,
  species,
  pval_threshold = DEP_PVAL_THRESHOLD,
  fc_threshold = DE_FC_THRESHOLD
) {
    treatments_present <- unique(metadata$Treatment)

    if (!"Ctrl" %in% treatments_present) {
        warning("No Ctrl samples found for ", species, " — skipping DEP.")
        return(NULL)
    }

    stress_present <- intersect(STRESS_TREATMENTS_PROT, treatments_present)
    if (length(stress_present) == 0) {
        warning("No stress treatments found for ", species, " — skipping DEP.")
        return(NULL)
    }

    metadata <- metadata %>%
        mutate(
            Treatment = factor(Treatment, levels = c("Ctrl", stress_present)),
            Replicate = factor(Replicate)
        )

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
            Species   = species,
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

#' Plot DEG vs DEP Percentage Bar Chart — All Species Combined
#'
#' @param deg_df      Tibble. call_deg() output with Species column.
#' @param dep_df      Tibble. call_dep_all() output with Species column.
#' @param output_path Chr. Output directory.
#' @return NULL invisibly.
plot_de_barplot_all <- function(deg_df, dep_df, output_path) {
    plots <- bind_rows(
        deg_df %>% mutate(type = "DEGs"),
        dep_df %>% mutate(type = "DEPs")
    ) %>%
        filter(Treatment %in% STRESS_TREATMENTS_PROT) %>%
        mutate(
            Treatment = factor(Treatment, levels = TREATMENT_ORDER),
            type      = factor(type, levels = c("DEGs", "DEPs"))
        ) %>%
        group_by(Species) %>%
        group_map(function(sp_df, sp_key) {
            species <- sp_key$Species

            x_max <- ceiling(max(sp_df$pct_de, na.rm = TRUE) / 10) * 10
            x_max <- max(x_max, 20)

            ggplot(sp_df, aes(x = pct_de, y = Treatment, fill = type)) +
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
                    title = format_species_title(species),
                    x     = "Percentage",
                    y     = "Stress Conditions"
                ) +
                theme_bw(base_size = 11) +
                theme(
                    plot.title         = element_text(hjust = 0.5, size = 11),
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
        })

    # Save per species
    walk2(plots, unique(deg_df$Species), function(p, species) {
        output_file <- file.path(output_path, species, "DEG_DEP_barplot")
        save_plot(
            plot = p,
            filepath = output_file,
            width = 16, height = 5
        )
        message("  Plot saved: ", output_file, MSG_FIG_FORMAT)
    })

    plots
}

#' Plot % DE vs mRNA-Protein Correlation — All Species
#'
#' @param deg_df      Tibble. call_deg() output.
#' @param dep_df      Tibble. call_dep_all() output.
#' @param corr_df     Tibble. calculate_correlation() output with Species column.
#' @param protq_name  Chr. Proteomics method label.
#' @param output_path Chr. Output directory.
#' @return NULL invisibly.
plot_de_vs_correlation_all <- function(
  deg_df, dep_df, corr_df, protq_name, output_path
) {
    corr_summary <- corr_df %>% distinct(Species, Treatment, R)

    build_plot_df <- function(de_df) {
        de_df %>%
            dplyr::select(Species, Treatment, pct_de) %>%
            inner_join(corr_summary, by = c("Species", "Treatment")) %>%
            filter(is.finite(pct_de), is.finite(R))
    }

    make_scatter <- function(plot_df, de_type) {
        if (nrow(plot_df) < 3) {
            warning("Too few points for ", de_type, " vs correlation — skipping.")
            return(invisible(NULL))
        }

        pearson_test <- cor.test(plot_df$pct_de, plot_df$R, method = "pearson")
        label <- paste0(
            sprintf("R = %.2f", pearson_test$estimate),
            "\np = ", format.pval(pearson_test$p.value,
                digits = 2,
                eps = .Machine$double.xmin
            )
        )

        species_levels <- sort(unique(plot_df$Species))

        p <- ggplot(plot_df, aes(
            x      = R,
            y      = pct_de,
            colour = Treatment
        )) +
            geom_smooth(
                aes(group = 1),
                method    = "lm",
                formula   = y ~ x,
                se        = TRUE,
                colour    = "grey50",
                fill      = "grey85",
                linewidth = 0.7
            ) +
            geom_point(aes(shape = Species), size = 3) +
            annotate(
                "label",
                x = max(plot_df$R),
                y = max(plot_df$pct_de),
                label = label,
                hjust = 1, vjust = 1,
                size = 3,
                fill = "white",
                label.padding = unit(0.2, "lines"),
                fontface = "italic"
            ) +
            scale_shape_manual(
                values = setNames(seq_along(species_levels) + 14, species_levels),
                labels = setNames(
                    lapply(species_levels, format_species_title),
                    species_levels
                )
            ) +
            labs(
                x = bquote("mRNA-protein correlation" ~ italic(R) ~
                    "(" * .(protq_name) * ")"),
                y = paste0("% ", de_type),
                colour = "Treatment",
                shape = "Species"
            ) +
            theme_bw(base_size = 11) +
            theme(
                axis.title       = element_text(size = 10),
                axis.text        = element_text(size = 9),
                panel.grid.minor = element_blank(),
                legend.text      = element_text(size = 8),
                legend.title     = element_text(size = 9),
                legend.key.size  = unit(0.45, "cm")
            )

        save_plot(
            plot = p,
            filepath = file.path(
                output_path,
                paste0(de_type, "_vs_correlation_", protq_name)
            ),
            width = 16, height = 5
        )
        message("  Plot saved: ", de_type, "_vs_correlation_", protq_name, MSG_FIG_FORMAT)
        p
    }

    return(
        list(
            DEG = make_scatter(build_plot_df(deg_df), "DEGs"),
            DEP = make_scatter(build_plot_df(dep_df), "DEPs")
        )
    )
}

# ============================================================================
# OUTPUT HELPERS
# ============================================================================

write_de_table <- function(de_df, output_path, filename) {
    de_df %>%
        dplyr::select(Species, Treatment, n_de, n_total, pct_de) %>%
        write.table(
            file.path(output_path, filename),
            sep = "\t", row.names = FALSE, quote = FALSE
        )
}

# ============================================================================
# MAIN
# ============================================================================

#' Run Full DEG / DEP Pipeline Across All Species
#'
#' @param tpm_data       Tibble. All-species long-format TPM data with Species,
#'   Protein_id, Treatment, TPM columns.
#' @param protq_data     Tibble. All-species long-format proteomics data with
#'   Species, Protein_id, Treatment, and quantification column.
#' @param protq_name     Chr. Name of the quantification column.
#' @param counts_path    Chr. Directory containing read counts Excel files.
#' @param norm_base      Chr. Base path to NormalyzerDE output directory.
#' @param norm_method    Chr. Normalization method name (e.g. "VSN").
#' @param output_path    Chr. Root output directory.
#' @param common_ids_only Logical. Restrict to IDs present in both DEG and
#'   DEP datasets. Default TRUE.
#' @param pval_threshold Num. FDR/BH threshold. Default DEG_PVAL_THRESHOLD.
#' @param fc_threshold   Num. |log2FC| threshold. Default DE_FC_THRESHOLD.
#'
#' @return List with deg_results, dep_results, corr_df (invisibly).
main_analysis_de <- function(
  tpm_data,
  protq_data,
  protq_name,
  counts_path,
  norm_base,
  norm_method,
  output_path,
  common_ids_only = TRUE,
  pval_threshold = DEG_PVAL_THRESHOLD,
  fc_threshold = DE_FC_THRESHOLD
) {
    output_de <- file.path(output_path, "DE")
    dir.create(output_de, recursive = TRUE, showWarnings = FALSE)

    # ---- Load DEG data for all species --------------------------------------
    message("Loading DEG data...")
    deg_data_all <- read_deg_data_all(counts_path)

    # ---- Load normalised proteomics for all species -------------------------
    message("Loading normalised proteomics...")
    species_list <- unique(tpm_data$Species)

    norm_data_all <- map_dfr(species_list, function(species) {
        norm_path <- file.path(
            norm_base, species, protq_name,
            paste0(norm_method, "-normalized.txt")
        )
        if (!file.exists(norm_path)) {
            warning("Normalised file not found for ", species, ": ", norm_path)
            return(NULL)
        }
        read_norm_proteomics(norm_path) %>% mutate(Species = species)
    })

    # ---- Optionally restrict to common IDs ----------------------------------
    if (common_ids_only) {
        common_ids <- intersect(
            unique(deg_data_all$Protein_id),
            unique(norm_data_all$Protein_id)
        )
        message("Common IDs (DEG ∩ DEP): ", length(common_ids))
        deg_data_all <- deg_data_all %>% filter(Protein_id %in% common_ids)
        norm_data_all <- norm_data_all %>% filter(Protein_id %in% common_ids)
        tpm_data <- tpm_data %>% filter(Protein_id %in% common_ids)
        protq_data <- protq_data %>% filter(Protein_id %in% common_ids)
    }

    # ---- Call DEGs and DEPs -------------------------------------------------
    message("Calling DEGs...")
    deg_results <- call_deg(deg_data_all, pval_threshold, fc_threshold)
    write_de_table(deg_results, output_de, "AllSpecies_DEG_results.tsv")

    message("Calling DEPs via limma...")
    dep_results <- call_dep_all(norm_data_all, pval_threshold, fc_threshold)
    write_de_table(dep_results, output_de, "AllSpecies_DEP_results.tsv")

    # ---- Compute mRNA-protein correlation -----------------------------------
    message("Computing mRNA-protein correlation...")
    tpm_log2 <- tpm_data %>%
        log2_transform(value_col = "TPM", pseudocount = TRUE) %>%
        dplyr::rename(mean_log2_TPM = mean_log2)

    mean_protq_col <- paste0("mean_log2_", protq_name)
    corr_df <- calculate_correlation(tpm_log2, protq_data, mean_protq_col)

    # ---- Plots --------------------------------------------------------------
    message("Generating plots...")
    plot_de_barplot_all(deg_results, dep_results, output_de)

    plot_de_vs_correlation_all(
        deg_df       = deg_results,
        dep_df       = dep_results,
        corr_df      = corr_df,
        protq_name   = protq_name,
        output_path  = output_de
    )

    message("\nDE analysis complete. Results written to: ", output_de)
    invisible(list(
        deg_results = deg_results,
        dep_results = dep_results,
        corr_df     = corr_df
    ))
}
