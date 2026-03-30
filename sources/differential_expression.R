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
plot_de_barplot_all <- function(deg_df, dep_df, output_path, prefix = NULL) {
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
                theme_publication()
        })

    # Save per species
    walk2(plots, unique(deg_df$Species), function(p, species) {
        output_file <- file.path(
            output_path, species, paste0(prefix, "_", "DEG_DEP_barplot")
        )
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
  deg_df, dep_df, corr_df, protq_name, output_path, prefix = NULL
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
            " - p = ", format.pval(pearson_test$p.value,
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
                size = 8,
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
            theme_publication()
        out_file <- file.path(
            output_path,
            paste(prefix, de_type, "vs_correlation", protq_name, sep = "_")
        )
        save_plot(
            plot = p,
            filepath = out_file,
            width = 16, height = 5
        )
        message("  Plot saved: ", out_file, MSG_FIG_FORMAT)
        p
    }

    return(
        list(
            DEG = make_scatter(build_plot_df(deg_df), "DEGs"),
            DEP = make_scatter(build_plot_df(dep_df), "DEPs")
        )
    )
}