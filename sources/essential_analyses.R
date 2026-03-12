# ============================================================================
# Essential genes analyses
# ============================================================================

#' Species Abbreviation Mapping
#' Maps full species names to abbreviations for data parsing
species_essential2tpm_map <- c(
    "Senterica"           = "Salmonella_enterica",
    "Saureus"             = "Staphylococcus_aureus",
    "Ypseudotuberculosis" = "Yersinia_pseudotuberculosis"
)

#' Plot Essential vs Non-Essential Gene Expression Distribution
#'
#' @description
#' Overlapping frequency histogram comparing global Log2 mean expression
#' between essential and non-essential genes, with a two-sided Wilcoxon
#' p-value matching the original analysis behaviour.
#'
#' @param all_df Tibble. All genes with a global Log2 mean column
#'   (output of log2_global_mean()).
#' @param essential_df Tibble. IDs of essential genes per species.
#' @param log2_col Chr. Name of the Log2 mean column in all_df.
#' @param x_lab Chr. Short label for the data type (e.g. "TPM", "iBAQ").
#' @param output_path Chr. Directory for saving output files.
#'
#' @return NULL invisibly.
plot_essential_distribution <- function(
  all_df, essential_df, log2_col, x_lab, output_path
) {
    condition_df <- all_df %>%
        left_join(
            essential_df %>% mutate(group = "Essential"),
            by = c("Species", "Protein_id")
        ) %>%
        mutate(
            group = factor(
                replace_na(group, "Non-Essential"),
                levels = c("Non-Essential", "Essential")
            )
        ) %>%
        filter(is.finite(.data[[log2_col]]))

    x_label <- bquote(
        .(ifelse(x_lab == "TPM", "mRNA", "Protein")) ~
            "- Mean" ~ .(x_lab) ~ (Log[2])
    )

    plots <- condition_df %>%
        group_by(Species) %>%
        group_map(function(sp_df, sp_key) {
            species <- sp_key$Species

            p_value <- wilcox.test(
                as.formula(paste(log2_col, "~ group")),
                data = sp_df
            )$p.value

            p <- ggplot(
                sp_df,
                aes(x = .data[[log2_col]], y = after_stat(density), fill = group)
            ) +
                geom_histogram(position = "identity", alpha = 0.5, bins = 40) +
                scale_fill_manual(
                    values = c("Non-Essential" = "#F4A8A0", "Essential" = "#80CEC8")
                ) +
                annotate(
                    "text",
                    x = Inf, y = Inf,
                    label = format.pval(p_value, digits = 2, eps = .Machine$double.xmin),
                    hjust = 1.1, vjust = 1.5, size = 3.5, fontface = "italic"
                ) +
                labs(
                    title = bquote(italic(.(gsub("_", " ", species)))),
                    x     = x_label,
                    y     = "Frequency",
                    fill  = NULL
                ) +
                theme_minimal() +
                theme(
                    plot.title      = element_text(face = "italic"),
                    legend.position = c(0.85, 0.85)
                )

            save_plot(p, file.path(
                output_path, species, paste0(x_lab, "_essential_distribution")
            ))
            p
        })

    # Combined figure — all species side by side
    wrap_plots(plots, nrow = 1) +
        plot_layout(guides = "collect") &
        theme(legend.position = "top")
}

plot_essential_vs_all_correlation <- function(
  corr_all_df, corr_essential_df, measurement, output_path
) {
    # Join both correlation sets by treatment
    plot_df <- inner_join(
        corr_all_df %>% distinct(Species, Treatment, R) %>% dplyr::rename(R_all = R),
        corr_essential_df %>% distinct(Species, Treatment, R) %>% dplyr::rename(R_essential = R),
        by = c("Species", "Treatment")
    )

    plots <- plot_df %>%
        group_by(Species) %>%
        group_map(function(sp_df, sp_key) {
            species <- sp_key$Species

            p <- ggplot(sp_df, aes(x = R_all, y = R_essential, color = Treatment)) +
                geom_abline(
                    slope     = 1,
                    intercept = 0,
                    linetype  = "dashed",
                    color     = "grey60",
                    linewidth = 0.6
                ) +
                geom_point(size = 4, alpha = 0.9) +
                scale_color_manual(values = setNames(
                    scales::hue_pal()(nrow(sp_df)),
                    sp_df$Treatment
                )) +
                coord_fixed(
                    ratio = 1,
                    xlim  = c(0.4, 1),
                    ylim  = c(0.4, 1)
                ) +
                labs(
                    title = bquote(italic(.(gsub("_", " ", species)))),
                    x     = "All genes'\nmRNA-protein level correlation",
                    y     = "Essential genes'\nmRNA-protein level correlation",
                    color = "Treatment"
                ) +
                theme_bw() +
                theme(
                    plot.title = element_text(face = "italic", hjust = 0.5),
                    axis.title = element_text(size = 10),
                    legend.position = "right"
                )

            save_plot(
                plot     = p,
                filepath = file.path(output_path, species, paste0(measurement, "_essential_vs_all_correlation")),
                width    = 5,
                height   = 5
            )
            p
        })

    wrap_plots(plots, nrow = 1) +
        plot_layout(guides = "collect") &
        theme(legend.position = "right")
}

#' Plot Essential vs Non-Essential SD Distribution
#'
#' @description
#' Frequency polygon comparing log2 expression standard deviation between
#' essential and non-essential genes, with a two-sided Wilcoxon p-value.
#'
#' @param all_df Tibble. All genes with SD column (output of log2_global_sd()).
#' @param essential_vec Chr vector. IDs of essential genes.
#' @param sd_col Chr. Name of the SD column in all_df.
#' @param x_lab Chr. Short label for the data type (e.g. "TPM", "iBAQ").
#' @param species Chr. Species name for plot title.
#' @param output_path Chr. Directory for saving output files.
#'
#' @return NULL invisibly.
plot_sd_distribution <- function(
  all_df, essential_df, sd_col, x_lab, output_path
) {
    condition_df <- all_df %>%
        left_join(
            essential_df %>% mutate(group = "Essential"),
            by = c("Species", "Protein_id")
        ) %>%
        mutate(
            group = factor(
                replace_na(group, "Non-Essential"),
                levels = c("Non-Essential", "Essential")
            )
        ) %>%
        filter(is.finite(.data[[sd_col]]))

    x_label <- if (x_lab == "TPM") {
        bquote("Standard deviation of mRNA - " ~ log[2](mean ~ TPM))
    } else {
        bquote("Standard deviation of protein - " ~ log[2](mean ~ .(x_lab)))
    }

    plots <- condition_df %>%
        group_by(Species) %>%
        group_map(function(sp_df, sp_key) {
            species <- sp_key$Species
            p_value <- wilcox.test(
                as.formula(paste(sd_col, "~ group")),
                data = sp_df
            )$p.value
            p <- ggplot(
                sp_df,
                aes(x = .data[[sd_col]], color = group)
            ) +
                geom_step(
                    aes(y = after_stat(density)),
                    stat = "bin",
                    binwidth = 0.05,
                    linewidth = 0.8,
                    direction = "mid" # centers the step on each bin
                ) +
                scale_color_manual(
                    values = c("Non-Essential" = "#F4A8A0", "Essential" = "#80CEC8")
                ) +
                annotate(
                    "text",
                    x = Inf, y = Inf,
                    label = sprintf("p-value: %.2e", p_value),
                    hjust = 1.1, vjust = 1.5,
                    size = 3.5,
                    fontface = "italic"
                ) +
                labs(
                    title = bquote(italic(.(gsub("_", " ", species)))),
                    x     = x_label,
                    y     = "Frequency",
                    color = NULL
                ) +
                theme_bw() +
                theme(
                    plot.title        = element_text(hjust = 0.5),
                    legend.position   = c(0.75, 0.85),
                    legend.background = element_rect(fill = "white", color = NA)
                )

            save_plot(
                plot     = p,
                filepath = file.path(output_path, species, paste0(x_lab, "_sd_distribution")),
                width    = 5,
                height   = 5
            )
            p
        })

    # Combined figure — all species side by side
    wrap_plots(plots, nrow = 1) +
        plot_layout(guides = "collect") &
        theme(legend.position = "top")
}
