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

library(ggstar)


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
    df <- inner_join(
        tpm_data, protQ_data, by = c("Species", "Protein_id", "Treatment")
        )

    stats_df <- df %>%
        group_by(Species, Treatment) %>%
        summarise(
            R = cor(
                mean_log2_TPM,
                .data[[mean_protQ_col]],
                method = "pearson",
                use    = "complete.obs"
            ),
            p_value = tryCatch(
                cor.test(
                    mean_log2_TPM,
                    .data[[mean_protQ_col]],
                    method = "pearson"
                )$p.value,
                error = function(e) NA_real_
            ),
            .groups = "drop"
        ) %>%
        mutate(
            label = sprintf("R: %.2f\np-value: %.2e", R, p_value)
        )
    left_join(df, stats_df, by = c("Species", "Treatment"))
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
#' @param output_path Character. Directory path for saving plots
#'
#' @return NULL (invisibly). Saves PDF and PNG files to output_path
plot_correlation <- function(
    correlation_df, mean_protQ_col, protQ_name, output_path
    ) {
    correlation_df %>%
        group_by(Species) %>%
        group_walk(function(sp_df, sp_key) {
            species <- sp_key$Species
            message("Plotting correlation for: ", species)

            # Compute point density per treatment
            sp_df <- sp_df %>%
                nest(data = -Treatment) %>%
                mutate(data = map(data, function(df) {
                    df %>% mutate(
                        density = compute_point_density(mean_log2_TPM, !!sym(mean_protQ_col))
                    )
                })) %>%
                unnest(data)

            # Correlation labels per treatment
            label_df <- sp_df %>%
                group_by(Treatment) %>%
                summarise(
                    x     = min(mean_log2_TPM, na.rm = TRUE),
                    y     = max(!!sym(mean_protQ_col), na.rm = TRUE),
                    label = dplyr::first(label),
                    .groups = "drop"
                )

            # Order treatments: Control first, then alphabetical
            treatments         <- unique(sp_df$Treatment)
            ordered_treatments <- c("Ctrl", sort(treatments[treatments != "Ctrl"]))
            sp_df    <- sp_df    %>% mutate(Treatment = factor(Treatment, levels = ordered_treatments))
            label_df <- label_df %>% mutate(Treatment = factor(Treatment, levels = ordered_treatments))

            p <- ggplot(sp_df, aes(x = mean_log2_TPM, y = !!sym(mean_protQ_col))) +
                geom_point(aes(colour = density), alpha = 0.6, size = 1.2) +
                scale_colour_gradientn(
                    colours = c(
                        "#c8d8e8", "#a8d8ea", "#8B7BB5", "#6B3A6B",
                        "#8B2020", "#C0392B", "#e8c8cc", "#c8d8e8"
                    ),
                    name = "counts"
                ) +
                geom_smooth(method = "lm", se = FALSE, colour = "red", linewidth = 0.8) +
                geom_text(
                    data = label_df,
                    aes(x = x, y = y, label = label),
                    hjust = 0, vjust = 1, size = 2.8, family = "mono"
                ) +
                facet_wrap(~Treatment, scales = "free", ncol = 6) +
                labs(
                    x     = "Mean TPM (log2)",
                    y     = paste("Mean", protQ_name, "(log2)"),
                    title = bquote(italic(.(gsub("_", " ", species))))
                ) +
                theme_bw(base_size = 10) +
                theme(
                    plot.title       = element_text(face = "italic"),
                    strip.background = element_blank(),
                    strip.text       = element_text(hjust = 0, face = "bold"),
                    panel.grid.minor = element_blank(),
                    legend.position  = "right"
                )

            filename <- file.path(
                output_path, 
                species,
                paste("TPM", protQ_name, "correlation_panels", sep = "_")
            )

            print(p)
            
            ggsave(paste0(filename, ".pdf"), plot = p, width = 24, height = 6, units = "in", dpi = 300)
            message("Saved: ", paste0(filename, ".pdf"))
            ggsave(paste0(filename, ".png"), plot = p, width = 24, height = 6, units = "in", dpi = 300)
            message("Saved: ", paste0(filename, ".png"))
        })

    invisible(NULL)
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
                    species     = sp$species_name,
                    measurement = meas_name
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
        mutate(label = sprintf("R: %.2f\np: %.3f", meta_R, meta_pval))
}

#' Build Label Position for Top-Right Annotation
#'
#' @description
#' Computes x/y coordinates for placing the correlation label at the
#' top-right corner of the plot, in the expanded zone above the data.
#'
#' @param df Tibble. Filtered to a single measurement. Must contain
#'   columns: exposure_label (factor), R (dbl).
#' @param stats Tibble. Single-row output of compute_exposure_correlations()
#'   filtered to the same measurement. Must contain: label (chr).
#'
#' @return Single-row tibble with columns: x (chr), y (dbl), label (chr)
build_label_position <- function(df, stats) {
    y_max <- max(df$R, na.rm = TRUE)
    y_range <- diff(range(df$R, na.rm = TRUE))

    tibble(
        x = last(levels(df$exposure_label)),
        y = y_max + y_range * 0.15
    ) %>%
        bind_cols(stats %>% dplyr::select(label))
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
        geom_star(aes(starshape = Treatment), size = 3, alpha = 0.8) +
        scale_starshape_manual(values = TREATMENT_SHAPES) +
        geom_label(
            data = label_pos,
            aes(x = x, y = y, label = label),
            hjust = 1,
            vjust = 1,
            size = 3,
            family = "mono",
            color = "black",
            fill = "white",
            label.color = "black",
            label.size = 0.4,
            label.padding = unit(0.4, "lines"),
            inherit.aes = FALSE
        ) +
        scale_y_continuous(
            expand = expansion(mult = c(0.05, 0.25))
        ) +
        labs(
            title = title,
            x     = "Exposure time",
            y     = "Pearson R",
            fill  = "Species",
            shape = "Treatment"
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
