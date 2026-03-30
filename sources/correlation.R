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
    tpm_data, protQ_data,
    by = c("Species", "Protein_id", "Treatment")
  )

  stats_df <- df %>%
    group_by(Species, Treatment) %>%
    summarise(
      cor_test = list(tryCatch(
        cor.test(mean_log2_TPM, .data[[mean_protQ_col]], method = "pearson"),
        error = function(e) NULL
      )),
      .groups = "drop"
    ) %>%
    mutate(
      R = map_dbl(cor_test, ~if (is.null(.x)) NA_real_ else .x$estimate),
      p_value = map_dbl(cor_test, ~if (is.null(.x)) NA_real_ else .x$p.value),
      label = sprintf("R: %.2f\np: %.2e", R, p_value)
    ) %>%
    dplyr::select(-cor_test)
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
  plots <- correlation_df %>%
    group_by(Species) %>%
    group_map(function(sp_df, sp_key) {
      species <- sp_key$Species
      message("Plotting correlation for: ", species)

      display <- PLOT_SPECIES_NAMES[species]
      display <- ifelse(is.na(display), gsub("_", " ", species), display)
      parts <- strsplit(display, " ")[[1]]
      title_sp <- format_species_title(species)
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
          x = min(mean_log2_TPM, na.rm = TRUE),
          y = max(!!sym(mean_protQ_col), na.rm = TRUE),
          label = dplyr::first(label),
          .groups = "drop"
        )

      # Order treatments: Control first, then alphabetical
      treatments <- unique(sp_df$Treatment)
      ordered_treatments <- c("Ctrl", sort(treatments[treatments != "Ctrl"]))
      sp_df <- sp_df %>% mutate(Treatment = factor(Treatment, levels = ordered_treatments))
      label_df <- label_df %>% mutate(Treatment = factor(Treatment, levels = ordered_treatments))

      p <- ggplot(sp_df, aes(x = mean_log2_TPM, y = !!sym(mean_protQ_col))) +
        geom_point(aes(colour = density), alpha = 0.5, size = 0.6) +
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
          aes(x = x, y = y, label = label), lineheight = 0.5,
          hjust = 0, vjust = 1, size = 5, fontface = "bold", family = "Arial"
        ) +
        facet_wrap(~Treatment, scales = "free_x", ncol = 6) +
        labs(
          x = "Mean TPM (log2)",
          y = paste("Mean", protQ_name, "(log2)"),
          title = title_sp
        ) +
        theme_publication(base_size = 24, legend_position = "right")

      filename <- file.path(
        output_path,
        species,
        paste("TPM", protQ_name, "correlation_panels", sep = "_")
      )

      save_plot(p, filename, width = 12, height = 6, units = "in", dpi = 300, out_format = FIGURE_FORMAT)
      message("mRNA-protein correlatio for", species, "saved: ", filename, MSG_FIG_FORMAT)

      p$layers[[3]] <- NULL
      p +
        geom_text(
          data = label_df,
          aes(x = x, y = y, label = label), lineheight = 0.5,
          hjust = 0, vjust = 1, size = 8, fontface = "bold", family = "Arial"
        )
    })

  # Combined: one species per row, shared x axis
  wrap_plots(plots, ncol = 1) +
    plot_layout(guides = "keep", axes = "collect_x") &
    theme(
      axis.title = element_text(size = 40),
      axis.text = element_text(size = 36),
      strip.text = element_text(size = 40),
      legend.position = "right"
    )
}

# ============================================================================
# CORRELATION VS EXPOSURE TIME ANALYSIS
# ============================================================================

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
    summarise(
      meta_R = cor(exposure_time, R, method = "spearman", use = "complete.obs"),
      meta_pval = tryCatch(
        cor.test(exposure_time, R, method = "spearman")$p.value,
        error = function(e) NA_real_
      )
    ) %>%
    mutate(label = sprintf("R: %.2f - p: %.3f", meta_R, meta_pval))
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
    y = y_max + y_range * 0.05
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
#'
#' @return ggplot object
build_exposure_plot <- function(df, label_pos) {
  ggplot(df, aes(x = exposure_label, y = R, fill = Species)) +
    geom_star(aes(starshape = Treatment), size = 3, alpha = 0.8) +
    scale_starshape_manual(values = TREATMENT_SHAPES) +
    geom_label(
      data = label_pos,
      aes(x = x, y = y, label = label),
      hjust = 4,
      vjust = 1,
      linewidth = 0.25,
      family = "Arial",
      color = "black",
      fill = "white",
      label.size = 0.4,
      label.padding = unit(0.4, "lines"),
      inherit.aes = FALSE
    ) +
    # scale_y_continuous(
    #     expand = expansion(mult = c(0.05, 0.25))
    # ) +
    scale_fill_discrete(labels = PLOT_SPECIES_NAMES) +
    labs(
      x = "Exposure time",
      y = "mRNA-protein Pearson correlation",
      fill = "Species",
      starshape = "Treatment" # was: shape = "Treatment"
    ) +
    theme_publication(legend_position = "right",)
}
