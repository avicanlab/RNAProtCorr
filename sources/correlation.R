# mRNA-Protein Correlation Analysis -------------------------------------------
#
# Description:
#   Functions for computing and visualising Pearson correlations between
#   transcript abundance (TPM) and protein quantification across bacterial
#   species and experimental conditions.
#
# Requirements:
#   library(ggplot2)
#   library(ggstar)
#   library(dplyr)
#   library(purrr)
#   library(tidyr)

library(ggstar)

# ── Correlation analysis ───────────────────────────────────────────────────────

#' Calculate Pearson Correlations Between TPM and Protein Abundance
#'
#' @description
#' Joins TPM and protein data, then computes per-species and per-treatment
#' Pearson correlation coefficients and p-values. Results are joined back onto
#' the input data for downstream use in scatter plots.
#'
#' @param tpm_data       Tibble. Log2-transformed TPM data with a
#'   `mean_log2_TPM` column.
#' @param protQ_data     Tibble. Protein quantification data with treatment
#'   groups.
#' @param mean_protQ_col Character. Column name for mean protein quantification
#'   values.
#'
#' @return Tibble with the original joined data plus per-treatment columns:
#'   `R` (Pearson estimate), `p_value`, and `label`.
calculate_correlation <- function(tpm_data, protQ_data, mean_protQ_col) {
  df <- dplyr::inner_join(
    tpm_data, protQ_data,
    by = c("Species", "Protein_id", "Treatment")
  )

  stats_df <- df |>
    dplyr::group_by(Species, Treatment) |>
    dplyr::summarise(
      cor_test = list(tryCatch(
        stats::cor.test(mean_log2_TPM, .data[[mean_protQ_col]], method = "pearson"),
        error = function(e) NULL
      )),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      R       = purrr::map_dbl(cor_test, ~ if (is.null(.x)) NA_real_ else .x$estimate),
      p_value = purrr::map_dbl(cor_test, ~ if (is.null(.x)) NA_real_ else .x$p.value),
      label   = sprintf("R: %.2f\np: %.2e\n%s", R, p_value, Treatment)
    ) |>
    dplyr::select(-cor_test)

  dplyr::left_join(df, stats_df, by = c("Species", "Treatment"))
}

# ── Visualisation ──────────────────────────────────────────────────────────────

#' Compute 2-D Bin Density for Points
#'
#' @description
#' Assigns each point a count equal to the number of other points sharing its
#' 2-D bin. This is a fast, dependency-free alternative to kernel density
#' estimation for colouring dense scatter plots.
#'
#' @param x      Numeric vector. x-coordinates.
#' @param y      Numeric vector. y-coordinates.
#' @param n_bins Integer. Number of bins per axis. Default 50.
#'
#' @return Numeric vector of bin counts (same length as `x`); NA for
#'   non-finite input values.
compute_point_density <- function(x, y, n_bins = 50) {
  keep    <- is.finite(x) & is.finite(y)
  x_clean <- x[keep]
  y_clean <- y[keep]

  x_breaks <- seq(min(x_clean), max(x_clean), length.out = n_bins + 1L)
  y_breaks <- seq(min(y_clean), max(y_clean), length.out = n_bins + 1L)

  ix <- pmin(.bincode(x_clean, x_breaks, include.lowest = TRUE), n_bins)
  iy <- pmin(.bincode(y_clean, y_breaks, include.lowest = TRUE), n_bins)

  bin_key    <- paste(ix, iy, sep = "_")
  bin_counts <- table(bin_key)

  density_out          <- rep(NA_real_, length(x))
  density_out[keep]    <- as.integer(bin_counts[bin_key])
  density_out
}

#' Plot mRNA-Protein Correlation Panels per Treatment and Species
#'
#' @description
#' Creates density-coloured scatter plots faceted by treatment, one figure per
#' species. A linear regression line and per-treatment correlation statistics
#' are annotated. Plots are saved to disk and returned as a list.
#'
#' @param correlation_df  Tibble. Output of `calculate_correlation()`. Must
#'   contain `mean_log2_TPM`, `Treatment`, `Species`, `label`, and the column
#'   named by `mean_protQ_col`.
#' @param mean_protQ_col  Character. Column name for mean protein quantification.
#' @param protQ_name      Character. Display name for the quantification measure
#'   (e.g. `"iBAQ"`).
#' @param output_path     Character. Directory for saving plot files.
#'
#' @return List of ggplot objects (one per species).
plot_correlation <- function(
  correlation_df, mean_protQ_col, protQ_name, output_path
) {
  plots <- correlation_df |>
    dplyr::group_by(Species) |>
    dplyr::group_map(function(sp_df, sp_key) {
      species  <- sp_key$Species
      message("Plotting correlation for: ", species)
      title_sp <- format_species_title(species)

      # Compute per-treatment point density
      sp_df <- sp_df |>
        tidyr::nest(data = -Treatment) |>
        dplyr::mutate(data = purrr::map(data, function(df) {
          dplyr::mutate(
            df,
            density = compute_point_density(mean_log2_TPM, .data[[mean_protQ_col]])
          )
        })) |>
        tidyr::unnest(data)

      # Annotation position per treatment
      label_df <- sp_df |>
        dplyr::group_by(Treatment) |>
        dplyr::summarise(
          x     = min(mean_log2_TPM, na.rm = TRUE),
          y     = max(.data[[mean_protQ_col]], na.rm = TRUE),
          label = dplyr::first(label),
          .groups = "drop"
        )

      # Ordered treatments: Control first, then alphabetical
      treatments         <- unique(sp_df$Treatment)
      ordered_treatments <- c("Ctrl", sort(treatments[treatments != "Ctrl"]))
      sp_df    <- dplyr::mutate(sp_df,    Treatment = factor(Treatment, levels = ordered_treatments))
      label_df <- dplyr::mutate(label_df, Treatment = factor(Treatment, levels = ordered_treatments))

      p <- ggplot2::ggplot(sp_df, ggplot2::aes(x = mean_log2_TPM, y = .data[[mean_protQ_col]])) +
        ggplot2::geom_point(ggplot2::aes(colour = density), alpha = 0.5, size = 0.6) +
        ggplot2::scale_colour_gradientn(
          colours = c(
            "#c8d8e8", "#a8d8ea", "#8B7BB5", "#6B3A6B",
            "#8B2020", "#C0392B", "#e8c8cc", "#c8d8e8"
          ),
          name = "counts"
        ) +
        ggplot2::geom_smooth(method = "lm", se = FALSE, colour = "red", linewidth = 0.8) +
        ggplot2::geom_text(
          data         = label_df,
          ggplot2::aes(x = x, y = y, label = label),
          lineheight   = 0.5,
          hjust        = 0, vjust = 1,
          size         = 5,
          fontface     = "bold",
          family       = "Arial"
        ) +
        ggplot2::facet_wrap(~Treatment, scales = "free_x", ncol = 6) +
        ggplot2::labs(
          x     = "Mean TPM (log2)",
          y     = paste("Mean", protQ_name, "(log2)"),
          title = title_sp
        ) +
        theme_publication(base_size = 24, legend_position = "right") +
        ggplot2::theme(strip.text = ggplot2::element_blank())

      filename <- file.path(
        output_path, species,
        paste("TPM", protQ_name, "correlation_panels", sep = "_")
      )
      save_plot(p, filename, width = 12, height = 6, units = "in", dpi = 300,
                out_format = FIGURE_FORMAT)
      message("mRNA-protein correlation for ", species, " saved: ", filename, MSG_FIG_FORMAT)

      # Return a clean version without the small annotation text for compositing
      p$layers[[3]] <- NULL
      p +
        ggplot2::geom_text(
          data       = label_df,
          ggplot2::aes(x = x, y = y, label = label),
          lineheight = 0.5,
          hjust      = 0, vjust = 1,
          size       = 8, fontface = "bold", family = "Arial"
        )
    })

  # Combined: one species per row, shared x-axis
  patchwork::wrap_plots(plots, ncol = 1) +
    patchwork::plot_layout(guides = "keep", axes = "collect_x") &
    ggplot2::theme(
      axis.title    = ggplot2::element_text(size = 40),
      axis.text     = ggplot2::element_text(size = 36),
      strip.text    = ggplot2::element_blank(),
      legend.position = "right"
    )
}

# ── Correlation vs exposure time ───────────────────────────────────────────────

#' Compute Spearman Meta-Correlation Between Exposure Time and R
#'
#' @description
#' Computes a global Spearman correlation between numeric exposure time and
#' mRNA-protein R values across all species and treatments.
#'
#' @param df Tibble. Must contain columns: `exposure_time` (dbl), `R` (dbl).
#'
#' @return Single-row tibble with columns: `meta_R`, `meta_pval`, `label`.
compute_exposure_correlations <- function(df) {
  df |>
    dplyr::summarise(
      meta_R    = stats::cor(exposure_time, R, method = "spearman", use = "complete.obs"),
      meta_pval = tryCatch(
        stats::cor.test(exposure_time, R, method = "spearman")$p.value,
        error = function(e) NA_real_
      )
    ) |>
    dplyr::mutate(label = sprintf("R: %.2f - p: %.3f", meta_R, meta_pval))
}

#' Build Label Position for Top-Right Annotation
#'
#' @description
#' Computes the x/y coordinates for placing the meta-correlation label at the
#' top-right corner of an exposure-time scatter plot.
#'
#' @param df    Tibble. Filtered to one measurement. Must contain columns:
#'   `exposure_label` (factor), `R` (dbl).
#' @param stats Single-row tibble from `compute_exposure_correlations()`.
#'   Must contain: `label` (chr).
#'
#' @return Single-row tibble with columns: `x` (chr), `y` (dbl), `label` (chr).
build_label_position <- function(df, stats) {
  y_max   <- max(df$R, na.rm = TRUE)
  y_range <- diff(range(df$R, na.rm = TRUE))

  tibble::tibble(
    x = dplyr::last(levels(df$exposure_label)),
    y = y_max + y_range * 0.05
  ) |>
    dplyr::bind_cols(dplyr::select(stats, label))
}

#' Build a Single Exposure-Time Scatter Plot
#'
#' @description
#' Creates a scatter plot of Pearson R vs exposure time for one measurement,
#' coloured by species and shaped by treatment (via `ggstar`). A top-right
#' annotation shows the global meta-correlation.
#'
#' @param df        Tibble. Data for one measurement. Output of
#'   `build_exposure_time_df()` filtered to a single `measurement` value.
#' @param label_pos Single-row tibble from `build_label_position()`.
#'
#' @return ggplot object.
build_exposure_plot <- function(df, label_pos) {
  ggplot2::ggplot(df, ggplot2::aes(x = exposure_label, y = R, fill = Species)) +
    ggstar::geom_star(ggplot2::aes(starshape = Treatment), size = 3, alpha = 0.8) +
    ggstar::scale_starshape_manual(values = TREATMENT_SHAPES) +
    ggplot2::geom_label(
      data          = label_pos,
      ggplot2::aes(x = x, y = y, label = label),
      hjust         = 4,
      vjust         = 1,
      linewidth     = 0.25,
      family        = "Arial",
      color         = "black",
      fill          = "white",
      label.size    = 0.4,
      label.padding = ggplot2::unit(0.4, "lines"),
      inherit.aes   = FALSE
    ) +
    ggplot2::scale_fill_discrete(labels = PLOT_SPECIES_NAMES) +
    ggplot2::labs(
      x         = "Exposure time",
      y         = "mRNA-protein Pearson correlation",
      fill      = "Species",
      starshape = "Treatment"
    ) +
    theme_publication(legend_position = "right")
}
