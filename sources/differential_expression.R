# Differential Expression Analysis — DEG (TPM) and DEP (Proteomics) -----------
#
# Description:
#   Functions for running limma-based differential expression analysis on both
#   transcriptomic (DEG) and proteomic (DEP) data, and for visualising results.
#
# Requirements:
#   library(limma)
#   library(dplyr)
#   library(tidyr)
#   library(stringr)
#   library(ggplot2)
#   library(purrr)

library(limma)

# ── Configuration ──────────────────────────────────────────────────────────────

#' Stress treatment display order (reversed alphabetical for y-axis readability)
TREATMENT_ORDER <- rev(sort(STRESS_TREATMENTS_PROT))

#' Fill colours for DEG vs DEP bar charts
DE_COLOURS <- c("DEGs" = "#89CFF0", "DEPs" = "#F5A623")

# ── limma helpers ──────────────────────────────────────────────────────────────

#' Build a Wide Expression Matrix from Long-Format Data
#'
#' @description
#' Pivots a long-format tibble to a wide numeric matrix suitable for limma,
#' with one column per sample (`Replicate_Treatment`) and one row per protein.
#' Rows with any missing value are removed.
#'
#' @param data       Tibble with columns: `Protein_id`, `Replicate`,
#'   `Treatment`, `value`.
#' @param treatments Character vector. Treatments to include; NULL keeps all.
#'   Default NULL.
#'
#' @return Named list:
#'   - `mat`      : numeric matrix (proteins × samples).
#'   - `metadata` : tibble with columns `sample`, `Replicate`, `Treatment`.
build_expression_matrix <- function(data, treatments = NULL) {
  if (!is.null(treatments)) data <- dplyr::filter(data, Treatment %in% treatments)

  wide <- data |>
    dplyr::mutate(sample = paste(Replicate, Treatment, sep = "_")) |>
    dplyr::select(Protein_id, sample, value) |>
    dplyr::group_by(Protein_id, sample) |>
    dplyr::summarise(value = mean(value, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = sample, values_from = value) |>
    tibble::column_to_rownames("Protein_id")

  wide <- wide[stats::complete.cases(wide), ]
  mat  <- as.matrix(wide)

  metadata <- tibble::tibble(sample = colnames(mat)) |>
    dplyr::mutate(
      Replicate = stringr::str_extract(sample, "^R\\d+"),
      Treatment = stringr::str_extract(sample, "[A-Za-z]+$")
    )

  list(mat = mat, metadata = metadata)
}

#' Run limma DEP Analysis for One Species
#'
#' @description
#' Fits a linear model with treatment and replicate as factors using limma,
#' applies empirical Bayes moderation, and extracts results for each stress
#' treatment vs control.
#'
#' @param mat            Numeric matrix. Rows = proteins, cols = samples.
#' @param metadata       Tibble with columns: `sample`, `Replicate`, `Treatment`.
#' @param species        Character. Species name (added to output).
#' @param pval_threshold Numeric. BH adj. p-value threshold. Default
#'   `DEP_PVAL_THRESHOLD`.
#' @param fc_threshold   Numeric. |log2FC| threshold. Default `DE_FC_THRESHOLD`.
#'
#' @return Tibble with columns: `Species`, `Treatment`, `n_de`, `n_total`,
#'   `pct_de`, `de_ids` — or NULL when no appropriate samples exist.
run_limma_dep <- function(
  mat,
  metadata,
  species,
  pval_threshold = DEP_PVAL_THRESHOLD,
  fc_threshold   = DE_FC_THRESHOLD
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

  metadata <- metadata |>
    dplyr::mutate(
      Treatment = factor(Treatment, levels = c("Ctrl", stress_present)),
      Replicate = factor(Replicate)
    )

  design       <- stats::model.matrix(~0 + Treatment + Replicate, data = metadata)
  colnames(design) <- stringr::str_replace(colnames(design), "^Treatment", "")

  fit          <- limma::lmFit(mat, design)
  contrast_mat <- limma::makeContrasts(
    contrasts = paste0(stress_present, " - Ctrl"),
    levels    = design
  )
  colnames(contrast_mat) <- stress_present
  fit2 <- limma::eBayes(limma::contrasts.fit(fit, contrast_mat))

  purrr::map_dfr(stress_present, function(trt) {
    tt <- limma::topTable(
      fit2,
      coef         = trt,
      number       = Inf,
      sort.by      = "none",
      adjust.method = "BH"
    )
    sig_mask <- !is.na(tt$adj.P.Val) &
      tt$adj.P.Val < pval_threshold &
      abs(tt$logFC) > fc_threshold

    tibble::tibble(
      Species   = species,
      Treatment = trt,
      n_de      = sum(sig_mask),
      n_total   = nrow(tt),
      pct_de    = 100 * sum(sig_mask) / nrow(tt),
      de_ids    = list(rownames(tt)[sig_mask])
    )
  })
}

# ── Visualisation ──────────────────────────────────────────────────────────────

#' Plot DEG vs DEP Percentage Bar Chart — All Species Combined
#'
#' @description
#' Horizontal grouped bar chart showing the percentage of differentially
#' expressed genes (DEGs) and proteins (DEPs) per stress condition, with one
#' panel per species. Plots are saved per species and returned as a list.
#'
#' @param deg_df      Tibble. DEG summary with columns: `Species`, `Treatment`,
#'   `pct_de`. Output of the DEG calling step.
#' @param dep_df      Tibble. DEP summary with the same structure. Output of
#'   `run_limma_dep()`.
#' @param output_path Character. Output directory.
#' @param prefix      Character or NULL. Optional filename prefix.
#'
#' @return List of ggplot objects (one per species).
plot_de_barplot_all <- function(deg_df, dep_df, output_path, prefix = NULL) {
  plots <- dplyr::bind_rows(
    dplyr::mutate(deg_df, type = "DEGs"),
    dplyr::mutate(dep_df, type = "DEPs")
  ) |>
    dplyr::filter(Treatment %in% STRESS_TREATMENTS_PROT) |>
    dplyr::mutate(
      Treatment = factor(Treatment, levels = TREATMENT_ORDER),
      type      = factor(type, levels = c("DEGs", "DEPs"))
    ) |>
    dplyr::group_by(Species) |>
    dplyr::group_map(function(sp_df, sp_key) {
      species <- sp_key$Species
      x_max   <- max(ceiling(max(sp_df$pct_de, na.rm = TRUE) / 10) * 10, 20)

      ggplot2::ggplot(sp_df, ggplot2::aes(x = pct_de, y = Treatment, fill = type)) +
        ggplot2::geom_col(
          position = ggplot2::position_dodge(width = 0.7, reverse = TRUE),
          width    = 0.6,
          colour   = NA
        ) +
        ggplot2::scale_fill_manual(values = DE_COLOURS, name = NULL) +
        ggplot2::scale_x_continuous(
          limits = c(0, x_max),
          expand = ggplot2::expansion(mult = c(0, 0.02)),
          breaks = seq(0, x_max, by = 20)
        ) +
        ggplot2::scale_y_discrete(drop = FALSE) +
        ggplot2::labs(
          title = format_species_title(species),
          x     = "Percentage",
          y     = "Stress Conditions"
        ) +
        theme_publication()
    })

  purrr::walk2(plots, unique(deg_df$Species), function(p, species) {
    output_file <- file.path(
      output_path,
      species,
      paste0(prefix, "_DEG_DEP_barplot")
    )
    save_plot(p, output_file, width = 16, height = 5)
    message("  Plot saved: ", output_file, MSG_FIG_FORMAT)
  })

  plots
}

#' Plot Percentage DE vs mRNA-Protein Correlation — All Species
#'
#' @description
#' Scatter plots (one for DEGs, one for DEPs) showing the relationship between
#' the percentage of differentially expressed genes or proteins and the
#' mRNA-protein Pearson R per condition and species. A linear trend is fitted.
#'
#' @param deg_df      Tibble. DEG summary (Species, Treatment, pct_de).
#' @param dep_df      Tibble. DEP summary (Species, Treatment, pct_de).
#' @param corr_df     Tibble. Correlation results (Species, Treatment, R).
#'   Output of `calculate_correlation()`.
#' @param protq_name  Character. Proteomics method label (e.g. `"iBAQ"`).
#' @param output_path Character. Output directory.
#' @param prefix      Character or NULL. Optional filename prefix.
#'
#' @return Named list: `DEG` and `DEP` ggplot objects.
plot_de_vs_correlation_all <- function(
  deg_df, dep_df, corr_df, protq_name, output_path, prefix = NULL
) {
  corr_summary <- dplyr::distinct(corr_df, Species, Treatment, R)

  # Helper: join DE data with correlation summary
  build_plot_df <- function(de_df) {
    de_df |>
      dplyr::select(Species, Treatment, pct_de) |>
      dplyr::inner_join(corr_summary, by = c("Species", "Treatment")) |>
      dplyr::filter(is.finite(pct_de), is.finite(R))
  }

  # Helper: build one scatter for a given DE type
  make_scatter <- function(plot_df, de_type) {
    if (nrow(plot_df) < 3) {
      warning("Too few points for ", de_type, " vs correlation — skipping.")
      return(invisible(NULL))
    }

    pearson_test  <- stats::cor.test(plot_df$pct_de, plot_df$R, method = "pearson")
    ann_label     <- paste0(
      sprintf("R = %.2f", pearson_test$estimate),
      " - p = ",
      format.pval(pearson_test$p.value, digits = 2, eps = .Machine$double.xmin)
    )
    species_levels <- sort(unique(plot_df$Species))

    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = R, y = pct_de, colour = Treatment)) +
      ggplot2::geom_smooth(
        ggplot2::aes(group = 1),
        method    = "lm",
        formula   = y ~ x,
        se        = TRUE,
        colour    = "grey50",
        fill      = "grey85",
        linewidth = 0.7
      ) +
      ggplot2::geom_point(ggplot2::aes(shape = Species), size = 3) +
      ggplot2::annotate(
        "label",
        x             = max(plot_df$R),
        y             = max(plot_df$pct_de),
        label         = ann_label,
        hjust         = 1, vjust = 1, size = 8,
        fill          = "white",
        label.padding = ggplot2::unit(0.2, "lines"),
        fontface      = "italic"
      ) +
      ggplot2::scale_shape_manual(
        values = setNames(seq_along(species_levels) + 14, species_levels),
        labels = setNames(lapply(species_levels, format_species_title), species_levels)
      ) +
      ggplot2::labs(
        x      = bquote("mRNA-protein correlation" ~ italic(R) ~ "(" * .(protq_name) * ")"),
        y      = paste0("% ", de_type),
        colour = "Treatment",
        shape  = "Species"
      ) +
      theme_publication()

    out_file <- file.path(
      output_path,
      paste(prefix, de_type, "vs_correlation", protq_name, sep = "_")
    )
    save_plot(p, out_file, width = 16, height = 5)
    message("  Plot saved: ", out_file, MSG_FIG_FORMAT)
    p
  }

  list(
    DEG = make_scatter(build_plot_df(deg_df), "DEGs"),
    DEP = make_scatter(build_plot_df(dep_df), "DEPs")
  )
}
