# Discordant mRNA-Protein regulation patterns ---------------------------------
#
# Description:
#   Identifies proteins with discordant mRNA-protein regulation patterns
#   across stress conditions using either logFC (limma) or Z-score approaches.
#   Also provides enrichment analysis and visualisation functions.
#
# Inputs:
#   deg_raw    : Raw DEG results (Species, Protein_id, Treatment, log2FC, fdr)
#   dep_raw    : Raw DEP results (Species, Protein_id, Treatment, logFC, adj.P.Val)
#   tpm_data   : TPM expression data (for Z-score approach)
#   protQ_data : Protein quantification data (for Z-score approach)
#
# Requirements:
#   library(limma)
#   library(clusterProfiler)

library(clusterProfiler)

# ── Cluster constants ──────────────────────────────────────────────────────────

CLUSTER_COLOURS <- c(
  "Cluster 1" = "#E8A838",
  "Cluster 2" = "#E87070",
  "Cluster 3" = "#C47EC4",
  "Cluster 4" = "#4BBFBF",
  "Cluster 5" = "#4BAF6E",
  "Cluster 6" = "#E8A070"
)

#' Human-readable labels for each discordance cluster
#' Format: "mRNA direction – protein direction"
CLUSTER_LABELS <- c(
  "Cluster 1" = "down \u2013 no change",
  "Cluster 2" = "down \u2013 up",
  "Cluster 3" = "no change \u2013 down",
  "Cluster 4" = "no change \u2013 up",
  "Cluster 5" = "up \u2013 down",
  "Cluster 6" = "up \u2013 no change"
)

# ── Shared helpers ─────────────────────────────────────────────────────────────

#' Tag Regulation Direction From a Numeric Value and Thresholds
#'
#' @param value          Numeric vector.
#' @param threshold      Numeric. Absolute value threshold.
#' @param pval           Numeric vector of p-values. NULL uses the Z-score
#'   approach (threshold only). Default NULL.
#' @param pval_threshold Numeric. P-value threshold; ignored when `pval` is
#'   NULL. Default NULL.
#'
#' @return Character vector: `"down"`, `"up"`, or `"unchanged"`.
tag_direction <- function (value, threshold, pval = NULL, pval_threshold = NULL) {
  if (!is.null(pval)) {
    dplyr::case_when(
      pval < pval_threshold & value < -threshold ~ "down",
      pval < pval_threshold & value > threshold ~ "up",
      TRUE ~ "unchanged"
    )
  } else {
    dplyr::case_when(
      value < -threshold ~ "down",
      value > threshold ~ "up",
      TRUE ~ "unchanged"
    )
  }
}

#' Assign Discordance Clusters to Per-Protein Regulation Patterns
#'
#' @description
#' Summarises per-treatment regulation flags into per-protein boolean
#' patterns, then maps those patterns onto the six pre-defined discordance
#' clusters. Proteins that do not match any cluster are dropped.
#'
#' Cluster definitions (mRNA – protein):
#'   1: down only    – unchanged in all
#'   2: down only    – up in ≥1 condition
#'   3: unchanged    – down in ≥1 condition
#'   4: unchanged    – up in ≥1 condition
#'   5: up only      – down in ≥1 condition
#'   6: up only      – unchanged in all
#'
#' @param joined Tibble. Must have columns: `Species`, `Protein_id`,
#'   `Treatment`, `rna_sig`, `prot_sig`.
#'
#' @return Tibble with columns: `Species`, `Protein_id`, `Treatment`, `cluster`.
assign_clusters <- function (joined) {
  joined |>
    dplyr::group_by(Species, Protein_id, Treatment) |>
    dplyr::summarise(
      rna_ever_down = any(rna_sig == "down"),
      rna_ever_up = any(rna_sig == "up"),
      rna_ever_unchanged = any(rna_sig == "unchanged"),
      prot_ever_down = any(prot_sig == "down"),
      prot_ever_up = any(prot_sig == "up"),
      prot_ever_unchanged = any(prot_sig == "unchanged"),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      cluster = dplyr::case_when(
        rna_ever_down &
          !rna_ever_up &
          !rna_ever_unchanged &
          !prot_ever_down &
          !prot_ever_up &
          prot_ever_unchanged ~ "Cluster 1",
        rna_ever_down &
          !rna_ever_up &
          !rna_ever_unchanged &
          !prot_ever_down &
          prot_ever_up &
          !prot_ever_unchanged ~ "Cluster 2",
        !rna_ever_down &
          !rna_ever_up &
          rna_ever_unchanged &
          prot_ever_down &
          !prot_ever_up &
          !prot_ever_unchanged ~ "Cluster 3",
        !rna_ever_down &
          !rna_ever_up &
          rna_ever_unchanged &
          !prot_ever_down &
          prot_ever_up &
          !prot_ever_unchanged ~ "Cluster 4",
        !rna_ever_down &
          rna_ever_up &
          !rna_ever_unchanged &
          prot_ever_down &
          !prot_ever_up &
          !prot_ever_unchanged ~ "Cluster 5",
        !rna_ever_down &
          rna_ever_up &
          !rna_ever_unchanged &
          !prot_ever_down &
          !prot_ever_up &
          prot_ever_unchanged ~ "Cluster 6",
        TRUE ~ NA_character_
      )
    ) |>
    dplyr::filter(!is.na(cluster)) |>
    dplyr::select(Species, Protein_id, Treatment, cluster)
}

#' Colour Cluster Strip Backgrounds in a gtable
#'
#' @description
#' Iterates over top-strip grobs in a `ggplot_gtable()` output and sets their
#' fill colours to match `CLUSTER_COLOURS`.
#'
#' @param g A gtable object from `ggplot_gtable(ggplot_build(p))`.
#'
#' @return The modified gtable, invisibly.
colour_cluster_strips <- function (g) {
  strip_idx <- which(grepl("strip-t", g$layout$name))
  cluster_names <- names(CLUSTER_COLOURS)
  for (i in seq_along(strip_idx)) {
    col <- CLUSTER_COLOURS[cluster_names[((i - 1L) %% length(cluster_names)) + 1L]]
    g$grobs[[strip_idx[i]]]$grobs[[1]]$children[[1]]$
      gp$
      fill <- col
  }
  g
}

# ── logFC approach ─────────────────────────────────────────────────────────────

#' Build Discordance Clusters from logFC (limma DEG / DEP)
#'
#' @description
#' Joins RNA and protein logFC per protein × treatment, tags regulation
#' direction, and assigns each protein to one of six discordance clusters
#' based on its pattern across all stress conditions.
#'
#' @param deg_raw        Tibble. Columns: `Species`, `Protein_id`, `Treatment`,
#'   `log2FC`, `fdr`.
#' @param dep_raw        Tibble. Columns: `Species`, `Protein_id`, `Treatment`,
#'   `logFC`, `adj.P.Val`.
#' @param fc_threshold   Numeric. |log2FC| threshold. Default `DE_FC_THRESHOLD`.
#' @param pval_threshold Numeric. Adjusted p-value threshold. Default
#'   `DEP_PVAL_THRESHOLD`.
#'
#' @return Tibble: `Species`, `Protein_id`, `Treatment`, `rna_log2FC`,
#'   `prot_log2FC`, `rna_sig`, `prot_sig`, `cluster`.
build_discordance_clusters <- function (
  deg_raw,
  dep_raw,
  fc_threshold = DE_FC_THRESHOLD,
  pval_threshold = DEP_PVAL_THRESHOLD
) {
  rna_df <- deg_raw |>
    dplyr::filter(!is.na(log2FC), !is.na(fdr)) |>
    dplyr::mutate(rna_sig = tag_direction(log2FC, fc_threshold, fdr, pval_threshold)) |>
    dplyr::select(Species, Protein_id, Treatment, rna_log2FC = log2FC, rna_sig)

  prot_df <- dep_raw |>
    dplyr::mutate(prot_sig = tag_direction(logFC, fc_threshold, adj.P.Val, pval_threshold)) |>
    dplyr::select(Species, Protein_id, Treatment, prot_log2FC = logFC, prot_sig)

  joined <- dplyr::inner_join(rna_df, prot_df, by = c("Species", "Protein_id", "Treatment"))
  protein_pattern <- assign_clusters(joined)

  joined |>
    dplyr::inner_join(protein_pattern, by = c("Species", "Protein_id")) |>
    dplyr::select(Species, Protein_id, Treatment,
                  rna_log2FC, prot_log2FC, rna_sig, prot_sig, cluster) |>
    dplyr::arrange(Species, cluster, Protein_id, Treatment)
}

#' Extract Raw Per-Protein DEP Results from limma
#'
#' @description
#' Fits a limma model and returns a tidy tibble of per-gene results for every
#' stress treatment vs control. Unlike `run_limma_dep()`, this function returns
#' raw statistics (logFC, adj.P.Val) for every protein without applying
#' significance thresholds.
#'
#' @param mat            Numeric matrix. Rows = proteins, cols = samples.
#' @param metadata       Tibble with columns: `sample`, `Replicate`, `Treatment`.
#' @param species        Character. Species name added to the output.
#'
#' @return Tibble: `Species`, `Protein_id`, `Treatment`, `logFC`, `adj.P.Val`.
extract_dep_raw <- function (mat, metadata, species) {
  stress_present <- intersect(STRESS_TREATMENTS_PROT, unique(metadata$Treatment))

  metadata <- metadata |>
    dplyr::mutate(
      Treatment = factor(Treatment, levels = c("Ctrl", stress_present)),
      Replicate = factor(Replicate)
    )

  design <- stats::model.matrix(~0 + Treatment + Replicate, data = metadata)
  colnames(design) <- stringr::str_replace(colnames(design), "^Treatment", "")
  contrast_mat <- limma::makeContrasts(
    contrasts = paste0(stress_present, " - Ctrl"),
    levels = design
  )
  colnames(contrast_mat) <- stress_present
  fit2 <- limma::eBayes(limma::contrasts.fit(limma::lmFit(mat, design), contrast_mat))

  purrr::map_dfr(stress_present, function (trt) {
    limma::topTable(fit2, coef = trt, number = Inf,
                    sort.by = "none", adjust.method = "BH") |>
      tibble::rownames_to_column("Protein_id") |>
      dplyr::mutate(Species = species, Treatment = trt) |>
      dplyr::select(Species, Protein_id, Treatment, logFC, adj.P.Val)
  })
}

# ── Z-score approach ───────────────────────────────────────────────────────────

#' Z-Score Transform Expression Data Relative to Control
#'
#' @description
#' Computes per-protein Z-scores using control mean and SD:
#'   `z = (X − mean_ctrl) / sd_ctrl`
#' Proteins with non-finite control mean or any infinite stress values are
#' excluded. Control samples are excluded from the output.
#'
#' @param df           Tibble. Must contain `Species`, `Protein_id`,
#'   `Treatment`, and the value column.
#' @param col_name     Character. Name of the value column to transform.
#' @param fudge_factor Numeric. Added to `sd_ctrl` for numerical stability.
#'   Default 0.
#'
#' @return Tibble: `Species`, `Protein_id`, `Treatment`, `value`, `mean_ctrl`,
#'   `sd_ctrl`, `zscore`.
zscore_transform <- function (df, col_name, fudge_factor = 0) {
  mean_sd_df <- df |>
    dplyr::filter(is.finite(.data[[col_name]])) |>
    dplyr::group_by(Species) |>
    dplyr::mutate(
      mean_ctrl = mean(.data[[col_name]][Treatment == "Ctrl"], na.rm = TRUE),
      sd_ctrl = stats::sd(.data[[col_name]][Treatment == "Ctrl"], na.rm = TRUE)
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(Treatment != "Ctrl")

  # Remove proteins with any infinite stress values
  infinite_ids <- mean_sd_df |>
    dplyr::filter(is.infinite(.data[[col_name]])) |>
    dplyr::distinct(Species, Protein_id)

  finite_df <- mean_sd_df |> dplyr::anti_join(infinite_ids, by = c("Species", "Protein_id"))
  global_median_sd <- stats::median(finite_df$sd_ctrl, na.rm = TRUE) + fudge_factor

  finite_df |>
    dplyr::mutate(zscore = (.data[[col_name]] - mean_ctrl) / sd_ctrl) |>
    dplyr::select(Species, Protein_id, Treatment,
                  value = dplyr::all_of(col_name), mean_ctrl, sd_ctrl, zscore)
}

#' Build Discordance Clusters from Z-Scores
#'
#' @description
#' Same cluster definitions as `build_discordance_clusters()`, but uses
#' Z-score thresholds rather than logFC + p-value.
#'
#' @param deg_zscore  Tibble. Output of `zscore_transform()` on TPM data.
#'   Must contain: `Species`, `Protein_id`, `Treatment`, `zscore`.
#' @param dep_zscore  Tibble. Output of `zscore_transform()` on protein data.
#' @param z_threshold Numeric. |Z-score| threshold. Default 1.96.
#'
#' @return Tibble: `Species`, `Protein_id`, `Treatment`, `rna_z`, `prot_z`,
#'   `rna_sig`, `prot_sig`, `cluster`.
build_discordance_clusters_zscore <- function (
  deg_zscore,
  dep_zscore,
  z_threshold = 1.96
) {
  rna_df <- deg_zscore |>
    dplyr::mutate(rna_sig = tag_direction(zscore, z_threshold)) |>
    dplyr::select(Species, Protein_id, Treatment, rna_z = zscore, rna_sig)

  prot_df <- dep_zscore |>
    dplyr::mutate(prot_sig = tag_direction(zscore, z_threshold)) |>
    dplyr::select(Species, Protein_id, Treatment, prot_z = zscore, prot_sig)

  joined <- dplyr::inner_join(rna_df, prot_df, by = c("Species", "Protein_id", "Treatment")) |>
    dplyr::distinct()
  protein_pattern <- assign_clusters(joined)

  joined |>
    dplyr::inner_join(protein_pattern, by = c("Species", "Protein_id", "Treatment")) |>
    dplyr::select(Species, Protein_id, Treatment,
                  rna_z, prot_z, rna_sig, prot_sig, cluster) |>
    dplyr::arrange(Species, cluster, Protein_id, Treatment)
}

# ── Plotting — scatter ─────────────────────────────────────────────────────────

#' Build One Discordance Scatter Plot for One Species (Internal)
#'
#' @description
#' Shared helper used by both `plot_discordance_scatter()` and
#' `plot_discordance_scatter_zscore()`. Renders background (grey) points and
#' foreground (cluster-coloured) star glyphs.
#'
#' @param bg            Tibble. Background (concordant) data for one species.
#' @param fg            Tibble. Foreground (discordant) data for one species.
#' @param x_col         Character. Column name for x-axis values.
#' @param y_col         Character. Column name for y-axis values.
#' @param x_lim         Numeric vector c(min, max). x-axis limits.
#' @param y_lim         Numeric vector c(min, max). y-axis limits.
#' @param species       Character. Species name for the plot title.
#' @param x_lab         Character. x-axis label.
#' @param y_lab         Character. y-axis label.
#' @param cluster_order Character vector. Draw order for clusters (largest first
#'   so smaller clusters are not hidden).
#'
#' @return ggplot object.
.build_scatter_plot <- function (
  bg, fg, x_col, y_col, x_lim, y_lim, species,
  x_lab, y_lab, cluster_order
) {
  n_discordant <- dplyr::n_distinct(fg$Protein_id)
  count_label <- paste0("n = ", n_discordant, " discordant")

  fg <- dplyr::mutate(fg, cluster = factor(cluster, levels = cluster_order))

  p <- ggplot2::ggplot() +
    ggplot2::geom_blank(
      data = tibble::tibble(x = x_lim, y = y_lim),
      ggplot2::aes(x = x, y = y)
    ) +
    ggplot2::geom_point(
      data = bg,
      ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]]),
      colour = "grey75", alpha = 0.4, size = 0.8
    )

  for (cl in cluster_order) {
    p <- p + ggplot2::geom_point(
      data = dplyr::filter(fg, cluster == cl),
      ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]], colour = cluster),
      alpha = 0.7, size = 1.5
    )
  }

  p +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30", linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30", linewidth = 0.5) +
    ggplot2::annotate(
      "text",
      x = x_lim[1], y = y_lim[2],
      label = count_label,
      hjust = -0.25, vjust = 2.5, size = 8,
      colour = "grey30", fontface = "italic"
    ) +
    ggplot2::scale_color_manual(
      values = CLUSTER_COLOURS,
      labels = CLUSTER_LABELS,
      name = "mRNA – Protein",
      breaks = names(CLUSTER_COLOURS),
      drop = FALSE
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        override.aes = list(shape = 16, size = 3, alpha = 1)
      )
    ) +
    ggplot2::scale_x_continuous(
      limits = x_lim,
      breaks = scales::pretty_breaks(n = 5),
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      limits = y_lim,
      breaks = scales::pretty_breaks(n = 5),
      expand = c(0, 0)
    ) +
    ggplot2::labs(title = format_species_title(species), x = x_lab, y = y_lab) +
    theme_publication(base_size = 10, legend_position = "right") +
    ggplot2::theme(
      legend.spacing = ggplot2::unit(0.2, "cm"),
      legend.margin = ggplot2::margin(0.25, 0.25, 0.25, 0.25),
      legend.box.spacing = ggplot2::unit(0.2, "cm")
    )
}

#' Plot mRNA vs Protein logFC Scatter — Discordance Clusters
#'
#' @param discordance_df   Tibble. Output of `build_discordance_clusters()`.
#' @param dep_raw          Tibble. Full DEP raw data (background points).
#' @param deg_raw          Tibble. Full DEG raw data (background points).
#' @param output_path      Character. Output directory.
#' @param xlim             Numeric vector or NULL. Shared x limits.
#' @param ylim             Numeric vector or NULL. Shared y limits.
#' @param save_per_species Logical. Save individual species plots. Default FALSE.
#' @param prefix           Character or NULL. Optional filename prefix.
#'
#' @return List of ggplot objects (one per species).
plot_discordance_scatter <- function (
  discordance_df, dep_raw, deg_raw, output_path,
  xlim = NULL, ylim = NULL,
  save_per_species = FALSE, prefix = NULL
) {
  all_logfc <- dplyr::inner_join(
    dplyr::select(deg_raw, Species, Protein_id, Treatment, rna_log2FC = log2FC),
    dplyr::select(dep_raw, Species, Protein_id, Treatment, prot_log2FC = log2FC),
    by = c("Species", "Protein_id", "Treatment")
  )

  background_df <- all_logfc |>
    dplyr::anti_join(dplyr::distinct(discordance_df, Species, Protein_id),
                     by = c("Species", "Protein_id")) |>
    dplyr::filter(is.finite(rna_log2FC), is.finite(prot_log2FC))

  foreground_df <- dplyr::filter(discordance_df, is.finite(rna_log2FC), is.finite(prot_log2FC))

  cluster_order <- foreground_df |>
    dplyr::count(cluster) |>
    dplyr::arrange(dplyr::desc(n)) |>
    dplyr::pull(cluster)

  purrr::map(sort(unique(foreground_df$Species)), function (species) {
    bg <- dplyr::filter(background_df, Species == species)
    fg <- dplyr::filter(foreground_df, Species == species)

    x_lim <- xlim %||% {
      r <- range(c(bg$rna_log2FC, fg$rna_log2FC), na.rm = TRUE)
      c(floor(r[1]), ceiling(r[2]))
    }
    y_lim <- ylim %||% {
      r <- range(c(bg$prot_log2FC, fg$prot_log2FC), na.rm = TRUE)
      c(floor(r[1]), ceiling(r[2]))
    }

    p <- .build_scatter_plot(bg, fg, "rna_log2FC", "prot_log2FC",
                             x_lim, y_lim, species,
                             "logFC(mRNA)", "logFC(Protein)", cluster_order)

    if (save_per_species) {
      out_file <- file.path(
        output_path, species,
        paste0(if (!is.null(prefix)) paste0(prefix, "_") else "", "discordance_scatter_logFC")
      )
      save_plot(p, out_file, width = 6, height = 5)
      message("Discordance scatter saved: ", out_file, MSG_FIG_FORMAT)
    }
    p
  })
}

#' Plot mRNA vs Protein Z-Score Scatter — Discordance Clusters
#'
#' @param discordance_df   Tibble. Output of `build_discordance_clusters_zscore()`.
#' @param dep_zscore       Tibble. Full DEP Z-score data (background points).
#' @param deg_zscore       Tibble. Full DEG Z-score data (background points).
#' @param output_path      Character. Output directory.
#' @param zlim             Numeric vector or NULL. Shared axis limits.
#' @param save_per_species Logical. Default FALSE.
#' @param prefix           Character or NULL. Optional filename prefix.
#'
#' @return List of ggplot objects (one per species).
plot_discordance_scatter_zscore <- function (
  discordance_df, dep_zscore, deg_zscore, output_path,
  zlim = NULL, save_per_species = FALSE, prefix = NULL
) {
  all_z <- dplyr::inner_join(
    dplyr::select(deg_zscore, Species, Protein_id, Treatment, rna_z = zscore),
    dplyr::select(dep_zscore, Species, Protein_id, Treatment, prot_z = zscore),
    by = c("Species", "Protein_id", "Treatment")
  )

  background_df <- all_z |>
    dplyr::anti_join(dplyr::distinct(discordance_df, Species, Protein_id),
                     by = c("Species", "Protein_id")) |>
    dplyr::filter(is.finite(rna_z), is.finite(prot_z))

  foreground_df <- dplyr::filter(discordance_df, is.finite(rna_z), is.finite(prot_z))

  cluster_order <- foreground_df |>
    dplyr::count(cluster) |>
    dplyr::arrange(dplyr::desc(n)) |>
    dplyr::pull(cluster)

  all_rna_z <- c(background_df$rna_z, foreground_df$rna_z)
  all_prot_z <- c(background_df$prot_z, foreground_df$prot_z)
  default_lim <- c(
    floor(min(c(all_rna_z, all_prot_z), na.rm = TRUE)),
    ceiling(max(c(all_rna_z, all_prot_z), na.rm = TRUE))
  )

  purrr::map(sort(unique(foreground_df$Species)), function (species) {
    bg <- dplyr::filter(background_df, Species == species)
    fg <- dplyr::filter(foreground_df, Species == species)
    lim <- zlim %||% default_lim

    p <- .build_scatter_plot(bg, fg, "rna_z", "prot_z",
                             lim, lim, species,
                             "Z-score (mRNA)", "Z-score (Protein)", cluster_order)

    if (save_per_species) {
      out_file <- file.path(
        output_path, species,
        paste0(if (!is.null(prefix)) paste0(prefix, "_") else "", "discordance_scatter")
      )
      save_plot(p, out_file, width = 6, height = 5)
      message("Z-score scatter saved: ", out_file, MSG_FIG_FORMAT)
    }
    p
  })
}

# ── Plotting — profiles ────────────────────────────────────────────────────────

#' Build a Discordance Profile Line Plot (Internal)
#'
#' @description
#' Shared helper for `plot_discordance_profiles()`. Renders a spaghetti-style
#' line plot of mRNA and protein values per cluster and optionally per species.
#'
#' @param df               Tibble. Long-format data with columns: `Protein_id`,
#'   `measure`, `Treatment`, `cluster`, and the value column (`y_col`).
#' @param count_df         Tibble. Protein counts per cluster and species.
#' @param y_col            Character. Value column to plot on y-axis.
#' @param y_lab            Character. y-axis label.
#' @param treatment_order  Character vector. Treatment display order.
#' @param species_filter   Character or NULL. If not NULL, show only this species.
#' @param species_labeller Named vector for `as_labeller()`.
#' @param threshold_band   Numeric or NULL. If provided, draws a shaded band at
#'   ±`threshold_band` (useful for Z-score plots). Default NULL.
#' @param cluster_spacing  Numeric. Horizontal spacing between cluster panels.
#'   Default 0.2.
#' @param species_spacing  Numeric. Vertical spacing between species rows.
#'   Default 0.5.
#'
#' @return ggplot object.
.build_profile_plot <- function (
  df, count_df, y_col, y_lab,
  treatment_order, species_filter, species_labeller,
  threshold_band = NULL,
  cluster_spacing = 0.2,
  species_spacing = 0.5
) {
  facet <- if (is.null(species_filter)) {
    ggplot2::facet_grid(
      Species ~ cluster,
      scales = "free_y",
      space = "fixed",
      labeller = ggplot2::labeller(
        Species = species_labeller,
        cluster = ggplot2::as_labeller(CLUSTER_LABELS)
      )
    )
  } else {
    ggplot2::facet_wrap(
      ~cluster, nrow = 1, scales = "free_y",
      labeller = ggplot2::as_labeller(CLUSTER_LABELS)
    )
  }

  base_size <- 24

  p <- ggplot2::ggplot(df, ggplot2::aes(
    x = Treatment,
    y = .data[[y_col]],
    group = interaction(Protein_id, measure),
    colour = measure
  )) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        colour = "grey40", linewidth = 0.4)

  if (!is.null(threshold_band)) {
    p <- p + ggplot2::annotate(
      "rect",
      xmin = -Inf, xmax = Inf,
      ymin = -threshold_band, ymax = threshold_band,
      alpha = 0.1, fill = "grey50"
    )
  }

  p +
    ggplot2::geom_line(alpha = 0.3, linewidth = 0.3) +
    ggplot2::geom_text(
      data = count_df,
      ggplot2::aes(
        x = length(treatment_order),
        y = Inf,
        label = paste0("n=", n)
      ),
      hjust = 1.1, vjust = 1.5, size = 8,
      colour = "grey30", fontface = "italic",
      inherit.aes = FALSE
    ) +
    facet +
    ggplot2::scale_colour_manual(
      values = c("mRNA" = "#E8C84A", "Protein" = "#7BBFDF"),
      name = NULL
    ) +
    ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = 45)) +
    ggplot2::labs(x = NULL, y = y_lab) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      strip.text.x = ggplot2::element_text(size = base_size + 2, face = "bold", colour = "white"),
      strip.text.y = ggplot2::element_text(size = base_size + 2, face = "italic", angle = -90),
      strip.background.x = ggplot2::element_rect(fill = "grey40"),
      strip.background.y = ggplot2::element_rect(fill = "white"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.spacing.x = ggplot2::unit(cluster_spacing, "cm"),
      panel.spacing.y = ggplot2::unit(species_spacing, "cm"),
      legend.position = "right",
      legend.text = ggplot2::element_text(size = base_size + 2),
      legend.key.size = ggplot2::unit(0.4, "cm"),
      axis.text.x = ggplot2::element_text(size = base_size),
      axis.text.y = ggplot2::element_text(size = base_size),
      axis.title.y = ggplot2::element_text(size = base_size + 2)
    )
}

#' Build Species Labeller for Profile Plots (Internal)
.species_labeller <- function (species_levels) {
  ggplot2::as_labeller(setNames(
    sapply(species_levels, function (s) {
      display <- PLOT_SPECIES_NAMES[s]
      ifelse(is.na(display), gsub("_", " ", s), display)
    }),
    species_levels
  ))
}

#' Plot Discordance Profiles per Cluster and Species
#'
#' @description
#' Spaghetti-style line plot showing per-protein mRNA and protein values across
#' stress treatments, faceted by discordance cluster. Supports both logFC
#' (from `build_discordance_clusters()`) and Z-score
#' (from `build_discordance_clusters_zscore()`) inputs via the `value_col`
#' argument. Optionally draws a shaded threshold band for Z-score plots.
#'
#' @param discordance_df Tibble. Output of either cluster-building function.
#'   Must contain columns: `Species`, `Protein_id`, `Treatment`, `cluster`,
#'   plus the two measure columns identified by `rna_col` and `prot_col`.
#' @param output_path    Character. Output directory (currently unused but
#'   retained for API consistency).
#' @param rna_col        Character. Column name for mRNA values.
#'   Default `"rna_log2FC"`.
#' @param prot_col       Character. Column name for protein values.
#'   Default `"prot_log2FC"`.
#' @param y_lab          Character. y-axis label. Default `"Log2 FC"`.
#' @param species_filter Character or NULL. If not NULL, plot only this species.
#' @param threshold_band Numeric or NULL. If provided, draws ±band shading
#'   (useful for Z-score approach). Default NULL.
#' @param cluster_spacing Numeric. Panel spacing. Default 0.2.
#' @param species_spacing Numeric. Row spacing. Default 0.5.
#'
#' @return A `wrap_elements` gtable with cluster-coloured strip backgrounds.
plot_discordance_profiles <- function (
  discordance_df,
  output_path,
  rna_col = "rna_z",
  prot_col = "prot_z",
  y_lab = "Log2 FC",
  species_filter = NULL,
  threshold_band = NULL,
  cluster_spacing = 0.2,
  species_spacing = 0.5
) {
  treatment_order <- sort(unique(
    discordance_df$Treatment[discordance_df$Treatment != "Ctrl"]
  ))

  filter_species <- function (df) {
    if (!is.null(species_filter)) dplyr::filter(df, Species == species_filter) else df
  }

  df <- discordance_df |>
    dplyr::filter(Treatment != "Ctrl") |>
    filter_species() |>
    dplyr::mutate(
      Treatment = factor(Treatment, levels = treatment_order),
      cluster = factor(cluster, levels = names(CLUSTER_COLOURS))
    ) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(c(rna_col, prot_col)),
      names_to = "measure",
      values_to = "value"
    ) |>
    dplyr::mutate(
      measure = factor(
        ifelse(measure == rna_col, "mRNA", "Protein"),
        levels = c("mRNA", "Protein")
      )
    )

  count_df <- discordance_df |>
    dplyr::filter(Treatment != "Ctrl") |>
    filter_species() |>
    dplyr::distinct(Species, Protein_id, cluster) |>
    dplyr::count(Species, cluster, name = "n") |>
    dplyr::mutate(cluster = factor(cluster, levels = names(CLUSTER_COLOURS)))

  species_levels <- sort(unique(df$Species))

  p <- .build_profile_plot(
    df, count_df, "value", y_lab,
    treatment_order, species_filter,
    .species_labeller(species_levels),
    threshold_band = threshold_band,
    cluster_spacing = cluster_spacing,
    species_spacing = species_spacing
  )

  g <- colour_cluster_strips(ggplot2::ggplot_gtable(ggplot2::ggplot_build(p)))
  patchwork::wrap_elements(g)
}

# ── Plotting — DE per cluster ──────────────────────────────────────────────────

#' Plot DEG vs DEP Percentage Bar Chart per Discordance Cluster
#'
#' @param deg_df      Tibble. DEG summary with `Species`, `Treatment`,
#'   `cluster`, `pct_de`.
#' @param dep_df      Tibble. DEP summary with the same structure.
#' @param output_path Character. Output directory.
#' @param prefix      Character or NULL. Optional filename prefix.
#'
#' @return List of ggplot objects (one per species).
plot_de_barplot_per_cluster <- function (deg_df, dep_df, output_path, prefix = NULL) {
  cluster_labels <- setNames(
    CLUSTER_LABELS[unique(deg_df$cluster)],
    unique(deg_df$cluster)
  )

  plots <- dplyr::bind_rows(
    dplyr::mutate(deg_df, type = "DEGs"),
    dplyr::mutate(dep_df, type = "DEPs")
  ) |>
    dplyr::filter(Treatment %in% STRESS_TREATMENTS_PROT) |>
    dplyr::mutate(
      Treatment = factor(Treatment, levels = TREATMENT_ORDER),
      type = factor(type, levels = c("DEGs", "DEPs"))
    ) |>
    dplyr::group_by(Species) |>
    dplyr::group_map(function (sp_df, sp_key) {
      species <- sp_key$Species
      x_max <- max(ceiling(max(sp_df$pct_de, na.rm = TRUE) / 10) * 10, 20)

      ggplot2::ggplot(sp_df, ggplot2::aes(x = pct_de, y = Treatment, fill = type)) +
        ggplot2::facet_wrap(
          ~cluster, nrow = 2,
          labeller = ggplot2::as_labeller(cluster_labels)
        ) +
        ggplot2::geom_col(
          position = ggplot2::position_dodge(width = 0.7, reverse = TRUE),
          width = 0.6,
          colour = NA
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
          x = "Percentage",
          y = "Stress Conditions"
        ) +
        theme_publication()
    })

  purrr::walk2(plots, unique(deg_df$Species), function (p, species) {
    out_file <- file.path(
      output_path, species,
      paste0(prefix, "discordance_DEG_DEP_barplot")
    )
    save_plot(p, out_file, width = 16, height = 5)
    message("  Plot saved: ", out_file, MSG_FIG_FORMAT)
  })

  plots
}

# ── Enrichment ─────────────────────────────────────────────────────────────────

#' Load Enrichment Annotation Files for One Species
#'
#' @param species Character. Species name.
#'
#' @return Named list: `term2gene`, `universe`.
load_enrichment_data <- function (species) {
  abbv <- switch(
    species,
    Salmonella_enterica = "SE",
    Staphylococcus_aureus = "SA",
    Yersinia_pseudotuberculosis = "YP",
    stop("No enrichment files configured for: ", species)
  )

  list(
    term2gene = readxl::read_excel(
      file.path(ENRICHMENT, paste0("TERM2GENE_", abbv, ".xlsx")),
      col_names = TRUE
    ),
    universe = readr::read_csv(
      file.path(ENRICHMENT, paste0("universe_", abbv, ".csv")),
      col_names = TRUE,
      show_col_types = FALSE
    ) |> dplyr::pull(1)
  )
}

#' Run clusterProfiler Enrichment for One Gene Set
#'
#' @param gene_ids    Character vector. Protein IDs to test.
#' @param term2gene   Tibble. TERM / GENE columns.
#' @param term2name   Tibble. TERM / ONTOLOGY columns.
#' @param universe    Character vector. Background gene set.
#' @param pval_cutoff Numeric. p-value cutoff. Default 0.05.
#' @param qval_cutoff Numeric. q-value cutoff. Default 0.2.
#' @param min_gs_size Integer. Minimum gene-set size. Default 5.
#' @param max_gs_size Integer. Maximum gene-set size. Default 500.
#'
#' @return Tibble of enrichment results, or NULL if nothing significant.
run_enricher <- function (
  gene_ids,
  term2gene,
  term2name,
  universe,
  pval_cutoff = 0.05,
  qval_cutoff = 0.2,
  min_gs_size = 5,
  max_gs_size = 500
) {
  if (length(gene_ids) == 0) {
    message("    No genes to enrich — skipping.")
    return(NULL)
  }

  result <- tryCatch(
    clusterProfiler::enricher(
      gene = gene_ids,
      TERM2GENE = term2gene,
      TERM2NAME = term2name,
      universe = universe,
      pvalueCutoff = pval_cutoff,
      qvalueCutoff = qval_cutoff,
      minGSSize = min_gs_size,
      maxGSSize = max_gs_size
    ),
    error = function (e) { message("    enricher error: ", e$message); NULL }
  )

  if (is.null(result) || nrow(result@result) == 0) {
    message("    No significant terms found.")
    return(NULL)
  }
  result@result
}

#' Run Enrichment Analysis per Cluster and Species
#'
#' @param discordance_df Tibble. Output of `build_discordance_clusters()`.
#' @param pval_cutoff    Numeric. Default 0.05.
#' @param qval_cutoff    Numeric. Default 0.2.
#' @param min_gs_size    Integer. Default 5.
#' @param max_gs_size    Integer. Default 500.
#'
#' @return Tibble of all enrichment results with `Species`, `Cluster`,
#'   `Protein_ID`, and GO-term annotation columns.
run_enrichment_all <- function (
  discordance_df,
  pval_cutoff = 0.05,
  qval_cutoff = 0.2,
  min_gs_size = 5,
  max_gs_size = 500
) {
  species_list <- unique(discordance_df$Species)
  clusters <- sort(unique(discordance_df$cluster))

  term2name <- readxl::read_excel(
    file.path(ENRICHMENT, "GO_terms_from_database.xlsx")
  ) |>
    dplyr::filter(Obsolete == FALSE)

  purrr::map_dfr(species_list, function (species) {
    message("\nRunning enrichment for: ", species)
    annot <- load_enrichment_data(species)

    purrr::map_dfr(clusters, function (cl) {
      message("  Cluster: ", cl)
      gene_ids <- discordance_df |>
        dplyr::filter(Species == species, cluster == cl) |>
        dplyr::distinct(Protein_id) |>
        dplyr::pull(Protein_id)

      message("    Genes: ", length(gene_ids))

      run_enricher(
        gene_ids = gene_ids,
        term2gene = annot$term2gene,
        term2name = term2name,
        universe = annot$universe,
        pval_cutoff = pval_cutoff,
        qval_cutoff = qval_cutoff,
        min_gs_size = min_gs_size,
        max_gs_size = max_gs_size
      ) |>
        dplyr::mutate(Cluster = cl)
    }) |>
      dplyr::mutate(Species = species)
  }) |>
    dplyr::inner_join(term2name, by = "ID") |>
    dplyr::mutate(Protein_ID = geneID)
}

#' Plot Enrichment Dot Plot — All Clusters and Species
#'
#' @param enrichment_results Tibble. Output of `run_enrichment_all()`.
#' @param fdr_threshold      Numeric. FDR threshold for displaying terms.
#'   Default 0.05.
#' @param top_n              Integer or NULL. Max terms per cluster; NULL = all.
#' @param output_path        Character or NULL. Directory for saving the plot.
#' @param cluster_spacing    Numeric. Horizontal spacing between clusters.
#'   Default 0.2.
#'
#' @return `wrap_elements` gtable object (invisibly).
plot_enrichment_dotplot <- function (
  enrichment_results,
  fdr_threshold = 0.05,
  top_n = NULL,
  output_path = NULL,
  cluster_spacing = 0.2
) {
  plot_df <- enrichment_results |>
    dplyr::filter(namespace == "Biological Process") |>
    dplyr::mutate(
      neg_log10_fdr = -log10(p.adjust),
      cluster = factor(Cluster, levels = names(CLUSTER_COLOURS)),
      Species = factor(Species, levels = sort(unique(Species)))
    )

  if (nrow(plot_df) == 0) {
    warning("No significant enrichment results to plot.")
    return(invisible(NULL))
  }

  if (!is.null(top_n)) {
    plot_df <- plot_df |>
      dplyr::group_by(cluster, Description) |>
      dplyr::summarise(min_fdr = min(p.adjust), .groups = "drop") |>
      dplyr::group_by(cluster) |>
      dplyr::slice_min(min_fdr, n = top_n) |>
      dplyr::select(cluster, Description) |>
      dplyr::inner_join(plot_df, by = c("cluster", "Description"))
  }

  sig_terms <- dplyr::distinct(dplyr::filter(plot_df, p.adjust < fdr_threshold), Description)
  plot_df <- dplyr::semi_join(plot_df, sig_terms, by = "Description")

  term_order <- sort(unique(plot_df$Description))
  plot_df <- dplyr::mutate(plot_df, Description = factor(Description, levels = rev(term_order)))

  species_levels <- levels(plot_df$Species)
  species_labels <- setNames(
    lapply(species_levels, format_species_title, linebreak = TRUE),
    species_levels
  )

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = Species, y = Description, size = Count, colour = neg_log10_fdr)
  ) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::geom_hline(
      yintercept = seq(0.5, length(term_order) + 0.5, 1),
      colour = "grey92", linewidth = 0.3
    ) +
    ggplot2::facet_wrap(
      ~cluster, nrow = 1,
      labeller = ggplot2::as_labeller(CLUSTER_LABELS)
    ) +
    ggplot2::scale_colour_gradient(
      low = "lightsalmon",
      high = "darkred",
      name = "-log10(FDR)",
      guide = ggplot2::guide_colourbar(barwidth = 0.8, barheight = 4, ticks = FALSE)
    ) +
    ggplot2::scale_size_continuous(
      name = "Count",
      range = c(3, 12),
      breaks = c(15, 30, 45, 60, 75, 90)
    ) +
    ggplot2::scale_x_discrete(labels = species_labels) +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_publication(base_size = 20, legend_position = "right") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      strip.background = ggplot2::element_rect(fill = "grey40"),
      panel.spacing.x = ggplot2::unit(cluster_spacing, "cm")
    )

  g <- colour_cluster_strips(ggplot2::ggplot_gtable(ggplot2::ggplot_build(p)))
  final <- patchwork::wrap_elements(g)

  if (!is.null(output_path)) {
    n_terms <- dplyr::n_distinct(plot_df$Description)
    n_clusters <- dplyr::n_distinct(plot_df$cluster)
    save_plot(
      final,
      filepath = file.path(output_path, "enrichment_dotplot"),
      width = max(24, n_clusters * 3),
      height = max(12, n_terms * 0.35 + 2)
    )
    message(
      "Enrichment dotplot saved: ",
      file.path(output_path, "enrichment_dotplot"),
      MSG_FIG_FORMAT
    )
  }

  invisible(final)
}
