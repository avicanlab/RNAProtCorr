# Discordant mRNA-Protein regulation pattern ---------------------------------
#
# Description:
#   Identifies proteins with discordant mRNA-protein regulation patterns
#   across stress conditions using either logFC (limma) or Z-score approach.
#
# Inputs:
#   - deg_raw    : Raw DEG results (Species, Protein_id, Treatment, log2FC, fdr)
#   - dep_raw    : Raw DEP results (Species, Protein_id, Treatment, logFC, adj.P.Val)
#   - tpm_data   : TPM expression data for Z-score approach
#   - protQ_data : Protein quantification data for Z-score approach
#
# Requirements:
#   library(ggstar)
#   library(clusterProfiler)

library(ggstar)
library(clusterProfiler)

CLUSTER_COLOURS <- c(
  "Cluster 1" = "#E8A838",
  "Cluster 2" = "#E87070",
  "Cluster 3" = "#C47EC4",
  "Cluster 4" = "#4BBFBF",
  "Cluster 5" = "#4BAF6E",
  "Cluster 6" = "#E8A070"
)

CLUSTER_LABELS <- c(
  "Cluster 1" = "down \u2013 no change",
  "Cluster 2" = "down \u2013 up",
  "Cluster 3" = "no change \u2013 down",
  "Cluster 4" = "no change \u2013 up",
  "Cluster 5" = "up \u2013 down",
  "Cluster 6" = "up \u2013 no change"
)

# ============================================================================
# SHARED HELPERS
# ============================================================================

#' Tag regulation direction from a numeric value and thresholds
#'
#' @param value      Numeric vector.
#' @param threshold  Num. Absolute threshold.
#' @param pval       Numeric vector of p-values (NULL for Z-score approach).
#' @param pval_threshold Num. P-value threshold (ignored if pval is NULL).
#' @return Chr vector: "down", "up", or "unchanged"
tag_direction <- function(value, threshold, pval = NULL, pval_threshold = NULL) {
  if (!is.null(pval)) {
    case_when(
      pval < pval_threshold & value < -threshold ~ "down",
      pval < pval_threshold & value > threshold ~ "up",
      TRUE ~ "unchanged"
    )
  } else {
    case_when(
      value < -threshold ~ "down",
      value > threshold ~ "up",
      TRUE ~ "unchanged"
    )
  }
}

#' Assign discordance cluster from per-protein regulation pattern
#'
#' @param joined Tibble with rna_sig and prot_sig columns, grouped by protein.
#' @param direction_method Chr. "strict" or "dominant".
#' @return Tibble with cluster column added.
assign_clusters <- function(joined) {
  protein_pattern <- joined %>%
    group_by(Species, Protein_id, Treatment) %>%
    summarise(
      rna_ever_down = any(rna_sig == "down"),
      rna_ever_up = any(rna_sig == "up"),
      rna_ever_unchanged = any(rna_sig == "unchanged"),
      prot_ever_down = any(prot_sig == "down"),
      prot_ever_up = any(prot_sig == "up"),
      prot_ever_unchanged = any(prot_sig == "unchanged"),
      .groups = "drop"
    )
  protein_pattern <- protein_pattern %>%
    mutate(
      cluster = case_when(
        rna_ever_down & !rna_ever_up & !rna_ever_unchanged & !prot_ever_down & !prot_ever_up & prot_ever_unchanged ~ "Cluster 1",
        rna_ever_down & !rna_ever_up & !rna_ever_unchanged & !prot_ever_down & prot_ever_up & !prot_ever_unchanged ~ "Cluster 2",
        !rna_ever_down & !rna_ever_up & rna_ever_unchanged & prot_ever_down & !prot_ever_up & !prot_ever_unchanged ~ "Cluster 3",
        !rna_ever_down & !rna_ever_up & rna_ever_unchanged & !prot_ever_down & prot_ever_up & !prot_ever_unchanged ~ "Cluster 4",
        !rna_ever_down & rna_ever_up & !rna_ever_unchanged & prot_ever_down & !prot_ever_up & !prot_ever_unchanged ~ "Cluster 5",
        !rna_ever_down & rna_ever_up & !rna_ever_unchanged & !prot_ever_down & !prot_ever_up & prot_ever_unchanged ~ "Cluster 6",
        TRUE ~ NA_character_
      )
    )
  protein_pattern <- protein_pattern %>%
    filter(!is.na(cluster)) %>%
    dplyr::select(Species, Protein_id, Treatment, cluster)
}

#' Colour cluster strip backgrounds in a gtable
colour_cluster_strips <- function(g) {
  strip_idx <- which(grepl("strip-t", g$layout$name))
  cluster_names <- names(CLUSTER_COLOURS)
  for (i in seq_along(strip_idx)) {
    col <- CLUSTER_COLOURS[cluster_names[((i - 1) %% length(cluster_names)) + 1]]
    g$grobs[[strip_idx[i]]]$grobs[[1]]$children[[1]]$
      gp$
      fill <- col
  }
  g
}

# ============================================================================
# LOGFC APPROACH
# ============================================================================

#' Build Discordance Clusters from logFC (limma DEG/DEP)
#'
#' @description
#' Joins RNA and protein logFC per protein x treatment, tags regulation
#' direction, and assigns each protein to a discordance cluster based on
#' its pattern across all stress conditions.
#'
#' Cluster definitions:
#'   1: mRNA down (only) in >= 1 condition, protein unchanged in all
#'   2: mRNA down (only) in >= 1 condition, protein up in >= 1 condition
#'   3: mRNA unchanged in all conditions, protein down in >= 1 condition
#'   4: mRNA unchanged in all conditions, protein up in >= 1 condition
#'   5: mRNA up (only) in >= 1 condition, protein down in >= 1 condition
#'   6: mRNA up (only) in >= 1 condition, protein unchanged in all
#'
#' @param deg_raw        Tibble. Species, Protein_id, Treatment, log2FC, fdr.
#' @param dep_raw        Tibble. Species, Protein_id, Treatment, logFC, adj.P.Val.
#' @param fc_threshold   Num. |log2FC| threshold. Default DE_FC_THRESHOLD.
#' @param pval_threshold Num. Adjusted p-value threshold. Default DEP_PVAL_THRESHOLD.
#' @param direction_method Chr. "strict" (default) or "dominant".
#'
#' @return Tibble: Species, Protein_id, Treatment, rna_log2FC, prot_log2FC,
#'   rna_sig, prot_sig, cluster
build_discordance_clusters <- function(
  deg_raw,
  dep_raw,
  fc_threshold = DE_FC_THRESHOLD,
  pval_threshold = DEP_PVAL_THRESHOLD
) {

  rna_df <- deg_raw %>%
    filter(!is.na(log2FC), !is.na(fdr)) %>%
    mutate(rna_sig = tag_direction(log2FC, fc_threshold, fdr, pval_threshold)) %>%
    dplyr::select(Species, Protein_id, Treatment, rna_log2FC = log2FC, rna_sig)

  prot_df <- dep_raw %>%
    mutate(prot_sig = tag_direction(logFC, fc_threshold, adj.P.Val, pval_threshold)) %>%
    dplyr::select(Species, Protein_id, Treatment, prot_log2FC = logFC, prot_sig)

  joined <- inner_join(rna_df, prot_df, by = c("Species", "Protein_id", "Treatment"))
  protein_pattern <- assign_clusters(joined)

  joined %>%
    inner_join(protein_pattern, by = c("Species", "Protein_id")) %>%
    dplyr::select(Species, Protein_id, Treatment,
                  rna_log2FC, prot_log2FC, rna_sig, prot_sig, cluster) %>%
    arrange(Species, cluster, Protein_id, Treatment)
}

#' Extract Raw Per-Protein DEP Results from limma
#'
#' @param mat            Numeric matrix. Rows = proteins, cols = samples.
#' @param metadata       Tibble. sample, Replicate, Treatment.
#' @param species        Chr. Species name.
#' @param pval_threshold Num. BH p-value threshold.
#' @param fc_threshold   Num. |logFC| threshold.
#'
#' @return Tibble: Species, Protein_id, Treatment, logFC, adj.P.Val
extract_dep_raw <- function(
  mat,
  metadata,
  species,
  pval_threshold = DEP_PVAL_THRESHOLD,
  fc_threshold = DE_FC_THRESHOLD
) {
  stress_present <- intersect(STRESS_TREATMENTS_PROT, unique(metadata$Treatment))

  metadata <- metadata %>%
    mutate(
      Treatment = factor(Treatment, levels = c("Ctrl", stress_present)),
      Replicate = factor(Replicate)
    )

  design <- model.matrix(~0 + Treatment + Replicate, data = metadata)
  colnames(design) <- str_replace(colnames(design), "^Treatment", "")
  contrast_mat <- makeContrasts(
    contrasts = paste0(stress_present, " - Ctrl"),
    levels = design
  )
  colnames(contrast_mat) <- stress_present
  fit2 <- eBayes(contrasts.fit(lmFit(mat, design), contrast_mat))

  map_dfr(stress_present, function(trt) {
    topTable(fit2, coef = trt, number = Inf,
             sort.by = "none", adjust.method = "BH") %>%
      rownames_to_column("Protein_id") %>%
      mutate(Species = species, Treatment = trt) %>%
      dplyr::select(Species, Protein_id, Treatment, logFC, adj.P.Val)
  })
}

# ============================================================================
# ZSCORE APPROACH
# ============================================================================

#' Z-Score Transform Expression Data Relative to Control
#'
#' @description
#' Computes per-protein Z-score using control mean and global median SD:
#'   z = (X - mean_ctrl) / global_median_sd
#' Proteins with zero/infinite control mean or SD are excluded.
#' Control samples are excluded from the output.
#'
#' @param df         Tibble. Must contain Species, Protein_id, Treatment,
#'   and the value column.
#' @param col_name   Chr. Name of the value column to transform.
#' @param fudge_factor Num. Added to global_median_sd for stability. Default 0.
#'
#' @return Tibble: Species, Protein_id, Treatment, value, mean_ctrl,
#'   sd_ctrl, zscore
zscore_transform <- function(df, col_name, fudge_factor = 0) {
  mean_sd_df <- df %>%
    filter(is.finite(.data[[col_name]])) %>%
    group_by(Species) %>%
    mutate(
      mean_ctrl = mean(.data[[col_name]][Treatment == "Ctrl"], na.rm = TRUE),
      sd_ctrl = sd(.data[[col_name]][Treatment == "Ctrl"], na.rm = TRUE)
    ) %>%
    ungroup() %>%
    filter(Treatment != "Ctrl")

  # Remove proteins with any infinite stress values
  infinite_ids <- mean_sd_df %>%
    filter(is.infinite(.data[[col_name]])) %>%
    distinct(Species, Protein_id)

  finite_df <- mean_sd_df %>% anti_join(infinite_ids, by = c("Species", "Protein_id"))
  global_median_sd <- median(finite_df$sd_ctrl, na.rm = TRUE) + fudge_factor

  finite_df %>%
    mutate(zscore = (.data[[col_name]] - mean_ctrl) / sd_ctrl) %>%
    dplyr::select(Species, Protein_id, Treatment,
                  value = all_of(col_name), mean_ctrl, sd_ctrl, zscore)
}

#' Build Discordance Clusters from Z-Scores
#'
#' @description
#' Assigns proteins to discordance clusters using Z-score thresholds
#' instead of logFC + p-value. Same cluster definitions as
#' build_discordance_clusters().
#'
#' @param deg_zscore     Tibble. Output of zscore_transform() on TPM data.
#'   Must contain: Species, Protein_id, Treatment, zscore.
#' @param dep_zscore     Tibble. Output of zscore_transform() on protein data.
#' @param z_threshold    Num. |Z-score| threshold. Default 1.96 (~95% CI).
#'
#' @return Tibble: Species, Protein_id, Treatment, rna_z, prot_z,
#'   rna_sig, prot_sig, cluster
build_discordance_clusters_zscore <- function(
  deg_zscore,
  dep_zscore,
  z_threshold = 1.96
) {
  rna_df <- deg_zscore %>%
    mutate(rna_sig = tag_direction(zscore, z_threshold)) %>%
    dplyr::select(Species, Protein_id, Treatment, rna_z = zscore, rna_sig)

  prot_df <- dep_zscore %>%
    mutate(prot_sig = tag_direction(zscore, z_threshold)) %>%
    dplyr::select(Species, Protein_id, Treatment, prot_z = zscore, prot_sig)

  joined <- inner_join(
    rna_df, prot_df, by = c("Species", "Protein_id", "Treatment")
  ) %>%
    distinct()
  protein_pattern <- assign_clusters(joined)

  joined <- joined %>%
    inner_join(protein_pattern, by = c("Species", "Protein_id", "Treatment"))
  joined %>%
    dplyr::select(Species, Protein_id, Treatment,
                  rna_z, prot_z, rna_sig, prot_sig, cluster) %>%
    arrange(Species, cluster, Protein_id, Treatment)
}

# ============================================================================
# PLOTTING — SCATTER
# ============================================================================

#' Build one discordance scatter plot for one species
#'
#' Internal helper used by both logFC and Z-score scatter functions.
.build_scatter_plot <- function(
  bg, fg, x_col, y_col, x_lim, y_lim, species,
  x_lab, y_lab, cluster_order
) {
  n_discordant <- n_distinct(fg$Protein_id)
  count_label <- paste0("n = ", n_discordant, " discordant")

  fg <- fg %>% mutate(cluster = factor(cluster, levels = cluster_order))

  p <- ggplot() +
    geom_blank(
      data = tibble(x = x_lim, y = y_lim),
      aes(x = x, y = y)
    ) +
    geom_point(
      data = bg,
      aes(x = .data[[x_col]], y = .data[[y_col]]),
      colour = "grey75", alpha = 0.4, size = 0.8
    )

  for (cl in cluster_order) {
    p <- p + geom_star(
      data = fg %>% filter(cluster == cl),
      aes(x = .data[[x_col]], y = .data[[y_col]],
          fill = cluster, starshape = Treatment),
      alpha = 0.7, size = 1.5, colour = NA
    )
  }

  p +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30", linewidth = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30", linewidth = 0.5) +
    annotate("text", x = x_lim[1], y = y_lim[2], label = count_label,
             hjust = -0.25, vjust = 2.5, size = 8, colour = "grey30", fontface = "italic") +
    scale_fill_manual(values = CLUSTER_COLOURS, labels = CLUSTER_LABELS,
                      name = "mRNA \u2013 Protein", breaks = names(CLUSTER_COLOURS), drop = FALSE) +
    scale_starshape_manual(values = TREATMENT_SHAPES, name = "Treatment") +
    guides(
      starshape = guide_legend(override.aes = list(fill = "grey40", colour = "grey40", size = 3, alpha = 1)),
      fill = guide_legend(override.aes = list(starshape = 15, colour = "white", size = 3, alpha = 1))
    ) +
    scale_x_continuous(limits = x_lim, breaks = scales::pretty_breaks(n = 5), expand = c(0, 0)) +
    scale_y_continuous(limits = y_lim, breaks = scales::pretty_breaks(n = 5), expand = c(0, 0)) +
    labs(title = format_species_title(species), x = x_lab, y = y_lab) +
    theme_publication(base_size = 10, legend_position = "right") +
    theme(
      legend.spacing = unit(0.2, "cm"),      # space between legend keys
      legend.margin = margin(0.25, 0.25, 0.25, 0.25), # margin around each legend
      legend.box.spacing = unit(0.2, "cm")   # space between the two legends
    )
}

#' Plot mRNA vs Protein LogFC Scatter — Discordance Clusters
#'
#' @param discordance_df Tibble. Output of build_discordance_clusters().
#' @param dep_raw        Tibble. Full DEP raw data (for background points).
#' @param deg_raw        Tibble. Full DEG raw data (for background points).
#' @param output_path    Chr. Output directory.
#' @param xlim           Num vector c(min, max). Shared x limits. NULL = per-species.
#' @param ylim           Num vector c(min, max). Shared y limits. NULL = per-species.
#' @param save_per_species Logical. Save individual species plots. Default FALSE.
#' @param prefix         Chr. Optional filename prefix for per-species plots.
#'
#' @return List of ggplot objects (one per species).
plot_discordance_scatter <- function(
  discordance_df, dep_raw, deg_raw, output_path,
  xlim = NULL, ylim = NULL,
  save_per_species = FALSE, prefix = NULL
) {
  all_logfc <- inner_join(
    deg_raw %>% dplyr::select(Species, Protein_id, Treatment, rna_log2FC = log2FC),
    dep_raw %>% dplyr::select(Species, Protein_id, Treatment, prot_log2FC = logFC),
    by = c("Species", "Protein_id", "Treatment")
  )

  background_df <- all_logfc %>%
    anti_join(discordance_df %>% distinct(Species, Protein_id), by = c("Species", "Protein_id")) %>%
    filter(is.finite(rna_log2FC), is.finite(prot_log2FC))

  foreground_df <- discordance_df %>%
    filter(is.finite(rna_log2FC), is.finite(prot_log2FC))

  cluster_order <- foreground_df %>%
    count(cluster) %>%
    arrange(desc(n)) %>%
    pull(cluster)

  species_list <- sort(unique(foreground_df$Species))

  map(species_list, function(species) {
    bg <- background_df %>% filter(Species == species)
    fg <- foreground_df %>% filter(Species == species)

    x_lim <- xlim %||% { rna_range <- range(c(bg$rna_log2FC, fg$rna_log2FC), na.rm = TRUE); c(floor(rna_range[1]), ceiling(rna_range[2])) }
    y_lim <- ylim %||% { prot_range <- range(c(bg$prot_log2FC, fg$prot_log2FC), na.rm = TRUE); c(floor(prot_range[1]), ceiling(prot_range[2])) }

    p <- .build_scatter_plot(bg, fg, "rna_log2FC", "prot_log2FC",
                             x_lim, y_lim, species,
                             "logFC(mRNA)", "logFC(Protein)", cluster_order)

    if (save_per_species) {
      output_file <- file.path(output_path, species,
                               paste0(if (!is.null(prefix)) paste0(prefix, "_") else "", "discordance_scatter_logFC"))
      save_plot(p, output_file, width = 6, height = 5)
      message("Discordance scatter saved: ", output_file, MSG_FIG_FORMAT)
    }
    p
  })
}

#' Plot mRNA vs Protein Z-Score Scatter — Discordance Clusters
#'
#' @param discordance_df Tibble. Output of build_discordance_clusters_zscore().
#' @param dep_zscore     Tibble. Full DEP Z-score data (for background points).
#' @param deg_zscore     Tibble. Full DEG Z-score data (for background points).
#' @param output_path    Chr. Output directory.
#' @param zlim           Num vector c(min, max). Shared limits. NULL = per-species.
#' @param save_per_species Logical. Default FALSE.
#' @param prefix         Chr. Optional filename prefix.
#'
#' @return List of ggplot objects (one per species).
plot_discordance_scatter_zscore <- function(
  discordance_df, dep_zscore, deg_zscore, output_path,
  zlim = NULL, save_per_species = FALSE, prefix = NULL
) {

  all_z <- inner_join(
    deg_zscore %>% dplyr::select(Species, Protein_id, Treatment, rna_z = zscore),
    dep_zscore %>% dplyr::select(Species, Protein_id, Treatment, prot_z = zscore),
    by = c("Species", "Protein_id", "Treatment")
  )

  background_df <- all_z %>%
    anti_join(discordance_df %>% distinct(Species, Protein_id), by = c("Species", "Protein_id")) %>%
    filter(is.finite(rna_z), is.finite(prot_z))

  foreground_df <- discordance_df %>%
    filter(is.finite(rna_z), is.finite(prot_z))

  cluster_order <- foreground_df %>%
    count(cluster) %>%
    arrange(desc(n)) %>%
    pull(cluster)

  # Shared limits from all data combined
  all_rna_z <- c(background_df$rna_z, foreground_df$rna_z)
  all_prot_z <- c(background_df$prot_z, foreground_df$prot_z)
  default_lim <- c(
    floor(min(c(all_rna_z, all_prot_z), na.rm = TRUE)),
    ceiling(max(c(all_rna_z, all_prot_z), na.rm = TRUE))
  )

  species_list <- sort(unique(foreground_df$Species))

  map(species_list, function(species) {
    bg <- background_df %>% filter(Species == species)
    fg <- foreground_df %>% filter(Species == species)

    lim <- zlim %||% default_lim

    p <- .build_scatter_plot(bg, fg, "rna_z", "prot_z",
                             lim, lim, species,
                             "Z-score (mRNA)", "Z-score (Protein)", cluster_order)

    if (save_per_species) {
      output_file <- file.path(output_path, species,
                               paste0(if (!is.null(prefix)) paste0(prefix, "_") else "", "discordance_scatter"))
      save_plot(p, output_file, width = 6, height = 5)
      message("Z-score scatter saved: ", output_file, MSG_FIG_FORMAT)
    }
    p
  })
}

# ============================================================================
# PLOTTING — PROFILES
# ============================================================================

#' Build discordance profile plot (line plot per cluster)
#'
#' Internal helper shared by logFC and Z-score profile functions.
.build_profile_plot <- function(
  df, count_df, y_col, y_lab,
  treatment_order, species_filter,
  species_labeller,
  threshold_band = NULL,
  cluster_spacing = 0.2,
  species_spacing = 0.5
) {
  facet <- if (is.null(species_filter)) {
    facet_grid(
      Species ~ cluster,
      scales = "free_y",
      space = "fixed",   # equal column widths
      labeller = labeller(
        Species = species_labeller,
        cluster = as_labeller(CLUSTER_LABELS)
      )
    )
  } else {
    facet_wrap(
      ~cluster, nrow = 1, scales = "free_y",
      labeller = as_labeller(CLUSTER_LABELS)
    )
  }

  base_size <- 24

  p <- ggplot(df, aes(
    x = Treatment,
    y = .data[[y_col]],
    group = interaction(Protein_id, measure),
    colour = measure
  )) +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "grey40", linewidth = 0.4)

  if (!is.null(threshold_band)) {
    p <- p + annotate("rect",
                      xmin = -Inf, xmax = Inf,
                      ymin = -threshold_band, ymax = threshold_band,
                      alpha = 0.1, fill = "grey50")
  }

  p +
    geom_line(alpha = 0.3, linewidth = 0.3) +
    geom_text(
      data = count_df,
      aes(x = length(treatment_order), y = Inf, label = paste0("n=", n)),
      hjust = 1.1, vjust = 1.5, size = 8,
      colour = "grey30", fontface = "italic", inherit.aes = FALSE
    ) +
    facet +
    scale_colour_manual(
      values = c("mRNA" = "#E8C84A", "Protein" = "#7BBFDF"),
      name = NULL
    ) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    labs(x = NULL, y = y_lab) +
    theme_bw(base_size = base_size) +
    theme(
      strip.text.x = element_text(size = base_size + 2, face = "bold",
                                  colour = "white"),
      strip.text.y = element_text(size = base_size + 2, face = "italic",
                                  angle = -90),
      strip.background.x = element_rect(fill = "grey40"),
      strip.background.y = element_rect(fill = "white"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.spacing.x = unit(cluster_spacing, "cm"),  # space between clusters
      panel.spacing.y = unit(species_spacing, "cm"),  # space between species
      legend.position = "right",
      legend.text = element_text(size = base_size + 2),
      legend.key.size = unit(0.4, "cm"),
      axis.text.x = element_text(size = base_size),
      axis.text.y = element_text(size = base_size),
      axis.title.y = element_text(size = base_size + 2)
    )
}

#' Plot LogFC Profiles per Cluster and Species
#'
#' @param discordance_df Tibble. Output of build_discordance_clusters().
#' @param output_path    Chr. Output directory.
#' @param species_filter Chr. Single species to plot. NULL = all species.
#'
#' @return wrap_elements gtable object.
plot_discordance_profiles <- function(discordance_df, output_path,
                                      species_filter = NULL) {
  treatment_order <- sort(unique(discordance_df$Treatment[discordance_df$Treatment != "Ctrl"]))

  filter_species <- function(df) {
    if (!is.null(species_filter)) filter(df, Species == species_filter) else df
  }

  df <- discordance_df %>%
    filter(Treatment != "Ctrl") %>%
    filter_species() %>%
    mutate(
      Treatment = factor(Treatment, levels = treatment_order),
      cluster = factor(cluster, levels = names(CLUSTER_COLOURS))
    ) %>%
    pivot_longer(cols = c(rna_log2FC, prot_log2FC),
                 names_to = "measure", values_to = "log2FC") %>%
    mutate(measure = factor(ifelse(measure == "rna_log2FC", "mRNA", "Protein"),
                            levels = c("mRNA", "Protein")))

  count_df <- discordance_df %>%
    filter(Treatment != "Ctrl") %>%
    filter_species() %>%
    distinct(Species, Protein_id, cluster) %>%
    count(Species, cluster, name = "n") %>%
    mutate(cluster = factor(cluster, levels = names(CLUSTER_COLOURS)))

  species_levels <- sort(unique(df$Species))
  species_labeller <- as_labeller(setNames(
    sapply(species_levels, function(s) {
      display <- PLOT_SPECIES_NAMES[s]
      ifelse(is.na(display), gsub("_", " ", s), display)
    }),
    species_levels
  ))

  p <- .build_profile_plot(df, count_df, "log2FC", "Log2 FC",
                           treatment_order, species_filter, species_labeller)

  g <- colour_cluster_strips(ggplot_gtable(ggplot_build(p)))
  wrap_elements(g)
}

#' Plot Z-Score Profiles per Cluster and Species
#'
#' @param discordance_df Tibble. Output of build_discordance_clusters_zscore().
#' @param output_path    Chr. Output directory.
#' @param species_filter Chr. Single species to plot. NULL = all species.
#' @param z_threshold    Num. Z-score threshold shown as shaded band. Default 1.96.
#'
#' @return wrap_elements gtable object.
plot_discordance_profiles <- function(
  discordance_df, output_path,
  species_filter = NULL,
  z_threshold = 1.96,
  cluster_spacing = 0.2,
  species_spacing = 0.5
) {
  treatment_order <- sort(unique(
    discordance_df$Treatment[discordance_df$Treatment != "Ctrl"]
  ))

  filter_species <- function(df) {
    if (!is.null(species_filter)) filter(df, Species == species_filter) else df
  }

  df <- discordance_df %>%
    filter(Treatment != "Ctrl") %>%
    filter_species() %>%
    mutate(
      Treatment = factor(Treatment, levels = treatment_order),
      cluster = factor(cluster, levels = names(CLUSTER_COLOURS))
    ) %>%
    pivot_longer(
      cols = c(rna_z, prot_z),
      names_to = "measure",
      values_to = "zscore"
    ) %>%
    mutate(
      measure = factor(
        ifelse(measure == "rna_z", "mRNA", "Protein"),
        levels = c("mRNA", "Protein")
      )
    )

  count_df <- discordance_df %>%
    filter(Treatment != "Ctrl") %>%
    filter_species() %>%
    distinct(Species, Protein_id, cluster) %>%
    count(Species, cluster, name = "n") %>%
    mutate(cluster = factor(cluster, levels = names(CLUSTER_COLOURS)))

  species_levels <- sort(unique(df$Species))
  species_labeller <- as_labeller(setNames(
    sapply(species_levels, function(s) {
      display <- PLOT_SPECIES_NAMES[s]
      ifelse(is.na(display), gsub("_", " ", s), display)
    }),
    species_levels
  ))

  p <- .build_profile_plot(
    df, count_df, "zscore", "Z-score",
    treatment_order, species_filter, species_labeller,
    cluster_spacing = cluster_spacing,
    species_spacing = species_spacing,
    threshold_band = z_threshold
  )

  g <- colour_cluster_strips(ggplot_gtable(ggplot_build(p)))
  wrap_elements(g)
}

# ============================================================================
# ENRICHMENT
# ============================================================================

#' Load Enrichment Annotation Files for One Species
#'
#' @param species Chr. Species name.
#' @return List: term2gene, term2name, universe.
load_enrichment_data <- function(species) {
  abbv <- switch(species,
                 Salmonella_enterica = "SE",
                 Staphylococcus_aureus = "SA",
                 Yersinia_pseudotuberculosis = "YP",
                 stop("No enrichment files configured for: ", species)
  )

  list(
    term2gene = read_excel(file.path(ENRICHMENT, paste0("TERM2GENE_", abbv, ".xlsx")),
                           col_names = TRUE),
    term2name = read_excel(file.path(ENRICHMENT, paste0("TERM2NAME_", abbv, ".xlsx")),
                           col_names = TRUE),
    universe = read_csv(file.path(ENRICHMENT, paste0("universe_", abbv, ".csv")),
                        col_names = TRUE, show_col_types = FALSE) %>% pull(1)
  )
}

#' Run clusterProfiler enricher for One Gene Set
#'
#' @param gene_ids       Chr vector. Protein IDs to test.
#' @param term2gene      Tibble. TERM, GENE columns.
#' @param term2name      Tibble. TERM, ONTOLOGY columns.
#' @param universe       Chr vector. Background gene set.
#' @param pval_cutoff    Num. p-value cutoff. Default 0.05.
#' @param qval_cutoff    Num. q-value cutoff. Default 0.2.
#' @param min_gs_size    Int. Minimum gene set size. Default 5.
#' @param max_gs_size    Int. Maximum gene set size. Default 500.
#'
#' @return enrichResult object or NULL.
run_enricher <- function(
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
    enricher(
      gene = gene_ids,
      TERM2GENE = term2gene,
      TERM2NAME = term2name,
      universe = universe,
      pvalueCutoff = pval_cutoff,
      qvalueCutoff = qval_cutoff,
      minGSSize = min_gs_size,
      maxGSSize = max_gs_size
    ),
    error = function(e) { message("    enricher error: ", e$message); NULL }
  )

  if (is.null(result) || nrow(result@result) == 0) {
    message("    No significant terms found.")
    return(NULL)
  }
  result
}

#' Run Enrichment Analysis per Cluster and Species
#'
#' @param discordance_df Tibble. Output of build_discordance_clusters().
#' @param pval_cutoff    Num. p-value cutoff. Default 0.05.
#' @param qval_cutoff    Num. q-value cutoff. Default 0.2.
#' @param min_gs_size    Int. Minimum gene set size. Default 5.
#' @param max_gs_size    Int. Maximum gene set size. Default 500.
#' @param output_path    Chr. Directory for saving TSV results.
#'
#' @return Nested list: results[[species]][[cluster]]
run_enrichment_all <- function(
  discordance_df,
  pval_cutoff = 0.05,
  qval_cutoff = 0.2,
  min_gs_size = 5,
  max_gs_size = 500,
  output_path = NULL
) {
  species_list <- unique(discordance_df$Species)
  clusters <- sort(unique(discordance_df$cluster))

  map(species_list, function(species) {
    message("\nRunning enrichment for: ", species)
    annot <- load_enrichment_data(species)

    map(clusters, function(cl) {
      message("  Cluster: ", cl)

      gene_ids <- discordance_df %>%
        filter(Species == species, cluster == cl) %>%
        distinct(Protein_id) %>%
        pull(Protein_id)

      message("    Genes: ", length(gene_ids))

      result <- run_enricher(
        gene_ids = gene_ids,
        term2gene = annot$term2gene,
        term2name = annot$term2name,
        universe = annot$universe,
        pval_cutoff = pval_cutoff,
        qval_cutoff = qval_cutoff,
        min_gs_size = min_gs_size,
        max_gs_size = max_gs_size
      )

      if (!is.null(result) && !is.null(output_path)) {
        out_file <- file.path(output_path, species,
                              paste0(species, "_", cl, "_enrichment.tsv"))
        write.table(as.data.frame(result), out_file,
                    sep = "\t", row.names = FALSE, quote = FALSE)
        message("    Saved: ", out_file)
      }
      result
    }) %>% setNames(clusters)
  }) %>% setNames(species_list)
}

#' Plot Enrichment Dot Plot — All Clusters and Species
#'
#' @param enrichment_results Nested list. Output of run_enrichment_all().
#' @param fdr_threshold      Num. FDR threshold. Default 0.05.
#' @param top_n              Int. Max terms per cluster. NULL = all.
#' @param output_path        Chr. Output directory.
#'
#' @return wrap_elements gtable object invisibly.
plot_enrichment_dotplot <- function(
  enrichment_results,
  fdr_threshold = 0.05,
  top_n = NULL,
  output_path = NULL,
  cluster_spacing = 0.2
) {
  plot_df <- map_dfr(names(enrichment_results), function(species) {
    map_dfr(names(enrichment_results[[species]]), function(cl) {
      res <- enrichment_results[[species]][[cl]]
      if (is.null(res)) return(NULL)
      as.data.frame(res) %>%
        as_tibble() %>%
        mutate(Species = species, cluster = cl)
    })
  }) %>%
    mutate(
      neg_log10_fdr = -log10(p.adjust),
      cluster = factor(cluster, levels = names(CLUSTER_COLOURS)),
      Species = factor(Species, levels = sort(unique(Species)))
    )

  if (nrow(plot_df) == 0) {
    warning("No significant enrichment results to plot.")
    return(invisible(NULL))
  }

  if (!is.null(top_n)) {
    plot_df <- plot_df %>%
      group_by(cluster, Description) %>%
      summarise(min_fdr = min(p.adjust), .groups = "drop") %>%
      group_by(cluster) %>%
      slice_min(min_fdr, n = top_n) %>%
      dplyr::select(cluster, Description) %>%
      inner_join(plot_df, by = c("cluster", "Description"))
  }

  sig_terms <- plot_df %>%
    filter(p.adjust < fdr_threshold) %>%
    distinct(Description)
  plot_df <- plot_df %>% semi_join(sig_terms, by = "Description")

  term_order <- sort(unique(plot_df$Description))
  plot_df <- plot_df %>%
    mutate(Description = factor(Description, levels = rev(term_order)))

  species_levels <- levels(plot_df$Species)
  species_labels <- setNames(
    lapply(species_levels, format_species_title, linebreak = TRUE),
    species_levels
  )

  p <- ggplot(plot_df, aes(x = Species, y = Description,
                           size = Count, colour = neg_log10_fdr)) +
    geom_point(alpha = 0.9) +
    geom_hline(yintercept = seq(0.5, length(term_order) + 0.5, 1),
               colour = "grey92", linewidth = 0.3) +
    facet_wrap(~cluster, nrow = 1, labeller = as_labeller(CLUSTER_LABELS)) +
    scale_colour_gradient(low = "lightsalmon", high = "darkred",
                          name = "-log10(FDR)",
                          guide = guide_colourbar(barwidth = 0.8, barheight = 4, ticks = FALSE)) +
    scale_size_continuous(name = "Count", range = c(1, 8),
                          breaks = c(15, 30, 45, 60, 75, 90)) +
    scale_x_discrete(labels = species_labels) +
    labs(x = NULL, y = NULL) +
    theme_publication(base_size = 32, legend_position = "right") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "grey40"),
      panel.spacing.x = unit(cluster_spacing, "cm")
    )

  g <- colour_cluster_strips(ggplot_gtable(ggplot_build(p)))
  final <- wrap_elements(g)

  if (!is.null(output_path)) {
    n_terms <- n_distinct(plot_df$Description)
    n_clusters <- n_distinct(plot_df$cluster)
    out_file <- file.path(output_path, "enrichment_dotplot")
    save_plot(final,
              filepath = file.path(output_path, "enrichment_dotplot"),
              width = max(24, n_clusters * 3),
              height = max(8, n_terms * 0.35 + 2))
    message("Enrichment dotplot saved: ",
            file.path(output_path, "enrichment_dotplot"), MSG_FIG_FORMAT)
  }

  invisible(final)
}

# Filter all three dataframes to one treatment (and optionally one species)
filter_treatment <- function(treat, species = NULL) {
  list(
    discordance = discordance_df %>%
      filter(Treatment == treat,
             if (!is.null(species)) Species == species else TRUE),
    prot = zscore_protQ_df %>%
      filter(Treatment == treat,
             if (!is.null(species)) Species == species else TRUE),
    rna = zscore_tpm_df %>%
      filter(Treatment == treat,
             if (!is.null(species)) Species == species else TRUE)
  )
}

# Build one combined-species plot for one treatment (no legend, title = treat)
make_treatment_row <- function(treat, show_legend = FALSE) {
  d <- filter_treatment(treat)
  plots <- plot_discordance_scatter_zscore(
    d$discordance, d$prot, d$rna, SUPP_OUTPUT_PATH
  )
  wrap_plots(plots, nrow = 1) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = treat,
      theme = theme(plot.title = element_text(face = "bold", size = 10))
    ) &
    theme(legend.position = if (show_legend) "right" else "none")
}

# Extract cluster-fill legend only (no starshape) from a plot
extract_cluster_legend <- function(treat) {
  d <- filter_treatment(treat)
  p <- plot_discordance_scatter_zscore(
    d$discordance, d$prot, d$rna, SUPP_OUTPUT_PATH
  )[[1]] +
    theme(legend.position = "right") +
    guides(
      starshape = "none",   # remove shape legend
      fill = guide_legend(
        override.aes = list(starshape = 15, colour = "white",
                            size = 3, alpha = 1)
      )
    )
  cowplot::get_legend(p) %>% wrap_elements()
}