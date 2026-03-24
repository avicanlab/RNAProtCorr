# ============================================================================
# Protein quantification from TMT MS data
# ============================================================================
#
# Description:
#   This script processes TMT (Tandem Mass Tag) proteomics data. It:
#     1. Filters TMT results per species and replicate
#     2. Performs in silico trypsin digestion from FASTA files
#     3. Computes relative intensity, iBAQ, and riBAQ quantifications
#
# Requirements:
#   install.packages("tidyverse")
#   install.packages("readxl")
#   if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#   BiocManager::install("Biostrings")
#   BiocManager::install("cleaver")
# ============================================================================

library(Biostrings)
library(cleaver)
library(purrr)
library(ggdist)       # stat_halfeye / raincloud helpers

# ============================================================================
# CONFIGURATION
# ============================================================================

# Column names in TMT results files
col_protein <- "Protein"
col_protein_id <- "Protein ID"
col_entry_name <- "Entry Name"
col_gene <- "Gene"
col_organism <- "Organism"
col_unique_peptides <- "Unique Peptides"
col_razor_peptides <- "Razor Peptides"
col_total_intensity <- "Total Intensity"

keep_cols <- c(
  col_protein,
  col_protein_id,
  col_entry_name,
  col_gene,
  col_organism,
  col_unique_peptides,
  col_razor_peptides,
  col_total_intensity
)

# ============================================================================
# IN SILICO TRYPSIN DIGEST FROM FASTA
# ============================================================================
#' Filter peptides to those within the configured length range
#'
#' @param peps A tibble with at least a `peptide` column
#' @param min_peptide_len Minimum peptide length (default 6)
#' @param max_peptide_len Maximum peptide length (default 30)
#' @return A filtered tibble
filter_peptides <- function(peps, min_peptide_len = 6, max_peptide_len = 30) {
  lengths <- nchar(peps$peptide)
  peps[lengths >= min_peptide_len & lengths <= max_peptide_len,]
}

#' Run in silico trypsin digestion on one or more FASTA files
#'
#' Digests proteins using a single enzyme configuration, filters peptides by
#' length, and writes summary tables to output_path.
#'
#' @param fasta_paths      A character vector of paths to FASTA files
#' @param output_path      Path to the output directory (must end with "/")
#' @param enzym            Enzyme to use for digestion. One of:
#'   - "trypsin"      : cleaves after K/R, not before P (standard)
#'   - "trypsin-high" : cleaves after K/R including before P (high specificity)
#'   - "trypsin-low"  : cleaves after K/R, more permissive rules (low specificity)
#' @param missed_cleavages Integer or integer vector of allowed missed cleavages
#'   (e.g. 0 for none, 2 for up to 0:2). Default: 0.
#' @param min_peptide_len  Minimum peptide length to retain. Default: 6.
#' @param max_peptide_len  Maximum peptide length to retain. Default: 30.
#'
#' @return A tibble with Species, Protein_id and observable peptide count
process_tryptic <- function(fasta_paths, output_path,
                            enzym = "trypsin",
                            missed_cleavages = 0,
                            min_peptide_len = 6,
                            max_peptide_len = 30) {

  all_fasta_ids <- c()
  all_raw_digests <- tibble()
  all_filtered_digests <- tibble()

  for (fasta_path in fasta_paths) {
    message("\nProcessing:", fasta_path, "\n")

    species <- tools::file_path_sans_ext(
      str_extract(basename(fasta_path), "^[A-Za-z]+_[A-Za-z]+")
    )
    fasta <- readAAStringSet(fasta_path)
    fasta_ids <- sub(" .*", "", names(fasta))

    # Warn on duplicate protein IDs across FASTA files (not expected)
    duplicated_ids <- intersect(fasta_ids, all_fasta_ids)
    if (length(duplicated_ids) > 0) {
      warning(
        length(duplicated_ids), " protein ID(s) in ", basename(fasta_path),
        " already present in a previous FASTA file.\n",
        "Duplicated IDs: ", paste(head(duplicated_ids, 5), collapse = ", "),
        if (length(duplicated_ids) > 5) " ..." else ""
      )
    }

    all_fasta_ids <- c(all_fasta_ids, fasta_ids)

    peptides_per_protein <- cleave(fasta, enzym = enzym, missedCleavages = missed_cleavages)
    raw_digests <- map_dfr(names(peptides_per_protein), function(protein_name) {
      tibble(
        Species = species,
        Protein_id = sub(" .*", "", protein_name),
        peptide = as.character(peptides_per_protein[[protein_name]])
      )
    })

    filtered_digests <- filter_peptides(raw_digests, min_peptide_len, max_peptide_len)

    all_raw_digests <- bind_rows(all_raw_digests, raw_digests)
    all_filtered_digests <- bind_rows(all_filtered_digests, filtered_digests)
  }

  # Peptide filtering comparison table — one row per species
  comparison_table <- all_raw_digests %>%
    count(Species, name = "n_before") %>%
    left_join(
      all_filtered_digests %>% count(Species, name = "n_after"),
      by = "Species"
    ) %>%
    mutate(
      n_removed = n_before - n_after,
      pct_kept = round(100 * n_after / n_before, 1)
    )

  message("\nPeptide filtering summary (length ", min_peptide_len, "-", max_peptide_len, " aa):\n")
  output_file <- file.path(output_path, "peptide_filter_comparison.tsv")
  write.table(
    comparison_table,
    output_file,
    sep = "\t", row.names = FALSE, quote = FALSE
  )

  # Observable peptide counts per protein
  observable_peptides <- all_filtered_digests %>%
    count(Species, Protein_id, name = "n_peptides")

  # Per-species summary statistics
  summary_table <- observable_peptides %>%
    group_by(Species) %>%
    summarise(
      Min = min(n_peptides),
      Q1 = quantile(n_peptides, 0.25),
      Median = median(n_peptides),
      Mean = round(mean(n_peptides), 2),
      Q3 = quantile(n_peptides, 0.75),
      Max = max(n_peptides),
      .groups = "drop"
    )

  message("\nSummary of observable peptide counts after filtering:\n")
  output_file <- file.path(output_path, "peptide_summary.tsv")
  write.table(
    summary_table,
    output_file,
    sep = "\t", row.names = FALSE, quote = FALSE
  )

  return(list(all_raw_digests = all_raw_digests, observable_peptides = observable_peptides))
}

# ── Shared theme ──────────────────────────────────────────────────────────────
theme_publication <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = "sans"),
      axis.title = element_text(size = base_size, colour = "#2c2c2a"),
      axis.text = element_text(size = base_size - 1, colour = "#444441"),
      axis.line = element_line(linewidth = 0.4, colour = "#888780"),
      axis.ticks = element_line(linewidth = 0.3, colour = "#888780"),
      panel.grid.major.y = element_line(colour = "#d3d1c7", linewidth = 0.3,
                                        linetype = "dashed"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.title = element_text(size = base_size + 1, face = "bold",
                                colour = "#2c2c2a", margin = margin(b = 4)),
      plot.subtitle = element_text(size = base_size - 1, colour = "#5f5e5a",
                                   margin = margin(b = 8)),
      plot.caption = element_text(size = base_size - 2, colour = "#888780",
                                  hjust = 0, margin = margin(t = 6)),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold",
                                colour = "#2c2c2a"),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(12, 16, 12, 12)
    )
}

# Colour palette — one tint per species (expand as needed)
species_colours <- c(
  "#1D9E75",   # teal-400
  "#7F77DD",   # purple-400
  "#D85A30",   # coral-400
  "#378ADD",   # blue-400
  "#639922",   # green-400
  "#BA7517"    # amber-400
)

# Treatment colour palette — cycles if more than 8 conditions
treatment_colours <- c(
  "#1D9E75",  # teal
  "#7F77DD",  # purple
  "#D85A30",  # coral
  "#378ADD",  # blue
  "#BA7517",  # amber
  "#639922",  # green
  "#D4537E",  # pink
  "#888780"   # gray
)
# =============================================================================
# HELPER: derive protein-length data from the FASTA objects already read by
# process_tryptic.  Call this BEFORE process_tryptic so you still have `fasta`.
# Alternatively pass `fasta_lengths_df` in directly if already available.
# -----------------------------------------------------------------------------
# Expected input:  a tibble with columns  Species  |  protein_length
#   e.g. built like:
#     fasta_lengths_df <- bind_rows(lapply(fasta_files, function(fp) {
#       sp  <- tools::file_path_sans_ext(
#                stringr::str_extract(basename(fp), "^[A-Za-z]+_[A-Za-z]+"))
#       fa  <- Biostrings::readAAStringSet(fp)
#       tibble::tibble(Species = sp, protein_length = Biostrings::width(fa))
#     }))
# =============================================================================


# =============================================================================
# PLOT 1 — Protein length diversity (raincloud)
# =============================================================================
plot_protein_length <- function(fasta_lengths_df,
                                log_scale = TRUE,
                                point_alpha = 0.25,
                                point_size = 0.6,
                                jitter_width = 0.08,
                                title = NULL) {

  species_lvls <- sort(unique(fasta_lengths_df$Species))
  n_sp <- length(species_lvls)
  pal <- setNames(species_colours[seq_len(n_sp)], species_lvls)

  df <- fasta_lengths_df %>%
    filter(!is.na(protein_length), protein_length > 0) %>%
    mutate(Species = factor(Species, levels = species_lvls))

  # Median labels for annotation
  medians <- df %>%
    group_by(Species) %>%
    summarise(med = median(protein_length), .groups = "drop")

  p <- ggplot(df, aes(x = Species, y = protein_length, fill = Species,
                      colour = Species)) +

    # ── Half-violin (right side) ──
    ggdist::stat_halfeye(
      aes(fill = Species),
      adjust = 0.8,
      width = 0.5,
      .width = 0,          # suppress credible-interval slab
      point_colour = NA,
      alpha = 0.55,
      justification = -0.25,
      normalize = "groups"
    ) +

    # ── Raw points (left side jitter) ──
    geom_jitter(
      width = jitter_width,
      height = 0,
      size = point_size,
      alpha = point_alpha
    ) +

    # ── Boxplot overlay ──
    geom_boxplot(
      aes(fill = Species),
      width = 0.12,
      outlier.shape = NA,
      alpha = 0.7,
      colour = "white",
      linewidth = 0.5
    ) +

    # ── Median text annotation ──
    geom_text(
      data = medians,
      aes(x = Species, y = med, label = scales::comma(round(med))),
      vjust = 0.8,
      hjust = 1.3,
      size = 3.5,
      face = "bold",
      colour = "white",
      inherit.aes = FALSE
    ) +

    scale_fill_manual(values = pal) +
    scale_colour_manual(values = pal) +
    scale_y_continuous(
      labels = scales::comma,
      trans = if (log_scale) "log10" else "identity",
      name = if (log_scale) "Protein length (aa, log\u2081\u2080)"
      else           "Protein length (aa)"
    ) +
    coord_flip() +
    labs(
      title = title,
      x = NULL,
    ) +
    theme_publication()

  p
}


# =============================================================================
# PLOT 2 — Peptide size distribution after in-silico digestion
# =============================================================================
plot_peptide_length <- function(all_raw_digests,
                                min_peptide_len = 6,
                                max_peptide_len = 30) {

  df <- all_raw_digests %>%
    mutate(
      pep_len = nchar(peptide),
      filtered = pep_len >= min_peptide_len & pep_len <= max_peptide_len
    )

  # Counts for caption
  n_total <- nrow(df)
  n_kept <- sum(df$filtered)
  pct_kept <- round(100 * n_kept / n_total, 1)

  # Truncate x-axis display slightly beyond filter window for readability
  x_max_disp <- max_peptide_len + 15

  # Per-length count table (all peptides, capped at x_max_disp)
  counts <- df %>%
    filter(pep_len <= x_max_disp) %>%
    count(pep_len, filtered)

  # Shaded filter-window ribbon (full height — will be behind bars)
  window_df <- data.frame(
    xmin = min_peptide_len - 0.5,
    xmax = max_peptide_len + 0.5
  )

  p <- ggplot(counts, aes(x = pep_len, y = n, fill = filtered)) +

    # ── Filter window band ──
    geom_rect(
      data = window_df,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = "#E1F5EE",   # teal-50
      alpha = 0.55
    ) +

    # ── Histogram bars ──
    geom_col(
      aes(fill = filtered),
      width = 0.85,
      alpha = 0.85
    ) +

    # ── Filter window border lines ──
    geom_vline(xintercept = min_peptide_len - 0.5,
               colour = "#0F6E56", linewidth = 0.5, linetype = "dashed") +
    geom_vline(xintercept = max_peptide_len + 0.5,
               colour = "#0F6E56", linewidth = 0.5, linetype = "dashed") +

    # ── Window label ──
    annotate(
      "text",
      x = (min_peptide_len + max_peptide_len) / 2,
      y = Inf,
      label = sprintf("Retained window\n%d\u2013%d aa  (%s%%)",
                      min_peptide_len, max_peptide_len, pct_kept),
      vjust = 1.25,
      size = 3,
      colour = "#0F6E56",
      family = "sans"
    ) +

    scale_fill_manual(
      values = c("TRUE" = "#1D9E75",   # teal-400  (kept)
                 "FALSE" = "#B4B2A9"),  # gray-200  (discarded)
      labels = c("TRUE" = "Retained",
                 "FALSE" = "Filtered out"),
      name = NULL
    ) +

    scale_x_continuous(
      breaks = c(seq(0, x_max_disp, by = 5)),
      limits = c(0, x_max_disp + 1),
      name = "Peptide length (aa)"
    ) +

    scale_y_continuous(
      labels = scales::comma,
      expand = expansion(mult = c(0, 0.12)),
      name = "Number of peptides"
    ) +

    labs(
      title = "Tryptic peptide length distribution",
      subtitle = sprintf(
        "In silico digest  \u2022  total = %s peptides  \u2022  %s%% retained after length filter",
        scales::comma(n_total), pct_kept
      )
    ) +

    theme_publication() +
    theme(
      legend.position = c(0.82, 0.88),
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.55, "lines"),
      legend.background = element_rect(fill = "white", colour = NA)
    )
  p
}


# =============================================================================
# COMBINED figure (stacked, suitable for a methods section)
# =============================================================================
make_ibaq_figure <- function(fasta_lengths_df,
                             all_raw_digests,
                             min_peptide_len = 6,
                             max_peptide_len = 30,
                             output_file = NULL,
                             width = 9,
                             height = 10) {

  p1 <- plot_protein_length(fasta_lengths_df)
  p2 <- plot_peptide_length(all_raw_digests,
                            min_peptide_len = min_peptide_len,
                            max_peptide_len = max_peptide_len)

  fig <- p1 / p2 +
    patchwork::plot_annotation(
      title = "Justification for iBAQ normalisation",
      subtitle = paste0(
        "Protein length heterogeneity (top) directly determines the number of ",
        "observable tryptic peptides (bottom),\n",
        "biasing raw spectral counts toward longer proteins. ",
        "iBAQ divides intensity by the observable peptide count to remove this bias."
      ),
      caption = "In silico trypsin digestion  |  0 missed cleavages  |  observable window 6\u201330 aa",
      theme = theme_publication() +
        theme(
          plot.title = element_text(size = 13, face = "bold"),
          plot.subtitle = element_text(size = 10, colour = "#5f5e5a",
                                       lineheight = 1.4)
        )
    ) &
    theme(plot.background = element_rect(fill = "white", colour = NA))

  if (!is.null(output_file)) {
    ggsave(output_file, fig, width = width, height = height, dpi = 300,
           bg = "white")
    message("Saved: ", output_file)
  }

  fig
}


# =============================================================================
# USAGE EXAMPLE
# Replace the mock data below with your real objects after running
# process_tryptic().
# =============================================================================
if (FALSE) {

  # ── 1. Build protein-length data frame ──────────────────────────────────────
  fasta_lengths_df <- bind_rows(lapply(fasta_files, function(fp) {
    sp <- tools::file_path_sans_ext(
      stringr::str_extract(basename(fp), "^[A-Za-z]+_[A-Za-z]+"))
    fa <- Biostrings::readAAStringSet(fp)
    tibble::tibble(Species = sp, protein_length = Biostrings::width(fa))
  }))

  # ── 2. Run digestion (returns observable peptide counts, but we also need   ──
  #       the raw digest tibble — uncomment the internal all_raw_digests export ──
  #       in process_tryptic, or rebuild it here) ───────────────────────────────
  protein2tryptic <- process_tryptic(
    fasta_paths = fasta_files,
    output_path = OUTPUT_PATH,
    min_peptide_len = MIN_PEPTIDE_LEN,
    max_peptide_len = MAX_PEPTIDE_LEN,
    missed_cleavages = 0
  )
  # If all_raw_digests isn't exported, rebuild quickly:
  # all_raw_digests <- bind_rows(lapply(fasta_files, function(fp) {
  #   sp  <- tools::file_path_sans_ext(
  #              stringr::str_extract(basename(fp), "^[A-Za-z]+_[A-Za-z]+"))
  #   fa  <- Biostrings::readAAStringSet(fp)
  #   peps <- cleave(fa, enzym = "trypsin", missedCleavages = 0)
  #   purrr::map_dfr(names(peps), ~tibble::tibble(
  #     Species    = sp,
  #     Protein_id = sub(" .*", "", .x),
  #     peptide    = as.character(peps[[.x]])
  #   ))
  # }))

  # ── 3. Render and save ───────────────────────────────────────────────────────
  fig <- make_ibaq_figure(
    fasta_lengths_df = fasta_lengths_df,
    all_raw_digests = all_raw_digests,
    min_peptide_len = MIN_PEPTIDE_LEN,
    max_peptide_len = MAX_PEPTIDE_LEN,
    output_file = file.path(OUTPUT_PATH, "ibaq_justification.pdf")
  )

  print(fig)
}
# ============================================================================
# QUANTIFICATION
# ============================================================================

#' Compute relative intensity per protein per sample
#'
#' Divides total intensity by the number of razor peptides, then normalises
#' within each Species/Replicate group. Adds log2, log10, and parts-per-billion
#' columns.
#'
#' @param tmt_data A unified TMT tibble with Species and Replicate columns
#' @return The input tibble with additional TopN quantification columns
Intensity_quantification <- function(tmt_data) {
  message("Computing intensity quantification for each sample...\n")

  # Compute per-row channel sum (across all channel_ columns)
  channel_cols <- grep("^channel_", colnames(tmt_data), value = TRUE)

  tmt_data %>%
    mutate(
      channel_sum = rowSums(across(all_of(channel_cols)), na.rm = TRUE),
    ) %>%
    pivot_longer(
      cols = all_of(channel_cols),
      names_to = "channel_col",
      values_to = "channel_intensity"
    ) %>%
    # Keep only the channel matching this row's treatment
    filter(channel_col == paste0("channel_", TMT_label)) %>%
    mutate(
      Intensity = .data[[col_total_intensity]] * channel_intensity / channel_sum,
    ) %>%
    filter(if_all(Intensity, ~. > 0 & !is.na(.))) %>%
    group_by(Species, Protein_id, Treatment) %>%
    mutate(
      Intensity_log2 = log2(Intensity),
      Intensity_meanlog2 = mean(Intensity_log2, na.rm = TRUE),
      Intensity_log2ratio = Intensity_log2 - Intensity_meanlog2,
      Intensity_norm = Intensity / sum(Intensity, na.rm = TRUE),
      Intensity_log = 10 + log10(Intensity_norm),
      Intensity_PpB = Intensity_norm * 1e9
    ) %>%
    ungroup() %>%
    select(-channel_sum, -channel_intensity)
}

#' Compute relative intensity per protein per sample
#'
#' Divides total intensity by the number of razor peptides, then normalises
#' within each Species/Replicate group. Adds log2, log10, and parts-per-billion
#' columns.
#'
#' @param tmt_data A unified TMT tibble with Species and Replicate columns
#' @return The input tibble with additional TopN quantification columns
RI_quantification <- function(tmt_data) {
  message("Computing relative intensity quantification for each sample...\n")

  # Compute per-row channel sum (across all channel_ columns)
  channel_cols <- grep("^channel_", colnames(tmt_data), value = TRUE)

  tmt_data %>%
    mutate(
      channel_sum = rowSums(across(all_of(channel_cols)), na.rm = TRUE),
      RI_total = (.data[[col_total_intensity]] / .data[[col_razor_peptides]])
    ) %>%
    pivot_longer(
      cols = all_of(channel_cols),
      names_to = "channel_col",
      values_to = "channel_intensity"
    ) %>%
    # Keep only the channel matching this row's treatment
    filter(channel_col == paste0("channel_", TMT_label)) %>%
    mutate(
      RI = RI_total * (channel_intensity / channel_sum)
    ) %>%
    filter(if_all(RI, ~. > 0 & !is.na(.))) %>%
    group_by(Species, Protein_id, Treatment) %>%
    mutate(
      RI_log2 = log2(RI),
      RI_meanlog2 = mean(RI_log2, na.rm = TRUE),
      RI_log2ratio = RI_log2 - RI_meanlog2,
      RI_norm = RI / sum(RI, na.rm = TRUE),
      RI_log = 10 + log10(RI_norm),
      RI_PpB = RI_norm * 1e9
    ) %>%
    ungroup() %>%
    select(-channel_sum, -channel_intensity)
}

#' Compute iBAQ and riBAQ intensities per protein per sample
#'
#' Joins observable peptide counts and divides total intensity by the number
#' of observable peptides (iBAQ). Then normalises within each Species/Replicate
#' group (riBAQ) and adds log-transformed columns.
#'
#' @param tmt_data          A unified TMT tibble with Species and Replicate columns
#' @param observable_peptides A tibble from process_tryptic() with protein_id and n_* columns
#' @param observable_col    Name of the column to use as denominator (default: "n_tryptic")
#' @return The input tibble joined with iBAQ quantification columns
ibaq_quantification <- function(
  tmt_data,
  observable_peptides,
  observable_col = "n_tryptic",
  output_path = NULL
) {
  message("Computing iBAQ for each sample...\n")

  channel_cols <- grep("^channel_", colnames(tmt_data), value = TRUE)

  tmt_data %>%
    left_join(observable_peptides, by = c("Species", "Protein_id")) %>%
    mutate(
      unmatched = is.na(.data[[observable_col]])
    ) %>%
  {
    unmatched_proteins <- distinct(filter(., unmatched), Species, Protein_id)
    n_unmatched <- nrow(unmatched_proteins)

    if (n_unmatched > 0) {
      warning(n_unmatched, " proteins had no match in FASTA and will have iBAQ = NA")

      if (!is.null(output_path)) {
        unmatched_path <- file.path(output_path, "iBAQ_unmatched_proteins.tsv")
        write.table(
          unmatched_proteins %>% arrange(Species, Protein_id),
          unmatched_path,
          sep = "\t", row.names = FALSE, quote = FALSE
        )
        message("Unmatched proteins written to: ", unmatched_path)
      }
    }
    .
  } %>%
    mutate(
      iBAQ_total = total_intensity / .data[[observable_col]],
      channel_sum = rowSums(across(all_of(channel_cols)), na.rm = TRUE)
    ) %>%
    pivot_longer(
      cols = all_of(channel_cols),
      names_to = "channel_col",
      values_to = "channel_intensity"
    ) %>%
    filter(channel_col == paste0("channel_", TMT_label)) %>%
    mutate(iBAQ = iBAQ_total * (channel_intensity / channel_sum)) %>%
    group_by(Species, Protein_id, Treatment) %>%
    mutate(
      iBAQ_log2 = log2(iBAQ),
      iBAQ_meanlog2 = mean(iBAQ_log2, na.rm = TRUE),
      iBAQ_log2ratio = log2(iBAQ) - iBAQ_meanlog2,
      riBAQ = iBAQ / sum(iBAQ, na.rm = TRUE),
      iBAQ_log = 10 + log10(riBAQ),
      iBAQ_PpB = riBAQ * 1e9
    ) %>%
    ungroup()
}

# =============================================================================
# CORE: one scatter panel for a single replicate pair within one condition
# -----------------------------------------------------------------------------
# pair_df  : two-column tibble  rep_x | rep_y  (log2 iBAQ, one row per protein)
# rep_x/y  : replicate labels (character), used for axis titles
# condition: character label for the strip / subtitle
# colour   : hex fill for points and regression line
# =============================================================================
.one_panel <- function(df_wide, sample_x, sample_y, colour) {

  d <- df_wide %>%
    select(x = all_of(sample_x), y = all_of(sample_y)) %>%
    filter(is.finite(x), is.finite(y))

  if (nrow(d) < 3) {
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5,
                      label = "insufficient data", colour = "#888780") +
             theme_void())
  }

  # ── Pearson correlation ───────────────────────────────────────────────────
  ct    <- cor.test(d$x, d$y, method = "pearson")
  R_val <- round(ct$estimate, 2)
  p_val <- ct$p.value
  p_fmt <- if (p_val < 0.001) "0.00e+00" else sprintf("%.3f", p_val)

  label_txt <- sprintf("Pearson correlation: %.2f\np-value: %s", R_val, p_fmt)

  # ── Axis limits (common, padded) ──────────────────────────────────────────
  lo <- floor(min(c(d$x, d$y)))
  hi <- ceiling(max(c(d$x, d$y)))

  # Annotation position — top-left, clear of points
  ann_x <- lo + 0.03 * (hi - lo)
  ann_y <- hi - 0.03 * (hi - lo)

  ggplot(d, aes(x = x, y = y)) +

    # Identity diagonal
    geom_abline(slope = 1, intercept = 0,
                colour = "#cccccc", linewidth = 0.35, linetype = "solid") +

    # Points
    geom_point(colour = colour, alpha = 0.35, size = 0.55, shape = 16) +

    # Stats label
    annotate("text",
             x = ann_x, y = ann_y,
             label      = label_txt,
             hjust      = 0, vjust = 1,
             size       = 2.4,
             colour     = "#2c2c2a",
             lineheight = 1.35,
             family     = "sans") +

    scale_x_continuous(limits = c(lo, hi),
                       breaks = scales::pretty_breaks(n = 4)) +
    scale_y_continuous(limits = c(lo, hi),
                       breaks = scales::pretty_breaks(n = 4)) +

    coord_fixed() +

    labs(
      title = sprintf("%s vs %s", sample_x, sample_y),
      x     = sample_x,
      y     = sample_y
    ) +
    theme_publication(base_size=10)
}


# =============================================================================
# BUILD: all pairwise panels for one species  →  one PDF
# -----------------------------------------------------------------------------
# df           : full data frame filtered to one species
#                columns: Species, Protein_id, Treatment, Replicate, iBAQ, iBAQ_log2
# species_name : character, used in file name and figure title
# output_dir   : directory to write the PDF into
# colour_map   : named vector treatment → hex colour (auto-generated if NULL)
# page_ncol    : number of pair-panels per row on each PDF page (default = 3,
#                one per replicate pair for 3 replicates)
# width/height : PDF dimensions in inches
# =============================================================================
plot_replicate_correlations <- function(df,
                                        species_name,
                                        sample_col  = NULL,
                                        colour_map  = NULL) {

  # ── Build a full sample label if not already a column ─────────────────────
  if (!is.null(sample_col)) {
    df <- df %>% rename(.sample = all_of(sample_col))
  } else {
    df <- df %>%
      mutate(.sample = paste0(Treatment, "_", Replicate))
  }

  treatments <- sort(unique(df$Treatment))
  n_tr       <- length(treatments)

  # Auto colour map
  if (is.null(colour_map)) {
    pal        <- treatment_colours[((seq_len(n_tr) - 1) %% length(treatment_colours)) + 1]
    colour_map <- setNames(pal, treatments)
  }

  # ── Pivot to wide: one column per sample label ────────────────────────────
  wide <- df %>%
    select(Protein_id, .sample, iBAQ_log2) %>%
    pivot_wider(names_from  = .sample,
                values_from = iBAQ_log2,
                values_fn   = mean)   # collapse duplicates if any

  # ── Build panel grid: treatments (rows) × pairs (cols) ───────────────────
  all_rows <- map(treatments, function(trt) {

    trt_df  <- df %>% filter(Treatment == trt)
    samples <- sort(unique(trt_df$.sample))

    if (length(samples) < 2) {
      message("  Skipping ", trt, ": fewer than 2 replicates.")
      return(NULL)
    }

    pairs <- combn(samples, 2, simplify = FALSE)

    panels <- map(pairs, function(pr) {
      .one_panel(wide, pr[1], pr[2], colour_map[[trt]])
    })

    # Force exactly 3 columns (one per pair for 3 replicates)
    wrap_plots(panels, ncol = 3)
  })

  all_rows <- compact(all_rows)

  if (length(all_rows) == 0) {
    warning("No panels generated for: ", species_name)
    return(invisible(NULL))
  }

  # ── Assemble full figure (treatments stacked) ─────────────────────────────
  fig <- wrap_plots(all_rows, ncol = 2) +
    plot_annotation(
      title   = species_name,
      caption = "log\u2082(iBAQ)  \u2022  Pearson correlation  \u2022  identity line shown",
      theme   = theme(
        plot.title      = element_text(size = 11, face = "bold", colour = "#2c2c2a",
                                       family = "sans"),
        plot.caption    = element_text(size = 7,  colour = "#888780", family = "sans"),
        plot.background = element_rect(fill = "white", colour = NA)
      )
    )

fig
}