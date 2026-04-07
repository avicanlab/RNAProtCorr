# Protein Normalization --------------------------------------------------------
#
# Description:
#   Functions for comparing multiple NormalyzerDE normalization methods across
#   species and quantification types, and for running a single chosen method
#   to produce tidy normalized data for downstream analyses.
#
# Requirements:
#   library(NormalyzerDE)
#   library(dplyr)
#   library(tidyr)
#   library(stringr)
#   library(purrr)
#   library(ggplot2)
#   library(patchwork)

library(NormalyzerDE)

# ── Constants ──────────────────────────────────────────────────────────────────

#' Normalization methods available in NormalyzerDE
NORM_METHODS <- c("Log2", "VSN", "CycLoess", "Quantile", "Median", "Mean", "GI", "RLR")

# ── Internal helpers ───────────────────────────────────────────────────────────

#' Write a Wide-Format Protein Matrix to a Temporary TSV File
#'
#' @description
#' Reshapes long-format data into a wide sample × protein matrix required by
#' NormalyzerDE and writes it to a temporary file. Sample columns are sorted
#' alphabetically for reproducibility.
#'
#' @param data   Tibble. Long-format data with columns: `Protein_id`,
#'   `Replicate`, `Treatment`, and the quantification column.
#' @param protQ  Character. Name of the quantification column.
#' @param tmp_dir Character. Temporary directory for the output file.
#'   Default `tempdir()`.
#'
#' @return Character. Path to the temporary TSV file.
.write_wide_matrix <- function(data, protQ, tmp_dir = tempdir()) {
  data_wide <- data |>
    dplyr::mutate(sample = paste0(Replicate, "_", Treatment)) |>
    tidyr::pivot_wider(
      id_cols     = Protein_id,
      names_from  = sample,
      values_from = dplyr::all_of(protQ)
    ) %>%
    dplyr::select(Protein_id, dplyr::all_of(sort(colnames(dplyr::select(., -Protein_id)))))

  temp_file <- tempfile(
    pattern = paste0("matrix_", protQ, "_"),
    tmpdir  = tmp_dir,
    fileext = ".tsv"
  )
  utils::write.table(data_wide, temp_file, sep = "\t", row.names = FALSE, quote = FALSE)
  temp_file
}

#' Write a NormalyzerDE Design File to a Temporary TSV
#'
#' @description
#' Reads column names from an existing matrix TSV and writes a two-column
#' design file (`sample`, `group`) where `group` is the treatment label
#' (everything after the replicate prefix `R<n>_`).
#'
#' @param data_path Character. Path to a matrix TSV produced by
#'   `.write_wide_matrix()`.
#' @param tmp_dir   Character. Temporary directory. Default `tempdir()`.
#'
#' @return Character. Path to the temporary design TSV file.
.write_design_file <- function(data_path, tmp_dir = tempdir()) {
  sample_cols <- utils::read.table(
    data_path, sep = "\t", header = TRUE,
    check.names = FALSE, nrows = 1
  ) |>
    dplyr::select(-Protein_id) |>
    colnames()

  design_df <- tibble::tibble(
    sample = sample_cols,
    group  = stringr::str_remove(sample_cols, "^R\\d+_")
  )

  temp_file <- tempfile(pattern = "design_", tmpdir = tmp_dir, fileext = ".tsv")
  utils::write.table(design_df, temp_file, sep = "\t", row.names = FALSE, quote = FALSE)
  temp_file
}

#' Read a NormalyzerDE Output File and Return Tidy Long-Format Data
#'
#' @description
#' Reads a `*-normalized.txt` file, removes auxiliary columns (`CV`,
#' `anovaP`), and pivots to long format with `Replicate` and `Treatment`
#' parsed from the sample column.
#'
#' @param norm_file Character. Full path to the normalized file.
#' @param species   Character. Species label to add to the output.
#'
#' @return Tibble with columns: `Protein_id`, `Species`, `Replicate`,
#'   `Treatment`, `value`.
.read_norm_file <- function(norm_file, species) {
  utils::read.table(norm_file, sep = "\t", header = TRUE, check.names = FALSE) |>
    dplyr::select(-dplyr::any_of(c("CV", "anovaP"))) |>
    tidyr::pivot_longer(
      cols      = -Protein_id,
      names_to  = "sample",
      values_to = "value"
    ) |>
    dplyr::mutate(
      Species   = species,
      Replicate = stringr::str_extract(sample, "^R\\d+"),
      Treatment = stringr::str_remove(sample,  "^R\\d+_")
    ) |>
    dplyr::select(-sample) |>
    dplyr::select(Species, Treatment, Replicate, Protein_id, value)
}

#' Compute PCA from a Normalized Matrix
#'
#' @description
#' Reads a normalized matrix file, removes non-numeric auxiliary columns,
#' drops incomplete rows, runs PCA (scaled and centred), and joins the first
#' two principal components with design metadata.
#'
#' @param norm_matrix  Data frame or tibble. Wide normalized matrix with
#'   a `Protein_id` row identifier.
#' @param design_table Character. Path to the design TSV file.
#' @param method_name  Character. Method label added to the plot title.
#'
#' @return Named list:
#'   - `df`   : tibble of PCA scores with sample metadata.
#'   - `plot` : ggplot PCA scatter plot.
.compute_pca <- function(norm_matrix, design_table, method_name) {
  mat <- norm_matrix |>
    dplyr::select(-dplyr::any_of(c("CV", "anovaP"))) |>
    tibble::column_to_rownames("Protein_id") |>
    dplyr::select(dplyr::where(is.numeric)) |>
    as.matrix()
  mat <- mat[stats::complete.cases(mat), ]

  pca     <- stats::prcomp(t(mat), scale. = TRUE, center = TRUE)
  var_exp <- round(summary(pca)$importance[2, 1:2] * 100, 1)

  design_df <- utils::read.table(design_table, header = TRUE, sep = "\t") |>
    tibble::as_tibble() |>
    dplyr::mutate(
      Treatment = group,
      Replicate = stringr::str_remove(sample, "_[a-zA-Z]+$")
    )

  pca_df <- as.data.frame(pca$x[, 1:2]) |>
    tibble::rownames_to_column("sample") |>
    dplyr::left_join(design_df, by = "sample")

  pca_plot <- ggplot2::ggplot(
    pca_df,
    ggplot2::aes(
      x     = PC1,
      y     = PC2,
      color = Treatment,
      shape = factor(Replicate),
      label = sample
    )
  ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(
      title = paste("PCA \u2014", method_name),
      x     = paste0("PC1 (", var_exp[1], "%)"),
      y     = paste0("PC2 (", var_exp[2], "%)"),
      shape = "Replicate"
    ) +
    theme_publication()

  list(df = pca_df, plot = pca_plot)
}

#' Compute, Save, and Return a PCA from a Normalized File
#'
#' @param norm_path   Character. Path to the normalized matrix file.
#' @param design_path Character. Path to the design TSV file.
#' @param output_path Character. Directory for saving the PCA coordinate TSV
#'   and PNG plot.
#' @param method_name Character. Method label used in file names and titles.
#'
#' @return NULL invisibly.
.process_plot_pca <- function(norm_path, design_path, output_path, method_name) {
  norm_matrix <- utils::read.table(norm_path, sep = "\t", header = TRUE, check.names = FALSE)
  pca_res     <- .compute_pca(norm_matrix, design_path, method_name)

  tsv_out <- file.path(output_path, paste0(method_name, "_PCA_coordinates.tsv"))
  utils::write.table(pca_res$df, tsv_out, sep = "\t", row.names = FALSE, quote = FALSE)
  message("PCA coordinates saved to: ", tsv_out)

  png_out <- file.path(output_path, paste0(method_name, ".png"))
  ggplot2::ggsave(filename = png_out, plot = pca_res$plot, width = 6, height = 5)
  message("PCA plot saved to: ", png_out)

  invisible(NULL)
}

# ── Main functions ─────────────────────────────────────────────────────────────

#' Run Normalization Method Comparison Across All Species
#'
#' @description
#' For each species, prepares a wide matrix and design file, runs all
#' normalization methods via `NormalyzerDE::normalyzer()`, and generates PCA
#' plots for each method. Output is written to `<output_path>/<species>/`.
#'
#' @param protQ_data  Tibble. Long-format data with columns: `Species`,
#'   `Replicate`, `Treatment`, `Protein_id`, and the quantification column.
#' @param protQ       Character. Name of the quantification column.
#' @param output_path Character. Root output directory.
#'
#' @return NULL invisibly.
run_normalization_comparison <- function(protQ_data, protQ, output_path) {
  tmp_dir      <- tempdir(check = FALSE)
  species_list <- unique(protQ_data$Species)

  for (species in species_list) {
    message("Normalizing ", species, " data ...")

    sp_data     <- dplyr::filter(protQ_data, Species == species)
    data_path   <- .write_wide_matrix(sp_data, protQ, tmp_dir)
    design_path <- .write_design_file(data_path, tmp_dir)

    sp_dir <- file.path(output_path, species)
    dir.create(sp_dir, showWarnings = FALSE)

    NormalyzerDE::normalyzer(
      jobName          = "Normalization",
      dataPath         = data_path,
      designPath       = design_path,
      outputDir        = sp_dir,
      writeReportAsPngs = TRUE,
      zeroToNA         = TRUE
    )

    pca_dir <- file.path(sp_dir, "Normalization", "PCA")
    dir.create(pca_dir, showWarnings = FALSE)

    for (norm_method in NORM_METHODS) {
      message("  Generating PCA for: ", norm_method)

      norm_file <- file.path(sp_dir, "Normalization", paste0(norm_method, "-normalized.txt"))
      if (!file.exists(norm_file)) {
        warning("Normalized file not found for method ", norm_method, ": ", norm_file)
        next
      }

      .process_plot_pca(
        norm_path   = norm_file,
        design_path = design_path,
        output_path = pca_dir,
        method_name = paste0(protQ, "_", norm_method)
      )
    }

    message("Species ", species, " complete. Outputs saved to: ", sp_dir)
  }

  invisible(NULL)
}

#' Plot PCA Comparison Grid for All Normalization Methods
#'
#' @description
#' Generates a 2-column grid of PCA plots — one panel per normalization
#' method, each faceted by species. Reads the normalized files produced by
#' `run_normalization_comparison()`.
#'
#' @param protQ_data  Tibble. Long-format data with `Species`, `Replicate`,
#'   `Treatment`, `Protein_id`, and the quantification column.
#' @param output_path Character. Root NormalyzerDE output directory (produced
#'   by `run_normalization_comparison()`).
#'
#' @return NULL invisibly.
plot_pca_comparison_grid <- function(protQ_data, output_path) {
  species_list <- unique(protQ_data$Species)

  design_map <- purrr::map(species_list, function(species) {
    protQ_data |>
      dplyr::filter(Species == species) |>
      dplyr::distinct(Replicate, Treatment) |>
      dplyr::mutate(
        sample = paste0(Replicate, "_", Treatment),
        group  = Treatment
      )
  }) |>
    setNames(species_list)

  pca_all <- purrr::map_dfr(NORM_METHODS, function(norm_method) {
    purrr::map_dfr(species_list, function(species) {
      norm_file <- file.path(
        output_path, species, "Normalization",
        paste0(norm_method, "-normalized.txt")
      )
      if (!file.exists(norm_file)) {
        warning("Missing: ", norm_file)
        return(NULL)
      }

      norm_matrix <- utils::read.table(norm_file, sep = "\t", header = TRUE, check.names = FALSE)

      mat <- norm_matrix |>
        dplyr::select(-dplyr::any_of(c("CV", "anovaP"))) |>
        tibble::column_to_rownames("Protein_id") |>
        dplyr::select(dplyr::where(is.numeric)) |>
        as.matrix()
      mat <- mat[stats::complete.cases(mat), ]

      pca     <- stats::prcomp(t(mat), scale. = TRUE, center = TRUE)
      var_exp <- round(summary(pca)$importance[2, 1:2] * 100, 1)

      as.data.frame(pca$x[, 1:2]) |>
        tibble::rownames_to_column("sample") |>
        dplyr::left_join(design_map[[species]], by = "sample") |>
        dplyr::mutate(
          Species      = species,
          method       = norm_method,
          var_exp_PC1  = var_exp[1],
          var_exp_PC2  = var_exp[2]
        )
    })
  })

  # Build one faceted plot per normalization method
  method_plots <- purrr::map(NORM_METHODS, function(norm_method) {
    df <- dplyr::filter(pca_all, method == norm_method)

    # Species facet labels using PLOT_SPECIES_NAMES
    species_labels <- setNames(
      sapply(unique(df$Species), function(s) {
        display <- PLOT_SPECIES_NAMES[s]
        ifelse(is.na(display), gsub("_", " ", s), display)
      }),
      unique(df$Species)
    )

    # Use the first species to derive shared axis labels
    axis_labels <- df |>
      dplyr::distinct(Species, var_exp_PC1, var_exp_PC2) |>
      dplyr::slice(1)

    ggplot2::ggplot(
      df,
      ggplot2::aes(x = PC1, y = PC2, color = Treatment, shape = factor(Replicate))
    ) +
      ggplot2::geom_point(size = 3.5, alpha = 0.9) +
      ggplot2::facet_wrap(
        ~Species, nrow = 1, scales = "free",
        labeller = ggplot2::as_labeller(species_labels)
      ) +
      ggplot2::labs(
        x     = paste0("PC1 (", axis_labels$var_exp_PC1, "%)"),
        y     = paste0("PC2 (", axis_labels$var_exp_PC2, "%)"),
        color = "Treatment",
        shape = "Replicate"
      ) +
      theme_publication(base_size = 50) +
      ggplot2::theme(
        plot.title       = ggplot2::element_text(face = "bold", hjust = 0),
        axis.text        = ggplot2::element_blank(),
        axis.ticks       = ggplot2::element_blank(),
        strip.background = ggplot2::element_rect(fill = "#8FBC8F"),
        strip.text       = ggplot2::element_text(face = "italic"),
        panel.grid.minor = ggplot2::element_blank(),
        legend.position  = "right"
      )
  })

  combined <- patchwork::wrap_plots(method_plots, ncol = 2) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(tag_levels = "A") &
    ggplot2::theme(
      legend.position = "right",
      plot.tag        = ggplot2::element_text(size = 60)
    )

  out_file <- file.path(output_path, "normalization_PCA_comparison")
  save_plot(combined, out_file, width = 24, height = 24, out_format = FIGURE_FORMAT)
  message("PCA comparison grid saved: ", out_file, MSG_FIG_FORMAT)

  invisible(NULL)
}

#' Normalize Protein Data Using a Single NormalyzerDE Method
#'
#' @description
#' Normalizes long-format proteomics data using one method from NormalyzerDE.
#' Runs normalization per species then recombines into a single tidy tibble.
#'
#' @param protQ_data  Tibble. Long-format data with columns: `Species`,
#'   `Replicate`, `Treatment`, `Protein_id`, and the quantification column.
#' @param protQ       Character. Name of the quantification column.
#' @param method      Character. One of `NORM_METHODS`. Default `"VSN"`.
#' @param output_path Character or NULL. If provided, the normalized data is
#'   saved as a TSV to this directory. Default NULL.
#'
#' @return Tibble in long format with columns: `Protein_id`, `Species`,
#'   `Replicate`, `Treatment`, `value` (normalized intensities).
run_normalization <- function(protQ_data, protQ, method = "VSN", output_path = NULL) {
  if (!method %in% NORM_METHODS) {
    stop(
      "Unknown normalization method: '", method, "'.\n",
      "Valid options: ", paste(NORM_METHODS, collapse = ", ")
    )
  }
  message("Running ", method, " normalization for ", protQ, "...")

  tmp_dir <- tempdir(check = FALSE)

  normalized_data <- purrr::map_dfr(unique(protQ_data$Species), function(species) {
    message("  Processing species: ", species)

    sp_data     <- dplyr::filter(protQ_data, Species == species)
    data_path   <- .write_wide_matrix(sp_data, protQ, tmp_dir)
    design_path <- .write_design_file(data_path, tmp_dir)

    job_name <- paste0("norm_", species)
    NormalyzerDE::normalyzer(
      jobName      = job_name,
      dataPath     = data_path,
      designPath   = design_path,
      outputDir    = tmp_dir,
      skipAnalysis = TRUE,
      zeroToNA     = TRUE,
      quiet        = TRUE
    )

    norm_file <- file.path(tmp_dir, job_name, paste0(method, "-normalized.txt"))
    if (!file.exists(norm_file)) {
      stop(
        "Normalized file not found: ", norm_file,
        "\nCheck that method '", method, "' is valid."
      )
    }

    .read_norm_file(norm_file, species)
  })

  if (!is.null(output_path)) {
    out_file <- file.path(output_path, paste0(method, "_", protQ, "_normalized.tsv"))
    utils::write.table(normalized_data, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
    message(method, " normalized data saved to: ", out_file)
  }

  normalized_data
}
