# Utilities function for all analyses -----------------------------------------
#
# Description:
#   Shared utility functions used across all analysis scripts.
#
# Requirements:
#   library(fs)
#   library(purrr)
#   library(ggplot2)

# ── Directory helpers ──────────────────────────────────────────────────────────

#' Check If a Directory Exists
#'
#' @param path         Character or fs_path. Path to the directory to check.
#' @param return_error Logical. If TRUE, stops execution with an error when the
#'   directory does not exist. Default FALSE.
#' @param show_warning Logical. If TRUE, emits a warning when the directory does
#'   not exist (only used when `return_error = FALSE`). Default FALSE.
#'
#' @return Logical. TRUE if the directory exists, FALSE otherwise.
check_dir <- function(path, return_error = FALSE, show_warning = FALSE) {
  exists <- fs::dir_exists(path)

  if (!exists) {
    if (return_error) {
      stop(path, " does not exist")
    } else {
      if (show_warning) warning(path, " does not exist")
      return(FALSE)
    }
  }

  TRUE
}

#' Create a Directory
#'
#' Mimics Python's `os.makedirs(path, exist_ok = ...)` behaviour.
#'
#' @param path         Character or fs_path. Path of the directory to create.
#' @param recursive    Logical. If TRUE, creates all intermediate directories
#'   (like `mkdir -p`). Default TRUE.
#' @param exist_ok     Logical. If FALSE, stops with an error when the directory
#'   already exists. Default FALSE.
#' @param show_warning Logical. If TRUE, emits a warning when the directory
#'   already exists and `exist_ok = TRUE`. Default FALSE.
#'
#' @return Invisible fs_path of the created directory.
create_dir <- function(path, recursive = TRUE, exist_ok = FALSE, show_warning = FALSE) {
  if (check_dir(path)) {
    if (exist_ok) {
      if (show_warning) warning(path, " already exists, files could be overwritten")
    } else {
      stop(path, " already exists")
    }
  } else {
    fs::dir_create(path, recurse = recursive)
  }

  invisible(fs::path(path))
}

# ── Species helpers ────────────────────────────────────────────────────────────

#' Extract Species Name from a TPM File Path
#'
#' @description
#' Extracts the first two-word species identifier (e.g. `"Salmonella_enterica"`)
#' from the basename of a TPM Excel file.
#'
#' @param tpm_file Character. Full path to a TPM `.xlsx` file.
#'
#' @return Character. Species name, e.g. `"Salmonella_enterica"`.
get_species_tpm <- function(tpm_file) {
  tools::file_path_sans_ext(
    stringr::str_extract(basename(tpm_file), "^[A-Za-z]+_[A-Za-z]+")
  )
}

#' Format Species Name as Italic Plot Title Expression
#'
#' @description
#' Converts a raw species key (e.g. `"Salmonella_enterica"`) into a plotmath
#' expression for use as a ggplot title. Genus and epithet are italicised;
#' strain names (third word) are rendered in plain text. Display names are
#' looked up from `PLOT_SPECIES_NAMES`, falling back to underscore-to-space
#' conversion when the key is not found.
#'
#' @param species   Character. Raw species key, e.g. `"Salmonella_enterica"`.
#' @param linebreak Logical. If TRUE, the strain name is placed on a second line
#'   via `atop()` — useful for axis labels where horizontal space is limited.
#'   If FALSE (default), the strain name follows on the same line separated by
#'   a thin space.
#'
#' @return A plotmath `call` object suitable for `labs(title = ...)` or
#'   `scale_*_discrete(labels = ...)`.
format_species_title <- function(species, linebreak = FALSE) {
  display <- PLOT_SPECIES_NAMES[species]
  display <- ifelse(is.na(display), gsub("_", " ", species), display)

  parts <- strsplit(display, " ")[[1]]

  if (length(parts) == 3) {
    genus_epithet <- paste(parts[1:2], collapse = " ")
    strain        <- parts[3]

    if (linebreak) {
      bquote(atop(italic(.(genus_epithet)), plain(.(strain))))
    } else {
      bquote(italic(.(genus_epithet)) ~ plain(.(strain)))
    }
  } else {
    bquote(italic(.(display)))
  }
}

# ── Data helpers ───────────────────────────────────────────────────────────────

#' Extract Protein IDs Common to Both RNA and Protein Datasets
#'
#' @param rna_data     Tibble. Must contain a `Protein_id` column.
#' @param protein_data Tibble. Must contain a `Protein_id` column.
#'
#' @return Character vector of shared `Protein_id` values.
extract_common_ids <- function(rna_data, protein_data) {
  intersect(
    dplyr::pull(rna_data,     Protein_id) |> unique(),
    dplyr::pull(protein_data, Protein_id) |> unique()
  )
}

#' Log2 Transform by Group
#'
#' @description
#' Computes mean log2-transformed values grouped by specified columns.
#' Zero values are excluded before transformation (treated as not-detected).
#' Use `pseudocount = TRUE` for TPM/RNA data to retain zero-expressed genes.
#'
#' @param data        Tibble with expression data.
#' @param value_col   Character. Name of the column containing raw values.
#' @param group_cols  Character vector. Columns to group by.
#'   Default `c("Protein_id", "Treatment")`.
#' @param pseudocount Logical. If TRUE, adds +1 before log2 (useful for TPM).
#'   Default FALSE.
#'
#' @return Tibble with a `mean_log2` column per group combination.
log2_transform <- function(
  data,
  value_col,
  group_cols  = c("Protein_id", "Treatment"),
  pseudocount = FALSE
) {
  data |>
    dplyr::filter(.data[[value_col]] > 0) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(
      mean_log2 = mean(
        log2(.data[[value_col]] + if (pseudocount) 1 else 0),
        na.rm = TRUE
      ),
      .groups = "drop"
    )
}

#' Compute SD of Log2 Treatment Means per Gene
#'
#' @description
#' Computes mean expression per treatment, log2-transforms, then computes the
#' standard deviation across treatments. This matches the label
#' "Standard deviation of log2(mean TPM)" used in figures.
#'
#' @param data        Tibble. Long-format data with `Protein_id`, `Treatment`,
#'   and a value column.
#' @param value_col   Character. Name of the column containing raw values.
#' @param out_col     Character. Name for the output SD column.
#' @param group_cols  Character vector. Outer grouping columns.
#'   Default `c("Protein_id", "Species")`.
#' @param pseudocount Logical. If TRUE, adds +1 before log2. Default FALSE.
#'
#' @return Tibble with columns from `group_cols` plus `<out_col>`.
log2_mean_sd <- function(
  data,
  value_col,
  out_col,
  group_cols  = c("Protein_id", "Species"),
  pseudocount = FALSE
) {
  data |>
    dplyr::filter(.data[[value_col]] > 0) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(union(group_cols, "Treatment")))) |>
    dplyr::summarise(
      mean_val = mean(.data[[value_col]], na.rm = TRUE),
      .groups  = "drop"
    ) |>
    dplyr::mutate(
      log2_mean = log2(mean_val + if (pseudocount) 1 else 0)
    ) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(
      !!out_col := sd(log2_mean, na.rm = TRUE),
      .groups = "drop"
    )
}

# ── Plot helpers ───────────────────────────────────────────────────────────────

#' Save a ggplot to One or More File Formats
#'
#' @description
#' Convenience wrapper around `ggplot2::ggsave()` that saves a plot to all
#' formats supplied in `out_format`, appending the appropriate extension to
#' `filepath`.
#'
#' @param plot       ggplot object. The plot to save.
#' @param filepath   Character. Full path **without** extension.
#' @param out_format Character vector. File extensions, e.g. `c("png", "pdf")`.
#'   Defaults to the global `FIGURE_FORMAT` variable.
#' @param width      Numeric. Width in inches. Default 10.
#' @param height     Numeric. Height in inches. Default 6.
#' @param units      Character. Size units passed to `ggsave()`. Default `"in"`.
#' @param dpi        Numeric. Resolution for raster formats. Default 300.
#'
#' @return NULL invisibly.
save_plot <- function(
  plot,
  filepath,
  out_format = FIGURE_FORMAT,
  width      = 10,
  height     = 6,
  units      = "in",
  dpi        = 300
) {
  purrr::walk(out_format, function(ext) {
    path <- paste0(filepath, ".", ext)
    ggplot2::ggsave(
      filename = path,
      plot     = plot,
      width    = width,
      height   = height,
      dpi      = dpi,
      units    = units
    )
    message("Saved: ", path)
  })
  invisible(NULL)
}

#' Publication-Ready ggplot Theme
#'
#' @description
#' A clean, minimal `theme_classic()`-based theme suited for print-quality
#' figures. Passes `base_size` and `base_family` to `theme_classic()`.
#'
#' @param base_size      Numeric. Base font size in points. Default 20.
#' @param base_family    Character. Base font family. Default `"Arial"`.
#' @param legend_position Character or NULL. Passed to `legend.position`.
#'   NULL uses ggplot2 default (right). Default NULL.
#'
#' @return A ggplot2 `theme` object.
theme_publication <- function(
  base_size       = 20,
  base_family     = "Arial",
  legend_position = NULL
) {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      axis.title     = ggplot2::element_text(size = base_size + 2, colour = "#2c2c2a"),
      axis.text      = ggplot2::element_text(size = base_size,     colour = "#444441"),
      axis.line      = ggplot2::element_line(linewidth = 0.4, colour = "#888780"),
      axis.ticks     = ggplot2::element_line(linewidth = 0.3, colour = "#888780"),

      legend.position = legend_position,
      legend.text     = ggplot2::element_text(size = base_size + 2),

      panel.grid.major.y = ggplot2::element_line(
        colour    = "#d3d1c7",
        linewidth = 0.3,
        linetype  = "dashed"
      ),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      panel.border       = ggplot2::element_rect(colour = "grey70"),
      panel.background   = ggplot2::element_rect(fill = "white", colour = NA),

      plot.title    = ggplot2::element_text(
        size   = base_size + 4,
        face   = "bold",
        colour = "#2c2c2a",
        margin = ggplot2::margin(b = 4)
      ),
      plot.subtitle = ggplot2::element_text(
        size   = base_size + 2,
        colour = "#5f5e5a",
        margin = ggplot2::margin(b = 8)
      ),
      plot.caption  = ggplot2::element_text(
        size   = base_size,
        colour = "#888780",
        hjust  = 0,
        margin = ggplot2::margin(t = 6)
      ),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),

      strip.background = ggplot2::element_blank(),
      strip.text       = ggplot2::element_text(
        size   = base_size,
        face   = "bold",
        colour = "#2c2c2a"
      )
    )
}

# ── Colour palettes ────────────────────────────────────────────────────────────

#' Species colour palette — one tint per species (extend as needed)
species_colours <- c(
  "#1D9E75", # teal-400
  "#7F77DD", # purple-400
  "#D85A30", # coral-400
  "#378ADD", # blue-400
  "#639922", # green-400
  "#BA7517"  # amber-400
)

#' Treatment colour palette — cycles if more than 8 conditions
treatment_colours <- c(
  "#1D9E75", # teal
  "#7F77DD", # purple
  "#D85A30", # coral
  "#378ADD", # blue
  "#BA7517", # amber
  "#639922", # green
  "#D4537E", # pink
  "#888780"  # gray
)
