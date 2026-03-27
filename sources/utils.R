# Utilities function for all analyses -----------------------------------------
#
# Description:
#   Utilities functions to process all analyses
#
# Inputs:
#
# Outputs (per species):
#
# Requirements:
#   library(fs)

#' Check if a directory exists
#'
#' @param path          Character or fs_path. Path to the directory to check.
#' @param return_error  Logical. If TRUE, stops execution with an error when the
#'                      directory does not exist. Default FALSE.
#' @param show_warning  Logical. If TRUE, emits a warning when the directory does
#'                      not exist (only used when return_error = FALSE). Default FALSE.
#' @return Logical. TRUE if the directory exists, FALSE otherwise.
check_dir <- function(path, return_error = FALSE, show_warning = FALSE) {
    exists <- dir_exists(path) # fs::dir_exists() instead of base dir.exists()

    if (!exists) {
        if (return_error) {
            stop(path, " does not exist")
        } else {
            if (show_warning) warning(path, " does not exist")
            return(FALSE)
        }
    }

    return(TRUE)
}


#' Create a directory
#'
#' Mimics Python's os.makedirs(path, exist_ok=...) behaviour.
#'
#' @param path          Character or fs_path. Path of the directory to create.
#' @param recursive     Logical. If TRUE, creates all intermediate directories
#'                      (like mkdir -p). Default TRUE.
#' @param exist_ok      Logical. If FALSE, stops with an error when the directory
#'                      already exists. Default FALSE.
#' @param show_warning  Logical. If TRUE, emits a warning when the directory
#'                      already exists and exist_ok = TRUE. Default FALSE.
#' @return Invisible fs_path of the created directory.
create_dir <- function(path, recursive = TRUE, exist_ok = FALSE, show_warning = FALSE) {
    if (check_dir(path)) {
        if (exist_ok) {
            if (show_warning) warning(path, " already exists, files could be overwritten")
        } else {
            stop(path, " already exists")
        }
    } else {
        # fs::dir_create() is recursive by default and never warns on existing dirs
        dir_create(path, recurse = recursive)
    }

    return(invisible(path(path))) # return a proper fs_path object
}

get_species_tpm <- function(tpm_file) {
    tools::file_path_sans_ext(
        str_extract(basename(tpm_file), "^[A-Za-z]+_[A-Za-z]+")
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
save_plot <- function(
  plot,
  filepath,
  out_format = FIGURE_FORMAT,
  width = 10,
  height = 6,
  units = "in",
  dpi = 300
) {
    walk(out_format, function(ext) {
        path <- paste0(filepath, ".", ext)
        print(path)
        ggsave(
            path,
            plot = plot,
            width = width,
            height = height,
            dpi = 300,
            units = units
        )
        message("Saved: ", path)
    })
    invisible(NULL)
}


#' Log2 Transform by Group
#'
#' @description
#' Computes mean log2-transformed values grouped by specified columns.
#' Zeros are excluded before transformation (not detected).
#' Use pseudocount = TRUE for TPM/RNA data to keep zero-expressed genes.
#'
#' @param data        Tibble with expression data
#' @param value_col   Character. Name of column containing raw values
#' @param group_cols  Character vector. Columns to group by.
#'                    Default c("Protein_id", "Treatment").
#' @param pseudocount Logical. If TRUE adds +1 before log2 (for TPM). Default FALSE.
#'
#' @return Tibble with mean_log2 values per group combination
log2_transform <- function(
  data,
  value_col,
  group_cols = c("Protein_id", "Treatment"),
  pseudocount = FALSE
) {
    data %>%
        filter(.data[[value_col]] > 0) %>%
        group_by(across(all_of(group_cols))) %>%
        summarise(
            mean_log2 = mean(
                log2(!!sym(value_col) + if (pseudocount) 1 else 0),
                na.rm = TRUE
            ),
            .groups = "drop"
        )
}

#' Compute SD of log2 Treatment Means per Gene
#'
#' @description
#' For TPM: computes mean per treatment first, then log2, then SD across
#' treatments. Matches "Standard deviation log2(mean TPM)" label.
#'
#' @param data Tibble. Long-format data with Protein_id, Treatment, value column.
#' @param value_col Chr. Name of column containing raw values.
#' @param out_col Chr. Name of the output SD column.
#' @param pseudocount Logical. If TRUE adds +1 before log2 (for TPM).
#'
#' @return Tibble with columns: Protein_id, <out_col>
log2_mean_sd <- function(
  data, value_col,
  out_col,
  group_cols = c("Protein_id", "Species"),
  pseudocount = FALSE
) {
    data %>%
        filter(.data[[value_col]] > 0) %>%
        group_by(across(all_of(union(group_cols, "Treatment")))) %>%
        summarise(
            mean_val = mean(!!sym(value_col), na.rm = TRUE),
            .groups = "drop"
        ) %>%
        mutate(
            log2_mean = log2(mean_val + if (pseudocount) 1 else 0)
        ) %>%
        group_by(across(all_of(group_cols))) %>%
        summarise(
            !!out_col := sd(log2_mean, na.rm = TRUE),
            .groups = "drop"
        )
}

extract_common_ids <- function(rna_data, protein_data) {
    prot_ids <- protein_data %>%
        pull(Protein_id) %>%
        unique()
    intersect(
        rna_data %>% pull(Protein_id) %>% unique(),
        prot_ids
    )
}

#' Format Species Name as Italic Plot Title Expression
#'
#' @description
#' Converts a raw species key (e.g. "Salmonella_enterica") into a plotmath
#' expression for use as a ggplot title. Genus and epithet are italicised;
#' strain names (third word) are rendered in plain text. Display names are
#' looked up from PLOT_SPECIES_NAMES, falling back to underscore-to-space
#' conversion if the species is not found in the map.
#'
#' @param species   Chr. Raw species key, e.g. "Salmonella_enterica".
#' @param linebreak Logical. If TRUE, strain name is placed on a second line
#'   via atop() — useful for axis labels where horizontal space is limited.
#'   If FALSE (default), strain name follows on the same line separated by
#'   a thin space.
#'
#' @return A plotmath call object suitable for use in labs(title = ...) or
#'   scale_*_discrete(labels = ...).
format_species_title <- function(species, linebreak = FALSE) {
    display <- PLOT_SPECIES_NAMES[species]
    display <- ifelse(is.na(display), gsub("_", " ", species), display)

    parts <- strsplit(display, " ")[[1]]

    if (length(parts) == 3) {
        genus_epithet <- paste(parts[1:2], collapse = " ")
        strain <- parts[3]

        if (linebreak) {
            bquote(atop(italic(.(genus_epithet)), plain(.(strain))))
        } else {
            bquote(italic(.(genus_epithet)) ~ plain(.(strain)))
        }
    } else {
        bquote(italic(.(display)))
    }
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