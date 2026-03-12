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
save_plot <- function(plot, filepath, width = 10, height = 6) {
    walk(c(".pdf", ".png"), function(ext) {
        path <- paste0(filepath, ext)
        ggsave(path, plot = plot, width = width, height = height, dpi = 300)
        message("Saved: ", path)
    })
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
