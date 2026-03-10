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
