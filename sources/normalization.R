# ============================================================================
# Protein Normalization Method Comparison
# ============================================================================
library(NormalyzerDE)

NORM_METHODS <- c("Log2", "VSN", "CycLoess", "Quantile", "Median", "Mean", "GI", "RLR")

# ============================================================================
# HELPERS
# ============================================================================

create_design_for_matrix <- function(data_path, tmp_dir = tempdir()) {
    sample_cols <- read.table(data_path,
        sep = "\t", header = TRUE,
        check.names = FALSE, nrows = 1
    ) %>%
        select(-Protein_id) %>%
        colnames()

    design_df <- tibble(
        sample = sample_cols,
        group  = str_remove(sample_cols, "^R\\d+_")
    )

    temp_file <- tempfile(pattern = "design_", tmpdir = tmp_dir, fileext = ".tsv")
    write.table(design_df, temp_file, sep = "\t", row.names = FALSE, quote = FALSE)
    temp_file
}

transform_data <- function(data, protQ, tmp_dir = tempdir()) {
    message("Transforming ", protQ, " data...")

    data_wide <- data %>%
        mutate(sample = paste0(Replicate, "_", Treatment)) %>%
        pivot_wider(
            id_cols     = Protein_id,
            names_from  = sample,
            values_from = all_of(protQ)
        ) %>%
        select(Protein_id, sort(colnames(.)[-1]))

    temp_file <- tempfile(
        pattern = paste0("transformed_", protQ, "_"),
        tmpdir  = tmp_dir,
        fileext = ".tsv"
    )
    write.table(data_wide, temp_file, sep = "\t", row.names = FALSE, quote = FALSE)
    temp_file
}

normalize_data <- function(data_path, design_path, job_name, output_path) {
    message("Normalizing ", job_name, " data...")
    dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

    normalyzer(
        jobName = job_name,
        dataPath = data_path,
        designPath = design_path,
        outputDir = output_path,
        writeReportAsPngs = TRUE,
        zeroToNA = TRUE
    )

    invisible(NULL)
}

plot_pca <- function(norm_matrix, design_table, method_name) {
    mat <- norm_matrix %>%
        select(-any_of(c("CV", "anovaP"))) %>%
        column_to_rownames("Protein_id") %>%
        select(where(is.numeric)) %>%
        as.matrix()
    mat <- mat[complete.cases(mat), ]

    pca <- prcomp(t(mat), scale. = TRUE, center = TRUE)
    var_exp <- round(summary(pca)$importance[2, 1:2] * 100, 1)

    design_df <- read.table(design_table, header = TRUE, sep = "\t") %>%
        as_tibble() %>%
        mutate(
            Treatment = group,
            Replicate = str_remove(sample, "_[a-zA-Z]+$")
        )

    pca_df <- as.data.frame(pca$x[, 1:2]) %>%
        rownames_to_column("sample") %>%
        left_join(design_df, by = "sample")

    pca_plot <- ggplot(pca_df, aes(
        x = PC1, y = PC2,
        color = Treatment,
        shape = factor(Replicate),
        label = sample
    )) +
        geom_point(size = 3) +
        labs(
            title = paste("PCA \u2014", method_name),
            x     = paste0("PC1 (", var_exp[1], "%)"),
            y     = paste0("PC2 (", var_exp[2], "%)"),
            shape = "Replicate"
        ) +
        theme_bw()

    list(df = pca_df, plot = pca_plot)
}

process_plot_pca <- function(norm_path, design_path, output_path, method_name) {
    norm_matrix <- read.table(norm_path, sep = "\t", header = TRUE, check.names = FALSE)
    pca_res <- plot_pca(norm_matrix, design_path, method_name)

    tsv_out <- file.path(output_path, paste0(method_name, "_PCA_coordinates.tsv"))
    write.table(pca_res$df, file = tsv_out, sep = "\t", row.names = FALSE, quote = FALSE)
    message("PCA coordinates saved to: ", tsv_out)

    png_out <- file.path(output_path, paste0(method_name, ".png"))
    ggsave(filename = png_out, plot = pca_res$plot, width = 6, height = 5)
    message("PCA plot saved to: ", png_out)

    invisible(NULL)
}

#' Plot PCA Comparison Grid for All Normalization Methods
#'
#' @description
#' Generates a grid of PCA plots — one panel per normalization method,
#' each showing all species side by side as facets. Layout is 2 columns
#' with ceiling(n_methods / 2) rows, matching the reference figure style.
#'
#' @param protQ_data  Tibble. Long-format data with Species, Replicate,
#'   Treatment, Protein_id, and quantification column.
#' @param norm_dir    Chr. Root NormalyzerDE output directory (from
#'   run_normalization_comparison).
#' @param output_path Chr. Directory for saving the combined figure.
#' @return NULL invisibly.
plot_pca_comparison_grid <- function(protQ_data, norm_dir, output_path) {
    species_list <- unique(protQ_data$Species)

    design_map <- map(species_list, function(species) {
        protQ_data %>%
            filter(Species == species) %>%
            distinct(Replicate, Treatment) %>%
            mutate(
                sample    = paste0(Replicate, "_", Treatment),
                group     = Treatment
            )
    }) %>% setNames(species_list)

    pca_all <- map_dfr(NORM_METHODS, function(norm_method) {
        map_dfr(species_list, function(species) {
            norm_file <- file.path(
                output_path,
                species,
                "Normalization",
                paste0(norm_method, "-normalized.txt")
            )

            if (!file.exists(norm_file)) {
                warning("Missing: ", norm_file)
                return(NULL)
            }

            norm_matrix <- read.table(norm_file,
                sep = "\t",
                header = TRUE, check.names = FALSE
            )

            mat <- norm_matrix %>%
                select(-any_of(c("CV", "anovaP"))) %>%
                column_to_rownames("Protein_id") %>%
                select(where(is.numeric)) %>%
                as.matrix()
            mat <- mat[complete.cases(mat), ]

            pca <- prcomp(t(mat), scale. = TRUE, center = TRUE)
            var_exp <- round(summary(pca)$importance[2, 1:2] * 100, 1)

            as.data.frame(pca$x[, 1:2]) %>%
                rownames_to_column("sample") %>%
                left_join(design_map[[species]], by = "sample") %>%
                mutate(
                    Species     = species,
                    method      = norm_method,
                    var_exp_PC1 = var_exp[1],
                    var_exp_PC2 = var_exp[2]
                )
        })
    })

    # Build one plot per method
    method_plots <- map(NORM_METHODS, function(norm_method) {
        df <- pca_all %>% filter(method == norm_method)

        # Per-species variance explained labels (take first row per species)
        axis_labels <- df %>%
            distinct(Species, var_exp_PC1, var_exp_PC2) %>%
            mutate(
                x_lab = paste0("PC1 (", var_exp_PC1, "%)"),
                y_lab = paste0("PC2 (", var_exp_PC2, "%)")
            )

        # Use first species x/y labels for shared axes
        x_lab <- axis_labels$x_lab[1]
        y_lab <- axis_labels$y_lab[1]

        # Species facet labels using PLOT_SPECIES_NAMES
        species_labels <- setNames(
            sapply(unique(df$Species), function(s) {
                display <- PLOT_SPECIES_NAMES[s]
                ifelse(is.na(display), gsub("_", " ", s), display)
            }),
            unique(df$Species)
        )

        ggplot(df, aes(x = PC1, y = PC2, color = Treatment, shape = factor(Replicate))) +
            geom_point(size = 2.5, alpha = 0.9) +
            facet_wrap(~Species,
                nrow = 1, scales = "free",
                labeller = as_labeller(species_labels)
            ) +
            labs(
                title  = norm_method,
                x      = x_lab,
                y      = y_lab,
                color  = "Treatment",
                shape  = "Replicate"
            ) +
            theme_bw(base_size = 10) +
            theme(
                plot.title       = element_text(face = "bold", hjust = 0),
                strip.background = element_rect(fill = "#8FBC8F"),
                strip.text       = element_text(face = "italic", size = 9),
                panel.grid.minor = element_blank(),
                legend.position  = "right"
            )
    })

    # Arrange in 2-column grid with method labels A, B, C...
    combined <- wrap_plots(method_plots, ncol = 2) +
        plot_layout(guides = "collect") +
        plot_annotation(tag_levels = "A") &
        theme(legend.position = "right")

    out_file <- file.path(output_path, paste0("normalization_PCA_comparison"))
    save_plot(combined, out_file,
        width = 14,
        height = 5 * ceiling(length(NORM_METHODS) / 2),
        out_format = FIGURE_FORMAT
    )
    message("PCA comparison grid saved: ", out_file, MSG_FIG_FORMAT)

    invisible(NULL)
}
# ============================================================================
# MAIN
# ============================================================================

#' Run normalization method comparison across all species and quantification types
#'
#' @param protQ_data   Tibble. Long-format data with Species, Replicate,
#'   Treatment, Protein_id, and quantification column.
#' @param protQ        Chr. Name of the quantification column.
#' @param output_path  Chr. Root output directory.
#' @return NULL invisibly.
run_normalization_comparison <- function(protQ_data, protQ, output_path) {
    tmp_dir <- tempdir(check = FALSE)
    # output_dir <- file.path(output_path, "Normalization")
    # dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    species_list <- unique(protQ_data$Species)

    for (species in species_list) {
        message("Normalizing ", species, " data ...")

        sp_data <- protQ_data %>% filter(Species == species)

        data_path <- transform_data(sp_data, protQ, tmp_dir)
        design_path <- create_design_for_matrix(data_path, tmp_dir)

        sp_dir <- file.path(output_path, species)
        dir.create(sp_dir, showWarnings = FALSE)

        normalyzer(
            jobName = "Normalization",
            dataPath = data_path,
            designPath = design_path,
            outputDir = sp_dir,
            writeReportAsPngs = TRUE,
            zeroToNA = TRUE
        )

        pca_dir <- file.path(sp_dir, "Normalization", "PCA")
        dir.create(pca_dir, showWarnings = FALSE)

        for (norm_method in NORM_METHODS) {
            message("  Generating PCA for: ", norm_method)

            norm_file <- file.path(
                sp_dir,
                "Normalization",
                paste0(norm_method, "-normalized.txt")
            )
            if (!file.exists(norm_file)) {
                warning("Normalized file not found for method ", norm_method, ": ", norm_file)
                next
            }

            process_plot_pca(
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

# ============================================================================
# Protein Normalization - Single Method
# ============================================================================

#' Normalize protein data using a single NormalyzerDE method
#'
#' @description
#' Normalizes long-format proteomics data using one normalization method
#' from NormalyzerDE. Returns a tidy long-format tibble with a new
#' normalized intensity column.
#'
#' @param protQ_data  Tibble. Long-format data with columns: Species,
#'   Replicate, Treatment, Protein_id, and the quantification column.
#' @param protQ       Chr. Name of the quantification column to normalize.
#' @param method      Chr. Normalization method. One of: "Log2", "VSN",
#'   "CycLoess", "Quantile", "Median", "Mean", "GI", "RLR".
#' @param output_path Chr. Optional directory to save normalized TSV.
#'
#' @return Tibble in long format with additional column
#'   `<method>_<protQ>` containing normalized values.
run_normalization <- function(protQ_data, protQ, method = "VSN", output_path = NULL) {
    message("Running ", method, " normalization for ", protQ, "...\n")

    tmp_dir <- tempdir(check = FALSE)

    # Normalize per species then recombine
    normalized_data <- map_dfr(unique(protQ_data$Species), function(species) {
        message("  Processing species: ", species)

        sp_data <- protQ_data %>% filter(Species == species)
        data_path <- transform_data(sp_data, protQ, tmp_dir)
        design_path <- create_design_for_matrix(data_path, tmp_dir)

        job_name <- paste0("norm_", species)
        normalyzer(
            jobName = job_name,
            dataPath = data_path,
            designPath = design_path,
            outputDir = tmp_dir,
            skipAnalysis = TRUE,
            zeroToNA = TRUE,
            quiet = TRUE
        )
        norm_file <- file.path(tmp_dir, job_name, paste0(method, "-normalized.txt"))
        if (!file.exists(norm_file)) {
            stop(
                "Normalized file not found: ", norm_file,
                "\nCheck that method '", method, "' is valid."
            )
        }

        out_col <- paste0(method, "_", protQ)

        read.table(norm_file, sep = "\t", header = TRUE, check.names = FALSE) %>%
            select(-any_of(c("CV", "anovaP"))) %>%
            pivot_longer(
                cols      = -Protein_id,
                names_to  = "sample",
                values_to = out_col
            ) %>%
            mutate(
                Species   = species,
                Replicate = str_extract(sample, "^R\\d+"),
                Treatment = str_remove(sample, "^R\\d+_")
            ) %>%
            select(-sample)
    })

    if (!is.null(output_path)) {
        out_file <- file.path(output_path, paste0(method, "_", protQ, "_normalized.tsv"))
        write.table(normalized_data, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
        message(method, " normalized data saved to: ", out_file)
    }

    normalized_data
}
