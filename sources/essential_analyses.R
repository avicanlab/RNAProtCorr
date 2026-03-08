# ============================================================================
# Essential genes analyses
# ============================================================================
#
# Description:
# ============================================================================

SOURCED_AS_MODULE <- TRUE
source("sources/correlation.R")

library(tidyverse)
library(readxl)
library(ggstar)

# ============================================================================
# CONFIGURATION
# ============================================================================

#' Species Abbreviation Mapping
#' Maps full species names to abbreviations for data parsing
species_essential2tpm_map <- c(
    "Senterica"           = "Salmonella_enterica",
    "Saureus"             = "Staphylococcus_aureus",
    "Ypseudotuberculosis" = "Yersinia_pseudotuberculosis"
)

#' Species Abbreviation Mapping
#' Maps full species names to abbreviations for data parsing
species_abbv_map <- c(
    "Salmonella_enterica"         = "SALMT",
    "Staphylococcus_aureus"       = "MSSA",
    "Yersinia_pseudotuberculosis" = "YPSTB"
)

#' Treatment List
#' Control and stress treatment conditions
treatment_list <- c(
    "Ctrl", "As", "Bs", "Mig", "Li", "Nd", "Ns", "Oss", "Oxs", "Sp", "Tm"
)

# ============================================================================
# DATA LOADING
# ============================================================================

#' Read TPM Data from Excel
#'
#' @description
#' Loads TPM expression data, reshapes to long format with Treatment
#' and Replicate columns.
#'
#' @param dataset_path Chr. Path to Excel file containing TPM data.
#' @param sp_abbv Chr. Species abbreviation for column pattern matching.
#'
#' @return Tibble with columns: Species, Protein_id, Treatment, Replicate, TPM
read_tpm <- function(dataset_path, sp_abbv) {
    select_columns <- c("Species", "New_locus_tag", "Old_locus_tag")

    pattern <- paste0(
        "^", sp_abbv, "_(",
        paste(treatment_list, collapse = "|"),
        ")_(\\d+) \\(GE\\) - TPM$"
    )

    read_excel(dataset_path, col_names = TRUE) %>%
        select(
            all_of(select_columns),
            matches(paste0(
                "(", paste(treatment_list, collapse = "|"), ").*TPM|",
                "TPM.*(", paste(treatment_list, collapse = "|"), ")"
            ))
        ) %>%
        mutate(
            New_locus_tag = if_else(
                is.na(New_locus_tag) | New_locus_tag == "N/A",
                Old_locus_tag,
                New_locus_tag
            )
        ) %>%
        select(-Old_locus_tag) %>%
        rename(Protein_id = New_locus_tag) %>%
        pivot_longer(
            cols      = -c(Species, Protein_id),
            names_to  = "col_name",
            values_to = "TPM"
        ) %>%
        mutate(
            Treatment = str_match(col_name, pattern)[, 2],
            Replicate = str_match(col_name, pattern)[, 3]
        ) %>%
        select(Species, Protein_id, Treatment, Replicate, TPM)
}

#' Read Essential Genes Vector from Excel
#'
#' @description
#' Loads an essential genes Excel file and returns a flat character vector
#' of gene IDs.
#'
#' @param filepath Chr. Path to essential genes Excel file.
#'
#' @return Character vector of essential gene IDs.
read_essential_vec <- function(filepath) {
    read_excel(filepath) %>%
        as.list() %>%
        unlist() %>%
        as.vector()
}

# ============================================================================
# DATA TRANSFORMATION
# ============================================================================

#' Compute Global Log2 Mean per Gene
#'
#' @description
#' Computes a single mean log2(value) per Protein_id across ALL replicates
#' and treatments combined. Zeros are excluded (not detected).
#' Use pseudocount = TRUE only for TPM/RNA data.
#'
#' @param data Tibble. Long-format data with Protein_id and value column.
#' @param value_col Chr. Name of column containing raw values.
#' @param out_col Chr. Name of the output mean column.
#' @param pseudocount Logical. If TRUE adds +1 before log2 (for TPM). Default FALSE.
#'
#' @return Tibble with columns: Protein_id, <out_col>
log2_global_mean <- function(data, value_col, out_col, pseudocount = FALSE) {
    data %>%
        filter(.data[[value_col]] > 0) %>%
        group_by(Protein_id) %>%
        summarise(
            !!out_col := mean(
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
log2_mean_sd <- function(data, value_col, out_col, pseudocount = FALSE) {
    data %>%
        filter(.data[[value_col]] > 0) %>%
        group_by(Protein_id, Treatment) %>%
        summarise(
            mean_val = mean(!!sym(value_col), na.rm = TRUE),
            .groups = "drop"
        ) %>%
        mutate(
            log2_mean = log2(mean_val + if (pseudocount) 1 else 0)
        ) %>%
        group_by(Protein_id) %>%
        summarise(
            !!out_col := sd(log2_mean, na.rm = TRUE),
            .groups = "drop"
        )
}

# ============================================================================
# VISUALIZATION
# ============================================================================

#' Save Plot to PDF and PNG
#'
#' @description
#' Saves a ggplot object to both PDF and PNG at the given filepath (no ext).
#'
#' @param plot ggplot object.
#' @param filepath Chr. Full path without file extension.
#' @param width Num. Width in inches (default 10).
#' @param height Num. Height in inches (default 6).
#'
#' @return NULL invisibly.
save_plot <- function(plot, filepath, width = 10, height = 6) {
    walk(c(".pdf", ".png"), function(ext) {
        path <- paste0(filepath, ext)
        ggsave(path, plot = plot, width = width, height = height, dpi = 300)
        message("Saved: ", path)
    })
}

#' Plot Essential vs Non-Essential Gene Expression Distribution
#'
#' @description
#' Overlapping frequency histogram comparing global Log2 mean expression
#' between essential and non-essential genes, with a two-sided Wilcoxon
#' p-value matching the original analysis behaviour.
#'
#' @param all_df Tibble. All genes with a global Log2 mean column
#'   (output of log2_global_mean()).
#' @param essential_vec Chr vector. IDs of essential genes.
#' @param log2_col Chr. Name of the Log2 mean column in all_df.
#' @param x_lab Chr. Short label for the data type (e.g. "TPM", "iBAQ").
#' @param species Chr. Species name for plot title.
#' @param output_path Chr. Directory for saving output files.
#'
#' @return NULL invisibly.
plot_essential_distribution <- function(
  all_df, essential_vec, log2_col, x_lab, species, output_path
) {
    # Tag each gene and drop non-finite values
    condition_df <- all_df %>%
        mutate(
            group = factor(
                ifelse(Protein_id %in% essential_vec, "Essential", "Non-Essential"),
                levels = c("Non-Essential", "Essential")
            )
        ) %>%
        filter(is.finite(.data[[log2_col]]))

    # Two-sided Wilcoxon (matches original wilcox.test(Mean ~ Condition))
    p_value <- wilcox.test(
        as.formula(paste(log2_col, "~ group")),
        data = condition_df
    )$p.value

    p_label <- sprintf("p-val: %.2e", p_value)

    x_label <- bquote(
        .(ifelse(x_lab == "TPM", "mRNA", "Protein")) ~
            "- Mean" ~ .(x_lab) ~ (Log[2])
    )

    p <- ggplot(
        condition_df,
        aes(x = .data[[log2_col]], y = after_stat(density), fill = group)
    ) +
        geom_histogram(position = "identity", alpha = 0.5, bins = 40) +
        scale_fill_manual(
            values = c("Non-Essential" = "#F4A8A0", "Essential" = "#80CEC8")
        ) +
        annotate(
            "text",
            x = Inf, y = Inf, label = p_label,
            hjust = 1.1, vjust = 1.5, size = 3.5, fontface = "italic"
        ) +
        labs(
            title = species,
            x     = x_label,
            y     = "Frequency",
            fill  = NULL
        ) +
        theme_minimal() +
        theme(
            plot.title      = element_text(face = "italic"),
            legend.position = c(0.85, 0.85)
        )

    save_plot(p, file.path(
        output_path, paste0(species, "_", x_lab, "_essential_distribution")
    ))

    invisible(NULL)
}

get_correlation <- function(
  tpm_df, protQ_df, essential_vec, log2_col
) {
    essential_tpm_df <- tpm_df %>%
        filter(Protein_id %in% essential_vec) %>%
        filter(is.finite(mean_Log2_TPM))

    essential_prot_df <- protQ_df %>%
        filter(Protein_id %in% essential_vec) %>%
        filter(is.finite(.data[[log2_col]]))

    corr_all <- calculate_correlation(tpm_df, protQ_df, log2_col)
    corr_essential <- calculate_correlation(
        essential_tpm_df, essential_prot_df, log2_col
    )
    list(corr_all = corr_all, corr_essential = corr_essential)
}

plot_essential_vs_all_correlation <- function(
  corr_all_df, corr_essential_df, species, measurement, output_path
) {
    # Join both correlation sets by treatment
    plot_df <- inner_join(
        corr_all_df %>% distinct(Treatment, R) %>% rename(R_all = R),
        corr_essential_df %>% distinct(Treatment, R) %>% rename(R_essential = R),
        by = "Treatment"
    )

    # Axis range: shared limits with a little padding, for the diagonal
    axis_min <- floor(min(plot_df$R_all, plot_df$R_essential) * 10) / 10 - 0.05
    axis_max <- ceiling(max(plot_df$R_all, plot_df$R_essential) * 10) / 10 + 0.05

    p <- ggplot(plot_df, aes(x = R_all, y = R_essential, color = Treatment)) +
        # Diagonal reference line (essential == all)
        geom_abline(
            slope     = 1,
            intercept = 0,
            linetype  = "dashed",
            color     = "grey60",
            linewidth = 0.6
        ) +
        geom_point(size = 4, alpha = 0.9) +
        scale_color_manual(values = setNames(
            scales::hue_pal()(nrow(plot_df)),
            plot_df$Treatment
        )) +
        coord_fixed(
            ratio  = 1,
            xlim   = c(0.4, 1),
            ylim   = c(0.4, 1)
        ) +
        labs(
            title = bquote(italic(.(gsub("_", " ", species)))),
            x     = "All genes'\nmRNA-protein level correlation",
            y     = "Essential genes'\nmRNA-protein level correlation",
            color = "Treatment"
        ) +
        theme_bw() +
        theme(
            plot.title      = element_text(face = "italic", hjust = 0.5),
            axis.title      = element_text(size = 10),
            legend.position = "right" # matches the figure (no legend shown)
        )

    save_plot(
        plot     = p,
        filepath = file.path(output_path, paste0(species, "_", measurement, "_essential_vs_all_correlation")),
        width    = 5,
        height   = 5
    )

    invisible(NULL)
}

#' Plot Essential vs Non-Essential SD Distribution
#'
#' @description
#' Frequency polygon comparing log2 expression standard deviation between
#' essential and non-essential genes, with a two-sided Wilcoxon p-value.
#'
#' @param all_df Tibble. All genes with SD column (output of log2_global_sd()).
#' @param essential_vec Chr vector. IDs of essential genes.
#' @param sd_col Chr. Name of the SD column in all_df.
#' @param x_lab Chr. Short label for the data type (e.g. "TPM", "iBAQ").
#' @param species Chr. Species name for plot title.
#' @param output_path Chr. Directory for saving output files.
#'
#' @return NULL invisibly.
plot_sd_distribution <- function(
  all_df, essential_vec, sd_col, x_lab, species, output_path, is_tpm = FALSE
) {
    condition_df <- all_df %>%
        mutate(
            group = factor(
                ifelse(Protein_id %in% essential_vec, "Essential", "Non-essential"),
                levels = c("Non-essential", "Essential")
            )
        ) %>%
        filter(is.finite(.data[[sd_col]]))

    p_value <- wilcox.test(
        as.formula(paste(sd_col, "~ group")),
        data = condition_df
    )$p.value

    x_label <- if (is_tpm) {
        bquote("Standard deviation of mRNA - " ~ log[2](mean ~ TPM))
    } else {
        bquote("Standard deviation of protein - " ~ log[2](mean ~ .(x_lab)))
    }

    p <- ggplot(
        condition_df,
        aes(x = .data[[sd_col]], color = group)
    ) +
        geom_step(
            aes(y = after_stat(density)),
            stat = "bin",
            binwidth = 0.05,
            linewidth = 0.8,
            direction = "mid" # centers the step on each bin
        ) +
        scale_color_manual(
            values = c("Non-essential" = "#F4A8A0", "Essential" = "#80CEC8")
        ) +
        annotate(
            "text",
            x = Inf, y = Inf,
            label = sprintf("p-value: %.2e", p_value),
            hjust = 1.1, vjust = 1.5,
            size = 3.5,
            fontface = "italic"
        ) +
        labs(
            title = bquote(italic(.(gsub("_", " ", species)))),
            x     = x_label,
            y     = "Frequency",
            color = NULL
        ) +
        theme_bw() +
        theme(
            plot.title        = element_text(hjust = 0.5),
            legend.position   = c(0.75, 0.85),
            legend.background = element_rect(fill = "white", color = NA)
        )

    save_plot(
        plot     = p,
        filepath = file.path(output_path, paste0(species, "_", x_lab, "_sd_distribution")),
        width    = 5,
        height   = 5
    )

    invisible(NULL)
}

# ============================================================================
# MAIN ORCHESTRATION
# ============================================================================

#' Process Single Species Essential Gene Analysis
#'
#' @description
#' Computes global log2 means for all quantification types and generates
#' essential vs non-essential distribution plots for each.
#'
#' @param essential_vec Chr vector. Essential gene IDs for this species.
#' @param tpm_species Tibble. Long-format raw TPM data (output of read_tpm()).
#' @param intensity_species Tibble. Long-format raw intensity data for this species.
#' @param RI_species Tibble. Long-format raw RI data for this species.
#' @param iBAQ_species Tibble. Long-format raw iBAQ data for this species.
#' @param iBAQ_mc_species Tibble. Long-format raw iBAQ MC data for this species.
#' @param species Chr. Species name.
#' @param output_path Chr. Directory for saving plots.
#'
#' @return NULL invisibly.
process_species <- function(
  essential_vec,
  tpm_species,
  intensity_species,
  RI_species,
  iBAQ_species,
  iBAQ_mc_species,
  species,
  output_path,
  pseudocount = FALSE, # set TRUE to reproduce old results
  common_ids_only = FALSE # set TRUE to restrict to genes in both TPM and proteomics
) {
    # --- Optionally filter to common IDs across all datasets ----------------
    if (common_ids_only) {
        common_ids <- Reduce(intersect, list(
            tpm_species %>% pull(Protein_id) %>% unique(),
            intensity_species %>% pull(Protein_id) %>% unique(),
            iBAQ_species %>% pull(Protein_id) %>% unique(),
            iBAQ_mc_species %>% pull(Protein_id) %>% unique()
        ))
        message("  Common IDs across all datasets: ", length(common_ids))

        tpm_species <- tpm_species %>% filter(Protein_id %in% common_ids)
        intensity_species <- intensity_species %>% filter(Protein_id %in% common_ids)
        RI_species <- RI_species %>% filter(Protein_id %in% common_ids)
        iBAQ_species <- iBAQ_species %>% filter(Protein_id %in% common_ids)
        iBAQ_mc_species <- iBAQ_mc_species %>% filter(Protein_id %in% common_ids)
    }

    # --- Compute global means -----------------------------------------------
    tpm_mean <- log2_global_mean(tpm_species, "TPM", "mean_Log2_TPM", pseudocount = TRUE)
    intensity_mean <- log2_global_mean(intensity_species, "Intensity", "mean_Log2_intensity", pseudocount = pseudocount)
    RI_mean <- log2_global_mean(RI_species, "RI", "mean_Log2_RI", pseudocount = pseudocount)
    iBAQ_mean <- log2_global_mean(iBAQ_species, "iBAQ", "mean_Log2_iBAQ", pseudocount = pseudocount)
    iBAQ_mc_mean <- log2_global_mean(iBAQ_mc_species, "iBAQ", "mean_Log2_iBAQ", pseudocount = pseudocount)

    # --- Compute global SDs -------------------------------------------------
    tpm_sd <- log2_mean_sd(tpm_species, "TPM", "sd_Log2_TPM", pseudocount = TRUE)
    intensity_sd <- log2_mean_sd(intensity_species, "Intensity", "sd_Log2_intensity", pseudocount = pseudocount)
    RI_sd <- log2_mean_sd(RI_species, "RI", "sd_Log2_RI", pseudocount = pseudocount)
    iBAQ_sd <- log2_mean_sd(iBAQ_species, "iBAQ", "sd_Log2_iBAQ", pseudocount = pseudocount)
    iBAQ_mc_sd <- log2_mean_sd(iBAQ_mc_species, "iBAQ", "sd_Log2_iBAQ", pseudocount = pseudocount)

    # --- Distribution plots (mean) ------------------------------------------
    list(
        list(df = tpm_mean, col = "mean_Log2_TPM", lab = "TPM", is_tpm = TRUE),
        list(df = intensity_mean, col = "mean_Log2_intensity", lab = "Intensity", is_tpm = FALSE),
        list(df = RI_mean, col = "mean_Log2_RI", lab = "RI", is_tpm = FALSE),
        list(df = iBAQ_mean, col = "mean_Log2_iBAQ", lab = "iBAQ", is_tpm = FALSE),
        list(df = iBAQ_mc_mean, col = "mean_Log2_iBAQ", lab = "iBAQ_mc", is_tpm = FALSE)
    ) %>%
        walk(~ plot_essential_distribution(
            all_df        = .x$df,
            essential_vec = essential_vec,
            log2_col      = .x$col,
            x_lab         = .x$lab,
            species       = species,
            output_path   = output_path
        ))

    # --- SD distribution plots ----------------------------------------------
    list(
        list(df = tpm_sd, col = "sd_Log2_TPM", lab = "TPM", is_tpm = TRUE),
        list(df = intensity_sd, col = "sd_Log2_intensity", lab = "Intensity", is_tpm = FALSE),
        list(df = RI_sd, col = "sd_Log2_RI", lab = "RI", is_tpm = FALSE),
        list(df = iBAQ_sd, col = "sd_Log2_iBAQ", lab = "iBAQ", is_tpm = FALSE),
        list(df = iBAQ_mc_sd, col = "sd_Log2_iBAQ", lab = "iBAQ_mc", is_tpm = FALSE)
    ) %>%
        walk(~ plot_sd_distribution(
            all_df        = .x$df,
            essential_vec = essential_vec,
            sd_col        = .x$col,
            x_lab         = .x$lab,
            is_tpm        = .x$is_tpm,
            species       = species,
            output_path   = output_path
        ))

    # --- Essential vs all correlation scatter -------------------------------
    tpm_log2 <- tpm_species %>%
        log2_transform(value_col = "TPM") %>%
        rename(mean_Log2_TPM = mean_log2)

    list(
        list(df = intensity_species, col = "Intensity_meanlog2", lab = "Intensity"),
        list(df = RI_species, col = "RI_meanlog2", lab = "RI"),
        list(df = iBAQ_species, col = "iBAQ_meanlog2", lab = "iBAQ"),
        list(df = iBAQ_mc_species, col = "iBAQ_meanlog2", lab = "iBAQ_mc")
    ) %>%
        walk(function(item) {
            pair_corr <- get_correlation(
                tpm_log2, item$df, essential_vec, item$col
            )
            plot_essential_vs_all_correlation(
                pair_corr$corr_all,
                pair_corr$corr_essential,
                species,
                item$lab,
                output_path
            )
        })

    invisible(NULL)
}

#' Run Full Essential Genes Analysis
#'
#' @description
#' Loads all data, pairs TPM and essential files per species, and runs
#' the complete distribution analysis for each species.
#'
#' @param essential_path Chr. Directory containing essential genes Excel files.
#' @param tpm_path Chr. Directory containing TPM Excel files.
#' @param RI_path Chr. Path to RI TSV file.
#' @param iBAQ_path Chr. Path to iBAQ TSV file.
#' @param iBAQ_mc_path Chr. Path to iBAQ MC TSV file.
#' @param output_path Chr. Root output directory.
#'
#' @return NULL invisibly.
main_analysis <- function(
  essential_path,
  tpm_path,
  intensity_path,
  RI_path,
  iBAQ_path,
  iBAQ_mc_path,
  output_path,
  pseudocount = FALSE, # set TRUE to reproduce old results
  common_ids_only = FALSE # set TRUE to restrict to common genes
) {
    # Load protein quantification tables (all species, long format)
    intensity_data <- read_tsv(
        intensity_path,
        col_names = TRUE,
        show_col_types = FALSE,
        skip_empty_rows = TRUE
    )
    RI_data <- read_tsv(
        RI_path,
        col_names = TRUE,
        show_col_types = FALSE,
        skip_empty_rows = TRUE
    )
    iBAQ_data <- read_tsv(
        iBAQ_path,
        col_names = TRUE,
        show_col_types = FALSE,
        skip_empty_rows = TRUE
    )
    iBAQ_mc_data <- read_tsv(
        iBAQ_mc_path,
        col_names = TRUE,
        show_col_types = FALSE,
        skip_empty_rows = TRUE
    )

    # File path helpers
    extract_species_tpm <- function(filepath) {
        tools::file_path_sans_ext(str_extract(basename(filepath), "^[A-Za-z]+_[A-Za-z]+"))
    }
    extract_species_essential <- function(filepath) {
        key <- tools::file_path_sans_ext(
            str_extract(basename(filepath), "Essential_genes_([A-Za-z]+)", group = 1)
        )
        species_essential2tpm_map[key]
    }

    # Build paired file list keyed by species name
    tpm_files <- list.files(tpm_path, pattern = "\\.xlsx$", full.names = TRUE)
    essential_files <- list.files(essential_path, pattern = "\\.xlsx$", full.names = TRUE)

    tpm_species_names <- map_chr(tpm_files, extract_species_tpm)
    essential_species_names <- map_chr(essential_files, extract_species_essential)

    species_files <- set_names(tpm_species_names) %>%
        map(function(species) {
            list(
                tpm       = tpm_files[[match(species, tpm_species_names)]],
                essential = essential_files[[match(species, essential_species_names)]]
            )
        })

    # Process each species
    iwalk(species_files, function(files, species) {
        message("Processing species: ", species)

        abbv <- species_abbv_map[species]
        if (is.na(abbv)) {
            warning("No abbreviation found for species '", species, "', skipping.")
            return(invisible(NULL))
        }

        # Load raw long-format TPM (no averaging yet)
        tpm_species <- read_tpm(files$tpm, abbv) %>%
            mutate(
                Species    = species,
                Treatment  = ifelse(Treatment == "Mig", "Hyp", Treatment),
                Protein_id = sub("^pi", "PI", Protein_id, ignore.case = FALSE)
            ) %>%
            filter(!is.na(Protein_id), !is.na(TPM))

        essential_vec <- read_essential_vec(files$essential)

        # Filter quantification data to this species (keep raw for global mean)
        intensity_species <- intensity_data %>% filter(Species == species)
        RI_species <- RI_data %>% filter(Species == species)
        iBAQ_species <- iBAQ_data %>% filter(Species == species)
        iBAQ_mc_species <- iBAQ_mc_data %>% filter(Species == species)

        process_species(
            essential_vec     = essential_vec,
            tpm_species       = tpm_species,
            intensity_species = intensity_species,
            RI_species        = RI_species,
            iBAQ_species      = iBAQ_species,
            iBAQ_mc_species   = iBAQ_mc_species,
            species           = species,
            output_path       = output_path,
            pseudocount       = pseudocount, # passed through
            common_ids_only   = common_ids_only # passed through
        )
    })

    message("All species processed.")
}

# ============================================================================
# ENTRY POINT
# ============================================================================

#' Run analysis with default development paths
#' @return NULL invisibly.
test_analysis <- function() {
    main_analysis(
        essential_path  = "DATA/Essential genes",
        tpm_path        = "DATA/Read_counts/",
        intensity_path  = "Analyses/Intensity_data.tsv",
        RI_path         = "Analyses/RI_data.tsv",
        iBAQ_path       = "Analyses/iBAQ_data.tsv",
        iBAQ_mc_path    = "Analyses/iBAQ_mc_data.tsv",
        output_path     = "Analyses",
        pseudocount     = FALSE, # match old +1 behaviour
        common_ids_only = TRUE # match old SA_common_IDs filter
    )
}

if (interactive()) {
    test_analysis()
}
