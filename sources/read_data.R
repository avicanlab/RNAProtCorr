# Function to read transcriptomic and proteomic data --------------------------
# Developed by Sena Gizem Suer, Yi Yang Lim, Jérôme Arnoux
#
# Description:
#
# Inputs:
#
# Outputs:
#
# Requirements:
#

## Configuration --------------------------------------------------------------

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

## Read transcriptomic data ---------------------------------------------------

#' Read TPM Data from Excel
#'
#' @description
#' Loads Transcript Per Million (TPM) expression data from Excel file,
#' extracts relevant columns, and reshapes to long format with treatment
#' and replicate information.
#'
#' @param dataset_path Character. Path to Excel file containing TPM data
#' @param sp_abbv Character. Species abbreviation for column pattern matching
#'
#' @return Tibble with columns: Species, Protein_id, Treatment, Replicate, TPM
read_tpm <- function(dataset_path, sp_abbv) {
    select_columns <- c("Species", "New_locus_tag", "Old_locus_tag")

    # Build regex pattern: ABBV_Treatment_Replicate (GE) - TPM
    pattern <- paste0(
        "^", sp_abbv, "_(",
        paste(treatment_list, collapse = "|"),
        ")_(\\d+) \\(GE\\) - TPM$"
    )

    data <- read_excel(dataset_path, col_names = TRUE, na = "N/A") %>%
        dplyr::select(
            all_of(select_columns),
            matches(
                paste0(
                    "(",
                    paste(treatment_list, collapse = "|"),
                    ").*TPM|TPM.*(", paste(treatment_list, collapse = "|"), ")"
                )
            )
        ) %>%
        mutate(
            New_locus_tag = if_else(
                is.na(New_locus_tag),
                Old_locus_tag,
                New_locus_tag
            )
        )

    data %>%
        dplyr::select(-Old_locus_tag) %>%
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
        dplyr::select(Species, Protein_id, Treatment, Replicate, TPM)
}


## Read proteomic data --------------------------------------------------------
#' Read a TMT result TSV file
#'
#' @param path Path to the TSV file
#' @return A tibble with the relevant columns
read_tmt <- function(path) {
    readr::read_tsv(
        path,
        col_names = TRUE,
        show_col_types = FALSE,
        col_select = c(all_of(keep_cols), starts_with("channel_")),
        skip_empty_rows = TRUE
    )
}

#' Filter TMT data for a given species and replicate
#'
#' Applies sequential filters:
#'   1. Keep rows matching the species
#'   2. Remove contaminants
#'   3. Keep proteins with > 1 unique peptide
#'   4. Remove rows with missing Gene
#'   5. Deduplicate by keeping the row with highest unique peptides and intensity per gene
#'
#' @param tmt_df     A tibble from read_tmt()
#' @param curr_species   Species name (as it appears in the Organism column)
#' @param curr_replicate Replicate identifier (e.g. "R1")
#' @return A named list with:
#'   - `filtered_tmt`: filtered tibble
#'   - `filter_log`:   named vector of row counts at each step
filter_tmt <- function(tmt_df, curr_species, curr_replicate) {
    step0 <- tmt_df %>%
        mutate(
            total_intensity = as.numeric(.data[[col_total_intensity]]),
            total_intensity = na_if(total_intensity, 0),
            Protein_id      = .data[[col_protein_id]],
            Organism        = str_replace(.data[[col_organism]], " ", "_")
        )

    step1 <- step0 %>%
        filter(grepl(curr_species, .data[[col_organism]]))

    step2 <- step1 %>%
        filter(!grepl("contam_sp", .data[[col_protein]]))

    step3 <- step2 %>%
        filter(.data[[col_unique_peptides]] > 1)

    step4 <- step3 %>%
        filter(total_intensity > 0)

    step5 <- step4 %>%
        filter(!is.na(.data[[col_gene]]))

    step6 <- step5 %>%
        group_by(.data[[col_gene]]) %>%
        filter(
            .data[[col_unique_peptides]] == max(.data[[col_unique_peptides]]),
            .data[[col_total_intensity]] == max(.data[[col_total_intensity]])
        ) %>%
        ungroup()

    filter_log <- c(
        species = curr_species,
        replicate = curr_replicate,
        initial = nrow(step0),
        remove_non_sp = nrow(step1),
        remove_contam = nrow(step2),
        remove_unique_pep = nrow(step3),
        remove_intensity_0 = nrow(step4),
        remove_na_gene = nrow(step5),
        remove_dedup = nrow(step6)
    )

    list(filtered_tmt = step6, filter_log = filter_log)
}

#' Process all TMT TSV files in a directory
#'
#' Reads, filters, and combines all TSV files found recursively under tmt_path.
#' Expects filenames of the form: <Species>_R<n>.tsv
#' Writes a filter log TSV to output_path.
#'
#' @param tmt_path    Path to the directory containing TSV files
#' @param output_path Path to the output directory (must end with "/")
#' @return A unified tibble with Species and Replicate columns prepended
process_tmt <- function(tmt_path, output_path) {
    tsv_files <- list.files(
        tmt_path,
        pattern    = "\\.tsv$",
        recursive  = TRUE,
        full.names = TRUE
    )

    filter_log <- list()
    filtered_tmt <- list()

    for (tsv in tsv_files) {
        curr_species <- sub(".*/(\\w+)_R\\d+\\.tsv$", "\\1", tsv)
        curr_replicate <- sub(".*/\\w+_(R\\d+)\\.tsv$", "\\1", tsv)

        tmt_data <- read_tmt(tsv)
        filter_result <- filter_tmt(tmt_data, curr_species, curr_replicate)

        filtered_tmt[[curr_species]][[curr_replicate]] <- filter_result$filtered_tmt
        filter_log[[curr_species]][[curr_replicate]] <- filter_result$filter_log
    }

    # Combine all species/replicates into a single tibble
    unified_tmt <- map_dfr(names(filtered_tmt), function(species) {
        map_dfr(names(filtered_tmt[[species]]), function(rep) {
            filtered_tmt[[species]][[rep]] %>%
                mutate(Species = species, Replicate = rep)
        })
    }) %>%
        select(Species, Replicate, everything())

    # Build and write filter log
    filter_log_df <- map_dfr(filter_log, ~ bind_rows(.x)) %>%
        mutate(across(-c("species", "replicate"), as.integer))

    message("TMT filtering summary:\n")
    print(filter_log_df)

    write.table(
        filter_log_df,
        paste0(output_path, "filter_log.tsv"),
        sep       = "\t",
        row.names = FALSE,
        quote     = FALSE
    )

    unified_tmt
}

#' Remove proteins with incomplete data across replicates and treatments
#'
#' Keeps only proteins that have at least one value for every unique
#' Species-Replicate-Treatment combination in the dataset.
#'
#' @param tmt_data A tibble with Species, Replicate, and Treatment columns
#' @return Filtered tibble with complete proteins only
filter_complete_proteins <- function(tmt_data) {
    # Get the total number of unique Species-Replicate-Treatment combinations
    n_conditions <- tmt_data %>%
        distinct(Replicate, Treatment) %>%
        nrow()

    all_proteins <- tmt_data %>%
        distinct(Species, Protein_id)

    n_before <- nrow(all_proteins)

    # Keep only proteins that appear in all conditions
    filtered <- tmt_data %>%
        group_by(Species, Protein_id) %>%
        filter(n_distinct(paste(Replicate, Treatment, sep = "_")) == n_conditions) %>%
        ungroup()

    kept_proteins <- filtered %>%
        distinct(Species, Protein_id)

    n_after <- nrow(kept_proteins)
    n_removed <- n_before - n_after

    # Find removed proteins
    removed_proteins <- all_proteins %>%
        anti_join(kept_proteins, by = c("Species", "Protein_id")) %>%
        # Get additional info for removed proteins from original data
        left_join(
            tmt_data,
            by = c("Species", "Protein_id")
        )

    message("Protein completeness filtering:\n")
    message("  Total unique conditions: ", n_conditions, "\n")
    message("  Proteins before filtering: ", n_before, "\n")
    message("  Proteins after filtering: ", n_after, "\n")
    message("  Proteins removed: ", n_removed, "\n\n")

    list(filtered_data = filtered, removed_proteins = removed_proteins)
}
