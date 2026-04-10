# ORA Analysis — MOBILE Framework Output --------------------------------------
#
# Description:
#   Processes MOBILE framework output to identify CTRL- and OSS-dependent
#   mRNA-protein interaction sets, performs Over-Representation Analysis (ORA)
#   using clusterProfiler, and generates dot-plot visualisations.
#
# Expected inputs (in the ORA directory):
#   - connections_ALL_<ORG>.csv  : MOBILE output per organism
#   - <ORG>_TERM2GENE.txt        : Term-to-gene mapping
#   - <ORG>_TERM2NAME.txt        : Term-to-name mapping
#
# Outputs (in ORA_OUTPUT directory):
#   - <ORG>_CTRL/OSS_dependent_matrix.csv
#   - unique_genes_LOGO_*_<ORG>.csv
#   - Updated_<ORG>_geneset.gmt
#   - ORA_dependent/specific_*.csv
#   - ORA_*_combined_CTRL_OSS.csv
#   - ORA_*_dotplot.png
#
# Requirements:
#   library(clusterProfiler)
#   library(dplyr)
#   library(purrr)
#   library(ggplot2)

# ── Constants ──────────────────────────────────────────────────────────────────

#' Mapping from organism abbreviation to full species name
ORG_SPECIES_MAP <- c(
  SE = "Salmonella_enterica",
  SA = "Staphylococcus_aureus",
  YP = "Yersinia_pseudotuberculosis"
)

#' ORA conditions to analyse
ORA_CONDITIONS <- c("CTRL dependent", "CTRL specific", "OSS dependent", "OSS specific")

#' Minimum absolute coefficient threshold for interaction filtering
ABS_COEFF_THRESHOLD <- 0.2

#' ORA parameters
ORA_MIN_GS <- 5
ORA_MAX_GS <- 500
ORA_PVAL <- 0.05
ORA_QVAL <- 0.25

#' Dot-plot display parameters
DOTPLOT_TOP_N <- 20
DOTPLOT_WIDTH <- 5
DOTPLOT_HEIGHT <- 8

# ── Data processing ────────────────────────────────────────────────────────────

#' Extract CTRL- and OSS-Dependent Interaction Matrices for One Organism
#'
#' @description
#' Filters MOBILE output to identify interactions gained or lost in the
#' CTRL- and OSS-specific leave-one-group-out (LOGO) runs compared to the
#' full-data run. Writes interaction matrices and gene lists to disk.
#'
#' @param data        Tibble. MOBILE connections output with columns:
#'   `presence_in_LOGO_CTRL`, `presence_in_LOGO_OSS`, `presence_in_FULL`,
#'   `ID_TPM`, `ID_iBAQ`, `FULLdata_coeff`, `abs`, `LOGO_CTRL_coeff`,
#'   `abs.3`, `LOGO_OSS_coeff`, `abs.8`.
#' @param org         Character. Organism abbreviation (e.g. `"SE"`).
#' @param output_path Character. Directory for saving output files.
#'
#' @return Named list:
#'   - `ctrl_matrix`           : tibble of CTRL-dependent interactions
#'   - `oss_matrix`            : tibble of OSS-dependent interactions
#'   - `unique_genes_ctrl`     : character vector of CTRL-dependent gene IDs
#'   - `unique_genes_oss`      : character vector of OSS-dependent gene IDs
#'   - `unique_genes_ctrl_only`: genes unique to CTRL (not in OSS)
#'   - `unique_genes_oss_only` : genes unique to OSS (not in CTRL)
process_organism <- function (data, org, output_path) {
  message("Processing organism: ", org)

  # Internal helper: extract gained/lost interactions for one LOGO condition
  extract_condition_matrix <- function (presence_col, coeff_col, abs_col) {
    data |>
      dplyr::filter(
        (.data[[presence_col]] == 1 & presence_in_FULL == 0) |
          (.data[[presence_col]] == 0 & presence_in_FULL == 1)
      ) |>
      dplyr::select(
        ID_TPM, ID_iBAQ,
        FULLdata_coeff, abs,
        presence_in_FULL,
        dplyr::all_of(c(coeff_col, abs_col, presence_col))
      ) |>
      dplyr::mutate(
        # Fall back to full-data abs when condition-specific abs is zero
        !!abs_col := ifelse(.data[[abs_col]] == 0, abs, .data[[abs_col]])
      ) |>
      dplyr::filter(.data[[abs_col]] > ABS_COEFF_THRESHOLD)
  }

  ctrl_matrix <- extract_condition_matrix("presence_in_LOGO_CTRL", "LOGO_CTRL_coeff", "abs.3")
  oss_matrix <- extract_condition_matrix("presence_in_LOGO_OSS", "LOGO_OSS_coeff", "abs.8")

  unique_genes_ctrl <- unique(c(ctrl_matrix$ID_TPM, ctrl_matrix$ID_iBAQ))
  unique_genes_oss <- unique(c(oss_matrix$ID_TPM, oss_matrix$ID_iBAQ))
  unique_genes_ctrl_only <- setdiff(unique_genes_ctrl, unique_genes_oss)
  unique_genes_oss_only <- setdiff(unique_genes_oss, unique_genes_ctrl)

  write_outputs <- list(
    list(data = ctrl_matrix, file = paste0(org, "_CTRL_dependent_matrix.csv"), sep = ";"),
    list(data = oss_matrix, file = paste0(org, "_OSS_dependent_matrix.csv"), sep = ";"),
    list(data = unique_genes_ctrl, file = paste0("unique_genes_LOGO_CTRL_", org, ".csv"), sep = ","),
    list(data = unique_genes_oss, file = paste0("unique_genes_LOGO_OSS_", org, ".csv"), sep = ","),
    list(data = unique_genes_ctrl_only, file = paste0("unique_genes_CTRL_only_not_in_OSS_", org, ".csv"), sep = ","),
    list(data = unique_genes_oss_only, file = paste0("unique_genes_OSS_not_in_CTRL_", org, ".csv"), sep = ",")
  )

  purrr::walk(write_outputs, function (x) {
    utils::write.table(x$data, file.path(output_path, x$file),
                       sep = x$sep, row.names = FALSE, quote = FALSE)
  })

  message("  CTRL-dependent genes        : ", length(unique_genes_ctrl))
  message("  OSS-dependent genes         : ", length(unique_genes_oss))
  message("  CTRL-specific (not OSS)     : ", length(unique_genes_ctrl_only))
  message("  OSS-specific  (not CTRL)    : ", length(unique_genes_oss_only))

  list(
    ctrl_matrix = ctrl_matrix,
    oss_matrix = oss_matrix,
    unique_genes_ctrl = unique_genes_ctrl,
    unique_genes_oss = unique_genes_oss,
    unique_genes_ctrl_only = unique_genes_ctrl_only,
    unique_genes_oss_only = unique_genes_oss_only
  )
}

# ── GMT file creation ──────────────────────────────────────────────────────────

#' Build a GMT Gene-Set File from TERM2GENE and TERM2NAME Files
#'
#' @description
#' Joins term-to-gene and term-to-name mappings, then writes a GMT-format
#' file where each line is: `TERM <TAB> ONTOLOGY <TAB> gene1 <TAB> gene2 ...`
#'
#' @param term2gene_file  Character. Path to TERM2GENE TSV (columns: `TERM`,
#'   `GENE`).
#' @param term2name_file  Character. Path to TERM2NAME TSV (columns: `TERM`,
#'   `ONTOLOGY`).
#' @param output_gmt_file Character. Path for the GMT output file.
#'
#' @return NULL invisibly.
create_gmt <- function (term2gene_file, term2name_file, output_gmt_file) {
  term2gene <- utils::read.delim(term2gene_file, header = TRUE)
  term2name <- utils::read.delim(term2name_file, header = TRUE)

  gmt_lines <- term2gene |>
    dplyr::left_join(term2name, by = "TERM") |>
    dplyr::group_by(TERM, ONTOLOGY) |>
    dplyr::summarise(genes = paste(GENE, collapse = "\t"), .groups = "drop") |>
    dplyr::mutate(line = paste(TERM, ONTOLOGY, genes, sep = "\t")) |>
    dplyr::pull(line)

  writeLines(gmt_lines, output_gmt_file)
  message("GMT file saved: ", output_gmt_file, " (", length(gmt_lines), " gene sets)")

  invisible(NULL)
}

# ── ORA ────────────────────────────────────────────────────────────────────────

#' Run Over-Representation Analysis for One Gene Set
#'
#' @description
#' Wraps `clusterProfiler::enricher()` and joins term-name annotation.
#' Returns an empty tibble if no enrichment is found.
#'
#' @param gene_vector  Character vector. Genes to test.
#' @param label        Character. Condition label added to the output
#'   (e.g. `"CTRL dependent"`).
#' @param gmt_pathways Tibble. TERM2GENE table from `clusterProfiler::read.gmt()`.
#' @param term2name    Tibble. Two-column table: `TERM`, `ONTOLOGY`.
#'
#' @return Tibble with ORA results and a `Condition` column, or an empty
#'   tibble when no enrichment is found.
run_ora <- function (gene_vector, label, gmt_pathways, term2name) {
  if (length(gene_vector) == 0) {
    message("  No genes for condition: ", label, " — skipping.")
    return(tibble::tibble())
  }

  ora <- tryCatch(
    clusterProfiler::enricher(
      gene = gene_vector,
      TERM2GENE = gmt_pathways,
      minGSSize = ORA_MIN_GS,
      maxGSSize = ORA_MAX_GS,
      pvalueCutoff = ORA_PVAL,
      qvalueCutoff = ORA_QVAL
    ),
    error = function (e) { message("  enricher error: ", e$message); NULL }
  )

  if (is.null(ora) || nrow(ora@result) == 0) {
    message("  No enrichment found for: ", label)
    return(tibble::tibble())
  }

  ora |>
    as.data.frame() |>
    tibble::as_tibble() |>
    dplyr::left_join(term2name, by = c("ID" = "TERM")) |>
    dplyr::mutate(
      Description = ONTOLOGY,
      Condition = label
    ) |>
    dplyr::filter(!is.na(qvalue))
}

#' Run ORA for CTRL and OSS Conditions for One Organism
#'
#' @description
#' Runs four ORA analyses per organism (CTRL / OSS × dependent / specific),
#' writes results to disk, and returns both result tibbles.
#'
#' @param org          Character. Organism abbreviation.
#' @param results      List. Output of `process_organism()` for this organism.
#' @param gmt_pathways Tibble. Gene sets from `clusterProfiler::read.gmt()`.
#' @param term2name    Tibble. TERM to ONTOLOGY mapping.
#' @param output_path  Character. Output directory.
#'
#' @return Named list: `dependent` and `specific` (tibbles of ORA results).
run_ora_organism <- function (org, results, gmt_pathways, term2name, output_path) {
  message("\nRunning ORA for: ", org)

  ora_ctrl_dep <- run_ora(results$unique_genes_ctrl, "CTRL dependent", gmt_pathways, term2name)
  ora_oss_dep <- run_ora(results$unique_genes_oss, "OSS dependent", gmt_pathways, term2name)
  ora_ctrl_spec <- run_ora(results$unique_genes_ctrl_only, "CTRL specific", gmt_pathways, term2name)
  ora_oss_spec <- run_ora(results$unique_genes_oss_only, "OSS specific", gmt_pathways, term2name)

  dependent <- dplyr::bind_rows(ora_ctrl_dep, ora_oss_dep)
  specific <- dplyr::bind_rows(ora_ctrl_spec, ora_oss_spec) |>
    dplyr::filter(!is.na(Description))

  utils::write.csv(dependent, file.path(output_path, paste0("ORA_dependent_CTRL_OSS_", org, ".csv")), row.names = FALSE)
  utils::write.csv(specific, file.path(output_path, paste0("ORA_specific_CTRL_OSS_", org, ".csv")), row.names = FALSE)

  list(dependent = dependent, specific = specific)
}

#' Combine ORA Results Across All Organisms
#'
#' @description
#' Row-binds per-organism ORA result tibbles, adds a `Dataset` column, and
#' coerces numeric columns that may have been read as character.
#'
#' @param ora_list Named list. Each element is a tibble of ORA results, named
#'   by organism abbreviation (e.g. `list(SE = ..., SA = ..., YP = ...)`).
#'
#' @return Combined tibble with a `Dataset` column.
combine_ora_results <- function (ora_list) {
  combined <- purrr::map_dfr(names(ora_list), function (org) {
    df <- ora_list[[org]]
    if (is.null(df) || nrow(df) == 0) return(NULL)
    dplyr::mutate(df, Dataset = org)
  })

  if (nrow(combined) == 0) stop("No non-empty ORA results to combine.")

  combined |>
    dplyr::filter(!is.na(Description)) |>
    dplyr::mutate(dplyr::across(c(pvalue, p.adjust, qvalue, Count), as.numeric)) |>
    dplyr::mutate(dplyr::across(c(pvalue, p.adjust, qvalue), ~tidyr::replace_na(.x, 1)))
}

# ── Visualisation ──────────────────────────────────────────────────────────────

#' Plot ORA Dot Plot for One Condition
#'
#' @description
#' Dot plot showing the top enriched pathways for one ORA condition across two
#' organisms (SE and YP). Dot size encodes gene count; dot colour encodes
#' −log10(FDR). The top `top_n` pathways by adjusted p-value are shown.
#'
#' @param combined   Tibble. Output of `combine_ora_results()`.
#' @param condition  Character. One of the `ORA_CONDITIONS` strings.
#' @param gap        Numeric. Horizontal gap between organism columns.
#'   Default 0.35.
#' @param top_n      Integer. Number of top pathways to display.
#'   Default `DOTPLOT_TOP_N`.
#'
#' @return ggplot object.
plot_ora_dotplot <- function (
  combined,
  condition = ORA_CONDITIONS,
  gap = 0.35,
  top_n = DOTPLOT_TOP_N
) {
  condition <- match.arg(condition)
  # term2name <- readxl::read_excel(
  #   file.path(ENRICHMENT, "GO_terms_from_database.xlsx")
  # ) |>
  #   dplyr::filter(Obsolete == FALSE, namespace == "Biological Process")
  # combined <- combined |>
  #   dplyr::filter(ID %in% term2name$ID)

  required_cols <- c("Condition", "Dataset", "Description", "Count", "p.adjust")
  missing_cols <- setdiff(required_cols, colnames(combined))
  if (length(missing_cols) > 0) {
    stop("Missing column(s): ", paste(missing_cols, collapse = ", "))
  }

  plot_data <- combined |>
    dplyr::filter(Condition == condition) |>
    dplyr::mutate(
      Dataset = trimws(as.character(Dataset)),
      Description = as.character(Description)
    )

  if (nrow(plot_data) == 0) stop("No rows found for Condition = ", condition)

  # Keep top_n pathways by best adjusted p-value across organisms
  top_paths <- plot_data |>
    dplyr::arrange(p.adjust) |>
    dplyr::distinct(Description) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::pull(Description)

  plot_data <- plot_data |>
    dplyr::filter(Description %in% top_paths) |>
    dplyr::mutate(
      neglog10 = -log10(p.adjust),
      pathway = factor(Description, levels = rev(top_paths)),
      Dataset = factor(Dataset, levels = c("SE", "SA", "YP"))
    )

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Dataset, y = pathway, size = Count, fill = neglog10)) +
    ggplot2::geom_point(shape = 21, colour = "black", alpha = 0.85) +
    ggplot2::scale_x_discrete(
      limits = c("SE", "SA", "YP"),
      drop = FALSE
    ) +
    ggplot2::scale_fill_gradient(
      low = "#fedad0",
      high = "#a72604",
      name = "-log10(FDR)"
    ) +
    ggplot2::scale_size_continuous(name = "Gene Count") +
    ggplot2::labs(
      title = paste(condition, "\u2014 ORA"),
      x = "Dataset",
      y = "Pathway"
    ) +
    theme_publication() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      panel.grid = ggplot2::element_blank()
    )
}
