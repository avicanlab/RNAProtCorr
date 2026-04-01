# ============================================================================
# ORA Analysis — MOBILE Framework Output
# ============================================================================
#
# Description:
#   Processes MOBILE framework output to identify CTRL- and OSS-dependent
#   mRNA-protein interaction sets, performs Over-Representation Analysis (ORA)
#   using clusterProfiler, and generates dotplot visualizations.
#
# Expected inputs (in ORA directory):
#   - connections_ALL_<ORG>.csv  : MOBILE output per organism
#   - <ORG>_TERM2GENE.txt        : Term-to-gene mapping
#   - <ORG>_TERM2NAME.txt        : Term-to-name mapping
#
# Outputs (in ORA_OUTPUT directory):
#   - <ORG>_CTRL/OSS_dependent_matrix.csv : Interaction matrices
#   - unique_genes_LOGO_*_<ORG>.csv        : Gene lists per condition
#   - Updated_<ORG>_geneset.gmt            : GMT gene set files
#   - ORA_dependent/specific_*.csv         : ORA results per organism
#   - ORA_*_combined_CTRL_OSS.csv          : Combined ORA results
#   - ORA_*_dotplot.png                    : Dotplot figures
#
# Requirements:
#   library(clusterProfiler)
#   library(tidyverse)
# ============================================================================

# Organism abbreviation to full species name mapping
ORG_SPECIES_MAP <- c(
  SE = "Salmonella_enterica",
  SA = "Staphylococcus_aureus",
  YP = "Yersinia_pseudotuberculosis"
)

# Conditions to analyse
ORA_CONDITIONS <- c("CTRL dependent", "CTRL specific",
                    "OSS dependent", "OSS specific")

# Minimum absolute coefficient threshold for interaction filtering
ABS_COEFF_THRESHOLD <- 0.01

# ORA parameters
ORA_MIN_GS <- 5
ORA_MAX_GS <- 500
ORA_PVAL <- 0.05
ORA_QVAL <- 0.25

# Dotplot parameters
DOTPLOT_TOP_N <- 20
DOTPLOT_WIDTH <- 5
DOTPLOT_HEIGHT <- 8

# ============================================================================
# DATA PROCESSING
# ============================================================================

#' Extract CTRL- and OSS-Dependent Interaction Matrices for One Organism
#'
#' @description
#' Filters MOBILE output to identify interactions gained or lost in the
#' CTRL- and OSS-specific runs (LOGO) compared to the full data run.
#' Writes interaction matrices and gene lists to disk.
#'
#' @param data        Tibble. MOBILE connections output with columns:
#'   presence_in_LOGO_CTRL, presence_in_LOGO_OSS, presence_in_FULL,
#'   ID_TPM, ID_iBAQ, FULLdata_coeff, abs, LOGO_CTRL_coeff, abs.3,
#'   LOGO_OSS_coeff, abs.8.
#' @param org         Chr. Organism abbreviation (e.g. "SE", "SA", "YP").
#' @param output_path Chr. Directory for saving output files.
#'
#' @return Named list:
#'   ctrl_matrix, oss_matrix,
#'   unique_genes_ctrl, unique_genes_oss,
#'   unique_genes_ctrl_only, unique_genes_oss_only
process_organism <- function(data, org, output_path) {
  message("Processing organism: ", org)

  # Helper: extract gained/lost interactions for one condition
  extract_condition_matrix <- function(presence_col, coeff_col, abs_col) {
    data %>%
      filter(
        (.data[[presence_col]] == 1 & presence_in_FULL == 0) |
          (.data[[presence_col]] == 0 & presence_in_FULL == 1)
      ) %>%
      dplyr::select(
        ID_TPM, ID_iBAQ,
        FULLdata_coeff, abs,
        presence_in_FULL,
        all_of(c(coeff_col, abs_col, presence_col))
      ) %>%
      mutate(
        # Use full-data abs if condition-specific abs is zero
        !!abs_col := ifelse(.data[[abs_col]] == 0, abs, .data[[abs_col]])
      ) %>%
      filter(.data[[abs_col]] > ABS_COEFF_THRESHOLD)
  }

  ctrl_matrix <- extract_condition_matrix(
    "presence_in_LOGO_CTRL", "LOGO_CTRL_coeff", "abs.3"
  )
  oss_matrix <- extract_condition_matrix(
    "presence_in_LOGO_OSS", "LOGO_OSS_coeff", "abs.8"
  )

  unique_genes_ctrl <- unique(c(ctrl_matrix$ID_TPM, ctrl_matrix$ID_iBAQ))
  unique_genes_oss <- unique(c(oss_matrix$ID_TPM, oss_matrix$ID_iBAQ))

  # Condition-specific gene sets (present in one condition but not the other)
  unique_genes_ctrl_only <- setdiff(unique_genes_ctrl, unique_genes_oss)
  unique_genes_oss_only <- setdiff(unique_genes_oss, unique_genes_ctrl)

  # Write all outputs
  write_outputs <- list(
    list(data = ctrl_matrix, file = paste0(org, "_CTRL_dependent_matrix.csv"), sep = ";"),
    list(data = oss_matrix, file = paste0(org, "_OSS_dependent_matrix.csv"), sep = ";"),
    list(data = unique_genes_ctrl, file = paste0("unique_genes_LOGO_CTRL_", org, ".csv"), sep = ","),
    list(data = unique_genes_oss, file = paste0("unique_genes_LOGO_OSS_", org, ".csv"), sep = ","),
    list(data = unique_genes_ctrl_only, file = paste0("unique_genes_CTRL_only_not_in_OSS_", org, ".csv"), sep = ","),
    list(data = unique_genes_oss_only, file = paste0("unique_genes_OSS_not_in_CTRL_", org, ".csv"), sep = ",")
  )

  walk(write_outputs, function(x) {
    write.table(x$data, file.path(output_path, x$file),
                sep = x$sep, row.names = FALSE, quote = FALSE)
  })

  message("  CTRL-dependent genes: ", length(unique_genes_ctrl))
  message("  OSS-dependent genes:  ", length(unique_genes_oss))
  message("  CTRL-specific (not OSS): ", length(unique_genes_ctrl_only))
  message("  OSS-specific (not CTRL): ", length(unique_genes_oss_only))

  list(
    ctrl_matrix = ctrl_matrix,
    oss_matrix = oss_matrix,
    unique_genes_ctrl = unique_genes_ctrl,
    unique_genes_oss = unique_genes_oss,
    unique_genes_ctrl_only = unique_genes_ctrl_only,
    unique_genes_oss_only = unique_genes_oss_only
  )
}

# ============================================================================
# GMT FILE CREATION
# ============================================================================

#' Build GMT Gene Set File from TERM2GENE and TERM2NAME Files
#'
#' @description
#' Joins term-to-gene and term-to-name mappings, then writes a GMT-format
#' file where each line is: TERM <TAB> ONTOLOGY <TAB> gene1 <TAB> gene2 ...
#'
#' @param term2gene_file   Chr. Path to TERM2GENE TSV (columns: TERM, GENE).
#' @param term2name_file   Chr. Path to TERM2NAME TSV (columns: TERM, ONTOLOGY).
#' @param output_gmt_file  Chr. Path to write the GMT output file.
#'
#' @return NULL invisibly.
create_gmt <- function(term2gene_file, term2name_file, output_gmt_file) {
  term2gene <- read.delim(term2gene_file, header = TRUE)
  term2name <- read.delim(term2name_file, header = TRUE)

  gmt_lines <- term2gene %>%
    left_join(term2name, by = "TERM") %>%
    group_by(TERM, ONTOLOGY) %>%
    summarise(genes = paste(GENE, collapse = "\t"), .groups = "drop") %>%
    mutate(line = paste(TERM, ONTOLOGY, genes, sep = "\t")) %>%
    pull(line)

  writeLines(gmt_lines, output_gmt_file)
  message("GMT file saved: ", output_gmt_file,
          " (", length(gmt_lines), " gene sets)")

  invisible(NULL)
}

# ============================================================================
# ORA
# ============================================================================

#' Run Over-Representation Analysis for One Gene Set
#'
#' @description
#' Wraps clusterProfiler::enricher() with term2name annotation joining.
#' Returns empty data frame if no enrichment is found.
#'
#' @param gene_vector  Chr vector. Genes to test.
#' @param label        Chr. Condition label added to output (e.g. "CTRL dependent").
#' @param gmt_pathways Tibble. Two-column TERM2GENE table from read.gmt().
#' @param term2name    Tibble. Two-column table: TERM, ONTOLOGY.
#'
#' @return Tibble with ORA results and Condition column, or empty tibble.
run_ora <- function(gene_vector, label, gmt_pathways, term2name) {
  if (length(gene_vector) == 0) {
    message("  No genes for condition: ", label, " — skipping.")
    return(tibble())
  }

  ora <- tryCatch(
    enricher(
      gene = gene_vector,
      TERM2GENE = gmt_pathways,
      minGSSize = ORA_MIN_GS,
      maxGSSize = ORA_MAX_GS,
      pvalueCutoff = ORA_PVAL,
      qvalueCutoff = ORA_QVAL
    ),
    error = function(e) { message("  enricher error: ", e$message); NULL }
  )

  if (is.null(ora) || nrow(ora@result) == 0) {
    message("  No enrichment found for: ", label)
    return(tibble())
  }

  ora %>%
    as.data.frame() %>%
    as_tibble() %>%
    left_join(term2name, by = c("ID" = "TERM")) %>%
    mutate(
      Description = ONTOLOGY,
      Condition = label
    ) %>%
    filter(!is.na(qvalue))
}

#' Run ORA for CTRL and OSS Conditions for One Organism
#'
#' @description
#' Runs four ORA analyses per organism:
#'   - CTRL dependent / OSS dependent (all genes in each condition)
#'   - CTRL specific  / OSS specific  (genes unique to each condition)
#' Writes results to disk and returns both result tibbles.
#'
#' @param org         Chr. Organism abbreviation.
#' @param results     List. Output of process_organism() for this organism.
#' @param gmt_pathways Tibble. Gene sets from read.gmt().
#' @param term2name   Tibble. TERM to ONTOLOGY mapping.
#' @param output_path Chr. Output directory.
#'
#' @return Named list: dependent, specific
run_ora_organism <- function(org, results, gmt_pathways, term2name, output_path) {
  message("\nRunning ORA for: ", org)

  ora_ctrl_dep <- run_ora(results$unique_genes_ctrl, "CTRL dependent", gmt_pathways, term2name)
  ora_oss_dep <- run_ora(results$unique_genes_oss, "OSS dependent", gmt_pathways, term2name)
  ora_ctrl_spec <- run_ora(results$unique_genes_ctrl_only, "CTRL specific", gmt_pathways, term2name)
  ora_oss_spec <- run_ora(results$unique_genes_oss_only, "OSS specific", gmt_pathways, term2name)

  dependent <- bind_rows(ora_ctrl_dep, ora_oss_dep)
  specific <- bind_rows(ora_ctrl_spec, ora_oss_spec) %>%
    filter(!is.na(Description))

  write.csv(dependent, file.path(output_path, paste0("ORA_dependent_CTRL_OSS_", org, ".csv")), row.names = FALSE)
  write.csv(specific, file.path(output_path, paste0("ORA_specific_CTRL_OSS_", org, ".csv")), row.names = FALSE)

  list(dependent = dependent, specific = specific)
}

#' Combine ORA Results Across All Organisms
#'
#' @description
#' Row-binds per-organism ORA results, adds Dataset column, and coerces
#' numeric columns that may have been read as character.
#'
#' @param ora_list Named list. Each element is a tibble of ORA results,
#'   named by organism abbreviation (e.g. list(SE = ..., SA = ..., YP = ...)).
#'
#' @return Combined tibble with Dataset column.
combine_ora_results <- function(ora_list) {
  combined <- map_dfr(names(ora_list), function(org) {
    df <- ora_list[[org]]
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df %>% mutate(Dataset = org)
  })

  if (nrow(combined) == 0) stop("No non-empty ORA results to combine.")

  combined %>%
    filter(!is.na(Description)) %>%
    mutate(across(c(pvalue, p.adjust, qvalue, Count), as.numeric)) %>%
    mutate(across(c(pvalue, p.adjust, qvalue), ~replace_na(.x, 1)))
}

# ============================================================================
# VISUALIZATION
# ============================================================================

#' Plot ORA Dot Plot for One Condition
#'
#' @description
#' Dot plot showing top enriched pathways for one condition (e.g. "OSS specific")
#' across two organisms (SE and YP). Dot size = gene count,
#' dot colour = -log10(FDR). Top 20 pathways by adjusted p-value are shown.
#'
#' @param combined   Tibble. Output of combine_ora_results().
#' @param condition  Chr. One of: "OSS specific", "OSS dependent",
#'   "CTRL specific", "CTRL dependent".
#' @param gap        Num. Horizontal gap between organism columns. Default 0.35.
#' @param top_n      Int. Number of top pathways to show. Default DOTPLOT_TOP_N.
#'
#' @return ggplot object.
plot_ora_dotplot <- function(
  combined,
  condition = c("OSS specific", "OSS dependent",
                "CTRL specific", "CTRL dependent"),
  gap = 0.35,
  top_n = DOTPLOT_TOP_N
) {
  condition <- match.arg(condition)

  required_cols <- c("Condition", "Dataset", "Description", "Count", "p.adjust")
  missing_cols <- setdiff(required_cols, colnames(combined))
  if (length(missing_cols) > 0) {
    stop("Missing column(s): ", paste(missing_cols, collapse = ", "))
  }

  plot_data <- combined %>%
    filter(Condition == condition) %>%
    mutate(
      Dataset = trimws(as.character(Dataset)),
      Description = as.character(Description)
    )

  if (nrow(plot_data) == 0) stop("No rows found for Condition = ", condition)

  # Keep top_n pathways by best adjusted p-value across organisms
  top_paths <- plot_data %>%
    arrange(p.adjust) %>%
    distinct(Description) %>%
    slice_head(n = top_n) %>%
    pull(Description)

  plot_data <- plot_data %>%
    filter(Description %in% top_paths) %>%
    mutate(
      neglog10 = -log10(p.adjust),
      pathway = factor(Description, levels = rev(top_paths)),
      x = case_when(
        Dataset == "SE" ~ 1,
        Dataset == "YP" ~ 1 + gap,
        TRUE ~ NA_real_
      )
    ) %>%
    filter(!is.na(x))

  ggplot(plot_data, aes(x = x, y = pathway, size = Count, fill = neglog10)) +
    geom_point(shape = 21, colour = "black", alpha = 0.85) +
    scale_x_continuous(
      breaks = c(1, 1 + gap),
      labels = c("SE", "YP"),
      limits = c(1 - gap / 2, 1 + gap + gap / 2)
    ) +
    scale_fill_gradient(
      low = "#fedad0",
      high = "#a72604",
      name = "-log10(FDR)"
    ) +
    scale_size_continuous(name = "Gene Count") +
    labs(
      title = paste(condition, "— ORA"),
      x = "Dataset",
      y = "Pathway"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
      axis.text.y = element_text(size = 8),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid = element_blank()
    )
}