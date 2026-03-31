# Process MOBILE output for one organism:
# get CTRL- and OSS-dependent interaction sets,
# extract gene lists, and compare condition-specific genes.

format_pathway <- function(x) {
  x
}

process_organism <- function(data, org, output_path) {

  # CTRL-dependent interactions
  ctrl_gained_lost <- data %>%
    filter(
      (presence_in_LOGO_CTRL == 1 & presence_in_FULL == 0) |
      (presence_in_LOGO_CTRL == 0 & presence_in_FULL == 1)
    )

  ctrl_matrix <- ctrl_gained_lost %>%
    select(
      ID_TPM,
      ID_iBAQ,
      FULLdata_coeff,
      abs,
      presence_in_FULL,
      LOGO_CTRL_coeff,
      abs.3,
      presence_in_LOGO_CTRL
    ) %>%
    mutate(abs.3 = ifelse(abs.3 == 0, abs, abs.3)) %>%
    filter(abs.3 > 0.01)

  write.table(
    ctrl_matrix,
    file = file.path(output_path, paste0(org, "_CTRL_dependent_matrix.csv")),
    sep = ";",
    row.names = FALSE,
    quote = FALSE
  )

  unique_genes_ctrl <- unique(c(ctrl_matrix$ID_TPM, ctrl_matrix$ID_iBAQ))

  write.csv(
    unique_genes_ctrl,
     file.path(output_path, paste0("unique_genes_LOGO_CTRL_", org, ".csv")),
    row.names = FALSE
  )

  # OSS-dependent interactions
  oss_gained_lost <- data %>%
    filter(
      (presence_in_LOGO_OSS == 1 & presence_in_FULL == 0) |
      (presence_in_LOGO_OSS == 0 & presence_in_FULL == 1)
    )

  oss_matrix <- oss_gained_lost %>%
    select(
      ID_TPM,
      ID_iBAQ,
      FULLdata_coeff,
      abs,
      presence_in_FULL,
      LOGO_OSS_coeff,
      abs.8,
      presence_in_LOGO_OSS
    ) %>%
    mutate(abs.8 = ifelse(abs.8 == 0, abs, abs.8)) %>%
    filter(abs.8 > 0.01)

  write.table(
    oss_matrix,
    file =  file.path(output_path, paste0(org, "_OSS_dependent_matrix.csv")),
    sep = ";",
    row.names = FALSE,
    quote = FALSE
  )

  unique_genes_oss <- unique(c(oss_matrix$ID_TPM, oss_matrix$ID_iBAQ))

  write.csv(
    unique_genes_oss,
     file.path(output_path, paste0("unique_genes_LOGO_OSS_", org, ".csv")),
    row.names = FALSE
  )

  # CTRL vs OSS comparison
  unique_genes_ctrl_only_not_oss <- setdiff(unique_genes_ctrl, unique_genes_oss)
  unique_genes_oss_only <- setdiff(unique_genes_oss, unique_genes_ctrl)

  write.csv(
    unique_genes_ctrl_only_not_oss,
     file.path(output_path, paste0("unique_genes_CTRL_only_not_in_OSS_", org, ".csv")),
    row.names = FALSE
  )

  write.csv(
    unique_genes_oss_only,
    file.path(output_path, paste0("unique_genes_OSS_not_in_CTRL_", org, ".csv")),
    row.names = FALSE
  )

  return(list(
    ctrl_matrix = ctrl_matrix,
    oss_matrix = oss_matrix,
    unique_genes_ctrl = unique_genes_ctrl,
    unique_genes_oss = unique_genes_oss,
    unique_genes_ctrl_only_not_oss = unique_genes_ctrl_only_not_oss,
    unique_genes_oss_only = unique_genes_oss_only
  ))
}

create_gmt <- function(term2gene_file, term2name_file, output_gmt_file) {
  term2gene <- read.delim(term2gene_file, header = TRUE)
  term2name <- read.delim(term2name_file, header = TRUE)

  gmt_data <- term2gene %>%
    left_join(term2name, by = "TERM") %>%
    group_by(TERM, ONTOLOGY) %>%
    summarize(Genes = paste(GENE, collapse = "\t"), .groups = "drop") %>%
    mutate(GMT_Line = paste(TERM, ONTOLOGY, Genes, sep = "\t"))

  writeLines(gmt_data$GMT_Line, output_gmt_file)

  cat("Updated GMT file saved as:", output_gmt_file, "\n")
}

runORA <- function(gene_vector, label, gmt_pathways, term2name) {
  ora <- enricher(
    gene          = gene_vector,
    TERM2GENE     = gmt_pathways,
    minGSSize     = 5,
    maxGSSize     = 500,
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.25
  )

  if (is.null(ora)) {
    message("No enrichment found for: ", label)
    return(data.frame())
  }

  ora_df <- ora %>%
    as.data.frame() %>%
    left_join(term2name, by = c("ID" = "TERM")) %>%
    mutate(
      Description = ONTOLOGY,
      Condition   = label
    ) %>%
    filter(!is.na(qvalue))

  return(ora_df)
}

combine_ORA_objects <- function(type = c("dependent", "specific")) {

  type <- match.arg(type)

  if (type == "dependent") {
    ora_list <- list(
      SE = ora_dependent_SE,
      SA = ora_dependent_SA,
      YP = ora_dependent_YP
    )
  } else {
    ora_list <- list(
      SE = ora_specific_SE,
      SA = ora_specific_SA,
      YP = ora_specific_YP
    )
  }

  combined <- bind_rows(
    lapply(names(ora_list), function(org) {
      df <- ora_list[[org]]

      if (is.null(df) || nrow(df) == 0) {
        return(NULL)
      }

      df$Dataset <- org
      df
    })
  )

  if (is.null(combined) || nrow(combined) == 0) {
    stop("No non-empty ORA objects found in the environment.")
  }

  combined <- combined %>%
    filter(!is.na(Description)) %>%
    mutate(
      pvalue   = as.numeric(pvalue),
      p.adjust = as.numeric(p.adjust),
      qvalue   = as.numeric(qvalue),
      Count    = as.numeric(Count)
    ) %>%
    mutate(
      pvalue   = ifelse(is.na(pvalue), 1, pvalue),
      p.adjust = ifelse(is.na(p.adjust), 1, p.adjust),
      qvalue   = ifelse(is.na(qvalue), 1, qvalue)
    )

  return(combined)
}

plot_ORA_dotplot <- function(combined,
                             condition = c("OSS specific", "OSS dependent",
                                           "CTRL specific", "CTRL dependent"),
                             gap = 0.35) {

  condition <- match.arg(condition)

  required_cols <- c("Condition", "Dataset", "Description", "Count", "p.adjust")
  missing_cols <- setdiff(required_cols, colnames(combined))

  if (length(missing_cols) > 0) {
    stop("Missing column(s): ", paste(missing_cols, collapse = ", "))
  }

  plot_data <- combined %>%
    dplyr::filter(Condition == condition) %>%
    dplyr::mutate(
      Dataset = trimws(as.character(Dataset)),
      Description = as.character(Description)
    )

  if (nrow(plot_data) == 0) {
    stop("No rows found for Condition = ", condition)
  }

  top_paths <- plot_data %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = 20) %>%
    dplyr::pull(Description) %>%
    unique()

  plot_data <- plot_data %>%
    dplyr::filter(Description %in% top_paths) %>%
    dplyr::mutate(
      neglog10 = -log10(p.adjust),
      pathway = format_pathway(Description)
    )

  plot_data$pathway <- factor(plot_data$pathway, levels = rev(unique(plot_data$pathway)))

  plot_data <- plot_data %>%
    dplyr::mutate(
      Dataset = trimws(as.character(Dataset)),
      x = dplyr::case_when(
        Dataset == "SE" ~ 1,
        Dataset == "YP" ~ 1 + gap,
        TRUE ~ NA_real_
      )
    ) %>%
    dplyr::filter(!is.na(x))

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = x, y = pathway, size = Count, fill = neglog10)
  ) +
    ggplot2::geom_point(shape = 21, color = "black", alpha = 0.85) +
    ggplot2::scale_x_continuous(
      breaks = c(1, 1 + gap),
      labels = c("SE", "YP"),
      limits = c(1 - gap/2, 1 + gap + gap/2)
    ) +
    ggplot2::scale_fill_gradient(
      low = "#a72604",
      high = "#fedad0",
      name = "-log10(FDR)"
    ) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = paste(condition, "ORA Dotplot"),
      x = "Dataset",
      y = "Pathway",
      size = "Gene Count"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, size = 10),
      axis.text.y = ggplot2::element_text(size = 8),
      plot.title  = ggplot2::element_text(hjust = 0.5, face = "bold"),
      panel.grid = ggplot2::element_blank()
    )

  return(p)
}