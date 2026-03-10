#' Compute Detected CDS Counts and Percentages per Species
#'
#' @description
#' For each species in annotation_data, counts the number of detected CDS
#' (unique Protein_id) in RNA and proteomics datasets relative to the total
#' CDS in the annotation reference.
#'
#' @param annotation_data Named list. Output of the annotation loading block.
#'   Each element is named by species and contains a tibble with at least
#'   a Protein_id column. Protein_id must already be resolved via
#'   New_locus_tag -> Old_locus_tag fallback before passing here.
#' @param rna_data        Tibble. Detected genes with at least columns Species
#'   and Protein_id. Typically the long-format TPM data from read_tpm().
#' @param prot_data       Tibble. Detected proteins with at least columns
#'   Species and Protein_id. Typically the filtered proteomics data from
#'   protein_quantification.R.
#'
#' @return Tibble with one row per species and columns:
#'   - Species       : species name
#'   - total_cds     : total annotated CDS in reference
#'   - rna_detected  : number of CDS detected in RNA data
#'   - prot_detected : number of CDS detected in proteomics data
#'   - pct_rna       : percentage of CDS detected in RNA
#'   - pct_prot      : percentage of CDS detected in proteomics
get_detected_cds <- function(annotation_data, rna_data, prot_data) {
    map_dfr(names(annotation_data), function(sp) {
        total_cds <- annotation_data[[sp]] %>%
            distinct(Protein_id) %>%
            nrow()

        rna_detected <- rna_data[[sp]] %>%
            distinct(Protein_id) %>%
            nrow()

        prot_detected <- prot_data[[sp]] %>%
            distinct(Protein_id) %>%
            nrow()

        tibble(
            Species       = sp,
            total_cds     = total_cds,
            rna_detected  = rna_detected,
            prot_detected = prot_detected,
            pct_rna       = 100 * rna_detected / total_cds,
            pct_prot      = 100 * prot_detected / total_cds
        )
    })
}

#' Plot Percentage of Detected CDS per Species
#'
#' @description
#' Generates a bar chart showing the percentage of annotated CDS detected
#' in RNA and proteomics data per species. Bars can be displayed side by
#' side or overlapping depending on the position argument.
#'
#' In overlay mode, the RNA bar is drawn wider than the protein bar so both
#' remain clearly visible even when protein coverage is substantially lower
#' than RNA coverage, which is typically the case.
#'
#' @param detected_cds_data Tibble. Output of get_detected_cds(). Must contain
#'   columns: Species, pct_rna, pct_prot.
#' @param position          Chr. Bar positioning mode:
#'   - "dodge"   : bars placed side by side (default)
#'   - "overlay" : bars overlapping, RNA behind and wider, protein in front
#'     and narrower
#'
#' @return A ggplot object. Save with save_plot() at the desired dimensions.
plot_barplot_detected_cds <- function(
  detected_cds_data, position = "overlay", output_path
) {
    plot_df <- detected_cds_data %>%
        pivot_longer(
            cols      = c(pct_rna, pct_prot),
            names_to  = "type",
            values_to = "pct"
        ) %>%
        mutate(
            type = factor(type,
                levels = c("pct_rna", "pct_prot"),
                labels = c("Transcriptome", "Proteome")
            ),
            Species = factor(Species)
        )

    species_levels <- levels(plot_df$Species)
    x_labels <- setNames(
        lapply(species_levels, function(s) bquote(italic(.(gsub("_", " ", s))))),
        species_levels
    )

    if (position == "overlay") {
        p <- ggplot(plot_df, aes(x = Species, y = pct, fill = type)) +
            geom_col(
                data     = ~ filter(.x, type == "Transcriptome"),
                width    = 0.6,
                alpha    = 0.5,
                position = "identity"
            ) +
            geom_col(
                data     = ~ filter(.x, type == "Proteome"),
                width    = 0.6,
                alpha    = 0.5,
                position = "identity"
            )
    } else {
        p <- ggplot(plot_df, aes(x = Species, y = pct, fill = type)) +
            geom_col(
                position = position_dodge(width = 0.75),
                width    = 0.65,
                alpha    = 0.85
            )
    }

    p <- p +
        scale_fill_manual(
            values = c("Transcriptome" = "#89CFF0", "Proteome" = "#f5dd23"),
            name   = NULL
        ) +
        scale_x_discrete(labels = x_labels) +
        scale_y_continuous(
            limits = c(0, 100),
            expand = expansion(mult = c(0, 0.02)),
            breaks = seq(0, 100, by = 20)
        ) +
        labs(x = NULL, y = "CDS detected (%)") +
        theme_bw(base_size = 11) +
        theme(
            axis.text.x = element_text(size = 9, angle = 15, hjust = 1),
            axis.text.y = element_text(size = 9),
            axis.title.y = element_text(size = 10),
            panel.grid.major.x = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "grey70"),
            legend.justification = "top",
            legend.border = element_rect(colour = "grey70"),
            legend.text = element_text(size = 9),
            legend.key.size = unit(0.4, "cm")
        )

    print(p)

    save_plot(
        plot     = p,
        filepath = output_path,
        width    = 5,
        height   = 5
    )
}
