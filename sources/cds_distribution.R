#' Compute Detected CDS Counts and Percentages per Species
#'
#' @description
#' For each species in annotation_data, counts the number of detected CDS
#' (unique Protein_id) in RNA and proteomics datasets relative to the total
#' CDS in the annotation reference.
#'
#' @param annotation_data Tibble. Output of the annotation loading block.
#'   Must contain at least columns Species and Protein_id.
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
    species_list <- unique(annotation_data$Species)

    print(annotation_data)
    print(rna_data)
    print(prot_data %>% select(c("Species", "Replicate", "Protein ID", "Total Intensity")))
    map_dfr(species_list, function(sp) {
        total_cds <- annotation_data %>%
            filter(Species == sp) %>%
            distinct(Protein_id) %>%
            nrow()

        rna_detected <- rna_data %>%
            filter(Species == sp) %>%
            distinct(Protein_id) %>%
            nrow()

        prot_detected <- prot_data %>%
            filter(Species == sp) %>%
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
            values = c("Transcriptome" = "#C2E3F2", "Proteome" = "#FFD892"),
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

#' Extract Common and Unique Protein IDs Between TPM and Proteomics
#'
#' @description
#' For each species, splits genes into three groups matching the figure:
#'   - mRNA        : all genes detected in TPM (full set)
#'   - Protein     : genes detected in both TPM and proteomics
#'   - Not protein : genes detected in TPM but absent from proteomics
#'
#' @param tpm_data   Tibble. Must contain columns Species, Protein_id,
#'   and mean_log2_TPM.
#' @param prot_data  Tibble. Must contain columns Species and Protein_id.
#' @param sp         Chr. Species name to filter on.
#'
#' @return Tibble with columns: Protein_id, mean_log2_TPM, group
#'   where group is a factor with levels: "mRNA", "Not protein", "Protein"
classify_tpm_groups <- function(tpm_data, prot_data, sp) {
    tpm_sp  <- tpm_data  %>% filter(Species == sp)
    prot_sp <- prot_data %>% filter(Species == sp)

    common_ids <- extract_common_ids(tpm_sp, prot_sp)

    tpm_unique <- tpm_sp %>%
        filter(is.finite(mean_log2_TPM)) %>%
        group_by(Protein_id) %>%
        summarise(mean_log2_TPM = mean(mean_log2_TPM, na.rm = TRUE), .groups = "drop")

    tpm_unique %>%
        mutate(
            group = factor(
                ifelse(Protein_id %in% common_ids, "Protein", "Not protein"),
                levels = c("mRNA", "Not protein", "Protein")
            )
        ) %>%
        bind_rows(
            tpm_unique %>%
                mutate(group = factor("mRNA", levels = c("mRNA", "Not protein", "Protein")))
        )
}

#' Plot TPM Distribution by Detection Group for One Species
#'
#' @description
#' Overlapping histogram with three layers (mRNA, Not protein, Protein)
#' reproducing the reference figure. Layers are drawn back-to-front so
#' smaller distributions remain visible.
#'
#' @param classified_df  Tibble. Output of classify_tpm_groups() for one species.
#' @param species        Chr. Species name for italic panel title.
#' @param x_max          Num. Maximum x-axis value. If NULL inferred from data.
#'
#' @return A ggplot object.
plot_tpm_detection_species <- function(classified_df, species, x_max = NULL) {
    if (is.null(x_max)) {
        x_max <- ceiling(max(classified_df$mean_log2_TPM, na.rm = TRUE))
    }

    # Species title: handle "Salmonella_enterica" -> "S. enterica Typhimurium"
    # style formatting — replace underscore with space, italicise genus + epithet
    title_str <- gsub("_", " ", species)

    # Draw order: mRNA (back), Not protein (middle), Protein (front)
    layer_order <- c("mRNA", "Not protein", "Protein")

    group_colours <- c(
        "mRNA"        = "#C2E3F2", # light blue
        "Not protein" = "#A3BECC", # dark blue-grey
        "Protein"     = "#FFD892" # orange/wheat
    )

    ggplot(
        classified_df %>%
            mutate(group = factor(group, levels = layer_order)),
        aes(x = mean_log2_TPM, fill = group)
    ) +
        geom_histogram(
            data     = ~ filter(.x, group == "mRNA"),
            binwidth = 0.5,
            alpha    = 0.85,
            colour   = NA
        ) +
        geom_histogram(
            data     = ~ filter(.x, group == "Not protein"),
            binwidth = 0.5,
            alpha    = 0.85,
            colour   = NA
        ) +
        geom_histogram(
            data     = ~ filter(.x, group == "Protein"),
            binwidth = 0.5,
            alpha    = 0.85,
            colour   = NA
        ) +
        scale_fill_manual(
            values = group_colours,
            breaks = layer_order, # legend order: mRNA, Protein, Not protein
            name   = NULL
        ) +
        scale_x_continuous(
            limits = c(0, x_max),
            expand = expansion(mult = c(0, 0.02)),
            breaks = seq(0, x_max, by = 2.5)
        ) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
        labs(
            title = bquote(italic(.(title_str))),
            x     = "Mean TPM (log2)",
            y     = "Number of genes"
        ) +
        theme_bw(base_size = 11) +
        theme(
            plot.title = element_text(face = "italic", hjust = 0.5, size = 11),
            axis.title = element_text(size = 10),
            axis.text = element_text(size = 9),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(colour = "grey92"),
            legend.position = c(0.72, 0.82),
            legend.background = element_rect(fill = "white", colour = NA),
            legend.key.size = unit(0.4, "cm"),
            legend.text = element_text(size = 9)
        )
}
