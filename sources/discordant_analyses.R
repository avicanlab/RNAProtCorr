# Discordant mRNA-Protein regulation pattern ---------------------------------
#
# Description:
#
# Inputs:
#
# Outputs:
#
# Requirements:
#

library(ggstar)

CLUSTER_COLOURS <- c(
    "Cluster 1" = "#E8A838",
    "Cluster 2" = "#E87070",
    "Cluster 3" = "#C47EC4",
    "Cluster 4" = "#4BBFBF",
    "Cluster 5" = "#4BAF6E",
    "Cluster 6" = "#E8A070"
)

CLUSTER_LABELS <- c(
    "Cluster 1" = "down \u2013 no change",
    "Cluster 2" = "down \u2013 up",
    "Cluster 3" = "no change \u2013 down",
    "Cluster 4" = "no change \u2013 up",
    "Cluster 5" = "up \u2013 down",
    "Cluster 6" = "up \u2013 no change"
)

#' Build LogFC Join Table and Assign Discordance Clusters
#'
#' @description
#' Joins RNA (DEG) and protein (DEP) log2FC values per protein and treatment,
#' then assigns each protein to a discordance cluster based on the direction
#' of regulation across all stress conditions:
#'
#'   Cluster 1: mRNA down in >= 1 condition, protein unchanged
#'   Cluster 2: mRNA down in >= 1 condition, protein up in >= 1 condition
#'   Cluster 3: mRNA unchanged in all conditions, protein down in >= 1 condition
#'   Cluster 4: mRNA unchanged in all conditions, protein up in >= 1 condition
#'   Cluster 5: mRNA up in >= 1 condition, protein down in >= 1 condition
#'   Cluster 6: mRNA up in >= 1 condition, protein unchanged
#'
#' Proteins not matching any cluster are filtered out.
#'
#' @param deg_raw    Tibble. Raw DEG data with columns: Species, Protein_id,
#'   Treatment, log2FC, fdr.
#' @param dep_raw    Tibble. Raw DEP data with columns: Species, Protein_id,
#'   Treatment, logFC, adj.P.Val.
#' @param fc_threshold   Num. |log2FC| threshold for calling regulation.
#' @param pval_threshold Num. Adjusted p-value threshold.
#'
#' @return Tibble with columns: Species, Protein_id, Treatment,
#'   rna_log2FC, prot_log2FC, rna_sig, prot_sig, cluster
build_discordance_clusters <- function(
  deg_raw,
  dep_raw,
  fc_threshold = DE_FC_THRESHOLD,
  pval_threshold = DEP_PVAL_THRESHOLD
) {
    # Tag RNA regulation direction per protein x treatment
    rna_df <- deg_raw %>%
        filter(!is.na(log2FC), !is.na(fdr)) %>%
        mutate(
            rna_sig = case_when(
                fdr < pval_threshold & log2FC < -fc_threshold ~ "down",
                fdr < pval_threshold & log2FC > fc_threshold ~ "up",
                TRUE ~ "unchanged"
            )
        ) %>%
        dplyr::select(Species, Protein_id, Treatment, rna_log2FC = log2FC, rna_sig)

    # Tag protein regulation direction per protein x treatment
    prot_df <- dep_raw %>%
        mutate(
            prot_sig = case_when(
                adj.P.Val < pval_threshold & logFC < -fc_threshold ~ "down",
                adj.P.Val < pval_threshold & logFC > fc_threshold ~ "up",
                TRUE ~ "unchanged"
            )
        ) %>%
        dplyr::select(Species, Protein_id, Treatment, prot_log2FC = logFC, prot_sig)

    joined <- inner_join(rna_df, prot_df, by = c("Species", "Protein_id", "Treatment"))

    protein_pattern <- joined %>%
        group_by(Species, Protein_id, Treatment) %>%
        summarise(
            rna_ever_down = any(rna_sig == "down"),
            rna_ever_up = any(rna_sig == "up"),
            rna_ever_changed = any(rna_sig != "unchanged"),
            prot_ever_down = any(prot_sig == "down"),
            prot_ever_up = any(prot_sig == "up"),
            prot_ever_changed = any(prot_sig != "unchanged"),
            .groups = "drop"
        ) %>%
        mutate(
            cluster = case_when(
                rna_ever_down & !rna_ever_up & !prot_ever_changed ~ "Cluster 1",
                rna_ever_down & !rna_ever_up & prot_ever_up ~ "Cluster 2",
                !rna_ever_changed & prot_ever_down ~ "Cluster 3",
                !rna_ever_changed & prot_ever_up ~ "Cluster 4",
                rna_ever_up & !rna_ever_down & prot_ever_down ~ "Cluster 5",
                rna_ever_up & !rna_ever_down & !prot_ever_changed ~ "Cluster 6",
                TRUE ~ NA_character_
            )
        )

    protein_pattern <- protein_pattern %>%
        filter(!is.na(cluster)) %>%
        dplyr::select(Species, Protein_id, cluster)

    joined %>%
        inner_join(protein_pattern, by = c("Species", "Protein_id")) %>%
        dplyr::select(
            Species, Protein_id, Treatment,
            rna_log2FC, prot_log2FC,
            rna_sig, prot_sig, cluster
        ) %>%
        arrange(Species, cluster, Protein_id, Treatment)
}

#' Extract Raw Per-Protein DEP Results from limma
#'
#' @return Tibble: Species, Protein_id, Treatment, logFC, adj.P.Val
extract_dep_raw <- function(mat, metadata, species,
                            pval_threshold = DEP_PVAL_THRESHOLD,
                            fc_threshold = DE_FC_THRESHOLD) {
    stress_present <- intersect(STRESS_TREATMENTS_PROT, unique(metadata$Treatment))

    metadata <- metadata %>%
        mutate(
            Treatment = factor(Treatment, levels = c("Ctrl", stress_present)),
            Replicate = factor(Replicate)
        )

    design <- model.matrix(~ 0 + Treatment + Replicate, data = metadata)
    colnames(design) <- str_replace(colnames(design), "^Treatment", "")
    contrast_mat <- makeContrasts(
        contrasts = paste0(stress_present, " - Ctrl"),
        levels    = design
    )
    colnames(contrast_mat) <- stress_present
    fit2 <- eBayes(contrasts.fit(lmFit(mat, design), contrast_mat))

    map_dfr(stress_present, function(trt) {
        topTable(fit2,
            coef = trt, number = Inf,
            sort.by = "none", adjust.method = "BH"
        ) %>%
            rownames_to_column("Protein_id") %>%
            mutate(Species = species, Treatment = trt) %>%
            dplyr::select(Species, Protein_id, Treatment, logFC, adj.P.Val)
    })
}

#' Plot mRNA vs Protein LogFC Scatter — Discordance Clusters
#'
#' @description
#' Scatter plot with rna_log2FC on x-axis and prot_log2FC on y-axis.
#' Points are coloured by discordance cluster, grey for non-discordant.
#' One facet per species. Dashed reference lines at x=0 and y=0.
#'
#' @param discordance_df Tibble. Output of build_discordance_clusters().
#'   Must contain: Species, Protein_id, Treatment, rna_log2FC, prot_log2FC, cluster.
#' @param dep_raw        Tibble. Full DEP raw data to include non-discordant
#'   proteins as grey background points.
#' @param deg_raw        Tibble. Full DEG raw data.
#' @param output_path    Chr. Directory for saving output files.
#'
#' @return ggplot object invisibly.
plot_discordance_scatter <- function(
  discordance_df, dep_raw, deg_raw, output_path,
  xlim = NULL, ylim = NULL, save_per_species = FALSE, prefix = NULL
) {
    all_logfc <- inner_join(
        deg_raw %>% dplyr::select(Species, Protein_id, Treatment, rna_log2FC = log2FC),
        dep_raw %>% dplyr::select(Species, Protein_id, Treatment, prot_log2FC = logFC),
        by = c("Species", "Protein_id", "Treatment")
    )

    background_df <- all_logfc %>%
        anti_join(
            discordance_df %>% distinct(Species, Protein_id),
            by = c("Species", "Protein_id")
        ) %>%
        filter(is.finite(rna_log2FC), is.finite(prot_log2FC))

    foreground_df <- discordance_df %>%
        filter(is.finite(rna_log2FC), is.finite(prot_log2FC)) %>%
        mutate(cluster = factor(cluster, levels = names(CLUSTER_COLOURS)))

    # Order clusters by size descending so smaller clusters are drawn on top
    cluster_order <- foreground_df %>%
        count(cluster) %>%
        arrange(desc(n)) %>%
        pull(cluster)

    foreground_df <- foreground_df %>%
        mutate(cluster = factor(cluster, levels = cluster_order))

    species_list <- sort(unique(foreground_df$Species))

    plots <- map(species_list, function(species) {
        bg <- background_df %>% filter(Species == species)
        fg <- foreground_df %>% filter(Species == species)

        n_discordant <- n_distinct(fg$Protein_id)
        count_label <- paste0("n = ", n_discordant, " discordant")

        x_lim <- if (!is.null(xlim)) {
            xlim
        } else {
            rna_range <- range(c(bg$rna_log2FC, fg$rna_log2FC), na.rm = TRUE)
            c(floor(rna_range[1]), ceiling(rna_range[2]))
        }
        y_lim <- if (!is.null(ylim)) {
            ylim
        } else {
            prot_range <- range(c(bg$prot_log2FC, fg$prot_log2FC), na.rm = TRUE)
            c(floor(prot_range[1]), ceiling(prot_range[2]))
        }

        # Draw clusters from largest to smallest so smaller ones appear on top
        p <- ggplot() +
            geom_blank(
                data = tibble(rna_log2FC = x_lim, prot_log2FC = y_lim),
                aes(x = rna_log2FC, y = prot_log2FC)
            ) +
            geom_point(
                data = bg,
                aes(x = rna_log2FC, y = prot_log2FC),
                colour = "grey75", alpha = 0.4, size = 0.8
            )

        # Add clusters one by one from largest to smallest
        for (cl in levels(fg$cluster)) {
            p <- p +
                geom_star(
                    data = fg %>% filter(cluster == cl),
                    aes(
                        x = rna_log2FC, y = prot_log2FC,
                        fill = cluster, starshape = Treatment
                    ),
                    alpha = 0.7, size = 1.5, colour = NA
                )
        }

        p <- p +
            geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30", linewidth = 0.5) +
            geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30", linewidth = 0.5) +
            annotate(
                "text",
                x = x_lim[1], y = y_lim[2],
                label = count_label,
                hjust = 0, vjust = 1,
                size = 3, colour = "grey30", fontface = "italic"
            ) +
            scale_fill_manual(
                values = CLUSTER_COLOURS,
                labels = CLUSTER_LABELS,
                name   = "mRNA \u2013 Protein",
                breaks = names(CLUSTER_COLOURS),
                drop   = FALSE
            ) +
            scale_starshape_manual(
                values = TREATMENT_SHAPES,
                name   = "Treatment"
            ) +
            guides(
                starshape = guide_legend(
                    override.aes = list(
                        fill   = "grey40",
                        colour = "grey40",
                        size   = 3,
                        alpha  = 1
                    )
                ),
                fill = guide_legend(
                    override.aes = list(
                        starshape = 15, # solid square for cluster legend
                        colour    = "white",
                        size      = 3,
                        alpha     = 1
                    )
                )
            ) +
            scale_x_continuous(
                limits = x_lim,
                breaks = scales::pretty_breaks(n = 5),
                expand = c(0, 0)
            ) +
            scale_y_continuous(
                limits = y_lim,
                breaks = scales::pretty_breaks(n = 5),
                expand = c(0, 0)
            ) +
            labs(
                title = format_species_title(species),
                x     = "logFC(mRNA)",
                y     = "logFC(Protein)"
            ) +
            theme_bw(base_size = 11) +
            theme(
                plot.title       = element_text(hjust = 0.5, size = 10),
                panel.grid.minor = element_blank(),
                panel.grid.major = element_line(colour = "grey92"),
                legend.position  = "right",
                legend.title     = element_text(size = 9),
                legend.text      = element_text(size = 8),
                legend.key.size  = unit(0.4, "cm"),
                axis.title       = element_text(size = 10),
                axis.text        = element_text(size = 9)
            )

        if (save_per_species) {
            output_file <- file.path(
                output_path,
                species,
                paste0(
                    if (!is.null(prefix)) paste0(prefix, "_") else "",
                    "discordance_scatter"
                )
            )
            save_plot(
                plot = p,
                filepath = output_file,
                width = 6, height = 5
            )
            message("Discordance scatter saved: ", output_file, MSG_FIG_FORMAT)
        }

        p
    })

    plots
}

plot_discordance_profiles <- function(discordance_df, output_path,
                                      species_filter = NULL) {
    treatment_order <- sort(unique(discordance_df$Treatment[discordance_df$Treatment != "Ctrl"]))

    df <- discordance_df %>%
        filter(Treatment != "Ctrl") %>%
        {
            if (!is.null(species_filter)) filter(., Species == species_filter) else .
        } %>%
        mutate(
            Treatment = factor(Treatment, levels = treatment_order),
            cluster   = factor(cluster, levels = names(CLUSTER_COLOURS))
        ) %>%
        pivot_longer(
            cols      = c(rna_log2FC, prot_log2FC),
            names_to  = "measure",
            values_to = "log2FC"
        ) %>%
        mutate(
            measure = factor(
                ifelse(measure == "rna_log2FC", "mRNA", "Protein"),
                levels = c("mRNA", "Protein")
            )
        )

    count_df <- discordance_df %>%
        filter(Treatment != "Ctrl") %>%
        {
            if (!is.null(species_filter)) filter(., Species == species_filter) else .
        } %>%
        distinct(Species, Protein_id, cluster) %>%
        count(Species, cluster, name = "n") %>%
        mutate(cluster = factor(cluster, levels = names(CLUSTER_COLOURS)))

    species_levels <- sort(unique(df$Species))
    species_labeller <- as_labeller(setNames(
        sapply(species_levels, function(s) {
            display <- PLOT_SPECIES_NAMES[s]
            ifelse(is.na(display), gsub("_", " ", s), display)
        }),
        species_levels
    ))

    facet <- if (is.null(species_filter)) {
        facet_grid(
            Species ~ cluster,
            scales = "free_y",
            labeller = labeller(
                Species = species_labeller,
                cluster = as_labeller(CLUSTER_LABELS)
            )
        )
    } else {
        facet_wrap(
            ~cluster,
            nrow     = 1,
            scales   = "free_y",
            labeller = as_labeller(CLUSTER_LABELS)
        )
    }

    p <- ggplot(
        df,
        aes(
            x = Treatment, y = log2FC,
            group = interaction(Protein_id, measure),
            colour = measure
        )
    ) +
        geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
        geom_line(alpha = 0.3, linewidth = 0.3) +
        geom_text(
            data = count_df,
            aes(x = length(treatment_order), y = Inf, label = paste0("n=", n)),
            hjust = 1.1, vjust = 1.5,
            size = 2.8, colour = "grey30", fontface = "italic",
            inherit.aes = FALSE
        ) +
        facet +
        scale_colour_manual(
            values = c("mRNA" = "#E8C84A", "Protein" = "#7BBFDF"),
            name   = NULL
        ) +
        scale_x_discrete(guide = guide_axis(angle = 45)) +
        labs(x = NULL, y = "Log2 FC") +
        theme_bw(base_size = 10) +
        theme(
            strip.text.x       = element_text(size = 9, face = "bold", colour = "white"),
            strip.text.y       = element_text(size = 9, face = "italic", angle = -90),
            strip.background.x = element_rect(fill = "grey40"),
            strip.background.y = element_rect(fill = "white"),
            panel.grid.minor   = element_blank(),
            panel.grid.major.x = element_blank(),
            legend.position    = "right",
            legend.text        = element_text(size = 9),
            legend.key.size    = unit(0.4, "cm"),
            axis.text.x        = element_text(size = 7),
            axis.text.y        = element_text(size = 7),
            axis.title.y       = element_text(size = 9)
        )

    # Colour cluster strip backgrounds
    g <- ggplot_gtable(ggplot_build(p))
    strip_idx <- which(grepl("strip-t", g$layout$name))
    cluster_names <- names(CLUSTER_COLOURS)
    for (i in seq_along(strip_idx)) {
        col <- CLUSTER_COLOURS[cluster_names[((i - 1) %% length(cluster_names)) + 1]]
        g$grobs[[strip_idx[i]]]$grobs[[1]]$children[[1]]$gp$fill <- col
    }

    wrap_elements(g)
}
