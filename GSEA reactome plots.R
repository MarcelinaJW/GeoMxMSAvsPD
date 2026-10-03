################################################################################
### GSEA Analysis — GO, KEGG, SynGO & Reactome - better plots, more refined
################################################################################

suppressPackageStartupMessages({
  library(fgsea)
  library(ReactomePA)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(GOSemSim)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(patchwork)
  library(scales)
})

cat("\n===GSEA Analysis ===\n")
cat("    GO | KEGG | SynGO | Reactome\n")
cat(rep("=", 50), "\n", sep = "")

################################################################################
###Shared theme
################################################################################

theme_gsea <- theme_classic(base_size = 11, base_family = "sans") +
  theme(
    # Panel
    panel.border      = element_rect(colour = "black", fill = NA,
                                     linewidth = 0.7),
    panel.background  = element_rect(fill = "white"),
    panel.grid.major.x = element_line(colour = "grey95", linewidth = 0.3),
    panel.grid.minor  = element_blank(),
    
    # Strip
    strip.background  = element_rect(fill = "grey95", colour = "black",
                                     linewidth = 0.5),
    strip.text        = element_text(face = "bold", size = 10,
                                     margin = margin(3, 0, 3, 0)),
    
    # Axes
    axis.line         = element_line(colour = "black", linewidth = 0.5),
    axis.ticks        = element_line(colour = "black", linewidth = 0.5),
    axis.text.y       = element_text(size = 9,  colour = "black"),
    axis.text.x       = element_text(size = 9,  colour = "black"),
    axis.title        = element_text(face = "bold", size = 10),
    axis.title.x      = element_text(margin = margin(t = 10)),
    
    # Titles
    plot.title        = element_text(face = "bold", size = 12, hjust = 0),
    plot.subtitle     = element_text(size = 10, colour = "grey30", hjust = 0),
    plot.caption      = element_text(size = 8,  colour = "grey40",
                                     hjust = 0, margin = margin(t = 8)),
    plot.tag          = element_text(size = 14, face = "bold"),
    
    # Legend
    legend.position   = "right",
    legend.background = element_rect(fill = "white", colour = NA),
    legend.key        = element_rect(fill = "white", colour = NA),
    legend.title      = element_text(size = 9, face = "bold"),
    legend.text       = element_text(size = 8),
    legend.spacing.y  = unit(0.2, "cm"),
    
    plot.margin       = margin(10, 20, 10, 20)
  )

################################################################################
### Reactome GSEA runner
################################################################################

run_reactome_gsea <- function(de_results,
                              comparison_name,
                              minSize     = 10,
                              maxSize     = 500,
                              out_dir     = OUTPUT_DIR) {
  
  cat("\n--- Reactome GSEA:", comparison_name, "---\n")
  
  # ── Ranked gene list (t-statistic) ──────────────────────────────────────────
  gene_list <- de_results$t
  symbols   <- rownames(de_results)
  
  ids <- mapIds(
    org.Hs.eg.db,
    keys      = symbols,
    column    = "ENTREZID",
    keytype   = "SYMBOL",
    multiVals = "first"
  )
  
  valid_idx <- !is.na(ids) & !duplicated(ids)
  gene_list <- gene_list[valid_idx]
  names(gene_list) <- ids[valid_idx]
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  cat("  Genes in ranked list:", length(gene_list), "\n")
  
  # ── Reactome pathways ────────────────────────────────────────────────────────
  pathways_raw      <- as.list(reactomePATHID2EXTID)
  pathways_filtered <- pathways_raw[vapply(pathways_raw, function(x) {
    len <- length(intersect(x, names(gene_list)))
    len >= minSize && len <= maxSize
  }, logical(1))]
  
  cat("  Pathways tested:", length(pathways_filtered), "\n")
  
  # ── fgsea ────────────────────────────────────────────────────────────────────
  fgsea_res <- fgseaMultilevel(
    pathways = pathways_filtered,
    stats    = gene_list,
    eps      = 0,
    minSize  = minSize,
    maxSize  = maxSize
  )
  
  # ── Add descriptions & standardise column names ───────────────────────────
  path_names <- as.list(reactomePATHID2NAME)
  
  fgsea_res <- fgsea_res %>%
    mutate(
      leadingEdge = vapply(leadingEdge, paste, character(1), collapse = ";"),
      Description = unlist(path_names[pathway]),
      Description = Description %>%
        str_remove("^Homo sapiens[: -]+") %>%
        str_remove(" \\(Homo sapiens\\)$") %>%
        str_trim(),
      # Standardise to enrichResult column names
      ID      = pathway,
      qvalue  = padj,
      setSize = size
    ) %>%
    dplyr::select(ID, Description, pval, padj, qvalue,
                  ES, NES, setSize, size, leadingEdge) %>%
    arrange(padj)
  
  cat("  Significant (padj < 0.05):",
      sum(fgsea_res$padj < 0.05, na.rm = TRUE), "\n")
  
  # ── Save CSV ─────────────────────────────────────────────────────────────────
  out_path <- file.path(
    out_dir,
    paste0("GSEA_Reactome_", gsub("[^A-Za-z0-9]", "_", comparison_name), ".csv")
  )
  write.csv(fgsea_res, out_path, row.names = FALSE)
  cat("  Saved:", basename(out_path), "\n")
  
  # ── Wrap as enrichResult for compatibility with plot_gsea_publication() ──────
  reactome_er <- new(
    "enrichResult",
    result        = as.data.frame(fgsea_res),
    pvalueCutoff  = 1,
    pAdjustMethod = "BH",
    qvalueCutoff  = 1,
    organism      = "human",
    ontology      = "Reactome",
    gene          = character(0),
    keytype       = "ENTREZID"
  )
  
  return(list(
    enrichResult = reactome_er,
    results      = fgsea_res,
    gene_list    = gene_list
  ))
}

################################################################################
### Core dotplot function
################################################################################

plot_gsea_publication <- function(gsea_res,
                                  comparison_name,
                                  ontology_name,
                                  padj_cutoff       = 0.25,
                                  top_n             = 10,
                                  remove_redundancy = FALSE) {
  
  # ── Input checks ─────────────────────────────────────────────────────────────
  if (is.null(gsea_res)) {
    cat("  [SKIP]", ontology_name, ": NULL\n")
    return(NULL)
  }
  
  if (nrow(gsea_res@result) == 0) {
    cat("  [SKIP]", ontology_name, ": No results\n")
    return(NULL)
  }
 
  # ── Filter & annotate ─────────────────────────────────────────────────────────
  df <- gsea_res@result %>%
    filter(!is.na(qvalue), qvalue < padj_cutoff) %>%
    mutate(
      Direction   = ifelse(NES > 0, "Up in MSA", "Up in PD"),
      Direction   = factor(Direction, levels = c("Up in MSA", "Up in PD")),
      abs_NES     = abs(NES),
      neg_log10_q = -log10(qvalue)
    )
  
  if (nrow(df) == 0) {
    cat("  [SKIP]", ontology_name,
        ": No significant terms at q <", padj_cutoff, "\n")
    return(NULL)
  }
  
  cat("  [OK]", ontology_name, ":", nrow(df), "significant terms\n")
  
  # ── Select top N per direction ────────────────────────────────────────────────
  df_top <- df %>%
    group_by(Direction) %>%
    arrange(desc(abs_NES)) %>%
    slice_head(n = top_n) %>%
    ungroup()
  
  # ── Clean descriptions ────────────────────────────────────────────────────────
  df_top <- df_top %>%
    mutate(
      Description = Description %>%
        gsub("\\s*\\(GO:\\d+\\)", "", .) %>%
        gsub("\\s*GO:\\d+",       "", .) %>%
        trimws() %>%
        str_wrap(width = 50)
    )
  
  # ── Order by NES within direction ─────────────────────────────────────────────
  df_top <- df_top %>%
    arrange(Direction, NES) %>%
    mutate(Description = factor(Description, levels = unique(Description)))
  
  # ── Subtitle counts ───────────────────────────────────────────────────────────
  n_msa <- sum(df$Direction == "Up in MSA")
  n_pd  <- sum(df$Direction == "Up in PD")
  
  # ── Plot ──────────────────────────────────────────────────────────────────────
  p <- ggplot(df_top, aes(x = NES, y = Description)) +
    
    geom_vline(xintercept = 0,
               linetype   = "solid",
               colour     = "grey30",
               linewidth  = 0.5) +
    
    geom_point(aes(size = setSize, fill = neg_log10_q),
               shape  = 21,
               colour = "black",
               alpha  = 0.85,
               stroke = 0.35) +
    
    facet_grid(Direction ~ ., scales = "free_y", space = "free_y") +
    
    scale_fill_viridis_c(
      option    = "C",
      direction = -1,
      name      = expression(-log[10]*"(q-value)"),
      guide     = guide_colorbar(
        barwidth        = 1,
        barheight       = 8,
        title.position  = "top",
        title.hjust     = 0.5
      )
    ) +
    
    scale_size_continuous(
      range = c(2, 8),
      name  = "Gene Set\nSize",
      guide = guide_legend(
        title.position = "top",
        title.hjust    = 0.5,
        override.aes   = list(fill   = "grey70",
                              colour = "black",
                              stroke = 0.3)
      )
    ) +
    
    scale_x_continuous(
      breaks = pretty_breaks(n = 6),
      expand = expansion(mult = c(0.1, 0.1))
    ) +
    
    coord_cartesian(clip = "off") +
    
    labs(
      title    = paste0(ontology_name, " Enrichment"),
      subtitle = paste0(
        comparison_name, "  |  ",
        n_msa + n_pd, " significant pathways  |  ",
        "MSA-up: ", n_msa, "  PD-up: ", n_pd
      ),
      x = "Normalised Enrichment Score (NES)",
      y = NULL
    ) +
    
    theme_gsea
  
  return(p)
}

################################################################################
### Combined SynGO plot
################################################################################

plot_syngo_combined <- function(gsea_bp, gsea_cc,
                                comparison_name,
                                padj_cutoff = 0.25,
                                top_n       = 8) {
  
  # ── Process BP ────────────────────────────────────────────────────────────────
  process_syngo <- function(gsea_res, ontology_label) {
    
    if (is.null(gsea_res)) return(NULL)
    if (nrow(gsea_res@result) == 0) return(NULL)
    
    df <- gsea_res@result %>%
      filter(!is.na(qvalue), qvalue < padj_cutoff) %>%
      mutate(
        Direction   = ifelse(NES > 0, "Up in MSA", "Up in PD"),
        Direction   = factor(Direction, levels = c("Up in MSA", "Up in PD")),
        abs_NES     = abs(NES),
        neg_log10_q = -log10(qvalue),
        Ontology    = ontology_label
      )
    
    if (nrow(df) == 0) return(NULL)
    
    # Top N per direction
    df %>%
      group_by(Direction) %>%
      arrange(desc(abs_NES)) %>%
      slice_head(n = top_n) %>%
      ungroup() %>%
      mutate(
        Description = Description %>%
          gsub("\\s*\\(GO:\\d+\\)", "", .) %>%
          gsub("\\s*GO:\\d+",       "", .) %>%
          trimws() %>%
          str_wrap(width = 40)
      )
  }
  
  df_bp <- process_syngo(gsea_bp, "Biological Process")
  df_cc <- process_syngo(gsea_cc, "Cellular Component")
  
  # ── Combine ───────────────────────────────────────────────────────────────────
  df_combined <- bind_rows(df_bp, df_cc)
  
  if (is.null(df_combined) || nrow(df_combined) == 0) {
    cat("  [SKIP] SynGO combined: no significant terms\n")
    return(NULL)
  }
  
  df_combined <- df_combined %>%
    arrange(Ontology, Direction, NES) %>%
    mutate(
      Description = factor(Description, levels = unique(Description)),
      Ontology    = factor(Ontology,
                           levels = c("Biological Process",
                                      "Cellular Component"))
    )
  
  n_total <- nrow(df_combined)
  n_msa   <- sum(df_combined$Direction == "Up in MSA")
  n_pd    <- sum(df_combined$Direction == "Up in PD")
  
  cat("  [OK] SynGO combined:", n_total, "terms",
      "(BP:", nrow(df_bp), "| CC:", nrow(df_cc), ")\n")
  
  # ── Plot ──────────────────────────────────────────────────────────────────────
  p <- ggplot(df_combined, aes(x = NES, y = Description)) +
    
    geom_vline(xintercept = 0,
               linetype   = "solid",
               colour     = "grey30",
               linewidth  = 0.5) +
    
    geom_point(aes(size = setSize, fill = neg_log10_q),
               shape  = 21,
               colour = "black",
               alpha  = 0.85,
               stroke = 0.35) +
    
    facet_grid(Direction ~ Ontology,
               scales = "free",
               space  = "free_y") +
    
    scale_fill_viridis_c(
      option    = "C",
      direction = -1,
      name      = expression(-log[10]*"(q-value)"),
      guide     = guide_colorbar(
        barwidth       = 1,
        barheight      = 8,
        title.position = "top",
        title.hjust    = 0.5
      )
    ) +
    
    scale_size_continuous(
      range = c(2, 8),
      name  = "Gene Set\nSize",
      guide = guide_legend(
        title.position = "top",
        title.hjust    = 0.5,
        override.aes   = list(fill   = "grey70",
                              colour = "black",
                              stroke = 0.3)
      )
    ) +
    
    scale_x_continuous(
      breaks = pretty_breaks(n = 5),
      expand = expansion(mult = c(0.1, 0.1))
    ) +
    
    coord_cartesian(clip = "off") +
    
    labs(
      title    = "SynGO Enrichment",
      subtitle = paste0(
        comparison_name, "  |  ",
        n_msa + n_pd, " significant terms (q < ", padj_cutoff, ")  |  ",
        "MSA-up: ", n_msa, "  PD-up: ", n_pd
      ),
      x = "Normalised Enrichment Score (NES)",
      y = NULL
    ) +
    
    theme_gsea +
    # Slightly wider strip for two-column facet
    theme(
      strip.text.x = element_text(face = "bold", size = 10),
      strip.text.y = element_text(face = "bold", size = 10)
    )
  
  return(p)
}

################################################################################
### Multi-panel creator — GO + KEGG + SynGO + Reactome
################################################################################
create_gsea_multipanel <- function(gsea_obj,
                                   comparison_name,
                                   padj_cutoff = 0.25,
                                   ontologies  = c("GO", "KEGG",
                                                   "SYNGO", "Reactome")) {
  
  cat("\n---", comparison_name, "---\n")
  
  plots <- list()
  
  # ── GO ────────────────────────────────────────────────────────────────────────
  if ("GO" %in% ontologies && !is.null(gsea_obj$GO)) {
    p <- plot_gsea_publication(
      gsea_obj$GO, comparison_name,
      ontology_name     = "GO Biological Process",
      padj_cutoff       = padj_cutoff,
      top_n             = 9,
      remove_redundancy = TRUE
    )
    if (!is.null(p)) plots[["GO"]] <- p
  }
  
  # ── KEGG ──────────────────────────────────────────────────────────────────────
  if ("KEGG" %in% ontologies && !is.null(gsea_obj$KEGG)) {
    p <- plot_gsea_publication(
      gsea_obj$KEGG, comparison_name,
      ontology_name     = "KEGG Pathways",
      padj_cutoff       = padj_cutoff,
      top_n             = 10,
      remove_redundancy = FALSE
    )
    if (!is.null(p)) plots[["KEGG"]] <- p
  }
  
  # ── SynGO ─────────────────────────────────
  if ("SYNGO" %in% ontologies &&
      (!is.null(gsea_obj$SYNGO_BP) || !is.null(gsea_obj$SYNGO_CC))) {
    
    p <- plot_syngo_combined(
      gsea_bp         = gsea_obj$SYNGO_BP,
      gsea_cc         = gsea_obj$SYNGO_CC,
      comparison_name = comparison_name,
      padj_cutoff     = padj_cutoff,
      top_n           = 8
    )
    if (!is.null(p)) plots[["SYNGO"]] <- p
  }
  
  # ── Reactome ──────────────────────────────────────────────────────────────────
  if ("Reactome" %in% ontologies && !is.null(gsea_obj$Reactome)) {
    p <- plot_gsea_publication(
      gsea_obj$Reactome$enrichResult,
      comparison_name,
      ontology_name     = "Reactome Pathways",
      padj_cutoff       = padj_cutoff,
      top_n             = 10,
      remove_redundancy = FALSE
    )
    if (!is.null(p)) plots[["Reactome"]] <- p
  }
  
  return(plots)
}
################################################################################
### Export tables
################################################################################

export_gsea_tables <- function(gsea_obj, prefix, out_dir = OUTPUT_DIR) {
  
  ont_map <- list(
    GO       = "GO",
    KEGG     = "KEGG",
    SYNGO_BP = "SYNGO_BP",
    SYNGO_CC = "SYNGO_CC"
  )
  
  # Standard enrichResult objects
  for (ont in names(ont_map)) {
    if (!is.null(gsea_obj[[ont]])) {
      out_path <- file.path(out_dir, paste0(prefix, "_", ont_map[[ont]], ".csv"))
      write.csv(gsea_obj[[ont]]@result, out_path, row.names = FALSE)
      cat("  Saved:", basename(out_path), "\n")
    }
  }
  
  # Reactome (fgsea tibble)
  if (!is.null(gsea_obj$Reactome)) {
    out_path <- file.path(out_dir, paste0(prefix, "_Reactome.csv"))
    write.csv(gsea_obj$Reactome$results, out_path, row.names = FALSE)
    cat("  Saved:", basename(out_path), "\n")
  }
}

################################################################################
### Summary statistics
################################################################################

create_gsea_summary <- function(gsea_obj,
                                comparison_name,
                                padj_cutoff = 0.25) {
  
  summary_rows <- list()
  
  # Standard enrichResult objects
  for (ont in c("GO", "KEGG", "SYNGO_BP", "SYNGO_CC")) {
    if (!is.null(gsea_obj[[ont]])) {
      df <- gsea_obj[[ont]]@result %>%
        filter(!is.na(qvalue), qvalue < padj_cutoff)
      
      summary_rows[[ont]] <- data.frame(
        Comparison        = comparison_name,
        Ontology          = ont,
        Total_Significant = nrow(df),
        Up_in_MSA         = sum(df$NES > 0, na.rm = TRUE),
        Up_in_PD          = sum(df$NES < 0, na.rm = TRUE)
      )
    }
  }
  
  # Reactome
  if (!is.null(gsea_obj$Reactome)) {
    df <- gsea_obj$Reactome$results %>%
      filter(!is.na(padj), padj < padj_cutoff)
    
    summary_rows[["Reactome"]] <- data.frame(
      Comparison        = comparison_name,
      Ontology          = "Reactome",
      Total_Significant = nrow(df),
      Up_in_MSA         = sum(df$NES > 0, na.rm = TRUE),
      Up_in_PD          = sum(df$NES < 0, na.rm = TRUE)
    )
  }
  
  bind_rows(summary_rows)
}

################################################################################
###Save
################################################################################
save_gsea_pdf <- function(plots, filename, out_dir = OUTPUT_DIR) {
  
  if (length(plots) == 0) {
    cat("  [SKIP] No plots to save for", filename, "\n")
    return(invisible(NULL))
  }
  
  panel_heights <- vapply(names(plots), function(nm) {
    if (nm == "SYNGO") 6 else 6 
  }, numeric(1))
  
  p_combined <- wrap_plots(plots, ncol = 1,
                           heights = panel_heights) +
    plot_annotation(
      tag_levels = "A",
      tag_suffix = ".",
      theme      = theme(
        plot.tag = element_text(size = 14, face = "bold")
      )
    )
  
  total_height <- sum(panel_heights)
  
  out_path <- file.path(out_dir, filename)
  
  ggsave(out_path, p_combined,
         width = 12, height = total_height,
         limitsize = FALSE)
  
  ggsave(sub("\\.pdf$", ".png", out_path), p_combined,
         width = 12, height = total_height,
         dpi = 600, bg = "white", limitsize = FALSE)
  
  cat("  Saved:", basename(out_path), "\n")
}
################################################################################
### Integrated enrichment summary barplot
################################################################################

create_integrated_summary <- function(summary_df, out_dir = OUTPUT_DIR) {
  
  # Reshape
  summary_long <- summary_df %>%
    pivot_longer(
      cols      = c(Up_in_MSA, Up_in_PD),
      names_to  = "Direction",
      values_to = "Count"
    ) %>%
    mutate(
      Direction = gsub("Up_in_", "", Direction),
      Direction = factor(Direction, levels = c("MSA", "PD")),
      Ontology  = factor(Ontology,
                         levels = c("GO", "KEGG",
                                    "SYNGO_BP", "SYNGO_CC", "Reactome"))
    )
  
  p <- ggplot(summary_long,
              aes(x = Ontology, y = Count, fill = Direction)) +
    
    geom_col(position  = "dodge",
             colour    = "black",
             linewidth = 0.3,
             alpha     = 0.85) +
    
    geom_text(aes(label = Count),
              position = position_dodge(width = 0.9),
              vjust    = -0.4,
              size     = 3,
              fontface = "bold") +
    
    facet_wrap(~ Comparison, ncol = 3) +
    
    scale_fill_manual(
      values = c("MSA" = "#D55E00", "PD" = "#0072B2"),
      name   = "Enriched in"
    ) +
    
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    
    labs(
      title   = "Pathway Enrichment Summary Across All Databases",
      caption = paste0(
        "GO = Gene Ontology Biological Process; ",
        "KEGG = Kyoto Encyclopaedia of Genes and Genomes; ",
        "SynGO = Synaptic Gene Ontologies; ",
        "Reactome = Reactome pathway database."
      ),
      x = NULL,
      y = "Number of Enriched Pathways"
    ) +
    
    theme_gsea +
    theme(
      axis.text.x  = element_text(angle = 40, hjust = 1),
      legend.position = "bottom",
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3)
    )
  
  # Save
  ggsave(file.path(out_dir, "07k_integrated_enrichment_summary.pdf"),
         p, width = 13, height = 6)
  ggsave(file.path(out_dir, "07k_integrated_enrichment_summary.png"),
         p, width = 13, height = 6, dpi = 600, bg = "white")
  
  write.csv(summary_df,
            file.path(out_dir, "Integrated_enrichment_summary.csv"),
            row.names = FALSE)
  
  cat("  Saved: 07k_integrated_enrichment_summary.pdf / .png\n")
  cat("  Saved: Integrated_enrichment_summary.csv\n")
  
  return(p)
}

################################################################################
### Reactome GSEA for all comparisons
################################################################################

cat("\n", rep("=", 50), "\n", sep = "")
cat("Running Reactome GSEA\n")
cat(rep("=", 50), "\n", sep = "")

reactome_neuron <- run_reactome_gsea(de_neuron, "MSA vs PD (Neurons)")
reactome_oligo  <- run_reactome_gsea(de_oligo,  "MSA vs PD (Oligodendrocytes)")
reactome_all    <- run_reactome_gsea(de_all,    "MSA vs PD (All ROIs)")

# ── Attach Reactome to existing GSEA list objects ─────────────────────────────
gsea_msa_vs_pd_neuron$Reactome <- reactome_neuron
gsea_msa_vs_pd_oligo$Reactome  <- reactome_oligo
gsea_msa_vs_pd_all$Reactome    <- reactome_all

################################################################################
### Multi-panel plots for all comparisons
################################################################################

cat("\n", rep("=", 50), "\n", sep = "")
cat("Creating multi-panel GSEA figures\n")
cat(rep("=", 50), "\n", sep = "")

PADJ <- 0.25   # single threshold — change here to update everywhere

plots_neuron <- create_gsea_multipanel(
  gsea_msa_vs_pd_neuron,
  "MSA vs PD (Neurons)",
  padj_cutoff = PADJ
)

plots_oligo <- create_gsea_multipanel(
  gsea_msa_vs_pd_oligo,
  "MSA vs PD (Oligodendrocytes)",
  padj_cutoff = PADJ
)

plots_all <- create_gsea_multipanel(
  gsea_msa_vs_pd_all,
  "MSA vs PD (All ROIs)",
  padj_cutoff = PADJ
)

################################################################################
### Save and export
################################################################################

cat("\n", rep("=", 50), "\n", sep = "")
cat("Saving figures\n")
cat(rep("=", 50), "\n", sep = "")

save_gsea_pdf(plots_neuron, "07a_GSEA_neurons_publication.pdf")
save_gsea_pdf(plots_oligo,  "07b_GSEA_oligos_publication.pdf")
save_gsea_pdf(plots_all,    "07c_GSEA_all_publication.pdf")


cat("\n", rep("=", 50), "\n", sep = "")
cat("Exporting results tables\n")
cat(rep("=", 50), "\n", sep = "")

export_gsea_tables(gsea_msa_vs_pd_neuron, "GSEA_neuron")
export_gsea_tables(gsea_msa_vs_pd_oligo,  "GSEA_oligo")
export_gsea_tables(gsea_msa_vs_pd_all,    "GSEA_all")

################################################################################
### Summary statistics & integrated barplot
################################################################################

cat("\n", rep("=", 50), "\n", sep = "")
cat("Creating summary statistics\n")
cat(rep("=", 50), "\n", sep = "")

summary_all <- bind_rows(
  create_gsea_summary(gsea_msa_vs_pd_neuron, "Neurons",          PADJ),
  create_gsea_summary(gsea_msa_vs_pd_oligo,  "Oligodendrocytes", PADJ),
  create_gsea_summary(gsea_msa_vs_pd_all,    "All ROIs",         PADJ)
)

write.csv(summary_all,
          file.path(OUTPUT_DIR, "GSEA_summary_statistics.csv"),
          row.names = FALSE)

cat("\n=== GSEA Summary ===\n")
print(summary_all)

p_integrated <- create_integrated_summary(summary_all)

cat("\n", rep("=", 50), "\n", sep = "")
cat("=== GSEA Analysis Complete ===\n")
cat(rep("=", 50), "\n", sep = "")
cat("\nFigures saved:\n")
cat("  07a  Neurons — all ontologies (stacked)\n")
cat("  07b  Oligodendrocytes — all ontologies (stacked)\n")
cat("  07c  All ROIs — all ontologies (stacked)\n")
cat("  07k  Integrated enrichment summary barplot\n")
cat("\nTables saved:\n")
cat("  GSEA_neuron / oligo / all — _GO / _KEGG / _SYNGO_BP / _SYNGO_CC / _Reactome .csv\n")
cat("  GSEA_summary_statistics.csv\n")
cat("  Integrated_enrichment_summary.csv\n")
cat("  GSEA_Reactome_*.csv (per comparison)\n")
cat("\nAll figures exported as PDF + 600 dpi PNG\n")