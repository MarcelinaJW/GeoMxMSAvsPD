###############################################################################
###ROI-Level Clustering Analysis with no Pseudobulking!!########################
###############################################################################
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(dplyr)
library(patchwork)
library(cluster)
library(RColorBrewer)
library(grid)

cat("\n=== ROI-Level Clustering Analysis ===\n")
###############################################################################
### Color schemes################################
###############################################################################
col_fun_expr <- colorRamp2(
  c(-2, -1, 0, 1, 2),
  c("#2166AC", "#92C5DE", "white", "#F4A582", "#B2182B")
)

###############################################################################
### ROI-level clustering function################################################
###############################################################################

run_roi_clustering <- function(spe_obj, cell_type_label,
                               n_var_genes = 500,
                               output_prefix = "10") {
  
  cat("\n--- ROI-level:", cell_type_label, "---\n")
  
  # Subset to cell type
  spe_sub <- spe_obj[, spe_obj$Cell == cell_type_label]
  cat("  ROIs:", ncol(spe_sub), "\n")
  cat("  Patients:", length(unique(spe_sub$Patient)), "\n")
  
  # ROI-level expression matrix (genes x ROIs) — NO aggregation
  roi_mat <- assay(spe_sub, "logcounts")
  
  # ROI-level metadata
  roi_groups  <- as.character(spe_sub$Group)
  roi_patient <- as.character(spe_sub$Patient)
  roi_ids     <- colnames(spe_sub)
  names(roi_groups)  <- roi_ids
  names(roi_patient) <- roi_ids
  
  # Select top variable genes (across ROIs, not patients)
  gene_vars <- apply(roi_mat, 1, var, na.rm = TRUE)
  n_var_use <- min(n_var_genes, sum(gene_vars > 0))
  top_var_genes <- names(sort(gene_vars, decreasing = TRUE)[1:n_var_use])
  
  cat("  Variable genes:", n_var_use, "\n")
  
  # Scale (Z-score) across ROIs
  roi_mat_scaled <- t(scale(t(roi_mat)))
  roi_mat_scaled[is.nan(roi_mat_scaled)] <- 0
  
  # Subset to top variable genes
  mat_plot <- roi_mat_scaled[top_var_genes, ]
  
  # Order ROIs by group (then patient, for readability within group blocks)
  roi_order <- order(roi_groups, roi_patient)
  mat_plot  <- mat_plot[, roi_order]
  roi_groups_ordered  <- roi_groups[roi_order]
  roi_patient_ordered <- roi_patient[roi_order]
  
  ###########################################################################
  ### ################---HEATMAP---########################################
  ###########################################################################
  
  ha_top <- HeatmapAnnotation(
    Group = roi_groups_ordered,
    
    col = list(Group = cols_group),
    
    annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
    annotation_name_side = "left",
    
    annotation_legend_param = list(
      Group = list(
        title = "Disease Group",
        title_gp = gpar(fontsize = 10, fontface = "bold"),
        labels_gp = gpar(fontsize = 9),
        grid_height = unit(4, "mm"),
        grid_width = unit(4, "mm")
      )
    ),
    
    show_legend = TRUE,
    height = unit(0.5, "cm")
  )
  
  ht <- Heatmap(
    mat_plot,
    name = "Z-score",
    
    col = col_fun_expr,
    
    cluster_columns = TRUE,
    cluster_rows = TRUE,
    clustering_distance_rows = "pearson",
    clustering_method_rows = "ward.D2",
    clustering_distance_columns = "pearson",
    clustering_method_columns = "ward.D2",
    
    show_row_dend = TRUE,
    show_column_dend = TRUE,
    row_dend_width = unit(2, "cm"),
    column_dend_height = unit(2, "cm"),
    
    top_annotation = ha_top,
    
    show_row_names = FALSE,
    show_column_names = FALSE,   
    
    column_split = NULL,
    column_gap = unit(2, "mm"),
    column_title_gp = gpar(fontsize = 11, fontface = "bold"),
    
    border = TRUE,
    rect_gp = gpar(col = NA),   
    
    heatmap_legend_param = list(
      title = "Expression\n(Z-score)",
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 9),
      grid_height = unit(4, "mm"),
      grid_width = unit(4, "mm"),
      legend_direction = "horizontal"
    ),
    
    width = unit(12, "cm"),
    height = unit(10, "cm")
  )
  
  draw(ht,
       column_title = paste0(
         "ROI-Level Expression Clustering: ", cell_type_label, "\n",
         "(", n_var_use, " most variable genes, n = ", ncol(mat_plot), " ROIs)"
       ),
       column_title_gp = gpar(fontsize = 12, fontface = "bold"),
       heatmap_legend_side = "bottom",
       merge_legend = TRUE,
       padding = unit(c(1, 1, 1, 10), "mm"))
  
  # Return data for PCA/silhouette
  return(list(
    roi_mat = roi_mat,
    roi_mat_scaled = roi_mat_scaled,
    top_var_genes = top_var_genes,
    roi_groups = roi_groups,
    roi_patient = roi_patient,
    roi_ids = roi_ids,
    n_var_genes = n_var_use
  ))
}


pdf(file.path(OUTPUT_DIR, "10a_ROI_clustering_publication.pdf"),
    width = 10, height = 8)

roi_neuron <- run_roi_clustering(
  spe_ruv,
  "Neuron",
  n_var_genes = 500
)

grid.newpage()

roi_oligo <- run_roi_clustering(
  spe_ruv,
  "Oligodendrocyte",
  n_var_genes = 500
)

dev.off()

cat("\n✓ Saved: 10a_ROI_clustering_publication.pdf\n")

###############################################################################
### PCA ################################################
###############################################################################

run_roi_pca_publication <- function(roi_results, cell_type_label) {
  
  mat <- roi_results$roi_mat[roi_results$top_var_genes, ]
  pca <- prcomp(t(mat), center = TRUE, scale. = TRUE)
  
  var_explained <- (pca$sdev^2 / sum(pca$sdev^2)) * 100
  
  pca_df <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    PC3 = pca$x[, 3],
    ROI = roi_results$roi_ids,
    Patient = roi_results$roi_patient,
    Group = roi_results$roi_groups,
    stringsAsFactors = FALSE
  )
  
  # Group centroids
  centroids <- pca_df %>%
    group_by(Group) %>%
    summarise(
      PC1_mean = mean(PC1),
      PC2_mean = mean(PC2),
      .groups = "drop"
    )
  
  p1 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, fill = Group)) +
    
    stat_ellipse(geom = "polygon", alpha = 0.1, level = 0.68,
                 show.legend = FALSE) +
    
    geom_point(size = 2, alpha = 0.5, shape = 21,
               color = "black", stroke = 0.3) +
    
    geom_point(data = centroids,
               aes(x = PC1_mean, y = PC2_mean),
               size = 6, shape = 18, color = "black") +
    
    scale_fill_manual(values = cols_group) +
    scale_color_manual(values = cols_group) +
    
    labs(
      x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
      y = paste0("PC2 (", round(var_explained[2], 1), "%)"),
      title = cell_type_label,
      subtitle = paste0("n = ", nrow(pca_df), " ROIs from ",
                        length(unique(pca_df$Patient)), " patients")
    ) +
    
    theme_classic(base_size = 11) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
      panel.grid.major = element_line(color = "grey95", linewidth = 0.3),
      
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.text = element_text(color = "black", size = 10),
      axis.title = element_text(face = "bold", size = 11),
      
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey30"),
      
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 10, face = "bold"),
      legend.key.size = unit(5, "mm"),
      
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  # Scree plot
  scree_df <- data.frame(
    PC = paste0("PC", 1:min(10, length(pca$sdev))),
    Variance = var_explained[1:min(10, length(pca$sdev))]
  )
  scree_df$PC <- factor(scree_df$PC, levels = scree_df$PC)
  
  p_scree <- ggplot(scree_df, aes(x = PC, y = Variance)) +
    geom_col(fill = "grey70", color = "black", width = 0.7) +
    geom_text(aes(label = paste0(round(Variance, 1), "%")),
              vjust = -0.5, size = 3) +
    labs(x = NULL, y = "Variance Explained (%)", title = "Scree Plot") +
    theme_classic(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )
  
  return(list(
    pca_plot = p1,
    scree_plot = p_scree,
    pca = pca,
    pca_df = pca_df,
    var_explained = var_explained
  ))
}

pca_roi_neuron <- run_roi_pca_publication(roi_neuron, "Neurons")
pca_roi_oligo  <- run_roi_pca_publication(roi_oligo, "Oligodendrocytes")

###############################################################################
###Combined PCA figure########################################################
###############################################################################

p_pca_combined <- (pca_roi_neuron$pca_plot | pca_roi_oligo$pca_plot) +
  plot_annotation(
    tag_levels = "A",
    tag_suffix = ".",
    theme = theme(plot.tag = element_text(size = 14, face = "bold"))
  )

p_pca_full <- (pca_roi_neuron$pca_plot | pca_roi_oligo$pca_plot) /
  (pca_roi_neuron$scree_plot | pca_roi_oligo$scree_plot) +
  plot_layout(heights = c(3, 1)) +
  plot_annotation(
    tag_levels = "A",
    tag_suffix = ".",
    theme = theme(plot.tag = element_text(size = 14, face = "bold"))
  )

ggsave(
  file.path(OUTPUT_DIR, "10b_ROI_PCA_publication.pdf"),
  p_pca_combined,
  width = 10, height = 5
)

ggsave(
  file.path(OUTPUT_DIR, "10b_ROI_PCA_with_scree.pdf"),
  p_pca_full,
  width = 10, height = 7
)

cat("\n✓ Saved: 10b_ROI_PCA_publication.pdf\n")
cat("✓ Saved: 10b_ROI_PCA_with_scree.pdf\n")

###############################################################################
### Silhouette analysis########################################
###############################################################################

compute_roi_silhouette_analysis <- function(roi_results, cell_type_label, n_pcs = 4) {
  
  cat("\n--- Silhouette analysis:", cell_type_label, "---\n")
  
  mat <- roi_results$roi_mat[roi_results$top_var_genes, ]
  pca <- prcomp(t(mat), center = TRUE, scale. = TRUE)
  
  d <- dist(pca$x[, 1:n_pcs, drop = FALSE])
  
  group_factor <- as.numeric(factor(roi_results$roi_groups, levels = c("PD", "MSA")))
  
  if (length(unique(group_factor)) >= 2 && min(table(group_factor)) >= 2) {
    sil <- cluster::silhouette(group_factor, d)
    sil_vals <- sil[, 3]
    sil_mean <- mean(sil_vals)
    sil_sd <- sd(sil_vals)
    
    sil_df <- data.frame(
      ROI     = roi_results$roi_ids,
      Patient = roi_results$roi_patient,
      Group   = roi_results$roi_groups,
      Silhouette = sil_vals,
      stringsAsFactors = FALSE
    ) %>%
      arrange(Group, desc(Silhouette))
    
    sil_df$ROI <- factor(sil_df$ROI, levels = sil_df$ROI)
    
    cat("  Mean silhouette:", round(sil_mean, 3), "\n")
    cat("  SD:", round(sil_sd, 3), "\n")
    cat("  Interpretation:\n")
    cat("    > 0.5:  Strong separation\n")
    cat("    0.25-0.5: Moderate separation\n")
    cat("    < 0.25: Weak/no separation\n")
    
    # ── Distribution plot ──────────────────────────
    p_sil <- ggplot(sil_df, aes(x = Group, y = Silhouette, fill = Group)) +
      geom_hline(yintercept = 0, linewidth = 0.5) +
      geom_hline(yintercept = c(0.25, 0.5), linetype = "dashed",
                 color = "grey50", linewidth = 0.4) +
      geom_violin(alpha = 0.2, color = NA, trim = TRUE) +
      geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7) +
      ggbeeswarm::geom_beeswarm(size = 1, alpha = 0.4, color = "black", cex = 1.2) +
      scale_fill_manual(values = cols_group) +
      labs(
        x = NULL,
        y = "Silhouette Coefficient",
        title = paste0("Clustering Quality: ", cell_type_label),
        subtitle = paste0("Mean = ", round(sil_mean, 3),
                          " (using first ", n_pcs, " PCs, n = ", nrow(sil_df), " ROIs)")
      ) +
      theme_classic(base_size = 11) +
      theme(
        axis.text.x = element_text(face = "bold", size = 10),
        axis.text.y = element_text(size = 9),
        axis.title = element_text(face = "bold", size = 10),
        plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5),
        legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
      )
    
    return(list(
      mean = sil_mean,
      sd = sil_sd,
      df = sil_df,
      plot = p_sil
    ))
    
  } else {
    cat("  Not enough samples for silhouette analysis\n")
    return(NULL)
  }
}

sil_roi_neuron <- compute_roi_silhouette_analysis(roi_neuron, "Neurons", n_pcs = 5)
sil_roi_oligo  <- compute_roi_silhouette_analysis(roi_oligo, "Oligodendrocytes", n_pcs = 5)

if (!is.null(sil_roi_neuron) && !is.null(sil_roi_oligo)) {
  
  p_sil_combined <- (sil_roi_neuron$plot | sil_roi_oligo$plot) +
    plot_annotation(
      tag_levels = "A",
      tag_suffix = ".",
      theme = theme(plot.tag = element_text(size = 14, face = "bold"))
    )
  
  ggsave(
    file.path(OUTPUT_DIR, "10c_ROI_silhouette_analysis_publication.pdf"),
    p_sil_combined,
    width = 10, height = 6
  )
  
  cat("\n✓ Saved: 10c_ROI_silhouette_analysis_publication.pdf\n")
}

# Summary table
sil_roi_summary <- data.frame(
  Cell_Type = c("Neurons", "Oligodendrocytes"),
  N_ROIs = c(
    if(!is.null(sil_roi_neuron)) nrow(sil_roi_neuron$df) else NA,
    if(!is.null(sil_roi_oligo))  nrow(sil_roi_oligo$df)  else NA
  ),
  Mean_Silhouette = c(
    if(!is.null(sil_roi_neuron)) round(sil_roi_neuron$mean, 3) else NA,
    if(!is.null(sil_roi_oligo))  round(sil_roi_oligo$mean, 3)  else NA
  ),
  SD_Silhouette = c(
    if(!is.null(sil_roi_neuron)) round(sil_roi_neuron$sd, 3) else NA,
    if(!is.null(sil_roi_oligo))  round(sil_roi_oligo$sd, 3)  else NA
  ),
  Interpretation = c(
    if(!is.null(sil_roi_neuron)) {
      ifelse(sil_roi_neuron$mean > 0.5, "Strong",
             ifelse(sil_roi_neuron$mean > 0.25, "Moderate", "Weak"))
    } else NA,
    if(!is.null(sil_roi_oligo)) {
      ifelse(sil_roi_oligo$mean > 0.5, "Strong",
             ifelse(sil_roi_oligo$mean > 0.25, "Moderate", "Weak"))
    } else NA
  )
)

write.csv(sil_roi_summary,
          file.path(OUTPUT_DIR, "10_ROI_silhouette_summary.csv"),
          row.names = FALSE)

cat("\n")
print(sil_roi_summary)

cat("\n=== ROI-Level Clustering Analysis Complete ===\n")
cat("Files created:\n")
cat("  10a_ROI_clustering_publication.pdf (heatmaps)\n")
cat("  10b_ROI_PCA_publication.pdf (PCA plots)\n")
cat("  10b_ROI_PCA_with_scree.pdf (PCA + scree plots)\n")
cat("  10c_ROI_silhouette_analysis_publication.pdf (clustering quality)\n")
cat("  10_ROI_silhouette_summary.csv (summary table)\n")