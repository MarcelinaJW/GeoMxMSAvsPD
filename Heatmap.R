###############################################################################
### DEG Heatmap
###############################################################################

library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(grid)
library(RColorBrewer)

cat("\n=== DEG Heatmap ===\n")


# Color schemes
cols_group <- c(PD = "#4C72B0", MSA = "#C44E52")
cols_cell  <- c(Neuron = "#E69F00", Oligodendrocyte = "#009E73")

patients <- unique(spe_qc$Patient)
cols_patient <- setNames(
  colorRampPalette(c("grey95", "grey30"))(length(patients)),
  patients
)

# Expression colors
col_fun_expr <- colorRamp2(
  c(-3, -1.5, 0, 1.5, 3),
  c("#2166AC", "#92C5DE", "white", "#F4A582", "#B2182B")
)

###############################################################################
### ALL ROIs HEATMAP
###############################################################################

cat("\n--- Creating ALL ROIs heatmap ---\n")

# Load DE results
de_all <- read.csv(file.path(OUTPUT_DIR, "DE_MSA_vs_PD_allROIs.csv"), row.names = 1)

# Select top DEGs
top_degs_all <- de_all %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::filter(adj.P.Val < 0.05) %>%
  dplyr::mutate(
    rank_score = -log10(adj.P.Val) * abs(logFC),
    direction = ifelse(logFC > 0, "MSA-up", "PD-up")
  ) %>%
  dplyr::arrange(desc(rank_score))

top_msa <- top_degs_all %>% dplyr::filter(direction == "MSA-up") %>% head(19)
top_pd  <- top_degs_all %>% dplyr::filter(direction == "PD-up") %>% head(19)

top_genes_all <- c(top_msa$gene, top_pd$gene)

cat("Selected genes (All ROIs):\n")
cat("  MSA-enriched:", nrow(top_msa), "\n")
cat("  PD-enriched:", nrow(top_pd), "\n")

# Expression data
expr_mat <- assay(spe_ruv, "logcounts")
expr_subset <- expr_mat[top_genes_all, ]

# Z-score
expr_z <- t(scale(t(expr_subset)))
expr_z[expr_z > 3] <- 3
expr_z[expr_z < -3] <- -3

# Order samples
sample_order <- order(spe_ruv$Group, spe_ruv$Patient, spe_ruv$Cell)
expr_z <- expr_z[, sample_order]

# Gene annotation (logFC + significance)
gene_annotation_all <- data.frame(
  gene = top_genes_all,
  logFC = de_all[top_genes_all, "logFC"],
  adj.P.Val = de_all[top_genes_all, "adj.P.Val"],
  stringsAsFactors = FALSE
)

###############################################################################
### Annotation — Effect Size + Significance
###############################################################################

ha_right_all <- rowAnnotation(
  `Effect Size` = anno_barplot(
    gene_annotation_all$logFC,
    gp = gpar(fill = ifelse(gene_annotation_all$logFC > 0, "#C44E52", "#4C72B0")),
    width = unit(2, "cm"),
    axis_param = list(gp = gpar(fontsize = 7))
  ),
  
  `Significance` = anno_barplot(
    -log10(gene_annotation_all$adj.P.Val),
    gp = gpar(fill = "grey40"),
    width = unit(1.5, "cm"),
    axis_param = list(gp = gpar(fontsize = 7))
  ),
  
  annotation_name_gp = gpar(fontsize = 9, fontface = "bold"),
  annotation_name_side = "top",
  annotation_name_rot = 0,
  
  gap = unit(2, "mm")
)

###############################################################################
### Top annotation
###############################################################################

ha_top_all <- HeatmapAnnotation(
  Group   = spe_ruv$Group[sample_order],
  Cell    = spe_ruv$Cell[sample_order],
  Patient = spe_ruv$Patient[sample_order],
  
  col = list(
    Group   = cols_group,
    Cell    = cols_cell,
    Patient = cols_patient
  ),
  
  annotation_name_gp = gpar(fontsize = 9, fontface = "bold"),
  annotation_name_side = "left",
  annotation_name_rot = 0,
  
  gap = unit(0.5, "mm"),
  
  annotation_legend_param = list(
    Group = list(
      title = "Group",
      title_gp = gpar(fontsize = 9, fontface = "bold"),
      labels_gp = gpar(fontsize = 8),
      grid_height = unit(3, "mm"),
      grid_width = unit(4, "mm")
    ),
    Cell = list(
      title = "Cell Type",
      title_gp = gpar(fontsize = 9, fontface = "bold"),
      labels_gp = gpar(fontsize = 8),
      grid_height = unit(3, "mm"),
      grid_width = unit(4, "mm")
    ),
    Patient = list(
      title = "Patient",
      title_gp = gpar(fontsize = 9, fontface = "bold"),
      labels_gp = gpar(fontsize = 7),
      grid_height = unit(2, "mm"),
      grid_width = unit(4, "mm"),
      nrow = 2
    )
  ),
  
  show_legend = c(TRUE, TRUE, FALSE)
)

###############################################################################
### Heatmap
###############################################################################

ht_all <- Heatmap(
  expr_z,
  name = "Expression\n(Z-score)",
  
  col = col_fun_expr,
  
  cluster_rows = TRUE,
  cluster_row_slices = TRUE,
  clustering_distance_rows = "pearson",
  clustering_method_rows = "ward.D2",
  show_row_dend = TRUE,
  row_dend_width = unit(1.5, "cm"),
  
  cluster_columns = TRUE,
  cluster_column_slices = FALSE,
  clustering_distance_columns = "euclidean",
  clustering_method_columns = "ward.D2",
  show_column_dend = FALSE,
  
  row_split = NULL,
  row_gap = unit(3, "mm"),
  
  right_annotation = ha_right_all,
  top_annotation = ha_top_all,
  
  row_names_gp = gpar(fontsize = 8, fontface = "italic"),
  row_names_side = "left",
  
  show_column_names = FALSE,
  
  row_title_gp = gpar(fontsize = 11, fontface = "bold"),
  row_title_rot = 0,
  row_title_side = "left",
  
  column_split = spe_ruv$Group[sample_order],  # visually group MSA | PD blocks
  column_gap = unit(1.5, "mm"),
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  
  border = TRUE,
  rect_gp = gpar(col = "white", lwd = 0.5),
  
  heatmap_legend_param = list(
    title = "Expression\n(Z-score)",
    title_gp = gpar(fontsize = 9, fontface = "bold"),
    labels_gp = gpar(fontsize = 8),
    grid_height = unit(4, "mm"),
    grid_width = unit(4, "mm"),
    legend_direction = "horizontal"
  ),
  
  width = unit(14, "cm"),
  height = unit(14, "cm"),
  
  row_names_max_width = unit(5, "cm")
)

# Save
pdf(file.path(OUTPUT_DIR, "09a_heatmap_top_DEGs_publication.pdf"),
    width = 12, height = 8)

draw(ht_all,
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom",
     merge_legend = TRUE,
     legend_gap = unit(0.5, "cm"),
     padding = unit(c(1, 1, 1, 10), "mm"))

grid.text(
  "Figure X. Differential gene expression in MSA vs PD substantia nigra",
  x = 0.5, y = 0.02,
  gp = gpar(fontsize = 9, fontface = "bold"),
  just = "center"
)

dev.off()

cat("✓ Saved: 09a_heatmap_top_DEGs_publication.pdf\n")
###############################################################################
### CELL-TYPE SPECIFIC HEATMAPS
###############################################################################

create_celltype_heatmap <- function(de_obj, spe_obj, celltype_name) {
  
  cat("\n--- Creating", celltype_name, "heatmap ---\n")
  
  ### DEG selection
  top_degs <- de_obj %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::filter(adj.P.Val < 0.05) %>%
    dplyr::mutate(
      rank_score = -log10(adj.P.Val) * abs(logFC),
      direction = ifelse(logFC > 0, "MSA-up", "PD-up")
    ) %>%
    dplyr::arrange(desc(rank_score))
  
  top_msa <- top_degs %>% dplyr::filter(direction == "MSA-up") %>% head(20)
  top_pd  <- top_degs %>% dplyr::filter(direction == "PD-up") %>% head(20)
  
  top_genes <- c(top_msa$gene, top_pd$gene)
  top_genes <- intersect(top_genes, rownames(spe_obj))
  
  cat("  Selected genes:", length(top_genes), "\n")
  
  ### Expression
  expr <- assay(spe_obj, "logcounts")[top_genes, ]
  expr_z <- t(scale(t(expr)))
  expr_z[expr_z > 3] <- 3
  expr_z[expr_z < -3] <- -3
  
  sample_order <- order(spe_obj$Group, spe_obj$Patient)
  expr_z <- expr_z[, sample_order]
  
  ### Gene annotation
  gene_annotation <- data.frame(
    gene = top_genes,
    logFC = de_obj[top_genes, "logFC"],
    adj.P.Val = de_obj[top_genes, "adj.P.Val"],
    stringsAsFactors = FALSE
  )
  
  ### Right annotation
  ha_right <- rowAnnotation(
    `Effect Size` = anno_barplot(
      gene_annotation$logFC,
      gp = gpar(fill = ifelse(gene_annotation$logFC > 0, "#C44E52", "#4C72B0")),
      width = unit(2, "cm"),
      axis_param = list(gp = gpar(fontsize = 7))
    ),
    
    `Significance` = anno_barplot(
      -log10(gene_annotation$adj.P.Val),
      gp = gpar(fill = "grey40"),
      width = unit(1.5, "cm"),
      axis_param = list(gp = gpar(fontsize = 7))
    ),
    
    annotation_name_gp = gpar(fontsize = 9, fontface = "bold"),
    annotation_name_side = "top",
    annotation_name_rot = 0,
    
    gap = unit(2, "mm")
  )
  
  ### Top annotation
  ha_top <- HeatmapAnnotation(
    Group   = spe_obj$Group[sample_order],
    Patient = spe_obj$Patient[sample_order],
    
    col = list(
      Group   = cols_group,
      Patient = cols_patient
    ),
    
    annotation_name_gp = gpar(fontsize = 9, fontface = "bold"),
    annotation_name_side = "left",
    annotation_name_rot = 0,
    
    gap = unit(0.5, "mm"),
    
    annotation_legend_param = list(
      Group = list(
        title = "Group",
        title_gp = gpar(fontsize = 9, fontface = "bold"),
        labels_gp = gpar(fontsize = 8),
        grid_height = unit(3, "mm"),
        grid_width = unit(4, "mm")
      ),
      Patient = list(
        title = "Patient",
        title_gp = gpar(fontsize = 9, fontface = "bold"),
        labels_gp = gpar(fontsize = 7),
        grid_height = unit(2, "mm"),
        grid_width = unit(4, "mm"),
        nrow = 2
      )
    ),
    
    show_legend = c(TRUE, FALSE)
  )
  
  ### Heatmap
  ht <- Heatmap(
    expr_z,
    name = "Expression\n(Z-score)",
    
    col = col_fun_expr,
    
    cluster_rows = TRUE,
    cluster_row_slices = TRUE,
    clustering_distance_rows = "pearson",
    clustering_method_rows = "ward.D2",
    show_row_dend = TRUE,
    row_dend_width = unit(1.5, "cm"),
    show_column_dend = FALSE,
    
    cluster_columns = TRUE,               # cluster WITHIN each Group block
    cluster_column_slices = FALSE,        
    clustering_distance_columns = "euclidean",
    clustering_method_columns = "ward.D2",
    
    row_split = NULL,
    row_gap = unit(3, "mm"),
    
    right_annotation = ha_right,
    top_annotation = ha_top,
    
    row_names_gp = gpar(fontsize = 8, fontface = "italic"),
    row_names_side = "left",
    
    show_column_names = FALSE,
    
    row_title_gp = gpar(fontsize = 11, fontface = "bold"),
    row_title_rot = 0,
    row_title_side = "left",
    
    column_split = spe_obj$Group[sample_order],   
    column_gap = unit(1.5, "mm"),
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    
    border = TRUE,
    rect_gp = gpar(col = "white", lwd = 0.5),
    
    heatmap_legend_param = list(
      title = "Expression\n(Z-score)",
      title_gp = gpar(fontsize = 9, fontface = "bold"),
      labels_gp = gpar(fontsize = 8),
      grid_height = unit(4, "mm"),
      grid_width = unit(4, "mm"),
      legend_direction = "horizontal"
    ),
    
    width = unit(14, "cm"),
    height = unit(14, "cm"),
    
    row_names_max_width = unit(5, "cm")
  )
  
  return(ht)
}

# Create cell-type heatmaps
ht_neuron <- create_celltype_heatmap(de_neuron, spe_neuron, "Neurons")
ht_oligo  <- create_celltype_heatmap(de_oligo, spe_oligo, "Oligodendrocytes")

# Save
pdf(file.path(OUTPUT_DIR, "09a_heatmap_DEGs_by_celltype.pdf"),
    width = 12, height = 8)

draw(ht_neuron,
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom",
     merge_legend = TRUE,
     legend_gap = unit(0.5, "cm"),
     padding = unit(c(1, 1, 1, 10), "mm"))

grid.text(
  "Differential gene expression in MSA vs PD substantia nigra (Neurons)",
  x = 0.5, y = 0.02,
  gp = gpar(fontsize = 9, fontface = "bold"),
  just = "center"
)

draw(ht_oligo,
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom",
     merge_legend = TRUE,
     legend_gap = unit(0.5, "cm"),
     padding = unit(c(1, 1, 1, 10), "mm"),
     newpage = TRUE)

grid.text(
  "Differential gene expression in MSA vs PD substantia nigra (Oligodendrocytes)",
  x = 0.5, y = 0.02,
  gp = gpar(fontsize = 9, fontface = "bold"),
  just = "center"
)

dev.off()

cat("✓ Saved: 09a_heatmap_DEGs_by_celltype.pdf\n")
cat("\n=== All publication heatmaps complete ===\n")