###############################################################################
### PD-Associated Genes and Cell-Type Enrichment
###############################################################################
cat("\n--- Creating PD gene cell-type enrichment plot ---\n")

pd_candidates_plot <- intersect(pd_candidates_present, rownames(de_celltype))

plot_df <- de_celltype[pd_candidates_plot, , drop = FALSE] %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  mutate(
    CellType = ifelse(logFC > 0, "Neuron-enriched", "Oligodendrocyte-enriched"),
    Gene     = reorder(Gene, logFC),
    sig      = p_to_stars(adj.P.Val),
    label_x  = logFC + sign(logFC) * (max(abs(logFC)) * 0.06)
  )

cols_cell_enrichment <- c(
  "Neuron-enriched"          = unname(cols_cell["Neuron"]),
  "Oligodendrocyte-enriched" = unname(cols_cell["Oligodendrocyte"])
)

x_max <- max(abs(plot_df$logFC)) * 1.25

p_pd_enrichment <- ggplot(plot_df, aes(x = logFC, y = Gene, fill = CellType)) +
  geom_vline(xintercept = 0, linewidth = 0.5, color = "grey40") +
  geom_col(width = 0.65, color = "black", linewidth = 0.3) +
  geom_text(aes(x = label_x, label = sig), fontface = "bold", size = 3.6, hjust = 0.5) +
  scale_fill_manual(values = cols_cell_enrichment, name = NULL) +   # FIXED
  scale_x_continuous(limits = c(-x_max, x_max), expand = expansion(mult = 0.02)) +
  labs(
    title    = "Cell-Type Enrichment of PD-Associated Genes",
    subtitle = expression("Differential expression: Neuron "*italic("vs")*" Oligodendrocyte"),
    x = expression(bold("log"[2]*" Fold Change (Neuron / Oligodendrocyte)")),
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 13, hjust = 0),
    plot.subtitle = element_text(size = 9.5, color = "grey30", hjust = 0, margin = margin(b = 10)),
    axis.title.x = element_text(face = "bold", size = 11, margin = margin(t = 8)),
    axis.text.x  = element_text(color = "black"),
    axis.text.y  = element_text(face = "bold.italic", size = 10.5, color = "black"),
    axis.line.y  = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "top",
    legend.justification = "left",
    legend.margin = margin(b = 4),
    legend.text = element_text(size = 9.5, face = "bold"),
    legend.key.size = unit(4, "mm"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.major.y = element_blank(),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave(
  file.path(OUTPUT_DIR, "08j_PD_gene_celltype_enrichment.pdf"),
  p_pd_enrichment, width = 6, height = 7, device = pdf
)

cat("  08j_PD_gene_celltype_enrichment.pdf\n")










###############################################################################
### Comparison Heatmap comparing All-ROIs DEGs across All ROIs / Neuron / Oligo contrasts
###############################################################################
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(grid)

cat("\n=== Creating cross-contrast logFC comparison heatmap ===\n")

###############################################################################
###Build logFC + significance matrix for the All-ROIs top DEGs
###############################################################################

sig_cut <- 0.05

top_genes_all <- de_all %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  filter(adj.P.Val < 0.05,
         abs(logFC) > 0.5) %>%
  arrange(desc(logFC)) %>%
  pull(gene)

# Pull logFC and adj.P.Val for the same gene set from all three DE results
comparison_df <- data.frame(gene = top_genes_all, stringsAsFactors = FALSE) %>%
  
  left_join(
    de_all %>% as.data.frame() %>% tibble::rownames_to_column("gene") %>%
      dplyr::select(gene, logFC_all = logFC, padj_all = adj.P.Val),
    by = "gene"
  ) %>%
  
  left_join(
    de_neuron %>% as.data.frame() %>% tibble::rownames_to_column("gene") %>%
      dplyr::select(gene, logFC_neuron = logFC, padj_neuron = adj.P.Val),
    by = "gene"
  ) %>%
  
  left_join(
    de_oligo %>% as.data.frame() %>% tibble::rownames_to_column("gene") %>%
      dplyr::select(gene, logFC_oligo = logFC, padj_oligo = adj.P.Val),
    by = "gene"
  )

# Order genes: MSA-up first (by original All-ROIs rank), then PD-up
comparison_df <- comparison_df %>%
  mutate(gene = factor(gene, levels = top_genes_all)) %>%
  arrange(gene)

###############################################################################
### Build logFC matrix (genes x 3 contrasts)
###############################################################################

logfc_mat <- comparison_df %>%
  dplyr::select(gene, logFC_all, logFC_neuron, logFC_oligo) %>%
  tibble::column_to_rownames("gene") %>%
  as.matrix()

colnames(logfc_mat) <- c("All ROIs", "Neuron", "Oligodendrocyte")

###############################################################################
### Build significance mask matrix (TRUE = significant, FALSE = not sig / NA)
###############################################################################

sig_mat <- comparison_df %>%
  dplyr::select(gene, padj_all, padj_neuron, padj_oligo) %>%
  tibble::column_to_rownames("gene") %>%
  as.matrix()

colnames(sig_mat) <- c("All ROIs", "Neuron", "Oligodendrocyte")

sig_mask <- !is.na(sig_mat) & sig_mat < sig_cut

# Keep original logFC for significant cells; set NA for non-significant
# (NA cells will be grey)
logfc_mat_masked <- logfc_mat
logfc_mat_masked[!sig_mask] <- NA

###############################################################################
###Colour function
###############################################################################

max_abs_lfc <- max(abs(logfc_mat), na.rm = TRUE)

col_fun_lfc <- colorRamp2(
  c(-1, -0.5, 0, 0.5, 1),
  c("#053061", "#92C5DE", "white",
    "#F4A582", "#67001F")
)

###############################################################################
### Cell function which adds significance stars and handles NA cells
###############################################################################

cell_fun_lfc <- function(j, i, x, y, width, height, fill) {
  
  val <- logfc_mat[i, j]
  is_sig <- sig_mask[i, j]
  
  if (is.na(val)) {
    grid.rect(x, y, width, height, gp = gpar(fill = "grey85", col = "white", lwd = 0.5))
  } else if (!is_sig) {
    # Non-significant — grey fill, no colour scale
    grid.rect(x, y, width, height, gp = gpar(fill = "grey85", col = "white", lwd = 0.5))
  } else {
    # Significant — coloured by logFC, with a star
    grid.rect(x, y, width, height, gp = gpar(fill = col_fun_lfc(val), col = "white", lwd = 0.5))
  }
}

###############################################################################
### Build heatmap
###############################################################################
ht_comparison <- Heatmap(
  logfc_mat,                      # full matrix (values used internally for clustering/colour scale range)
  name = "log2FC",
  
  col = col_fun_lfc,
  na_col = "grey85",
  
  cell_fun = cell_fun_lfc,
  
  # Rows, keep in original MSA-up/PD-up rank order
  cluster_rows = FALSE,
  
  # Columns with fixed order, no clustering
  cluster_columns = FALSE,
  
  row_names_gp = gpar(fontsize = 8, fontface = "italic"),
  row_names_side = "left",
  
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_rot = 0,
  column_names_centered = TRUE,
  
  border = TRUE,
  rect_gp = gpar(col = NA),   
  
  heatmap_legend_param = list(
    title = expression(log[2]*"FC"),
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 9),
    grid_height = unit(4, "mm"),
    grid_width = unit(4, "mm"),
    legend_direction = "vertical"
  ),
  
  width = unit(8, "cm"),
  height = unit(0.35 * nrow(logfc_mat), "cm"),
  
  row_names_max_width = unit(5, "cm")
)

###############################################################################
### Save
###############################################################################

pdf(file.path(OUTPUT_DIR, "09b_DEG_comparison_across_contrasts.pdf"),
    width = 6, height = 0.3 * nrow(logfc_mat) + 3)

draw(ht_comparison,
     column_title = "All-ROIs DEGs: Effect Size Across Contrasts",
     column_title_gp = gpar(fontsize = 12, fontface = "bold"),
     heatmap_legend_side = "right",
     padding = unit(c(2, 2, 2, 10), "mm"))

grid.text(
  "* FDR < 0.05 in that contrast. Grey = not significant.",
  x = 0.5, y = 0.01,
  gp = gpar(fontsize = 8, fontface = "italic"),
  just = "center"
)

dev.off()

cat("✓ Saved: 09b_DEG_comparison_across_contrasts.pdf\n")