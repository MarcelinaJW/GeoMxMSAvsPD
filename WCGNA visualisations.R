###############################################################################
### WGCNA Visualizations
###############################################################################
library(WGCNA)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggbeeswarm)
library(igraph)
library(scales)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
library(stringr)
library(ggraph)
library(tidygraph)

cat("\n=== Creating WGCNA Plots ===\n")
###############################################################################
### 1. Power Selection Plot
###############################################################################
cat("\n--- Creating power selection plot ---\n")

sft_df <- read.csv(file.path(OUTPUT_DIR, "WGCNA_power_selection.csv"))

p_power_fit <- ggplot(sft_df, aes(x = Power, y = SFT.R.sq)) +
  geom_hline(yintercept = c(0.6, 0.7, 0.8), 
             linetype = c("dotted", "dashed", "solid"),
             color = c("grey50", "orange", "red"), 
             linewidth = 0.5) +
  geom_line(color = "#0072B2", linewidth = 1) +
  geom_point(aes(size = selected), 
             color = ifelse(sft_df$selected, "#D55E00", "#0072B2"),
             alpha = 0.8) +
  scale_size_manual(values = c("FALSE" = 3, "TRUE" = 6), guide = "none") +
  annotate("text", x = max(sft_df$Power) * 0.7, y = 0.8,
           label = "R² = 0.8", vjust = -0.5, size = 3.5, color = "red") +
  labs(
    x = "Soft Threshold (Power)",
    y = expression(bold("Scale Free Topology Model Fit (signed R"^2*")")),
    title = "Scale-Free Topology Fit Index"
  ) +
  theme_classic(base_size = 11) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    panel.grid.major = element_line(color = "grey95", linewidth = 0.3),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5)
  )

p_power_conn <- ggplot(sft_df, aes(x = Power, y = MeanConnectivity)) +
  geom_line(color = "#D55E00", linewidth = 1) +
  geom_point(aes(size = selected),
             color = ifelse(sft_df$selected, "#D55E00", "#0072B2"),
             alpha = 0.8) +
  scale_size_manual(values = c("FALSE" = 3, "TRUE" = 6), guide = "none") +
  scale_y_log10(breaks = c(1, 10, 100, 1000)) +
  labs(
    x = "Soft Threshold (Power)",
    y = "Mean Connectivity (log scale)",
    title = "Mean Network Connectivity"
  ) +
  theme_classic(base_size = 11) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    panel.grid.major = element_line(color = "grey95", linewidth = 0.3),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5)
  )

p_power_combined <- p_power_fit + p_power_conn +
  plot_annotation(
    tag_levels = "A",
    tag_suffix = ".",
    theme = theme(plot.tag = element_text(size = 14, face = "bold"))
  )

ggsave(
  file.path(OUTPUT_DIR, "11a_WGCNA_power_selection.pdf"),
  p_power_combined,
  width = 10,
  height = 4.5
)

cat("✓ Saved: 11a_WGCNA_power_selection.pdf\n")

###############################################################################
### 2. Module Sizes Barplot
###############################################################################
cat("\n--- Creating module sizes plot ---\n")

module_sizes <- table(module_colors)
module_sizes <- sort(module_sizes[names(module_sizes) != "grey"], decreasing = TRUE)

df_sizes <- data.frame(
  Module = names(module_sizes),
  Count = as.numeric(module_sizes)
)
df_sizes$Module <- factor(df_sizes$Module, levels = df_sizes$Module)

# Color mapping
color_map <- setNames(as.character(df_sizes$Module), df_sizes$Module)

p_sizes <- ggplot(df_sizes, aes(x = Module, y = Count, fill = Module)) +
  geom_col(color = "black", linewidth = 0.3, width = 0.7) +
  geom_text(aes(label = Count), vjust = -0.5, size = 3, fontface = "bold") +
  scale_fill_manual(values = color_map) +
  labs(
    x = "Module",
    y = "Number of Genes",
    title = "Gene Module Sizes"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 9),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

ggsave(
  file.path(OUTPUT_DIR, "11b_WGCNA_module_sizes.pdf"),
  p_sizes,
  width = 10,
  height = 6
)

cat("✓ Saved: 11b_WGCNA_module_sizes.pdf\n")

###############################################################################
### 3. Module-Trait Correlation Heatmap
###############################################################################
cat("\n--- Creating module-trait heatmap ---\n")

# Load results
module_trait_results <- read.csv(file.path(OUTPUT_DIR, "WGCNA_module_trait_results.csv"))

# Prepare correlation matrix for heatmap, Use patient-level correlations
cor_mat <- module_trait_cor
pval_mat <- module_trait_pval

# Remove ME prefix
rownames(cor_mat) <- gsub("^ME", "", rownames(cor_mat))
rownames(pval_mat) <- gsub("^ME", "", rownames(pval_mat))

# Keep only Group columns
keep_cols <- c("Group_MSA", "Group_PD")
cor_mat <- cor_mat[, keep_cols, drop = FALSE]
pval_mat <- pval_mat[, keep_cols, drop = FALSE]

colnames(cor_mat) <- gsub("Group_", "", colnames(cor_mat))
colnames(pval_mat) <- gsub("Group_", "", colnames(pval_mat))

# Create text matrix
text_mat <- matrix(
  paste0(
    sprintf("%.2f", cor_mat),
    "\n",
    ifelse(pval_mat < 0.001, "***",
           ifelse(pval_mat < 0.01, "**",
                  ifelse(pval_mat < 0.05, "*", "")))
  ),
  nrow = nrow(cor_mat),
  ncol = ncol(cor_mat)
)

rownames(text_mat) <- rownames(cor_mat)
colnames(text_mat) <- colnames(cor_mat)

# Color function
col_fun <- colorRamp2(
  seq(-0.8, 0.8, length.out = 9),
  c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", 
    "white", 
    "#FDDBC7", "#F4A582", "#D6604D", "#B2182B")
)

# Create heatmap
ht_trait <- Heatmap(
  cor_mat,
  name = "Correlation",
  
  col = col_fun,
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(text_mat[i, j], x, y, 
              gp = gpar(fontsize = 10, fontface = "bold"))
  },
  
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "ward.D2",
  show_row_dend = TRUE,
  row_dend_width = unit(2, "cm"),
  
  row_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_gp = gpar(fontsize = 11, fontface = "bold"),
  column_names_rot = 0,
  column_names_centered = TRUE,
  
  border = TRUE,
  rect_gp = gpar(col = "white", lwd = 1),
  
  width = unit(5, "cm"),
  height = unit(12, "cm"),
  
  # Legend
  heatmap_legend_param = list(
    title = "Correlation",
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 9),
    grid_height = unit(4, "mm"),
    grid_width = unit(4, "mm"),
    legend_direction = "vertical",
    title_position = "topcenter",
    at = seq(-0.8, 0.8, 0.4)
  )
)

pdf(file.path(OUTPUT_DIR, "11c_WGCNA_module_trait_heatmap.pdf"),
    width = 6, height = 10)

draw(ht_trait,
     column_title = "Module-Disease Group Correlations",
     column_title_gp = gpar(fontsize = 13, fontface = "bold"),
     heatmap_legend_side = "right",
     padding = unit(c(2, 2, 2, 10), "mm"))

dev.off()

cat("✓ Saved: 11c_WGCNA_module_trait_heatmap.pdf\n")

###############################################################################
### 4. Module Eigengene Expression Plots
###############################################################################
cat("\n--- Creating module eigengene plots ---\n")

# Get top 6 significant modules
top_modules <- head(module_trait_results$Module, 6)
top_modules <- gsub("^ME", "", top_modules)

ME_long <- MEs %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Patient") %>%
  tidyr::pivot_longer(
    cols = starts_with("ME"),
    names_to = "Module",
    values_to = "Eigengene"
  ) %>%
  mutate(Module = gsub("^ME", "", Module)) %>%
  left_join(patient_meta, by = "Patient")

# Filter to top modules
ME_plot_data <- ME_long %>%
  filter(Module %in% top_modules) %>%
  mutate(Module = factor(Module, levels = top_modules))

# Summary stats
ME_summary <- ME_plot_data %>%
  group_by(Module, Group) %>%
  summarise(
    mean = mean(Eigengene, na.rm = TRUE),
    sem = sd(Eigengene, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# Get p-values
ME_pvals <- module_trait_results %>%
  mutate(Module = gsub("^ME", "", Module)) %>%
  filter(Module %in% top_modules) %>%
  select(Module, P.Value) %>%
  mutate(
    p_label = case_when(
      P.Value < 0.001 ~ "***",
      P.Value < 0.01 ~ "**",
      P.Value < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    y_pos = Inf
  )

# Plot
p_eigengenes <- ggplot(ME_plot_data, aes(x = Group, y = Eigengene, fill = Group)) +
  
  geom_violin(alpha = 0.2, color = NA, scale = "width", trim = TRUE) +
  
  geom_boxplot(width = 0.3, outlier.shape = NA, alpha = 0.7) +
  
  ggbeeswarm::geom_beeswarm(alpha = 0.6, size = 1.5, cex = 2) +
  
  geom_errorbar(data = ME_summary,
                aes(x = Group, y = mean, ymin = mean - sem, ymax = mean + sem),
                width = 0.2, linewidth = 0.7, color = "black",
                inherit.aes = FALSE) +
  geom_point(data = ME_summary,
             aes(x = Group, y = mean),
             size = 2.5, shape = 18, color = "black",
             inherit.aes = FALSE) +
  
  # Statistical annotation
  geom_text(data = ME_pvals,
            aes(x = 1.5, y = y_pos, label = p_label),
            vjust = 1.2, size = 5, fontface = "bold",
            inherit.aes = FALSE) +
  
  facet_wrap(~Module, scales = "fixed", ncol = 6) +
  
  scale_fill_manual(values = cols_group) +
  
  labs(
    x = NULL,
    y = "Module Eigengene",
    title = "Module Eigengene Expression by Disease"
  ) +
  
  # Theme
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey95", color = "black", linewidth = 0.5),
    strip.text = element_text(face = "bold", size = 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    axis.text = element_text(color = "black", size = 9),
    axis.title = element_text(face = "bold", size = 11),
    axis.text.x = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    legend.position = "none"
  )

ggsave(
  file.path(OUTPUT_DIR, "11d_WGCNA_module_eigengenes.pdf"),
  p_eigengenes,
  width = 10,
  height = 3
)

cat("✓ Saved: 11d_WGCNA_module_eigengenes.pdf\n")

###############################################################################
### 5. Module Eigengene Correlation Network
###############################################################################
cat("\n--- Creating eigengene network ---\n")

MEs_ordered <- orderMEs(MEs)
ME_cor <- cor(MEs_ordered, use = "pairwise.complete.obs", method = "pearson")

# Remove ME prefix
colnames(ME_cor) <- gsub("^ME", "", colnames(ME_cor))
rownames(ME_cor) <- gsub("^ME", "", rownames(ME_cor))

# Tree for clustering
ME_tree <- hclust(as.dist(1 - ME_cor), method = "average")

col_fun_cor <- colorRamp2(
  seq(-1, 1, length.out = 9),
  c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
    "white",
    "#FDDBC7", "#F4A582", "#D6604D", "#B2182B")
)

# Heatmap
ht_cor <- Heatmap(
  ME_cor,
  name = "Pearson r",
  
  col = col_fun_cor,
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (abs(ME_cor[i, j]) > 0.3) {  # Only label strong correlations
      grid.text(sprintf("%.2f", ME_cor[i, j]), x, y,
                gp = gpar(fontsize = 7))
    }
  },
  
  cluster_rows = ME_tree,
  cluster_columns = ME_tree,
  show_row_dend = TRUE,
  show_column_dend = TRUE,
  row_dend_width = unit(2, "cm"),
  column_dend_height = unit(2, "cm"),
  
  row_names_gp = gpar(fontsize = 9, fontface = "bold"),
  column_names_gp = gpar(fontsize = 9, fontface = "bold"),
  
  border = TRUE,
  rect_gp = gpar(col = "white", lwd = 0.5),
  
  heatmap_legend_param = list(
    title = "Pearson r",
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 9),
    grid_height = unit(4, "mm"),
    grid_width = unit(4, "mm"),
    legend_direction = "vertical",
    at = seq(-1, 1, 0.5)
  )
)

pdf(file.path(OUTPUT_DIR, "11e_WGCNA_eigengene_network.pdf"),
    width = 9, height = 8)

draw(ht_cor,
     column_title = "Module Eigengene Correlations",
     column_title_gp = gpar(fontsize = 13, fontface = "bold"),
     heatmap_legend_side = "right",
     padding = unit(c(2, 2, 2, 10), "mm"))

dev.off()

cat("✓ Saved: 11e_WGCNA_eigengene_network.pdf\n")

###############################################################################
### 6. Hub Gene Heatmaps (Top 6 Modules)
###############################################################################
cat("\n--- Creating hub gene heatmaps ---\n")

# Load hub genes must have SYMBOL column
top_6_modules <- head(gsub("^ME", "", module_trait_results$Module), 16)

pdf(file.path(OUTPUT_DIR, "11f_WGCNA_hub_gene_heatmaps.pdf"),
    width = 10, height = 8)

for (mod in top_6_modules) {
  
  hub_file <- file.path(OUTPUT_DIR, paste0("WGCNA_hub_genes_", mod, ".csv"))
  
  if (!file.exists(hub_file)) {
    cat("  Hub gene file not found for module", mod, "\n")
    next
  }
  
  hub_genes <- read.csv(hub_file)
  top_hubs <- head(hub_genes, 25)
  
  cat("  Creating heatmap for module", mod, "-", nrow(top_hubs), "hub genes\n")
  
  # Get expression
  genes_use <- intersect(top_hubs$ENSEMBL, rownames(pseudobulk_mat))
  
  if (length(genes_use) < 5) {
    cat("    Not enough genes found\n")
    next
  }
  
  expr_subset <- pseudobulk_mat[genes_use, ]
  
  # Z-score
  expr_z <- t(scale(t(expr_subset)))
  expr_z[expr_z > 3] <- 3
  expr_z[expr_z < -3] <- -3
  
  # Match symbols to rows
  symbol_map <- setNames(top_hubs$SYMBOL, top_hubs$ENSEMBL)
  rownames(expr_z) <- symbol_map[rownames(expr_z)]
  
  sample_order <- order(patient_meta$Group)
  expr_z <- expr_z[, patient_meta$Patient[sample_order]]
  
  ha_top <- HeatmapAnnotation(
    Group = patient_meta$Group[sample_order],
    col = list(Group = cols_group),
    annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
    annotation_name_side = "left",
    show_legend = TRUE,
    height = unit(0.5, "cm")
  )
  
  # MM barplot
  mm_values <- setNames(top_hubs$MM, top_hubs$SYMBOL)[rownames(expr_z)]
  
  ha_left <- rowAnnotation(
    MM = anno_barplot(
      mm_values,
      gp = gpar(fill = mod, col = "black"),
      width = unit(2, "cm"),
      axis_param = list(gp = gpar(fontsize = 8))
    ),
    annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
    annotation_name_rot = 0
  )
  
  # Heatmap
  ht <- Heatmap(
    expr_z,
    name = "Z-score",
    
    col = colorRamp2(
      c(-2, -1, 0, 1, 2),
      c("#2166AC", "#92C5DE", "white", "#F4A582", "#B2182B")
    ),
    
    top_annotation = ha_top,
    left_annotation = ha_left,
    
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    clustering_distance_rows = "pearson",
    clustering_method_rows = "ward.D2",
    show_row_dend = TRUE,
    row_dend_width = unit(1.5, "cm"),
    
    row_names_gp = gpar(fontsize = 8, fontface = "italic"),
    column_names_gp = gpar(fontsize = 7),
    show_column_names = TRUE,
    
    border = TRUE,
    rect_gp = gpar(col = "white", lwd = 0.3),
    
    heatmap_legend_param = list(
      title = "Expression\n(Z-score)",
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 9),
      legend_direction = "vertical"
    ),
    
    width = unit(10, "cm"),
    height = unit(12, "cm")
  )
  
  draw(ht,
       column_title = paste0("Top 25 Hub Genes: ", mod, " Module"),
       column_title_gp = gpar(fontsize = 12, fontface = "bold"),
       heatmap_legend_side = "right",
       padding = unit(c(2, 2, 2, 10), "mm"))
  
  if (mod != top_6_modules[length(top_6_modules)]) {
    grid.newpage()
  }
}

dev.off()

cat("✓ Saved: 11f_WGCNA_hub_gene_heatmaps.pdf\n")

###############################################################################
### 7. Module-Trait Association Forest Plot
###############################################################################
cat("\n--- Creating forest plot ---\n")

# Get significant modules
forest_data <- module_trait_results %>%
  mutate(
    Module = gsub("^ME", "", Module),
    CI_lower = Beta - 1.96 * SE,
    CI_upper = Beta + 1.96 * SE,
    Significant = ifelse(adj.P.Val < 0.05, "Yes", "No")
  ) %>%
  filter(P.Value < 0.1) %>%  # Show nominally significant
  arrange(Beta)

forest_data$Module <- factor(forest_data$Module, levels = forest_data$Module)

p_forest <- ggplot(forest_data, aes(x = Beta, y = Module)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey50", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper, color = Significant),
                 height = 0.3, linewidth = 0.8) +
  geom_point(aes(fill = Significant), size = 4, shape = 21, color = "black", stroke = 0.5) +
  scale_color_manual(values = c("Yes" = "#D55E00", "No" = "grey50")) +
  scale_fill_manual(values = c("Yes" = "#D55E00", "No" = "grey70")) +
  labs(
    x = "Effect Size (Beta: MSA - PD)",
    y = "Module",
    title = "Module-Disease Group Associations",
    subtitle = "Error bars = 95% CI"
  ) +
  theme_classic(base_size = 11) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    panel.grid.major.x = element_line(color = "grey95", linewidth = 0.3),
    axis.text.y = element_text(face = "bold", size = 10),
    axis.text.x = element_text(size = 10),
    axis.title = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 9)
  )

ggsave(
  file.path(OUTPUT_DIR, "11g_WGCNA_forest_plot.pdf"),
  p_forest,
  width = 8,
  height = 6
)

cat("✓ Saved: 11g_WGCNA_forest_plot.pdf\n")

###############################################################################
### 8. Summary Figure
###############################################################################
cat("\n--- Creating summary figure ---\n")

p_summary <- (p_power_combined) /
  (p_sizes + p_forest) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    tag_levels = "A",
    tag_suffix = ".",
    title = "WGCNA Network Analysis Summary",
    theme = theme(
      plot.tag = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  )

ggsave(
  file.path(OUTPUT_DIR, "11h_WGCNA_summary_figure.pdf"),
  p_summary,
  width = 14,
  height = 10
)

cat("✓ Saved: 11h_WGCNA_summary_figure.pdf\n")

cat("\n=== WGCNA Visualization Complete ===\n")
cat("Files created:\n")
cat("  11a: Power selection\n")
cat("  11b: Module sizes\n")
cat("  11c: Module-trait correlations\n")
cat("  11d: Module eigengenes by group\n")
cat("  11e: Eigengene correlation network\n")
cat("  11f: Hub gene heatmaps (top 3 modules)\n")
cat("  11g: Forest plot of associations\n")
cat("  11h: Multi-panel summary figure\n")






###############################################################################
### WGCNA Follow-up Analysis
###############################################################################
cat("\n--- Loading data ---\n")

gene_modules <- read.csv(file.path(OUTPUT_DIR, "WGCNA_gene_modules_with_symbols.csv"))
trait_tests <- read.csv(file.path(OUTPUT_DIR, "WGCNA_module_trait_results.csv"))

sig_modules <- trait_tests %>%
  filter(P.Value < 0.97) %>%
  mutate(module_color = gsub("^ME", "", Module)) %>%
  pull(module_color)

cat("Significant modules (P < 0.05):", paste(sig_modules, collapse = ", "), "\n")

if (length(sig_modules) == 0) {
  sig_modules <- trait_tests %>%
    arrange(P.Value) %>%
    slice(1:16) %>%
    mutate(module_color = gsub("^ME", "", Module)) %>%
    pull(module_color)
  
  cat("Using top 6 modules:", paste(sig_modules, collapse = ", "), "\n")
}

###############################################################################
### 1. Improved Module Enrichment Plots
###############################################################################
cat("\n--- Creating enrichment plots ---\n")

cols_group <- c("PD" = "#4C72B0", "MSA" = "#C44E52")

plot_module_enrichment <- function(enrich_res, module_color,
                                               ontology = "GO",
                                               padj_cutoff = 0.25,
                                               top_n = 10) {
  
  res <- enrich_res[[ontology]]
  
  if (is.null(res) || nrow(res@result) == 0) {
    cat("  ", module_color, "-", ontology, ": No results\n")
    return(NULL)
  }
  
  # Filter significant
  df <- res@result %>%
    filter(!is.na(qvalue), qvalue < padj_cutoff)
  
  if (nrow(df) == 0) {
    cat("  ", module_color, "-", ontology, ": No significant terms\n")
    return(NULL)
  }
  
  cat("  ", module_color, "-", ontology, ":", nrow(df), "significant terms\n")
  
  # Top N
  df_top <- df %>%
    arrange(qvalue) %>%
    head(top_n) %>%
    mutate(
      Description = str_wrap(Description, width = 50),
      neg_log10_q = -log10(qvalue)
    )
  
  # Order by significance
  df_top$Description <- factor(df_top$Description,
                               levels = df_top$Description[order(df_top$neg_log10_q)])
  
  # Plot
  p <- ggplot(df_top, aes(x = neg_log10_q, y = Description)) +
    
    geom_point(aes(size = Count, fill = neg_log10_q),
               shape = 21, color = "black", alpha = 0.8, stroke = 0.3) +
    
    scale_fill_viridis_c(
      option = "C",
      direction = -1,
      name = expression(-log[10]*"(q)"),
      guide = guide_colorbar(
        barwidth = 1,
        barheight = 6,
        title.position = "top",
        title.hjust = 0.5
      )
    ) +
    
    scale_size_continuous(
      range = c(2, 8),
      name = "Gene\nCount",
      guide = guide_legend(
        title.position = "top",
        title.hjust = 0.5,
        override.aes = list(fill = "grey70", color = "black", stroke = 0.3)
      )
    ) +
    
    scale_x_continuous(
      breaks = pretty_breaks(n = 5),
      expand = expansion(mult = c(0.05, 0.1))
    ) +
    
    labs(
      x = expression(bold("-log"[10]*" q-value")),
      y = NULL,
      title = paste0(ontology, " Enrichment: ", module_color, " Module"),
      subtitle = paste0(nrow(df), " significant pathways")
    ) +
    
    # Theme
    theme_classic(base_size = 11) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
      panel.background = element_rect(fill = "white"),
      panel.grid.major.x = element_line(color = "grey95", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.text.y = element_text(size = 9, color = "black"),
      axis.text.x = element_text(size = 9, color = "black"),
      axis.title = element_text(face = "bold", size = 10),
      axis.title.x = element_text(margin = margin(t = 10)),
      
      plot.title = element_text(face = "bold", size = 12, hjust = 0),
      plot.subtitle = element_text(size = 10, color = "grey30", hjust = 0),
      
      legend.position = "right",
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8),
      
      plot.margin = margin(10, 10, 10, 10)
    )
  
  return(p)
}

###############################################################################
### 2. Generate enrichment plots
###############################################################################
for (mod in sig_modules) {
  
  cat("\n--- Module:", mod, "---\n")
  
  enrich_res <- list()
  
  # GO
  go_file <- file.path(OUTPUT_DIR, paste0("WGCNA_GO_", mod, ".csv"))
  if (file.exists(go_file)) {
    go_df <- read.csv(go_file)
    
    # Create enrichResult object
    go_res <- new("enrichResult",
                  result = go_df,
                  pvalueCutoff = 1,
                  pAdjustMethod = "BH",
                  qvalueCutoff = 1,
                  organism = "hsa",
                  ontology = "BP",
                  gene = character(0),
                  keytype = "SYMBOL")
    
    enrich_res$GO <- go_res
  }
  
  # KEGG
  kegg_file <- file.path(OUTPUT_DIR, paste0("WGCNA_KEGG_", mod, ".csv"))
  if (file.exists(kegg_file)) {
    kegg_df <- read.csv(kegg_file)
    
    kegg_res <- new("enrichResult",
                    result = kegg_df,
                    pvalueCutoff = 1,
                    pAdjustMethod = "BH",
                    qvalueCutoff = 1,
                    organism = "hsa",
                    ontology = "KEGG",
                    gene = character(0),
                    keytype = "ENTREZID")
    
    enrich_res$KEGG <- kegg_res
  }
  
  # Reactome
  react_file <- file.path(OUTPUT_DIR, paste0("WGCNA_Reactome_", mod, ".csv"))
  if (file.exists(react_file)) {
    react_df <- read.csv(react_file)
    
    react_res <- new("enrichResult",
                     result = react_df,
                     pvalueCutoff = 1,
                     pAdjustMethod = "BH",
                     qvalueCutoff = 1,
                     organism = "human",
                     ontology = "Reactome",
                     gene = character(0),
                     keytype = "ENTREZID")
    
    enrich_res$Reactome <- react_res
  }
  
  # Create plots
  plots <- list()
  
  for (ont in c("GO", "KEGG", "Reactome")) {
    p <- plot_module_enrichment(
      enrich_res,
      mod,
      ontology = ont,
      padj_cutoff = 0.25,
      top_n = 12
    )
    
    if (!is.null(p)) {
      plots[[ont]] <- p
    }
  }
  
  if (length(plots) > 0) {
    
    p_combined <- wrap_plots(plots, ncol = 1) +
      plot_annotation(
        tag_levels = "A",
        tag_suffix = ".",
        theme = theme(plot.tag = element_text(size = 14, face = "bold"))
      )
    
    ggsave(
      file.path(OUTPUT_DIR, paste0("11i_WGCNA_enrichment_", mod, "_improved.pdf")),
      p_combined,
      width = 10,
      height = 4 * length(plots),
      limitsize = FALSE
    )
    
    cat("✓ Saved: 11i_WGCNA_enrichment_", mod, "_improved.pdf\n")
  }
}

###############################################################################
### 3. Hub Gene Expression Plots
###############################################################################
cat("\n--- Creating hub gene expression plots ---\n")

plot_hub_expression <- function(module_color, n_top = 16) {
  
  cat("  Module:", module_color, "\n")
  
  hub_file <- file.path(OUTPUT_DIR, paste0("WGCNA_hub_genes_", module_color, ".csv"))
  if (!file.exists(hub_file)) {
    cat("    Hub file not found\n")
    return(NULL)
  }
  
  hub_data <- read.csv(hub_file)
  top_hubs <- head(hub_data, n_top)
  
  cat("    Top", n_top, "hub genes:", paste(top_hubs$SYMBOL, collapse = ", "), "\n")
  
  # Get expression
  assay(spe_qc, "logcounts") <- cpm(
    assay(spe_qc, "counts"),
    log         = TRUE,
    prior.count = 1
  )
  logmat <- assay(spe_qc, "logcounts")
  
  # Match genes
  genes_present <- intersect(top_hubs$ENSEMBL, rownames(logmat))
  
  if (length(genes_present) == 0) {
    cat("    No hub genes found in expression data\n")
    return(NULL)
  }
  
  # Extract expression
  expr_df <- as.data.frame(t(logmat[genes_present, , drop = FALSE]))
  expr_df$ROI <- rownames(expr_df)
  
  # Add metadata
  meta_df <- as.data.frame(colData(spe_qc))
  meta_df$ROI <- rownames(meta_df)
  
  expr_long <- expr_df %>%
    left_join(meta_df, by = "ROI") %>%
    pivot_longer(
      cols = all_of(genes_present),
      names_to = "ENSEMBL",
      values_to = "logcounts"
    ) %>%
    left_join(top_hubs %>% select(ENSEMBL, SYMBOL, MM), by = "ENSEMBL")
  
  # Order genes by MM
  expr_long <- expr_long %>%
    mutate(SYMBOL = factor(SYMBOL, levels = top_hubs$SYMBOL[order(-top_hubs$MM)]))
  
  # Summary stats
  summary_stats <- expr_long %>%
    group_by(SYMBOL, Group) %>%
    summarise(
      mean = mean(logcounts, na.rm = TRUE),
      sem = sd(logcounts, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  # Statistical tests
  pvals <- expr_long %>%
    group_by(SYMBOL) %>%
    summarise(
      p = tryCatch({
        t.test(logcounts ~ Group)$p.value
      }, error = function(e) NA),
      .groups = "drop"
    ) %>%
    mutate(
      p_label = case_when(
        is.na(p) ~ "",
        p < 0.001 ~ "***",
        p < 0.01 ~ "**",
        p < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )
  
  # Plot
  p <- ggplot(expr_long, aes(x = Group, y = logcounts, fill = Group)) +
    
    geom_violin(alpha = 0.2, color = NA, scale = "width", trim = TRUE) +
    
    geom_boxplot(width = 0.3, outlier.shape = NA, alpha = 0.7) +
    
    ggbeeswarm::geom_beeswarm(alpha = 0.6, size = 1.5, cex = 2.5) +
    
    geom_errorbar(data = summary_stats,
                  aes(x = Group, y = mean, ymin = mean - sem, ymax = mean + sem),
                  width = 0.2, linewidth = 0.7, color = "black",
                  inherit.aes = FALSE) +
    geom_point(data = summary_stats,
               aes(x = Group, y = mean),
               size = 2.5, shape = 18, color = "black",
               inherit.aes = FALSE) +
    
    # Statistical annotation
    geom_text(data = pvals,
              aes(x = 1.5, y = Inf, label = p_label),
              vjust = 1.2, size = 5, fontface = "bold",
              inherit.aes = FALSE) +
    
    facet_wrap(~SYMBOL, scales = "free_y", ncol = 3) +
    
    scale_fill_manual(values = cols_group) +
    
    labs(
      x = NULL,
      y = expression(bold("Expression (log"[2]*" CPM)")),
      title = paste0("Hub Gene Expression: ", module_color, " Module")
    ) +
    
    # Theme
    theme_classic(base_size = 11) +
    theme(
      strip.background = element_rect(fill = "grey95", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold.italic", size = 10),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
      axis.text = element_text(color = "black", size = 9),
      axis.title = element_text(face = "bold", size = 11),
      axis.text.x = element_text(face = "bold", size = 10),
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      legend.position = "none"
    )
  
  return(p)
}

# Create plots for each module
for (mod in sig_modules) {
  if (mod == "grey") next
  
  p <- plot_hub_expression(mod, n_top = 6)
  
  if (!is.null(p)) {
    ggsave(
      file.path(OUTPUT_DIR, paste0("11j_WGCNA_hub_expression_", mod, "_improved.pdf")),
      p,
      width = 10,
      height = 7
    )
    
    cat("✓ Saved: 11j_WGCNA_hub_expression_", mod, "_improved.pdf\n")
  }
}

###############################################################################
### 4. Hub-DEGs Overlap
###############################################################################
cat("\n--- Creating hub-DEG overlap plots ---\n")

hub_de_file <- file.path(OUTPUT_DIR, "WGCNA_hub_genes_DE_overlap.csv")

if (file.exists(hub_de_file)) {
  
  hub_de <- read.csv(hub_de_file) %>%
    filter(!is.na(adj.P.Val))
  
  # Scatterplot MM vs logFC
  p_overlap <- ggplot(hub_de, aes(x = MM, y = logFC)) +
    
    geom_hline(yintercept = 0, linetype = "solid", color = "grey50", linewidth = 0.5) +
    geom_point(aes(color = in_DE_sig, size = -log10(adj.P.Val)),
               alpha = 0.7, shape = 16) +
    
    scale_color_manual(
      values = c("TRUE" = "#D55E00", "FALSE" = "grey70"),
      labels = c("TRUE" = "DE (FDR < 0.05)", "FALSE" = "Not DE"),
      name = NULL
    ) +
    
    scale_size_continuous(
      range = c(1, 6),
      name = expression(-log[10]*" FDR")
    ) +
    
    facet_wrap(~module, ncol = 2) +
    
    labs(
      x = "Module Membership (MM)",
      y = expression(bold("Effect Size (log"[2]*" FC)")),
      title = "Hub Gene Overlap with Differential Expression"
    ) +
    
    # Theme
    theme_classic(base_size = 11) +
    theme(
      strip.background = element_rect(fill = "grey95", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 10),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
      panel.grid.major = element_line(color = "grey95", linewidth = 0.3),
      axis.text = element_text(color = "black", size = 10),
      axis.title = element_text(face = "bold", size = 11),
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      legend.position = "bottom",
      legend.box = "vertical"
    )
  
  ggsave(
    file.path(OUTPUT_DIR, "11k_WGCNA_hub_DE_overlap.pdf"),
    p_overlap,
    width = 10,
    height = 15
  )
  
  cat("✓ Saved: 11k_WGCNA_hub_DE_overlap.pdf\n")
  
  # Summary barplot
  overlap_summary <- hub_de %>%
    group_by(module) %>%
    summarise(
      Total = n(),
      DE_sig = sum(in_DE_sig, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(cols = c(Total, DE_sig),
                 names_to = "Category",
                 values_to = "Count")
  
  p_summary <- ggplot(overlap_summary, aes(x = module, y = Count, fill = Category)) +
    geom_col(position = "dodge", color = "black", linewidth = 0.3) +
    geom_text(aes(label = Count), position = position_dodge(width = 0.9),
              vjust = -0.5, size = 3, fontface = "bold") +
    scale_fill_manual(
      values = c("Total" = "grey70", "DE_sig" = "#D55E00"),
      labels = c("Total" = "Hub Genes", "DE_sig" = "Also DE (FDR < 0.05)"),
      name = NULL
    ) +
    labs(
      x = "Module",
      y = "Number of Genes",
      title = "Hub Gene Overlap with Differential Expression"
    ) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x = element_text(face = "bold", size = 10),
      axis.title = element_text(face = "bold", size = 11),
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      legend.position = "bottom",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )
  
  ggsave(
    file.path(OUTPUT_DIR, "11l_WGCNA_hub_DE_summary.pdf"),
    p_summary,
    width = 12,
    height = 6
  )
  
  cat("✓ Saved: 11l_WGCNA_hub_DE_summary.pdf\n")
}

###############################################################################
### 5. Known PD Genes
###############################################################################

cat("\n--- Creating PD genes plot ---\n")

pd_genes_file <- file.path(OUTPUT_DIR, "WGCNA_PD_genes_modules.csv")

if (file.exists(pd_genes_file)) {
  
  # Load ALL genes
  pd_genes_df <- read.csv(pd_genes_file)
  
  cat("Total PD genes in CSV:", nrow(pd_genes_df), "\n")
  cat("Modules represented:", paste(unique(pd_genes_df$module), collapse = ", "), "\n")
  
  #Load all hub files and bind into one MM lookup table
  all_modules <- unique(pd_genes_df$module)
  
  mm_lookup <- lapply(all_modules, function(mod) {
    
    hub_file <- file.path(OUTPUT_DIR, paste0("WGCNA_hub_genes_", mod, ".csv"))
    
    if (!file.exists(hub_file)) {
      cat("  No hub file for module:", mod, "\n")
      return(NULL)
    }
    
    hub_data <- read.csv(hub_file)
    
    if (!all(c("SYMBOL", "MM") %in% colnames(hub_data))) {
      cat("  Hub file missing SYMBOL or MM column for module:", mod, "\n")
      return(NULL)
    }
    
    hub_data %>%
      dplyr::select(SYMBOL, MM) %>%
      mutate(module = mod)
    
  }) %>%
    bind_rows()
  
  cat("MM lookup table rows:", nrow(mm_lookup), "\n")
  
  pd_with_mm <- pd_genes_df %>%
    # Standardise case to avoid mismatch
    mutate(SYMBOL = trimws(SYMBOL),
           module = trimws(module)) %>%
    left_join(
      mm_lookup %>% mutate(SYMBOL = trimws(SYMBOL),
                           module = trimws(module)),
      by = c("SYMBOL", "module")
    ) %>%
    # Replace NA MM with median so all genes still appear
    mutate(MM = replace_na(MM, median(mm_lookup$MM, na.rm = TRUE)))
  
  cat("Genes after join:", nrow(pd_with_mm), "\n")
  cat("Genes with MM value:", sum(!is.na(pd_with_mm$MM)), "\n")
  
  # Exclude grey
  grey_genes <- pd_with_mm %>% filter(module == "grey")
  if (nrow(grey_genes) > 0) {
    cat("Grey module genes (unassigned, excluded from plot):",
        paste(grey_genes$SYMBOL, collapse = ", "), "\n")
  }
  
  pd_with_mm <- pd_with_mm %>% filter(module != "grey")
  
  cat("Genes in final plot:", nrow(pd_with_mm), "\n")
  
  if (nrow(pd_with_mm) > 0) {
    
    # Order genes by module then MM for cleaner layout
    pd_with_mm <- pd_with_mm %>%
      arrange(module, desc(MM)) %>%
      mutate(SYMBOL = factor(SYMBOL, levels = unique(SYMBOL)))
    
    # Plot
    p_pd <- ggplot(pd_with_mm, aes(x = module, y = SYMBOL)) +
      
      geom_point(
        aes(size = MM, fill = module),
        shape  = 21,
        colour = "black",
        stroke = 0.5,
        alpha  = 0.85
      ) +
      
      scale_fill_identity(guide = "none") +
      
      scale_size_continuous(
        range  = c(3, 10),
        name   = "Module\nMembership",
        breaks = c(0.3, 0.5, 0.7, 0.9)
      ) +
      
      geom_hline(
        yintercept = seq_along(levels(pd_with_mm$SYMBOL)),
        colour     = "grey93",
        linewidth  = 0.3
      ) +
      
      labs(
        x       = "WGCNA Module",
        y       = NULL,
        title   = "Known PD/Neurodegeneration Genes in Network Modules",
        caption = paste0(
          "Point size = module membership (MM). ",
          nrow(grey_genes), " gene(s) in grey (unassigned) module not shown."
        )
      ) +
      
      theme_classic(base_size = 11) +
      theme(
        axis.text.x      = element_text(face = "bold", angle = 45,
                                        hjust = 1, size = 10),
        axis.text.y      = element_text(face = "bold.italic", size = 10),
        axis.title       = element_text(face = "bold", size = 11),
        plot.title       = element_text(face = "bold", size = 12, hjust = 0.5),
        plot.caption     = element_text(size = 8, colour = "grey40", hjust = 0),
        panel.border     = element_rect(colour = "black", fill = NA,
                                        linewidth = 0.5),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
        panel.grid.major.y = element_blank(),   # replaced by geom_hline above
        legend.position  = "right"
      )
    
    # Dynamic height so all gene labels fit 
    plot_height <- max(5, nrow(pd_with_mm) * 0.35 + 2)
    
    ggsave(
      file.path(OUTPUT_DIR, "11m_WGCNA_PD_genes.pdf"),
      p_pd,
      width  = 4,
      height = plot_height,
      limitsize = FALSE
    )
    
    ggsave(
      file.path(OUTPUT_DIR, "11m_WGCNA_PD_genes.png"),
      p_pd,
      width     = 4,
      height    = plot_height,
      dpi       = 600,
      bg        = "white",
      limitsize = FALSE
    )
    
    cat("✓ Saved: 11m_WGCNA_PD_genes.pdf\n")
    cat("  Plot dimensions: 8 x", round(plot_height, 1), "inches\n")
  }
}
cat("\n=== WGCNA Follow-up Visualization Complete ===\n")
cat("Files created:\n")
cat("  11i: Module enrichment plots (GO/KEGG/Reactome per module)\n")
cat("  11j: Hub gene expression plots (per module)\n")
cat("  11k: Hub-DE overlap scatterplot\n")
cat("  11l: Hub-DE overlap summary barplot\n")
cat("  11m: Known PD genes in modules\n")











###############################################################################
### Improved Module-Clinical Trait Heatmap
###############################################################################

library(ComplexHeatmap)
library(circlize)
library(grid)
library(dplyr)
library(readxl)

cat("\n=== Creating Module-Clinical Trait Heatmap ===\n")

###############################################################################
### 1. Load data
###############################################################################
cat("\n--- Loading clinical data ---\n")

# Helper functions
strip_dcc <- function(x) sub("\\.dcc$", "", as.character(x), ignore.case = TRUE)
canon <- function(x) { x <- trimws(as.character(x)); toupper(x) }

trait_df <- readxl::read_excel(TRAIT_FILE) %>% as.data.frame()

cat("Trait file columns:", paste(colnames(trait_df), collapse = ", "), "\n")

# Clean ROI names
trait_df$ROI <- trimws(strip_dcc(trait_df$ROI))
rownames(trait_df) <- trait_df$ROI

# Align with MEs_roi
me_rois_raw <- rownames(MEs_roi)
me_rois_clean <- strip_dcc(me_rois_raw)

# Check alignment
missing <- setdiff(me_rois_clean, rownames(trait_df))
if (length(missing) > 0) {
  cat("WARNING:", length(missing), "ROIs missing from trait file\n")
}

trait_df <- trait_df[me_rois_clean, , drop = FALSE]
rownames(trait_df) <- me_rois_raw

###############################################################################
### 2. Create clean numeric trait matrix
###############################################################################
cat("\n--- Creating numeric trait matrix ---\n")

id_cols <- c("ROI", "Patient")
keep_binary <- c("Group", "Sex")
keep_numeric <- c("Age", "PMI", "Duration", "Tau", "AB")

keep_binary <- intersect(keep_binary, colnames(trait_df))
keep_numeric <- intersect(keep_numeric, colnames(trait_df))

cat("Binary traits:", paste(keep_binary, collapse = ", "), "\n")
cat("Numeric traits:", paste(keep_numeric, collapse = ", "), "\n")

# Create numeric matrix
trait_numeric <- data.frame(row.names = rownames(trait_df))

# Binary encodings
if ("Group" %in% keep_binary) {
  g <- trimws(toupper(as.character(trait_df$Group)))
  trait_numeric$Group_MSA <- as.numeric(factor(g, levels = c("PD", "MSA"))) - 1
}

if ("Sex" %in% keep_binary) {
  s <- trimws(toupper(as.character(trait_df$Sex)))
  s[s %in% c("MALE", "M")] <- "M"
  s[s %in% c("FEMALE", "F")] <- "F"
  
  ux <- sort(unique(na.omit(s)))
  if (length(ux) == 2) {
    trait_numeric$Sex_Female <- ifelse(s == "F", 1, ifelse(s == "M", 0, NA))
  }
}

# Numeric traits
numify <- function(x) suppressWarnings(as.numeric(x))

for (nm in keep_numeric) {
  v <- numify(trait_df[[nm]])
  if (sd(v, na.rm = TRUE) > 0) {
    trait_numeric[[nm]] <- v
  } else {
    cat("  Removing", nm, "(constant)\n")
  }
}

# Remove all-NA columns
trait_numeric <- trait_numeric[, colSums(!is.na(trait_numeric)) > 0, drop = FALSE]

cat("Final trait matrix:", nrow(trait_numeric), "ROIs ×", ncol(trait_numeric), "traits\n")
cat("Traits:", paste(colnames(trait_numeric), collapse = ", "), "\n")

###############################################################################
### 3. Calculate correlations
###############################################################################
cat("\n--- Calculating correlations ---\n")

# Order module eigengenes
MEs_roi_ord <- WGCNA::orderMEs(MEs_roi)

# Correlations
module_trait_cor_roi <- cor(MEs_roi_ord, trait_numeric, use = "pairwise.complete.obs")
module_trait_p_roi <- WGCNA::corPvalueStudent(module_trait_cor_roi, nrow(MEs_roi_ord))

# Remove ME prefix
rownames(module_trait_cor_roi) <- gsub("^ME", "", rownames(module_trait_cor_roi))
rownames(module_trait_p_roi) <- gsub("^ME", "", rownames(module_trait_p_roi))

cat("Correlation matrix:", nrow(module_trait_cor_roi), "modules ×", 
    ncol(module_trait_cor_roi), "traits\n")


cor_plot <- module_trait_cor_roi
p_plot <- module_trait_p_roi

###############################################################################
### 4. Create heatmap
###############################################################################
cat("\n--- Creating heatmap ---\n")

cor_plot <- as.matrix(cor_plot)
p_plot <- as.matrix(p_plot)

# Create text matrices
cor_text <- matrix(
  sprintf("%.2f", cor_plot),
  nrow = nrow(cor_plot),
  ncol = ncol(cor_plot),
  dimnames = dimnames(cor_plot)
)

sig_text <- matrix(
  ifelse(p_plot < 0.001, "***",
         ifelse(p_plot < 0.01, "**",
                ifelse(p_plot < 0.05, "*", ""))),
  nrow = nrow(p_plot),
  ncol = ncol(p_plot),
  dimnames = dimnames(p_plot)
)

# Combined text (correlation + significance)
cell_text <- matrix(
  paste0(cor_text, "\n", sig_text),
  nrow = nrow(cor_plot),
  ncol = ncol(cor_plot),
  dimnames = dimnames(cor_plot)
)

colnames_clean <- colnames(cor_plot)
colnames_clean <- gsub("_", " ", colnames_clean)
colnames_clean <- gsub("Group MSA", "Disease (MSA)", colnames_clean)
colnames_clean <- gsub("Sex Female", "Sex (Female)", colnames_clean)

max_abs <- max(abs(cor_plot), na.rm = TRUE)
max_abs <- min(max_abs, 1)  # Cap at 1

col_fun <- colorRamp2(
  seq(-max_abs, max_abs, length.out = 9),
  c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
    "white",
    "#FDDBC7", "#F4A582", "#D6604D", "#B2182B")
)

module_names <- rownames(cor_plot)

module_col_anno <- rowAnnotation(
  Module = module_names,
  col = list(Module = setNames(module_names, module_names)),
  show_annotation_name = FALSE,
  show_legend = FALSE,
  width = unit(0.5, "cm")
)

# Create heatmap
ht <- Heatmap(
  cor_plot,
  name = "Correlation",
  
  col = col_fun,
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    # Only show text for cells with p < 0.1 to reduce clutter
    if (p_plot[i, j] < 0.1) {
      grid.text(
        cell_text[i, j],
        x, y,
        gp = gpar(
          fontsize = ifelse(p_plot[i, j] < 0.05, 9, 7),
          fontface = ifelse(p_plot[i, j] < 0.05, "bold", "plain"),
          col = ifelse(abs(cor_plot[i, j]) > 0.4, "white", "black")
        )
      )
    }
  },
  
  # Clustering
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "ward.D2",
  show_row_dend = TRUE,
  row_dend_width = unit(2, "cm"),
  
  left_annotation = module_col_anno,
  
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_rot = 45,
  column_labels = colnames_clean,
  
  border = TRUE,
  rect_gp = gpar(col = "white", lwd = 1),
  
  row_title_gp = gpar(fontsize = 11, fontface = "bold"),
  
  width = unit(ncol(cor_plot) * 1.5, "cm"),
  height = unit(nrow(cor_plot) * 0.8, "cm"),
  
  heatmap_legend_param = list(
    title = "Pearson r",
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 9),
    grid_height = unit(4, "mm"),
    grid_width = unit(4, "mm"),
    legend_direction = "vertical",
    title_position = "topcenter",
    at = seq(-max_abs, max_abs, length.out = 5),
    labels = sprintf("%.2f", seq(-max_abs, max_abs, length.out = 5))
  )
)

fig_width <- 4 + ncol(cor_plot) * 1.5
fig_height <- 5 + nrow(cor_plot) * 0.8

pdf(
  file.path(OUTPUT_DIR, "11p_WGCNA_module_clinical_traits.pdf"),
  width = min(fig_width, 16),
  height = min(fig_height, 20)
)

draw(ht,
     column_title = "Module-Clinical Trait Associations",
     column_title_gp = gpar(fontsize = 14, fontface = "bold"),
     heatmap_legend_side = "right",
     padding = unit(c(2, 2, 2, 12), "mm"))

# Add significance legend
grid.text(
  "* P < 0.05\n** P < 0.01\n*** P < 0.001",
  x = 0.98,
  y = 0.15,
  just = c("right", "top"),
  gp = gpar(fontsize = 9, fontface = "bold")
)

dev.off()

cat("✓ Saved: 11p_WGCNA_module_clinical_traits.pdf\n")

###############################################################################
### 5. Create heatmap with only significant associations
###############################################################################
cat("\n--- Creating heatmap with significant associations only ---\n")

# Identify modules and traits with at least one significant correlation
sig_modules <- rownames(cor_plot)[apply(p_plot < 0.05, 1, any)]
sig_traits <- colnames(cor_plot)[apply(p_plot < 0.05, 2, any)]

if (length(sig_modules) > 0 && length(sig_traits) > 0) {
  
  cat("Significant associations:\n")
  cat("  Modules:", length(sig_modules), "-", paste(sig_modules, collapse = ", "), "\n")
  cat("  Traits:", length(sig_traits), "-", paste(sig_traits, collapse = ", "), "\n")
  
  # Subset
  cor_sig <- cor_plot[sig_modules, sig_traits, drop = FALSE]
  p_sig <- p_plot[sig_modules, sig_traits, drop = FALSE]
  
  cor_text_sig <- sprintf("%.2f", cor_sig)
  sig_text_sig <- ifelse(p_sig < 0.001, "***",
                         ifelse(p_sig < 0.01, "**",
                                ifelse(p_sig < 0.05, "*", "")))
  cell_text_sig <- paste0(cor_text_sig, "\n", sig_text_sig)
  dim(cell_text_sig) <- dim(cor_sig)
  
  module_col_anno_sig <- rowAnnotation(
    Module = sig_modules,
    col = list(Module = setNames(sig_modules, sig_modules)),
    show_annotation_name = FALSE,
    show_legend = FALSE,
    width = unit(0.5, "cm")
  )
  
  # Clean names
  colnames_clean_sig <- gsub("_", " ", colnames(cor_sig))
  colnames_clean_sig <- gsub("Group MSA", "Disease (MSA)", colnames_clean_sig)
  colnames_clean_sig <- gsub("Sex Female", "Sex (Female)", colnames_clean_sig)
  
  # Heatmap
  ht_sig <- Heatmap(
    cor_sig,
    name = "Correlation",
    
    col = col_fun,
    
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(
        cell_text_sig[i, j],
        x, y,
        gp = gpar(
          fontsize = 10,
          fontface = "bold",
          col = ifelse(abs(cor_sig[i, j]) > 0.4, "white", "black")
        )
      )
    },
    
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "ward.D2",
    show_row_dend = TRUE,
    row_dend_width = unit(2, "cm"),
    
    left_annotation = module_col_anno_sig,
    
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 11, fontface = "bold"),
    column_names_gp = gpar(fontsize = 11, fontface = "bold"),
    column_names_rot = 45,
    column_labels = colnames_clean_sig,
    
    border = TRUE,
    rect_gp = gpar(col = "white", lwd = 1),
    
    width = unit(ncol(cor_sig) * 2, "cm"),
    height = unit(nrow(cor_sig) * 1, "cm"),
    
    heatmap_legend_param = list(
      title = "Pearson r",
      title_gp = gpar(fontsize = 11, fontface = "bold"),
      labels_gp = gpar(fontsize = 10),
      grid_height = unit(5, "mm"),
      grid_width = unit(5, "mm"),
      legend_direction = "vertical",
      title_position = "topcenter",
      at = seq(-max_abs, max_abs, length.out = 5),
      labels = sprintf("%.2f", seq(-max_abs, max_abs, length.out = 5))
    )
  )
  
  pdf(
    file.path(OUTPUT_DIR, "11q_WGCNA_module_clinical_traits_significant_only.pdf"),
    width = 4 + ncol(cor_sig) * 2,
    height = 5 + nrow(cor_sig) * 1
  )
  
  draw(ht_sig,
       column_title = "Significant Module-Clinical Trait Associations",
       column_title_gp = gpar(fontsize = 14, fontface = "bold"),
       heatmap_legend_side = "right",
       padding = unit(c(2, 2, 2, 12), "mm"))
  
  grid.text(
    "* P < 0.05 | ** P < 0.01 | *** P < 0.001",
    x = 0.98,
    y = 0.1,
    just = c("right", "top"),
    gp = gpar(fontsize = 10, fontface = "bold")
  )
  
  dev.off()
  
  cat("✓ Saved: 11q_WGCNA_module_clinical_traits_significant_only.pdf\n")
  
} else {
  cat("No significant associations found\n")
}

###############################################################################
### 6. Export correlation table
###############################################################################
cat("\n--- Exporting correlation table ---\n")

# Create long-format table
cor_table <- cor_plot %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Module") %>%
  tidyr::pivot_longer(
    cols = -Module,
    names_to = "Trait",
    values_to = "Correlation"
  )

p_table <- p_plot %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Module") %>%
  tidyr::pivot_longer(
    cols = -Module,
    names_to = "Trait",
    values_to = "P.Value"
  )

cor_table <- cor_table %>%
  left_join(p_table, by = c("Module", "Trait")) %>%
  mutate(
    Significant = P.Value < 0.05,
    FDR = p.adjust(P.Value, method = "BH")
  ) %>%
  arrange(P.Value)

write.csv(
  cor_table,
  file.path(OUTPUT_DIR, "WGCNA_module_clinical_correlations.csv"),
  row.names = FALSE
)

cat("✓ Saved: WGCNA_module_clinical_correlations.csv\n")

# Summary
cat("\nSummary:\n")
cat("  Total correlations tested:", nrow(cor_table), "\n")
cat("  Significant (P < 0.05):", sum(cor_table$P.Value < 0.05), "\n")
cat("  FDR-significant (FDR < 0.05):", sum(cor_table$FDR < 0.05), "\n")

# Top associations
cat("\nTop 10 associations:\n")
print(head(cor_table %>% select(Module, Trait, Correlation, P.Value), 10))

cat("\n=== Module-Clinical Trait Visualization Complete ===\n")
cat("Files created:\n")
cat("  11p: Full heatmap (all modules × all traits)\n")
cat("  11q: Focused heatmap (only significant associations)\n")
cat("  WGCNA_module_clinical_correlations.csv: Complete correlation table\n")










###############################################################################
### TOM Network Visualization
###############################################################################
library(WGCNA)

# Ensure correct orientation: samples as rows, genes as columns, expr_use should be samples × genes
cat("expr_use dimensions:", nrow(expr_use), "samples x", ncol(expr_use), "genes\n")

# Run blockwiseModules with TOM saving
net <- blockwiseModules(
  expr_use,
  
  power             = 12,   # from soft threshold analysis
  networkType       = "signed",
  TOMType           = "signed",

  minModuleSize     = 30,
  mergeCutHeight    = 0.25,
  deepSplit         = 2,
  
  # Save TOM to disk
  saveTOMs          = TRUE,
  saveTOMFileBase   = file.path(OUTPUT_DIR, "WGCNA_TOM"),
  
  # Reproducibility
  numericLabels     = TRUE,
  randomSeed        = 42,
  
  # Parallelisation
  nThreads          = 4,
  
  verbose           = 3
)

cat("Modules found:", length(unique(net$colors)), "\n")
cat("Module sizes:\n")
print(table(net$colors))

saveRDS(net, file.path(OUTPUT_DIR, "WGCNA_net.rds"))
cat("Saved: WGCNA_net.rds\n")

tom_file <- file.path(OUTPUT_DIR, "WGCNA_TOM-block.1.RData")

cat("Loading TOM from:", tom_file, "\n")
cat("File exists:", file.exists(tom_file), "\n")

# Load into a temporary environment
tmp_env <- new.env()
load(tom_file, envir = tmp_env)

cat("Objects in TOM file:", paste(ls(tmp_env), collapse = ", "), "\n")

TOM <- as.matrix(tmp_env$TOM)

# Set gene names as row and column names
rownames(TOM) <- colnames(TOM) <- colnames(expr_use)

cat("TOM loaded successfully\n")
cat("TOM dimensions:", nrow(TOM), "x", ncol(TOM), "\n")

# Clean up temp environment
rm(tmp_env)
gc()

# ── Compute dissTOM and gene tree ─────────────────────────────────────────────
dissTOM  <- 1 - TOM
geneTree <- hclust(as.dist(dissTOM), method = "average")

# ── Module colours ────────────────────────────────────────────────────────────
moduleColors <- labels2colors(net$colors)

cat("dissTOM computed\n")
cat("Module colours:", paste(unique(moduleColors), collapse = ", "), "\n")

cat("\n=== Object check ===\n")
cat("expr_use:     ", nrow(expr_use),    "x", ncol(expr_use),    "\n")
cat("TOM:          ", nrow(TOM),         "x", ncol(TOM),         "\n")
cat("dissTOM:      ", nrow(dissTOM),     "x", ncol(dissTOM),     "\n")
cat("geneTree:     ", length(geneTree$order), "genes\n")
cat("moduleColors: ", length(moduleColors),   "genes\n")
cat("Modules:      ", length(unique(moduleColors)), "\n")


expr_use <- roi_expr_common

cat("Expression matrix:  ", nrow(expr_use), "samples ×", ncol(expr_use), "genes\n")
cat("TOM:                ", nrow(TOM),      "×",         ncol(TOM),      "\n")
cat("module_colors:      ", length(module_colors), "genes\n")

# ── Align all three objects to the same genes ─────────────────────────────────
common_genes <- rownames(TOM) %>%
  intersect(names(module_colors)) %>%
  intersect(colnames(expr_use))

cat("\nCommon genes across all objects:", length(common_genes), "\n")

if (length(common_genes) == 0) {
  stop("ERROR: No common genes — check your objects")
}

TOM           <- TOM[common_genes, common_genes, drop = FALSE]
module_colors <- module_colors[common_genes]
expr_use      <- expr_use[, common_genes, drop = FALSE]

cat("✓ TOM:           ", nrow(TOM),            "×", ncol(TOM),          "\n")
cat("✓ module_colors: ", length(module_colors),                           "\n")
cat("✓ expr_use:      ", nrow(expr_use),        "×", ncol(expr_use),     "\n")
cat("✓ Modules:       ", paste(unique(module_colors), collapse = ", "),   "\n")

###############################################################################
### Get significant modules to plot
###############################################################################
cat("\n--- Getting modules to plot ---\n")

if (!exists("sig_modules")) {
  cat("sig_modules not found — using all non-grey modules\n")
  sig_modules <- setdiff(unique(module_colors), "grey")
}

# NO grey
sig_modules <- sig_modules[sig_modules != "grey"]

cat("Modules to plot:", paste(sig_modules, collapse = ", "), "\n")
cat("Total modules:  ", length(sig_modules), "\n")

###############################################################################
### Get DEGs for highlighting
###############################################################################
cat("\n--- Part 3: Getting DEGs ---\n")

deg_genes <- character(0)

if (exists("de_neuron")) {
  
  de_df <- de_neuron %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene_id")
  
  cat("DE results: ", nrow(de_df), "genes\n")
  cat("DE ID type: ", de_df$gene_id[1], "\n")
  
  # Convert SYMBOL to ENSEMBL if needed
  if (!grepl("^ENSG", de_df$gene_id[1])) {
    cat("Converting DEG IDs: SYMBOL → ENSEMBL\n")
    
    de_ensembl <- mapIds(
      org.Hs.eg.db,
      keys      = de_df$gene_id,
      column    = "ENSEMBL",
      keytype   = "SYMBOL",
      multiVals = "first"
    )
    de_df$gene_id <- ifelse(!is.na(de_ensembl), de_ensembl, NA)
  }
  
  deg_genes <- de_df %>%
    filter(!is.na(gene_id))    %>%
    filter(adj.P.Val < 0.05)   %>%
    filter(abs(logFC)  > 0.5)  %>%
    pull(gene_id)               %>%
    unique()
  
  cat("✓ DEGs for highlighting:", length(deg_genes), "\n")
  
} else {
  cat("No DE results found — nodes will not be highlighted\n")
}

###############################################################################
### Network plotting function
###############################################################################

plot_tom_network <- function(module,
                             TOM,
                             moduleColors,
                             expr_use,
                             deg_genes   = character(0),
                             top_n_edges = 500,
                             min_tom     = 0.08) {
  
  cat("\n--- Module:", module, "---\n")
  
  # ── Get module genes ──────────────────────────────────────────────
  mod_genes <- names(moduleColors)[moduleColors == module] %>%
    intersect(rownames(TOM))       %>%
    intersect(colnames(expr_use))
  
  cat("  Genes:", length(mod_genes), "\n")
  
  if (length(mod_genes) < 10) {
    cat("  SKIP: fewer than 10 genes\n")
    return(NULL)
  }
  
  # ── Subset TOM ────────────────────────────────────────────────────
  TOM_mod       <- TOM[mod_genes, mod_genes, drop = FALSE]
  diag(TOM_mod) <- 0
  
  cat("  TOM range:", round(min(TOM_mod), 3), 
      "to",          round(max(TOM_mod), 3), "\n")
  
  # ── Build edges ───────────────────────────────────────────────────
  ut        <- which(upper.tri(TOM_mod), arr.ind = TRUE)
  edges_all <- tibble(
    from   = rownames(TOM_mod)[ut[, 1]],
    to     = colnames(TOM_mod)[ut[, 2]],
    weight = TOM_mod[ut]
  ) %>% arrange(desc(weight))
  
  # Try min_tom threshold first
  edges_use <- edges_all %>% filter(weight > min_tom)
  
  # If too sparse, drop to top 5% of weights
  if (nrow(edges_use) < 10) {
    adaptive  <- quantile(edges_all$weight, 0.95)
    cat("  Only", nrow(edges_use), "edges at TOM >", min_tom,
        "— using adaptive threshold:", round(adaptive, 4), "\n")
    edges_use <- edges_all %>% filter(weight > adaptive)
  }
  
  if (nrow(edges_use) < 2) {
    cat("  SKIP: not enough edges\n")
    return(NULL)
  }
  
  edges <- head(edges_use, top_n_edges)
  cat("  Edges:", nrow(edges), "\n")
  
  # ── Compute eigengene + kME ──────────────────────────────────────
  MEs      <- moduleEigengenes(expr_use, moduleColors)$eigengenes
  eig_name <- paste0("ME", module)
  
  if (!eig_name %in% colnames(MEs)) {
    cat("  SKIP: eigengene", eig_name, "not found\n")
    return(NULL)
  }
  
  kME        <- cor(
    expr_use[, mod_genes, drop = FALSE],
    MEs[, eig_name,       drop = FALSE],
    use = "pairwise.complete.obs"
  )
  kME_vec        <- as.numeric(kME)
  names(kME_vec) <- mod_genes
  
  cat("  kME range:", round(min(kME_vec, na.rm = TRUE), 3),
      "to",          round(max(kME_vec, na.rm = TRUE), 3), "\n")
  
  # ── Convert ENSEMBL → SYMBOL for labels ───────────────────────────
  gene_labels <- mapIds(
    org.Hs.eg.db,
    keys      = mod_genes,
    column    = "SYMBOL",
    keytype   = "ENSEMBL",
    multiVals = "first"
  )
  gene_labels[is.na(gene_labels)] <- names(gene_labels)[is.na(gene_labels)]
  
  # ── Build node table ──────────────────────────────────────────────
  keep_nodes <- unique(c(edges$from, edges$to)) %>%
    intersect(names(kME_vec))
  
  if (length(keep_nodes) < 2) {
    cat("  SKIP: not enough nodes\n")
    return(NULL)
  }
  
  hub_cutoff <- quantile(abs(kME_vec[keep_nodes]), 0.90, na.rm = TRUE)
  
  nodes <- tibble(
    name      = keep_nodes,
    label     = gene_labels[keep_nodes],
    node_size = rescale(abs(kME_vec[keep_nodes]), to = c(2, 10)),
    hub       = abs(kME_vec[keep_nodes]) > hub_cutoff,
    DEG       = keep_nodes %in% deg_genes
  )
  
  cat("  Nodes:", nrow(nodes),
      "| Hubs:", sum(nodes$hub),
      "| DEGs:", sum(nodes$DEG), "\n")
  
  # ── Step 7: Build graph ───────────────────────────────────────────────────
  g     <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
  g_tbl <- as_tbl_graph(g)
  
  # ── Step 8: Layout ────────────────────────────────────────────────────────
  set.seed(42)
  lay <- create_layout(g_tbl, layout = "fr")
  
  # ── Step 9: Plot ──────────────────────────────────────────────────────────
  p <- ggraph(lay) +
    
    # Edges
    geom_edge_link(
      aes(width = weight, alpha = weight),
      color = "grey50",
      show.legend = FALSE
    ) +
    scale_edge_width_continuous(range = c(0.2, 1.5)) +
    scale_edge_alpha_continuous(range = c(0.2, 0.8)) +
    
    # Regular nodes (not hub, not DEG)
    geom_node_point(
      data  = function(x) subset(x, !DEG & !hub),
      aes(size = node_size),
      fill  = "grey85",
      color = "grey40",
      shape = 21, stroke = 0.3, alpha = 0.8
    ) +
    
    # Hub nodes
    geom_node_point(
      data  = function(x) subset(x, !DEG & hub),
      aes(size = node_size),
      fill  = module,
      color = "black",
      shape = 21, stroke = 0.5, alpha = 0.9
    ) +
    
    # DEG nodes — only added if DEGs exist in this module
    {if (sum(nodes$DEG) > 0) list(
      
      # DEG, not hub
      geom_node_point(
        data  = function(x) subset(x, DEG & !hub),
        aes(size = node_size),
        fill  = "#D55E00",
        color = "black",
        shape = 21, stroke = 0.5, alpha = 0.9
      ),
      
      # DEG + hub (diamond shape)
      geom_node_point(
        data  = function(x) subset(x, DEG & hub),
        aes(size = node_size),
        fill  = "#D55E00",
        color = "black",
        shape = 23, stroke = 0.8, alpha = 1.0
      )
    )} +
    
    scale_size_identity() +
    
    # Labels on hub genes only
    geom_node_text(
      aes(label = ifelse(hub, label, "")),
      size          = 2.5,
      fontface      = "italic",
      repel         = TRUE,
      max.overlaps  = Inf,
      box.padding   = 0.5,
      point.padding = 0.3,
      segment.size  = 0.2,
      segment.color = "grey30"
    ) +
    
    labs(
      title    = paste0(module, " Module Network"),
      subtitle = paste0(
        nrow(nodes),    " genes  |  ",
        nrow(edges),    " edges  |  ",
        sum(nodes$hub), " hubs"
      )
    ) +
    
    theme_void() +
    theme(
      plot.title      = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle   = element_text(size = 10,               hjust = 0.5),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  return(p)
}

###############################################################################
###Generate and save 
###############################################################################

cat("\n--- Generating and saving plots---\n")

if (!exists("OUTPUT_DIR")) {
  OUTPUT_DIR <- getwd()
  cat("OUTPUT_DIR not set — saving to working directory:", OUTPUT_DIR, "\n")
}

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

plots_saved  <- 0
plots_failed <- 0

for (mod in sig_modules) {
  
  # Generate the plot
  p <- tryCatch({
    plot_tom_network(
      module       = mod,
      TOM          = TOM,
      moduleColors = module_colors,
      expr_use     = expr_use,
      deg_genes    = deg_genes,
      top_n_edges  = 500,
      min_tom      = 0.05
    )
  }, error = function(e) {
    cat("  ERROR in module", mod, ":", e$message, "\n")
    return(NULL)
  })
  
  if (!is.null(p)) {
    
    out_file <- file.path(
      OUTPUT_DIR,
      paste0("WGCNA_network_", mod, ".pdf")
    )
    
    ggsave(
      filename = out_file,
      plot     = p,
      width    = 10,
      height   = 10,
      device   = "pdf"
    )
    
    cat("  ✓ Saved:", out_file, "\n")
    plots_saved <- plots_saved + 1
    
  } else {
    plots_failed <- plots_failed + 1
  }
}

cat("\n=== Done ===\n")
cat("✓ Plots saved: ", plots_saved,  "\n")
cat("✗ Plots failed:", plots_failed, "\n")
cat("  Output dir:  ", OUTPUT_DIR,   "\n")
















