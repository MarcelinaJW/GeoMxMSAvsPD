###############################################################################
### WGCNA####################################################################
###############################################################################
#purpose= find groups of genes (modules) that behave similarly across samples + relate those modules to biological traits (MSA vs PD, PMI, Duration, Sex, Age).

suppressPackageStartupMessages({
  library(WGCNA)
  library(edgeR)
  library(ggplot2)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)
  library(readxl)
  library(grid)
})

cat("\n=== WGCNA Analysis (Neurons - Pseudobulk) ===\n")
allowWGCNAThreads()

set.seed(12345)
###############################################################################
### Create Pseudobulk Data
###############################################################################
cat("\n---Creating Pseudobulk ---\n")

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

# Helpers
strip_version <- function(x) sub("\\.\\d+$", "", x)
canon <- function(x) { x <- trimws(as.character(x)); toupper(x) }
looks_ensembl <- function(x) any(grepl("^ENSG", x))

spe_neuronWGCNA <- spe_ruv[, spe_ruv$Cell == "Neuron"]
cat("ROI-level data:", ncol(spe_neuronWGCNA), "oligodendrocyte ROIs\n")

# Get expression
if ("logcounts" %in% assayNames(spe_neuronWGCNA)) {
  expr0 <- assay(spe_neuronWGCNA, "logcounts")
} else {
  expr0 <- cpm(assay(spe_neuronWGCNA, "counts"), log = TRUE, prior.count = 1)
}
expr0 <- as.matrix(expr0)

gene_ids_raw <- rownames(expr0)

###############################################################################
### ID Standardization
###############################################################################
cat("Standardizing gene IDs to ENSEMBL (no version)...\n")

if (looks_ensembl(gene_ids_raw)) {
  # Already ENSEMBL so strip version and average duplicates
  ensg <- strip_version(gene_ids_raw)
  
  summed <- rowsum(expr0, group = ensg, reorder = FALSE)
  counts <- table(ensg)
  denom <- as.numeric(counts[match(rownames(summed), names(counts))])
  logcpm <- summed / denom
  rownames(logcpm) <- rownames(summed)
  
} else {
  # Convert SYMBOL → ENSEMBL
  symbols <- canon(gene_ids_raw)
  
  map <- AnnotationDbi::select(org.Hs.eg.db,
                               keys = unique(symbols),
                               keytype = "SYMBOL",
                               columns = "ENSEMBL")
  map <- map[!is.na(map$ENSEMBL), ]
  map <- map[!duplicated(map$SYMBOL), c("SYMBOL", "ENSEMBL")]
  rownames(map) <- canon(map$SYMBOL)
  
  ensg_vec <- strip_version(map[symbols, "ENSEMBL"])
  keep <- !is.na(ensg_vec)
  
  expr1 <- expr0[keep, , drop = FALSE]
  expr1 <- as.matrix(expr1)
  ensg_vec <- ensg_vec[keep]
  
  summed <- rowsum(expr1, group = ensg_vec, reorder = FALSE)
  counts <- table(ensg_vec)
  denom <- as.numeric(counts[match(rownames(summed), names(counts))])
  logcpm <- summed / denom
  rownames(logcpm) <- rownames(summed)
}

stopifnot(!any(duplicated(rownames(logcpm))))
cat("Standardized:", nrow(logcpm), "genes (ENSEMBL) x", ncol(logcpm), "ROIs\n")

###############################################################################
### Create pseudobulk from standardized logcpm
###############################################################################

patients <- spe_neuronWGCNA$Patient
unique_patients <- unique(patients)

pseudobulk_mat <- matrix(0, nrow = nrow(logcpm), ncol = length(unique_patients))
rownames(pseudobulk_mat) <- rownames(logcpm)  
colnames(pseudobulk_mat) <- unique_patients

for (pat in unique_patients) {
  roi_idx <- which(patients == pat)
  pseudobulk_mat[, pat] <- rowMeans(logcpm[, roi_idx, drop = FALSE])
}

cat("Pseudobulk:", ncol(pseudobulk_mat), "patients\n")

# Patient metadata
patient_meta <- data.frame(
  Patient = unique_patients,
  Group = sapply(unique_patients, function(p) {
    unique(spe_neuronWGCNA$Group[spe_neuronWGCNA$Patient == p])
  }),
  stringsAsFactors = FALSE
)

cat("  MSA:", sum(patient_meta$Group == "MSA"), "patients\n")
cat("  PD:", sum(patient_meta$Group == "PD"), "patients\n")

###############################################################################
###Filter to variable genes
###############################################################################
cat("\n--- Filtering genes ---\n")

gene_mad <- apply(pseudobulk_mat, 1, mad)
n_keep <- min(3500, nrow(pseudobulk_mat))
top_genes <- names(sort(gene_mad, decreasing = TRUE)[1:n_keep])

expr_mat <- t(pseudobulk_mat[top_genes, ])  # Rows = patients, cols = genes

cat("WGCNA input:", nrow(expr_mat), "patients ×", ncol(expr_mat), "genes\n")

# QC step
gsg <- goodSamplesGenes(expr_mat, verbose = 3)
if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0) {
    cat("Removing", sum(!gsg$goodGenes), "problematic genes\n")
    expr_mat <- expr_mat[, gsg$goodGenes]
  }
  if (sum(!gsg$goodSamples) > 0) {
    cat("WARNING: Removing", sum(!gsg$goodSamples), "problematic samples\n")
    expr_mat <- expr_mat[gsg$goodSamples, , drop = FALSE]
    patient_meta <- patient_meta[match(rownames(expr_mat), patient_meta$Patient), , drop = FALSE]
  }
}

if (!gsg$allOK) {
  expr_mat <- expr_mat[, gsg$goodGenes]
  
  module_colors <- module_colors[gsg$goodGenes]
  

  names(module_colors) <- colnames(expr_mat)
}

cat("\n--- Sample clustering ---\n")

sampleTree <- hclust(dist(expr_mat), method = "average")

pdf(file.path(OUTPUT_DIR, "11a_WGCNA_sampleTree_pseudobulk.pdf"), width = 12, height = 6)
par(mar = c(2, 5, 2, 1))
plot(sampleTree, main = "Sample clustering (Neurons - Pseudobulk)",
     sub = "", xlab = "", cex = 0.7)
abline(h = 100, col = "red", lty = 2)
dev.off()

###############################################################################
###Power selection
###############################################################################
cat("\n--- Selecting soft-thresholding power ---\n")

powers <- c(seq(1, 10, by = 1), seq(12, 20, by = 2), seq(22, 30, by = 2))

sft <- pickSoftThreshold(
  expr_mat,
  powerVector = powers,
  networkType = "signed",
  corFnc = "bicor",
  verbose = 5
)

sft_df <- data.frame(
  Power = sft$fitIndices[, 1],
  SFT.R.sq = -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  Slope = sft$fitIndices[, 3],
  MeanConnectivity = sft$fitIndices[, 5],
  MedianConnectivity = sft$fitIndices[, 6]
)

cat("\nPower selection table:\n")
print(sft_df)

power_at_08 <- sft_df$Power[which(sft_df$SFT.R.sq >= 0.8)[1]]
power_at_07 <- sft_df$Power[which(sft_df$SFT.R.sq >= 0.7)[1]]
power_at_06 <- sft_df$Power[which(sft_df$SFT.R.sq >= 0.6)[1]]

max_rsq <- max(sft_df$SFT.R.sq, na.rm = TRUE)
power_at_max <- sft_df$Power[which.max(sft_df$SFT.R.sq)]

cat("\nThresholds:\n")
cat("  R² ≥ 0.8:", if (is.na(power_at_08)) "NONE" else power_at_08, "\n")
cat("  R² ≥ 0.7:", if (is.na(power_at_07)) "NONE" else power_at_07, "\n")
cat("  R² ≥ 0.6:", if (is.na(power_at_06)) "NONE" else power_at_06, "\n")
cat("  Maximum R²:", round(max_rsq, 3), "at power", power_at_max, "\n")

if (!is.na(power_at_08)) {
  power_use <- power_at_08
  cat("\nUsing power", power_use, "(R² ≥ 0.8)\n")
} else if (!is.na(power_at_07)) {
  power_use <- power_at_07
  cat("\nWARNING: Using power", power_use, "(R² ≥ 0.7, never reached 0.8)\n")
} else if (max_rsq >= 0.5) {
  power_use <- power_at_max
  cat("\nWARNING: Using power", power_use, "(max R² =", round(max_rsq, 3), ")\n")
} else {
  power_use <- 6
  cat("\nWARNING: Max R² =", round(max_rsq, 3), "< 0.5\n")
  cat("Data may not be suitable for WGCNA. Using conservative power =", power_use, "\n")
}

mean_conn <- sft_df$MeanConnectivity[sft_df$Power == power_use]
cat("Mean connectivity at selected power:", round(mean_conn, 1), "\n")

if (mean_conn < 5) {
  cat("WARNING: Connectivity very low. Adjusting power...\n")
  candidates <- sft_df[sft_df$MeanConnectivity >= 10 & sft_df$MeanConnectivity <= 50, ]
  if (nrow(candidates) > 0) {
    power_use <- candidates$Power[which.max(candidates$SFT.R.sq)]
    mean_conn <- sft_df$MeanConnectivity[sft_df$Power == power_use]
    cat("Adjusted to power", power_use, "(connectivity =", round(mean_conn, 1), ")\n")
  }
}

# Plots
sft_df$selected <- sft_df$Power == power_use

p1 <- ggplot(sft_df[order(sft_df$Power), ], aes(x = Power, y = SFT.R.sq)) +
  geom_point(aes(color = selected, size = selected)) +
  geom_line(color = "dodgerblue") +
  geom_hline(yintercept = c(0.8, 0.7), linetype = c("solid", "dashed"),
             color = c("red", "orange"), alpha = 0.5) +
  scale_color_manual(values = c("FALSE" = "dodgerblue", "TRUE" = "red")) +
  scale_size_manual(values = c("FALSE" = 3, "TRUE" = 6)) +
  labs(title = "Scale-free topology fit",
       subtitle = paste0("Selected: power = ", power_use, ", R² = ",
                         round(sft_df$SFT.R.sq[sft_df$Power == power_use], 3)),
       x = "Soft-thresholding power", y = "Signed R²") +
  theme_classic(base_size = 12) + theme(legend.position = "none")

p2 <- ggplot(sft_df[order(sft_df$Power), ], aes(x = Power, y = MeanConnectivity)) +
  geom_point(aes(color = selected, size = selected)) +
  geom_line(color = "darkorange") +
  scale_color_manual(values = c("FALSE" = "darkorange", "TRUE" = "red")) +
  scale_size_manual(values = c("FALSE" = 3, "TRUE" = 6)) +
  scale_y_log10() +
  labs(title = "Mean connectivity",
       subtitle = paste0("Selected: ", round(mean_conn, 1), " connections"),
       x = "Soft-thresholding power", y = "Mean connectivity (log)") +
  theme_classic(base_size = 12) + theme(legend.position = "none")

ggsave(file.path(OUTPUT_DIR, "11b_WGCNA_power_selection.pdf"),
       p1 | p2, width = 12, height = 5)

write.csv(sft_df, file.path(OUTPUT_DIR, "WGCNA_power_selection.csv"), row.names = FALSE)

if (max_rsq < 0.5) {
  cat("\n", paste(rep("=", 70), collapse = ""), "\n", sep = "")
  cat("CRITICAL: Scale-free topology R² < 0.5 - data not suitable for WGCNA\n")
  cat(paste(rep("=", 70), collapse = ""), "\n", sep = "")
  cat("\nRecommendation: Focus on individual DEG analysis instead.\n")
  cat("This is a valid scientific finding - not all datasets have modular structure.\n\n")
}

###############################################################################
###Build network
###############################################################################
cat("\n--- Network construction ---\n")

net <- blockwiseModules(
  expr_mat,
  power = power_use,
  networkType = "signed",
  corType = "bicor",
  
  deepSplit = 3,           #3 = more granular (more, smaller modules), 1 = conservative, fewer large modules Neurons =3
  minModuleSize = 55,      # is this too high? Neurons =55 oligo =40
  mergeCutHeight = 0.25,     #mergeCutHeight = 0.25 → only merge very similar modules, 0.50 → merges many modules (can merge biology away) Neurons =0.25 oligo =0.2
  pamRespectsDendro = TRUE,   #N = true
  maxBlockSize = ncol(expr_mat),
  saveTOMs = FALSE,       #N=false
  
  numericLabels = TRUE,
  verbose = 3
)

cat("\nNetwork built successfully\n")
cat("  Genes:", length(net$colors), "\n")
cat("  Module numbers:", paste(sort(unique(net$colors)), collapse = ", "), "\n")

# Store gene names BEFORE conversion
gene_names <- colnames(expr_mat)

cat("\nGene names (first 10):", paste(head(gene_names, 10), collapse = ", "), "\n")

# Verify alignment
if (length(gene_names) != length(net$colors)) {
  stop("ERROR: Mismatch between gene names (", length(gene_names), 
       ") and module assignments (", length(net$colors), ")")
}

# Convert numeric labels to colors
module_colors <- labels2colors(net$colors)

# Re-assign names  as labels2colors drops them
names(module_colors) <- gene_names

cat("\nModule color assignment:\n")
cat("  module_colors length:", length(module_colors), "\n")
cat("  module_colors are named:", !is.null(names(module_colors)), "\n")
cat("  First 10 genes:", paste(head(names(module_colors), 10), collapse = ", "), "\n")
cat("  Unique colors:", paste(sort(unique(module_colors)), collapse = ", "), "\n")
cat("\nModule sizes:\n")
print(table(module_colors))

# Verify no names were lost
stopifnot(!is.null(names(module_colors)))
stopifnot(length(names(module_colors)) == length(module_colors))

###############################################################################
### Save module assignments
###############################################################################

module_df <- data.frame(
  gene = names(module_colors), 
  module = module_colors,
  stringsAsFactors = FALSE
)

cat("\nModule dataframe:\n")
cat("  Rows:", nrow(module_df), "\n")
cat("  Columns:", ncol(module_df), "\n")
cat("  First 5 rows:\n")
print(head(module_df, 5))

write.csv(module_df, 
          file.path(OUTPUT_DIR, "WGCNA_gene_modules_QCd.csv"), 
          row.names = FALSE)

cat("\nSaved module assignments to WGCNA_gene_modules_QCd.csv\n")

###############################################################################
### Dendrogram
###############################################################################
pdf(file.path(OUTPUT_DIR, "11c_WGCNA_dendrogram.pdf"), width = 14, height = 6)
plotDendroAndColors(net$dendrograms[[1]], 
                    module_colors[net$blockGenes[[1]]],
                    "Modules",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram (Neurons)")
dev.off()

cat("\n=== Module Detection Complete ===\n")

###############################################################################
### Module eigengenes and disease association
###############################################################################
cat("\n---  Module-disease associations ---\n")

MEs <- moduleEigengenes(expr_mat, module_colors)$eigengenes
MEs <- orderMEs(MEs)

# Trait matrix for visualisation
trait_data <- data.frame(
  Group_MSA = as.numeric(patient_meta$Group == "MSA"),
  Group_PD  = as.numeric(patient_meta$Group == "PD"),
  row.names = patient_meta$Patient
)

# Correlations (patient-level)
module_trait_cor <- cor(MEs, trait_data, use = "pairwise.complete.obs")
module_trait_pval <- corPvalueStudent(module_trait_cor, nrow(expr_mat))

###############################################################################
### ROI-level eigengenes
###############################################################################
cat("\n---  Computing ROI-level eigengenes ---\n")

stopifnot(exists("logcpm"), exists("module_colors"))

roi_expr <- t(logcpm)  # Rows = ROIs, cols = genes (ENSEMBL)
roi_genes <- colnames(roi_expr)
net_genes <- names(module_colors)

cat("  ROI genes:", length(roi_genes), "\n")
cat("  Network genes:", length(net_genes), "\n")

common_genes <- intersect(roi_genes, net_genes)
cat("  Common genes:", length(common_genes), "\n")

if (length(common_genes) < 100) {
  warning("Very few genes matched (", length(common_genes), 
          "). Check if IDs were standardized correctly.")
}

# Subset to common genes and assign colors
roi_expr_common <- roi_expr[, common_genes, drop = FALSE]
colors_roi <- module_colors[common_genes]

# Check color distribution
cat("  Module assignments:\n")
print(table(colors_roi))

if (sum(colors_roi != "grey") < 100) {
  warning("Most genes assigned to grey module. Network may be weak.")
}

# Compute ROI-level eigengenes
MEs_roi <- moduleEigengenes(roi_expr_common, colors = colors_roi)$eigengenes
MEs_roi <- orderMEs(MEs_roi)

# Remove modules with too few genes
module_sizes <- table(colors_roi[colors_roi != "grey"])
min_genes <- 30
small_modules <- names(module_sizes)[module_sizes < min_genes]

if (length(small_modules) > 0) {
  cat("  Removing", length(small_modules), "small modules (<", min_genes, "genes)\n")
  keep_cols <- setdiff(colnames(MEs_roi), c("MEgrey", paste0("ME", small_modules)))
  MEs_roi <- MEs_roi[, keep_cols, drop = FALSE]
}

cat("  Final ROI-level eigengenes:", ncol(MEs_roi), "modules x", nrow(MEs_roi), "ROIs\n")

###############################################################################
### Statistical testing
###############################################################################
cat("\n---  Statistical testing ---\n")

module_trait_results <- data.frame()

for (module_name in names(MEs)) {
  test_data <- data.frame(
    Eigengene = MEs[[module_name]],
    Group = factor(patient_meta$Group, levels = c("PD", "MSA"))
  )
  
  test <- t.test(Eigengene ~ Group, data = test_data)  # Welch
  
  pd_vals  <- test_data$Eigengene[test_data$Group == "PD"]
  msa_vals <- test_data$Eigengene[test_data$Group == "MSA"]
  
  # Manual Welch SE
  SE <- sqrt(var(pd_vals) / length(pd_vals) + var(msa_vals) / length(msa_vals))
  
  result <- data.frame(
    Module  = module_name,
    Beta    = mean(msa_vals) - mean(pd_vals),
    SE      = SE,
    t       = as.numeric(test$statistic),
    P.Value = test$p.value,
    PD_mean = mean(pd_vals),
    MSA_mean= mean(msa_vals),
    PD_n    = length(pd_vals),
    MSA_n   = length(msa_vals)
  )
  
  module_trait_results <- rbind(module_trait_results, result)
}

module_trait_results$adj.P.Val <- p.adjust(module_trait_results$P.Value, method = "BH")
module_trait_results <- module_trait_results[order(module_trait_results$P.Value), ]

write.csv(module_trait_results, 
          file.path(OUTPUT_DIR, "WGCNA_module_trait_results.csv"), 
          row.names = FALSE)

cat("\n=== Module-Disease Results (patient-level) ===\n")
print(module_trait_results[, c("Module", "Beta", "PD_mean", "MSA_mean", "P.Value", "adj.P.Val")])

cat("\nSummary:\n")
cat("  Total modules:", nrow(module_trait_results), "\n")
cat("  FDR-significant (adj.P < 0.05):", sum(module_trait_results$adj.P.Val < 0.05), "\n")
cat("  Nominally significant (P < 0.05):", sum(module_trait_results$P.Value < 0.05), "\n")

###############################################################################
###Visualizations
###############################################################################
cat("\n---Creating visualizations ---\n")

# Heatmap
cor_mat  <- module_trait_cor
pval_mat <- module_trait_pval

if (!is.matrix(cor_mat))  cor_mat  <- as.matrix(cor_mat)
if (!is.matrix(pval_mat)) pval_mat <- as.matrix(pval_mat)

# Module names remove 'ME' prefix
rownames(cor_mat)  <- gsub("^ME", "", rownames(cor_mat))
rownames(pval_mat) <- gsub("^ME", "", rownames(pval_mat))


sig_labels <- matrix(
  ifelse(pval_mat < 0.001, "***",
         ifelse(pval_mat < 0.01,  "**",
                ifelse(pval_mat < 0.05, "*", ""))),
  nrow = nrow(pval_mat), ncol = ncol(pval_mat),
  dimnames = dimnames(pval_mat)
)

cor_labels <- matrix(
  sprintf("%.2f", cor_mat),
  nrow = nrow(cor_mat), ncol = ncol(cor_mat),
  dimnames = dimnames(cor_mat)
)

pdf(file.path(OUTPUT_DIR, "11d_WGCNA_module_trait_heatmap_patient.pdf"), width = 6, height = 8)

col_fun <- colorRamp2(c(-0.6, 0, 0.6), c("#3B5BA5", "white", "#B40426"))

ht <- Heatmap(
  cor_mat,
  name = "Correlation",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_gp = gpar(fontsize = 11, fontface = "bold"),
  column_names_rot = 0,
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.text(paste0(cor_labels[i, j], "\n", sig_labels[i, j]),
              x, y, gp = gpar(fontsize = 8))
  },
  column_title = "Module–Group Associations\n(Neurons - Pseudobulk, patient-level)",
  column_title_gp = gpar(fontsize = 12, fontface = "bold")
)

draw(ht, heatmap_legend_side = "right")
dev.off()

# Eigengene network (patient-level)
ME_cor <- cor(MEs, use = "pairwise.complete.obs")
ME_tree <- hclust(as.dist(1 - ME_cor))

colnames(ME_cor) <- gsub("^ME", "", colnames(ME_cor))
rownames(ME_cor) <- gsub("^ME", "", rownames(ME_cor))

pdf(file.path(OUTPUT_DIR, "11e_WGCNA_eigengene_network.pdf"), width = 7, height = 7)
ht2 <- Heatmap(
  ME_cor,
  name = "Pearson r",
  col = colorRamp2(c(-1, 0, 1), c("#3B5BA5", "white", "#B40426")),
  cluster_rows = ME_tree,
  cluster_columns = ME_tree,
  row_names_gp = gpar(fontsize = 9, fontface = "bold"),
  column_names_gp = gpar(fontsize = 9, fontface = "bold"),
  column_title = "Module eigengene correlations (patient-level)"
)
draw(ht2)
dev.off()

# Boxplots for top modules by P
top_modules <- head(module_trait_results$Module, 4)

pdf(file.path(OUTPUT_DIR, "11f_WGCNA_eigengene_boxplots.pdf"), width = 12, height = 8)
par(mfrow = c(2, 2))
for (mod in top_modules) {
  eigengene <- MEs[[mod]]
  group <- patient_meta$Group
  mod_res <- module_trait_results[module_trait_results$Module == mod, ]
  
  boxplot(eigengene ~ group,
          main = paste0(mod, "\nP = ", format.pval(mod_res$P.Value, digits = 3),
                        ", Beta = ", round(mod_res$Beta, 3)),
          xlab = "Disease group", ylab = "Module eigengene",
          col = c("PD" = "#4575B4", "MSA" = "#D73027"),
          outline = FALSE)
  
  stripchart(eigengene ~ group, vertical = TRUE, method = "jitter",
             add = TRUE, pch = 16, cex = 1.2, col = "black")
}
dev.off()

###############################################################################
### WGCNA Eigengene Correlation Network
###############################################################################

MEs_ordered <- orderMEs(MEs)
ME_cor_pub <- cor(MEs_ordered, use = "pairwise.complete.obs", method = "pearson")
ME_tree_pub <- hclust(as.dist(1 - ME_cor_pub))

colnames(ME_cor_pub) <- gsub("^ME", "", colnames(ME_cor_pub))
rownames(ME_cor_pub) <- gsub("^ME", "", rownames(ME_cor_pub))

max_abs <- max(abs(ME_cor_pub))
col_fun_pub <- circlize::colorRamp2(c(-max_abs, 0, max_abs), c("#3B5BA5", "white", "#B40426"))

pdf(file.path(OUTPUT_DIR, "11g_WGCNA_eigengene_network_publication.pdf"),
    width = 7, height = 7)
ht_pub <- ComplexHeatmap::Heatmap(
  ME_cor_pub,
  name = "Pearson r",
  col = col_fun_pub,
  cluster_rows = ME_tree_pub,
  cluster_columns = ME_tree_pub,
  row_names_gp = grid::gpar(fontsize = 9, fontface = "bold"),
  column_names_gp = grid::gpar(fontsize = 9, fontface = "bold"),
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid::grid.text(sprintf("%.2f", ME_cor_pub[i, j]),
                    x, y, gp = grid::gpar(fontsize = 6))
  },
  heatmap_legend_param = list(
    title = "Pearson r",
    title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
    labels_gp = grid::gpar(fontsize = 8)
  )
)
ComplexHeatmap::draw(ht_pub)
dev.off()
###############################################################################
### CORRELATION FIGURE CLINICAL TRAITS
###############################################################################
cat("\n--- ROI-level trait correlations (clinical only) ---\n")
stopifnot(exists("MEs_roi")) 

strip_dcc <- function(x) sub("\\.dcc$", "", as.character(x), ignore.case = TRUE)

# Read ROI-level traits
trait_df <- readxl::read_excel(TRAIT_FILE) |> as.data.frame()
stopifnot(all(c("ROI","Patient","Group","Age","Sex", "Tau", "AB") %in% colnames(trait_df)))

trait_df$ROI <- trimws(strip_dcc(trait_df$ROI))
rownames(trait_df) <- trait_df$ROI

# Align traits to ROI-level eigengenes (rows = ROIs)
me_rois_raw   <- rownames(MEs_roi)
me_rois_clean <- strip_dcc(me_rois_raw)

missing <- setdiff(me_rois_clean, rownames(trait_df))
if (length(missing) > 0) {
  stop("ROIs missing from trait file after stripping .dcc: ",
       paste(head(missing, 10), collapse = ", "),
       if (length(missing) > 10) paste0(" ... +", length(missing) - 10, " more") else "")
}
trait_df <- trait_df[me_rois_clean, , drop = FALSE]
rownames(trait_df) <- me_rois_raw
stopifnot(all(rownames(trait_df) == rownames(MEs_roi)))

# Build a  numeric trait matrix

id_cols      <- c("ROI","Patient")  # anything identifying
keep_binary  <- c("Group","Sex")  # will map to 0/1 columns
keep_numeric <- c(
  "Age", "PMI", "Duration",
  "Duration", "Tau", "AB"
)

# Collect
keep_binary  <- intersect(keep_binary,  colnames(trait_df))
keep_numeric <- intersect(keep_numeric, colnames(trait_df))

# Numeric matrix
trait_numeric <- data.frame(row.names = rownames(trait_df))

# Add binary encodings (0/1)
if ("Group" %in% keep_binary) {
  # Code PD=0, MSA=1
  g <- as.character(trait_df$Group)
  g <- trimws(toupper(g))
  trait_numeric$Group_MSA <- as.numeric(factor(g, levels = c("PD","MSA"))) - 1
}

if ("Sex" %in% keep_binary) {
  # Code M=0, F=1
  s <- as.character(trait_df$Sex)
  s <- trimws(toupper(s))
  # Try to harmonize common forms
  s[s %in% c("MALE","M")]   <- "M"
  s[s %in% c("FEMALE","F")] <- "F"
  # Only encode if truly binary
  ux <- sort(unique(na.omit(s)))
  if (length(ux) == 2) {
    trait_numeric$Sex_Female <- as.numeric(factor(s, levels = c("M","F"))) - 0  
    trait_numeric$Sex_Female <- ifelse(s == "F", 1, ifelse(s == "M", 0, NA))
  }
}

# Add numeric traits
numify <- function(x) suppressWarnings(as.numeric(x))
for (nm in keep_numeric) {
  v <- numify(trait_df[[nm]])
  if (sd(v, na.rm = TRUE) > 0) {
    trait_numeric[[nm]] <- v
  }
}

# Remove any all-NA columns
if (ncol(trait_numeric) == 0 || all(colSums(!is.na(trait_numeric)) == 0)) {
  stop("No usable clinical traits found after filtering. Check column names and content in TRAIT_FILE.")
}
trait_numeric <- trait_numeric[, colSums(!is.na(trait_numeric)) > 0, drop = FALSE]

# Module–trait correlations
MEs_roi_ord <- WGCNA::orderMEs(MEs_roi)
module_trait_cor_roi <- cor(MEs_roi_ord, trait_numeric, use = "pairwise.complete.obs")
module_trait_p_roi   <- WGCNA::corPvalueStudent(module_trait_cor_roi, nSamples = nrow(MEs_roi_ord))

rownames(module_trait_cor_roi) <- gsub("^ME", "", rownames(module_trait_cor_roi))
rownames(module_trait_p_roi)   <- gsub("^ME", "", rownames(module_trait_p_roi))

top_k <- min(10, ncol(module_trait_cor_roi))
score_cols <- apply(abs(module_trait_cor_roi), 2, max, na.rm = TRUE)
top_cols <- names(sort(score_cols, decreasing = TRUE))[1:top_k]

cor_plot <- module_trait_cor_roi[, top_cols, drop = FALSE]
p_plot   <- module_trait_p_roi[,   top_cols, drop = FALSE]

# --- Make sure these are matrices with dimnames ---
cor_plot <- as.matrix(cor_plot)
p_plot   <- as.matrix(p_plot)

# Label matrices (preserve dims + dimnames)
cor_lab_mat <- matrix(
  sprintf("%.2f", as.numeric(cor_plot)),
  nrow = nrow(cor_plot),
  ncol = ncol(cor_plot),
  dimnames = dimnames(cor_plot)
)

p_lab_mat <- matrix(
  ifelse(p_plot < 0.001, "***",
         ifelse(p_plot < 0.01,  "**",
                ifelse(p_plot < 0.05, "*", ""))),
  nrow = nrow(p_plot),
  ncol = ncol(p_plot),
  dimnames = dimnames(p_plot)
)

# Color scale
max_abs_trait <- max(abs(cor_plot), na.rm = TRUE)
col_fun_trait <- circlize::colorRamp2(c(-max_abs_trait, 0, max_abs_trait),
                                      c("#3B5BA5", "white", "#B40426"))

# Build heatmap; note
ht_trait <- ComplexHeatmap::Heatmap(
  cor_plot,
  name = "Module–trait r",
  col = col_fun_trait,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_gp = grid::gpar(fontsize = 9, fontface = "bold"),
  column_names_gp = grid::gpar(fontsize = 9, fontface = "bold"),
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid::grid.text(
      paste0(cor_lab_mat[i, j], "\n", p_lab_mat[i, j]),
      x, y, gp = grid::gpar(fontsize = 7)
    )
  }
)

pdf(file.path(OUTPUT_DIR, "11h_WGCNA_module_trait_heatmap_ROI_CLINICAL.pdf"), width = 7, height = 6)
ComplexHeatmap::draw(ht_trait, heatmap_legend_side = "right")
dev.off()

cat("Clinical traits included: ", paste(colnames(cor_plot), collapse = ", "), "\n")

###############################################################################
###Hub genes with SYMBOL conversion
###############################################################################

cat("\n---Exporting hub genes (with gene symbols) ---\n")

library(org.Hs.eg.db)

sig_modules <- gsub("^ME", "", module_trait_results$Module[module_trait_results$P.Value < 0.97])

if (length(sig_modules) > 0) {
  cat("Exporting hub genes for", length(sig_modules), "modules (P < 0.99)\n")
  
  for (module_col in sig_modules) {
    
    # Get genes in module
    module_genes <- module_df$gene[module_df$module == module_col]
    if (length(module_genes) == 0) next
    
    cat("  ", module_col, ":", length(module_genes), "genes\n")
    
    # Calculate module membership
    eigengene <- MEs[[paste0("ME", module_col)]]
    gene_expr <- expr_mat[, module_genes, drop = FALSE]
    
    mm <- cor(gene_expr, eigengene, use = "pairwise.complete.obs")
    
    # Create hub genes dataframe
    hub_genes <- data.frame(
      ENSEMBL = module_genes,
      module = module_col,
      MM = mm[, 1],
      stringsAsFactors = FALSE
    )
    
    ###########################################################################
    ### Convert ENSEMBL to SYMBOL
    ###########################################################################
    
    # Map ENSEMBL → SYMBOL
    gene_map <- mapIds(
      org.Hs.eg.db,
      keys = hub_genes$ENSEMBL,
      column = "SYMBOL",
      keytype = "ENSEMBL",
      multiVals = "first"
    )
    
    hub_genes$SYMBOL <- gene_map[hub_genes$ENSEMBL]
    
    # For unmapped genes, keep ENSEMBL as symbol
    unmapped <- is.na(hub_genes$SYMBOL)
    if (any(unmapped)) {
      hub_genes$SYMBOL[unmapped] <- hub_genes$ENSEMBL[unmapped]
      cat("    ", sum(unmapped), "genes could not be mapped to symbols\n")
    }
    
    # Reorder columns: SYMBOL first, then ENSEMBL, then MM
    hub_genes <- hub_genes[, c("SYMBOL", "ENSEMBL", "module", "MM")]
    
    # Sort by absolute MM (most correlated with module)
    hub_genes <- hub_genes[order(-abs(hub_genes$MM)), ]
    
    # Save top 50
    write.csv(head(hub_genes, 50),
              file.path(OUTPUT_DIR, paste0("WGCNA_hub_genes_", module_col, ".csv")),
              row.names = FALSE)
    
    cat("    Top 5 hub genes:", paste(head(hub_genes$SYMBOL, 5), collapse = ", "), "\n")
  }
  
  cat("\nHub gene files saved with both SYMBOL and ENSEMBL IDs\n")
}

###############################################################################
### Update the main module assignment file with symbols
###############################################################################

cat("\n--- Adding gene symbols to module assignment file ---\n")

# Convert all module genes to symbols
all_gene_map <- mapIds(
  org.Hs.eg.db,
  keys = module_df$gene,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

module_df$SYMBOL <- all_gene_map[module_df$gene]

# For unmapped, keep ENSEMBL
unmapped_all <- is.na(module_df$SYMBOL)
if (any(unmapped_all)) {
  module_df$SYMBOL[unmapped_all] <- module_df$gene[unmapped_all]
  cat("  ", sum(unmapped_all), "of", nrow(module_df), 
      "genes could not be mapped to symbols\n")
}

# Reorder columns: SYMBOL, ENSEMBL, module
module_df_export <- data.frame(
  SYMBOL = module_df$SYMBOL,
  ENSEMBL = module_df$gene,
  module = module_df$module,
  stringsAsFactors = FALSE
)

write.csv(module_df_export,
          file.path(OUTPUT_DIR, "WGCNA_gene_modules_with_symbols.csv"),
          row.names = FALSE)

cat("Saved full module assignments with gene symbols\n")

###############################################################################
### Final Summary
###############################################################################

cat("\n=== WGCNA Analysis Complete ===\n")
cat("Power used:", power_use, "(R² =", round(sft_df$SFT.R.sq[sft_df$Power == power_use], 3), ")\n")
cat("Modules detected:", length(unique(module_colors)), "\n")
cat("Significant modules (FDR < 0.05):", sum(module_trait_results$adj.P.Val < 0.05), "\n")
cat("Nominally significant (P < 0.05):", sum(module_trait_results$P.Value < 0.05), "\n\n")
















###############################################################################
### WGCNA Follow-up Analysis
###############################################################################

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(dplyr)
  library(enrichplot)
  library(ggplot2)
  library(tidyr)
  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(scales)
})

cat("\n=== WGCNA Follow-up Analysis ===\n")

###############################################################################
### Load data and verify
###############################################################################

cat("\n--- Loading data ---\n")

# Load module assignments (with SYMBOLS)
gene_modules <- read.csv(file.path(OUTPUT_DIR, "WGCNA_gene_modules_with_symbols.csv"))
cat("Gene modules loaded:", nrow(gene_modules), "genes\n")
cat("Columns:", paste(colnames(gene_modules), collapse = ", "), "\n")

# Load module-trait results
trait_tests <- read.csv(file.path(OUTPUT_DIR, "WGCNA_module_trait_results.csv"))
cat("Module-trait results loaded:", nrow(trait_tests), "modules\n")

# Identify significant modules (use P < 0.10 for follow-up)
sig_modules <- trait_tests %>%
  filter(P.Value < 0.96) %>%
  mutate(module_color = gsub("^ME", "", Module)) %>%
  pull(module_color)

cat("\nModules for follow-up (P < 0.05):", paste(sig_modules, collapse = ", "), "\n")
cat("  FDR-significant (adj.P < 0.05):", sum(trait_tests$adj.P.Val < 0.05), "\n")

if (length(sig_modules) == 0) {
  cat("\nWARNING: No modules with P < 0.10\n")
  cat("Consider:\n")
  cat("  1. Using all modules for exploratory analysis\n")
  cat("  2. Focusing on top 3 modules by P-value\n")
  cat("  3. Prioritizing hub gene analysis over pathway enrichment\n\n")
  
  sig_modules <- trait_tests %>%
    arrange(P.Value) %>%
    slice(1:min(16, nrow(trait_tests))) %>%
    mutate(module_color = gsub("^ME", "", Module)) %>%
    pull(module_color)
  
  cat("Proceeding with top", length(sig_modules), "modules:", 
      paste(sig_modules, collapse = ", "), "\n")
}

###############################################################################
### Pathway Enrichment (GO, KEGG, Reactome)
###############################################################################

cat("\n---  Pathway enrichment ---\n")

# Background: all genes in network
all_genes_symbol <- unique(gene_modules$SYMBOL)
all_genes_entrez <- mapIds(
  org.Hs.eg.db,
  keys = all_genes_symbol,
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)
bg_entrez <- na.omit(unique(all_genes_entrez))

cat("Background: ", length(all_genes_symbol), " symbols, ", 
    length(bg_entrez), " ENTREZ IDs\n", sep = "")

# Function: pathway enrichment for one module
enrich_module <- function(module_color, top_n = 50) {
  cat("\n  Enriching module:", module_color, "\n")
  
  # Get hub genes for this module
  hub_file <- file.path(OUTPUT_DIR, paste0("WGCNA_hub_genes_", module_color, ".csv"))
  
  if (!file.exists(hub_file)) {
    cat("    Hub gene file not found, using all module genes\n")
    mod_genes_symbol <- gene_modules %>%
      filter(module == module_color) %>%
      pull(SYMBOL) %>%
      head(top_n)
  } else {
    hub_data <- read.csv(hub_file)
    mod_genes_symbol <- gene_modules %>%
      filter(module == module_color) %>%
      pull(SYMBOL)
    
  }
  
  cat("    Genes:", length(mod_genes_symbol), "\n")
  
  if (length(mod_genes_symbol) < 5) {
    cat("    Too few genes, skipping\n")
    return(NULL)
  }
  
  # GO enrichment
  go_res <- tryCatch({
    enrichGO(
      gene = mod_genes_symbol,
      universe = all_genes_symbol,
      OrgDb = org.Hs.eg.db,
      keyType = "SYMBOL",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1
    )
  }, error = function(e) {
    cat("    GO error:", e$message, "\n")
    NULL
  })
  
  # Convert to ENTREZ for KEGG/Reactome
  mod_genes_entrez <- mapIds(
    org.Hs.eg.db,
    keys = mod_genes_symbol,
    column = "ENTREZID",
    keytype = "SYMBOL",
    multiVals = "first"
  )
  mod_genes_entrez <- na.omit(unique(mod_genes_entrez))
  
  cat("    Mapped to", length(mod_genes_entrez), "ENTREZ IDs\n")
  
  # KEGG enrichment
  kegg_res <- tryCatch({
    enrichKEGG(
      gene = mod_genes_entrez,
      universe = bg_entrez,
      organism = "hsa",
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1
    )
  }, error = function(e) {
    cat("    KEGG error:", e$message, "\n")
    NULL
  })
  
  # Reactome enrichment
  reactome_res <- tryCatch({
    enrichPathway(
      gene = mod_genes_entrez,
      universe = bg_entrez,
      organism = "human",
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1
    )
  }, error = function(e) {
    cat("    Reactome error:", e$message, "\n")
    NULL
  })
  
  list(GO = go_res, KEGG = kegg_res, Reactome = reactome_res)
}

# Run enrichment
enrichment_results <- lapply(sig_modules, enrich_module)
names(enrichment_results) <- sig_modules

# Save top results
for (mod in sig_modules) {
  res <- enrichment_results[[mod]]
  if (is.null(res)) next
  
  # GO
  if (!is.null(res$GO) && nrow(res$GO@result) > 0) {
    go_top <- res$GO@result %>%
      arrange(pvalue) %>%
      head(50) %>%
      select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count)
    
    write.csv(go_top, 
              file.path(OUTPUT_DIR, paste0("WGCNA_GO_", mod, ".csv")),
              row.names = FALSE)
    
    cat("\n", mod, "- GO top 5:\n")
    print(head(go_top[, c("Description", "pvalue", "p.adjust")], 5))
  }
  
  # KEGG
  if (!is.null(res$KEGG) && nrow(res$KEGG@result) > 0) {
    kegg_top <- res$KEGG@result %>%
      arrange(pvalue) %>%
      head(50)
    
    write.csv(kegg_top,
              file.path(OUTPUT_DIR, paste0("WGCNA_KEGG_", mod, ".csv")),
              row.names = FALSE)
    
    cat("\n", mod, "- KEGG top 5:\n")
    print(head(kegg_top[, c("Description", "pvalue", "p.adjust")], 5))
  }
  
  # Reactome
  if (!is.null(res$Reactome) && nrow(res$Reactome@result) > 0) {
    react_top <- res$Reactome@result %>%
      arrange(pvalue) %>%
      head(50)
    
    write.csv(react_top,
              file.path(OUTPUT_DIR, paste0("WGCNA_Reactome_", mod, ".csv")),
              row.names = FALSE)
  }
}

###############################################################################
### Enrichment dotplots
###############################################################################

cat("\n---Creating enrichment plots ---\n")

# Simplify GO terms (remove redundancy)
prune_go_terms <- function(go_res, cutoff = 0.7) {
  if (is.null(go_res) || nrow(go_res@result) == 0) return(go_res)
  
  tryCatch({
    clusterProfiler::simplify(
      go_res,
      cutoff = cutoff,
      by = "qvalue",
      select_fun = min,
      measure = "Wang"
    )
  }, error = function(e) {
    message("GO pruning failed: ", e$message)
    go_res
  })
}

# Plotting function
plot_module_enrichment <- function(enrich_res, module_color, 
                                   padj_cutoff = 0.25, top_n = 10) {
  
  plots <- list()
  
  for (ont in c("GO", "KEGG", "Reactome")) {
    res <- enrich_res[[ont]]
    if (is.null(res) || nrow(res@result) == 0) next
    
    
    # Prune but limit GO terms beforehand
    if (ont == "GO") {
      if (nrow(res@result) > 150) {
        res@result <- res@result %>% arrange(qvalue) %>% slice_head(n = 150)
      }
      res <- prune_go_terms(res)
    }
    
    
    df <- res@result %>%
      filter(!is.na(qvalue) & qvalue < padj_cutoff)
    
    if (nrow(df) == 0) next
    
    # Top N terms
    df_top <- df %>%
      slice_min(order_by = qvalue, n = top_n)
    
    df_top$Description <- stringr::str_wrap(df_top$Description, width = 45)
    
    p <- ggplot(df_top, aes(
      x = -log10(qvalue),
      y = reorder(Description, -log10(qvalue)),
      size = Count,
      color = qvalue
    )) +
      geom_point(alpha = 0.9) +
      scale_color_viridis_c(option = "C", direction = -1, name = "q-value") +
      scale_size_continuous(range = c(2.5, 7), name = "Gene count") +
      theme_classic(base_size = 12) +
      theme(
        axis.text.y = element_text(size = 10),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 14)
      ) +
      labs(
        title = paste0(ont, " Enrichment: ", module_color, " module"),
        x = "-log10(q-value)",
        y = NULL
      )
    
    plots[[ont]] <- p
  }
  
  plots
}

# Generate plots
pdf(file.path(OUTPUT_DIR, "11i_WGCNA_module_enrichment.pdf"),
    width = 10, height = 8)

for (mod in sig_modules) {
  enrich_res <- enrichment_results[[mod]]
  if (is.null(enrich_res)) next
  
  module_plots <- plot_module_enrichment(enrich_res, mod, 
                                         padj_cutoff = 0.25, top_n = 10)
  
  for (p in module_plots) {
    if (!is.null(p)) print(p)
  }
}

dev.off()

cat("Saved: 11i_WGCNA_module_enrichment.pdf\n")

###############################################################################
### Hub genes vs DE results
###############################################################################

cat("\n---  Hub genes vs DE overlap ---\n")

# Check which DE results exist
de_available <- c()
if (exists("de_neuron")) de_available <- c(de_available, "neuron")
if (exists("de_oligo")) de_available <- c(de_available, "oligo")
if (exists("de_all")) de_available <- c(de_available, "all")

if (length(de_available) == 0) {
  cat("WARNING: No DE results found (de_neuron, de_oligo, de_all)\n")
  cat("Skipping DE overlap analysis\n")
} else {
  cat("DE results available:", paste(de_available, collapse = ", "), "\n")
  
  # Use neuron DE by default (since WGCNA was on neurons)
  de_use <- if("neuron" %in% de_available) de_neuron else 
    if("all" %in% de_available) de_all else de_oligo
  
  de_use_name <- if("neuron" %in% de_available) "neuron" else
    if("all" %in% de_available) "all" else "oligo"
  
  cat("Using DE results from:", de_use_name, "\n")
  
  # Convert DE rownames to have SYMBOL column
  de_df <- de_use %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene_id")
  
  # If DE uses ENSEMBL, convert to SYMBOL
  if (grepl("^ENSG", de_df$gene_id[1])) {
    de_symbols <- mapIds(org.Hs.eg.db,
                         keys = de_df$gene_id,
                         column = "SYMBOL",
                         keytype = "ENSEMBL",
                         multiVals = "first")
    de_df$SYMBOL <- de_symbols
  } else {
    de_df$SYMBOL <- de_df$gene_id
  }
  
  de_df <- de_df %>% filter(!is.na(SYMBOL))
  
  # Hub gene overlap
  hub_gene_overlap <- data.frame()
  
  for (mod in sig_modules) {
    hub_file <- file.path(OUTPUT_DIR, paste0("WGCNA_hub_genes_", mod, ".csv"))
    
    if (!file.exists(hub_file)) {
      cat("  ", mod, ": hub file not found\n")
      next
    }
    
    hub_genes <- read.csv(hub_file)
    top_hubs <- head(hub_genes, 20)
    
    # Match by SYMBOL
    top_hubs <- top_hubs %>%
      left_join(de_df %>% select(SYMBOL, logFC, adj.P.Val, t),
                by = "SYMBOL",
                suffix = c("", "_DE"))
    
    top_hubs$in_DE_sig <- !is.na(top_hubs$adj.P.Val) & top_hubs$adj.P.Val < 0.05
    top_hubs$module <- mod
    
    hub_gene_overlap <- rbind(hub_gene_overlap, top_hubs)
  }
  
  # Save
  write.csv(hub_gene_overlap,
            file.path(OUTPUT_DIR, "WGCNA_hub_genes_DE_overlap.csv"),
            row.names = FALSE)
  
  # Summary
  cat("\nHub genes that are also DE (adj.P < 0.05):\n")
  overlap_summary <- hub_gene_overlap %>%
    filter(in_DE_sig) %>%
    select(module, SYMBOL, MM, logFC, adj.P.Val) %>%
    arrange(module, adj.P.Val)
  
  print(overlap_summary)
  cat("\nTotal overlapping:", nrow(overlap_summary), "\n")
}

###############################################################################
### Known Neurodegenerative genes in modules
###############################################################################

cat("\n---  Known Neurodegenerative genes ---\n")

pd_genes <- c(
  "SNCA", "PARK7", "LRRK2", "PINK1", "PRKN", "GBA", "VPS35",
  "MAPT", "UCHL1", "ATP13A2", "PLA2G6", "FBXO7", "DNAJC6",
  "SYNJ1", "GCH1", "DJ-1", "APOE", "APP", "PSEN1", "PSEN2", "TREM2", "CLU", "PICALM", 
  "HTT", "TARDBP", "C9orf72", "FUS", "SOD1", "COQ2", "ALS2", "PRNP" 
)

pd_genes_present <- pd_genes[pd_genes %in% gene_modules$SYMBOL]

cat("Known PD genes in network:", length(pd_genes_present), "/", 
    length(pd_genes), "\n")
if (length(pd_genes_present) > 0) {
  cat("  ", paste(pd_genes_present, collapse = ", "), "\n")
}

# Module assignments
if (length(pd_genes_present) > 0) {
  pd_gene_modules <- gene_modules %>%
    filter(SYMBOL %in% pd_genes_present) %>%
    arrange(module, SYMBOL)
  
  cat("\nPD gene module assignments:\n")
  print(pd_gene_modules)
  
  write.csv(pd_gene_modules,
            file.path(OUTPUT_DIR, "WGCNA_PD_genes_modules.csv"),
            row.names = FALSE)
  
  # Check if any are hubs
  pd_hubs <- data.frame()
  
  for (mod in unique(pd_gene_modules$module)) {
    if (mod == "grey") next
    
    hub_file <- file.path(OUTPUT_DIR, paste0("WGCNA_hub_genes_", mod, ".csv"))
    
    if (file.exists(hub_file)) {
      hub_genes <- read.csv(hub_file)
      
      pd_in_hubs <- hub_genes %>%
        filter(SYMBOL %in% pd_genes_present) %>%
        mutate(is_hub = MM > 0.7)
      
      if (nrow(pd_in_hubs) > 0) {
        pd_hubs <- rbind(pd_hubs, pd_in_hubs)
      }
    }
  }
  
  if (nrow(pd_hubs) > 0) {
    cat("\nPD genes that are hubs (MM > 0.7):\n")
    print(pd_hubs)
  } else {
    cat("\nNo PD genes are hubs\n")
  }
}

###############################################################################
### Hub gene expression plots
###############################################################################

cat("\n--- Hub gene expression ---\n")

plot_hub_expression <- function(module_color, n_top = 6, spe_obj) {
  
  cat("  Plotting", module_color, "module\n")
  
  # Load hub genes
  hub_file <- file.path(OUTPUT_DIR, paste0("WGCNA_hub_genes_", module_color, ".csv"))
  if (!file.exists(hub_file)) {
    cat("    Hub file not found\n")
    return(NULL)
  }
  
  hub_data <- read.csv(hub_file)
  top_symbols <- head(hub_data$SYMBOL, n_top)
  
  # Get expression
  if ("logcounts" %in% assayNames(spe_obj)) {
    logmat <- assay(spe_obj, "logcounts")
  } else {
    logmat <- cpm(assay(spe_obj, "counts"), log = TRUE, prior.count = 1)
  }
  
  # Match genes (handle ENSEMBL vs SYMBOL)
  if (grepl("^ENSG", rownames(logmat)[1])) {
    # Convert ENSEMBL to SYMBOL
    gene_map <- mapIds(org.Hs.eg.db,
                       keys = rownames(logmat),
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")
    
    # Find indices
    symbol_to_ensembl <- names(gene_map)[match(top_symbols, gene_map)]
    genes_present <- symbol_to_ensembl[!is.na(symbol_to_ensembl)]
  } else {
    genes_present <- top_symbols[top_symbols %in% rownames(logmat)]
  }
  
  if (length(genes_present) == 0) {
    cat("    No hub genes found in expression matrix\n")
    return(NULL)
  }
  
  # Extract expression
  expr_df <- as.data.frame(t(logmat[genes_present, , drop = FALSE]))
  expr_df$ROI <- rownames(expr_df)
  
  # Add metadata
  meta_df <- as.data.frame(colData(spe_obj))
  meta_df$ROI <- rownames(meta_df)
  
  expr_long <- expr_df %>%
    left_join(meta_df, by = "ROI") %>%
    pivot_longer(cols = all_of(genes_present),
                 names_to = "gene_id",
                 values_to = "logcounts")
  
  # Convert gene IDs to symbols for plot
  if (grepl("^ENSG", expr_long$gene_id[1])) {
    gene_map <- mapIds(org.Hs.eg.db,
                       keys = expr_long$gene_id,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")
    expr_long$gene <- gene_map[expr_long$gene_id]
  } else {
    expr_long$gene <- expr_long$gene_id
  }
  
  # Plot
  p <- ggplot(expr_long, aes(x = Group, y = logcounts, color = Group)) +
    geom_violin(alpha = 0.3, fill = "grey90", color = NA) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
    stat_summary(fun = mean, geom = "point", size = 3, shape = 18, color = "black") +
    facet_wrap(~gene, scales = "free_y", ncol = 3) +
    scale_color_manual(values = c("PD" = "#4575B4", "MSA" = "#D73027")) +
    labs(title = paste0("Hub genes: ", module_color, " module"),
         x = NULL, y = "log2 CPM") +
    theme_classic(base_size = 11) +
    theme(
      strip.text = element_text(face = "bold.italic"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
  
  ggsave(file.path(OUTPUT_DIR, paste0("11j_WGCNA_", module_color, "_hub_expression.pdf")),
         p, width = 10, height = 8)
  
  cat("    Saved hub expression plot\n")
  
  return(p)
}

# Plot for each significant module
for (mod in sig_modules) {
  if (mod == "grey") next
  plot_hub_expression(mod, n_top = 6, spe_obj = spe_qc)
}
###############################################################################
### Summary
###############################################################################

cat("\n=== WGCNA Follow-up Complete ===\n")
cat("Modules analyzed:", length(sig_modules), "\n")
if (exists("overlap_summary")) {
  cat("Hub-DE overlap:", nrow(overlap_summary), "genes\n")
}
cat("Known PD genes:", length(pd_genes_present), "\n")
cat("\nFiles created:\n")
cat("  - WGCNA_GO_*.csv (pathway enrichment)\n")
cat("  - WGCNA_KEGG_*.csv (KEGG enrichment)\n")
cat("  - 11i_WGCNA_module_enrichment.pdf (dotplots)\n")
cat("  - WGCNA_hub_genes_DE_overlap.csv\n")
cat("  - WGCNA_PD_genes_modules.csv\n")
cat("  - 11j_WGCNA_*_hub_expression.pdf (per module)\n")
cat("  - 11k_WGCNA_network_*.pdf (per module)\n")
