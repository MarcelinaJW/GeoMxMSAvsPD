###GeoMx Spatial Transcriptomics Analysis
#Project: Transcrptomic changes in the SN of PD and MSA patients 2026
#Tissue: Substantia nigra (block 15 FFPE) PDUK Tissue bank
#Platform: NanoString GeoMx Digital Spatial Profiler (now Bruker)
#Panel: Whole Transcriptome Atlas (WTA)
#Experimental groups: MSA (n=11) and PD (n=10)
#Technical factors: SlideName, Batch
#Biological factors: Patient, Group, Cell

#Loading packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# The following initializes most up to date version of Bioc
BiocManager::install()

suppressPackageStartupMessages({
  library(NanoStringNCTools)
  library(GeomxTools)
  library(GeoMxWorkflows)
  library(standR)
  library(SpatialExperiment)
  library(readxl)
  library(tidyverse)
  library(here)
  library(scater)
  library(edgeR)
  library(limma)
  library(ggpubr)
  library(patchwork)
  library(pheatmap)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
  library(DOSE)
  library(cluster)
  library(ggVennDiagram)
  library(pheatmap)
  library(scales)
  library(WGCNA)
  library(ReactomePA)
  library(stringr)
  library(enrichplot)
  library(GOSemSim)
  library(grid)
  library(gridExtra)
  library(ggplot2)
  library(cowplot)
  library(ggrepel)
  library(ggvenn)
})
set.seed(1234)

###Configuration----
datadir <- file.path("/Users/mw4217/Desktop/Geomx_data/R")
DCC_DIR <- file.path(datadir, "dccs")
PKC_FILE <- file.path(datadir, "pkcs", "Hs_R_NGS_WTA_v1.0.pkc")
LAB_WORKSHEET <- dir(file.path(datadir, "extras"), pattern = ".txt$",
                     full.names = TRUE, recursive = TRUE)
ANNOT_FILE <-   dir(file.path(datadir, "annotationv2"), pattern = ".xlsx$",
                    full.names = TRUE, recursive = TRUE)
ANNOT_FILE <- ANNOT_FILE[1]
OUTPUT_DIR <- "/Users/mw4217/Desktop/Geomx_data/R/results"
SN_REF_FILE <- file.path(datadir, "results", "spatialdecon", "E-HCAD-25_reference_matrix.rds")
TRAIT_FILE <-file.path(datadir, "metadata", "geomx_metadata.xlsx")
SEQ_QC_FILES <- dir(file.path(datadir, "seqqc"), pattern = "\\.txt$", full.names = TRUE, recursive = TRUE
)
# Keep only real text files (non-empty)
SEQ_QC_FILES <- SEQ_QC_FILES[file.info(SEQ_QC_FILES)$size > 0]

#Theme
theme_pub <- theme_classic(base_size = 11) +
  theme(axis.text = element_text(colour = "black"))

group_labels <- c(MSA_Neuron = "MSA Neuron",
                  MSA_Oligodendrocyte = "MSA Oligodendrocyte",
                  PD_Neuron = "PD Neuron",
                  PD_Oligodendrocyte = "PD Oligodendrocyte")

format_p <- function(p) {
  if (p < 0.001) return("< 0.001")
  as.character(signif(p, 3))
}
 
#Color palettes
cols_group <- c(PD = "#4C72B0", MSA = "#C44E52")
cols_batch <- c(Batch1 = "#DD8452", Batch2 = "#55A868")
cols_cell <- c(Neuron ="#E69F00", Oligodendrocyte = "#009E73")
cols_groupcell <- c(MSA_Neuron ="#009E73", MSA_Oligodendrocyte = "#6e00b2", PD_Neuron = "#D55E00", PD_Oligodendrocyte = "#b200b2")

###QC thresholds----------------------------------------------------------
#Based on NanoString standard parameters
QC <- list(
  #Probe QC
  minProbeRatio = 0.1,        #Flag probes with geomean ratio < 0.1 DO NOT CHANGE
  percentFailGrubbs = 20,     #Flag probes failing Grubbs in >= 20% segments DO NOT CHANGE
  #Segment QC (sequencing)
  min_aligned = 5000,      
  pct_trimmed_min = 80,       
  pct_stitched_min = 80,     
  pct_aligned_min = 75,      
  #Segment QC (tissue)
  min_area = 3000,       
  #Gene QC (LOQ-based) - using geometric mean method
  loq_sd = 1.75,   #SGR used 2, i did 1.5 before   nature paper does 1   #LOQ = NegGeoMean * NegGeoSD^loq_sd
  loq_min = 2,   #SGR used 2 
  gene_detection_pct = 0.05,    #SGR used 0.1, i did 0.05 before, nat paper does 0.1
  segment_detection_pct = 0.01  #SGR used 0.01, i also did 0.01 before, nat paper did 0.03
)

###Build annotation from LabWorksheet-----------------
lab_lines <- readLines(LAB_WORKSHEET)
annot_start <- which(grepl("^Annotations", lab_lines)) + 2
lab_data <- lab_lines[annot_start:length(lab_lines)] %>%
  .[. != ""] %>%
  map_dfr(function(line) {
    parts <- str_split(line, "\t")[[1]]
    if (!str_starts(parts[1], "DSP-")) return(NULL)
    is_ntc <- grepl("No Template Control", parts[2])
    tibble(
      Sample_ID = parts[1],
      SlideName = if (is_ntc) "NTC" else parts[2],
      Scan_Name = if (is_ntc) NA_character_ else parts[3],
      ROI_Label = if (is_ntc) NA_integer_ else as.integer(str_extract(parts[5], "\\d+")),
      Area = if (is_ntc) NA_real_ else as.numeric(parts[8]),
      is_NTC = is_ntc
    )
  })

cat("Parsed", nrow(lab_data), "samples from LabWorksheet\n")

#Load ROI annotations 
roi_annot <- read_excel(ANNOT_FILE) %>%
  dplyr::mutate(
    ROI_Label = as.integer(ROI_Label),
  ) %>%
  dplyr::select(Scan_Name, ROI_Label, Patient, Group, Cell, Sample_ID, Slide, Case)

#Merge and add batch information
annot_df <- lab_data %>%
  left_join(roi_annot, by = c("Scan_Name", "ROI_Label", "Sample_ID")) %>%
  dplyr::mutate(
    DCC_Filename = paste0(Sample_ID, ".dcc"),
    Patient = if_else(is_NTC, "NTC", Patient),
    Group = if_else(is_NTC, "NTC", Group),
    Cell = if_else(is_NTC, "NTC", Cell),
    Slide = if_else(is_NTC, "NTC", Slide),
    Case = if_else(is_NTC, "NTC", Case),
    #NTC slide name must contain "No Template Control" for GeomxTools to recognise it
    SlideName = if_else(is_NTC, "No Template Control", SlideName),
    Batch = case_when(
      SlideName %in% c("Slide1RepTHTPPP_271025", 
                       "Slide3THTPPP_271025",
                       "Slide4THTPPP_271025") ~ "Batch1",
      SlideName %in% c("Slide5THTPPP_261125", 
                       "Slide6THTPPP_261125",
                       "Slide7THTPPP_261125",
                       "Slide8THTPPP_261125") ~ "Batch2",
      TRUE ~ "NTC"
    )
  ) %>%
  dplyr::rename(`slide name` = SlideName) %>%
  dplyr::select(-is_NTC)

cat("\nAnnotation summary:\n")
print(table(annot_df$Batch, annot_df$Cell))
cat("\nBatch distribution:\n")
print(table(annot_df$Batch, annot_df$Group))


###Load sequencing QC-------------------------------------
parse_seq_qc <- function(filepath) {
  lines <- readLines(filepath)
  qc_lines <- lines[grepl("^File:", lines)]
  map_dfr(qc_lines, function(line) {
    file_match <- str_match(line, "File:\\s*([^,]+)")
    tibble(
      Sample_ID = trimws(file_match[1, 2]),
      Raw = as.numeric(str_match(line, "Raw:\\s*(\\d+)")[1, 2]),
      Trimmed = as.numeric(str_match(line, "Trimmed:\\s*(\\d+)")[1, 2]),
      Stitched = as.numeric(str_match(line, "Stitched:\\s*(\\d+)")[1, 2]),
      Aligned = as.numeric(str_match(line, "Aligned:\\s*(\\d+)")[1, 2]),
      Deduplicated = as.numeric(str_match(line, "Deduplicated:\\s*(\\d+)")[1, 2])
    )
  })
}

seq_qc <- map_dfr(SEQ_QC_FILES, parse_seq_qc) %>%
  dplyr::mutate(pct_trimmed = Trimmed / Raw * 100,
                pct_stitched = Stitched / Trimmed * 100,
                pct_aligned = Aligned / Stitched * 100)

cat("Sequencing QC:", nrow(seq_qc), "samples, aligned range:",
    comma(min(seq_qc$Aligned, na.rm = TRUE)), "-", comma(max(seq_qc$Aligned, na.rm = TRUE)), "\n")

###Load GeoMx data------------------------------
dcc_files <- list.files(DCC_DIR, pattern = "\\.dcc$", full.names = TRUE, ignore.case = TRUE)
cat("\nFound", length(dcc_files), "DCC files\n")

geomx_data <- readNanoStringGeoMxSet(
  dccFiles = dcc_files,
  pkcFiles = PKC_FILE,
  phenoData = annot_df,
  phenoDataDccColName = "DCC_Filename",
  phenoDataColNames = c("Patient", "Group", "Cell", "slide name", "Area", "Batch", "Case"),
  protocolDataColNames = c("Sample_ID", "Scan_Name", "ROI_Label")
)

cat("Loaded:", ncol(geomx_data), "samples,", nrow(geomx_data), "probes\n")
cat("Feature type:", featureType(geomx_data), "\n")

######################################################
###Create Sankey Diagram-------------------------------------
################################################
library(dplyr)
library(ggplot2)
library(ggalluvial)
library(stringr)

# ---- Figure size 
width_mm  <- 183
height_mm <- 140
width_in  <- width_mm / 25.4
height_in <- height_mm / 25.4

# ---- Prepare data
pdata_df <- as.data.frame(pData(geomx_data))

count_mat <- pdata_df %>%
  count(Slide, Group, Case, Cell)

count_mat$Slide <- gsub("THTPPP_", "", count_mat$Slide)
count_mat$Slide <- gsub("Slide", "S", count_mat$Slide)
count_mat$Cell <- recode(count_mat$Cell,
                         "Oligodendrocyte" = "Oligo")

# Get unique cases in their original form
u <- unique(count_mat$Case)

# Extract prefix (letters) and number
prefix <- str_to_upper(str_extract(u, "^[A-Za-z]+"))   # "M", "PD"
num    <- as.numeric(str_extract(u, "\\d+"))           # 1, 2, 10

ord <- order(prefix, num, na.last = TRUE)
levels_case <- u[ord]

count_mat$Case <- factor(count_mat$Case, levels = levels_case)

# ---- Plot
sank <- ggplot(count_mat,
               aes(axis1 = Slide,
                   axis2 = Group,
                   axis3 = Case,
                   axis4 = Cell,
                   y = n)) +
  
  geom_alluvium(aes(fill = Group),
                width = 0.18,
                alpha = 0.8,
                knot.pos = 0.4) +
  
  geom_stratum(width = 0.18,
               fill = "grey92",
               colour = "grey60",
               size = 0.25) +
  
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)),
            size = 2.8) +
  
  scale_fill_manual(values = cols_group) +
  
  scale_x_discrete(limits = c("Slide","Group","Case", "Cell"),
                   expand = c(0.05, 0.05)) +
  
  scale_y_continuous(expand = c(0, 0)) +
  
  theme_minimal(base_size = 8) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    plot.margin = margin(6, 8, 6, 8)
  )

# ---- Export
pdf(file.path(OUTPUT_DIR, "01_Sankey_qc.pdf"),
    width = width_in,
    height = height_in,
    useDingbats = FALSE)

print(sank)
dev.off()

################################################
###Extract negative probes----------
######################################################
#This must be done before shiftCountsOne changes the counts
neg_idx <- which(fData(geomx_data)$CodeClass == "Negative")
neg_counts_orig <- exprs(geomx_data)[neg_idx, , drop = FALSE]

neg_geomean <- apply(neg_counts_orig, 2, function(x) exp(mean(log(pmax(x, 1)))))
neg_geosd   <- apply(neg_counts_orig, 2, function(x) exp(sd(log(pmax(x, 1)))))

cat("Negative probes:", length(neg_idx),
    "| Geomean range:", round(min(neg_geomean), 2), "-", round(max(neg_geomean), 2), "\n")

############################################################
#####Probe-level QC (GeomxTools standard)
##################################################################
cat("\n=== Probe QC (GeomxTools) ===\n")

geomx_data <- shiftCountsOne(geomx_data, useDALogic = TRUE)
geomx_data <- setBioProbeQCFlags(geomx_data,
                                 qcCutoffs = list(minProbeRatio = QC$minProbeRatio,
                                                  percentFailGrubbs = QC$percentFailGrubbs),
                                 removeLocalOutliers = TRUE)

qc_flags <- fData(geomx_data)[["QCFlags"]]
probes_exclude <- which(qc_flags$LowProbeRatio | qc_flags$GlobalGrubbsOutlier)
cat("Probes excluded:", length(probes_exclude), "of", nrow(geomx_data), "\n")

if (length(probes_exclude) > 0) geomx_data <- geomx_data[-probes_exclude, ]

geomx_target <- aggregateCounts(geomx_data)
cat("Targets after aggregation:", nrow(geomx_target), "\n")

##################################################################
###Convert to SpatialExperiment####################################
############################################################
## Gene-level count matrix (already aggregated)
count_mat <- exprs(geomx_target)
rownames(count_mat) <- make.unique(fData(geomx_target)$TargetName, sep = "_")

## Sample annotation — MUST come from geomx_target
sample_anno <- pData(geomx_target) %>% as.data.frame() %>%
  dplyr::mutate(SlideName = `slide name`,
                Sample_ID_clean = gsub("\\.dcc$", "", rownames(.)))

## Match sequencing QC
seq_match <- match(sample_anno$Sample_ID_clean, seq_qc$Sample_ID)
neg_match <- match(rownames(sample_anno), names(neg_geomean))
sample_anno <- sample_anno %>%
  dplyr::mutate(Raw = seq_qc$Raw[seq_match], Trimmed = seq_qc$Trimmed[seq_match],
                Stitched = seq_qc$Stitched[seq_match], Aligned = seq_qc$Aligned[seq_match],
                Deduplicated = seq_qc$Deduplicated[seq_match],
                pct_trimmed = seq_qc$pct_trimmed[seq_match],
                pct_stitched = seq_qc$pct_stitched[seq_match],
                pct_aligned = seq_qc$pct_aligned[seq_match],
                NegGeoMean = neg_geomean[neg_match],
                NegGeoSD = neg_geosd[neg_match])
rownames(sample_anno) <- colnames(count_mat)


spe <- SpatialExperiment(
  assays = list(counts = count_mat),
  colData = sample_anno,
  rowData = fData(geomx_target) %>% as.data.frame() %>% `rownames<-`(rownames(count_mat))
)
metadata(spe) <- list(NegProbes = neg_counts_orig, NegGeoMean = neg_geomean, NegGeoSD = neg_geosd)

cat("SpatialExperiment:", nrow(spe), "genes x", ncol(spe), "ROIs\n")

#################################################################################
### Segment QC
#################################################################################
cat("\n=== Segment QC ===\n")

#######Derived metrics####################################
spe$LibrarySize    <- colSums(assay(spe, "counts"))
spe$DetectionRate  <- colSums(assay(spe, "counts") > 1) / nrow(spe)

######QC flags##############################
spe$QC_Fail <- (!is.na(spe$Aligned)      & spe$Aligned      < QC$min_aligned)     |
  (!is.na(spe$pct_trimmed)  & spe$pct_trimmed  < QC$pct_trimmed_min) |
  (!is.na(spe$pct_stitched) & spe$pct_stitched < QC$pct_stitched_min)|
  (!is.na(spe$pct_aligned)  & spe$pct_aligned  < QC$pct_aligned_min) |
  (spe$Area < QC$min_area)

cat("Segments failed QC:", sum(spe$QC_Fail), "of", ncol(spe), "\n")

######Column data##################################################################
cd         <- as.data.frame(colData(spe))
cd$Batch   <- as.factor(cd$Batch)
cd$QC_Fail <- spe$QC_Fail   # carry flag into plot data

######Colour palettes ################################################
fill_grp <- c(
  PD  = "#4C72B0",MSA = "#C44E52", NTC = "#999999"    
)

cols_batch <- setNames(
  colorspace::qualitative_hcl(length(levels(cd$Batch)),
                              palette = "Dark 3"),
  levels(cd$Batch)
)

############QC threshold label helper####################################
vline_label <- function(x, label,
                        hjust = -0.08, vjust = -0.5,
                        size  = 3.2) {
  list(
    geom_vline(xintercept = x,
               linetype   = "dashed",
               linewidth  = 0.6,
               colour     = "grey30"),
    annotate("text",
             x      = x,
             y      = Inf,
             label  = label,
             hjust  = hjust,
             vjust  = vjust,
             size   = size,
             colour = "grey30",
             fontface = "italic")
  )
}

hline_label <- function(y, label,
                        hjust = -0.08, vjust = -0.5,
                        size  = 3.2) {
  list(
    geom_hline(yintercept = y,
               linetype   = "dashed",
               linewidth  = 0.6,
               colour     = "grey30"),
    annotate("text",
             x      = Inf,
             y      = y,
             label  = label,
             hjust  = hjust,
             vjust  = vjust,
             size   = size,
             colour = "grey30",
             fontface = "italic")
  )
}

######Master theme######################################################
theme_thesis <- theme_classic(base_size = 11, base_family = "sans") +
  theme(
    # Axes
    axis.title       = element_text(size = 10, colour = "black"),
    axis.text        = element_text(size = 9,  colour = "black"),
    axis.line        = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks       = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks.length = unit(2.5, "pt"),
    
    # Legend
    legend.title     = element_blank(),
    legend.text      = element_text(size = 9),
    legend.key.size  = unit(0.45, "cm"),
    legend.spacing.y = unit(0.15, "cm"),
    legend.background = element_blank(),
    legend.position  = "right",
    
    # Panel
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border     = element_blank(),
    
    # Strip (facets, if used)
    strip.background = element_blank(),
    strip.text       = element_text(size = 9, face = "bold"),
    
    # Panel tag (A, B, C …) — set by patchwork
    plot.tag         = element_text(size = 11, face = "bold", colour = "black"),
    plot.tag.position = "topleft",
    
    # Titles
    plot.title       = element_blank(),   # panel letter handled by patchwork tags
    plot.margin      = margin(8, 10, 4, 6)
  )

#######Shared y-axis label################################################
ylab_rois <- "Number of ROIs"

#######Panel A — Aligned reads################################################
p1 <- ggplot(cd, aes(Aligned, fill = Group)) +
  geom_histogram(bins = 35, alpha = 0.85, colour = "white", linewidth = 0.2) +
  vline_label(QC$min_aligned, paste0("min = ", scales::comma(QC$min_aligned))) +
  scale_x_log10(
    labels = scales::comma_format(),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  scale_fill_manual(values = fill_grp) +
  labs(subtitle = "Aligned reads",
       x = "Aligned reads (log\u2081\u2080)",
       y = ylab_rois) +
  theme_thesis

######Panel B — Alignment rate######################################################
p2 <- ggplot(cd, aes(pct_aligned, fill = Group)) +
  geom_histogram(bins = 35, alpha = 0.85, colour = "white", linewidth = 0.2) +
  vline_label(QC$pct_aligned_min, paste0("min = ", QC$pct_aligned_min, "%")) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  scale_fill_manual(values = fill_grp) +
  labs(subtitle = "Alignment rate",
       x = "Alignment (%)",
       y = ylab_rois) +
  theme_thesis

#######Panel C — Aligned reads by slide##########################################
p3 <- ggplot(cd, aes(x = reorder(SlideName, Aligned, FUN = median),
                     y = Aligned,
                     fill = Batch)) +
  geom_boxplot(outlier.shape  = 21,
               outlier.size   = 1.2,
               outlier.stroke = 0.3,
               outlier.alpha  = 0.6,
               linewidth      = 0.4,
               alpha          = 0.85) +
  hline_label(QC$min_aligned, paste0("min = ", scales::comma(QC$min_aligned)),
              hjust = 1.02, vjust = -0.5) +
  scale_y_log10(
    labels = scales::comma_format(),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  scale_fill_manual(values = cols_batch) +
  labs(subtitle = "Aligned reads by slide",
       x = NULL,
       y = "Aligned reads (log\u2081\u2080)") +
  theme_thesis +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 8))

######Panel D — ROI area ######################################################
p4 <- ggplot(cd, aes(Area, fill = Group)) +
  geom_histogram(bins = 35, alpha = 0.85, colour = "white", linewidth = 0.2) +
  vline_label(QC$min_area, paste0("min = ", scales::comma(QC$min_area))) +
  scale_x_continuous(
    labels = scales::comma_format(),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  scale_fill_manual(values = fill_grp) +
  labs(subtitle = expression("ROI area"),
       x = expression(paste("Area (", mu, m^2, ")")),
       y = ylab_rois) +
  theme_thesis

######Panel E — Negative probe signal##########################################
p5 <- ggplot(cd, aes(NegGeoMean, fill = Group)) +
  geom_histogram(bins = 35, alpha = 0.85, colour = "white", linewidth = 0.2) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  scale_fill_manual(values = fill_grp) +
  labs(subtitle = "Negative probe geometric mean",
       x = "Negative probe signal (geometric mean)",
       y = ylab_rois) +
  theme_thesis

######Panel F — Library size vs detection rate ####################################
p6 <- ggplot(cd, aes(LibrarySize, DetectionRate)) +
  geom_point(aes(colour = Group, shape = QC_Fail),
             size = 1.8, alpha = 0.75, stroke = 0.3) +
  scale_x_log10(
    labels = scales::comma_format(),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  scale_colour_manual(values = fill_grp) +
  scale_shape_manual(values  = c(`FALSE` = 16, `TRUE` = 4),
                     labels  = c(`FALSE` = "Pass", `TRUE` = "Fail"),
                     name    = "QC") +
  labs(subtitle = "Library size vs. detection rate",
       x = "Library size (log\u2081\u2080)",
       y = "Detection rate") +
  guides(
    colour = guide_legend(order = 1, override.aes = list(size = 3)),
    shape  = guide_legend(order = 2, override.aes = list(size = 3))
  ) +
  theme_thesis

######Assemble with patchwork################################################
final_plot <- (p1 | p2) /
  (p3 | p4) /
  (p5 | p6) +
  plot_layout(guides = "collect") +
  plot_annotation(
    tag_levels  = "A",                          
    title       = "Segment-level quality control",
    caption     = paste0(
      "Dashed lines indicate QC thresholds. ",
      "Crosses (×) in panel F denote segments failing \u22651 QC criterion. ",
      "n = ", ncol(spe), " ROIs total; ",
      sum(spe$QC_Fail), " failed QC."
    ),
    theme = theme(
      plot.title   = element_text(size = 13, face = "bold",
                                  margin = margin(b = 6)),
      plot.caption = element_text(size = 8, colour = "grey40",
                                  hjust = 0, margin = margin(t = 6)),
      plot.margin  = margin(10, 10, 10, 10)
    )
  ) &
  theme(legend.position = "right")

######Export##################################################################
pdf(file.path(OUTPUT_DIR, "02a_segment_qc_thesis.pdf"),
    width = 14, height = 12, useDingbats = FALSE)
print(final_plot)
dev.off()

ggsave(
  filename = file.path(OUTPUT_DIR, "02a_segment_qc_thesis.png"),
  plot     = final_plot,
  width    = 14, height = 12,
  dpi      = 600,
  bg       = "white"
)

######QC summary CSV ##################################################################
write.csv(
  as.data.frame(colData(spe)) %>%
    dplyr::select(Group, Patient, SlideName, Batch,
                  Area, Raw:pct_aligned,
                  LibrarySize, DetectionRate, QC_Fail),
  file.path(OUTPUT_DIR, "segment_qc_summary.csv"),
  row.names = FALSE
)

# Remove failed segments
spe <- spe[, !spe$QC_Fail]
metadata(spe)$NegProbes  <- metadata(spe)$NegProbes[, colnames(spe), drop = FALSE]
metadata(spe)$NegGeoMean <- metadata(spe)$NegGeoMean[colnames(spe)]
metadata(spe)$NegGeoSD   <- metadata(spe)$NegGeoSD[colnames(spe)]

cat("After segment QC:", ncol(spe), "ROIs\n")

################################################################################
###################PCA BIPLOT OF QC METRICS ##############################
################################################################################

library(ggplot2)
library(ggrepel)
library(dplyr)
library(grid)
library(patchwork)

###### QC variables######################################################
qc_vars <- c(
  "Aligned", "pct_aligned", "pct_trimmed", "pct_stitched",
  "Area", "LibrarySize", "DetectionRate", "NegGeoMean"
)

qc_var_labels <- c(
  Aligned       = "Aligned reads",
  pct_aligned   = "Alignment (%)",
  pct_trimmed   = "Trimmed (%)",
  pct_stitched  = "Stitched (%)",
  Area          = "ROI area",
  LibrarySize   = "Library size",
  DetectionRate = "Detection rate",
  NegGeoMean    = "Neg. probe signal"
)

######Matrix prep ##################################################################
qc_mat    <- cd[, qc_vars] %>% as.data.frame()
qc_scaled <- scale(qc_mat)
keep      <- complete.cases(qc_scaled)
qc_scaled <- qc_scaled[keep, ]

######PCA ##################################################################
pca           <- prcomp(qc_scaled, center = TRUE, scale. = TRUE)
var_explained <- summary(pca)$importance[2, 1:2] * 100   # % variance PC1, PC2

######Scores ########################################################################
scores        <- as.data.frame(pca$x[, 1:2])
scores$Batch  <- droplevels(as.factor(cd$Batch[keep]))
scores$Group  <- droplevels(as.factor(cd$Group[keep]))
scores$QC_Fail <- cd$QC_Fail[keep]

######Loadings ########################################################################
loadings          <- as.data.frame(pca$rotation[, 1:2])
loadings$Variable <- qc_var_labels[rownames(loadings)]   # human-readable names

####### Proportional arrow scaling
arrow_scale <- min(
  diff(range(scores$PC1)) / diff(range(loadings$PC1)),
  diff(range(scores$PC2)) / diff(range(loadings$PC2))
) * 0.65

loadings <- loadings %>%
  mutate(
    PC1_arrow = PC1 * arrow_scale,
    PC2_arrow = PC2 * arrow_scale
  )

######Loading contribution magnitude########################
loadings <- loadings %>%
  mutate(loading_mag = sqrt(PC1^2 + PC2^2))

###### Axis limits with padding ################################################
pad    <- 0.12    # 12% padding
x_rng  <- range(c(scores$PC1, loadings$PC1_arrow))
y_rng  <- range(c(scores$PC2, loadings$PC2_arrow))
x_pad  <- diff(x_rng) * pad
y_pad  <- diff(y_rng) * pad
x_lims <- x_rng + c(-x_pad, x_pad)
y_lims <- y_rng + c(-y_pad, y_pad)

######Theme####################################
theme_thesis_pca <- theme_classic(base_size = 11, base_family = "sans") +
  theme(
    axis.title        = element_text(size = 10, colour = "black"),
    axis.text         = element_text(size = 9,  colour = "black"),
    axis.line         = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks        = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks.length = unit(2.5, "pt"),
    
    panel.grid.major  = element_line(colour = "grey92", linewidth = 0.3),
    panel.grid.minor  = element_blank(),
    
    legend.title      = element_text(size = 9, face = "bold"),
    legend.text       = element_text(size = 9),
    legend.key.size   = unit(0.45, "cm"),
    legend.spacing.y  = unit(0.15, "cm"),
    legend.background = element_blank(),
    legend.box        = "vertical",      # stack Batch + Group legends
    legend.position   = "right",
    
    plot.title        = element_text(size = 13, face = "bold",
                                     margin = margin(b = 4)),
    plot.subtitle     = element_text(size = 9, colour = "grey40",
                                     margin = margin(b = 10)),
    plot.caption      = element_text(size = 8, colour = "grey40",
                                     hjust = 0, margin = margin(t = 8)),
    plot.margin       = margin(10, 10, 10, 10)
  )

###### Build plot##################################################################
p_pca_biplot <- ggplot(scores, aes(PC1, PC2)) +
  
  # Zero reference lines
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70", linetype = "solid") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey70", linetype = "solid") +
  
  # 95% confidence ellipses per Batch (norm = Gaussian)
  stat_ellipse(
    aes(colour = Batch),
    type      = "norm",
    linewidth = 0.55,
    linetype  = "dashed",
    alpha     = 0.55,
    show.legend = FALSE
  ) +
  
  # Points
  geom_point(
    aes(colour = Batch,
        shape  = Group,
        size   = QC_Fail),     
    alpha  = 0.82,
    stroke = 0.35
  ) +
  scale_size_manual(
    values = c(`FALSE` = 2.2, `TRUE` = 3.2),
    guide  = "none"          
  ) +
  
  # Overlay crosses on QC-fail points
  geom_point(
    data   = subset(scores, QC_Fail),
    aes(PC1, PC2),
    shape  = 4,               # cross ×
    size   = 3.5,
    stroke = 0.7,
    colour = "black",
    inherit.aes = FALSE
  ) +
  
  # Loading arrows — shaded by magnitude
  geom_segment(
    data = loadings,
    aes(x = 0, y = 0,
        xend   = PC1_arrow,
        yend   = PC2_arrow,
        alpha  = loading_mag),
    inherit.aes = FALSE,
    arrow     = arrow(length = unit(0.20, "cm"),
                      type   = "closed"),
    linewidth = 0.55,
    colour    = "grey15",
    show.legend = FALSE
  ) +
  scale_alpha_continuous(range = c(0.45, 1), guide = "none") +
  
  geom_text_repel(
    data          = loadings,
    aes(PC1_arrow, PC2_arrow, label = Variable),
    inherit.aes   = FALSE,
    size          = 3.0,
    fontface      = "italic",
    colour        = "grey15",
    segment.size  = 0.25,
    segment.colour = "grey50",
    box.padding   = 0.45,
    point.padding = 0.2,
    max.overlaps  = Inf,
    seed          = 42        
  ) +
  
  # Scales
  scale_colour_manual(values = cols_batch, name = "Batch") +
  scale_shape_manual(
    values = c(PD = 16, MSA = 17, NTC = 15),   # circle, triangle, square
    name   = "Group"
  ) +
  
  # Axis labels with variance explained
  labs(
    title    = "PCA of segment-level QC metrics",
    subtitle = paste0(
      "Arrows indicate variable loadings; dashed ellipses = 95% CI per batch.\n",
      "Crosses (×) denote segments failing \u22651 QC criterion."
    ),
    x = paste0("PC1 (", round(var_explained[1], 1), "% variance explained)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "% variance explained)"),
    caption = paste0(
      "PCA performed on ", ncol(qc_scaled), " QC variables (z-scored). ",
      "n = ", nrow(scores), " ROIs; ",
      sum(scores$QC_Fail), " failed QC."
    )
  ) +
  
  # Equal aspect ratio — essential for biplots
  coord_equal(xlim = x_lims, ylim = y_lims, clip = "off") +
  
  # Guides — clean stacked legends
  guides(
    colour = guide_legend(order = 1,
                          override.aes = list(size = 3, shape = 16,
                                              alpha = 1)),
    shape  = guide_legend(order = 2,
                          override.aes = list(size = 3, colour = "grey30",
                                              alpha = 1))
  ) +
  
  theme_thesis_pca

######Export ############################################################
pdf(file.path(OUTPUT_DIR, "02b_segment_qc_PCAbiplot_thesis.pdf"),
    width = 7.5, height = 6.5, useDingbats = FALSE)
print(p_pca_biplot)
dev.off()

ggsave(
  filename = file.path(OUTPUT_DIR, "02b_segment_qc_PCAbiplot_thesis.png"),
  plot     = p_pca_biplot,
  width    = 7.5, height = 6.5,
  dpi      = 600,
  bg       = "white"
)
##PC1 (X%)   — typically driven by library size, aligned reads, detection rate  
##PC2 (Y%)   — often driven by background, alignment %, area













################################################################################
####### LOQ-BASED GENE FILTERING##############################
################################################################################
cat("\n=== LOQ Filtering (Geometric Mean) ===\n")

######LOQ calculation##################################################################
loq     <- pmax(
  metadata(spe)$NegGeoMean * (metadata(spe)$NegGeoSD ^ QC$loq_sd),
  QC$loq_min
)
loq_mat   <- matrix(rep(loq, each = nrow(spe)), nrow = nrow(spe))
above_loq <- assay(spe, "counts") > loq_mat

###### Segment filtering ##################################################################
seg_loq_pct <- colSums(above_loq) / nrow(spe)
low_seg     <- seg_loq_pct < QC$segment_detection_pct
cat("Segments below LOQ threshold:", sum(low_seg), "\n")

spe       <- spe[, !low_seg]
above_loq <- above_loq[, !low_seg]

# Sync metadata
metadata(spe)$NegProbes  <- metadata(spe)$NegProbes[, colnames(spe), drop = FALSE]
metadata(spe)$NegGeoMean <- metadata(spe)$NegGeoMean[colnames(spe)]
metadata(spe)$NegGeoSD   <- metadata(spe)$NegGeoSD[colnames(spe)]

######Gene detection per ROI######################################################
gene_loq_pct <- rowSums(above_loq) / ncol(spe)

###### Group factor for plotting ################################################
spe$Group_PDMSA_qc <- factor(
  case_when(
    spe$Group == "MSA" & spe$Cell == "Neuron"          ~ "MSA_Neuron",
    spe$Group == "PD"  & spe$Cell == "Neuron"          ~ "PD_Neuron",
    spe$Group == "MSA" & spe$Cell == "Oligodendrocyte" ~ "MSA_Oligodendrocyte",
    spe$Group == "PD"  & spe$Cell == "Oligodendrocyte" ~ "PD_Oligodendrocyte"
  ),
  levels = c("MSA_Neuron", "MSA_Oligodendrocyte",
             "PD_Neuron",  "PD_Oligodendrocyte")
)

######Plot data ##################################################################─
genes_detected_per_roi <- colSums(above_loq)

qc_genes <- data.frame(
  GenesDetected = genes_detected_per_roi,
  Group_PDMSA   = spe$Group_PDMSA_qc,
  Patient       = spe$Patient
)

######Shared theme ######################################################
theme_thesis <- theme_classic(base_size = 11, base_family = "sans") +
  theme(
    axis.title        = element_text(size = 10, colour = "black"),
    axis.text         = element_text(size = 9,  colour = "black"),
    axis.line         = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks        = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks.length = unit(2.5, "pt"),
    
    panel.grid.major  = element_line(colour = "grey92", linewidth = 0.3),
    panel.grid.minor  = element_blank(),
    
    legend.title      = element_text(size = 9, face = "bold"),
    legend.text       = element_text(size = 9),
    legend.key.size   = unit(0.45, "cm"),
    legend.background = element_blank(),
    legend.position   = "right",
    
    strip.background  = element_blank(),
    strip.text        = element_text(size = 9, face = "bold"),
    
    plot.tag          = element_text(size = 11, face = "bold"),
    plot.tag.position = "topleft",
    plot.title        = element_text(size = 13, face = "bold",
                                     margin = margin(b = 4)),
    plot.subtitle     = element_text(size = 9,  colour = "grey40",
                                     margin = margin(b = 8)),
    plot.caption      = element_text(size = 8,  colour = "grey40",
                                     hjust = 0, margin = margin(t = 8)),
    plot.margin       = margin(10, 12, 8, 10)
  )

######Summary statistics for annotation ##############################
group_stats <- qc_genes %>%
  group_by(Group_PDMSA) %>%
  summarise(
    med    = median(GenesDetected, na.rm = TRUE),
    n      = n(),
    .groups = "drop"
  )

###### Panel A — Genes detected per ROI ##########################################
p_genes_roi <- ggplot(qc_genes,
                      aes(x = Group_PDMSA, y = GenesDetected)) +
  
  # Violin
  geom_violin(
    aes(fill = Group_PDMSA),
    alpha     = 0.55,
    scale     = "width",
    linewidth = 0.35,
    colour    = "grey30",
    trim      = TRUE
  ) +
  
  # Boxplot inset
  geom_boxplot(
    aes(fill = Group_PDMSA),
    width         = 0.14,
    alpha         = 0.90,
    linewidth     = 0.4,
    outlier.shape = NA,    # outliers shown via jitter below
    colour        = "grey20"
  ) +
  
  # Jittered individual ROIs, coloured by patient
  geom_jitter(
    aes(colour = Patient),
    width  = 0.12,
    size   = 1.6,
    alpha  = 0.75,
    stroke = 0.2
  ) +
  
  # Median n label beneath each group
  geom_text(
    data  = group_stats,
    aes(x = Group_PDMSA, y = -Inf,
        label = paste0("n = ", n)),
    vjust    = -0.5,
    size     = 2.8,
    colour   = "grey40",
    fontface = "italic",
    inherit.aes = FALSE
  ) +
  
  scale_fill_manual(values = cols_groupcell, guide = "none") +
  scale_colour_discrete(name = "Patient") +
  scale_x_discrete(labels = group_labels) +
  scale_y_continuous(
    expand = expansion(mult = c(0.08, 0.06)),
    labels = scales::comma_format()
  ) +
  
  labs(
    subtitle = paste0(
      "After segment LOQ filtering; ",
      ncol(spe), " ROIs retained"
    ),
    x = NULL,
    y = "Genes detected per ROI (above LOQ)"
  ) +
  
  guides(
    colour = guide_legend(
      title    = "Patient",
      ncol     = 1,
      override.aes = list(size = 3, alpha = 1)
    )
  ) +
  
  theme_thesis +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 9)
  )

######Panel B — LOQ gene filtering histogram##########################################

# Threshold value as percentage
thresh_pct <- QC$gene_detection_pct * 100

loq_df <- data.frame(
  detection_pct = gene_loq_pct * 100,    # convert to % for axis
  retained      = gene_loq_pct >= QC$gene_detection_pct
)

# Counts for subtitle
n_retained <- sum(loq_df$retained)
n_removed  <- sum(!loq_df$retained)

p_loq_hist <- ggplot(loq_df, aes(detection_pct, fill = retained)) +
  
  geom_histogram(
    bins      = 50,
    alpha     = 0.85,
    colour    = "white",
    linewidth = 0.2
  ) +
  
  # Threshold line
  geom_vline(
    xintercept = thresh_pct,
    linetype   = "dashed",
    linewidth  = 0.65,
    colour     = "grey30"
  ) +
  
  # Threshold label
  annotate(
    "text",
    x        = thresh_pct,
    y        = Inf,
    label    = paste0("Threshold\n", thresh_pct, "%"),
    hjust    = -0.08,
    vjust    = 1.4,
    size     = 3.0,
    colour   = "grey20",
    fontface = "italic"
  ) +
  
  
  scale_fill_manual(
    values = c("TRUE"  = "#009E73", 
               "FALSE" = "grey72"),
    labels = c("TRUE"  = "Retained",
               "FALSE" = "Removed"),
    name   = NULL
  ) +
  
  scale_x_continuous(
    expand = expansion(mult = c(0.01, 0.05)),
    labels = function(x) paste0(x, "%")
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.18)),   # headroom for annotations
    labels = scales::comma_format()
  ) +
  
  labs(
    subtitle = paste0(
      n_retained, " genes retained; ",
      n_removed,  " genes removed at ",
      thresh_pct, "% LOQ threshold"
    ),
    x = "ROIs with gene expression above LOQ (%)",
    y = "Number of genes"
  ) +
  
  theme_thesis +
  theme(legend.position = "top")

######Assemble panels with patchwork ##########################################
combined_loq <- (p_genes_roi | p_loq_hist) +
  plot_annotation(
    title   = "LOQ-based quality control and gene filtering",
    caption = paste0(
      "LOQ = NegGeoMean \u00d7 NegGeoSD^", QC$loq_sd,
      " (minimum = ", QC$loq_min, "). ",
      "Genes retained if detected above LOQ in \u2265",
      thresh_pct, "% of ROIs. ",
      "Segments retained if \u2265",
      QC$segment_detection_pct * 100, "% of genes above LOQ."
    ),
    tag_levels = "A",
    theme = theme(
      plot.title   = element_text(size = 13, face = "bold",
                                  margin = margin(b = 4)),
      plot.caption = element_text(size = 8, colour = "grey40",
                                  hjust = 0, margin = margin(t = 8)),
      plot.margin  = margin(10, 10, 10, 10)
    )
  )

######Export ##################################################################
pdf(file.path(OUTPUT_DIR, "03_loq_filtering.pdf"),
    width = 12, height = 5.5, useDingbats = FALSE)
print(combined_loq)
dev.off()

ggsave(
  filename = file.path(OUTPUT_DIR, "03_loq_filtering.png"),
  plot     = combined_loq,
  width    = 12, height = 5.5,
  dpi      = 600,
  bg       = "white"
)

######Receptor detection rates##########################################
receptor_check <- c("TH", "TPPP", "SNCA")
cat("Receptor detection rates (before gene filtering):\n")
for (g in receptor_check) {
  if (g %in% names(gene_loq_pct))
    cat("  ", g, ":", round(gene_loq_pct[g] * 100, 1), "% of ROIs above LOQ\n")
  else
    cat("  ", g, ": not in probe panel\n")
}

######Gene filtering ############################################################
raw_counts_prefilter <- assay(spe, "counts")
spe <- spe[gene_loq_pct >= QC$gene_detection_pct, ]

cat("After LOQ:", nrow(spe), "genes x", ncol(spe), "ROIs\n")
cat("Genes removed by LOQ:", sum(gene_loq_pct < QC$gene_detection_pct),
    "of", length(gene_loq_pct), "\n")

###### Clean up ########################################################################
spe$Group_PDMSA_qc <- NULL








####################################################################################
###Create analysis factors---------------------------------------------
####################################################################################
spe$Group_PDMSA <- case_when(
  spe$Group == "MSA" & spe$Cell == "Neuron" ~ "MSA_Neuron",
  spe$Group == "PD" & spe$Cell == "Neuron" ~ "PD_Neuron",
  spe$Group == "MSA" & spe$Cell == "Oligodendrocyte" ~ "MSA_Oligodendrocyte",
  spe$Group == "PD" & spe$Cell == "Oligodendrocyte" ~ "PD_Oligodendrocyte"
)
spe$Group_PDMSA <- factor(spe$Group_PDMSA, levels = c("MSA_Neuron", "PD_Neuron", "MSA_Oligodendrocyte", "PD_Oligodendrocyte"))

cat("\nGroup_PDMSA distribution:\n")
print(table(spe$Group_PDMSA))







################################################################################
# ROI QUALITY CONTROL
################################################################################
cat("\n=== ROI Quality Control: Substantia Nigra ===\n\n")

spe_working <- spe

################################################################################
# Marker scores
################################################################################

markers <- list(
  TH_Core   = c("TH", "SLC6A3", "DDC", "SLC18A2", "MAP2", "SNAP25",
                "TUBB3", "ELAVL2", "BSN", "KCNJ6", "FOXA2", "ATP1A3"),
  TPPP_Core = c("TPPP", "OLIG1", "OLIG2", "CNP", "SOX10", "PLP1",
                "QKI", "MYRF", "MOG", "CLDN11", "OPALIN", "TMEM144",
                "MOBP", "MAG"),
  Vascular  = c("VWF", "PECAM1", "CD34", "FLT1", "EGFL7", "ANPEP",
                "CLDN5", "ITIH2", "ITIH3", "SEMA3G"),
  Immune    = c("PTPRC", "AIF1", "CD68", "ALDH1L1", "AQPR4",
                "C1QC", "GFAP", "SERPINA3")
)

logcpm <- cpm(assay(spe_working, "counts"), log = TRUE, prior.count = 1)

th_genes    <- intersect(markers$TH_Core,   rownames(spe_working))
tppp_genes  <- intersect(markers$TPPP_Core, rownames(spe_working))
vasc_genes  <- intersect(markers$Vascular,  rownames(spe_working))
immune_genes <- intersect(markers$Immune,   rownames(spe_working))

th_scores     <- colMeans(logcpm[th_genes,    , drop = FALSE])
tppp_scores   <- colMeans(logcpm[tppp_genes,  , drop = FALSE])
vasc_scores   <- if (length(vasc_genes) > 0)
  colMeans(logcpm[vasc_genes, , drop = FALSE]) else
    rep(0, ncol(spe_working))
immune_scores <- colMeans(logcpm[immune_genes, , drop = FALSE])
contam_score  <- (vasc_scores + immune_scores) / 2

################################################################################
#Ratio of wanted to unwanted markers
################################################################################

neuron_mask <- spe_working$Cell == "Neuron"
oligo_mask  <- spe_working$Cell == "Oligodendrocyte"

purity_ratio <- rep(NA_real_, ncol(spe_working))
purity_ratio[neuron_mask] <- th_scores[neuron_mask] /
  (th_scores[neuron_mask] + tppp_scores[neuron_mask])
purity_ratio[oligo_mask]  <- tppp_scores[oligo_mask] /
  (th_scores[oligo_mask] + tppp_scores[oligo_mask])

################################################################################
# Thresholds & flags
################################################################################

PURITY_THRESHOLD <- 0.496
CONTAM_THRESHOLD <- quantile(contam_score, 0.97)

flag_low_purity  <- replace(purity_ratio < PURITY_THRESHOLD,
                            is.na(purity_ratio), FALSE)
flag_high_contam <- contam_score > CONTAM_THRESHOLD
qc_pass          <- !(flag_low_purity | flag_high_contam)

################################################################################
# Results data frame
################################################################################

qc_results <- data.frame(
  ROI              = colnames(spe_working),
  Manual_Label     = spe_working$Cell,
  Group            = spe_working$Group,
  Patient          = spe_working$Patient,
  TH_Score         = th_scores,
  TPPP_Score       = tppp_scores,
  Purity_Ratio     = purity_ratio,
  Contamination    = contam_score,
  Flag_Low_Purity  = flag_low_purity,
  Flag_High_Contam = flag_high_contam,
  QC_Pass          = qc_pass,
  stringsAsFactors = FALSE
)

# ── Console summary ───────────────────────────────────────────────────────────
cat("Overall:\n")
cat("  Total ROIs:", nrow(qc_results), "\n")
cat("  Pass QC:", sum(qc_pass),
    sprintf("(%.1f%%)\n", 100 * mean(qc_pass)))
cat("  Fail QC:", sum(!qc_pass),
    sprintf("(%.1f%%)\n\n", 100 * mean(!qc_pass)))
cat("Failure breakdown:\n")
cat("  Low purity (<50%):", sum(flag_low_purity), "\n")
cat("  High contamination (>P97):", sum(flag_high_contam), "\n")
cat("  Both:", sum(flag_low_purity & flag_high_contam), "\n\n")

################################################################################
# Shared theme
################################################################################

theme_thesis <- theme_classic(base_size = 11, base_family = "sans") +
  theme(
    axis.title        = element_text(size = 10, colour = "black"),
    axis.text         = element_text(size = 9,  colour = "black"),
    axis.line         = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks        = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks.length = unit(2.5, "pt"),
    
    panel.grid.major  = element_line(colour = "grey92", linewidth = 0.3),
    panel.grid.minor  = element_blank(),
    
    legend.title      = element_text(size = 9, face = "bold"),
    legend.text       = element_text(size = 9),
    legend.key.size   = unit(0.45, "cm"),
    legend.background = element_blank(),
    legend.position   = "right",
    
    strip.background  = element_blank(),
    strip.text        = element_text(size = 9, face = "bold",
                                     margin = margin(b = 4)),
    
    plot.tag          = element_text(size = 11, face = "bold"),
    plot.tag.position = "topleft",
    plot.subtitle     = element_text(size = 8.5, colour = "grey40",
                                     margin = margin(b = 6)),
    plot.margin       = margin(8, 10, 6, 8)
  )

pal_cell <- c(
  Neuron          = "#0072B2",   
  Oligodendrocyte = "#E69F00"   
)


# QC shape encoding
shape_qc <- c("TRUE" = 16, "FALSE" = 4)  

# Threshold line style helper
thresh_line_discrete <- function(yint, label,
                                 x_pos   = -Inf,
                                 y_nudge = 0.02,
                                 hjust   = 0,
                                 size    = 2.9) {
  list(
    geom_hline(yintercept = yint,
               linetype   = "dashed",
               linewidth  = 0.55,
               colour     = "grey30"),
    annotate("text",
             x        = x_pos,
             y        = yint + y_nudge,
             label    = label,
             size     = size,
             colour   = "grey25",
             fontface = "italic",
             hjust    = hjust)
  )
}

thresh_line_continuous <- function(yint, label,
                                   x_pos   = -Inf,
                                   y_nudge = 0.02,
                                   hjust   = 0,
                                   size    = 2.9) {
  list(
    geom_hline(yintercept = yint,
               linetype   = "dashed",
               linewidth  = 0.55,
               colour     = "grey30"),
    annotate("text",
             x        = x_pos,
             y        = yint + y_nudge,
             label    = label,
             size     = size,
             colour   = "grey25",
             fontface = "italic",
             hjust    = hjust)
  )
}

################################################################################
# Panel A — TH vs TPPP marker scatter
################################################################################

# Density layer to show mass of points
p_a <- ggplot(qc_results, aes(TH_Score, TPPP_Score)) +
  
  
  # Identity line (equal marker expression)
  geom_abline(slope     = 1,
              intercept = 0,
              linetype  = "dashed",
              linewidth = 0.45,
              colour    = "grey60") +
  
  # Points
  geom_point(
    aes(colour = Manual_Label,
        shape  = QC_Pass),
    size   = 2.2,
    alpha  = 0.80,
    stroke = 0.35
  ) +
  
  # Region labels
  annotate("text",
           x = quantile(qc_results$TH_Score,   0.85),
           y = quantile(qc_results$TPPP_Score,  0.25),
           label    = "Neuron-enriched",
           colour   = pal_cell["Neuron"],
           size     = 3.0,
           fontface = "italic") +
  annotate("text",
           x = quantile(qc_results$TH_Score,   0.25),
           y = quantile(qc_results$TPPP_Score,  0.85),
           label    = "Oligodendrocyte-enriched",
           colour   = pal_cell["Oligodendrocyte"],
           size     = 3.0,
           fontface = "italic") +
  
  scale_colour_manual(values = pal_cell, name = "Cell type") +
  scale_shape_manual(
    values = shape_qc,
    labels = c("TRUE" = "Pass", "FALSE" = "Fail"),
    name   = "QC status"
  ) +
  
  labs(
    subtitle = paste0(
      "Dashed line = equal marker expression. ",
      "Contours show point density per cell type."
    ),
    x = "TH marker score (log\u2082 CPM)",
    y = "TPPP marker score (log\u2082 CPM)"
  ) +
  
  guides(
    colour = guide_legend(order = 1,
                          override.aes = list(size = 3, alpha = 1)),
    shape  = guide_legend(order = 2,
                          override.aes = list(size = 3, colour = "grey30"))
  ) +
  
  theme_thesis

################################################################################
# Panel B — Purity ratio
################################################################################
# Per-group per-cell-type median labels
purity_stats <- qc_results %>%
  group_by(Manual_Label, Group) %>%
  summarise(med = median(Purity_Ratio, na.rm = TRUE),
            n   = n(),
            .groups = "drop")

p_b <- ggplot(qc_results,
              aes(x = Manual_Label, y = Purity_Ratio)) +
  
  # ← use discrete-safe helper; x_pos = -Inf pins label to left edge
  thresh_line_discrete(
    PURITY_THRESHOLD,
    paste0("Threshold = ", round(PURITY_THRESHOLD * 100, 0), "%"),
    x_pos   = -Inf,
    y_nudge = 0.025,
    hjust   = -0.05
  ) +
  
  geom_violin(
    aes(fill = Manual_Label),
    alpha     = 0.50,
    trim      = TRUE,
    linewidth = 0.35,
    colour    = "grey30",
    scale     = "width"
  ) +
  
  geom_boxplot(
    aes(fill = Manual_Label),
    width         = 0.14,
    alpha         = 0.88,
    linewidth     = 0.40,
    outlier.shape = NA,
    colour        = "grey20"
  ) +
  
  geom_jitter(
    aes(shape = QC_Pass),
    colour = "grey25",
    width  = 0.12,
    size   = 1.5,
    alpha  = 0.65,
    stroke = 0.25
  ) +
  
  geom_text(
    data  = purity_stats,
    aes(x = Manual_Label,
        y = med,
        label = paste0(round(med * 100, 0), "%")),
    size        = 2.8,
    colour      = "white",
    fontface    = "bold",
    inherit.aes = FALSE
  ) +
  
  geom_text(
    data = purity_stats,
    aes(x = Manual_Label, y = -Inf,
        label = paste0("n = ", n)),
    vjust       = -0.5,
    size        = 2.7,
    colour      = "grey45",
    fontface    = "italic",
    inherit.aes = FALSE
  ) +
  
  facet_wrap(~ Group) +
  
  scale_fill_manual(values = pal_cell, guide = "none") +
  scale_shape_manual(
    values = shape_qc,
    labels = c("TRUE" = "Pass", "FALSE" = "Fail"),
    name   = "QC status",
    guide  = "none"
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0.10, 0.10))
  ) +
  
  labs(
    subtitle = paste0(
      "Neurons: ", round(mean(purity_ratio[neuron_mask], na.rm = TRUE) * 100, 0),
      "% TH | Oligodendrocytes: ",
      round(mean(purity_ratio[oligo_mask],  na.rm = TRUE) * 100, 0),
      "% TPPP (biologically realistic for SN)"
    ),
    x = NULL,
    y = "Purity ratio (target marker / total)"
  ) +
  
  theme_thesis +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

################################################################################
# Panel C — Contamination marker score
################################################################################
contam_stats <- qc_results %>%
  group_by(Manual_Label, Group) %>%
  summarise(med = median(Contamination, na.rm = TRUE),
            n   = n(),
            .groups = "drop")

p_c <- ggplot(qc_results,
              aes(x = Manual_Label, y = Contamination)) +
  
  thresh_line_discrete(
    CONTAM_THRESHOLD,
    paste0("97th percentile = ", round(CONTAM_THRESHOLD, 2)),
    x_pos   = -Inf,
    y_nudge = 0.03,
    hjust   = -0.05
  ) +
  
  geom_violin(
    aes(fill = Manual_Label),
    alpha     = 0.50,
    trim      = TRUE,
    linewidth = 0.35,
    colour    = "grey30",
    scale     = "width"
  ) +
  
  geom_boxplot(
    aes(fill = Manual_Label),
    width         = 0.14,
    alpha         = 0.88,
    linewidth     = 0.40,
    outlier.shape = NA,
    colour        = "grey20"
  ) +
  
  geom_jitter(
    aes(shape = QC_Pass),
    colour = "grey25",
    width  = 0.12,
    size   = 1.5,
    alpha  = 0.65,
    stroke = 0.25
  ) +
  
  geom_text(
    data = contam_stats,
    aes(x = Manual_Label, y = -Inf,
        label = paste0("n = ", n)),
    vjust       = -0.5,
    size        = 2.7,
    colour      = "grey45",
    fontface    = "italic",
    inherit.aes = FALSE
  ) +
  
  facet_wrap(~ Group) +
  
  scale_fill_manual(values = pal_cell, guide = "none") +
  scale_shape_manual(
    values = shape_qc,
    labels = c("TRUE" = "Pass", "FALSE" = "Fail"),
    name   = "QC status",
    guide  = "none"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.10, 0.12))
  ) +
  
  labs(
    subtitle = "Only extreme outliers (>97th percentile) flagged for removal",
    x = NULL,
    y = expression("Contamination score (mean log"[2]*" CPM)")
  ) +
  
  theme_thesis +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

################################################################################
# Panel D — Results by patient
################################################################################

p_d <- ggplot(qc_results,
              aes(x = Patient, y = Purity_Ratio,
                  colour = Manual_Label)) +
  
  geom_hline(
    yintercept = PURITY_THRESHOLD,
    linetype   = "dashed",
    linewidth  = 0.50,
    colour     = "grey30"
  ) +
  
  # Per-patient median bar
  stat_summary(
    aes(group = interaction(Patient, Manual_Label)),
    fun      = median,
    geom     = "crossbar",
    width    = 0.45,
    linewidth = 0.35,
    fatten   = 1.5,
    alpha    = 0.55
  ) +
  
  geom_jitter(
    aes(shape = QC_Pass),
    width  = 0.18,
    size   = 1.8,
    alpha  = 0.78,
    stroke = 0.30
  ) +
  
  facet_wrap(~ Group, scales = "free_x") +
  
  scale_colour_manual(values = pal_cell, name = "Cell type") +
  scale_shape_manual(
    values = shape_qc,
    labels = c("TRUE" = "Pass", "FALSE" = "Fail"),
    name   = "QC status"
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0.06, 0.08))
  ) +
  
  labs(
    subtitle = "Crossbars show per-patient median purity",
    x = "Patient ID",
    y = "Purity ratio"
  ) +
  
  guides(
    colour = guide_legend(order = 1,
                          override.aes = list(size = 3, alpha = 1)),
    shape  = guide_legend(order = 2,
                          override.aes = list(size = 3, colour = "grey30"))
  ) +
  
  theme_thesis +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 8))

################################################################################
# Assemble & export
################################################################################

final_roi_qc <- (p_a | p_b) / (p_c | p_d) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title   = "ROI quality control in the substantia nigra",
    caption = paste0(
      "Purity ratio: target marker score / (target + off-target marker score). ",
      "Contamination: mean log\u2082 CPM of vascular (n = ", length(vasc_genes),
      ") and immune (n = ", length(immune_genes), ") marker genes. ",
      "Thresholds: purity > ", round(PURITY_THRESHOLD * 100, 0),
      "%, contamination < 97th percentile. ",
      "Pass: ", sum(qc_pass), "/", length(qc_pass),
      " ROIs (", round(100 * mean(qc_pass), 1), "%)."
    ),
    tag_levels = "A",
    theme = theme(
      plot.title   = element_text(size = 13, face = "bold",
                                  margin = margin(b = 3)),
      plot.caption = element_text(size = 8,  colour = "grey40",
                                  hjust = 0, margin = margin(t = 8)),
      plot.margin  = margin(10, 10, 10, 10)
    )
  ) &
  theme(legend.position = "right")

pdf(file.path(OUTPUT_DIR, "03a_ROI_QC_substantia_nigra.pdf"),
    width = 14, height = 10, useDingbats = FALSE)
print(final_roi_qc)
dev.off()

ggsave(
  filename = file.path(OUTPUT_DIR, "03_ROI_QC_substantia_nigra.png"),
  plot     = final_roi_qc,
  width    = 14, height = 10,
  dpi      = 600,
  bg       = "white"
)

################################################################################
# Apply filter & save
################################################################################

write.csv(qc_results,
          file.path(OUTPUT_DIR, "ROI_QC_results.csv"),
          row.names = FALSE)

if (sum(!qc_pass) > 0) {
  write.csv(qc_results[!qc_pass, ],
            file.path(OUTPUT_DIR, "ROI_QC_failed.csv"),
            row.names = FALSE)
  spe_qc <- spe_working[, qc_pass]
} else {
  spe_qc <- spe_working
}

metadata(spe_qc)$NegProbes  <- metadata(spe_working)$NegProbes[, colnames(spe_qc), drop = FALSE]
metadata(spe_qc)$NegGeoMean <- metadata(spe_working)$NegGeoMean[colnames(spe_qc)]
metadata(spe_qc)$NegGeoSD   <- metadata(spe_working)$NegGeoSD[colnames(spe_qc)]

saveRDS(spe_qc, file.path(OUTPUT_DIR, "spe_qc.rds"))

cat("\n=== QC Complete ===\n")
cat("spe_qc:", ncol(spe_qc), "ROIs after QC\n")
cat("Removed:", ncol(spe_working) - ncol(spe_qc), "ROIs\n")







########################################################################
###Normalisation (TMM)################################################
#####################################################################
cat("\n=== TMM Normalisation ===\n")

dge_norm <- calcNormFactors(DGEList(counts = assay(spe_qc, "counts")), method = "TMM")
assay(spe_qc, "logcounts") <- cpm(dge_norm, log = TRUE, prior.count = 1)

#Housekeeping gene stability check, use counts not logcounts for cv or do sd on logcounts
hk_genes <- c("ACTB", "GAPDH", "TUBB", "RPL19", "RPLP0", "B2M")
hk_present <- intersect(hk_genes, rownames(spe))
if (length(hk_present) > 0) {
  hk_cv <- apply(assay(spe_qc, "logcounts")[hk_present, , drop = FALSE], 1, sd, )
                 #function(x) sd(x)/mean(x))
  cat("Housekeeping gene CV (lower = more stable):\n")
  print(round(sort(hk_cv), 3))
}

# PCA before batch correction (pre-RUV4)
spe_pre <- runPCA(spe_qc, exprs_values = "logcounts", ncomponents = 10)
pca_pre <- reducedDim(spe_pre, "PCA")
pca_pre_var <- attr(pca_pre, "percentVar")
if (is.null(pca_pre_var)) {
  pca_pre_prcomp <- prcomp(t(assay(spe, "logcounts")), center = TRUE, scale. = FALSE)
  pca_pre_var <- (pca_pre_prcomp$sdev^2 / sum(pca_pre_prcomp$sdev^2)) * 100
}

pca_pre_df <- data.frame(
  PC1 = pca_pre[, 1], PC2 = pca_pre[, 2],
  Group = spe_qc$Group_PDMSA, Batch = spe_qc$SlideName
)

p_pca_pre <- ggplot(pca_pre_df, aes(PC1, PC2, colour = Group)) +
  geom_point(size = 2.5, alpha = 0.7) +
  stat_ellipse(level = 0.95, linewidth = 0.6, linetype = "dashed", show.legend = FALSE) +
  scale_colour_manual(values = cols_groupcell, labels = group_labels) +
  labs(x = paste0("PC1 (", round(pca_pre_var[1], 1), "%)"),
       y = paste0("PC2 (", round(pca_pre_var[2], 1), "%)"),
       title = "Before batch correction",
       colour = NULL) +
  theme_pub + theme(legend.position = "top")

ggsave(file.path(OUTPUT_DIR, "02c_pca_pre_ruv4.pdf"), p_pca_pre, width = 6, height = 5)













###Batch correction (RUV4) with k selection-sgr used 500 
cat("\n=== Batch Correction (RUV4) ===\n")

spe_qc <- findNCGs(spe_qc, batch_name = "SlideName", top_n = 750)
cat("Found", length(metadata(spe_qc)$NCGs), "negative control genes\n")

#Evaluate k values
ruv_eval <- data.frame(k = 1:5, sil_group = NA, sil_batch = NA)

for (i in 1:5) {
  spe_k <- geomxBatchCorrection(spe_qc, factors = "Group_PDMSA", NCGs = metadata(spe_qc)$NCGs, k = i)
  spe_k <- runPCA(spe_k, exprs_values = "logcounts", ncomponents = 10)
  pca_coords <- reducedDim(spe_k, "PCA")[, 1:5]
  
  ruv_eval$sil_group[i] <- mean(silhouette(as.numeric(spe_k$Group_PDMSA), dist(pca_coords))[, 3])
  ruv_eval$sil_batch[i] <- mean(silhouette(as.numeric(factor(spe_k$SlideName)), dist(pca_coords))[, 3])
  
  cat("  k =", i, "| Group sil:", round(ruv_eval$sil_group[i], 3),
      "| Batch sil:", round(ruv_eval$sil_batch[i], 3), "\n")
}

# Select k by elbow: find where batch silhouette stops decreasing
# The elbow is the last k before the curve flattens or reverses
ruv_eval$batch_drop <- c(NA, -diff(ruv_eval$sil_batch))

# Elbow = last k with meaningful batch reduction before plateau/reversal
# A drop is meaningful if batch sil still decreases; reversal signals overcorrection
for (ek in 2:nrow(ruv_eval)) {
  if (ruv_eval$batch_drop[ek] <= 0) break  # Batch sil worsened â€” previous k was the elbow
}
best_k <- ruv_eval$k[ek - 1]

cat("Selected k =", best_k, "(elbow: last k before batch silhouette reversal)\n")
cat("Batch sil drops per k:", paste(round(ruv_eval$batch_drop[-1], 3), collapse = ", "), "\n")

write_csv(ruv_eval, file.path(OUTPUT_DIR, "ruv4_k_evaluation.csv"))

# Diagnostic plot: silhouette scores by k
ruv_long <- ruv_eval %>%
  dplyr::select(k, Group = sil_group, Batch = sil_batch) %>%
  pivot_longer(-k, names_to = "Metric", values_to = "Silhouette")

p_ruv <- ggplot(ruv_long, aes(k, Silhouette, colour = Metric)) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  geom_vline(xintercept = best_k, linetype = "dashed", colour = "black") +
  annotate("text", x = best_k + 0.15, y = max(ruv_eval$sil_group),
           label = paste0("k = ", best_k), hjust = 0, fontface = "bold") +
  scale_colour_manual(values = c(Group = "#C62828", Batch = "#1565C0")) +
  scale_x_continuous(breaks = 1:5) +
  labs(title = "RUV4 k Selection",
       subtitle = "Elbow: last k before batch silhouette reversal",
       x = "Number of RUV factors (k)", y = "Mean Silhouette Width") +
  theme_pub + theme(legend.position = "top", legend.title = element_blank())

ggsave(file.path(OUTPUT_DIR, "04b_ruv4_k_selection.pdf"), p_ruv, width = 6, height = 5)

#Select k=3 SGR DID 2
spe_ruv <- geomxBatchCorrection(spe_qc, factors = "Group_PDMSA", NCGs = metadata(spe_qc)$NCGs, k = 3)

spe_ruv$Group <- factor(spe_ruv$Group, levels = c("MSA", "PD"))
spe_ruv$Cell <- factor(spe_ruv$Cell, levels = c("Oligodendrocyte", "Neuron"))
spe_ruv$Group_PDMSA <- factor(spe_ruv$Group_PDMSA, levels = c("MSA_Neuron", "PD_Neuron", "MSA_Oligodendrocyte", "PD_Oligodendrocyte"))
spe_ruv$Patient <- factor(spe_ruv$Patient)
spe_ruv$SlideName <- factor(spe_ruv$SlideName)
spe_ruv$Batch <- factor(spe_ruv$Batch)

spe_ruv <- runPCA(spe_ruv, exprs_values = "logcounts", ncomponents = 10)
spe_ruv <- runUMAP(spe_ruv, dimred = "PCA")

saveRDS(spe_ruv, file.path(OUTPUT_DIR, "spe_ruv4.rds"))






###ROI-level heterogeneity##############################################################################
# Characterise within- vs between-patient variance and confirm batch correction. Informs the choice to use duplicateCorrelation in DE model.

cat("\n=== ROI-Level Heterogeneity ===\n")

pca_coords_full <- reducedDim(spe_ruv, "PCA")
pca_df <- data.frame(
  PC1 = pca_coords_full[, 1],
  PC2 = pca_coords_full[, 2],
  PC3 = pca_coords_full[, 3],
  Group = spe_ruv$Group_PDMSA,
  Patient = spe_ruv$Patient,
  Cell    = spe_ruv$Cell,
  Batch = spe_ruv$Batch
)

pca_var <- attr(reducedDim(spe_ruv, "PCA"), "percentVar")
if (is.null(pca_var)) {
  logcounts_mat <- assay(spe_ruv, "logcounts")
  pca_prcomp <- prcomp(t(logcounts_mat), center = TRUE, scale. = FALSE)
  pca_var <- (pca_prcomp$sdev^2 / sum(pca_prcomp$sdev^2)) * 100
}

hulls <- pca_df %>%
  group_by(Patient) %>%
  filter(n() >= 3) %>%
  slice(chull(PC1, PC2)) %>%
  ungroup()

# Within- vs between-patient distances (PCA space, first 5 PCs)
pca_5 <- pca_coords_full[, 1:min(5, ncol(pca_coords_full))]
dist_mat <- as.matrix(dist(pca_5))
patients <- as.character(spe_ruv$Patient)
n_rois <- length(patients)

within_dists <- c()
between_dists <- c()
for (i in 1:(n_rois - 1)) {
  for (j in (i + 1):n_rois) {
    if (patients[i] == patients[j]) {
      within_dists <- c(within_dists, dist_mat[i, j])
    } else {
      between_dists <- c(between_dists, dist_mat[i, j])
    }
  }
}

dist_df <- data.frame(
  Distance = c(within_dists, between_dists),
  Type = c(rep("Within-patient", length(within_dists)),
           rep("Between-patient", length(between_dists)))
)

wt_dist <- wilcox.test(within_dists, between_dists)
cat("Within-patient distance (median):", round(median(within_dists), 2), "\n")
cat("Between-patient distance (median):", round(median(between_dists), 2), "\n")
cat("Ratio (between/within):", round(median(between_dists) / median(within_dists), 2), "\n")
cat("Wilcoxon p:", signif(wt_dist$p.value, 3), "\n")

patient_stats <- pca_df %>%
  dplyr::group_by(Patient, Group, Cell) %>%
  dplyr::summarise(
    n_ROIs = dplyr::n(),
    PC1_sd = sd(PC1),
    PC2_sd = sd(PC2),
    spread = sqrt(PC1_sd^2 + PC2_sd^2),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Group, Patient) %>%
  dplyr::mutate(
    Patient = factor(Patient, levels = unique(Patient))   # <-- FIX
  )

group_centroids <- pca_df %>%
  group_by(Group, Cell) %>%
  summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop")

p_group <- ggplot(pca_df, aes(PC1, PC2, colour = Group, shape = Cell)) +
  geom_point(size = 2.5, alpha = 0.7) +
  stat_ellipse(level = 0.95, linewidth = 0.6, linetype = "dashed", show.legend = FALSE) +
  geom_point(data = group_centroids, aes(colour = Group),
             size = 4, shape = 3, stroke = 1.2, show.legend = FALSE) +
  scale_colour_manual(values = cols_groupcell, labels = group_labels) +
  labs(x = paste0("PC1 (", round(pca_var[1], 1), "%)"),
       y = paste0("PC2 (", round(pca_var[2], 1), "%)"),
       title = "After batch correction (RUV-4)",
       colour = NULL) +
  theme_pub + theme(legend.position = "top")

ggsave(file.path(OUTPUT_DIR, "02b_pca_post_ruv4.pdf"), p_group, width = 6, height = 5)

# 0) Ensure required columns are present in pca_df
pca_df$Group_PDMSA <- spe_ruv$Group_PDMSA   # <-- you map 'shape' to this
pca_df$Group       <- spe_ruv$Group         # if you still use Group anywhere

# 1) Ensure Group_PDMSA is a factor with the intended 4 levels
pca_df$Group_PDMSA <- factor(
  pca_df$Group_PDMSA,
  levels = c("MSA_Neuron","PD_Neuron","MSA_Oligodendrocyte","PD_Oligodendrocyte")
)

# 2) Define shapes that match those levels (you already have this—but keep here for clarity)
shape_vals <- c(
  "MSA_Neuron"            = 16,
  "PD_Neuron"             = 17,
  "MSA_Oligodendrocyte"   = 15,
  "PD_Oligodendrocyte"    = 18
)

# 3) Make a batch palette that matches the actual Batch levels 
pca_df$Batch <- droplevels(factor(pca_df$Batch))
cols_batch <- setNames(scales::hue_pal()(nlevels(pca_df$Batch)), levels(pca_df$Batch))


pca_df$Patient <- droplevels(factor(pca_df$Patient))
cols_patient <- setNames(scales::hue_pal()(nlevels(pca_df$Patient)),
                         levels(pca_df$Patient))

pdf(file.path(OUTPUT_DIR, "02_roi_heterogeneity.pdf"), width = 16, height = 12)


# colour by Patient, shape by Group_PDMSA
p1 <- ggplot(pca_df, aes(PC1, PC2)) +
  geom_polygon(data = hulls, aes(fill = Patient), alpha = 0.1, show.legend = FALSE) +
  geom_point(aes(colour = Patient, shape = Group_PDMSA), size = 3, alpha = 0.8) +
  scale_colour_manual(values = cols_patient) +
  scale_fill_manual(values = cols_patient) +
  scale_shape_manual(values = shape_vals) +
  labs(title = "PCA: Within-Patient ROI Clustering",
       subtitle = "Convex hulls show ROI spread per patient",
       x = paste0("PC1 (", round(pca_var[1], 1), "%)"),
       y = paste0("PC2 (", round(pca_var[2], 1), "%)")) +
  theme_pub

# colour by Batch, shape by Group_PDMSA
p2 <- ggplot(pca_df, aes(PC1, PC2)) +
  geom_point(aes(colour = Batch, shape = Group_PDMSA), size = 3, alpha = 0.8) +
  scale_colour_manual(values = cols_batch) +
  scale_shape_manual(values = shape_vals) +
  labs(title = "PCA: Batch Effects Post-RUV4",
       x = paste0("PC1 (", round(pca_var[1], 1), "%)"),
       y = paste0("PC2 (", round(pca_var[2], 1), "%)")) +
  theme_pub

p3 <- ggplot(dist_df, aes(Type, Distance, fill = Type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 0.5, alpha = 0.3) +
  scale_fill_manual(values = c("Within-patient" = "#4DAF4A", "Between-patient" = "#984EA3")) +
  annotate("text", x = 1.5, y = max(dist_df$Distance) * 0.95,
           label = paste0("p = ", signif(wt_dist$p.value, 3)), size = 4) +
  labs(title = "ROI Pairwise Distances (PC1-5)",
       subtitle = paste0("Median ratio: ",
                         round(median(between_dists) / median(within_dists), 2), "x"),
       x = NULL, y = "Euclidean Distance") +
  theme_pub + guides(fill = "none")

p4 <- ggplot(patient_stats, aes(Patient, n_ROIs, fill = Group)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = n_ROIs), vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = cols_group) +
  labs(title = "ROIs per Patient",
       subtitle = paste0("Total: ", ncol(spe_ruv), " ROIs from ",
                         nlevels(spe_ruv$Patient), " patients"),
       x = NULL, y = "Number of ROIs") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print((p1 + p2) / (p3 + p4) +
        plot_annotation(title = "ROI-Level Heterogeneity Assessment",
                        theme = theme(plot.title = element_text(size = 14, face = "bold"))))
invisible(dev.off())

write_csv(patient_stats, file.path(OUTPUT_DIR, "patient_roi_stats.csv"))
write_csv(dist_df, file.path(OUTPUT_DIR, "roi_pairwise_distances.csv"))



###Dimensionality reduction plots----
pdf(file.path(OUTPUT_DIR, "04_pca_umap.pdf"), width = 12, height = 5)
p1 <- plotDR(spe_ruv, dimred = "PCA", col = Group_PDMSA) + 
  scale_colour_manual(values = cols_groupcell) + labs(title = "PCA")
p2 <- plotDR(spe_ruv, dimred = "UMAP", col = Group_PDMSA) + 
  scale_colour_manual(values = cols_groupcell) + labs(title = "UMAP")
print(p1 + p2)
dev.off()


cat("\n=== Final Dataset ===\n")
cat("Genes:", nrow(spe_ruv), "\n")
cat("ROIs:", ncol(spe_ruv), "\n")
print(table(spe_ruv$Group_PDMSA, spe_ruv$SlideName))




###########################################################################################
### Differential Expression Analysis 
###########################################################################################

cat("\n=== Differential Expression Analysis ===\n")

library(limma)
library(edgeR)
library(dplyr)

cat("Using spe_qc:", ncol(spe_ruv), "ROIs\n")
cat("Sample distribution:\n")
print(table(spe_ruv$Group, spe_ruv$Cell))

###########################################################################################
### Helper function: Extract results with CI####################################
###########################################################################################

extract_de_with_ci <- function(fit, coef_name) {
  res <- topTable(fit, coef = coef_name, sort.by = "P", n = Inf)
  
  coef_idx <- which(colnames(fit$coefficients) == coef_name)
  se <- sqrt(fit$s2.post) * fit$stdev.unscaled[, coef_idx]
  
  res$SE <- se[rownames(res)]
  res$CI_lower <- res$logFC - 1.96 * res$SE
  res$CI_upper <- res$logFC + 1.96 * res$SE
  
  res %>%
    dplyr::select(logFC, CI_lower, CI_upper, SE, AveExpr, t, P.Value, adj.P.Val, B)
}

###########################################################################################
### Comparison 1: MSA vs PD (in All ROIs)##################################################################
###########################################################################################

cat("\n=== Comparison 1: MSA vs PD (All ROIs) ===\n")

# Prepare data
dge_all <- SE2DGEList(spe_ruv)

# Build design with RUV factors
ruv_cols <- grep("^ruv_W", colnames(colData(spe_ruv)), value = TRUE)
cat("RUV factors:", paste(ruv_cols, collapse = ", "), "\n")

design_all <- model.matrix(
  ~ 0 + Group + Cell + ruv_W1 + ruv_W2  + ruv_W3,
  data = colData(spe_ruv)
)
colnames(design_all) <- gsub("Group", "", colnames(design_all))
colnames(design_all) <- gsub("Cell", "Cell_", colnames(design_all))

# Contrast
contr_all <- makeContrasts(
  MSA_vs_PD = MSA - PD,
  levels = colnames(design_all)
)

# Filter low-expressed genes
keep <- filterByExpr(dge_all, design_all, min.count = 3)
cat("Genes retained:", sum(keep), "of", nrow(dge_all), "\n")
dge_all <- dge_all[keep, , keep.lib.sizes = FALSE]

# Estimate dispersion
dge_all <- estimateDisp(dge_all, design_all, robust = TRUE)

# Iterated voom-duplicateCorrelation
cat("Running iterated voom-duplicateCorrelation...\n")

# First iteration
v <- voom(dge_all, design_all, plot = FALSE)
corfit <- duplicateCorrelation(v, design_all, block = spe_ruv$Patient)

# Second iteration
v <- voom(dge_all, design_all, 
          block = spe_ruv$Patient,
          correlation = corfit$consensus.correlation,
          plot = FALSE)
corfit <- duplicateCorrelation(v, design_all, block = spe_ruv$Patient)

within_cor <- corfit$consensus.correlation

if (within_cor > 0.5)
  cat("NOTE: High within-patient correlation suggests substantial pseudoreplication.\n",
      "  Effective sample size is closer to n_patients than n_ROIs.\n",
      "  Interpret DEG counts with caution.\n")
cat("Within-patient correlation:", round(within_cor, 3), "\n")

if (within_cor > 0.5) {
  cat("NOTE: High within-patient correlation (", round(within_cor, 3), ")\n", sep="")
  cat("      Effective sample size closer to n_patients (", 
      length(unique(spe_ruv$Patient)), ") than n_ROIs (", ncol(spe_ruv), ")\n", sep="")
}

# Final voom with consensus correlation
v <- voom(dge_all, design_all, 
          block = spe_ruv$Patient, 
          correlation = within_cor,
          plot = FALSE)

# Fit model
fit <- lmFit(v, design_all, 
             block = spe_ruv$Patient,
             correlation = within_cor)
fit <- contrasts.fit(fit, contr_all)
fit <- eBayes(fit, robust = TRUE)

## Diagnostics: BCV, voom mean–variance, and eBayes trend
pdf(file.path(OUTPUT_DIR, "05a_diagnostics_All.pdf"), width = 12, height = 4)
op <- par(mfrow = c(1, 3), mar = c(5, 5, 3, 2))

# 1) BCV plot on the pre-voom DGEList
plotBCV(dge_all, main = "Biological Coefficient of Variation")

# 2) voom mean–variance trend
v_tmp <- voom(
  dge_all, design_all,
  block = spe_ruv$Patient,
  correlation = within_cor,
  plot = TRUE
)
title("Voom: Mean–variance trend", line = 1)

# 3) eBayes mean–variance trend after fitting
plotSA(fit, main = "eBayes: Mean–variance trend")

par(op)
dev.off()## Diagnostics: BCV, voom mean–variance, and eBayes trend

# Results
results_all <- decideTests(fit, adjust.method = "BH", p.value = 0.05)
cat("\nDEGs (FDR < 0.05):\n")
print(summary(results_all))

de_all <- extract_de_with_ci(fit, "MSA_vs_PD")
de_all_sig <- de_all %>% filter(adj.P.Val < 0.05)

cat("  Total DEGs:", nrow(de_all_sig), "\n")
cat("  Up in MSA:", sum(de_all_sig$logFC > 0), "\n")
cat("  Up in PD:", sum(de_all_sig$logFC < 0), "\n")

write.csv(de_all, file.path(OUTPUT_DIR, "DE_MSA_vs_PD_allROIs.csv"), row.names = TRUE)

###########################################################################################
### Comparison 2: MSA vs PD (Neurons only)######################################################
###########################################################################################

cat("\n=== Comparison 2: MSA vs PD (Neurons) ===\n")

spe_neuron <- spe_ruv[, spe_ruv$Cell == "Neuron"]
spe_neuron$Group <- droplevels(spe_neuron$Group)

cat("Neuron ROIs:", ncol(spe_neuron), "\n")
cat("  MSA:", sum(spe_neuron$Group == "MSA"), "\n")
cat("  PD:", sum(spe_neuron$Group == "PD"), "\n")

dge_neuron <- SE2DGEList(spe_neuron)

design_neuron <- model.matrix(
  ~ 0 + Group + ruv_W1 + ruv_W2 + ruv_W3,
  data = colData(spe_neuron)
)
colnames(design_neuron) <- gsub("Group", "", colnames(design_neuron))

contr_neuron <- makeContrasts(
  MSA_vs_PD = MSA - PD,
  levels = colnames(design_neuron)
)

keep <- filterByExpr(dge_neuron, design_neuron, min.count = 4)
cat("Genes retained:", sum(keep), "\n")
dge_neuron <- dge_neuron[keep, , keep.lib.sizes = FALSE]
dge_neuron <- estimateDisp(dge_neuron, design_neuron, robust = TRUE)

# Iterated voom-duplicateCorrelation
v <- voom(dge_neuron, design_neuron, plot = FALSE)
corfit <- duplicateCorrelation(v, design_neuron, block = spe_neuron$Patient)

v <- voom(dge_neuron, design_neuron,
          block = spe_neuron$Patient,
          correlation = corfit$consensus.correlation,
          plot = FALSE)
corfit <- duplicateCorrelation(v, design_neuron, block = spe_neuron$Patient)

within_cor <- corfit$consensus.correlation
cat("Within-patient correlation:", round(within_cor, 3), "\n")

v <- voom(dge_neuron, design_neuron,
          block = spe_neuron$Patient,
          correlation = within_cor,
          plot = FALSE)

fit <- lmFit(v, design_neuron,
             block = spe_neuron$Patient,
             correlation = within_cor)
fit <- contrasts.fit(fit, contr_neuron)
fit <- eBayes(fit, robust = TRUE)

## Diagnostics: BCV, voom mean–variance, and eBayes trend
pdf(file.path(OUTPUT_DIR, "05b_diagnostics_Neuron.pdf"), width = 12, height = 4)
op <- par(mfrow = c(1, 3), mar = c(5, 5, 3, 2))

# 1) BCV plot on the pre-voom DGEList
plotBCV(dge_neuron, main = "Biological Coefficient of Variation")

# 2) voom mean–variance trend
v_tmp <- voom(
  dge_neuron, design_neuron,
  block = spe_neuron$Patient,
  correlation = within_cor,
  plot = TRUE
)
title("Voom: Mean–variance trend", line = 1)

# 3) eBayes mean–variance trend after fitting
plotSA(fit, main = "eBayes: Mean–variance trend")

par(op)
dev.off()## Diagnostics: BCV, voom mean–variance, and eBayes trend

results_neuron <- decideTests(fit, adjust.method = "BH", p.value = 0.05)
cat("\nDEGs (FDR < 0.05):\n")
print(summary(results_neuron))

de_neuron <- extract_de_with_ci(fit, "MSA_vs_PD")
de_neuron_sig <- de_neuron %>% filter(adj.P.Val < 0.05)

cat("  Total DEGs:", nrow(de_neuron_sig), "\n")
cat("  Up in MSA:", sum(de_neuron_sig$logFC > 0), "\n")
cat("  Up in PD:", sum(de_neuron_sig$logFC < 0), "\n")

write.csv(de_neuron, file.path(OUTPUT_DIR, "DE_MSA_vs_PD_neuron.csv"), row.names = TRUE)

###########################################################################################
### Comparison 3: MSA vs PD (Oligodendrocytes only)##########################################
###########################################################################################

cat("\n=== Comparison 3: MSA vs PD (Oligodendrocytes) ===\n")

spe_oligo <- spe_ruv[, spe_ruv$Cell == "Oligodendrocyte"]
spe_oligo$Group <- droplevels(spe_oligo$Group)

cat("Oligo ROIs:", ncol(spe_oligo), "\n")
cat("  MSA:", sum(spe_oligo$Group == "MSA"), "\n")
cat("  PD:", sum(spe_oligo$Group == "PD"), "\n")

dge_oligo <- SE2DGEList(spe_oligo)

design_oligo <- model.matrix(
  ~ 0 + Group + ruv_W1 + ruv_W2 + ruv_W3,
  data = colData(spe_oligo)
)
colnames(design_oligo) <- gsub("Group", "", colnames(design_oligo))

contr_oligo <- makeContrasts(
  MSA_vs_PD = MSA - PD,
  levels = colnames(design_oligo)
)

keep <- filterByExpr(dge_oligo, design_oligo, min.count = 3)
cat("Genes retained:", sum(keep), "\n")
dge_oligo <- dge_oligo[keep, , keep.lib.sizes = FALSE]
dge_oligo <- estimateDisp(dge_oligo, design_oligo, robust = TRUE)

# Iterated voom-duplicateCorrelation
v <- voom(dge_oligo, design_oligo, plot = FALSE)
corfit <- duplicateCorrelation(v, design_oligo, block = spe_oligo$Patient)

v <- voom(dge_oligo, design_oligo,
          block = spe_oligo$Patient,
          correlation = corfit$consensus.correlation,
          plot = FALSE)
corfit <- duplicateCorrelation(v, design_oligo, block = spe_oligo$Patient)

within_cor <- corfit$consensus.correlation
cat("Within-patient correlation:", round(within_cor, 3), "\n")

v <- voom(dge_oligo, design_oligo,
          block = spe_oligo$Patient,
          correlation = within_cor,
          plot = FALSE)

fit <- lmFit(v, design_oligo,
             block = spe_oligo$Patient,
             correlation = within_cor)
fit <- contrasts.fit(fit, contr_oligo)
fit <- eBayes(fit, robust = TRUE)

## Diagnostics: BCV, voom mean–variance, and eBayes trend
pdf(file.path(OUTPUT_DIR, "05c_diagnostics_Olgo.pdf"), width = 12, height = 4)
op <- par(mfrow = c(1, 3), mar = c(5, 5, 3, 2))

# 1) BCV plot on the pre-voom DGEList
plotBCV(dge_oligo, main = "Biological Coefficient of Variation")

# 2) voom mean–variance trend
v_tmp <- voom(
  dge_oligo, design_oligo,
  block = spe_oligo$Patient,
  correlation = within_cor,
  plot = TRUE
)
title("Voom: Mean–variance trend", line = 1)

# 3) eBayes mean–variance trend after fitting
plotSA(fit, main = "eBayes: Mean–variance trend")

par(op)
dev.off()## Diagnostics: BCV, voom mean–variance, and eBayes trend

results_oligo <- decideTests(fit, adjust.method = "BH", p.value = 0.05)
cat("\nDEGs (FDR < 0.05):\n")
print(summary(results_oligo))

de_oligo <- extract_de_with_ci(fit, "MSA_vs_PD")
de_oligo_sig <- de_oligo %>% filter(adj.P.Val < 0.05)

cat("  Total DEGs:", nrow(de_oligo_sig), "\n")
cat("  Up in MSA:", sum(de_oligo_sig$logFC > 0), "\n")
cat("  Up in PD:", sum(de_oligo_sig$logFC < 0), "\n")

write.csv(de_oligo, file.path(OUTPUT_DIR, "DE_MSA_vs_PD_oligo.csv"), row.names = TRUE)

###########################################################################################
### Comparison 4: Neurons vs Oligodendrocytes (Cell Type Markers)########################
###########################################################################################

cat("\n=== Comparison 4: Neurons vs Oligodendrocytes ===\n")

spe_celltype <- spe_ruv[, spe_ruv$Cell %in% c("Neuron", "Oligodendrocyte")]
spe_celltype$Cell <- droplevels(spe_celltype$Cell)

cat("ROIs:", ncol(spe_celltype), "\n")
cat("  Neurons:", sum(spe_celltype$Cell == "Neuron"), "\n")
cat("  Oligos:", sum(spe_celltype$Cell == "Oligodendrocyte"), "\n")

dge_celltype <- SE2DGEList(spe_celltype)

# Design: Cell type as main factor, Group as covariate
design_celltype <- model.matrix(
  ~ 0 + Cell + Group + ruv_W1 + ruv_W2 + ruv_W3,
  data = colData(spe_celltype)
)
colnames(design_celltype) <- gsub("^Cell", "", colnames(design_celltype))
colnames(design_celltype) <- gsub("^Group", "Group_", colnames(design_celltype))

contr_celltype <- makeContrasts(
  Neuron_vs_Oligo = Neuron - Oligodendrocyte,
  levels = colnames(design_celltype)
)

keep <- filterByExpr(dge_celltype, design_celltype, min.count = 2)
cat("Genes retained:", sum(keep), "\n")
dge_celltype <- dge_celltype[keep, , keep.lib.sizes = FALSE]
dge_celltype <- estimateDisp(dge_celltype, design_celltype, robust = TRUE)

# Iterated voom-duplicateCorrelation
v <- voom(dge_celltype, design_celltype, plot = FALSE)
corfit <- duplicateCorrelation(v, design_celltype, block = spe_celltype$Patient)

v <- voom(dge_celltype, design_celltype,
          block = spe_celltype$Patient,
          correlation = corfit$consensus.correlation,
          plot = FALSE)
corfit <- duplicateCorrelation(v, design_celltype, block = spe_celltype$Patient)

within_cor <- corfit$consensus.correlation
cat("Within-patient correlation:", round(within_cor, 3), "\n")

v <- voom(dge_celltype, design_celltype,
          block = spe_celltype$Patient,
          correlation = within_cor,
          plot = FALSE)

fit <- lmFit(v, design_celltype,
             block = spe_celltype$Patient,
             correlation = within_cor)
fit <- contrasts.fit(fit, contr_celltype)
fit <- eBayes(fit, robust = TRUE)

## Diagnostics: BCV, voom mean–variance, and eBayes trend
pdf(file.path(OUTPUT_DIR, "05d_diagnostics_Celltype.pdf"), width = 12, height = 4)
op <- par(mfrow = c(1, 3), mar = c(5, 5, 3, 2))

# 1) BCV plot on the pre-voom DGEList
plotBCV(dge_celltype, main = "Biological Coefficient of Variation")

# 2) voom mean–variance trend
v_tmp <- voom(
  dge_celltype, design_celltype,
  block = spe_celltype$Patient,
  correlation = within_cor,
  plot = TRUE
)
title("Voom: Mean–variance trend", line = 1)

# 3) eBayes mean–variance trend after fitting
plotSA(fit, main = "eBayes: Mean–variance trend")

par(op)
dev.off()## Diagnostics: BCV, voom mean–variance, and eBayes trend

results_celltype <- decideTests(fit, adjust.method = "BH", p.value = 0.05)
cat("\nDEGs (FDR < 0.05):\n")
print(summary(results_celltype))

de_celltype <- extract_de_with_ci(fit, "Neuron_vs_Oligo")
de_celltype_sig <- de_celltype %>% filter(adj.P.Val < 0.05)

cat("  Total DEGs:", nrow(de_celltype_sig), "\n")
cat("  Neuron-enriched:", sum(de_celltype_sig$logFC > 0), "\n")
cat("  Oligo-enriched:", sum(de_celltype_sig$logFC < 0), "\n")

write.csv(de_celltype, file.path(OUTPUT_DIR, "DE_Neuron_vs_Oligodendrocyte.csv"), row.names = TRUE)

###########################################################################################
### Summary Table
###########################################################################################

cat("\n=== DE Analysis Summary ===\n")

summary_df <- data.frame(
  Comparison = c("MSA vs PD (All)", "MSA vs PD (Neurons)", 
                 "MSA vs PD (Oligos)", "Neurons vs Oligos"),
  N_ROIs = c(ncol(spe_ruv), ncol(spe_neuron), ncol(spe_oligo), ncol(spe_celltype)),
  N_Patients = c(length(unique(spe_ruv$Patient)),
                 length(unique(spe_neuron$Patient)),
                 length(unique(spe_oligo$Patient)),
                 length(unique(spe_celltype$Patient))),
  Within_Cor = c(
    round(duplicateCorrelation(voom(dge_all, design_all), design_all, block = spe_ruv$Patient)$consensus.correlation, 3),
    round(duplicateCorrelation(voom(dge_neuron, design_neuron), design_neuron, block = spe_neuron$Patient)$consensus.correlation, 3),
    round(duplicateCorrelation(voom(dge_oligo, design_oligo), design_oligo, block = spe_oligo$Patient)$consensus.correlation, 3),
    round(duplicateCorrelation(voom(dge_celltype, design_celltype), design_celltype, block = spe_celltype$Patient)$consensus.correlation, 3)
  ),
  Total_DEGs = c(nrow(de_all_sig), nrow(de_neuron_sig), nrow(de_oligo_sig), nrow(de_celltype_sig))
)

print(summary_df)
write.csv(summary_df, file.path(OUTPUT_DIR, "DE_summary.csv"), row.names = FALSE)

cat("\n=== DE Analysis Complete ===\n")




###############################################################################
###Volcano Plots
###############################################################################

library(ggplot2)
library(ggrepel)
library(dplyr)
library(tibble)
library(patchwork)

cat("\n=== Volcano Plots ===\n")

make_volcano <- function(
    de_results, 
    title,
    comparison_type = "disease",
    comparison_label = "MSA vs PD",
    n_label = 10,
    fdr_cutoff = 0.05,
    fc_cutoff = 0.5,
    label_top_by = "combined"
) {
  
  # Prepare data
  df <- de_results %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::mutate(
      neg_log10_p = -log10(P.Value),
      neg_log10_fdr = -log10(adj.P.Val),
      abs_logFC = abs(logFC),
      
      # Different categories based on comparison type
      category = if(comparison_type == "disease") {
        dplyr::case_when(
          adj.P.Val >= fdr_cutoff ~ "Not Significant",
          logFC > fc_cutoff ~ "Up in MSA",
          logFC < -fc_cutoff ~ "Up in PD",
          TRUE ~ "Significant\n(|logFC| < 0.5)"
        )
      } else {
        # Cell type comparison
        dplyr::case_when(
          adj.P.Val >= fdr_cutoff ~ "Not Significant",
          logFC > fc_cutoff ~ "Up in Neurons",
          logFC < -fc_cutoff ~ "Up in Oligodendrocytes",
          TRUE ~ "Significant\n(|logFC| < 0.5)"
        )
      },
      
      # Factor for plotting order
      category = if(comparison_type == "disease") {
        factor(
          category,
          levels = c("Not Significant", "Significant\n(|logFC| < 0.5)", 
                     "Up in PD", "Up in MSA")
        )
      } else {
        factor(
          category,
          levels = c("Not Significant", "Significant\n(|logFC| < 0.5)", 
                     "Up in Oligodendrocytes", "Up in Neurons")
        )
      },
      
      # Ranking score for labeling
      rank_score = case_when(
        label_top_by == "pvalue" ~ neg_log10_p,
        label_top_by == "logFC" ~ abs_logFC,
        label_top_by == "combined" ~ neg_log10_p * abs_logFC,
        TRUE ~ neg_log10_p * abs_logFC
      )
    ) %>%
    dplyr::arrange(desc(rank_score))
  
  # Count the significant genes
  n_total_sig <- sum(df$adj.P.Val < fdr_cutoff)
  
  if(comparison_type == "disease") {
    n_up_1 <- sum(df$adj.P.Val < fdr_cutoff & df$logFC > fc_cutoff)
    n_up_2 <- sum(df$adj.P.Val < fdr_cutoff & df$logFC < -fc_cutoff)
    label_1 <- "MSA"
    label_2 <- "PD"
  } else {
    n_up_1 <- sum(df$adj.P.Val < fdr_cutoff & df$logFC > fc_cutoff)
    n_up_2 <- sum(df$adj.P.Val < fdr_cutoff & df$logFC < -fc_cutoff)
    label_1 <- "Neurons"
    label_2 <- "Oligodendrocytes"
  }
  
  n_sig_small <- sum(df$adj.P.Val < fdr_cutoff & abs(df$logFC) < fc_cutoff)
  
  cat("\n", title, ":\n")
  cat("  Total DEGs (FDR <", fdr_cutoff, "):", n_total_sig, "\n")
  cat("    Up in", label_1, ":", n_up_1, "\n")
  cat("    Up in", label_2, ":", n_up_2, "\n")
  cat("    Significant but |logFC| <", fc_cutoff, ":", n_sig_small, "\n")
  
  # Select genes to label (top N by rank score, excluding NS)
  genes_to_label <- df %>%
    dplyr::filter(category != "Not Significant") %>%
    dplyr::slice_head(n = n_label)
  
  # Axis limits (symmetric for logFC)
  max_logfc <- max(abs(df$logFC), na.rm = TRUE)
  xlim_max <- ceiling(max_logfc * 1.1)
  ylim_max <- ceiling(max(df$neg_log10_p, na.rm = TRUE) * 1.1)
  
  if(comparison_type == "disease") {
    colors <- c(
      "Not Significant" = "grey80",
      "Significant\n(|logFC| < 0.5)" = "grey50",
      "Up in PD" = "#4C72B0",          # Blue
      "Up in MSA" = "#C44E52"          # red
    )
  } else {
    colors <- c(
      "Not Significant" = "grey80",
      "Significant\n(|logFC| < 0.5)" = "grey50",
      "Up in Oligodendrocytes" = "#009E73",  
      "Up in Neurons" = "#E69F00"            
    )
  }
  
  # Create plot
  p <- ggplot(df, aes(x = logFC, y = neg_log10_p)) +
    
    # Background grid (subtle)
    theme_classic(base_size = 11) +
    
    # Points (plot NS first so significant are on top)
    geom_point(
      data = df %>% filter(category == "Not Significant"),
      aes(color = category),
      size = 1.5,
      alpha = 0.3,
      shape = 16
    ) +
    geom_point(
      data = df %>% filter(category != "Not Significant"),
      aes(color = category),
      size = 2,
      alpha = 0.7,
      shape = 16
    ) +
    
    # Threshold lines
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed",
      color = "black",
      linewidth = 0.5,
      alpha = 0.7
    ) +
    geom_vline(
      xintercept = c(-fc_cutoff, fc_cutoff),
      linetype = "dashed",
      color = "black",
      linewidth = 0.5,
      alpha = 0.7
    ) +
    
    # Gene labels
    geom_text_repel(
      data = genes_to_label,
      aes(label = gene),
      size = 3,
      fontface = "italic",
      max.overlaps = Inf,
      min.segment.length = 0,
      segment.size = 0.3,
      segment.color = "grey30",
      segment.alpha = 0.6,
      box.padding = 0.4,
      point.padding = 0.3,
      force = 2,
      force_pull = 0.5,
      seed = 123
    ) +
    
    # Colors
    scale_color_manual(
      values = colors,
      name = NULL,
      labels = if(comparison_type == "disease") {
        c(
          "Not Significant" = "Not significant",
          "Significant\n(|logFC| < 0.5)" = paste0("FDR < 0.05, |logFC| < 0.5 (", n_sig_small, ")"),
          "Up in PD" = paste0("Up in PD (", n_up_2, ")"),
          "Up in MSA" = paste0("Up in MSA (", n_up_1, ")")
        )
      } else {
        c(
          "Not Significant" = "Not significant",
          "Significant\n(|logFC| < 0.5)" = paste0("FDR < 0.05, |logFC| < 0.5 (", n_sig_small, ")"),
          "Up in Oligodendrocytes" = paste0("Up in Oligodendrocytes (", n_up_2, ")"),
          "Up in Neurons" = paste0("Up in Neurons (", n_up_1, ")")
        )
      }
    ) +
    
    # Scales
    scale_x_continuous(
      limits = c(-xlim_max, xlim_max),
      expand = c(0.02, 0)
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.02, 0.05))
    ) +
    
    # Labels
    labs(
      x = if(comparison_type == "disease") {
        expression(bold("log"[2]*" fold change (MSA / PD)"))
      } else {
        expression(bold("log"[2]*" fold change (Neuron / Oligodendrocyte)"))
      },
      y = expression(bold("-log"[10]*" "*italic(P)*"-value")),
      title = title,
      subtitle = paste0(
        n_total_sig, " DEGs (FDR < ", fdr_cutoff, ") | ",
        n_up_1 + n_up_2, " with |logFC| > ", fc_cutoff
      )
    ) +
    
    theme(
      # Plot area
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white", color = NA),
      
      # Grid
      panel.grid.major = element_line(color = "grey95", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      
      # Axes
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.text = element_text(color = "black", size = 10),
      axis.title = element_text(face = "bold", size = 11),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10)),
      
      # Title
      plot.title = element_text(
        face = "bold",
        size = 13,
        hjust = 0.5,
        margin = margin(b = 5)
      ),
      plot.subtitle = element_text(
        size = 10,
        hjust = 0.5,
        color = "grey30",
        margin = margin(b = 10)
      ),
      
      # Legend
      legend.position = "bottom",
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      legend.text = element_text(size = 9),
      legend.margin = margin(t = 5),
      legend.box.spacing = unit(0.2, "cm")
    ) +
    
    # Legend formatting
    guides(
      color = guide_legend(
        override.aes = list(size = 3.5, alpha = 1),
        nrow = 2,
        byrow = TRUE
      )
    )
  
  return(p)
}

###############################################################################
### Create plots
###############################################################################

cat("\n--- Creating volcano plots ---\n")

# 1. MSA vs PD (All ROIs)
p_vol_all <- make_volcano(
  de_results = de_all,
  title = "Differential Expression: MSA vs PD",
  comparison_type = "disease",
  comparison_label = "MSA vs PD",
  n_label = 15,
  fdr_cutoff = 0.05,
  fc_cutoff = 0.5,
  label_top_by = "combined"
)

# 2. MSA vs PD (Neurons)
p_vol_neuron <- make_volcano(
  de_results = de_neuron,
  title = "Differential Expression in Neurons: MSA vs PD",
  comparison_type = "disease",
  comparison_label = "MSA vs PD",
  n_label = 15,
  fdr_cutoff = 0.05,
  fc_cutoff = 0.5,
  label_top_by = "combined"
)

# 3. MSA vs PD (Oligodendrocytes)
p_vol_oligo <- make_volcano(
  de_results = de_oligo,
  title = "Differential Expression in Oligodendrocytes: MSA vs PD",
  comparison_type = "disease",
  comparison_label = "MSA vs PD",
  n_label = 15,
  fdr_cutoff = 0.05,
  fc_cutoff = 0.5,
  label_top_by = "combined"
)

# 4. Neuron vs Oligo
if (exists("de_celltype")) {
  p_vol_celltype <- make_volcano(
    de_results = de_celltype,
    title = "Differential Expression: Neuron vs Oligodendrocyte",
    comparison_type = "celltype",  
    comparison_label = "Neuron vs Oligo",
    n_label = 15,
    fdr_cutoff = 0.05,
    fc_cutoff = 0.5,
    label_top_by = "combined"
  )
}

###############################################################################
### Save
###############################################################################

ggsave(
  file.path(OUTPUT_DIR, "06a_volcano_MSA_vs_PD_all.pdf"),
  p_vol_all,
  width = 8, height = 8
)

ggsave(
  file.path(OUTPUT_DIR, "06b_volcano_MSA_vs_PD_neurons.pdf"),
  p_vol_neuron,
  width = 8, height = 8
)

ggsave(
  file.path(OUTPUT_DIR, "06c_volcano_MSA_vs_PD_oligos.pdf"),
  p_vol_oligo,
  width = 8, height = 8
)

if (exists("de_celltype")) {
  ggsave(
    file.path(OUTPUT_DIR, "06d_volcano_neuron_vs_oligo.pdf"),
    p_vol_celltype,
    width = 8, height = 8
  )
}

cat("\n✓ Saved individual volcano plots\n")

cat("\n=== Volcano plots complete ===\n")
cat("Files created:\n")
cat("  06a-d: Individual volcano plots\n")
}






####################################################################################
###GSEA KEGG and SYNGO Analysis##########################################
##############################################################################
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)

# ---- Load Syngo Excel files ----
syngo_anno  <- read_excel(file.path(datadir, "syngo", "annotations.xlsx"))
syngo_genes <- read_excel(file.path(datadir, "syngo", "genes.xlsx"))
syngo_onto  <- read_excel(file.path(datadir, "syngo", "ontologies.xlsx"))

# ---- Standardize column names ----
syngo_anno  <- syngo_anno  %>% rename_with(tolower)
syngo_genes <- syngo_genes %>% rename_with(tolower)
syngo_onto  <- syngo_onto  %>% rename_with(tolower)

# ---- Clean & normalise ----
syngo_anno <- syngo_anno %>%
  mutate(
    hgnc_id     = trimws(as.character(hgnc_id)),   # already in HGNC:#### format
    hgnc_symbol = toupper(trimws(hgnc_symbol)),
    go_domain   = tolower(trimws(go_domain))
  )

syngo_genes <- syngo_genes %>%
  mutate(
    hgnc_id     = trimws(as.character(hgnc_id)),   # matches annotations
    hgnc_symbol = toupper(trimws(hgnc_symbol)),
    entrez_id   = trimws(as.character(entrez_id))
  )

# ---- JOIN BY HGNC ID ----
joined_id <- syngo_anno %>%
  left_join(
    syngo_genes %>% dplyr::select(hgnc_id, entrez_id),
    by = "hgnc_id"
  )

joined_id_ok <- joined_id %>%
  filter(!is.na(entrez_id)) %>%
  dplyr::select(go_name, go_domain, hgnc_symbol, hgnc_id, entrez_id)

# ---- Handle rows that failed to join by ID → join by SYMBOL ----
anno_missing <- joined_id %>%
  filter(is.na(entrez_id)) %>%
  dplyr::select(go_name, go_domain, hgnc_symbol) %>%
  mutate(
    hgnc_symbol = str_replace_all(hgnc_symbol, "[,\\/]", ";")
  ) %>%
  separate_rows(hgnc_symbol, sep = "\\s*;\\s*") %>%
  filter(hgnc_symbol != "")

recovered_by_symbol <- anno_missing %>%
  left_join(
    syngo_genes %>% dplyr::select(hgnc_symbol, entrez_id),
    by = "hgnc_symbol"
  ) %>%
  filter(!is.na(entrez_id)) %>%
  distinct(go_name, go_domain, hgnc_symbol, entrez_id)

# ---- Combine ID matching + symbol matching ----
syngo_join <- bind_rows(joined_id_ok, recovered_by_symbol) %>%
  distinct(go_name, go_domain, entrez_id)

# ---- Domain → BP/CC mapping ----
syngo_join <- syngo_join %>%
  mutate(
    go_domain = tolower(go_domain),   # you have "bp" and "cc"
    branch = case_when(
      go_domain == "bp" ~ "BP",
      go_domain == "cc" ~ "CC",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(branch))

# ---- Build TERM2GENE ----
syngo_bp_t2g <- syngo_join %>%
  filter(branch == "BP") %>%
  dplyr::select(term_name = go_name, entrez_id) %>%
  distinct()

syngo_cc_t2g <- syngo_join %>%
  filter(branch == "CC") %>%
  dplyr::select(term_name = go_name, entrez_id) %>%
  distinct()

# ---- Build TERM2NAME using ontology file ----
syngo_t2n <- syngo_onto %>%
  dplyr::select(term_name = name, definition = name) %>%
  distinct()


# ---- Diagnostics ----
cat("SynGO loaded:\n")
cat("  matched by ID:", nrow(joined_id_ok), "\n")
cat("  recovered by SYMBOL:", nrow(recovered_by_symbol), "\n")
cat("  final links:", nrow(syngo_join), "\n")
cat("  BP terms:", n_distinct(syngo_bp_t2g$term_name), "\n")
cat("  CC terms:", n_distinct(syngo_cc_t2g$term_name), "\n")


## ---- GSEA Analysis --------------------------------------------------------
cat("\n=== GSEA Analysis ===\n")

# Function to run GSEA for a given DE table
run_gsea <- function(de_results, comparison_name) {
  cat("\nRunning GSEA for:", comparison_name, "\n")
  
  if (!"t" %in% colnames(de_results)) {
    stop("Column 't' not found in de_results. Make sure this is limma output with a t-statistic.")
  }
  
  gene_list <- de_results$t
  symbols   <- rownames(de_results)
  
  # Map SYMBOL -> ENTREZ (deduplicate)
  ids <- mapIds(
    org.Hs.eg.db,
    keys     = symbols,
    column   = "ENTREZID",
    keytype  = "SYMBOL",
    multiVals = "first"
  )
  
  valid_idx <- !is.na(ids) & !duplicated(ids)
  gene_list        <- gene_list[valid_idx]
  names(gene_list) <- ids[valid_idx]
  
  gene_list <- sort(gene_list, decreasing = TRUE)
  cat("  Gene list size:", length(gene_list), "\n")
  
  # ---- GO BP ----
  gse_go <- tryCatch({
    gseGO(geneList = gene_list,
          ont = "BP",
          keyType = "ENTREZID",
          minGSSize = 15,
          maxGSSize = 500,
          pvalueCutoff = 0.5,
          OrgDb = org.Hs.eg.db,
          pAdjustMethod = "BH",
          nPermSimple = 100000,
          eps = 0,
          verbose = FALSE)
  }, error = function(e) { cat("  GO error:", e$message, "\n"); NULL })
  
  # ---- KEGG ----
  gse_kegg <- tryCatch({
    gseKEGG(geneList = gene_list,
            keyType = "ncbi-geneid",
            organism = "hsa",
            minGSSize = 5,
            maxGSSize = 500,
            pvalueCutoff = 0.5,
            pAdjustMethod = "BH",
            nPermSimple = 100000,
            eps = 0,
            verbose = FALSE)
  }, error = function(e) { cat("  KEGG error:", e$message, "\n"); NULL })
  
  # ---- SynGO BP ----
  gse_syngo_bp <- tryCatch({
    GSEA(
      geneList      = gene_list,
      TERM2GENE     = syngo_bp_t2g,
      TERM2NAME     = if (nrow(syngo_t2n)) syngo_t2n else NULL,
      pvalueCutoff  = 0.5,
      pAdjustMethod = "BH",
      minGSSize     = 5,
      maxGSSize     = 500,
      verbose       = FALSE
    )
  }, error = function(e) { cat("  SynGO BP error:", e$message, "\n"); NULL })
  
  # ---- SynGO CC ----
  gse_syngo_cc <- tryCatch({
    GSEA(
      geneList      = gene_list,
      TERM2GENE     = syngo_cc_t2g,
      TERM2NAME     = if (nrow(syngo_t2n)) syngo_t2n else NULL,
      pvalueCutoff  = 0.5,
      pAdjustMethod = "BH",
      minGSSize     = 5,
      maxGSSize     = 500,
      verbose       = FALSE
    )
  }, error = function(e) { cat("  SynGO CC error:", e$message, "\n"); NULL })
  
  # ---- Diagnostics ----
  if (!is.null(gse_go)) {
    n_na <- sum(is.na(gse_go@result$pvalue))
    cat("  GO pathways tested:", nrow(gse_go@result), "(", n_na, "with NA p-values)\n")
    cat("  GO qval < 0.05:", sum(gse_go@result$qvalue < 0.05, na.rm = TRUE), "\n")
    cat("  GO qval < 0.10:", sum(gse_go@result$qvalue < 0.10, na.rm = TRUE), "\n")
    cat("  GO qval < 0.25:", sum(gse_go@result$qvalue < 0.25, na.rm = TRUE), "\n")
  }
  if (!is.null(gse_kegg)) {
    n_na <- sum(is.na(gse_kegg@result$pvalue))
    cat("  KEGG pathways tested:", nrow(gse_kegg@result), "(", n_na, "with NA p-values)\n")
    cat("  KEGG qval < 0.05:", sum(gse_kegg@result$qvalue < 0.05, na.rm = TRUE), "\n")
    cat("  KEGG qval < 0.10:", sum(gse_kegg@result$qvalue < 0.10, na.rm = TRUE), "\n")
    cat("  KEGG qval < 0.25:", sum(gse_kegg@result$qvalue < 0.25, na.rm = TRUE), "\n")
  }
  if (!is.null(gse_syngo_bp)) {
    n_na <- sum(is.na(gse_syngo_bp@result$pvalue))
    cat("  SynGO BP tested:", nrow(gse_syngo_bp@result), "(", n_na, "with NA p-values)\n")
  }
  if (!is.null(gse_syngo_cc)) {
    n_na <- sum(is.na(gse_syngo_cc@result$pvalue))
    cat("  SynGO CC tested:", nrow(gse_syngo_cc@result), "(", n_na, "with NA p-values)\n")
  }
  
  return(list(
    GO        = gse_go,
    KEGG      = gse_kegg,
    SYNGO_BP  = gse_syngo_bp,
    SYNGO_CC  = gse_syngo_cc,
    gene_list = gene_list
  ))
}
gsea_msa_vs_pd_oligo <- run_gsea(de_oligo,  "MSA vs PD (Oligodendrocytes)")
gsea_msa_vs_pd_neuron <- run_gsea(de_neuron, "MSA vs PD (Neurons)")
gsea_msa_vs_pd_all <- run_gsea(de_all,       "MSA vs PD (All ROIs)")


saveRDS(gsea_msa_vs_pd_oligo,  file.path(OUTPUT_DIR, "gsea_oligo_full.rds"))
saveRDS(gsea_msa_vs_pd_neuron, file.path(OUTPUT_DIR, "gsea_neuron_full.rds"))
saveRDS(gsea_msa_vs_pd_all,    file.path(OUTPUT_DIR, "gsea_all_full.rds"))
cat("\n=== GSEA Complete ===\n")




# # # # # # # # # # # #
#REACTTOME
# # # # # # # #
library(ReactomePA)
library(reactome.db)
library(fgsea)
run_reactome_gsea <- function(de_results, comparison_name, minSize = 15, maxSize = 500) {
  
  cat("\nRunning Reactome GSEA for:", comparison_name, "\n")
  
  # Prepare ranked gene list (ENTREZ IDs)
  gene_list <- de_results$t
  symbols <- rownames(de_results)
  
  ids <- mapIds(
    org.Hs.eg.db,
    keys = symbols,
    column = "ENTREZID",
    keytype = "SYMBOL",
    multiVals = "first"
  )
  
  valid_idx <- !is.na(ids) & !duplicated(ids)
  gene_list <- gene_list[valid_idx]
  names(gene_list) <- ids[valid_idx]
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  pathways_raw <- as.list(reactomePATHID2EXTID)
  
  # Filter pathways by gene set size
  pathways_filtered <- pathways_raw[sapply(pathways_raw, function(x) {
    len <- length(intersect(x, names(gene_list)))
    len >= minSize && len <= maxSize
  })]
  
  cat("  Pathways tested:", length(pathways_filtered), "\n")
  
  # Run fgseaMultilevel
  fgsea_res <- fgseaMultilevel(
    pathways = pathways_filtered,
    stats = gene_list,
    eps = 0
  )
  
  # Convert leadingEdge to comma-separated strings & add Description
  path_names <- as.list(reactomePATHID2NAME)
  
  fgsea_res <- fgsea_res %>%
    mutate(
      leadingEdge = sapply(leadingEdge, function(x) paste(x, collapse = ";")),
      Description = unlist(path_names[pathway])
    ) %>%
    dplyr::select(Description, everything()) %>%
    arrange(padj)
  
  # Map Reactome IDs to names
  path_names <- as.list(reactomePATHID2NAME)
  
  fgsea_res <- fgsea_res %>%
    mutate(leadingEdge = sapply(leadingEdge, function(x) paste(x, collapse = ";")),
           Description = unlist(path_names[pathway])) %>%
    dplyr::select(Description, everything()) %>%
    arrange(padj)
  
  cat("  Reactome pathways found:", nrow(fgsea_res), "\n")
  
  write.csv(
    fgsea_res,
    file.path(OUTPUT_DIR, paste0("GSEA_Reactome_", gsub(" ", "_", comparison_name), ".csv")),
    row.names = FALSE
  )
  
  return(list(results = fgsea_res, gene_list = gene_list))
}

reactome_msa_vs_pd_oligo <- run_reactome_gsea(de_oligo, "MSA_vs_PD_Oligo")
reactome_msa_vs_pd_neuron <- run_reactome_gsea(de_neuron, "MSA_vs_PD_Neuron")
reactome_msa_vs_pd_all <- run_reactome_gsea(de_all, "MSA_vs_PD_All")

























###############################################################################
###Gene Expression Plots######################################################
###############################################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(patchwork)
library(ggbeeswarm)
library(ggsignif)

cat("\n=== Creating Gene Expression Plots ===\n")

plot_gene_expression <- function(
    genes_of_interest,
    plot_title,
    logmat,
    meta_df,
    cols_group,
    de_results,
    show_stats = TRUE,
    ncol = NULL) {
  
  genes_present <- genes_of_interest[genes_of_interest %in% rownames(logmat)]
  
  if (length(genes_present) == 0) {
    cat("No genes found in expression data\n")
    return(NULL)
  }
  
  cat("Plotting", length(genes_present), "genes\n")
  
  # Convert to long format
  expr_long <- as.data.frame(t(logmat[genes_present, , drop = FALSE])) %>%
    rownames_to_column("ROI") %>%
    pivot_longer(cols = all_of(genes_present),
                 names_to = "gene",
                 values_to = "logcounts") %>%
    left_join(meta_df %>% dplyr::select(ROI, Group, Cell, Patient), by = "ROI") %>%
    mutate(
      gene = factor(gene, levels = genes_present),
      Group = factor(Group, levels = c("PD", "MSA")),  # PD first (control)
      Cell = factor(Cell, levels = c("Neuron", "Oligodendrocyte"))
    )
  
  # Calculate summary statistics
  summary_stats <- expr_long %>%
    group_by(gene, Cell, Group) %>%
    summarise(
      mean = mean(logcounts, na.rm = TRUE),
      sem = sd(logcounts, na.rm = TRUE) / sqrt(n()),
      median = median(logcounts, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
  
  n_genes <- length(genes_present)
  if (is.null(ncol)) {
    ncol <- ifelse(n_genes <= 3, n_genes, 3)
  }
  
  # Create plot
  p <- ggplot(expr_long, aes(x = Group, y = logcounts, fill = Group)) +
    
    # Violin plot
    geom_violin(alpha = 0.2, color = NA, scale = "width", trim = TRUE) +
    # Boxplot
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.7,
                 position = position_dodge(0.9)) +
    
    # Individual points
    geom_beeswarm(alpha = 0.6, size = 1.2, cex = 2.5,
                  dodge.width = 0.9, color = "black") +
    
    # Mean ± SEM (error bars)
    geom_errorbar(data = summary_stats,
                  aes(x = Group, y = mean, ymin = mean - sem, ymax = mean + sem),
                  width = 0.15, linewidth = 0.8, color = "black",
                  position = position_dodge(0.9)) +
    
    # Mean point
    geom_point(data = summary_stats,
               aes(x = Group, y = mean),
               size = 3, shape = 18, color = "black",
               position = position_dodge(0.9)) +
    
    # Faceting
    facet_grid(Cell ~ gene, scales = "fixed", switch = "y") +
    
    # Colors
    scale_fill_manual(values = cols_group) +
    
    # Labels
    labs(
      x = NULL,
      y = expression(bold("Expression (log"[2]*" CPM)")),
      title = plot_title
    ) +
  
    # Theme
    theme_classic(base_size = 11) +
    theme(
      # Strip (facet labels)
      strip.background = element_rect(fill = "grey95", color = "black", linewidth = 0.5),
      strip.text.x = element_text(face = "bold.italic", margin = margin(3, 0, 3, 0)),
      strip.text.y = element_text(face = "bold", angle = 0, hjust = 0.5),
      strip.placement = "outside",
      
      # Axes
      axis.text.x = element_text(face = "bold", color = "black"),
      axis.text.y = element_text(color = "black"),
      axis.title.y = element_text(face = "bold", size = 11, margin = margin(0, 10, 0, 0)),
      axis.line = element_line(linewidth = 0.5, color = "black"),
      axis.ticks = element_line(linewidth = 0.5, color = "black"),
      
      # Title
      plot.title = element_text(face = "bold", size = 13, hjust = 0),
      
      # Legend
      legend.position = "none",
      
      # Panel
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.spacing = unit(0.3, "lines"),
      
      # Background
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  # Add DE statistics from limma-voom
  if (show_stats && n_genes <= 6) {
    
    stat_results <- de_results[genes_present, , drop = FALSE] %>%
      as.data.frame() %>%
      rownames_to_column("gene") %>%
      mutate(
        p_label = case_when(
          adj.P.Val < 0.0001 ~ "****",
          adj.P.Val < 0.001  ~ "***",
          adj.P.Val < 0.01   ~ "**",
          adj.P.Val < 0.05   ~ "*",
          TRUE               ~ "ns"
        )
      )
    
    # repeat labels for both cell rows
    stat_results <- tidyr::crossing(
      stat_results,
      Cell = levels(expr_long$Cell)
    )
    
    p <- p +
      geom_text(
        data = stat_results,
        aes(x = 1.5, y = Inf, label = p_label),
        vjust = 1.2,
        size = 4.2,
        inherit.aes = FALSE,
        fontface = "bold"
      )
  }
  
  return(p)
}

###############################################################################
### Prepare data
###############################################################################

logmat <- assay(spe_ruv, "logcounts")

# Metadata
meta_df <- as.data.frame(colData(spe_qc)) %>% 
  rownames_to_column("ROI")

cols_group <- c(PD = "#4C72B0", MSA = "#C44E52")

###############################################################################
### Gene sets
###############################################################################

# Top 6 DEGs from all ROIs analysis
top_degs_all <- de_all %>%
  as.data.frame() %>%
  filter(adj.P.Val < 0.05, logFC > 0) %>%
  arrange(adj.P.Val, desc(logFC)) %>%
  head(6) %>%
  rownames()

# Top 6 DEGs from neuron analysis
top_degs_neuron <- de_neuron %>%
  as.data.frame() %>%
  filter(adj.P.Val < 0.05, logFC > 0) %>%
  arrange(adj.P.Val, desc(logFC)) %>%
  head(6) %>%
  rownames()

# Top 6 DEGs from oligo analysis
top_degs_oligo <- de_oligo %>%
  as.data.frame() %>%
  filter(adj.P.Val < 0.05, logFC > 0) %>%
  arrange(adj.P.Val, desc(logFC)) %>%
  head(6) %>%
  rownames()


# PD/neurodegeneration candidate genes
pd_candidates <- c("SNCA", "PARK7", "LRRK2", "PINK1", "PRKN", "PARK2", "PARK8",
                  "GBA", "MAPT", "COQ2")
# Filter to genes present in data
pd_candidates_present <- intersect(pd_candidates, rownames(logmat))

cat("\nGene sets:\n")
cat("  Top DEGs (All):", length(top_degs_all), "\n")
cat("  Top DEGs (Neuron):", length(top_degs_neuron), "\n")
cat("  Top DEGs (Oligo):", length(top_degs_oligo), "\n")
cat("  PD candidates present:", length(pd_candidates_present), "/", 
    length(pd_candidates), "\n")

###############################################################################
### Create plots
###############################################################################

cat("\n--- Creating plots ---\n")

# Top 6 DEGs (All ROIs)
p_top_all <- plot_gene_expression(
  genes_of_interest = top_degs_all,
  plot_title = "Top Differentially Expressed Genes (MSA vs PD)",
  logmat = logmat,
  meta_df = meta_df,
  cols_group = cols_group,
  de_results = de_all,
  show_stats = TRUE
)

ggsave(
  file.path(OUTPUT_DIR, "08a_expression_top_DEGs.pdf"),
  p_top_all,
  width = 10, height = 6,
  device = pdf
)

# Top 6 DEGs (Neurons only)
meta_neuron <- meta_df %>% filter(Cell == "Neuron")

rois_neuron <- meta_neuron$ROI
logmat_neuron <- logmat[, rois_neuron]

p_top_neuron <- plot_gene_expression(
  genes_of_interest = top_degs_neuron,
  plot_title = "Top Differentially Expressed Genes in Neurons",
  logmat = logmat_neuron,
  meta_df = meta_neuron,
  cols_group = cols_group,
  de_results = de_neuron,
  show_stats = TRUE
)

# Modify to remove Cell faceting
p_top_neuron <- p_top_neuron + 
  facet_wrap(~gene, scales = "fixed", ncol = 3)
ggsave(
  file.path(OUTPUT_DIR, "08b_expression_neuron_DEGs.pdf"),
  p_top_neuron,
  width = 9, height = 6,
  device = pdf
)

# Top 6 DEGs (Oligos only)
meta_oligo <- meta_df %>% filter(Cell == "Oligodendrocyte")

rois_oligo <- meta_oligo$ROI
logmat_oligo <- logmat[, rois_oligo]

p_top_oligo <- plot_gene_expression(
  genes_of_interest = top_degs_oligo,
  plot_title = "Top Differentially Expressed Genes in Oligodendrocytes",
  logmat = logmat_oligo,
  meta_df = meta_oligo,
  cols_group = cols_group,
  de_results = de_oligo,
  show_stats = TRUE
)

p_top_oligo <- p_top_oligo + 
  # facet_wrap(~gene, scales = "free_y", ncol = 3)
  facet_wrap(~gene, scales = "fixed", ncol = 3)
ggsave(
  file.path(OUTPUT_DIR, "08c_expression_oligo_DEGs.pdf"),
  p_top_oligo,
  width = 9, height = 6,
  device = pdf
)

# PD candidate genes
if (length(pd_candidates_present) > 0) {
  
  p_pd <- plot_gene_expression(
    genes_of_interest = pd_candidates_present,
    plot_title = "Parkinson's Disease-Associated Genes",
    logmat = logmat,
    meta_df = meta_df,
    cols_group = cols_group,
    de_results = de_all,
    show_stats = TRUE
  )

  
  ggsave(
    file.path(OUTPUT_DIR, "08d_expression_PD_candidates.pdf"),
    p_pd,
    width = 9, height = 6,
    device = pdf
  )
}

cat("\n===plots complete ===\n")
cat("Files created:\n")
cat("  08a_expression_top_DEGs.pdf\n")
cat("  08b_expression_neuron_DEGs.pdf\n")
cat("  08c_expression_oligo_DEGs.pdf\n")
cat("  08d_expression_PD_candidates.pdf\n")


###############################################################################
### Cell Type Marker Genes (Neuron vs Oligodendrocyte)
###############################################################################

cat("\n--- Creating cell type marker plot ---\n")

# Top markers from Neuron vs Oligodendrocyte DE analysis
top_neuron_markers <- de_celltype %>%
  as.data.frame() %>%
  filter(adj.P.Val < 0.05, logFC > 0) %>%   # Up in Neuron
  arrange(adj.P.Val, desc(logFC)) %>%
  head(6) %>%
  rownames()

top_oligo_markers <- de_celltype %>%
  as.data.frame() %>%
  filter(adj.P.Val < 0.05, logFC < 0) %>%   
  arrange(adj.P.Val, (logFC)) %>%                         
  head(6) %>%
  rownames()

celltype_markers <- c(top_neuron_markers, top_oligo_markers)

cat("  Neuron markers:", paste(top_neuron_markers, collapse = ", "), "\n")
cat("  Oligo markers:  ", paste(top_oligo_markers, collapse = ", "), "\n")

stopifnot(!any(duplicated(celltype_markers)))

###############################################################################
### Prepare expression data
###############################################################################

expr_long <- as.data.frame(t(logmat[celltype_markers, , drop = FALSE])) %>%
  rownames_to_column("ROI") %>%
  pivot_longer(
    cols = all_of(celltype_markers),
    names_to = "gene",
    values_to = "logcounts"
  ) %>%
  left_join(
    meta_df %>% dplyr::select(ROI, Cell),
    by = "ROI"
  ) %>%
  mutate(
    gene = factor(gene, levels = celltype_markers),
    Cell = factor(Cell, levels = c("Neuron", "Oligodendrocyte"))
  )
###############################################################################
### Statistical testing — limma/voom FDRs from de_celltype
###############################################################################

global_max <- max(expr_long$logcounts, na.rm = TRUE)
global_min <- min(expr_long$logcounts, na.rm = TRUE)

star_y_pos <- global_max + (global_max - global_min) * 0.08

stat_results <- de_celltype[celltype_markers, , drop = FALSE] %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  mutate(
    p_label = case_when(
      adj.P.Val < 0.0001 ~ "****",
      adj.P.Val < 0.001  ~ "***",
      adj.P.Val < 0.01   ~ "**",
      adj.P.Val < 0.05   ~ "*",
      TRUE               ~ "ns"
    ),
    y_pos = star_y_pos
  )
###############################################################################
### Plot
###############################################################################

p_celltype_markers <- ggplot(
  expr_long,
  aes(x = Cell, y = logcounts, fill = Cell)
) +
  geom_violin(alpha = 0.2, color = NA) +
  geom_boxplot(width = 0.20, outlier.shape = NA, alpha = 0.7) +
  geom_beeswarm(size = 1.2, alpha = 0.6, color = "black") +
  
  # Significance stars — same y position across all panels
  geom_text(
    data = stat_results,
    aes(x = 1.5, y = y_pos, label = p_label),
    inherit.aes = FALSE,
    size = 4.2,
    fontface = "bold"
  ) +
  
  facet_wrap(~gene, ncol = 3, scales = "fixed") +  
  
  coord_cartesian(ylim = c(global_min, star_y_pos + (global_max - global_min) * 0.05), clip = "off") +
  
  scale_fill_manual(
    values = c("Neuron" = "#E69F00",
               "Oligodendrocyte" = "#009E73")
  ) +
  labs(
    title = "Top Neuron- and Oligodendrocyte-Enriched Marker Genes",
    x = NULL,
    y = expression(bold("Expression (log"[2]*" CPM)"))
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey95", color = "black"),
    strip.text = element_text(face = "bold.italic"),
    legend.position = "none",
    panel.border = element_rect(colour = "black", fill = NA),
    axis.text.x = element_text(face = "bold", color = "black"),
    plot.title = element_text(face = "bold", size = 13)
  )



ggsave(
  file.path(OUTPUT_DIR, "08h_expression_celltype_markers.pdf"),
  p_celltype_markers,
  width = 10,
  height = 8,
  device = pdf
)

cat("  08h_expression_celltype_markers.pdf\n")


###############################################################################
### PD-Associated Gene Set: Neuron vs Oligodendrocyte
###############################################################################

cat("\n--- Creating PD-associated gene expression plot (Neuron vs Oligo) ---\n")

# PD/neurodegeneration candidate genes
pd_candidates <- c("SNCA", "PARK7", "LRRK2", "PINK1", "PRKN", "PARK2", "PARK8","GBA", "MAPT", "COQ2")

# Filter to genes present in the expression matrix
pd_candidates_present <- intersect(pd_candidates, rownames(logmat))

cat("  PD candidates present:", length(pd_candidates_present), "/",
    length(pd_candidates), "\n")
cat("  Present:", paste(pd_candidates_present, collapse = ", "), "\n")

missing_genes <- setdiff(pd_candidates, pd_candidates_present)
if (length(missing_genes) > 0) {
  cat("  Missing from expression data:", paste(missing_genes, collapse = ", "), "\n")
}

###############################################################################
### Prepare expression data
###############################################################################

expr_long_pd <- as.data.frame(t(logmat[pd_candidates_present, , drop = FALSE])) %>%
  rownames_to_column("ROI") %>%
  pivot_longer(
    cols = all_of(pd_candidates_present),
    names_to = "gene",
    values_to = "logcounts"
  ) %>%
  left_join(
    meta_df %>% dplyr::select(ROI, Cell),
    by = "ROI"
  ) %>%
  mutate(
    gene = factor(gene, levels = pd_candidates_present),
    Cell = factor(Cell, levels = c("Neuron", "Oligodendrocyte"))
  )


###############################################################################
### Statistical testing — limma/voom FDRs from de_celltype
###############################################################################

global_max <- max(expr_long_pd$logcounts, na.rm = TRUE)
global_min <- min(expr_long_pd$logcounts, na.rm = TRUE)

star_y_pos <- global_max + (global_max - global_min) * 0.08

pd_candidates_plot <- intersect(
  pd_candidates_present,
  rownames(de_celltype)
)

pd_candidates_plot

stat_results_pd <- de_celltype[pd_candidates_plot, , drop = FALSE] %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  mutate(
    p_label = case_when(
      adj.P.Val < 0.0001 ~ "****",
      adj.P.Val < 0.001  ~ "***",
      adj.P.Val < 0.01   ~ "**",
      adj.P.Val < 0.05   ~ "*",
      TRUE               ~ "ns"
    ),
    y_pos = star_y_pos
  )
###############################################################################
### Plot
###############################################################################

p_pd_markers <- ggplot(
  expr_long_pd,
  aes(x = Cell, y = logcounts, fill = Cell)
) +
  geom_violin(alpha = 0.2, color = NA) +
  geom_boxplot(width = 0.20, outlier.shape = NA, alpha = 0.7) +
  geom_beeswarm(size = 1.2, alpha = 0.6, color = "black") +
  
  # Significance stars
  geom_text(
    data = stat_results_pd,
    aes(x = 1.5, y = y_pos, label = p_label),
    inherit.aes = FALSE,
    size = 4.2,
    fontface = "bold"
  ) +
  
  facet_wrap(~gene, ncol = 6, scales = "fixed") +  
  
  
  scale_fill_manual(
    values = c("Neuron" = "#E69F00",
               "Oligodendrocyte" = "#009E73")
  ) +
  labs(
    title = "Parkinson's Disease-Associated Gene Expression: Neuron vs Oligodendrocyte",
    x = NULL,
    y = expression(bold("Expression (log"[2]*" CPM)"))
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey95", color = "black"),
    strip.text = element_text(face = "bold.italic"),
    legend.position = "none",
    panel.border = element_rect(colour = "black", fill = NA),
    axis.text.x = element_text(face = "bold", color = "black"),
    plot.title = element_text(face = "bold", size = 13)
  )

ggsave(
  file.path(OUTPUT_DIR, "08i_expression_PD_genes_celltype.pdf"),
  p_pd_markers,
  width = 12,
  height = 4,
  device = pdf
)

cat("  08i_expression_PD_genes_celltype.pdf\n")













###############################################################################
###DEG Heatmap
###############################################################################

library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(grid)
library(RColorBrewer)

cat("\n=== Creating DEG Heatmap ===\n")

cols_group <- c(PD = "#4C72B0", MSA = "#C44E52")
cols_cell  <- c(Neuron = "#E69F00", Oligodendrocyte = "#009E73")

patients <- unique(spe_qc$Patient)
cols_patient <- setNames(
  colorRampPalette(c("grey95", "grey30"))(length(patients)),
  patients
)

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
### Right annotation (Effect Size + Significance)
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
  
  column_split = spe_ruv$Group[sample_order],  # visually group MSA and PD blocks
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
pdf(file.path(OUTPUT_DIR, "09a_heatmap_top_DEGs.pdf"),
    width = 12, height = 8)

draw(ht_all,
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom",
     merge_legend = TRUE,
     legend_gap = unit(0.5, "cm"),
     padding = unit(c(1, 1, 1, 10), "mm"))

grid.text(
  "Differential gene expression in MSA vs PD substantia nigra",
  x = 0.5, y = 0.02,
  gp = gpar(fontsize = 9, fontface = "bold"),
  just = "center"
)

dev.off()

cat("✓ Saved: 09a_heatmap_top_DEGs.pdf\n")
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
  
  ### Gene annotation (logFC + significance only)
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
cat("\n=== All heatmaps complete ===\n")