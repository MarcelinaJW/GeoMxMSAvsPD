# GeoMxMSAvsPD
GeoMx Spatial Transcriptomics Analysis Pipeline
Author: Marcelina Wojewska Imperial College London Dept of Brain Sciences
If using this repository, please cite: Marcelina's PhD thesis

**Overview**
This repository contains code for the analysis of NanoString GeoMx Digital Spatial Profiler (DSP) transcriptomic data. The pipeline performs quality control, data normalisation, differential expression analysis, pathway enrichment, and WGCNA of spatially resolved gene expression profiles.

**Project Aims**
Characterise spatial transcriptomic changes across PD vs MSA oligodendrocytes and neurons.
Identify differentially expressed genes.
Explore biological pathways associated with disease pathology.
Visualise spatial patterns of gene expression.
Identify co-expressed modules associated with disease, age, tau, ab .. etc.

**Samples**
Cohort: PD and MSA FFPE substantia nigra
Number of samples: 21 samples (11 MSA and 10 PD)
Number of ROIs/AOIs: 171 (3 neg)

**Platform**
NanoString GeoMx DSP with Whole Transcriptome Atlas (WTA) panel
Cell Types: TH-positive neurons and TPPP-positive oligodendrocytes

**Quality control**
The workflow integrates:
GeoMxTools
GeoMxWorkflows
standR
SpatialExperiment
edgeR
limma-voom
RUV4 batch correction

It is designed to address technical variability, ROI quality, batch effects, and pseudoreplication arising from multiple ROIs sampled from the same patient.

**1. Data Import**
Input data include:
GeoMx DCC count files
WTA PKC annotation file
Laboratory worksheet
ROI annotation spreadsheet
Sequencing QC reports

Metadata are merged to assign:
Patient
Disease group (PD/MSA)
Cell type sampled (Neuron/Oligodendrocyte)
Slide
Batch (1 or 2)
ROI area

**2. Probe-Level Quality Control**
Probe-level QC is performed using GeomxTools.
Probes are removed if they exhibit a Geometric mean probe ratio < 0.1or a Global Grubbs outlier failure in ≥20% of segments.
Remaining probes are aggregated to target-level gene counts.

**3. Segment-Level Quality Control**
ROIs are evaluated using sequencing and tissue metrics.
Aligned reads ≥ 5,000
Percent trimmed ≥ 80%
Percent stitched ≥ 80%
Percent aligned ≥ 75%
ROI area ≥ 3,000 µm²

Additional metrics:
Library size
Detection rate
Negative probe geometric mean

ROIs failing any criterion are excluded.

QC visualisations include:
Sequencing-depth distributions
Alignment-rate distributions
ROI-area distributions
Background signal distributions
PCA of quality control metrics

**4. LOQ-Based Filtering**
Background is estimated using negative control probes.
Segment filtering - ROIs are removed if fewer than 1% of genes exceed LOQ.
Gene filtering - Genes are retained if detected above LOQ in at least 5% of ROIs.

**5. Expression of Markers associated with dopaminergic neurons and oligodendrocytes**
To evaluate the cell-type specificity of the ROI the expression of Markers associated with dopaminergic neurons and oligodendrocytes, vascular cells and immune cells are calculated. Only extreme outliers removed. 

**6. Normalisation**
Library-size normalisation is performed using edgeR's and TMM

**7. Batch Correction**
Batch effects are corrected using RUV4.

**8. ROI-Level Heterogeneity Assessment**
To assess patient effects and ROI variability:
PCA clustering is performed
ROI convex hulls are generated per patient
Pairwise ROI distances are calculated
Within-patient and between-patient distances are compared

This analysis supports modelling patient-level dependence during differential expression analysis.

**Differential Expression Analysis**
Differential expression is performed using:
edgeR
limma-voom
duplicateCorrelation

**Workflow**
Counts -> filterByExpr() -> estimateDisp() -> voom() -> duplicateCorrelation() -> voom() -> lmFit() -> contrasts.fit() -> eBayes()

**Patient-level blocking**
Multiple ROIs are sampled from each patient.
To account for non-independence duplicateCorrelation is used to estimate within-patient correlation, reducing pseudoreplication.

Disease comparison (all ROIs)
~ 0 + Group + Cell + ruv_W1 + ruv_W2 + ruv_W3

Disease comparison (neurons or oligodendroytes)
~ 0 + Group + ruv_W1 + ruv_W2 + ruv_W3

Cell-type comparison
~ 0 + Cell + Group + ruv_W1 + ruv_W2 + ruv_W3

**Outputs**
Full differential expression tables
Log2 fold change
Standard error
95% confidence intervals
Raw P values
Benjamini-Hochberg FDR

**Diagnostic outputs include:**
BCV plots
Voom mean-variance trends
Empirical Bayes trends
PCA and UMAP visualisations
Volcano plots

**Key Output Objects**
spe	-> Raw SpatialExperiment
spe_qc	-> QC-filtered dataset
spe_ruv	-> Batch-corrected dataset
de_all	-> MSA vs PD (all ROIs)
de_neuron	-> MSA vs PD (neurons)
de_oligo	-> MSA vs PD (oligodendrocytes)
de_celltype	-> Neuron vs oligodendrocyte
