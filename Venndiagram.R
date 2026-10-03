################################################################################
### Venn Diagrams, DEG Overlap in  MSA vs PD
################################################################################
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

cat("\n=== Venn Diagram Analysis: DEG Overlap ===\n")
cat(rep("=", 50), "\n", sep = "")

################################################################################
### Configuration
################################################################################

padj_cut <- 0.05
fc_cut   <- 0.5
out_tag  <- paste0("padj", padj_cut, "_fc", fc_cut)

COL_MSA  <- "#C44E52"   
COL_PD   <- "#4C72B0"   

################################################################################
### Helper functions — DEG sets
################################################################################

strip_version <- function(x) sub("\\.\\d+$", "", as.character(x))

get_deg_sets <- function(df,
                         padj_col  = "adj.P.Val",
                         logfc_col = "logFC",
                         padj_cut  = 0.05,
                         fc_cut    = 0.5) {
  
  stopifnot(all(c(padj_col, logfc_col) %in% colnames(df)))
  
  g  <- strip_version(rownames(df))
  ok <- !is.na(df[[padj_col]]) & !is.na(df[[logfc_col]]) & !is.na(g)
  df <- df[ok, , drop = FALSE]
  g  <- g[ok]
  
  list(
    MSA_up = unique(g[df[[padj_col]] < padj_cut & df[[logfc_col]] >=  fc_cut]),
    PD_up  = unique(g[df[[padj_col]] < padj_cut & df[[logfc_col]] <= -fc_cut])
  )
}

################################################################################
### Build DEG sets
################################################################################

deg_neuron <- get_deg_sets(de_neuron, padj_cut = padj_cut, fc_cut = fc_cut)
deg_oligo  <- get_deg_sets(de_oligo,  padj_cut = padj_cut, fc_cut = fc_cut)

cat("\nDEG set sizes:\n")
cat("  MSA-up   Neuron:          ", length(deg_neuron$MSA_up), "\n")
cat("  PD-up    Neuron:          ", length(deg_neuron$PD_up),  "\n")
cat("  MSA-up   Oligodendrocyte: ", length(deg_oligo$MSA_up),  "\n")
cat("  PD-up    Oligodendrocyte: ", length(deg_oligo$PD_up),   "\n")

sets_msa_up <- list(
  "Neuron"          = deg_neuron$MSA_up,
  "Oligodendrocyte" = deg_oligo$MSA_up
)

sets_pd_up <- list(
  "Neuron"          = deg_neuron$PD_up,
  "Oligodendrocyte" = deg_oligo$PD_up
)

################################################################################
### Two-circle geometry
################################################################################

# ── Circle outline points ────────────────────────────────────────────────────
circle_points <- function(cx, cy, r, n = 200) {
  theta <- seq(0, 2 * pi, length.out = n)
  data.frame(X = cx + r * cos(theta), Y = cy + r * sin(theta))
}

# ── Intersection points of two equal circles ─────────────────────────────────
circle_intersections <- function(cx1, cy1, r1, cx2, cy2, r2) {
  d <- sqrt((cx2 - cx1)^2 + (cy2 - cy1)^2)
  a <- (r1^2 - r2^2 + d^2) / (2 * d)
  h <- sqrt(r1^2 - a^2)
  
  xm <- cx1 + a * (cx2 - cx1) / d
  ym <- cy1 + a * (cy2 - cy1) / d
  
  x3_1 <- xm + h * (cy2 - cy1) / d
  y3_1 <- ym - h * (cx2 - cx1) / d
  x3_2 <- xm - h * (cy2 - cy1) / d
  y3_2 <- ym + h * (cx2 - cx1) / d
  
  list(p1 = c(x3_1, y3_1), p2 = c(x3_2, y3_2))
}

# ── Lens (intersection)
lens_polygon <- function(cx1, cy1, r1, cx2, cy2, r2, n = 100) {
  
  ints <- circle_intersections(cx1, cy1, r1, cx2, cy2, r2)
  p1 <- ints$p1
  p2 <- ints$p2
  
  ang1_a <- atan2(p1[2] - cy1, p1[1] - cx1)
  ang1_b <- atan2(p2[2] - cy1, p2[1] - cx1)
  ang2_a <- atan2(p1[2] - cy2, p1[1] - cx2)
  ang2_b <- atan2(p2[2] - cy2, p2[1] - cx2)
  
  mid_x <- (cx1 + cx2) / 2
  
  arc1 <- seq(ang1_a, ang1_b, length.out = n)
  if (abs(diff(range(arc1))) > pi) arc1 <- seq(ang1_a, ang1_b + 2*pi, length.out = n)
  pts1 <- data.frame(X = cx1 + r1 * cos(arc1), Y = cy1 + r1 * sin(arc1))
  # Keep the arc that is on the circle-2 side (closer to mid_x, toward circle 2)
  if (mean(pts1$X) < mid_x - 0.001 && cx2 > cx1) {
    arc1 <- seq(ang1_b, ang1_a, length.out = n)
    if (abs(diff(range(arc1))) > pi) arc1 <- seq(ang1_b, ang1_a + 2*pi, length.out = n)
    pts1 <- data.frame(X = cx1 + r1 * cos(arc1), Y = cy1 + r1 * sin(arc1))
  }
  
  arc2 <- seq(ang2_b, ang2_a, length.out = n)
  if (abs(diff(range(arc2))) > pi) arc2 <- seq(ang2_b, ang2_a + 2*pi, length.out = n)
  pts2 <- data.frame(X = cx2 + r2 * cos(arc2), Y = cy2 + r2 * sin(arc2))
  if (mean(pts2$X) > mid_x + 0.001 && cx1 < cx2) {
    arc2 <- seq(ang2_a, ang2_b, length.out = n)
    if (abs(diff(range(arc2))) > pi) arc2 <- seq(ang2_a, ang2_b + 2*pi, length.out = n)
    pts2 <- data.frame(X = cx2 + r2 * cos(arc2), Y = cy2 + r2 * sin(arc2))
  }
  
  rbind(pts1, pts2)
}

################################################################################
### theme
################################################################################

theme_venn <- function() {
  theme_void(base_size = 11, base_family = "sans") +
    theme(
      plot.title = element_text(
        size = 12, face = "bold", hjust = 0.5, color = "black",
        margin = margin(b = 3)
      ),
      plot.subtitle = element_text(
        size = 9.5, hjust = 0.5, color = "grey30",
        margin = margin(b = 12)
      ),
      legend.position = "none",
      plot.margin = margin(12, 14, 10, 14),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

################################################################################
### Core Venn drawing function
################################################################################

draw_venn <- function(sets,
                             accent_col,
                             plot_title,
                             plot_subtitle,
                             out_path,
                             fig_width  = 5.5,
                             fig_height = 5.5) {
  
  set_names <- names(sets)
  n1 <- length(sets[[1]])
  n2 <- length(sets[[2]])
  n_overlap <- length(intersect(sets[[1]], sets[[2]]))
  n1_only <- n1 - n_overlap
  n2_only <- n2 - n_overlap
  total_union <- length(unique(c(sets[[1]], sets[[2]])))
  
  # ── Fixed geometry with two equal circles, horizontal ───────────────────
  r  <- 1.3
  d  <- 1.5   # distance between centers
  cx1 <- -d / 2; cy1 <- 0
  cx2 <-  d / 2; cy2 <- 0
  
  circ1 <- circle_points(cx1, cy1, r)
  circ2 <- circle_points(cx2, cy2, r)
  lens  <- lens_polygon(cx1, cy1, r, cx2, cy2, r)
  
  fill_pale  <- scales::alpha(accent_col, 0.22)
  fill_solid <- accent_col
  
  p <- ggplot() +
    
    # Pale unique-region fills (full circles, pale — overlap redrawn solid on top)
    geom_polygon(data = circ1, aes(X, Y), fill = fill_pale, color = NA) +
    geom_polygon(data = circ2, aes(X, Y), fill = fill_pale, color = NA) +
    
    # Solid overlap region on top
    geom_polygon(data = lens, aes(X, Y), fill = fill_solid, color = NA) +
    
    # Crisp black circle outlines on top of everything
    geom_path(data = circ1, aes(X, Y), color = "black", linewidth = 0.7) +
    geom_path(data = circ2, aes(X, Y), color = "black", linewidth = 0.7) +
    
    # Set name labels (above each circle)
    annotate("text", x = cx1, y = r + 0.28, label = set_names[1],
             fontface = "bold", size = 4.2, color = "grey15") +
    annotate("text", x = cx2, y = r + 0.28, label = set_names[2],
             fontface = "bold", size = 4.2, color = "grey15") +
    
    # Unique-region counts
    annotate("text", x = cx1 - 0.55, y = 0,
             label = paste0(n1_only, "\n(", round(100*n1_only/total_union,1), "%)"),
             size = 3.6, fontface = "bold", color = "grey10", lineheight = 0.85) +
    annotate("text", x = cx2 + 0.55, y = 0,
             label = paste0(n2_only, "\n(", round(100*n2_only/total_union,1), "%)"),
             size = 3.6, fontface = "bold", color = "grey10", lineheight = 0.85) +
    
    # Overlap count (white text on solid fill)
    annotate("text", x = 0, y = 0,
             label = paste0(n_overlap, "\n(", round(100*n_overlap/total_union,1), "%)"),
             size = 3.6, fontface = "bold", color = "white", lineheight = 0.85) +
    
    coord_equal(clip = "off", xlim = c(-2.6, 2.6), ylim = c(-1.8, 2.0)) +
    
    labs(title = plot_title, subtitle = plot_subtitle) +
    
    theme_venn()
  
  ggsave(out_path, p, width = fig_width, height = fig_height, device = pdf)
  
  png_path <- sub("\\.pdf$", ".png", out_path)
  ggsave(png_path, p, width = fig_width, height = fig_height,
         dpi = 600, bg = "white")
  
  cat("  Saved:", basename(out_path), "\n")
  cat("  Saved:", basename(png_path), "\n")
  
  invisible(p)
}


cat("\n--- Rendering Venn diagrams ---\n")

p_msa <- draw_venn(
  sets          = sets_msa_up,
  accent_col    = COL_MSA,
  plot_title    = "MSA-enriched DEGs",
  plot_subtitle = paste0(
    "Upregulated in MSA vs PD | FDR < ", padj_cut,
    " | log\u2082FC \u2265 ", fc_cut
  ),
  out_path      = file.path(OUTPUT_DIR,
                            paste0("12a_Venn_", out_tag,
                                   "_MSA_enriched.pdf"))
)

p_pd <- draw_venn(
  sets          = sets_pd_up,
  accent_col    = COL_PD,
  plot_title    = "PD-enriched DEGs",
  plot_subtitle = paste0(
    "Upregulated in PD vs MSA | FDR < ", padj_cut,
    " | log\u2082FC \u2265 ", fc_cut
  ),
  out_path      = file.path(OUTPUT_DIR,
                            paste0("12b_Venn_", out_tag,
                                   "_PD_enriched.pdf"))
)

################################################################################
### Combined figure
################################################################################

p_combined <- (p_msa | p_pd) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(plot.tag = element_text(size = 13, face = "bold"))
  )

ggsave(
  file.path(OUTPUT_DIR, paste0("12c_Venn_", out_tag, "_combined.pdf")),
  p_combined, width = 11, height = 5.5, device = pdf
)

ggsave(
  file.path(OUTPUT_DIR, paste0("12c_Venn_", out_tag, "_combined.png")),
  p_combined, width = 11, height = 5.5, dpi = 600, bg = "white"
)

cat("  Saved: 12c_Venn_combined.pdf/.png\n")

################################################################################
### Export overlap tables
################################################################################

cat("\n--- Exporting tables ---\n")

shared_msa <- intersect(deg_neuron$MSA_up, deg_oligo$MSA_up)
shared_pd  <- intersect(deg_neuron$PD_up,  deg_oligo$PD_up)

cat("  Shared MSA-enriched (both cell types):", length(shared_msa), "\n")
cat("  Shared PD-enriched  (both cell types):", length(shared_pd),  "\n")

writeLines(shared_msa,
           file.path(OUTPUT_DIR,
                     paste0("12_shared_MSA_enriched_", out_tag, ".txt")))
writeLines(shared_pd,
           file.path(OUTPUT_DIR,
                     paste0("12_shared_PD_enriched_", out_tag, ".txt")))

build_shared_table <- function(genes, de_n, de_o, direction) {
  if (length(genes) == 0) {
    cat("  No shared genes for:", direction, "\n")
    return(data.frame())
  }
  pull_cols <- function(de, genes, cols) {
    g   <- strip_version(rownames(de))
    idx <- match(genes, g)
    de[idx, cols, drop = FALSE]
  }
  n_vals <- pull_cols(de_n, genes, c("logFC", "adj.P.Val"))
  o_vals <- pull_cols(de_o, genes, c("logFC", "adj.P.Val"))
  data.frame(
    Gene          = genes,
    Direction     = direction,
    logFC_Neuron  = round(n_vals$logFC,      3),
    FDR_Neuron    = signif(n_vals$adj.P.Val, 3),
    logFC_Oligo   = round(o_vals$logFC,      3),
    FDR_Oligo     = signif(o_vals$adj.P.Val, 3),
    stringsAsFactors = FALSE
  ) %>% arrange(FDR_Neuron)
}

shared_table <- bind_rows(
  build_shared_table(shared_msa, de_neuron, de_oligo, "MSA-enriched"),
  build_shared_table(shared_pd,  de_neuron, de_oligo, "PD-enriched")
)

write.csv(
  shared_table,
  file.path(OUTPUT_DIR, paste0("12_shared_DEGs_annotated_", out_tag, ".csv")),
  row.names = FALSE
)

cat("\n=== Venn Diagram Analysis Complete ===\n")