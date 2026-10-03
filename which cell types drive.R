###############################################################################
### Which cell type drives the All-ROIs DEGs???????????
###############################################################################
# Get the All ROIs top DEGs
driver_check <- de_all[top_genes_all, ] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::select(gene, logFC_all = logFC, adj.P.Val_all = adj.P.Val) %>%
  
  # Join neuron-only results
  left_join(
    de_neuron %>% as.data.frame() %>% tibble::rownames_to_column("gene") %>%
      dplyr::select(gene, logFC_neuron = logFC, adj.P.Val_neuron = adj.P.Val),
    by = "gene"
  ) %>%
  
  # Join oligo-only results
  left_join(
    de_oligo %>% as.data.frame() %>% tibble::rownames_to_column("gene") %>%
      dplyr::select(gene, logFC_oligo = logFC, adj.P.Val_oligo = adj.P.Val),
    by = "gene"
  ) %>%
  
  mutate(
    sig_in_neuron = !is.na(adj.P.Val_neuron) & adj.P.Val_neuron < 0.05,
    sig_in_oligo  = !is.na(adj.P.Val_oligo)  & adj.P.Val_oligo  < 0.05,
    
    driven_by = case_when(
      sig_in_neuron & sig_in_oligo   ~ "Both cell types",
      sig_in_neuron & !sig_in_oligo  ~ "Neuron-driven",
      !sig_in_neuron & sig_in_oligo  ~ "Oligodendrocyte-driven",
      TRUE                           ~ "Neither (subthreshold in both)"
    )
  ) %>%
  arrange(adj.P.Val_all)

print(driver_check)

table(driver_check$driven_by)

write.csv(driver_check,
          file.path(OUTPUT_DIR, "AllROIs_DEG_celltype_driver_check.csv"),
          row.names = FALSE)