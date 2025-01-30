# loading required packages
library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(MAST)
library(ComplexHeatmap)
library(circlize)

# MDD-associated genes in cell types --------------------------------------------------------------------------------

# MDD-associated gene expression in major cell types
mdd_major <- FindAllMarkers(seurat, 
                            features = c("SLC6A2", "TPH2", "FKBP5", "HTR2A", "CRHR1", "BDNF", "TSPAN5", "ERICH3", "NR3C1", 
                                         "GRM7", "PCLO", "NEGR1"),
                            test.use = "MAST")

# reshape data containing log2FC values for each cluster for each gene
mdd_major_plot <- mdd_major %>%
  select(gene, cluster, avg_log2FC) %>% 
  pivot_wider(names_from = cluster, values_from = avg_log2FC) %>%
  column_to_rownames(var = "gene") %>%
  mutate(across(everything(), ~replace(., is.na(.), 0))) %>%
  as.matrix()

# create colour gradient
col_fun = colorRamp2(c(-5, 0, 5), c("#377eb8", "white", "#e41a1c"))

# plot heatmap
column_order <- c("Neuron", "NP", "NE", "Astrocyte", "OPC", "OL", "Microglia", "Ependymal", "Endothelial", "Mural", "VLMC")
row_order <- c("FKBP5", "HTR2A", "BDNF", "TSPAN5", "ERICH3", "NR3C1", "GRM7", "PCLO", "NEGR1")

mdd_major_plot <- Heatmap(mdd_major_plot,
                          name = "log2FC",
                          col = col_fun,
                          show_row_names = T,
                          show_column_names = T,
                          column_names_side = "top",
                          column_dend_side = "bottom",
                          column_order = column_order,
                          row_order = row_order,
                          row_names_gp = gpar(fontsize = 12),
                          column_names_gp = gpar(fontsize = 12),
                          width = unit(10, "cm"),
                          height = unit(8.5, "cm"))

# MDD-associated gene expression in neuronal cell types
mdd_neuron <- FindAllMarkers(neuron, 
                             features = c("SLC6A2", "TPH2", "FKBP5", "HTR2A", "CRHR1", "BDNF", "TSPAN5", "ERICH3", "NR3C1", 
                                          "GRM7", "PCLO", "NEGR1"),
                             test.use = "MAST")

# reshape data containing log2FC values for each cluster for each gene
mdd_neuron_plot <- mdd_neuron %>%
  select(gene, cluster, avg_log2FC) %>% 
  pivot_wider(names_from = cluster, values_from = avg_log2FC) %>%
  column_to_rownames(var = "gene") %>%
  mutate(across(everything(), ~replace(., is.na(.), 0))) %>%
  select(-Unknown) %>%
  as.matrix()

# plot heatmap
column_order <- c("Glutamatergic", "GABAergic", "Histaminergic", "Cholinergic")
row_order <- c("FKBP5", "HTR2A", "BDNF", "TSPAN5", "ERICH3", "NR3C1", "GRM7", "PCLO", "NEGR1")

mdd_neuron_plot <- Heatmap(mdd_neuron_plot,
                           name = "log2FC",
                           col = col_fun,
                           show_row_names = T,
                           show_column_names = T,
                           column_names_side = "top",
                           column_dend_side = "bottom",
                           column_order = column_order,
                           row_order = row_order,
                           row_names_gp = gpar(fontsize = 12),
                           column_names_gp = gpar(fontsize = 12),
                           width = unit(3.6, "cm"),
                           height = unit(8.5, "cm"))

# MDD-associated gene expression in hypothalamic neurons
hypothalamic <- SetIdent(hypothalamic, value = hypothalamic$nuclei) # set idents as nuclei not markers

mdd_hypothalamic <- FindAllMarkers(hypothalamic, 
                                   features = c("SLC6A2", "TPH2", "FKBP5", "HTR2A", "CRHR1", "BDNF", "TSPAN5", "ERICH3", "NR3C1", 
                                                "GRM7", "PCLO", "NEGR1"),
                                   test.use = "MAST")

# reshape data containing log2FC values for each cluster for each gene
mdd_hypothalamic_plot <- mdd_hypothalamic %>%
  select(gene, cluster, avg_log2FC) %>% 
  pivot_wider(names_from = cluster, values_from = avg_log2FC) %>%
  column_to_rownames(var = "gene") %>%
  mutate(across(everything(), ~replace(., is.na(.), 0))) %>%
  as.matrix()

# plot heatmap
col_fun = colorRamp2(c(-5, 0, 5), c("#377eb8", "white", "#e41a1c"))

column_order <- c("PVN/SON", "ARC", "SCN", "AH", "LH", "PH", "DMH", "VMH", "TH", "TT", "TE", "TMN", "MMN", "SMN/MMN", "SMN", "Thalamus", "Prethalamus")
row_order <- c("FKBP5", "HTR2A", "CRHR1", "BDNF", "TSPAN5", "ERICH3", "NR3C1", "GRM7", "PCLO", "NEGR1")

mdd_hypothalamic_plot <- Heatmap(mdd_hypothalamic_plot,
                                 name = "log2FC",
                                 col = col_fun,
                                 show_row_names = T,
                                 show_column_names = T,
                                 column_names_side = "top",
                                 column_dend_side = "bottom",
                                 column_order = column_order,
                                 row_order = row_order,
                                 row_names_gp = gpar(fontsize = 12),
                                 column_names_gp = gpar(fontsize = 12),
                                 width = unit(13, "cm"),
                                 height = unit(8.5, "cm"))

# ASD-associated genes in cell types --------------------------------------------------------------------------------

# ASD-associated gene expression in major cell types
asd_major <- FindAllMarkers(seurat, features = c("SCN2A", "CDH8", "SHANK3", "MECP2", "PTEN", "NRXN1", "DYRK1A", "MEF2C", "NLGN3",
                                                 "FOXP1", "FOXP2", "UBE3A", "FMR1", "GRIN2B", "SYNGAP1", "ADNP"),
                            test.use = "MAST")

# reshape data containing log2FC values for each cluster for each gene
asd_major_plot <- asd_major %>%
  select(gene, cluster, avg_log2FC) %>% 
  pivot_wider(names_from = cluster, values_from = avg_log2FC) %>%
  column_to_rownames(var = "gene") %>%
  mutate(across(everything(), ~replace(., is.na(.), 0))) %>%
  as.matrix()

# create colour gradient
col_fun = colorRamp2(c(-6, 0, 5), c("#377eb8", "white", "#e41a1c"))

# plot heatmap
column_order <- c("Neuron", "NP", "NE", "Astrocyte", "OPC", "OL", "Microglia", "Ependymal", "Endothelial", "Mural", "VLMC")
row_order <- c("SCN2A", "CDH8", "SHANK3", "MECP2", "PTEN", "NRXN1", "DYRK1A", "MEF2C", "NLGN3", "FOXP1", "FOXP2", "UBE3A", "FMR1", 
               "GRIN2B", "SYNGAP1", "ADNP")

asd_major_plot <- Heatmap(asd_major_plot,
                          name = "log2FC",
                          col = col_fun,
                          show_row_names = T,
                          show_column_names = T,
                          column_names_side = "top",
                          column_dend_side = "bottom", 
                          column_order = column_order,
                          row_order = row_order,
                          row_names_gp = gpar(fontsize = 12),
                          column_names_gp = gpar(fontsize = 12),
                          width = unit(10, "cm"),
                          height = unit(14, "cm"))

# ASD-associated gene expression in neuronal cell types
asd_neuron <- FindAllMarkers(neuron, 
                             features = c("SCN2A", "CDH8", "SHANK3", "MECP2", "PTEN", "NRXN1", "DYRK1A", "MEF2C", "NLGN3",
                                          "FOXP1", "FOXP2", "UBE3A", "FMR1", "GRIN2B", "SYNGAP1", "ADNP"),
                             test.use = "MAST")

# reshape data containing log2FC values for each cluster for each gene
asd_neuron_plot <- asd_neuron %>%
  select(gene, cluster, avg_log2FC) %>% 
  pivot_wider(names_from = cluster, values_from = avg_log2FC) %>%
  column_to_rownames(var = "gene") %>%
  mutate(across(everything(), ~replace(., is.na(.), 0))) %>%
  select(-Unknown) %>%
  as.matrix()

# create colour gradient
col_fun = colorRamp2(c(-5, 0, 5), c("#377eb8", "white", "#e41a1c"))

# plot heatmap
column_order <- c("Glutamatergic", "GABAergic", "Histaminergic", "Cholinergic")
row_order <- c("SCN2A", "CDH8", "SHANK3", "MECP2", "PTEN", "NRXN1", "MEF2C", "NLGN3", "FOXP1", "FOXP2", "UBE3A", "FMR1", 
               "GRIN2B", "SYNGAP1", "ADNP")

# plot heatmap
asd_neuron_plot <- Heatmap(asd_neuron_plot,
                           name = "log2FC",
                           col = col_fun,
                           show_row_names = T,
                           show_column_names = T,
                           column_names_side = "top",
                           column_dend_side = "bottom",
                           column_order = column_order,
                           row_order = row_order,
                           row_names_gp = gpar(fontsize = 12),
                           column_names_gp = gpar(fontsize = 12),
                           width = unit(3.7, "cm"),
                           height = unit(14, "cm"))

# ASD-associated gene expression in hypothalamic neurons
asd_hypothalamic <- FindAllMarkers(hypothalamic, 
                                   features = c("SCN2A", "CDH8", "SHANK3", "MECP2", "PTEN", "NRXN1", "DYRK1A", "MEF2C", "NLGN3",
                                                "FOXP1", "FOXP2", "UBE3A", "FMR1", "GRIN2B", "SYNGAP1", "ADNP"),
                                   test.use = "MAST")

# reshape data containing log2FC values for each cluster for each gene
asd_hypothalamic_plot <- asd_hypothalamic %>%
  select(gene, cluster, avg_log2FC) %>% 
  pivot_wider(names_from = cluster, values_from = avg_log2FC) %>%
  column_to_rownames(var = "gene") %>%
  mutate(across(everything(), ~replace(., is.na(.), 0))) %>%
  as.matrix()

# create colour gradient
col_fun = colorRamp2(c(-10, 0, 5), c("#377eb8", "white", "#e41a1c"))

# plot heatmap
column_order <- c("PVN/SON", "ARC", "SCN", "AH", "LH", "PH", "DMH", "VMH", "TH", "TT", "TE", "TMN", "MMN", "SMN/MMN", "SMN", "Thalamus", "Prethalamus")
row_order <- c("SCN2A", "CDH8", "SHANK3", "MECP2", "PTEN", "NRXN1", "DYRK1A", "MEF2C", "NLGN3", "FOXP1", "FOXP2", "UBE3A", "FMR1", 
               "GRIN2B", "SYNGAP1", "ADNP")

asd_hypothalamic_plot <- Heatmap(asd_hypothalamic_plot,
                                 name = "log2FC",
                                 col = col_fun,
                                 show_row_names = T,
                                 show_column_names = T,
                                 column_names_side = "top",
                                 column_dend_side = "bottom",
                                 column_order = column_order,
                                 row_order = row_order,
                                 row_names_gp = gpar(fontsize = 12),
                                 column_names_gp = gpar(fontsize = 12),
                                 width = unit(13, "cm"),
                                 height = unit(14, "cm"))
