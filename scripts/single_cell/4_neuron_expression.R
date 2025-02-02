# load required packages
library(Seurat)
library(tidyverse)
library(ggplot2)
library(MAST)
library(ComplexHeatmap)
library(circlize)

# neuron subtype expression --------------------------------------------------------------------------------

# neuropeptide expression in neuronal subtypes
neuronsubtype_neuropeptide <-  FindAllMarkers(neuron,
                                           features = c("CRH", "OXT", "AVP", "HCRT", "TRH", "POMC", "AGRP", "NPY", 
                                                        "TAC1", "SST", "KISS1", "GNRH1", "GHRH", "PMCH"),
                                           test.use = "MAST")


# reshape data containing log2FC values for each cluster for each gene
neuronsubtype_neuropeptide_plot <- neuronsubtype_neuropeptide %>%
  filter(cluster != "Unknown") %>%
  dplyr::select(avg_log2FC, cluster, gene) %>% 
  pivot_wider(names_from = cluster, values_from = avg_log2FC) %>%
  column_to_rownames(var = "gene") %>%
  mutate(across(everything(), ~replace(., is.na(.), 0))) %>%
  as.matrix()

neuronsubtype_neuropeptide_plot <- t(neuronsubtype_neuropeptide_plot)

# create colour gradient
col_fun = colorRamp2(c(-10, 0, 5), c("#377eb8", "white", "#e41a1c"))

# plot heatmap
row_order <- c("Glutamatergic", "GABAergic", "Histaminergic", "Cholinergic")
column_order <- c("OXT", "AVP", "HCRT", "TRH", "POMC", "AGRP", "NPY", "TAC1", "SST", "KISS1", "GNRH1", "GHRH", "PMCH")

neuronsubtype_neuropeptide_plot <- Heatmap(neuronsubtype_neuropeptide_plot,
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
                                           height = unit(12, "cm"))

# IEG expression in neuronal subtypes
neuronsubtype_IEG <-  FindAllMarkers(neuron, features = c("ARC", "DUSP1", "DUSP6", "DUSP14", "EGR1", "EGR2", "EGR3", "EGR4",
                                                           "NR4A1", "NR4A2", "NR4A3", "NPAS1", "NPAS2", "NPAS3", "NPAS4",
                                                           "JUN", "FOS", "FOSL2", "FOSB"),
                                              test.use = "MAST")

# create colour gradient
col_fun = colorRamp2(c(-2, 0, 2), c("#377eb8", "white", "#e41a1c"))

# reshape data containing log2FC values for each cluster for each gene
neuronsubtype_IEG_plot <- neuronsubtype_IEG %>%
  filter(cluster != "Unknown") %>%
  dplyr::select(avg_log2FC, cluster, gene) %>% 
  pivot_wider(names_from = cluster, values_from = avg_log2FC) %>%
  column_to_rownames(var = "gene") %>%
  mutate(across(everything(), ~replace(., is.na(.), 0))) %>%
  as.matrix()

neuronsubtype_IEG_plot <- t(neuronsubtype_IEG_plot)

# plot heatmap
column_order <- c("FOS", "FOSB", "JUN", "ARC", "DUSP1", "DUSP6", "DUSP14", "EGR1", "NR4A1", "NR4A2", "NR4A3", 
               "NPAS1", "NPAS2", "NPAS3", "NPAS4")

neuronsubtype_IEG_plot <- Heatmap(neuronsubtype_IEG_plot,
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
                                           height = unit(16, "cm"))

# neuropeptides amd IEG expression by gestational week in neurons --------------------------------------------------------------------------------

# set identity as gestational week for DEG analysis
neuron_GW <- SetIdent(neuron, value = neuron$gestational.week)

neuron_neuropeptide_GW <- FindAllMarkers(neuron_GW, features = c("CRH", "OXT", "AVP", "HCRT", "TRH", "POMC", "AGRP", "NPY", 
                                                        "TAC1", "SST", "KISS1", "GNRH1", "GHRH", "PMCH"),
                                         test.use = "MAST")

# reshape data containing log2FC values for each cluster for each gene
neuron_neuropeptide_GW_plot <- neuron_neuropeptide_GW %>%
  dplyr::select(gene, cluster, avg_log2FC) %>% 
  pivot_wider(names_from = cluster, values_from = avg_log2FC) %>%
  column_to_rownames(var = "gene") %>%
  mutate(across(everything(), ~replace(., is.na(.), 0))) %>%
  as.matrix()
neuron_neuropeptide_GW_plot <- t(neuron_neuropeptide_GW_plot)

# create colour gradient
col_fun = colorRamp2(c(-10, 0, 5), c("#377eb8", "white", "#e41a1c"))

# plot heatmap
row_order <- c("GW7", "GW8", "GW10", "GW12", "GW15", "GW18", "GW20")
column_order <- c("CRH", "OXT", "AVP", "HCRT", "TRH", "POMC", "AGRP", "NPY", "TAC1", "SST", "KISS1", "GNRH1", "GHRH", "PMCH")

neuron_neuropeptide_GW_plot <- Heatmap(neuron_neuropeptide_GW_plot,
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
                                  width = unit(6, "cm"),
                                  height = unit(12, "cm"))

# expression of IEGs in neurons by gestational week
neuron_IEG_GW <- FindAllMarkers(neuron_GW, features = c("ARC", "DUSP1", "DUSP6", "DUSP14", "EGR1", "EGR2", "EGR3", "EGR4",
                                                        "NR4A1", "NR4A2", "NR4A3", "NPAS1", "NPAS2", "NPAS3", "NPAS4",
                                                        "JUN", "FOS", "FOSL2", "FOSB"),
                                test.use = "MAST")

# reshape data containing log2FC values for each cluster for each gene
neuron_IEG_GW_plot <- neuron_IEG_GW %>%
  dplyr::select(gene, cluster, avg_log2FC) %>% 
  pivot_wider(names_from = cluster, values_from = avg_log2FC) %>%
  column_to_rownames(var = "gene") %>%
  mutate(across(everything(), ~replace(., is.na(.), 0))) %>%
  as.matrix()
neuron_IEG_GW_plot <- t(neuron_IEG_GW_plot)

# set colour gradient
col_fun = colorRamp2(c(-5, 0, 5), c("#377eb8", "white", "#e41a1c"))

# plot heatmap
row_order <- c("GW7", "GW8", "GW10", "GW12", "GW15", "GW18", "GW20")
column_order <- c("FOS", "FOSB", "JUN", "ARC", "DUSP1", "DUSP6", "DUSP14", "EGR1", "EGR2", "EGR3", "NR4A1", "NR4A2", "NR4A3",
               "NPAS1", "NPAS2", "NPAS3", "NPAS4")

neuron_IEG_GW_plot <- Heatmap(neuron_IEG_GW_plot,
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
                                  width = unit(6, "cm"),
                                  height = unit(15, "cm"))

