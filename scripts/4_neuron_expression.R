# loading required packages
library(Seurat)
library(tidyverse)
library(ggplot2)
library(MAST)
library(ComplexHeatmap)
library(circlize)

# neuropeptides amd IEG expression in major cell types --------------------------------------------------------------------------------

# log fold change in neuropeptide expression in all neuronal cells compared to other major cell types
neuron_neuropeptide <- FindMarkers(seurat,
                                  ident.1 = "Neuron",
                                  features = c("CRH", "OXT", "AVP", "HCRT", "TRH", "POMC", "AGRP", "NPY",
                                               "TAC1", "SST", "KISS1", "GNRH1", "GHRH", "PMCH"),
                                  test.use = "MAST")

# plot barplot
neuron_neuropeptide$gene <- factor(rownames(neuron_neuropeptide), levels = c("CRH", "OXT", "AVP", "HCRT", "TRH", "POMC", "AGRP", "NPY", 
                                                                             "TAC1", "SST", "KISS1", "GNRH1", "GHRH", "PMCH"))

ggplot(neuron_neuropeptide, aes(x = gene, y = avg_log2FC)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5), width = 0.5, fill = "#377eb8") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + 
  labs(x = "Neuropeptide", y = "Log2 Fold Change", title = "Neuropeptide expression in neurons") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_classic() +
  theme(axis.text.x = element_text(family = "Arial", size = 11, color = "black"),
        axis.text.y = element_text(family = "Arial", size = 11, color = "black"),
        axis.title = element_text(family = "Arial", size = 12, color = "black"),
        axis.title.x = element_blank(),
        plot.title = element_text(family = "Arial", size = 16, hjust = 0.5, color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(family = "Arial", size = 11, color = "black"),
        panel.grid = element_blank())

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

# IEG expression in neurons 
neuron_IEG <- FindMarkers(seurat,
                          ident.1 = "Neuron",
                          features = c("ARC", "DUSP1", "DUSP6", "DUSP14", "EGR1", "EGR2", "EGR3", "EGR4",
                                                      "NR4A1", "NR4A2", "NR4A3", "NPAS1", "NPAS2", "NPAS3", "NPAS4",
                                                      "JUN", "FOS", "FOSL2", "FOSB"),
                          test.use = "MAST")

# plot barplot
ggplot(neuron_IEG, aes(x = rownames(neuron_IEG), y = avg_log2FC)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5), width = 0.5, fill = "#377eb8") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + 
  labs(x = "Neuropeptide", y = "Log2 Fold Change", title = "IEG expression in neurons") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_classic() +
  theme(axis.text.x = element_text(family = "Arial", size = 11, color = "black"),
        axis.text.y = element_text(family = "Arial", size = 11, color = "black"),
        axis.title = element_text(family = "Arial", size = 12, color = "black"),
        axis.title.x = element_blank(),
        plot.title = element_text(family = "Arial", size = 16, hjust = 0.5, color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(family = "Arial", size = 11, color = "black"),
        panel.grid = element_blank())

# neuropeptide expression in neuronal subtypes
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

# neuropeptides amd IEG expression by gestational week --------------------------------------------------------------------------------

# set identity as gestational week
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

