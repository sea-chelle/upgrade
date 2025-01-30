# loading required packages
library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggtext)
library(MAST)
library(ComplexHeatmap)
library(circlize)

# sex differences in neuropeptide expression during development --------------------------------------------------------------------------------

# subset hypothalamic neurons in embryonic samples for sex differential analysis
embryonic <- subset(hypothalamic, subset = sample %in% c("GW7_sample1", "GW8_sample1")) # GW7 male, GW8 female


# subset hypothalamic neurons in fetal samples for sex differential analysis
fetal <- subset(hypothalamic, subset = sample %in% c("GW12_sample1", "GW12_sample2", "GW15_sample1")) # GW12 male, GW15 female

# set identities of seurat objects to sex to compare sex differences
embryonic <- SetIdent(embryonic, value = embryonic$sex)
fetal <- SetIdent(fetal, value = fetal$sex)

# sex differences in neuropeptide expression in neurons during embryonic development
embryonic_neuropeptide <- FindMarkers(embryonic,
                                       ident.1 = "male",
                                       ident.2 = "female",
                                       features = c("CRH", "OXT", "AVP", "HCRT", "TRH", "POMC", "AGRP", "NPY",
                                                    "TAC1", "SST", "KISS1", "GNRH1", "GHRH", "PMCH"),
                                       test.use = "MAST")

# sex differences in neuropeptide expression in neurons during embryonic development
fetal_neuropeptide <- FindMarkers(fetal,
                                   ident.1 = "male",
                                   ident.2 = "female",
                                   features = c("CRH", "OXT", "AVP", "HCRT", "TRH", "POMC", "AGRP", "NPY",
                                                "TAC1", "SST", "KISS1", "GNRH1", "GHRH", "PMCH"),
                                   test.use = "MAST")

# prepare data for plotting
embryonic_neuropeptide$gene <- rownames(embryonic_neuropeptide)
fetal_neuropeptide$gene <- rownames(fetal_neuropeptide)
embryonic_neuropeptide$stage <- "Embryonic"
fetal_neuropeptide$stage <- "Fetal"
developmental_neuropeptide <- bind_rows(embryonic_neuropeptide, fetal_neuropeptide)

developmental_neuropeptide <- developmental_neuropeptide %>%
  dplyr::select(gene, avg_log2FC, stage) %>%
  complete(gene, stage, fill = list(avg_log2FC = 0)) %>% 
  arrange(stage == "Embryonic", gene) %>% 
  mutate(stage = factor(stage, levels = c("Fetal", "Embryonic"))) %>%
  mutate(gene = factor(gene, levels = unique(gene))) 


# plot barplot
ggplot(developmental_neuropeptide, aes(x = gene, y = avg_log2FC, fill = stage)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5), width = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(values = c("Embryonic" = "#fc8d62", Fetal = "#66c2a5")) +  # Color by stage
  labs(x = "Gene", y = "Log2 Fold Change", fill = "Stage", title = "Male vs Female neuropeptide expression") +
  annotate("text", x = Inf, y = max(developmental_neuropeptide$avg_log2FC), label = "Upregulated", vjust = -0.5, hjust = 1.2, color = "black", size = 4) +
  annotate("text", x = Inf, y = min(developmental_neuropeptide$avg_log2FC), label = "Downregulated", vjust = 1.5, hjust = 1.2, color = "black", size = 4) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_classic() +
  theme(axis.text.x = element_text(family = "Arial", size = 11, color = "black"),
    axis.text.y = element_text(family = "Arial", size = 11, color = "black"),
    axis.title = element_text(family = "Arial", size = 12, color = "black"),
    axis.title.x = element_blank(),
    plot.title = element_text(family = "Arial", size = 16, hjust = 0.5, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(family = "Arial", size = 11, color = "black"),
    panel.grid = element_blank()) +
  coord_flip()

# sex differences in IEG expression during development ----------------------------------------

# sex differences in IEG expression in neurons during embryonic development
embryonic_IEG <- FindMarkers(embryonic,
                                      ident.1 = "male",
                                      ident.2 = "female",
                                      features = c("ARC", "DUSP1", "DUSP6", "DUSP14", "EGR1", "EGR2", "EGR3", "EGR4",
                                                   "NR4A1", "NR4A2", "NR4A3", "NPAS1", "NPAS2", "NPAS3", "NPAS4",
                                                   "JUN", "FOS", "FOSL2", "FOSB"),
                                      test.use = "MAST")

# sex differences in IEG expression in neurons during embryonic development
fetal_IEG <- FindMarkers(fetal,
                         ident.1 = "male",
                         ident.2 = "female",
                         features = c("ARC", "DUSP1", "DUSP6", "DUSP14", "EGR1", "EGR2", "EGR3", "EGR4",
                                               "NR4A1", "NR4A2", "NR4A3", "NPAS1", "NPAS2", "NPAS3", "NPAS4",
                                               "JUN", "FOS", "FOSL2", "FOSB"),
                         test.use = "MAST")

# prepare data for plotting
embryonic_IEG$gene <- rownames(embryonic_IEG)
fetal_IEG$gene <- rownames(fetal_IEG)
embryonic_IEG$stage <- "Embryonic"
fetal_IEG$stage <- "Fetal"

developmental_IEG <- bind_rows(embryonic_IEG, fetal_IEG)
developmental_IEG <- developmental_IEG%>%
  dplyr::select(gene, avg_log2FC, stage) %>%
  complete(gene, stage, fill = list(avg_log2FC = 0))

# plot barplot
ggplot(developmental_IEG, aes(x = gene, y = avg_log2FC, fill = stage)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5), width = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(values = c("Embryonic" = "#4daf4a", Fetal = "#984ea3")) +  # Color by stage
  labs(x = "Gene", y = "Log2 Fold Change", fill = "Stage", title = "Male vs Female IEG expression") +
  annotate("text", x = Inf, y = max(developmental_IEG$avg_log2FC), label = "Upregulated", vjust = -0.5, hjust = 1.2, color = "black", size = 4) +
  annotate("text", x = Inf, y = min(developmental_IEG$avg_log2FC), label = "Downregulated", vjust = 6, hjust = 1.2, color = "black", size = 4) +
  scale_y_continuous(expand = c(0, 0), limits = c(-3, 5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_bw() +
  theme(axis.text.x = element_text(family = "Arial", size = 11, color = "black"),
        axis.text.y = element_text(family = "Arial", size = 11, color = "black"),
        axis.title = element_text(family = "Arial", size = 12, color = "black"),
        axis.title.x = element_blank(),
        plot.title = element_text(family = "Arial", size = 16, hjust = 0.5, color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(family = "Arial", size = 11, color = "black"),
        panel.grid = element_blank())

# sex differences in MDD-associated gene expression during development ----------------------------------------

# sex differences in IEG expression in neurons during embryonic development
embryonic_MDD <- FindMarkers(embryonic,
                             ident.1 = "male",
                             ident.2 = "female",
                             features = c("SLC6A2", "TPH2", "FKBP5", "HTR2A", "CRHR1", "BDNF", "TSPAN5", "ERICH3", "NR3C1", 
                                          "GRM7", "PCLO", "NEGR1"),
                             test.use = "MAST")

# sex differences in IEG expression in neurons during embryonic development
fetal_MDD <- FindMarkers(fetal,
                       ident.1 = "male",
                       ident.2 = "female",
                       features = c("SLC6A2", "TPH2", "FKBP5", "HTR2A", "CRHR1", "BDNF", "TSPAN5", "ERICH3", "NR3C1", 
                                    "GRM7", "PCLO", "NEGR1"),
                       test.use = "MAST")

# prepare data for plotting
embryonic_MDD$gene <- rownames(embryonic_MDD)
fetal_MDD$gene <- rownames(fetal_MDD)
embryonic_MDD$stage <- "Embryonic"
fetal_MDD$stage <- "Fetal"

developmental_MDD <- bind_rows(embryonic_MDD, fetal_MDD)
developmental_MDD <- developmental_MDD%>%
  select(gene, avg_log2FC, stage) %>%
  complete(gene, stage, fill = list(avg_log2FC = 0))

# plot barplot
ggplot(developmental_MDD, aes(x = gene, y = avg_log2FC, fill = stage)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5), width = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(values = c("Embryonic" = "#4daf4a", Fetal = "#984ea3")) +  # Color by stage
  labs(x = "Gene", y = "Log2 Fold Change", fill = "Stage", title = "Male vs Female depression-associated gene expression") +
  annotate("text", x = Inf, y = max(developmental_MDD$avg_log2FC), label = "Upregulated", vjust = -3.5, hjust = 1.2, color = "black", size = 4) +
  annotate("text", x = Inf, y = min(developmental_MDD$avg_log2FC), label = "Downregulated", vjust = 4, hjust = 1.2, color = "black", size = 4) +
  scale_y_continuous(expand = c(0, 0), limits = c(-2.5, 2)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_bw() +
  theme(axis.text.x = element_text(family = "Arial", size = 11, color = "black"),
        axis.text.y = element_text(family = "Arial", size = 11, color = "black"),
        axis.title = element_text(family = "Arial", size = 12, color = "black"),
        axis.title.x = element_blank(),
        plot.title = element_text(family = "Arial", size = 16, hjust = 0.5, color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(family = "Arial", size = 11, color = "black"),
        panel.grid = element_blank())

# sex differences in ASD-associated gene expression during development ----------------------------------------

# sex differences in ASD expression in neurons during embryonic development
embryonic_ASD <- FindMarkers(embryonic,
                             ident.1 = "male",
                             ident.2 = "female",
                             features = c("SCN2A", "CDH8", "SHANK3", "MECP2", "PTEN", "NRXN1", "DYRK1A", "NLGN3",
                                          "FOXP1", "FOXP2", "UBE3A", "FMR1", "GRIN2B", "SYNGAP1", "ADNP"),
                             test.use = "MAST")

# sex differences in ASD expression in neurons during embryonic development
fetal_ASD <- FindMarkers(fetal,
                       ident.1 = "male",
                       ident.2 = "female",
                       features = c("SCN2A", "CDH8", "SHANK3", "MECP2", "PTEN", "NRXN1", "DYRK1A", "NLGN3",
                                          "FOXP1", "FOXP2", "UBE3A", "FMR1", "GRIN2B", "SYNGAP1", "ADNP"),
                       test.use = "MAST")

# prepare data for plotting
embryonic_ASD$gene <- rownames(embryonic_ASD)
fetal_ASD$gene <- rownames(fetal_ASD)
embryonic_ASD$stage <- "Embryonic"
fetal_ASD$stage <- "Fetal"

developmental_ASD <- bind_rows(embryonic_ASD, fetal_ASD)
developmental_ASD <- developmental_ASD %>%
  select(gene, avg_log2FC, stage) %>%
  complete(gene, stage, fill = list(avg_log2FC = 0))

# plot barplot
ggplot(developmental_ASD, aes(x = gene, y = avg_log2FC, fill = stage)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5), width = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(values = c("Embryonic" = "#4daf4a", Fetal = "#984ea3")) +  # Color by stage
  labs(x = "Gene", y = "Log2 Fold Change", fill = "Stage", title = "Male vs Female ASD-associated gene expression") +
  annotate("text", x = Inf, y = max(developmental_ASD$avg_log2FC), label = "Upregulated", vjust = -0.5, hjust = 1.2, color = "black", size = 4) +
  annotate("text", x = Inf, y = min(developmental_ASD$avg_log2FC), label = "Downregulated", vjust = 1.5, hjust = 1.2, color = "black", size = 4) +
  scale_y_continuous(expand = c(0, 0), limits = c(-4.5, 4)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_bw() +
  theme(axis.text.x = element_text(family = "Arial", size = 11, color = "black"),
        axis.text.y = element_text(family = "Arial", size = 11, color = "black"),
        axis.title = element_text(family = "Arial", size = 12, color = "black"),
        axis.title.x = element_blank(),
        plot.title = element_text(family = "Arial", size = 16, hjust = 0.5, color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(family = "Arial", size = 11, color = "black"),
        panel.grid = element_blank())

