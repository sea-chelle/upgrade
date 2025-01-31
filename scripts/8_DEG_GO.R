# loading required packages
library(Seurat)
library(tidyverse)
library(ggplot2)
library(edgeR)
library(EnhancedVolcano)
library(clusterProfiler)

# prepare data for edgeR DEG analysis --------------------------------------------------------------------------------

# aggregate counts to sample level
cts_embryonic <- (AggregateExpression(embryonic, 
                                      group.by = "sample",
                                      assays = "RNA",
                                      slot = "counts",
                                      return.seurat = F))
cts_embryonic <- as.data.frame(cts_embryonic$RNA)

cts_fetal <- (AggregateExpression(fetal, 
                                  group.by = "sample",
                                  assays = "RNA",
                                  slot = "counts",
                                  return.seurat = F))
cts_fetal <- as.data.frame(cts_fetal$RNA)

# perform DEG analysis on embryonic samples using edgeR --------------------------------------------------------------------------------

# convert counts dataframe to a DGEList object
dge_embryonic <- DGEList(counts = cts_embryonic)
dge_embryonic$samples$group <- factor(c("male", "female"))
dge_embryonic <- calcNormFactors(dge_embryonic)

# calculate CPM (counts per million)
cpm <- cpm(dge_embryonic)

# filter out genes with CPM > 1 in at least 1 samples
keep <- rowSums(cpm > 1) >= 1
dge_embryonic <- dge_embryonic[keep, ]

# assign group factors (conditions) to each sample - in this case sex
design_embryonic <- model.matrix(~group, data = dge_embryonic$samples)

# estimate dispersion by assuming a common dispersion 
dge_embryonic  <- estimateDisp(dge_embryonic)

# DEG
bcv = 0.1

gw7_vs_gw8 <- exactTest(dge_embryonic, pair =c("female", "male"), dispersion=bcv^2)
res_embryonic <- topTags(gw7_vs_gw8, n = Inf)
res_embryonic <- as.data.frame(res_embryonic)

# create column for significantly down or upregulated genes
res_embryonic <- res_embryonic %>% mutate(diffexpressed = case_when(
  logFC > 1 & FDR < 0.05 ~ "UP",
  logFC < 1 & FDR < 0.05 ~ "DOWN",
  PValue > 0.05 ~ "NO"
))

# plot volcano plot of DEG
EnhancedVolcano(res_embryonic, x = "logFC", y = "PValue", lab = rownames(res_embryonic),
                pCutoff = 0.05,
                FCcutoff = 1,
                gridlines.major = F,
                gridlines.minor = F)

# perform DEG analysis on fetal samples using edgeR ----------------------------------------

# create metadata for the samples (needed since there are more than 1 sample for some conditions)
metadata_fetal <- data.frame(Sample = c("GW12_sample1", "GW12_sample2", "GW15_sample1"),
                             Group = c("male", "male", "female"))

# convert counts dataframe to a DGEList object
dge_fetal <- DGEList(counts = cts_fetal)
dge_fetal$samples$group <- metadata_fetal$Group
dge_fetal <- calcNormFactors(dge_fetal)

# calculate CPM (counts per million)
cpm <- cpm(dge_fetal)

# filter out genes with CPM > 1 in at least 1 samples
keep <- rowSums(cpm > 1) >= 1
dge_fetal <- dge_fetal[keep, ]

# assign group factors (conditions) to each sample - in this case sex
design_fetal <- model.matrix(~group, data = dge_fetal$samples)

# estimate dispersion by assuming a common dispersion 
dge_fetal <- estimateDisp(dge_fetal)

# DEG
gw12_vs_gw15 <- exactTest(dge_fetal, pair =c("female", "male"), dispersion=bcv^2)
res_fetal <- topTags(gw12_vs_gw15, n = Inf)
res_fetal <- as.data.frame(res_fetal)

# create column for significantly down or upregulated genes
res_fetal <- res_fetal %>% mutate(diffexpressed = case_when(
  logFC > 1 & FDR < 0.05 ~ "UP",
  logFC < 1 & FDR < 0.05 ~ "DOWN",
  PValue > 0.05 ~ "NO"
))

# plot volcano plot of DEG
EnhancedVolcano(res_fetal, x = "logFC", y = "PValue", lab = rownames(res_fetal),
                pCutoff = 0.05,
                FCcutoff = 1,
                gridlines.major = F,
                gridlines.minor = F)

# GO and enrichment analysis on embryonic samples with clusterProfiler --------------------------------------------------------------------------------

# create table of significantly differentially expressed genes and find overlapping genes between time points
sigs_embryonic <- na.omit(res_embryonic)
sigs_embryonic <- sigs_embryonic[sigs_embryonic$diffexpressed != "NO",]

sigs_fetal <- na.omit(res_fetal)
sigs_fetal <- sigs_fetal[sigs_fetal$diffexpressed != "NO",]

overlap_genes <- intersect(rownames(sigs_embryonic), rownames(sigs_fetal))

# getting a list of upregulated and downregulated genes to test for GO enrichment
genes_to_test_embryonic <- split(sigs_embryonic, sigs_embryonic$diffexpressed)
upregulated_embryonic <- rownames(genes_to_test_embryonic$UP)
downregulated_embryonic <- rownames(genes_to_test_embryonic$DOWN)

# getting a list of top 500 upregulated and top 500 downregulated genes to test for GO enrichment
upregulated_embryonic <- sigs_embryonic[sigs_embryonic$diffexpressed == "UP", ]
downregulated_embryonic <- sigs_embryonic[sigs_embryonic$diffexpressed == "DOWN", ]

upregulated_embryonic <- upregulated_embryonic[order(upregulated_embryonic$PValue), ]
downregulated_embryonic<- downregulated_embryonic[order(downregulated_embryonic$PValue), ]

upregulated_embryonic <- rownames(upregulated_embryonic)
downregulated_embryonic <- rownames(downregulated_embryonic)

# running gene ontology for unregulated genes for embryonic samples
GO_upregulated_embryonic <- enrichGO(gene = upregulated_embryonic, 
                                     OrgDb = "org.Hs.eg.db",
                                     keyType = "SYMBOL",
                                     ont = "BP")

# plotting results in barplot
plot(barplot(GO_upregulated_embryonic, showCategory = 20))

# run gene ontology for downregulated genes
GO_downregulated_embryonic <- enrichGO(gene = downregulated_embryonic, 
                                       OrgDb = "org.Hs.eg.db",
                                       keyType = "SYMBOL",
                                       ont = "BP")

# plotting results in barplot
plot(barplot(GO_downregulated_embryonic, showCategory = 20))

# GO and enrichment analysis on fetal samples with clusterProfiler --------------------------------------------------------------------------------

# getting a list of upregulated and downregulated genes to test for GO enrichment
genes_to_test_fetal <- split(sigs_fetal, sigs_fetal$diffexpressed)
upregulated_fetal <- rownames(genes_to_test_fetal$UP)
downregulated_fetal <- rownames(genes_to_test_fetal$DOWN)

# getting a list of top 500 upregulated and top 500 downregulated genes to test for GO enrichment
upregulated_fetal <- sigs_fetal[sigs_fetal$diffexpressed == "UP", ]
downregulated_fetal <- sigs_fetal[sigs_fetal$diffexpressed == "DOWN", ]

upregulated_fetal<- upregulated_fetal[order(upregulated_fetal$PValue), ]
downregulated_fetal<- downregulated_fetal[order(downregulated_fetal$PValue), ]

upregulated_fetal <- rownames(upregulated_fetal)
downregulated_fetal <- rownames(downregulated_fetal)

# running gene ontology for unregulated genes for fetal samples
GO_upregulated_fetal <- enrichGO(gene = upregulated_fetal, 
                                     OrgDb = "org.Hs.eg.db",
                                     keyType = "SYMBOL",
                                     ont = "BP")

# plotting results in barplot
plot(barplot(GO_upregulated_fetal, showCategory = 20))

# run gene ontology for downregulated genes
GO_downregulated_fetal <- enrichGO(gene = downregulated_fetal, 
                                       OrgDb = "org.Hs.eg.db",
                                       keyType = "SYMBOL",
                                       ont = "BP")

# plotting results in barplot
plot(barplot(GO_downregulated_fetal, showCategory = 20))
