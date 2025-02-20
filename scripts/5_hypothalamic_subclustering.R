# load required packages
library(Seurat)
library(tidyverse)
library(ggplot2)

# subclustering hypothalamic neurons --------------------------------------------------------------------------------

# subset hypothalamic neurons
hypothalamic <- subset(neuron, idents = c("Glutamatergic", "GABAergic", "Histaminergic"))

# perform standard seurat workflow
hypothalamic <- NormalizeData(hypothalamic)
hypothalamic <- FindVariableFeatures(hypothalamic, 
                                     selection.method = "vst", 
                                     nfeatures = 2268)

top10_hypothalamic <- head(VariableFeatures(hypothalamic), 10)
plot_top10_hypothalamic <- VariableFeaturePlot(hypothalamic)
LabelPoints(plot = plot_top10_hypothalamic, points = top10_hypothalamic, repel = TRUE, xnudge = 0, ynudge = 0)

hypothalamic <- ScaleData(hypothalamic, 
                          vars.to.regress = c("S.Score", "G2M.Score", "percentage_mt"))

hypothalamic <- RunPCA(hypothalamic, features = VariableFeatures(object = hypothalamic))
ElbowPlot(hypothalamic)

# run UMAP
hypothalamic <- FindNeighbors(hypothalamic, dims = 1:35)
hypothalamic <- FindClusters(hypothalamic, resolution = 1.6)
hypothalamic <- RunUMAP(hypothalamic, dims = 1:35)

# plot UMAP to visualise clusters
DimPlot(hypothalamic, reduction = "umap", group.by = "seurat_clusters", label = T)

# hypothalamic nuclei identification using markers from paper --------------------------------------------------------------------------------------

# suprachiasmatic nucleus (SCN)
FeaturePlot(hypothalamic, 
            features = c("ONECUT1", "ONECUT2", "SIX6"))
VlnPlot(hypothalamic, 
        features = c("ONECUT1", "ONECUT2", "SIX6"), pt.size = 0)

# paraventricular nucleus/supraoptic nucleus (PVN/SON)
FeaturePlot(hypothalamic, 
            features = c("OTP", "SIM1"))
VlnPlot(hypothalamic, 
        features = c("OTP", "SIM1"), pt.size = 0)

# anterior hypohtlamaus(AH)
FeaturePlot(hypothalamic, 
            features = c("SIX3", "DLX1", "DLX2", "DLX5", "DLX6"))
VlnPlot(hypothalamic, 
        features = c("SIX3", "DLX1", "DLX2", "DLX5", "DLX6"), pt.size = 0)

# tuberal hypothalamus (TH)
FeaturePlot(hypothalamic, 
            features = c("ISL1", "MEIS2", "SST", "DLK1", "ASCL1", "TBX3"))
VlnPlot(hypothalamic, 
        features = c("ISL1", "MEIS2", "SST", "DLK1", "ASCL1", "TBX3"), pt.size = 0)

# anterior hypothalamus/dorsomedial hypothalamus (AH/DMH)
FeaturePlot(hypothalamic, 
            features = c("ACHE", "PCSK1N"))
VlnPlot(hypothalamic, 
        features = c("ACHE", "PCSK1N"), pt.size = 0)

# ventromedial hypothalamus (VMH)
FeaturePlot(hypothalamic, 
            features = c("FEZF1", "NPTX2"))
VlnPlot(hypothalamic, 
        features = c("FEZF1", "NPTX2"), pt.size = 0)

# arcuate nucelus (ARC)
FeaturePlot(hypothalamic, 
            features = c("GHRH", "GSX1", "SST", "PCP4", "OTP", "AGRP", "NPY", "NPY1R", "POMC", "CBLN1", "TBX3", "KISS1", "TAC3", "PTGER3"))
VlnPlot(hypothalamic, 
        features = c("GHRH", "GSX1", "SST", "PCP4", "OTP", "AGRP", "NPY", "NPY1R", "POMC", "CBLN1", "TBX3", "KISS1", "TAC3", "PTGER3"), pt.size = 0)

# lateral hypothalamus (LH)
FeaturePlot(hypothalamic, 
            features = c("PDYN", "HCRT"))
VlnPlot(hypothalamic, 
        features = c("PDYN", "HCRT"), pt.size = 0)

# posterior hypothalamus (PH)
FeaturePlot(hypothalamic, 
            features = c("CRABP1", "MEIS2", "CRABP2", "ONECUT1", "ONECUT2", "NHLH1", "IRX1", "IRX2", "IRX4", "LAMP5"))
VlnPlot(hypothalamic, 
        features = c("CRABP1", "MEIS2", "CRABP2", "ONECUT1", "ONECUT2", "NHLH1", "IRX1", "IRX2", "IRX4", "LAMP5"), pt.size = 0)

# premammillary nucleus (PMN)
FeaturePlot(hypothalamic, 
            features = c("ADCYAP1", "SSTR2", "ONECUT1", "PRPH"))
VlnPlot(hypothalamic, 
        features = c("ADCYAP1", "SSTR2", "ONECUT1", "PRPH"), pt.size = 0)

# tuberomammillary nucleus (TMN)
FeaturePlot(hypothalamic, 
            features = c("HDC", "SLC18A2", "GULP1", "CBLN4"))
VlnPlot(hypothalamic, 
        features = c("HDC", "SLC18A2", "GULP1", "CBLN4"), pt.size = 0)

# tuberomamillary terminal (TT)
FeaturePlot(hypothalamic, 
            features = c("SST", "LHX6", "ARX", "CRABP1", "CALB2", "NFIB"))
VlnPlot(hypothalamic, 
        features = c("SST", "LHX6", "ARX", "CRABP1", "CALB2", "NFIB"), pt.size = 0)

# supramammilary nucelus/medial mammillary nucleus (SMN/MMN)
FeaturePlot(hypothalamic, 
            features = c("ADCYAP1", "PITX2", "LHX5", "SIM1", "EMX2"))
VlnPlot(hypothalamic, 
        features = c("ADCYAP1", "PITX2", "LHX5", "SIM1", "EMX2"), pt.size = 0)

# medial mamillary nucleus (MMN)
FeaturePlot(hypothalamic, 
            features = c("CRABP1", "LHX5", "NHLH1"))
VlnPlot(hypothalamic, 
        features = c("CRABP1", "LHX5", "NHLH1"), pt.size = 0)

# supramamillary nucleus (SMN)
FeaturePlot(hypothalamic, 
            features = c("CALB2", "BARHL1", "BARHL2", "ADCYAP1", "IRX2", "IRX3", "IRX5", "NEUROD2"))
VlnPlot(hypothalamic, 
        features = c("CALB2", "BARHL1", "BARHL2", "ADCYAP1", "IRX2", "IRX3", "IRX5", "NEUROD2"), pt.size = 0)

# TE
FeaturePlot(hypothalamic, 
            features = c("TBR1", "EOMES", "EMX2", "BHLHE22"))
VlnPlot(hypothalamic, 
        features = c("TBR1", "EOMES", "EMX2", "BHLHE22"), pt.size = 0)

# prethalamus
FeaturePlot(hypothalamic, 
            features = c("SP9", "ARX"))
VlnPlot(hypothalamic, 
        features = c("SP9", "ARX"), pt.size = 0)

# thalamus
FeaturePlot(hypothalamic, 
            features = c("GBX2", "TCF7L2", "GATA2", "GATA3"))
VlnPlot(hypothalamic, 
        features = c("GBX2", "TCF7L2", "GATA2", "GATA3"), pt.size = 0)

# hypothalamic nuclei identification using differentially expressed genes --------------------------------------------------------------------------------------

# find top expressed genes from each cluster
hypothalamic_markers <- FindAllMarkers(hypothalamic, only.pos = T, 
                                       min.pct = 0.25,
                                       min.diff.pct = 0.25,
                                       logfc.threshold = 0.25)

# create function to get top 5 markers for each cluster
top10_markers_hypothalamic <- function(x){
  hypothalamic_markers[hypothalamic_markers$cluster == x, ] %>% head(n=10)}

# create a table of top 5 markers for each cluster
top10_markers_hypohtlamaic_df <- map_dfr(0:49, top10_markers_hypothalamic)

# renaming clusters --------------------------------------------------------------------------------------------

# assign names to clusters based on expression of neuropeptides and transcription factors
new_cluster_ids_hypothalamic <- c("HDC/SLC18A2/GULP1", "GBX2/TCF7L2", "PVN/SON", "GBX2/TCF7L2", "GBX2/TCF7L2", "SST/PCP4/OTP", "PITX2/LHX5/SIM1", "GBX2/TCF7L2", "GHRH/GSX1", "ADCYAP1/PITX2", "TBX3/ISL1",
                                  "MMN", "CRABP1/CALB2", "FEZF1/NPTX2", "KISS1/TAC3", "PH", "Unknown", "AH", "HDC/SLC18A2/GULP1", "Unknown", "HDC/SLC18A2/GULP1",
                                  "SST/LHX6/ARX", "SST/LHX6/ARX", "POMC/CLBN1", "AH", "PVN/SON", "CALB2/BARHL1/BARHL2", "TBX3/ISL1", "GHRH/GSX1", "ISL1/MEIS2/SST", "DMH",
                                  "HDC/SLC18A2/CBLN4", "GBX2/GATA2/GATA3", "ISL1/MEIS2/SST", "AGRP/NPY/NPY1R", "TE", "AH", "LH", "BARHL2/NEUROD2", "POMC/TBX3", "Prethalamus",
                                  "Unknown", "Prethalamus", "Prethalamus", "PVN/SON", "DLK1/ASCL1", "Unknown", "SCN", "GBX2/GATA2/GATA3", "GBX2/TCF7L2")

names(new_cluster_ids_hypothalamic) <- levels(hypothalamic)
hypothalamic <- RenameIdents(hypothalamic, new_cluster_ids_hypothalamic)

# visualise UMAP after assigning cluster names 
hypothalamic@meta.data$"markers" <- as.factor(hypothalamic@active.ident)
DimPlot(hypothalamic, reduction = "umap", group.by="markers") 

# map marker names to hypothalamic area
hypothalamic@meta.data <- hypothalamic@meta.data %>%
  mutate(nuclei = case_when(
    markers == "SCN" ~ "SCN",
    markers == "PVN/SON" ~ "PVN/SON",
    markers == "AH" ~ "AH",
    markers == "ISL1/MEIS2/SST" ~ "TH",
    markers == "DLK1/ASCL1" ~ "TH",
    markers == "TBX3/ISL1" ~ "TH",
    markers == "DMH" ~ "DMH",
    markers == "FEZF1/NPTX2" ~ "VMH",
    markers == "GHRH/GSX1" ~ "ARC",
    markers == "SST/PCP4/OTP" ~ "ARC",
    markers == "AGRP/NPY/NPY1R" ~ "ARC",
    markers == "POMC/CLBN1" ~ "ARC",
    markers == "POMC/TBX3" ~ "ARC",
    markers == "KISS1/TAC3" ~ "ARC",
    markers == "LH" ~ "LH",
    markers == "PH" ~ "PH",
    markers == "HDC/SLC18A2/GULP1" ~ "TMN",
    markers == "HDC/SLC18A2/CBLN4" ~ "TMN",
    markers == "ADCYAP1/PITX2" ~ "SMN/MMN",
    markers == "PITX2/LHX5/SIM1" ~ "SMN/MMN",
    markers == "CALB2/BARHL1/BARHL2" ~ "SMN",
    markers == "BARHL2/NEUROD2" ~ "SMN",
    markers == "TE" ~ "TE",
    markers == "Prethalamus" ~ "Prethalamus",
    markers == "MMN" ~ "MMN",
    markers == "GBX2/TCF7L2" ~ "Thalamus",
    markers == "GBX2/GATA2/GATA3" ~ "Thalamus",
    markers == "SST/LHX6/ARX" ~ "TT",
    markers == "CRABP1/CALB2" ~ "TT",
  ))

# visualise hypothalamic seurat object by hypothalamic area
DimPlot(hypothalamic, reduction = "umap", group.by= "nuclei", label = T) 

ummary_genes_nuclei <- tapply(hypothalamic$nFeature_RNA, hypothalamic$nuclei, summary)



