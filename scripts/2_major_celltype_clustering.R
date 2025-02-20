# load required packages
library(Seurat)
library(tidyverse)
library(ggplot2)

# UMAP dimensionality reduction for clustering --------------------------------------------------------------------------------

# run UMAP
seurat <- FindNeighbors(seurat, reduction = "harmony", dims = 1:50)
seurat <- FindClusters(seurat, resolution = 1.21)
seurat <- RunUMAP(seurat, reduction = "harmony", dims = 1:50)

# plot UMAP for cluster visualisation
DimPlot(seurat, reduction = "umap", group.by = "cell_type", label = T) 

# look at how many cells are in each cluster
table(seurat@meta.data$seurat_clusters) 

# major cell type identification using markers from paper --------------------------------------------------------------------------------

# neuroepithelial (NE)
FeaturePlot(seurat, 
            features = c("VIM", "HES1", "HMGA2", "ARHGAP28"))
VlnPlot(seurat, 
        features = c("VIM", "HES1", "HMGA2", "ARHGAP28"), pt.size = 0)

# neural progenitor (NP)
FeaturePlot(seurat, 
            features = c("NES", "ASCL1", "NBPF10", "NBPF15", "PLAGL1", "SFTA3"))
VlnPlot(seurat, 
        features = c("NES", "ASCL1", "NBPF10", "NBPF15", "PLAGL1", "SFTA3"), pt.size = 0)

# neuron
FeaturePlot(seurat, 
            features = c("STMN2", "SYT1"))
VlnPlot(seurat, 
        features = c("STMN2", "SYT1"), pt.size = 0)

# oligoprogenitor cell (OPC)
FeaturePlot(seurat, 
            features = c("PDGFRA", "APOD", "GPR75-ASB3"))
VlnPlot(seurat, 
        features = c("PDGFRA", "APOD", "GPR75-ASB3"), pt.size = 0)

# oligodendrocyte (OL)
FeaturePlot(seurat, 
            features = c("MBP", "MOG"))
VlnPlot(seurat, 
        features = c("MBP", "MOG"), pt.size = 0)

# astrocytes
FeaturePlot(seurat, 
            features = c("AQP4", "AGT", "FAM107A"))
VlnPlot(seurat, 
        features = c("AQP4", "AGT", "FAM107A"), pt.size = 0)

# ependymal
FeaturePlot(seurat, 
            features = c("CCDC153", "CCDC74A", "CCDC74B", "LRTOMT"))
VlnPlot(seurat, 
        features = c("CCDC153", "CCDC74A", "CCDC74B", "LRTOMT"), pt.size = 0)

# microglia
FeaturePlot(seurat, 
            features = c("AIF1", "CX3CR1"))
VlnPlot(seurat, 
        features = c("AIF1", "CX3CR1"), pt.size = 0)

# endothelial cells
FeaturePlot(seurat, 
            features = c("CLDN5", "FLT1", "SLC38A5"))
VlnPlot(seurat, 
        features = c("CLDN5", "FLT1", "SLC38A5"), pt.size = 0)

# mural
FeaturePlot(seurat, 
            features = c("NDUFA4L2", "PDGFRB"))
VlnPlot(seurat, 
        features = c("NDUFA4L2", "PDGFRB"), pt.size = 0)

# vascular and leptomeningeal cell (VLMC)
FeaturePlot(seurat, 
            features = c("PTGDS", "COL1A1"))
VlnPlot(seurat, 
        features = c("PTGDS", "COL1A1"), pt.size = 0)

# major cell type identification using canonical markers --------------------------------------------------------------------------------------

# NE
FeaturePlot(seurat, 
            features = c("NES", "SOX1", "SOX2", "PAX6", "HES1", "HES2", "NOTCH1", "OCLN", "CDH1", "VIM", "ZIC1", "ZIC3"))
VlnPlot(seurat, features = c("NES", "SOX1", "SOX2", "PAX6", "HES1", "HES2", "NOTCH1", "OCLN", "CDH1", "VIM", "ZIC1", "ZIC3"), pt.size = 0)
            
# NP
FeaturePlot(seurat, 
            features = c("NES", "SOX2", "PAX6", "MSI1", "PROM1", "DCX", "ASCL1"))
VlnPlot(seurat, 
        features = c("NES", "SOX2", "PAX6", "MSI1", "PROM1", "DCX", "ASCL1"), pt.size = 0)
            
# neuron
FeaturePlot(seurat, 
            features = c("TUBB3", "MAP2", "RBFOX3", "NEFL", "NEFM", "NEFH", "SYT1", "SYT2", "DLG1"))
VlnPlot(seurat, 
        features = c("TUBB3", "MAP2", "RBFOX3", "NEFL", "NEFM", "NEFH", "SYT1", "SYT2", "DLG1"), pt.size = 0)
                    
# OPC
FeaturePlot(seurat, 
            features = c("PDGFRA", "CSPG4", "OLIG2", "SOX10"))
VlnPlot(seurat,
        features = c("PDGFRA", "CSPG4", "OLIG2", "SOX10"), pt.size = 0)
                    
# OL
FeaturePlot(seurat, 
            features = c("OLIG1", "OLIG2", "OLIG3", "CLDN11", "MBP", "MOG", "SOX10"))
VlnPlot(seurat, 
        features = c("OLIG1", "OLIG2", "OLIG3", "CLDN11", "MBP", "MOG", "SOX10"), pt.size = 0)
                    
# astrocytes
FeaturePlot(seurat, 
            features = c("ALDH1L1", "AQP4", "GFAP", "S100B", "SLC1A2"))
VlnPlot(seurat, 
        features = c("ALDH1L1", "AQP4", "GFAP", "S100B", "SLC1A2"), pt.size = 0)
                    
                    
# microglia
FeaturePlot(seurat, 
            features = c("TMEM119", "TREM2", "ITGAM", "PTPRC", "CX3CR1", "P2RY12"))
VlnPlot(seurat, 
        features = c("TMEM119", "TREM2", "ITGAM", "PTPRC", "CX3CR1", "P2RY12"), pt.size = 0)
                
                    
# cell type identification using top DEGs --------------------------------------------------------------------------------------
                    
# find top expressed genes for each cluster
major_markers <- FindAllMarkers(seurat, only.pos = T, 
                                min.pct = 0.25,
                                min.diff.pct = 0.25,
                                logfc.threshold = 0.25)
                    
# create function to get top 5 markers for each cluster
top5_markers <- function(x){
major_markers[major_markers$cluster == x, ] %>% head(n=5)}
                    
# create a table of top 5 markers for each cluster
top5_markers_df <- map_dfr(0:41, top5_markers)
                    
# renaming clusters --------------------------------------------------------------------------------------
                    
# assign cell-type names to clusters
new_cluster_ids <- c("Astrocyte", "OPC", "Astrocyte", "NE", "OPC", "Neuron", "Neuron", "Neuron", "Neuron", "Microglia", "NP",
                     "Astrocyte", "Astrocyte", "Neuron", "Neuron", "NE", "OL", "Neuron", "Neuron", "Neuron", "NE",
                     "NE", "OL", "Ependymal", "OPC", "Neuron", "VLMC", "Neuron", "Mural", "NP", "NE",
                     "NP", "Endothelial", "Microglia", "Neuron", "Neuron", "NE", "Microglia", "Neuron", "Neuron", "Mural",
                     "Microglia")
                    
names(new_cluster_ids) <- levels(seurat)
seurat <- RenameIdents(seurat, new_cluster_ids)
                    
# visualise UMAP after cell-type identification
seurat@meta.data$"cell_type" <- as.factor(seurat@active.ident)
DimPlot(seurat, reduction = "umap", group.by="cell_type", label = TRUE) 
