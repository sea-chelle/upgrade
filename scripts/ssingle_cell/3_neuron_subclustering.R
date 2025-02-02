# load required packages
library(Seurat)
library(tidyverse)
library(ggplot2)

# UMAP for subclustering of neuronal cells ----------------------------------------

# subset neurons for subclustering
neuron <- subset(seurat, idents = "Neuron")

# perform standard seurat workflow
neuron <- NormalizeData(neuron)
neuron <- FindVariableFeatures(neuron, 
                               selection.method = "vst", 
                               nfeatures = 2268)

top10_neuron <- head(VariableFeatures(neuron), 10)
plot_top10_neuron <- VariableFeaturePlot(neuron)
LabelPoints(plot = plot_top10_neuron, points = top10_neuron, repel = TRUE, xnudge = 0, ynudge = 0)

neuron <- ScaleData(neuron, 
                    vars.to.regress = c("S.Score", "G2M.Score", "percentage_mt"))

neuron <- RunPCA(neuron, features = VariableFeatures(object = neuron))
ElbowPlot(neuron)

# run UMAP
neuron <- FindNeighbors(neuron, dims = 1:50)
neuron <- FindClusters(neuron, resolution = 1.3)
neuron <- RunUMAP(neuron, dims = 1:50)

# plot clusters
DimPlot(neuron, reduction = "umap", group.by = "seurat_clusters", label = T)

# neuron cell type identification using marker expression ----------------------------------------------

# neuron makrers from paper
FeaturePlot(neuron, 
            features = c("SLC32A1", "SLC17A6", "SLC18A2", "HDC", "ACHE", "SYT1"))
VlnPlot(neuron, 
        features = c("SLC32A1", "SLC17A6", "SLC18A2", "HDC", "ACHE", "SYT1"))

# GABAergic - canonical markers
FeaturePlot(neuron, 
            features = c("SLC32A1", "SLC6A1", "GAD1", "GAD2", "GABRA1", "GABRA2"))
VlnPlot(neuron, 
        features = c("SLC32A1", "SLC6A1", "GAD1", "GAD2", "GABRA1", "GABRA2"))

# glutaminergic - canonical markers
FeaturePlot(neuron, 
            features = c("SLC17A6", "SLC17A7", "GRIN1", "GRIN2B", "GLS"))
VlnPlot(neuron, 
        features = c("SLC17A6", "SLC17A7", "GRIN1", "GRIN2B", "GLS"))

# histaminergic - from paper
FeaturePlot(neuron, 
            features = "HDC", "SLC18A2")
VlnPlot(neuron, 
        features = "HDC", "SLC18A2")

# midbrain cholinergic = from paper
FeaturePlot(neuron, 
            features = "ACHE")
VlnPlot(neuron, 
        features = "ACHE")

# neuron subtype identification using differentially expressed genes ----------------------------------------------

# find top expressed genes from each cluster
neuron_markers <- FindAllMarkers(neuron, only.pos = T, 
                                 min.pct = 0.25,
                                 min.diff.pct = 0.25,
                                 logfc.threshold = 0.25)

# create a function to get top 5 markers for each cluster
top5_markers_neuron <- function(x){
  neuron_markers[neuron_markers$cluster == x, ] %>% head(n=5)}

# create a table of top 5 markers for each cluster
top5_markers_neuron_df <- map_dfr(0:49, top5_markers_neuron)

# renaming clusters ----------------------------------------------

# assign neuronal subtype names to each cluster
new_cluster_ids_neuron <- c("Glutamatergic","Histaminergic", "Glutamatergic", "Glutamatergic", "Cholinergic", "Unknown", "GABAergic", "GABAergic", "Glutamatergic", "Glutamatergic", "Glutamatergic",
                            "Glutamatergic", "Glutamatergic", "GABAergic", "GABAergic", "Glutamatergic", "Histaminergic", "GABAergic", "GABAergic", "Histaminergic", "GABAergic",
                            "GABAergic", "Glutamatergic", "GABAergic", "Glutamatergic", "GABAergic", "GABAergic", "Glutamatergic", "Glutamatergic", "Glutamatergic", "GABAergic",
                            "Histaminergic", "GABAergic", "GABAergic", "GABAergic", "GABAergic", "Unknown", "Glutamatergic", "GABAergic", "Unknown", "Glutamatergic",
                            "Glutamatergic", "Glutamatergic", "GABAergic", "GABAergic", "GABAergic", "Glutamatergic", "GABAergic", "GABAergic", "GABAergic")

names(new_cluster_ids_neuron) <- levels(neuron)
neuron <- RenameIdents(neuron, new_cluster_ids_neuron)

# visualise clusters after cell-type identification
neuron@meta.data$"cell_type" <- as.factor(neuron@active.ident)

DimPlot(neuron, reduction = "umap", group.by="cell_type", label = TRUE) 
