# loading required packages
library(Seurat)
library(tidyverse)
library(ggplot2)
library(SingleCellExperiment)
library(scDblFinder)
library(BiocParallel)
library(harmony)

# import count matrices and converting to seurat object --------------------------------------------------------------------------------

# get data location
dirs <- list.dirs(path = "./datasets/zhou_2022/matrix/", recursive = FALSE, full.names = FALSE)

# import files and create seurat object for each folder in dirs
for(x in dirs) {
  name <- gsub("_matrix", "", x)
  
  cts <- ReadMtx(mtx = paste0("./datasets/zhou_2022/matrix/", x, "/matrix.mtx.gz"),
                 features = paste0("./datasets/zhou_2022/matrix/", x, "/features.tsv.gz"),
                 cells = paste0("./datasets/zhou_2022/matrix/", x, "/barcodes.tsv.gz"))
  
  assign(name, CreateSeuratObject(counts = cts))
}

# merge seurat objects from the same biological replicate together
GW7_sample1_whole_male <- merge(GW7_sample1_lane1_whole_male, GW7_sample1_lane2_whole_male)
GW8_sample1_whole_female <- merge(GW8_sample1_lane1_whole_female, GW8_sample1_lane2_whole_female)

# remove unmerged files from environment
rm(list = c("GW7_sample1_lane1_whole_male", "GW7_sample1_lane2_whole_male", 
            "GW8_sample1_lane1_whole_female", "GW8_sample1_lane2_whole_female"))

# merge all biological replicates into one seurat object for QC
seurat <- merge(GW10_sample1_whole_male, y = c(GW12_sample1_whole_male, GW12_sample2_whole_male,
                                               GW15_sample1_anterior_female, GW15_sample1_medial_female, 
                                               GW15_sample1_posterior_female, GW18_sample1_anterior_male, 
                                               GW18_sample1_medial_male, GW18_sample1_posterior_male, 
                                               GW18_sample2_anterior_male, GW18_sample2_medial_male, 
                                               GW18_sample2_posterior_male, GW20_sample1_anterior_male, 
                                               GW20_sample1_medial_male, GW20_sample1_posterior_male, 
                                               GW7_sample1_whole_male, GW8_sample1_whole_female),
                add.cell.ids = ls()[3:19],
                project = "zhou_2022")

# create new column to add sample info to seurat metadata
seurat$sample <- rownames(seurat@meta.data)

# split sample column into sex, condition and barcode
seurat@meta.data <- separate(seurat@meta.data,
                             col = "sample",
                             into = c("gestational.week", "sample", "hypothalamus.location", "sex"),
                             sep = "_")

seurat@meta.data$sample <- paste(seurat@meta.data$gestational.week, seurat@meta.data$sample, sep = "_")

# check how many cells in dataset before QC
seurat # 295985 cells

# QC and filtering --------------------------------------------------------------------------------

# calculate percentage mitochondiral DNA for each cell and add to metadata
seurat$percentage_mt <- PercentageFeatureSet(seurat, pattern = "^MT-")

# calculate percentage ribosomal DNA for each cell and add to metadata
seurat$percentage_rb <- PercentageFeatureSet(seurat, pattern = "^RP[SL]")

# calculate percentage haemoglobin DNA for each cell and add to metadata
seurat$percentage_hb <- PercentageFeatureSet(seurat, pattern = "HB[^(P)]")

# filter low quality cells 
# expressing more than 200 but less than 4000 genes
# expressing less than 10% mtDNA
# expressing less than 10% rbDNA
# expressing less than 5% hbDNA
seurat <- subset(seurat, 
                 subset = nFeature_RNA > 200 &
                   nFeature_RNA < 4000 & 
                   percentage_mt < 10 &
                   percentage_rb < 10 &
                   percentage_hb < 5)

# check how many cells after initial QC
seurat # 111431 cells

# create SingleCellExperiment
seurat <- JoinLayers(seurat) 
sce <- as.SingleCellExperiment(seurat)

# doublet removal using scDblFinder
sce <- scDblFinder(sce, 
                   samples = sce@int_colData@rownames,
                   clusters = FALSE)

# convert back to seurat object
doublets <- as.Seurat(sce, counts = "counts", data = NULL)

# remove doublets from original seurat object
seurat <- subset(doublets, subset = scDblFinder.class  == "singlet")

# check how many cells after doublet removal
seurat # 102881 cells

# perform standard seurat workflow and harmony batch correction --------------------------------------------------------------------------------

# perform standard seurat workflow
seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat, 
                               selection.method = "vst", 
                               nfeatures = 2268)

# plot top 10 most variable genes amongst cells
top10 <- head(VariableFeatures(seurat), 10)
plot_top10 <- VariableFeaturePlot(seurat)
LabelPoints(plot = plot_top10, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)

# calculate cell cycle scores
seurat <- CellCycleScoring(object = seurat, 
                           g2m.features = cc.genes$g2m.genes,
                           s.features = cc.genes$s.genes)

# continue with standard seurat workflow - regress out cell cycle stores and mtDNA
seurat <- ScaleData(seurat, 
                    vars.to.regress = c("S.Score", "G2M.Score", "percentage_mt"))

seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))

# run harmony for batch correction
seurat <- RunHarmony(seurat, group.by.vars = "sample")
