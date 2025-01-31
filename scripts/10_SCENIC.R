# loading required packages
library(Seurat)
library(tidyverse)
library(ggplot2)
library(SCENIC)
library(AUCell)
library(RcisTarget)
library(GENIE3)
library(doRNG)
library(R2HTML)
library(pheatmap)
library(RColorBrewer)


# cell inforation and phenotype data ----------------------------------------

# subset hypothalamic neurons to only include age-matched sex difference samples
development <- subset(hypothalamic, subset = sample %in% c("GW7_sample1", "GW8_sample1",
                                                           "GW12_sample1", "GW12_sample2", "GW15_sample1"))
# download human-specific database for RcisTarget (motif ranking)
dbFiles <- c("https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.feather",
             "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.feather")

# extract raw counts as matrix from seurat object
exprMat <- as.matrix(GetAssayData(object = development, assay = "RNA", layer = "counts"))

# get metadata from seurat object
metadata_hypothalamic <- hypothalamic@meta.data

# get cell type information fpr each cell
development <- SetIdent(development, value = "nuclei")
cellInfo <- data.frame(CellType = Idents(development))

# get marker information for each cell
development <- SetIdent(development, value = "markers")
markerInfo <- data.frame(Marker = Idents(development))

# get sex for each cell
development <- SetIdent(development, value = "sex")
sexInfo <- data.frame(Sex = Idents(development))

# get development info for each cell
development@meta.data <- development@meta.data %>%
  mutate(stage = case_when(sample %in% c("GW7_sample1", "GW8_sample1") ~ "embryonic",
                           sample %in% c("GW12_sample1", "GW12_sample2", "GW15_sample1") ~ "fetal"))

development <- SetIdent(development, value = "stage")
stageInfo <- data.frame(Stage = Idents(development))

# removing cells labelled as N/A for downstream analysis
metadata_hypothalamic <- na.omit(metadata_hypothalamic)
cellInfo <- na.omit(cellInfo)
exprMat <- exprMat[, colnames(exprMat) %in% rownames(cellInfo)]

# colour to assign to the variables (each cell type)
colVars <- list(CellType=c("ARC" = "#9367bb",
                           "AH" = "#d42a2d",
                           "DMH" = "navy",
                           "LH" = "#fe802b",
                           "MMN" = "#1ebecc",
                           "TE" = "#7f7f7f",
                           "TH" = "#bcbd3c",
                           "TT" = "#1b9e77",
                           "PH" = "#2177b2",
                           "Prethalamus" = "#8c564c",
                           "PVN/SON" = "#30a139",
                           "SCN" = "#e278c0",
                           "SMN" = "#ffff33",
                           "SMN/MMN" = "#fb9a99",
                           "Thalamus" = "#e7298a",
                           "TMN" = "lightcyan"))

colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$CellType)]

# export exprMat, cellInfo and colVars
saveRDS(exprMat, file = "/Users/michellelam/int/exprMat.rds")
saveRDS(cellInfo, file = "/Users/michellelam/SCENIC/cellInfo.rds")
saveRDS(colVars, file= "/Users/michellelam/SCENIC/colVars.rds")

# plot legend 
plot.new()
legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))

# initialize scenic settings ----------------------------------------

# specify organism
org <- "hgnc" # for mouse, use mgi

# get database directory location - download from https://resources.aertslab.org/cistarget/databases/old/
dbDir = "/Users/michellelam/int"

# choose a name for the analysis
myDatasetTitle <- "SCENIC on the developing human hypothalamus"

# load scenic data
data("defaultDbNames")

# get database
dbs <- defaultDbNames[[org]]

# load in the motifannotation
data(list="motifAnnotations_hgnc_v9", package="RcisTarget")

# rename the motif annnotion by attributing it to the variable that is in the error
motifAnnotations_hgnc <- motifAnnotations_hgnc_v9

# initialise scenic
scenicOptions <- initializeScenic(org = org, 
                                  dbDir = dbDir, 
                                  dbs = dbs, 
                                  datasetTitle = myDatasetTitle,
                                  nCores = 8)

# add cellInfo and colVars into scenic options
scenicOptions@inputDatasetInfo$cellInfo <- "/Users/michellelam/int/cellInfo.rds"
scenicOptions@inputDatasetInfo$colVars <- "/Users/michellelam/int/colVars.rds"

# co-expression network ----------------------------------------

# filter out genes which are expressed at low levels or too few cells or those not available in RcisTarget database
genesKept <- geneFiltering(exprMat, scenicOptions= scenicOptions,
                           minCountsPerGene= 3*.01*ncol(exprMat),
                           minSamples= ncol(exprMat)*.01)
# 9815 genes remaning

# check if there are any genes of interest (neuropeptides, IEGs, mental health-associated genes or hypothalamic TFs) removed from filtering
interestingGenes <- c("CRH", "OXT", "AVP", "HCRT", "TRH", "POMC", "AGRP", "NPY",
                      "TAC1", "SST", "KISS1", "GNRH1", "GHRH", "PMCH", "ARC", "DUSP1", "DUSP6", "DUSP14", "EGR1", "EGR2", "EGR3", "EGR4",
                      "NR4A1", "NR4A2", "NR4A3", "NPAS1", "NPAS2", "NPAS3", "NPAS4",
                      "JUN", "FOS", "FOSL2", "FOSB", "SLC6A2", "TPH2", "FKBP5", "HTR2A", "CRHR1", "BDNF", "TSPAN5", "ERICH3", "NR3C1", 
                      "GRM7", "PCLO", "NEGR1", "SLC6A2", "TPH2", "FKBP5", "HTR2A", "CRHR1", "BDNF", "TSPAN5", "ERICH3", "NR3C1", 
                      "GRM7", "PCLO", "NEGR1", "PAX6", "LHX2", "LHX9", "SIM1", "TFAP2A", "OTX2", "NTS", "SIX3", "RFX1", "RFX3", "SF1",
                      "DBX1", "FOXG1", "FEZ1", "HES1", "ISL1", "GATA2", "GATA3", "BRN3")

genesToAdd <- interestingGenes[which(!interestingGenes %in% genesKept)]

# add any genes back into the list which were removed during filtering
genesKept <- c(genesKept, genesToAdd)

# filter expression matrix to containing only highly expressed and interesting genes
exprMat_filtered <- exprMat[rownames(exprMat) %in% genesKept, ]

# calculated correlation between expression matrix and TF database
runCorrelation(exprMat_filtered, scenicOptions)

# save scenicOptions to use later
saveRDS(scenicOptions, file = "/Users/michellelam/SCENIC/scenicOptions.rds")

# normalise expression matrix
exprMat_filtered <- log2(exprMat_filtered+1)

# save expression matrix 
saveRDS(exprMat_filtered, "/Users/michellelam/SCENIC/exprMat_filtered.rds")

# load in the motifannotation
data(list="motifAnnotations_hgnc_v9", package="RcisTarget")

# rename the motif annnotion by attributing it to the correct name
motifAnnotations_hgnc <- motifAnnotations_hgnc_v9

# register the number of cores
nCores <- 8
cl <- makeCluster(nCores)
registerDoParallel(cl)

# run GENIE3
runGenie3(exprMat_filtered, scenicOptions)

# build the gene regulatory network and identify cell states ----------------------------------------

# set parameters for running
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 123

# set up parallel backend
registerDoParallel(cores = 10)

# import correct expression matrix into scenicOptions
scenicOptions@inputDatasetInfo$int_01 <- exprMat_filtered

# get co-expression modules
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)

# get regulons using Rcistarget (TF motif analysis)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod c("top5perTarget"))

# score GRN in the cells with AUCell
exprMat_log <- log2(exprMat+1) # normalised unfiltered expression matrix
scenicOptions <- runSCENIC_3_scorecells(scenicOptions, exprMat_log)

# save the analysis
scenicOptions@settings$dbs

# binarise the network activity (regulon on/off) and cluster on regulon activity ----------------------------------------

# import new scenicOptions again
scenicOptions <- readRDS("/Users/michellelam/int/scenicOptions.rds")

# binarise the AUC
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)

# select PCs
nPcs <- c(30, 40, 50, 60, 70)

# run t-SNE with different settings
fileNames_AUC <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(20, 30, 40, 50, 60))
fileNames_AUC <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(20, 30, 40, 50, 60), onlyHighConf=TRUE, filePrefix="int/tSNE_oHC")

# plot all tSNE plots ot visualise them
png("tSNE_plots.png", width = 2400, height = 2400) 

# par(mar = c(4, 4, 2, 1))
par(mfrow=c(length(nPcs), 5))
fileNames_AUC <- paste0("int/",grep(".Rds", grep("tSNE_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames_AUC, scenicOptions, showLegend=FALSE, varName="CellType", cex=.5)
dev.off() # 50 PCs and prepl 40 looks best

# save default tsne parameters
scenicOptions@settings$defaultTsne$aucType <- "AUC"
scenicOptions@settings$defaultTsne$dims <- 50
scenicOptions@settings$defaultTsne$perpl <- 40
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

# check if right tsne map is being used
print(tsneFileName(scenicOptions))

# exploring and interpreting results ----------------------------------------

# improt tsne map 
tSNE_scenic <- readRDS(tsneFileName(scenicOptions))

# load aucell regulon results
aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")

# get regulon scores
regulon_matrix <- getAUC(aucell_regulonAUC)

# get tSNE data
tSNE_data <- data.frame(tSNE_scenic$Y)

# get celltype, marker, sex and development info
cellInfo$barcode <- rownames(cellInfo)
markerInfo$barcode <- rownames(markerInfo)
sexInfo$barcode <- rownames(sexInfo)
stageInfo$barcode <- rownames(stageInfo)
tSNE_data$barcode <- rownames(tSNE_data)

# merge dataframes together
merging <- merge(cellInfo, markerInfo, by = "barcode")
merging <- merge(merging, sexInfo, by = "barcode")
merging <- merge(merging, stageInfo, by = "barcode")
tSNE_data <- merge(tSNE_data, merging, by = "barcode")

# plot tsne plot by cell type
ggplot(tSNE_data, aes(x = tSNE_scenic$Y[, 1], y = tSNE_scenic$Y[, 2], color = CellType)) +
  geom_point(size = 1) +
  scale_color_manual(values = c("ARC" = "#A6CEE3",
                                "AH" = "#b3b3b3",
                                "DMH" = "#a6d854",
                                "LH" = "#bf812d",
                                "MMN" = "#fdbf6f",
                                "TE" = "#fb9a99",
                                "TH" = "#cab2d6",
                                "TT" = "#fc8d62",
                                "PH" = "#ffd92f",
                                "Prethalamus" = "#5aae61",
                                "PVN/SON" = "#e5c494",
                                "SCN" = "#b2df8a",
                                "SMN" = "#e46c6c",
                                "SMN/MMN" = "#e78ac3",
                                "Thalamus" = "#66c2a5",
                                "TMN" = "#8da0cb")) +
  theme_minimal() +
  labs(color = NULL) +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12, family = "Arial", face = "bold"))

# plot heatmap of regulon activity by hypothalamic nuclei

# plot DLX1 activity
tSNE_data$DLX1 <- regulon_matrix["DLX1_extended (25g)", ]

ggplot(tSNE_data, aes(x = tSNE_scenic$Y[, 1], y = tSNE_scenic$Y[, 2], color = DLX1)) +
  geom_point() +
  scale_color_gradient(low = "#377eb8", high = "#e41a1c") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12, family = "Arial", face = "bold")) +
  labs(title = "DLX")

# plot ARX activity
tSNE_data$ARX <- regulon_matrix["ARX_extended (62g)", ]

ggplot(tSNE_data, aes(x = tSNE_scenic$Y[, 1], y = tSNE_scenic$Y[, 2], color = ARX)) +
  geom_point() +
  scale_color_gradient(low = "#377eb8", high = "#e41a1c") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12, family = "Arial", face = "bold")) +
  labs(title = "ARX")

# plot DLX2 activity
tSNE_data$DLX2 <- regulon_matrix["DLX2_extended (21g)", ]

ggplot(tSNE_data, aes(x = tSNE_scenic$Y[, 1], y = tSNE_scenic$Y[, 2], color = DLX2)) +
  geom_point() +
  scale_color_gradient(low = "#377eb8", high = "#e41a1c") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12, family = "Arial", face = "bold")) +
  labs(title = "DLX2")

# plot DLX5 activity
tSNE_data$DLX5 <- regulon_matrix["DLX5_extended (28g)", ]

ggplot(tSNE_data, aes(x = tSNE_scenic$Y[, 1], y = tSNE_scenic$Y[, 2], color = DLX5)) +
  geom_point() +
  scale_color_gradient(low = "#377eb8", high = "#e41a1c") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12, family = "Arial", face = "bold")) +
  labs(title = "DLX5")

# plot LHX1 activity
tSNE_data$LHX1 <- regulon_matrix["LHX1_extended (24g)", ]

ggplot(tSNE_data, aes(x = tSNE_scenic$Y[, 1], y = tSNE_scenic$Y[, 2], color = LHX1)) +
  geom_point() +
  scale_color_gradient(low = "#377eb8", high = "#e41a1c") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12, family = "Arial", face = "bold")) +
  labs(title = "LHX1")

# plot ONECUT2 activity
tSNE_data$ONECUT2 <- regulon_matrix["ONECUT2 (19g)", ]

ggplot(tSNE_data, aes(x = tSNE_scenic$Y[, 1], y = tSNE_scenic$Y[, 2], color = ONECUT2)) +
  geom_point() +
  scale_color_gradient(low = "#377eb8", high = "#e41a1c") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12, family = "Arial", face = "bold")) +
  labs(title = "ONECUT2")

# plot NHL1H1 activity
tSNE_data$NHL1H1 <- regulon_matrix["NHLH1 (24g)", ]

ggplot(tSNE_data, aes(x = tSNE_scenic$Y[, 1], y = tSNE_scenic$Y[, 2], color = NHL1H1)) +
  geom_point() +
  scale_color_gradient(low = "#377eb8", high = "#e41a1c") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12, family = "Arial", face = "bold")) +
  labs(title = "NHLH1")

# plot EMX2 activity
tSNE_data$EMX2 <- regulon_matrix["EMX2_extended (34g)", ]

ggplot(tSNE_data, aes(x = tSNE_scenic$Y[, 1], y = tSNE_scenic$Y[, 2], color = EMX2)) +
  geom_point() +
  scale_color_gradient(low = "#377eb8", high = "#e41a1c") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12, family = "Arial", face = "bold")) +
  labs(title = "EMX2")

# plot FOXB1 activity
tSNE_data$FOXB1 <- regulon_matrix["FOXB1_extended (10g)", ]

ggplot(tSNE_data, aes(x = tSNE_scenic$Y[, 1], y = tSNE_scenic$Y[, 2], color = FOXB1)) +
  geom_point() +
  scale_color_gradient(low = "#377eb8", high = "#e41a1c") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12, family = "Arial", face = "bold")) +
  labs(title = "FOXB2")

# get regulon scores
tSNE_data$BARHL1 <- regulon_matrix["BARHL1_extended (23g)", ]

ggplot(tSNE_data, aes(x = tSNE_scenic$Y[, 1], y = tSNE_scenic$Y[, 2], color = BARHL1)) +
  geom_point() +
  scale_color_gradient(low = "#377eb8", high = "#e41a1c") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12, family = "Arial", face = "bold")) +
  labs(title = "BARHL1")

# get regulon scores
tSNE_data$LMX1A <- regulon_matrix["LMX1A_extended (30g)", ]

ggplot(tSNE_data, aes(x = tSNE_scenic$Y[, 1], y = tSNE_scenic$Y[, 2], color = LMX1A)) +
  geom_point() +
  scale_color_gradient(low = "#377eb8", high = "#e41a1c") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12, family = "Arial", face = "bold")) +
  labs(title = "LMX1A")

# plot GATA3 activity
tSNE_data$GATA3 <- regulon_matrix["GATA3 (29g)", ]

ggplot(tSNE_data, aes(x = tSNE_scenic$Y[, 1], y = tSNE_scenic$Y[, 2], color = GATA3)) +
  geom_point() +
  scale_color_gradient(low = "#377eb8", high = "#e41a1c") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12, family = "Arial", face = "bold")) +
  labs(title = "GATA3")

# plot OTX2 activity
tSNE_data$OTX2 <- regulon_matrix["OTX2_extended (12g)", ]

ggplot(tSNE_data, aes(x = tSNE_scenic$Y[, 1], y = tSNE_scenic$Y[, 2], color = OTX2)) +
  geom_point() +
  scale_color_gradient(low = "#377eb8", high = "#e41a1c") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12, family = "Arial", face = "bold")) +
  labs(title = "OTX2")

# plot MEIS1 activity
tSNE_data$MEIS1 <- regulon_matrix["MEIS1_extended (14g)", ]

ggplot(tSNE_data, aes(x = tSNE_scenic$Y[, 1], y = tSNE_scenic$Y[, 2], color = MEIS1)) +
  geom_point() +
  scale_color_gradient(low = "#377eb8", high = "#e41a1c") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12, family = "Arial", face = "bold")) +
  labs(title = "MEIS1")

# plot ONECUT1 activity
tSNE_data$ONECUT1 <- regulon_matrix["ONECUT1_extended (10g)", ]

ggplot(tSNE_data, aes(x = tSNE_scenic$Y[, 1], y = tSNE_scenic$Y[, 2], color = ONECUT1)) +
  geom_point() +
  scale_color_gradient(low = "#377eb8", high = "#e41a1c") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12, family = "Arial", face = "bold")) +
  labs(title = "ONECUT1")

# plot CUX1 activity
tSNE_data$CUX1 <- regulon_matrix["CUX1_extended (24g)", ]

ggplot(tSNE_data, aes(x = tSNE_scenic$Y[, 1], y = tSNE_scenic$Y[, 2], color = CUX1)) +
  geom_point() +
  scale_color_gradient(low = "#377eb8", high = "#e41a1c") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12, family = "Arial", face = "bold")) +
  labs(title = "CUX1")

# plot ZMAT4 activity
tSNE_data$ZMAT4 <- regulon_matrix["ZMAT4 (20g)", ]

ggplot(tSNE_data, aes(x = tSNE_scenic$Y[, 1], y = tSNE_scenic$Y[, 2], color = ZMAT4)) +
  geom_point() +
  scale_color_gradient(low = "#377eb8", high = "#e41a1c") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12, family = "Arial", face = "bold")) +
  labs(title = "ZMAT4")

# plot LEF1 activity
tSNE_data$LEF1 <- regulon_matrix["LEF1_extended (89g)", ]

ggplot(tSNE_data, aes(x = tSNE_scenic$Y[, 1], y = tSNE_scenic$Y[, 2], color = LEF1)) +
  geom_point() +
  scale_color_gradient(low = "#377eb8", high = "#e41a1c") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12, family = "Arial", face = "bold")) +
  labs(title = "LEF1")

# geet NKX2-1 activity
tSNE_data$NKX2.1 <- regulon_matrix["NKX2-1_extended (26g)", ]

ggplot(tSNE_data, aes(x = tSNE_scenic$Y[, 1], y = tSNE_scenic$Y[, 2], color = NKX2.1)) +
  geom_point() +
  scale_color_gradient(low = "#377eb8", high = "#e41a1c") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12, family = "Arial", face = "bold")) +
  labs(title = "NKX2-1")

# get BCL11A activity
tSNE_data$BCL11A <- regulon_matrix["BCL11A_extended (17g)", ]

ggplot(tSNE_data, aes(x = tSNE_scenic$Y[, 1], y = tSNE_scenic$Y[, 2], color = BCL11A)) +
  geom_point() +
  scale_color_gradient(low = "#377eb8", high = "#e41a1c") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12, family = "Arial", face = "bold")) +
  labs(title = "BCL11A")

# get SOX11 activity
tSNE_data$SOX11 <- regulon_matrix["SOX11_extended (73g)", ]

ggplot(tSNE_data, aes(x = tSNE_scenic$Y[, 1], y = tSNE_scenic$Y[, 2], color = SOX11)) +
  geom_point() +
  scale_color_gradient(low = "#377eb8", high = "#e41a1c") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12, family = "Arial", face = "bold")) +
  labs(title = "SOX11")

# get JUND activity
tSNE_data$JUND <- regulon_matrix["JUND_extended (275g)", ]

ggplot(tSNE_data, aes(x = tSNE_scenic$Y[, 1], y = tSNE_scenic$Y[, 2], color = JUND)) +
  geom_point() +
  scale_color_gradient(low = "#377eb8", high = "#e41a1c") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12, family = "Arial", face = "bold")) +
  labs(title = "JUND")

# plot ZEB1 activity
tSNE_data$ZEB1 <- regulon_matrix["ZEB1_extended (74g)", ]

ggplot(tSNE_data, aes(x = tSNE_scenic$Y[, 1], y = tSNE_scenic$Y[, 2], color = ZEB1)) +
  geom_point() +
  scale_color_gradient(low = "#377eb8", high = "#e41a1c") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12, family = "Arial", face = "bold")) +
  labs(title = "ZEB1")

# plot SIX6 activity
tSNE_data$SIX6 <- regulon_matrix["SIX6 (21g)", ]
ggplot(tSNE_data, aes(x = tSNE_scenic$Y[, 1], y = tSNE_scenic$Y[, 2], color = SIX6)) +
  geom_point() +
  scale_color_gradient(low = "#377eb8", high = "#e41a1c") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12, family = "Arial", face = "bold")) +
  labs(title = "SIX6")

# plot ASCL1 activity
tSNE_data$ASCL1 <- regulon_matrix["ASCL1_extended (12g)", ]

ggplot(tSNE_data, aes(x = tSNE_scenic$Y[, 1], y = tSNE_scenic$Y[, 2], color = ASCL1)) +
  geom_point() +
  scale_color_gradient(low = "#377eb8", high = "#e41a1c") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12, family = "Arial", face = "bold")) +
  labs(title = "ASCL1")


# plot LHX6 activity
tSNE_data$LHX6 <- regulon_matrix["LHX6_extended (14g)", ]

ggplot(tSNE_data, aes(x = tSNE_scenic$Y[, 1], y = tSNE_scenic$Y[, 2], color = LHX6)) +
  geom_point() +
  scale_color_gradient(low = "#377eb8", high = "#e41a1c") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12, family = "Arial", face = "bold")) +
  labs(title = "LHX6")

# plot SOX2 activity
tSNE_data$SOX2 <- regulon_matrix["SOX2_extended (14g)", ]

ggplot(tSNE_data, aes(x = tSNE_scenic$Y[, 1], y = tSNE_scenic$Y[, 2], color = SOX2)) +
  geom_point() +
  scale_color_gradient(low = "#377eb8", high = "#e41a1c") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12, family = "Arial", face = "bold")) +
  labs(title = "SOX2")

# plot JUN activity
tSNE_data$JUN <- regulon_matrix["JUN_extended (41g)", ]

ggplot(tSNE_data, aes(x = tSNE_scenic$Y[, 1], y = tSNE_scenic$Y[, 2], color = JUN)) +
  geom_point() +
  scale_color_gradient(low = "#377eb8", high = "#e41a1c") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12, family = "Arial", face = "bold")) +
  labs(title = "JUN")

# plot FOS activity
tSNE_data$FOS <- regulon_matrix["FOS_extended (37g)", ]

ggplot(tSNE_data, aes(x = tSNE_scenic$Y[, 1], y = tSNE_scenic$Y[, 2], color = FOS)) +
  geom_point() +
  scale_color_gradient(low = "#377eb8", high = "#e41a1c") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12, family = "Arial", face = "bold")) +
  labs(title = "FOS")

ggsave("FOS_regulon_tsne.png", height = 8, width = 8)

# gplot HDAC2 activity
tSNE_data$HDAC2 <- regulon_matrix["HDAC2_extended (42g)", ]

ggplot(tSNE_data, aes(x = tSNE_scenic$Y[, 1], y = tSNE_scenic$Y[, 2], color = HDAC2)) +
  geom_point() +
  scale_color_gradient(low = "#377eb8", high = "#e41a1c") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12, family = "Arial", face = "bold")) +
  labs(title = "HDAC2")

# plot GBX2 activity
tSNE_data$GBX2 <- regulon_matrix["GBX2_extended (28g)", ]

ggplot(tSNE_data, aes(x = tSNE_scenic$Y[, 1], y = tSNE_scenic$Y[, 2], color = GBX2)) +
  geom_point() +
  scale_color_gradient(low = "#377eb8", high = "#e41a1c") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12, family = "Arial", face = "bold")) +
  labs(title = "GBX2")

# plot LHX9 activity
tSNE_data$LHX9 <- regulon_matrix["LHX9_extended (115g)", ]

ggplot(tSNE_data, aes(x = tSNE_scenic$Y[, 1], y = tSNE_scenic$Y[, 2], color = LHX9)) +
  geom_point() +
  scale_color_gradient(low = "#377eb8", high = "#e41a1c") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12, family = "Arial", face = "bold")) +
  labs(title = "LHX9")

# ggsave("LHX9_regulon_tsne.png", height = 8, width = 8)

# get regulon scores
tSNE_data$FOXP2 <- regulon_matrix["FOXP2_extended (77g)", ]

ggplot(tSNE_data, aes(x = tSNE_scenic$Y[, 1], y = tSNE_scenic$Y[, 2], color = FOXP2)) +
  geom_point() +
  scale_color_gradient(low = "#377eb8", high = "#e41a1c") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12, family = "Arial", face = "bold")) +
  labs(title = "FOXP2")

# aggreagate regulon activity by cell type
aggregated_regulon <- tSNE_data %>%
  group_by(Marker) %>%
  summarise(across(-c(barcode, tsne1, tsne2, CellType, Sex, Stage), \(x) mean(as.numeric(x), na.rm = TRUE))) %>%
  column_to_rownames(var = "Marker") %>%
  as.matrix()

# exploring sex differences ----------------------------------------

# aggregate regulon activity by sex and stage
aggregated_regulon_sex <- tSNE_data %>%
  group_by(Sex, Stage) %>%
  summarise(across(-c(barcode, tsne1, tsne2, CellType, Marker), \(x) mean(as.numeric(x), na.rm = TRUE)))

# gene regulatory network visualisation
# load TF-gene file
regulons <- loadInt(scenicOptions, "regulons")

# convert list to a data frame
regulons <- bind_rows(lapply(names(regulons), function(tf) {
  tibble(TF = tf, Target = regulons[[tf]])
  }))

# pvn TFs
pvn_tfs <- c("ZMAT4", "LHX1", "ONECUT1", "ONECUT2", "SOX11", "FOXP2", "FOXP2_extended", "SOX11_extended",
             "ONECUT1_extended", "ONECUT2_extended", "ZMAT4_extended", "LHX1_extended")

# filter regulon data for PVN TFs
pvn_regulons <- regulons[regulons$TF %in% pvn_tfs, ]

# remove _extended suffix
pvn_regulons$TF <- gsub("_extended", "", pvn_regulons$TF)

# export PVN regulon df
write.csv(pvn_regulons, "PVN_regulons.csv", row.names = F, quote = F) # import to Cytoscape for visualisation



