# loading required packages
library(Seurat)
library(tidyverse)
library(ggplot2)
library(SeuratWrappers)
library(monocle3)
library(dtw)

# convert to cell dataset object and transfer clustering information --------------------------------------------------------------------------------

# convert seurat object to celldataset object
cds <- as.cell_data_set(seurat)

# assign partitions
partitions <- c(rep(1, length(cds@colData@rownames)))
names(partitions) <- cds@colData@rownames
partitions <- as.factor(partitions)
cds@clusters$UMAP$partitions <- partitions

# assign cluster info from seurat clustering
list_cluster <- seurat$seurat_clusters
cds@clusters$UMAP$clusters <- list_cluster

# assign UMAP coordinates from seurat clustering
cds@int_colData@listData$reducedDims$UMAP <- seurat@reductions$umap@cell.embeddings

# plot UMAP without cell type annotations
plot_cells(cds,
           color_cells_by = "cluster",
           label_groups_by_cluster = F,
           group_label_size = 5) +
  theme(legend.position = "right")

# plot UMAP with cell type annotations
plot_cells(cds,
           color_cells_by = "cell_type",
           label_groups_by_cluster = F,
           group_label_size = 5) +
  theme(legend.position = "right")

# trajectory graph and psuedotime analysis --------------------------------------------------------------------------------

# subset cds object to include samples used for sex differences
cds_sex <- cds[, colData(cds)$sample == c("GW7_sample1", "GW8_sample1",
                                          "GW12_sample1", "GW12_sample2", "GW15_sample1")]

# subset cds_sex object to include neurons, NE and NP cells only to study neuronal trajectory
cds_sex <- cds_sex[, colData(cds_sex)$cell_type == c("Neuron", "NP", "NE")]

# recompute clusters
cds_sex <- cluster_cells(cds_sex)

# learn trajectory graph
cds_sex <- learn_graph(cds_sex, use_partition = F)

# visualise trajectory graph
plot_cells(cds_sex,
           color_cells_by = "cell_type",
           label_groups_by_cluster = F,
           label_branch_points = F,
           label_roots = F,
           label_leaves = F,
           group_label_size = 5,
           show_trajectory_graph = F) +
  theme(legend.position = "right")

# order cells by psuedotime
cds_sex <- order_cells(cds_sex, 
                       reduction_method = "UMAP", 
                       root_cells = rownames(colData(cds_sex))[colData(cds_sex)$cell_type == "NE"]) 

# visualise pseudotime trajectory
plot_cells(cds_sex,
           color_cells_by = "pseudotime",
           label_groups_by_cluster = F,
           label_branch_points = F,
           label_roots = F,
           label_leaves = F) +
  theme(legend.position = "right")
# root cells are ones with smallest psuedotime (earlier states)

# store pseudotime data into metadata of cds object
cds_sex$pseudotime <- pseudotime(cds_sex)

# finding genes which change as a function of pseudotime --------------------------------------------------------------------------------

# run DEG on pseudotime trajectory
pseudo_deg <- graph_test(cds_sex, neighbor_graph = "principal_graph", cores = 4)

# filter for top differentially expressed genes
pseudo_deg %>%
  arrange(q_value) %>%
  filter(status == "OK") %>%
  head(n = 50)

# run DEG on pseudotime trajectory with sex as a covariate
pseudo_sex_deg <- fit_models(cds_sex, model_formula_str = "~sex + pseudotime")

# extract results from model
pseudo_sex_deg_summary <- pseudo_sex_deg$model_summary

# initialize an empty list to store results
sex_results_list <- list()

# loop through the summary_model column and extract relevant information for each gene
for (gene in names(pseudo_sex_deg_summary)) {
  model <- pseudo_sex_deg_summary[[gene]]  # get the model for this gene
  
  # Extract the coefficient for sex and its p-value
  sex_coef <- model$coefficients["sexmale", "Estimate"]  # Coefficient for sex
  sex_pvalue <- model$coefficients["sexmale", "Pr(>|t|)"]  # P-value for sex
  
  # Store the results in a list (instead of a data frame to avoid row issues)
  sex_results_list[[gene]] <- c(gene_id = gene,
                                coefficient = sex_coef,
                                pvalue = sex_pvalue)
}

# Convert the list to a data frame
pseudo_deg_sex_results <- do.call(rbind, sex_results_list)
pseudo_deg_sex_results <- as.data.frame(pseudo_deg_sex_results, stringsAsFactors = FALSE)

# Convert the columns to appropriate types
pseudo_deg_sex_results$sex_coefficient <- as.numeric(pseudo_deg_sex_results$coefficient)
pseudo_deg_sex_results$sex_pvalue <- as.numeric(pseudo_deg_sex_results$pvalue)

# perform pvalue adjustment for multiple comparisons
pseudo_deg_sex_results$pval_adj <- p.adjust(pseudo_deg_sex_results$pvalue, method = "BH") 

# extract significant results
pseudo_sex_deg_sig <- pseudo_deg_sex_results[pseudo_deg_sex_results$pval_adj < 0.05,]

# male trajectory graph and psuedotime analysis ----------------------------------------

# subset cds object to include male samples used for sex differences
cds_male <- cds_sex[, colData(cds_sex)$sex == "male"]

# recompute clusters
cds_male <- cluster_cells(cds_male)

# learn trajectory graph
cds_male <- learn_graph(cds_male, use_partition = F)

# visualise trajectory graph
plot_cells(cds_male,
           color_cells_by = "cell_type",
           label_groups_by_cluster = F,
           label_branch_points = F,
           label_roots = F,
           label_leaves = F,
           group_label_size = 5) +
  theme(legend.position = "right")

# order cells by psuedotime
cds_male <- order_cells(cds_male, 
                        reduction_method = "UMAP", 
                        root_cells = rownames(colData(cds_male))[colData(cds_male)$cell_type == "NE"]) 

# visualise pseudotime trajectory
plot_cells(cds_male,
           color_cells_by = "pseudotime",
           label_groups_by_cluster = F,
           label_branch_points = F,
           label_roots = F,
           label_leaves = F) +
  theme(legend.position = "right")
# root cells are ones with smallest psuedotime (earlier states)

# store pseudotime data into metadata of cds object
cds_male$pseudotime <- pseudotime(cds_male)

# run DEG on pseudotime trajectory
pseudo_male_deg <- graph_test(cds_male, neighbor_graph = "principal_graph", cores = 4)

# filter for top differentially expressed genes
pseudo_male_deg %>%
  arrange(q_value) %>%
  filter(status == "OK") %>%
  head(n = 50)

# female trajectory graph and psuedotime analysis --------------------------------------------------------------------------------


# subset cds object to include female samples used for sex differences
cds_female <- cds_sex[, colData(cds_sex)$sex == "female"]

# recompute clusters
cds_female <- cluster_cells(cds_female)

# learn trajectory graph
cds_female <- learn_graph(cds_female, use_partition = F)

# visualise trajectory graph
plot_cells(cds_female,
           color_cells_by = "cell_type",
           label_groups_by_cluster = F,
           label_branch_points = F,
           label_roots = F,
           label_leaves = F,
           group_label_size = 5) +
  theme(legend.position = "right")

# order cells by psuedotime
cds_female <- order_cells(cds_female, 
                          reduction_method = "UMAP", 
                          root_cells = rownames(colData(cds_female))[colData(cds_female)$cell_type == "NE"]) 

# visualise pseudotime trajectory
plot_cells(cds_female,
           color_cells_by = "pseudotime",
           label_groups_by_cluster = F,
           label_branch_points = F,
           label_roots = F,
           label_leaves = F) +
  theme(legend.position = "right")
# root cells are ones with smallest psuedotime (earlier states)

# store pseudotime data into metadata of cds object
cds_female$pseudotime <- pseudotime(cds_female)

# run DEG on pseudotime trajectory
pseudo_female_deg <- graph_test(cds_female, neighbor_graph = "principal_graph", cores = 4)

# filter for top differentially expressed genes
pseudo_female_deg %>%
  arrange(q_value) %>%
  filter(status == "OK") %>%
  head(n = 50)

# dynamic time wwraping --------------------------------------------------------------------------------

# extract pseudotime trajectories for male and female
pseudotime_male <- pseudotime(cds_male)
pseudotime_female <- pseudotime(cds_female)

# interpolation of the male pseudotime to match the female pseudotime length
pseudotime_male <- approx(pseudotime_male, n = length(pseudotime_female))$y

# align pseudotime trajecotries between male and female samples
dtw <- dtw(pseudotime_male, pseudotime_female, keep.internals = T)

# output the DTW distance and visualise the warping path
dtw$distance
plot(dtw)
# wraping path shows vertical line and then a sharp bend later on
# this suggests that male graph is more compressed at the start so neural development occurs faster here
# there is then a sharp bend which suggests sex-divergence here

# plot density graph of cells across pseudotime for each sex to investigate differences
pseudotime_values <- pseudotime(cds_sex)
pseudotime_df <- data.frame(cell = colnames(cds_sex), 
                            pseudotime = pseudotime_values, 
                            sex = colData(cds_sex)$sex)

ggplot(pseudotime_df, aes(x = pseudotime, fill = sex)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density of Cells Across Pseudotime by Sex",
       x = "Pseudotime",
       y = "Density") +
  theme_classic() +
  scale_fill_manual(values = c("#e78ac3", "#8da0cb"))

# extract gene names from seurat metadata nad add to monocle3 objects
gene_metadata <- data.frame(
  gene_id = rownames(seurat),  # Gene IDs (rownames of Seurat object)
  gene_short_name = rownames(seurat)  # Add a column for gene names (or use custom names if available)
)
rowData(cds_sex) <- gene_metadata
rowData(cds_male) <- gene_metadata
rowData(cds_female) <- gene_metadata

# subset cds object on interesting genes
pseudotime_subset <- rownames(subset(as.data.frame(rowData(cds_sex)), gene_short_name %in% c("CRMP1", "CSDE1")))
cds_subset <- cds_male[pseudotime_subset,]

# plot line graph of pseudotime values for genes
plot_genes_in_pseudotime(cds_subset, color_cells_by = "pseudotime")
cds_subset <- cds_female[pseudotime_subset,]
plot_genes_in_pseudotime(cds_subset, color_cells_by = "pseudotime")

# subset hypothalamic neuron data by sex subsetted samples
development <- subset(hypothalamic, subset = sample %in% c("GW7_sample1", "GW8_sample1", "GW12_sample1", "GW12_sample2", "GW15_sample1"))

# set ident as sex
development <- SetIdent(development, value = development$sex)

# plot heatmap of interesting genes by sex
pseudo_genes <- c("STMN1", "MLLT11", "CFL1", "CD24", "GAP43", "UCHL1", "PPIA", "ZBTB20", "MAP1B", "MARCKSL1", "TXN", "TCF7L2", "FTH1", "HNRNPK", "MIF",
                  "EIF4G2", "YWHAE", "MARCKS", "PRDX2", "GDI1", "YWHAZ", "HSP90AB1", "DLG2", "NSG1", "GPC2", "RPL34", "EEF1A1", "PPDPF", "HDAC2", "RACK1",
                  "CRMP1", "CSDE1")

DotPlot(development, features = pseudo_genes, group.by = "sex") +
  scale_color_gradient(low = "lightgrey", high = "#8da0cb")

# Extract average expression for selected genes grouped by sex
pseudo_genes_av <- AverageExpression(development, features = pseudo_genes, group.by = "sex")$RNA

