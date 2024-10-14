# Prep_Tan2020_AllCells.R

# This script processes single-cell RNA-seq data from Tan et al., 2020, performing clustering and generating cell type reference files for CIBERSORTx.
# It combines all samples (P7, P2, E18), performs clustering, and exports categorized cells from P7_Rep1 and P7_Rep2 for further deconvolution analysis.

# Reference:
# Tan K, Song HW, Wilkinson MF. Single-cell RNAseq analysis of testicular germ and somatic cell development during the perinatal period. 
# Development. 2020 Feb 3;147(3):dev183251. doi: 10.1242/dev.183251. PMID: 31964773; PMCID: PMC7033731.

# Load necessary libraries
library(Seurat)
library(Matrix)
library(ggplot2)

# Clear workspace and set a random seed for reproducibility
rm(list = ls())
set.seed(34)

# Define paths for saving plots and output files
path_to_plots <- "/Users/kuzel/NEL/deconvolution_scRNAseq/plots"
path_to_outputs <- "/Users/kuzel/NEL/deconvolution_scRNAseq/output_files"

# Function to load, process, and create a Seurat object from matrix files
process_matrix <- function(barcodes_path, features_path, matrix_path) {
  barcodes <- readLines(barcodes_path)
  features <- read.delim(features_path, header = FALSE, row.names = 1)
  matrix <- readMM(file = matrix_path)
  colnames(matrix) <- barcodes
  rownames(matrix) <- features$V2
  matrix[unique(rownames(matrix)), ] # Retain only unique genes
}

# Load data for all samples
# P7 samples
matrix_P7_Rep1 <- process_matrix("/Users/kuzel/NEL/deconvolution_scRNAseq/raw_files/Tan2020/GSM4185522_P7_Rep1_barcodes.tsv", 
                                 "/Users/kuzel/NEL/deconvolution_scRNAseq/raw_files/Tan2020/GSM4185522_P7_Rep1_features.tsv", 
                                 "/Users/kuzel/NEL/deconvolution_scRNAseq/raw_files/Tan2020/GSM4185522_P7_Rep1_matrix.mtx")
matrix_P7_Rep2 <- process_matrix("/Users/kuzel/NEL/deconvolution_scRNAseq/raw_files/Tan2020/GSM4185523_P7_Rep2_barcodes.tsv", 
                                 "/Users/kuzel/NEL/deconvolution_scRNAseq/raw_files/Tan2020/GSM4185523_P7_Rep2_features.tsv", 
                                 "/Users/kuzel/NEL/deconvolution_scRNAseq/raw_files/Tan2020/GSM4185523_P7_Rep2_matrix.mtx")

# P2 samples
matrix_P2_Rep1 <- process_matrix("/Users/kuzel/NEL/deconvolution_scRNAseq/raw_files/Tan2020/GSM3744443_P2_Rep1_barcodes.tsv", 
                                 "/Users/kuzel/NEL/deconvolution_scRNAseq/raw_files/Tan2020/GSM3744443_P2_Rep1_features.tsv", 
                                 "/Users/kuzel/NEL/deconvolution_scRNAseq/raw_files/Tan2020/GSM3744443_P2_Rep1_matrix.mtx")
matrix_P2_Rep2 <- process_matrix("/Users/kuzel/NEL/deconvolution_scRNAseq/raw_files/Tan2020/GSM3744444_P2_Rep2_barcodes.tsv", 
                                 "/Users/kuzel/NEL/deconvolution_scRNAseq/raw_files/Tan2020/GSM3744444_P2_Rep2_features.tsv", 
                                 "/Users/kuzel/NEL/deconvolution_scRNAseq/raw_files/Tan2020/GSM3744444_P2_Rep2_matrix.mtx")

# E18 samples
matrix_E18_Rep1 <- process_matrix("/Users/kuzel/NEL/deconvolution_scRNAseq/raw_files/Tan2020/GSM3744441_E18_Rep1_barcodes.tsv", 
                                  "/Users/kuzel/NEL/deconvolution_scRNAseq/raw_files/Tan2020/GSM3744441_E18_Rep1_features.tsv", 
                                  "/Users/kuzel/NEL/deconvolution_scRNAseq/raw_files/Tan2020/GSM3744441_E18_Rep1_matrix.mtx")
matrix_E18_Rep2 <- process_matrix("/Users/kuzel/NEL/deconvolution_scRNAseq/raw_files/Tan2020/GSM3744442_E18_Rep2_barcodes.tsv", 
                                  "/Users/kuzel/NEL/deconvolution_scRNAseq/raw_files/Tan2020/GSM3744442_E18_Rep2_features.tsv", 
                                  "/Users/kuzel/NEL/deconvolution_scRNAseq/raw_files/Tan2020/GSM3744442_E18_Rep2_matrix.mtx")

# Create Seurat objects for each sample
seurat_P7_Rep1 <- CreateSeuratObject(counts = matrix_P7_Rep1)
seurat_P7_Rep2 <- CreateSeuratObject(counts = matrix_P7_Rep2)
seurat_P2_Rep1 <- CreateSeuratObject(counts = matrix_P2_Rep1)
seurat_P2_Rep2 <- CreateSeuratObject(counts = matrix_P2_Rep2)
seurat_E18_Rep1 <- CreateSeuratObject(counts = matrix_E18_Rep1)
seurat_E18_Rep2 <- CreateSeuratObject(counts = matrix_E18_Rep2)

# Assign sample IDs to the Seurat objects for tracking
seurat_P7_Rep1$orig.ident <- "P7-1"
seurat_P7_Rep2$orig.ident <- "P7-2"
seurat_P2_Rep1$orig.ident <- "P2-1"
seurat_P2_Rep2$orig.ident <- "P2-2"
seurat_E18_Rep1$orig.ident <- "E18-1"
seurat_E18_Rep2$orig.ident <- "E18-2"

# Merge all samples into a combined Seurat object
P7 <- merge(seurat_P7_Rep1, y = seurat_P7_Rep2)
P2 <- merge(seurat_P2_Rep1, y = seurat_P2_Rep2)
E18 <- merge(seurat_E18_Rep1, y = seurat_E18_Rep2)

combined_seurat_obj <- merge(P7, y = c(P2, E18))

# Clean up the workspace
rm(list = setdiff(ls(), c("combined_seurat_obj", "path_to_plots", "path_to_outputs")))

# Calculate the percentage of mitochondrial gene expression (routine check, expecting none)
mitochondrial_genes <- grep("^MT-", rownames(combined_seurat_obj), value = TRUE)
combined_seurat_obj[["percent.mt"]] <- PercentageFeatureSet(combined_seurat_obj, features = mitochondrial_genes)

# Filter cells based on RNA features: remove cells with <200 or >6000 genes and cells with high mitochondrial gene expression
combined_seurat_obj <- subset(combined_seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 0.2)

# Normalize data and identify highly variable features for downstream analysis
combined_seurat_obj <- NormalizeData(combined_seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
combined_seurat_obj <- FindVariableFeatures(combined_seurat_obj, selection.method = "vst", nfeatures = 2000)

# Perform PCA for dimensionality reduction and then UMAP for visualization
combined_seurat_obj <- ScaleData(combined_seurat_obj)
combined_seurat_obj <- RunPCA(combined_seurat_obj, features = VariableFeatures(object = combined_seurat_obj))
combined_seurat_obj <- RunUMAP(combined_seurat_obj, dims = 1:10)

# Find clusters of cells using a low-resolution clustering method
combined_seurat_obj <- FindNeighbors(combined_seurat_obj, dims = 1:10)
combined_seurat_obj <- FindClusters(combined_seurat_obj, resolution = 0.025)

# Example feature plot showing germ cells to recapitulate cell types and clusters in Tan et al., 2020
FeaturePlot(combined_seurat_obj, features = c("Dazl", "Gfra1", "Kit"), ncol = 3)

# UMAP plot of clusters
DimPlot(combined_seurat_obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE)

# Identify cells in each cluster corresponding to major cell types
germ_cells <- WhichCells(combined_seurat_obj, idents = 2)
sertoli_cells <- WhichCells(combined_seurat_obj, idents = 0)
leydig_cells <- WhichCells(combined_seurat_obj, idents = 3)
ptm_cells <- WhichCells(combined_seurat_obj, idents = 6)
stroma_cells <- WhichCells(combined_seurat_obj, idents = 1)

# Reload raw P7_Rep1 and P7_Rep2 matrices to categorize cells
barcodes_P7_Rep1 <- readLines("/Users/kuzel/NEL/deconvolution_scRNAseq/raw_files/Tan2020/GSM4185522_P7_Rep1_barcodes.tsv")
features_P7_Rep1 <- read.delim("/Users/kuzel/NEL/deconvolution_scRNAseq/raw_files/Tan2020/GSM4185522_P7_Rep1_features.tsv", header = FALSE, row.names = 1)
matrix_P7_Rep1 <- readMM(file = "/Users/kuzel/NEL/deconvolution_scRNAseq/raw_files/Tan2020/GSM4185522_P7_Rep1_matrix.mtx")
colnames(matrix_P7_Rep1) <- barcodes_P7_Rep1
rownames(matrix_P7_Rep1) <- features_P7_Rep1$V2  

barcodes_P7_Rep2 <- readLines("/Users/kuzel/NEL/deconvolution_scRNAseq/raw_files/Tan2020/GSM4185523_P7_Rep2_barcodes.tsv")
features_P7_Rep2 <- read.delim("/Users/kuzel/NEL/deconvolution_scRNAseq/raw_files/Tan2020/GSM4185523_P7_Rep2_features.tsv", header = FALSE, row.names = 1)
matrix_P7_Rep2 <- readMM(file = "/Users/kuzel/NEL/deconvolution_scRNAseq/raw_files/Tan2020/GSM4185523_P7_Rep2_matrix.mtx")
colnames(matrix_P7_Rep2) <- barcodes_P7_Rep2
rownames(matrix_P7_Rep2) <- features_P7_Rep2$V2  

# Update column names of matrices with sample tags
colnames(matrix_P7_Rep1) <- paste0(barcodes_P7_Rep1, "_1_1")
colnames(matrix_P7_Rep2) <- paste0(barcodes_P7_Rep2, "_2_1")

# Categorize each cell according to clusters
categorize_cell <- function(cell_id) {
  if (cell_id %in% germ_cells) {
    return("Germ_cells")
  } else if (cell_id %in% sertoli_cells) {
    return("Sertoli_cells")
  } else if (cell_id %in% leydig_cells) {
    return("Leydig_cells")
  } else if (cell_id %in% ptm_cells) {
    return("PTM_cells")
  } else if (cell_id %in% stroma_cells) {
    return("Stroma_cells")
  } else {
    return("Other")
  }
}

# Categorize cells for P7_Rep1
new_colnames_1 <- vector("character", length = length(colnames(matrix_P7_Rep1)))
for (i in seq_along(new_colnames_1)) {
  new_colnames_1[i] <- categorize_cell(colnames(matrix_P7_Rep1)[i])
}

# Rename the columns of the matrix
colnames(matrix_P7_Rep1) <- new_colnames_1

# Get gene names
gene_names_1<-rownames((matrix_P7_Rep1))

# Create a data frame that includes GeneSymbol as the first column
matrix_df_Rep1 <- data.frame(GeneSymbol = gene_names_1, matrix_P7_Rep1, check.names = FALSE)

# Write P7_Rep1 data frame to a text file for further use with CIBERSORTx
write.table(matrix_df_Rep1, file = paste(path_to_outputs,"/","Tan2020_AllCells_P7_Rep1",".txt",sep=""), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


# Categorize cells for P7_Rep2
new_colnames_2 <- vector("character", length = length(colnames(matrix_P7_Rep2)))
for (i in seq_along(new_colnames_2)) {
  new_colnames_2[i] <- categorize_cell(colnames(matrix_P7_Rep2)[i])
}

# Rename the columns of the matrix for P7_Rep2
colnames(matrix_P7_Rep2) <- new_colnames_2

# Get gene names
gene_names_2 <- rownames((matrix_P7_Rep2))

# Create a data frame for P7_Rep2 including GeneSymbol as the first column
matrix_df_Rep2 <- data.frame(GeneSymbol = gene_names_2, matrix_P7_Rep2, check.names = FALSE)

# Write P7_Rep2 data frame to a text file for further use with CIBERSORTx
write.table(matrix_df_Rep2, file = paste(path_to_outputs,"/","Tan2020_AllCells_P7_Rep2",".txt",sep=""), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Below is the code for reclustering of germ cells (Figure S2E and F)
######################################################################

#clear the workspace
rm(list = setdiff(ls(), c("combined_seurat_obj","germ_cells","path_to_plots","path_to_outputs")))

# Ensure that 'cells_in_cluster_2' is a character vector of cell barcodes
germ_cells <- as.character(germ_cells)

# Subset 'combined_seurat_obj' to only include cells from 'germ_cells'
subset_seurat_obj <- subset(combined_seurat_obj, cells = germ_cells)

# Calculate mitochondrial gene percentage
mitochondrial_genes <- grep("^MT-", rownames(subset_seurat_obj), value = TRUE)
subset_seurat_obj[["percent.mt"]] <- PercentageFeatureSet(subset_seurat_obj, features = mitochondrial_genes)

# Basic filtering
subset_seurat_obj <- subset(subset_seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 0.2)

# Normalize Data and Identify Highly Variable Genes
subset_seurat_obj <- NormalizeData(subset_seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
subset_seurat_obj <- FindVariableFeatures(subset_seurat_obj, selection.method = "vst", nfeatures = 2000)


# Perform PCA for dimensionality reduction and then UMAP for visualization
subset_seurat_obj <- ScaleData(subset_seurat_obj)
subset_seurat_obj <- RunPCA(subset_seurat_obj, features = VariableFeatures(object = subset_seurat_obj))
subset_seurat_obj <- RunUMAP(subset_seurat_obj, dims = 1:10)

# UMAP plot of samples
DimPlot(subset_seurat_obj, reduction = "umap",group.by = "orig.ident")

# Feature plots for cell types of interest
# SSCs
FeaturePlot(subset_seurat_obj, features = c("Id4","Gfra1"),ncol=2,order = TRUE)
# DiffSPGs
FeaturePlot(subset_seurat_obj, features = c("Kit","Dmrtb1"),ncol=2,order = TRUE)


# Find clusters of cells 
subset_seurat_obj <- FindNeighbors(subset_seurat_obj, dims = 1:15)
subset_seurat_obj <- FindClusters(subset_seurat_obj, resolution = 0.1)

# UMAP plot of clusters
DimPlot(subset_seurat_obj, reduction = "umap", group.by = "seurat_clusters",label = TRUE)

# finding cell IDs for each cluster
GermReclustered_SSCs <- WhichCells(subset_seurat_obj, idents = 0)
GermReclustered_DiffSPGs <- WhichCells(subset_seurat_obj, idents = 2)

# Reload raw P7_Rep1 and P7_Rep2 matrices to categorize cells
barcodes_P7_Rep1 <- readLines("/Users/kuzel/NEL/deconvolution_scRNAseq/raw_files/Tan2020/GSM4185522_P7_Rep1_barcodes.tsv")
features_P7_Rep1 <- read.delim("/Users/kuzel/NEL/deconvolution_scRNAseq/raw_files/Tan2020/GSM4185522_P7_Rep1_features.tsv", header = FALSE, row.names = 1)
matrix_P7_Rep1 <- readMM(file = "/Users/kuzel/NEL/deconvolution_scRNAseq/raw_files/Tan2020/GSM4185522_P7_Rep1_matrix.mtx")
colnames(matrix_P7_Rep1) <- barcodes_P7_Rep1
rownames(matrix_P7_Rep1) <- features_P7_Rep1$V2  

barcodes_P7_Rep2 <- readLines("/Users/kuzel/NEL/deconvolution_scRNAseq/raw_files/Tan2020/GSM4185523_P7_Rep2_barcodes.tsv")
features_P7_Rep2 <- read.delim("/Users/kuzel/NEL/deconvolution_scRNAseq/raw_files/Tan2020/GSM4185523_P7_Rep2_features.tsv", header = FALSE, row.names = 1)
matrix_P7_Rep2 <- readMM(file = "/Users/kuzel/NEL/deconvolution_scRNAseq/raw_files/Tan2020/GSM4185523_P7_Rep2_matrix.mtx")
colnames(matrix_P7_Rep2) <- barcodes_P7_Rep2
rownames(matrix_P7_Rep2) <- features_P7_Rep2$V2  

# Update column names of matrices with sample tags
colnames(matrix_P7_Rep1) <- paste0(barcodes_P7_Rep1, "_1_1")
colnames(matrix_P7_Rep2) <- paste0(barcodes_P7_Rep2, "_2_1")


# Subset matrix_P7_Rep1
matrix_P7_Rep1_subset <- matrix_P7_Rep1[, colnames(matrix_P7_Rep1) %in% germ_cells]
matrix_P7_Rep1<-matrix_P7_Rep1_subset

# Subset matrix_P7_Rep2
matrix_P7_Rep2_subset <- matrix_P7_Rep2[, colnames(matrix_P7_Rep2) %in% germ_cells]
matrix_P7_Rep2<-matrix_P7_Rep2_subset

# Define a function to categorize each cell for new_clust
categorize_cell <- function(cell_id) {
  if (cell_id %in% GermReclustered_SSCs) {
    return("SSCs")
  } else if (cell_id %in% GermReclustered_DiffSPGs) {
    return("DiffSPGs")
  } else {
    return("Other")
  }
}


# Categorize cells for P7_Rep1
new_colnames_1 <- vector("character", length = length(colnames(matrix_P7_Rep1)))
for (i in seq_along(new_colnames_1)) {
  new_colnames_1[i] <- categorize_cell(colnames(matrix_P7_Rep1)[i])
}

# Rename the columns of the matrix
colnames(matrix_P7_Rep1) <- new_colnames_1
# Get gene names
gene_names_1<-rownames((matrix_P7_Rep1))

# Create a data frame that includes GeneSymbol as the first column
matrix_df_Rep1 <- data.frame(GeneSymbol = gene_names_1, matrix_P7_Rep1, check.names = FALSE)

# Write P7_Rep1 data frame to a text file
write.table(matrix_df_Rep1, file = paste(path_to_outputs,"/","Tan2020_GermCellsReclustered_P7_Rep1",".txt",sep=""), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


# Categorize cells for P7_Rep2
new_colnames_2 <- vector("character", length = length(colnames(matrix_P7_Rep2)))
for (i in seq_along(new_colnames_2)) {
  new_colnames_2[i] <- categorize_cell(colnames(matrix_P7_Rep2)[i])
}

# Rename the columns of the matrix for P7_Rep2
colnames(matrix_P7_Rep2) <- new_colnames_2
# Get gene names
gene_names_2 <- rownames((matrix_P7_Rep2))

# Create a data frame for P7_Rep2 including GeneSymbol as the first column
matrix_df_Rep2 <- data.frame(GeneSymbol = gene_names_2, matrix_P7_Rep2, check.names = FALSE)

# Write P7_Rep2 data frame to a text file
write.table(matrix_df_Rep2, file = paste(path_to_outputs,"/","Tan2020_GermCellsReclustered_P7_Rep2",".txt",sep=""), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

