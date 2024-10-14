# Prep_Green2018.R
# This script processes single-cell RNA-seq data from Green et al., 2018, performing clustering and generating cell type reference files for CIBERSORTx.
# Filtering steps are based on the methods described in Green et al., 2018.

# Reference:
# Green CD, Ma Q, Manske GL, et al. A Comprehensive Roadmap of Murine Spermatogenesis Defined by Single-Cell RNA-Seq. 
# Developmental Cell. 2018;46(5):651-667.e10. doi: 10.1016/j.devcel.2018.07.025

# Load required libraries
library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)
library(data.table)
library(ggplot2)

# Clear workspace and set random seed for reproducibility
rm(list = ls())
set.seed(34)

# Define paths for input and output data
home <- "/Users/kuzel/NEL/deconvolution_scRNAseq"
path_to_outputs <- "/Users/kuzel/NEL/deconvolution_scRNAseq/output_files"

# Load the merged gene expression matrix for the 25 ST batches
file_path <- paste0(home, "GSE112393_MergedAdultMouseST25_DGE.txt")
mouse_data <- fread(file_path)

# Exclude the gene name column for downstream calculations
dgedata <- mouse_data[, .SD, .SDcols = -"V1"]

# Filter genes based on thresholds for UMI counts and cells
nUMIperGene <- rowSums(dgedata)
nCellperGene <- rowSums(dgedata > 0)

# Define thresholds for gene filtering
thres.nUMIperGene <- 20
thres.nCellperGene <- 15

# Convert the data table to a standard data frame and set gene names as rownames
dgedata <- as.data.frame(dgedata)
rownames(dgedata) <- mouse_data$V1

# Identify genes passing the thresholds
gene1 <- rownames(dgedata)[nUMIperGene > thres.nUMIperGene & nCellperGene > thres.nCellperGene]

# Add back 24 known genes of interest that do not meet the filtering thresholds
gene2 <- c("Ngf", "T", "Tbx20", "Tert", "Olig2", "Prl", "Hmx2", "Itgb6", "Mageb1", "Neurod1", 
           "Neurod2", "Ptchd4", "Prl6a1", "Fgf17", "Hpx", "Npy4r", "Fgf6", "Foxf2", "Prokr2", 
           "Sry", "Ccl3", "Prf1", "Pax8", "Dcx", "Ly6g5c", "Ereg", "Pgr", "Sycp2l", "Prnd", 
           "Cysltr1", "Pdf")
gene2 <- gene2[!(gene2 %in% gene1)]  # Exclude duplicates already in gene1

# Combine filtered genes with genes of interest
genekeep <- c(gene1, gene2)
dgedata_filtered <- dgedata[genekeep, ]

# Create a Seurat object using the filtered gene expression data
dge <- CreateSeuratObject(counts = dgedata_filtered)

# Normalize the data using log normalization
dge <- NormalizeData(dge, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
dge <- FindVariableFeatures(dge, selection.method = "vst", nfeatures = 2000)

# Scale the data and perform PCA for dimensionality reduction
dge <- ScaleData(dge)
dge <- RunPCA(dge, features = VariableFeatures(object = dge))

# Perform UMAP for visualization
dge <- RunUMAP(dge, dims = 1:10)

# Generate a UMAP plot showing cluster assignments
DimPlot(dge, reduction = "umap", group.by = "orig.ident")

# Generate feature plots for germ cell markers
features1 <- c("Dazl", "Ddx4", "Zbtb16", "Gfra1", "Kit", "Stra8", "Id4", "Uchl1")
fplot1 <- FeaturePlot(dge, features = features1, ncol = 4, order = TRUE)

# Perform clustering with a resolution of 0.11
dge <- FindNeighbors(dge, dims = 1:20)
dge <- FindClusters(dge, resolution = 0.11)

# Generate a UMAP plot with cluster labels
DimPlot(dge, reduction = "umap", group.by = "seurat_clusters", label = TRUE)

# Identify cell types based on cluster assignments
germ_cells <- WhichCells(dge, idents = 4)
sertoli_cells <- WhichCells(dge, idents = 8)
leydig_cells <- WhichCells(dge, idents = 9)
ptm_cells <- WhichCells(dge, idents = 10)
stroma_cells <- WhichCells(dge, idents = 6)

# Define a function to categorize cells by type based on cluster identity
categorize_cell <- function(cell_id) {
  if (cell_id %in% germ_cells) {
    return("Germ_cells")
  } else if (cell_id %in% sertoli_cells) {
    return("Sertoli_cells")
  } else if (cell_id %in% leydig_cells) {
    return("Leydig_cells")
  } else if (cell_id %in% ptm_cells) {
    return("Ptm_cells")
  } else if (cell_id %in% stroma_cells) {
    return("Stroma_cells")
  } else {
    return("Other")
  }
}

# List of library names used for the deconvolution analysis
library_names <- c("ST1", "ST2", "ST3", "ST4", "ST5", "ST6", "ST7", "ST8", "SPG1", "SPG2", "SPG3",
                   "SER1", "SER2", "SER3", "SER4", "SER5", "SER6", "SER7", "SER8")

# Loop through each library and generate a categorized file for each
for (temp_name in library_names) {
  # Select columns matching the current library name
  selected_columns <- grep(paste0("^", temp_name, "_"), names(mouse_data), value = TRUE)
  selected_data <- mouse_data[, selected_columns, with = FALSE]
  
  # Categorize cells based on their cluster identities
  new_colnames <- sapply(colnames(selected_data), categorize_cell)
  colnames(selected_data) <- new_colnames
  
  # Create a data frame with gene symbols and categorized cell types
  matrix_df <- data.frame(GeneSymbol = mouse_data$V1, selected_data, check.names = FALSE)
  
  # Write the categorized data to a text file for further use in CIBERSORTx
  output_filename <- paste(path_to_outputs, "/Green2018_", temp_name, ".txt", sep = "")
  write.table(matrix_df, file = output_filename, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  # Clean up the workspace, keeping only necessary objects
  rm(list = setdiff(ls(), c("mouse_data", "path_to_outputs", "categorize_cell", "library_names", "temp_name",
                            "germ_cells", "sertoli_cells", "leydig_cells", "ptm_cells", "stroma_cells")))
}
