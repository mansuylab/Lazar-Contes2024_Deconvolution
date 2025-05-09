# Deconvolution.R
# This script runs deconvolution of bulk RNA-seq data using CIBERSORTx, with pre-prepared reference datasets.
# It processes multiple reference datasets and uses them to estimate cell type proportions in bulk RNA-seq samples in Lazar-Contes et al 2024.

# Reference:
# Newman AM, Steen CB, Liu CL, Gentles AJ, Chaudhuri AA, Scherer F, Khodadoust MS, Esfahani MS, Luca BA, Steiner D, Diehn M, Alizadeh AA. 
# Determining cell type abundance and expression from bulk tissues with digital cytometry. 
# Nat Biotechnol. 2019 Jul;37(7):773-782. doi: 10.1038/s41587-019-0114-2. Epub 2019 May 6. PMID: 31061481; PMCID: PMC6610714.


# Clear the workspace
rm(list = ls())

# Set working directory (adjust as needed)
setwd("path")

# CIBERSORTx login details
username <- "username"
token <- "token"

# Input bulk RNA-seq dataset for deconvolution
mixture_name <- "LazarContes2024_Samples.txt"

# CIBERSORTx Singularity image file
sif_file <- "./cibersortx-fractions.sif"

# Directory containing reference datasets
reference_dir <- "path/for/references"

# Define the specific reference datasets to be used in this analysis
# Could be generated from PrepTan2020.R, PrepGreen2018.R and PrepHermann2018.R
reference_names <- c(
  "Tan2020_GermCellsReclustered_P7_Rep1.txt",
  "Tan2020_GermCellsReclustered_P7_Rep2.txt",
  "Tan2020_AllCells_P7_Rep1.txt",
  "Tan2020_AllCells_P7_Rep2.txt",
  "Hermann2018_P6.txt",
  "Hermann2018_Adu_ID4sorted.txt",
  paste0("Green2018_", c("ST1", "ST2", "ST3", "ST4", "ST5", "ST6", "ST7", "ST8", 
                         "SPG1", "SPG2", "SPG3", "SER1", "SER2", "SER3", "SER4", 
                         "SER5", "SER6", "SER7", "SER8"), ".txt")
)

# Tag to be appended to output folders
folder_tag <- "LazarContes2024_Deconvolution_"

# Generate folder names based on reference datasets
folder_names <- paste0(folder_tag, sub(".txt$", "", reference_names))

# Loop through each reference dataset and run the deconvolution
for (i in seq_along(folder_names)) {
  folder_name <- folder_names[i]
  reference_name <- reference_names[i]
  
  # Create the output directory if it doesn't exist
  output_dir <- paste0("./outputs/", folder_name)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Construct the CIBERSORTx command
  cmd <- sprintf(
    "apptainer exec --bind ./data:/src/data --bind ./outputs/%s:/src/outdir %s /src/CIBERSORTxFractions --username %s --token %s --single_cell TRUE --refsample %s --mixture %s", 
    folder_name, sif_file, username, token, reference_name, mixture_name
  )
  
  # Execute the CIBERSORTx command
  system(cmd)
}
