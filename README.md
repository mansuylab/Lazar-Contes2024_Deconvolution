This repository contains scripts used for the deconvolution analysis in Lazar-Contes et al., eLife 2024, Dynamics of transcriptional programs and chromatin accessibility in mouse spermatogonial cells from early postnatal to adult life (https://doi.org/10.7554/eLife.91528.3).

## Scripts

**Prep_.R scripts**  
These scripts prepare the reference datasets, by processing single-cell RNA-seq data, performing clustering (if cluster IDs are not provided by the original publication), and generating reference files for CIBERSORTx.

**Deconvolution.R**  
This script runs the deconvolution analysis (CIBERSORTx) using the bulk RNA-seq dataset in Lazar-Contes et al., 2024 and all of the prepared reference datasets.
