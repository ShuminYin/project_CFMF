
# CFMF

[![R-CMD-check](https://github.com/bioklab/CFMF/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ShuminYin/project_CFMF.git)
The goal of CFMF (Clustering-Free Cell Marker Finder) is to identify cell marker 
genes from single-cell RNA sequencing (scRNA-seq) data without relying on prior cell clustering. 
By leveraging local neighborhood enrichment and global compactness in the transcriptomic space, 
CFMF avoids the limitations of resolution parameter selection and effectively detects markers for rare cell populations and continuous biological states.



## Installation

You can install the development version of CFMF from [GitHub](https://github.com/ShuminYin/project_CFMF.git) with:

``` r
# install.packages("devtools")
install_github("ShuminYin/project_CFMF", subdir = "CFMF")



## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(CFMF)
library(CFMF)
library(Seurat)
library(dplyr)

# 1. Load your Seurat object (Must have PCA pre-calculated)
# seu <- Read10X(...) %>% CreateSeuratObject() %>% ... %>% RunPCA()

# 2. Run the CFMF algorithm
# neighbor = 50 is the recommended default
# ncores = 10 sets the number of threads for parallel processing
cfmf_results <- RunCFMF(seu, neighbor = 50, ncores = 10)

# 3. Filter for significant markers
# Recommended thresholds: -log10(FDR) > 200 and Silhouette > 0
final_markers <- cfmf_results %>%
  mutate(logFDR = -log10(FDR)) %>%
  filter(logFDR > 200 & silhouette > 0) %>%
  arrange(desc(silhouette))

# View results
head(final_markers)
