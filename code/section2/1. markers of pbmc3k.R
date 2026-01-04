library(Seurat)
library(tidyverse)
library(foreach)
library(patchwork)
library(scales)

setwd("/home/yinshumin/pro_cfmf/")

seu       <- readRDS("./input/pbmc/pbmc3k.rds")
seu       <- NormalizeData(seu)
seu       <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seu)
seu       <- ScaleData(seu, features = all.genes)
seu       <- RunPCA(seu, features = VariableFeatures(seu))
seu       <- FindNeighbors(seu, dims = 1:30)
seu       <- FindClusters(seu, resolution = 1)
seu       <- RunUMAP(seu, dims = 1:30)


####1. FDR and silhouette of pbmc3k####
cutoff    <- 0 
gene.pct  <- data.frame(gene = rownames(seu), pct = rowMeans(seu[["RNA"]]@data > cutoff), 
                       count = rowSums(seu[["RNA"]]@data > cutoff))  
gene.pct  <- gene.pct %>% filter(count > 3, pct < 1)

source("./code/section1/FDR.function.R")
purity     <- calculate.purity(seu, gene.pct = gene.pct, neighbor = 50, ncores = 40)

source("./code/section1/Silhouette.function.R")
silhouette <- calculate.silhouette(seu, purity, ncores = 40)

save(list = c("silhouette", "seu"), file = "./output/section2/FDR_and_silhouette_of_pbmc3k.RData")



####2. DEGs#####
load("./output/section2/FDR_and_silhouette_of_pbmc3k.RData")

Idents(seu) <- "cell_type"
DEG         <- FindAllMarkers(seu)
save("DEG", file = "./output/section2/GEDs_of_celltype.RData")


