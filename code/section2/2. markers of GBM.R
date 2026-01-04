library(Seurat)
library(foreach)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(scales) 
library(harmony)

setwd("/home/yinshumin/pro_cfmf/")


#######1. GBM CFMF########
seu         <- readRDS("./input/GBM/Data_Brain_Neftel2019_SmartSeq2_TPM.rds")
seu.ls      <- SplitObject(seu, split.by = "sample")
seu.ma      <- subset(seu, malignant == "yes")
seu.ma      <- seu.ma %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>%  RunUMAP(dims = 1:30)
seu         <- seu.ma %>% RunHarmony(group.by.vars = "sample") %>%  RunUMAP(reduction = "harmony", dims = 1:30)

pct.ls <- lapply(seu.ls, function(x) {
  cutoff    <- 0 
  gene.pct  <- data.frame(gene = rownames(x), pct = rowMeans(x[["RNA"]]@data > cutoff), 
                          count = rowSums(x[["RNA"]]@data > cutoff))  
  gene.pct  <- gene.pct %>% filter(count > 3, pct < 1)
})
names(pct.ls) <- names(seu.ls)


source("./code/section1/FDR.function.R")
purity.ls <- foreach(i = seu.ls, j = pct.ls) %do% {
  purity     <- calculate.purity(i, gene.pct = j, neighbor = 50, ncores = 40)
}
names(purity.ls) <- names(seu.ls)


source("./code/section1/Silhouette.function.R")
silhouette.ls <- foreach(i = seu.ls, j = purity.ls) %do% {
  silhouette <- calculate.silhouette(i, purity = j, ncores = 40)
}
names(silhouette.ls) <- names(seu.ls)

save(list = c("seu","silhouette.ls","seu.ls"),
     file = "./output/section2/seu and silhouette of Data_Brain_Neftel2019_SmartSeq2_TPM.RData")








