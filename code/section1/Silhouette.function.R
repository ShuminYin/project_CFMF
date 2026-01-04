library(Seurat)
library(parallel)
library(doParallel)
library(foreach)
library(tidyverse)


calculate.silhouette <- function(seu, purity, ncores = 10) {
  cl    <- makeCluster(ncores)
  registerDoParallel(cl)
  
  data  <- seu@reductions$pca@cell.embeddings |> as.matrix() |> dist()
  genes <- purity$gene
  
  df    <- foreach(i = genes, .combine = "rbind", .packages = c("Seurat", "cluster","dplyr")) %dopar% {
    
    expr              <- seu@assays$RNA@data[i, ]
    clusters          <- ifelse(expr > 0, 1, 2)
    
    silhouette.scores <- cluster::silhouette(x = as.numeric(clusters),dist = data)
    
    silhouette.df     <- as.data.frame(silhouette.scores) |> filter(cluster == 1)
    
    pos.silhouette    <- median(silhouette.df$sil_width)
    
    data.frame(gene = i, silhouette = pos.silhouette)
  }
   
  silhouette.df       <- merge(df, purity, by = "gene")
  
  stopCluster(cl)
  return(silhouette.df)
}
