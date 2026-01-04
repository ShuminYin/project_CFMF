library(tidyverse)
library(parallel)
library(foreach)
library(doParallel)
library(Seurat)

calculate.purity <- function(seu, gene.pct, neighbor,ncores = 10) {
  cl       <- makeCluster(ncores)
  registerDoParallel(cl)
  
  data     <- seu@reductions$pca@cell.embeddings %>% as.matrix()
  neighbor <- BiocNeighbors::findKNN(data, k = neighbor, get.distance=TRUE)
  
  cutoff   <- 0
  genes    <- gene.pct$gene
  
  gene.purity <- foreach(i = genes, .combine = "rbind", .packages = c('Seurat', 'doParallel')) %dopar% {
    gene.exp  <- seu[["RNA"]]@data[i, ]
    clusters  <- gene.exp > cutoff
    
    cells     <- 1:ncol(seu)
    pos.cells <- which(clusters)
    
    df        <- foreach(j = pos.cells, .combine = 'rbind') %do% {
      other.cells <- cells[-j]
      is.neighbor <- other.cells %in% neighbor$index[j, ]  
      is.pos      <- other.cells %in% pos.cells
      
      fit         <- fisher.test(x = is.neighbor, y = is.pos)
      p.val       <- fit$p.value
      
      data.frame(cell = j, p.val, row.names = NULL)
    }
    
    chisq_stat   <- -2 * sum(log(df$p.val))
    v            <- 2 * length(df$p.val)
    meta.p       <- pchisq(chisq_stat, v, lower.tail = F)
    
    
    data.frame(gene = i, meta.p)
  }
  
  gene.purity <- gene.purity |> mutate(FDR = p.adjust(meta.p, method = "BH"))
  gene.purity <- merge(gene.purity, gene.pct, by = "gene")
  
  stopCluster(cl)
  return(gene.purity)
}


