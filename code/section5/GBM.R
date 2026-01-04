.libPaths("/home/yinshumin/R/x86_64-pc-linux-gnu-library/4.2/")  ## necessary for finding the packages used in the script.
library(Seurat)
library(tidyverse)
library(foreach)
library(ggplot2)
library(patchwork)
library(MatrixGenerics) 
library(scran)
library(matrixStats)
library(RaceID)
library(M3Drop)


setwd("/home/yinshumin/pro_cfmf/")


##############AUC of GBM dataset#####################################
GBM.names <- dir("./input/GBM/")
GBM.names <- GBM.names[c(1,2,4)]
GBM.names <- gsub(".rds", "", GBM.names)
GBM.names <- gsub("Data_Brain_", "", GBM.names)

lapply(GBM.names, function(GBM) {
  seu.all <- readRDS(paste0("./input/GBM/Data_Brain_", GBM, ".rds"))
  seu.ls  <- SplitObject(seu.all, split.by = "sample")
  save("seu.ls",file = paste0("./output/section5/",GBM,".seu.ls.RData"))
})

#####genes of 7 methods

foreach(GBM = GBM.names) %do% {
  load(paste0("./output/section5/",GBM,".seu.ls.RData"))
  
  gene.methods7.ls <- foreach(j = seu.ls) %do% {
    seu           <- j
    
    counts        <- seu@assays$RNA@counts %>% as.matrix()
    sda           <- t(counts)
    
    seu           <- seu %>% NormalizeData() %>% FindVariableFeatures(nfeatures = 2000)
    hvg.genes     <- VariableFeatures(seu)
    
    
    library("M3Drop")
    norm          <- M3DropConvertData(counts, is.counts=TRUE) |> as("dgCMatrix")
    genes         <- M3DropFeatureSelection(norm, suppress.plot = T, mt_threshold = 2)
    m3drop.genes  <- genes %>% arrange(p.value, desc(effect.size)) %>% head(2000) %>% pull(Gene)
    
    
    ######
    sobj          <- SCTransform(seu, verbose = FALSE, variable.features.n = 2000)
    SCT.genes     <- VariableFeatures(sobj)
    
    
    FanoFactor_fun <- function(expr){
      Fano <- apply(expr,2,function(x) var(x)/mean(x))
      Fano <- Fano %>% 
        as.data.frame() %>%
        tibble::rownames_to_column(var = "Gene")
      
      colnames(Fano) <- c("Gene","fano")
      Fano <- Fano %>% dplyr::filter(!is.na(fano)) %>% dplyr::arrange(desc(fano))
      
      Fano <- Fano %>% 
        dplyr::mutate(
          p.value = 1-pnorm(
            Fano$fano,
            mean = mean(Fano$fano),
            sd = sd(Fano$fano)))
      
      return(Fano)
    }
    Fano.genes     <- FanoFactor_fun(sda)
    Fano.genes     <- Fano.genes %>% arrange(p.value, desc(fano)) %>% head(2000) %>% pull(Gene)
    
    
    raceid_fun <- function(x, mthr = -1){
      uvar  <- function(x,fit){
        err <- coef(summary(fit))[, "Std. Error"]
        2**(coef(fit)[1] + err[1] + log2(x)*(coef(fit)[2] + err[2]) + (coef(fit)[3] + err[3]) * log2(x)**2)
      }
      
      x <- t(x)
      m <- apply(x, 1, mean)
      v <- apply(x, 1, var)
      ml <- log2(m)
      vl <- log2(v)
      f <- ml > -Inf & vl > -Inf
      ml <- ml[f]
      vl <- vl[f]
      mm <- -8
      repeat {
        fit <- lm(vl ~ ml + I(ml^2))
        if (coef(fit)[3] >= 0 | mm >= mthr) {
          break
        }
        mm <- mm + 0.5
        f <- ml > mm
        ml <- ml[f]
        vl <- vl[f]
      }
      vln <- log2(v) - log2(sapply(m, FUN = uvar, fit = fit))
      tibble(
        Gene = names(vln),
        vln = vln
      ) %>%
        dplyr::filter(!is.na(vln)) %>%
        dplyr::arrange(desc(vln)) -> race.res
      
      race.res <- race.res %>% 
        dplyr::mutate(
          p.value = 1-pnorm(
            race.res$vln,
            mean = mean(race.res$vln),
            sd = sd(race.res$vln)))
      return(race.res)
    }
    RaceID.genes  <- raceid_fun(sda)
    RaceID.genes  <- RaceID.genes %>% arrange(p.value, desc(vln)) %>% head(2000) %>% pull(Gene)
    
    
    library(ROGUE)
    ent.res       <- SE_fun(counts)
    SE.genes      <- ent.res %>% arrange(p.value, desc(entropy)) %>% head(2000) %>% pull(Gene)
    
    
    ####CFMF###
    seu <- seu %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>%  RunUMAP(dims = 1:30)
    cutoff    <- 0 
    gene.pct  <- data.frame(gene = rownames(seu), pct = rowMeans(seu[["RNA"]]@data > cutoff), 
                            count = rowSums(seu[["RNA"]]@data > cutoff))  
    gene.pct  <- gene.pct %>% filter(count > 3, pct < 1)
    
    source("./code/section1/FDR.function.R")
    purity     <- calculate.purity(seu, gene.pct = gene.pct, neighbor = 50, ncores = 30)
    
    source("./code/section1/Silhouette.function.R")
    silhouette <- calculate.silhouette(seu, purity = purity, ncores = 30)
    
    
    CFMF.genes    <- silhouette %>%
      arrange(FDR, desc(silhouette)) %>%
      head(2000) %>%
      pull(gene) %>%
      as.character()
    

    gene.methods <- list(hvg.genes, m3drop.genes, SCT.genes, Fano.genes, RaceID.genes, SE.genes, CFMF.genes)
    names(gene.methods) <- c('hvg.genes', 'm3drop.genes', 'SCT.genes', 'Fano.genes', 'RaceID.genes', 'SE.genes', 'CFMF.genes')
    
    return(gene.methods)
  }
  
  names(gene.methods7.ls) <- names(seu.ls)
  
  save("gene.methods7.ls", file = paste0("./output/section5/", GBM, "_genes_of_7methods.RData"))
  
}



#####AUC of genes
Abdelfattah2022_GBM <- c("rGBM-01-A","ndGBM-11-C","ndGBM-02-5", "rGBM-03-3", "rGBM-05-3" , "ndGBM-07")
Wang2019            <- c("SF10127","SF11979sn","SF12090","SF4400","SF6996","SF9259S")
Neftel2019_10x      <- c("102","114","115","118","124","125","126","105")


foreach(GBM = GBM.names) %do% {
  load(paste0("./output/section5/",GBM,".seu.ls.RData"))
  load(paste0("./output/section5/", GBM, "_genes_of_7methods.RData"))
  samples <- get(GBM)
  seu.ls  <- seu.ls[samples]
  
  auc.sample.ls <- foreach(j = names(seu.ls)) %do% {
    seu     <- seu.ls[[j]]
    gene.ls <- gene.methods7.ls[[j]]
    
    AUC.method.ls <- foreach(g = gene.ls) %do% {
      cell.type   <- unique(seu$cell_type)
      cell.type   <- setdiff(cell.type, "")
      
      AUC.celltype.df <- foreach(c = cell.type, .combine = "cbind") %do% {
        seu$test         <- ifelse(seu$cell_type == c, c, "other")
        
        out              <- scran::scoreMarkers(x = seu[["RNA"]]@data, groups = seu$test)
        markers.AUC      <- out[[c]] %>% as.data.frame() %>% dplyr::select(min.AUC) %>% rownames_to_column(var = "gene")
        AUC.celltype.df  <- markers.AUC %>% filter(gene %in% g) %>% column_to_rownames(var = "gene")
        
      }
      
      rownames(AUC.celltype.df) <- NULL
      AUC.max.df                <- apply(AUC.celltype.df, 1, max) %>% as.data.frame() %>% setNames("AUC")
      
      return(AUC.max.df)
    }
    
    
    names(AUC.method.ls) <- names(gene.ls)
 
    
    return(AUC.method.ls)
  }
  
  names(auc.sample.ls) <- samples
  
  save("auc.sample.ls",file = paste0("./output/section5/", GBM, "_AUC_of_7methods.RData"))
}


