.libPaths("/home/yinshumin/R/x86_64-pc-linux-gnu-library/4.2/")  ## necessary for finding the packages used in the script.
library(Seurat)
library(foreach)
library(tibble)
library(magrittr)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(purrr)
library(reldist)
library(entropy)

setwd("/home/yinshumin/pro_cfmf/")

############################genes of 7 method in pbmc3k###########################
seu           <- readRDS("./input/pbmc/pbmc3k.rds")

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


CFMF.genes <- silhouette %>%
  arrange(FDR, desc(silhouette)) %>%
  head(2000) %>%
  pull(gene) %>%
  as.character()



gene.methods <- list(hvg.genes, m3drop.genes, SCT.genes, Fano.genes, RaceID.genes, SE.genes, CFMF.genes)
names(gene.methods) <- c('hvg.genes', 'm3drop.genes', 'SCT.genes', 'Fano.genes', 'RaceID.genes', 'SE.genes', 'CFMF.genes')

save("gene.methods", file =  "./output/section5/pbmc3k: gene of 7 mothed.RData")


AUC.ls <- foreach(i = gene.methods) %do% {
  genes.name  <- i
  cell.type   <- unique(seu$cell_type)
  
  AUC.df <- foreach(j = cell.type, .combine = "cbind") %do% {
    seu$test    <- ifelse(seu$cell_type == j, j, "other")
    
    out         <- scran::scoreMarkers(x = seu[["RNA"]]@data, groups = seu$test)
    markers.AUC <- out[[j]] %>% as.data.frame() %>% dplyr::select(min.AUC) %>% rownames_to_column(var = "gene")
    AUC.df      <- markers.AUC %>% filter(gene %in% genes.name) %>% column_to_rownames(var = "gene")
    
  }
  
  return(AUC.df)
}

names(AUC.ls) <- names(gene.methods)


save("AUC.ls", file =  "./output/section5/pbmc3k:AUC.ls.RData")




########################replacement:silhouette in pbmc3k########################
load("./output/section5/pbmc3k: gene of 7 mothed.RData")
seu <- readRDS("./input/pbmc/pbmc3k.rds")

seu.ls <- lapply(gene.methods, function(x) {
  seu <- seu |>  NormalizeData() |> ScaleData(features = rownames(seu)) |> RunPCA(features = x) |> RunUMAP(dims = 1:30)
})
names(seu.ls) <- names(gene.methods)


seu.feature.gene.ls <- lapply(seu.ls, function(i) {
  Idents(i)              <- "cell_type"
  data                   <- i@reductions$pca@cell.embeddings |> as.matrix() |> dist()
  clusters               <- as.integer(as.factor(Idents(i)))
  silhouette.scores      <- cluster::silhouette(x = clusters, dist = data)
  silhouette.df          <- as.data.frame(silhouette.scores)[, c("cluster", "sil_width")]
  rownames(silhouette.df)<- rownames(i@meta.data)
  i$silhouette           <- silhouette.df$sil_width
  return(i)
})

save("seu.feature.gene.ls", file = "./output/section5/seu.ls_of_7_methods_with_silhouette.RData")


##############################singler:sampling of pbmc3k########################
seu <- readRDS("./input/pbmc/pbmc3k.rds")

set.seed(1234)
seeds         <- sample(1000, 50)

seu.train.ls <- foreach(i = seeds, .packages = c("Seurat","dplyr")) %dopar% {
  set.seed(i)
  sample.size  <- floor(0.7 * length(colnames(seu)))
  sample.cells <- sample(colnames(seu), size = sample.size)
  seu.train <- subset(seu, cells = sample.cells)
  seu.train <- seu.train %>% NormalizeData() %>% FindVariableFeatures()
}

seu.test.ls <- lapply(seu.train.ls, function(seu.train) {
  sample.cells <- colnames(seu.train)
  seu.test     <- subset(seu, cells = base::setdiff(colnames(seu), sample.cells))
})

save(list = c("seu.train.ls", "seu.test.ls"), file = "./output/section5/50 times sample seu.RData")

###cfmf###
CFMF.silhouette.ls <- foreach(seu.train = seu.train.ls) %dopar% {
  cutoff    <- 0 
  gene.pct  <- data.frame(gene = rownames(seu.train), pct = rowMeans(seu.train[["RNA"]]@data > cutoff), 
                          count = rowSums(seu.train[["RNA"]]@data > cutoff))  
  gene.pct  <- gene.pct %>% filter(count > 3, pct < 1)
  
  
  source("/home/yinshumin/pro_cfmf/code/1. formula/FDR.function.R")
  purity     <- calculate.purity(seu.train, gene.pct = gene.pct, neighbor = 50, ncores = 20)
  
  source("/home/yinshumin/pro_cfmf/code/1. formula/Silhouette.function.R")
  silhouette <- calculate.silhouette(seu.train, purity, ncores = 20)
}

save(list = c("CFMF.silhouette.ls"), file = "./output/section5/50 times sample silhouette.RData")


###genes of 7 mothod###
method.gene.ls <- foreach(i = 1:50, .packages = c("Seurat","dplyr")) %do% {
  seu.train  <- seu.train.ls[[i]]
  silhouette <- CFMF.silhouette.ls[[i]]
  counts     <- seu.train@assays$RNA@counts %>% as.matrix()
  sda        <- t(counts)
  
  #######1. HVG########
  seu.hvg   <- seu.train %>% NormalizeData() %>% FindVariableFeatures(nfeatures = 2000)
  hvg.genes <- VariableFeatures(seu.hvg)
  
  #######2. M3Drop########
  norm          <- M3DropConvertData(counts, is.counts=TRUE) |> as("dgCMatrix")
  genes         <- M3DropFeatureSelection(norm, suppress.plot = T, mt_threshold = 2)
  m3drop.genes  <- genes %>% arrange(p.value, desc(effect.size)) %>% head(2000) %>% pull(Gene)
  
  
  #######3. SCTransform########
  sobj          <- SCTransform(seu.train, verbose = FALSE, variable.features.n = 2000)
  SCT.genes     <- VariableFeatures(sobj)
  
  #######4. Fano Factor########
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
  
  #######5. CFMF########
  CFMF.genes    <- silhouette %>% arrange(FDR, desc(silhouette)) %>% head(2000) %>% pull(gene) %>% as.character()
  
  #######6. SE###########
  ent.res       <- ROGUE::SE_fun(counts)
  SE.genes      <- ent.res %>% arrange(p.value, desc(entropy)) %>% head(2000) %>% pull(Gene)
  
  
  #######7. RaceID3###########
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
  
  
  df <- data.frame(hvg = hvg.genes, m3drop = m3drop.genes, SCT = SCT.genes,
                   Fano = Fano.genes, CFMF = CFMF.genes, SE = SE.genes, RaceID3 = RaceID.genes)
}


save(list = c("method.gene.ls"), file = "./output/section5/feature_genes_of_7method_of_50pbmc.RData")



accuracy.ls <- foreach(i = 1:50, .packages = c("Seurat","dplyr","foreach")) %do% {
  
  method = c('hvg', 'm3drop', 'SCT', 'Fano', 'CFMF', 'SE', 'RaceID3')
  
  gene.7.method <- method.gene.ls[[i]]
  seu.train     <- seu.train.ls[[i]]
  seu.test      <- seu.test.ls[[i]]
  
  gene.accuracy.df <- foreach(j = method, .packages = c("Seurat","dplyr"), .combine = "cbind") %do% {
    
    genes         <- gene.7.method[[j]]
    
    ref           <- AggregateExpression(seu.train, group.by = "cell_type")$RNA
    colnames(ref) <- as.character(colnames(ref))
    ref           <- apply(ref, 2, \(x) log2(1e6 * x / sum(x) + 1))
    
    pred.HVGs     <- SingleR::SingleR(test = seu.test[["RNA"]]@counts, ref = ref, labels = colnames(ref), restrict = genes)
    pct.HVGs      <- sum(seu.test$cell_type == pred.HVGs$pruned.labels, na.rm = TRUE) / ncol(seu.test)
    
    data.frame(method = pct.HVGs)
  }
  colnames(gene.accuracy.df) <- c('hvg', 'm3drop', 'SCT', 'Fano', 'CFMF', 'SE', 'RaceID3')
  
  return(gene.accuracy.df)
}


save(list = c("accuracy.ls"), file = "./output/section5/accuracy_of_50pbmc_compared_7methods.RData")




##################################hotspot#######################################
library(sceasy)
library(reticulate)
library(readr)
load("./output/section5/pbmc3k: gene of 7 mothed.RData")
seu <- readRDS("./input/pbmc/pbmc3k.rds")

# #####hotspot in python#############
# conda create -n hotspot_env python=3.9 -y


use_condaenv("hotspot_env", required = TRUE)

genes <- unique(unlist(gene.methods))

convertFormat(seu[genes, ], from = "seurat", to = "anndata", main_layer = "data", transfer_layers = c("counts"), drop_single_values = FALSE, 
              outFile = "./output/section5/adata.h5ad")


# conda activate hotspot_env
# python
# import scanpy as sc
# import hotspot
# adata_path = "adata.h5ad"
# adata = sc.read_h5ad(adata_path)
# hs = hotspot.Hotspot(adata, model='normal', latent_obsm_key='X_pca')
# hs.create_knn_graph(n_neighbors=30, weighted_graph=False)
# hs_results = hs.compute_autocorrelations()
# output_file = "./output/section5/hotspot_all_genes_results.csv"
# hs_results.to_csv(output_file)





