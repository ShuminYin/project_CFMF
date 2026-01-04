library(Seurat)
library(harmony)
library(tidyverse)
library(parallel)
library(doParallel)
library(foreach)
library(ggplot2)
library(patchwork)

setwd("/home/yinshumin/pro_cfmf/")

#################1.CFMF of CRC dataset###############################
CRC.name <- dir("./input/CRC/")
CRC.name <- CRC.name[1:6]
CRC.name <- gsub(".rds","",CRC.name)

lapply(CRC.name, function(CRC) {
  seu <- readRDS(paste0("./input/CRC/",CRC,".rds"))
  
  seu.ls        <- foreach(i = seu.ls) %do% {
    seu <- i %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% RunUMAP(dims = 1:30)
  }
  names(seu.ls) <- seu.name
  
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
  
  save(list = c("seu.ls","silhouette.ls"),file = paste0("./output/section4/",CRC,"_seu_and_silhouette_ls.RData"))
  
})

#############2. integeration of dataset#######
load("./output/section4/Data_Colorectal_Johannes2023.RData")
seu.Joh.ls        <- seu.ls
silhouette.Joh.ls <- silhouette.ls


load("./output/section4/Data_Colorectal_Lee2020.RData")
seu.Lee.ls        <- seu.ls
silhouette.Lee.ls <- silhouette.ls

load("./output/section4/Data_Colorectal_Qian2020_normalized.RData")
seu.Qia.ls        <- seu.ls
silhouette.Qia.ls <- silhouette.ls


load("./output/section4/Data_Colorectal_Pelka2021.RData")
seu.Pel.ls        <- seu.ls
silhouette.Pel.ls <- silhouette.ls


load("./output/section4/Data_Colorectal_Uhlitz2021_normalized.RData")
seu.Uhl.ls        <- seu.ls
silhouette.Uhl.ls <- silhouette.ls


load("./output/section4/Data_Colorectal_Che2021.RData")
seu.GSE.ls <- seu.ls
silhouette.GSE.ls <- silhouette.ls


seu.ls.all <- c(seu.meta.ls, seu.lee.ls, seu.Pel.ls, seu.Qia.ls, seu.Uhl.ls, seu.GSE.ls)
silhouette.ls.all <- c(silhouette.meta.ls, silhouette.lee.ls, silhouette.Pel.ls, silhouette.Qia.ls, silhouette.Uhl.ls, silhouette.GSE.ls)

save(list = c("seu.ls.all", "silhouette.ls.all"),
     file = "./output/section4/silhouette and seu of all sample.RData")


seu.malig.ls <- Filter(function(seu) {
  "sample" %in% colnames(seu@meta.data) && 
    any(seu$sample_origin %in% c("metastasis", "primary","metastasis (liver)"))
}, seu.ls.all)

silhouette.malig.ls <- silhouette.ls.all[names(seu.malig.ls)]

seu  <- Reduce(function(x, y) merge(x, y), seu.malig.ls)

seu$pri_meta <- ifelse(seu$sample_origin == "primary","primary","metastasis")


save(list = c("seu","silhouette.malig.ls","seu.malig.ls"), 
     file = "./output/section4/merge seu of malignant sample and malignant silhouette.ls.RData")


#############3. HVG:  harmony of dataset#######
seu.malignant <- subset(seu, malignant == "yes")
seu.malignant@meta.data <- seu.malignant@meta.data[c("orig.ident", "nCount_RNA", "nFeature_RNA", "cell_type", "sample", "sample_origin", "malignant","pri_meta")]
seu.malignant <- seu.malignant %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  RunHarmony(group.by.vars = "sample") %>% RunUMAP(reduction = "harmony", dims = 1:30) %>% FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters()

save("seu.malignant", 
     file = "./output/section4/harmony seu.malignant of tumor sample.RData")

rm(list = ls());gc()


#############4. CFMF:  harmony of dataset#######
load("./output/section4/harmony seu.malignant of tumor sample.RData")
load("./output/section4/merge seu of malignant sample and malignant silhouette.ls.RData")

malig.pct.ls <- foreach(seu = seu.malig.ls, silhouette = silhouette.malig.ls) %do% {
  genes          <- silhouette %>% filter(-log10(FDR) > 200, silhouette > 0) %>% pull(gene)
  
  if (length(genes) == 0) {
    return(NULL)
  }
  
  data           <- GetAssayData(seu, slot = "data")[genes, ] %>% as.matrix()
  expr.mt        <- data > 0
  cell.type      <- seu$malignant
  
  cancer.gene.df   <- data.frame(gene = character(), malig.in.pos = numeric(), malig.in.can = numeric())
  cancer.count     <- ncol(subset(seu, malignant == "yes"))
  
  for (gene in genes) {
    pos.cells      <- which(expr.mt[gene, ] == 1)
    pos.cell.types <- cell.type[pos.cells]
    cancer.ratio   <- sum(pos.cell.types == "yes") / length(pos.cells)
    cancer.pct     <- sum(pos.cell.types == "yes") / cancer.count
    cancer.gene.df <- rbind(cancer.gene.df, data.frame(gene = gene, malig.in.pos = cancer.ratio,
                                                       malig.in.can = cancer.pct))
  }
  
  cancer.gene.df <- merge(cancer.gene.df, silhouette, by = "gene")
  
  return(cancer.gene.df)
}
names(malig.pct.ls) <- names(seu.malig.ls)

save("malig.pct.ls", file = "./output/section4/malig.pct.ls.RData")


malig.pct.ls <- Filter(Negate(is.null), malig.pct.ls)

density.top.ls <- lapply(malig.pct.ls, function(x) {
  density.df <- x %>% filter(malig.in.pos > 0.5, malig.in.can < 0.5)
})

desity.top1000.ls                 <- lapply(density.top.ls, function(x) {gene <- x %>% arrange(desc(silhouette)) %>% head(1000) %>% pull(gene)})
desity.top1000.count.df           <- as.data.frame(table(unlist(desity.top1000.ls)))
colnames(desity.top1000.count.df) <- c("gene", "count")

save(list = c("desity.top1000.ls","desity.top1000.count.df"),
     file = "./output/section4/top1000 gene count.RData")




load("./output/section4/harmony seu.malignant of tumor sample.RData")

HVGs <- desity.top1000.count.df %>% arrange(desc(count)) %>% head(2000) %>% pull(gene) %>% as.character()
seu  <- seu.malignant %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData(features = rownames(seu.malignant)) %>% RunPCA(features = HVGs) %>% 
  RunHarmony(group.by.vars = "sample") %>% RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% FindClusters()

save("seu", 
     file = "./output/section4/CFMF markers harmony seu.malignant of tumor sample.RData")



#######5. low quality clusters######
seu$log10_nCount   <- log10(seu$nCount_RNA)
seu$log10_nFeature <- log10(seu$nFeature_RNA)
clusters           <- seu$seurat_clusters
avg.nCount         <- tapply(seu$log10_nCount, clusters, mean)
avg.nFeature       <- tapply(seu$log10_nFeature, clusters, mean)
avg.nCount.df      <- data.frame(Cluster      = names(avg.nCount),
                                 Avg.nCount   = as.numeric(avg.nCount),
                                 Avg.nFeature = as.numeric(avg.nFeature))

Idents(seu) <- "seurat_clusters"
p1 <- FeaturePlot(seu, c("nCount_RNA", "nFeature_RNA"),ncol = 1,label = T, raster = F) & scale_color_gradient(trans = "log10", low = "lightgrey", high = "blue")
p2 <- VlnPlot(seu, features = c("nCount_RNA", "nFeature_RNA"),ncol = 1, pt.size = 0, raster = F)
p1|p2

seu$quality <- ifelse(seu$seurat_clusters %in% c(3,4,5,17,18,19), "low", "high")
seu.high    <- subset(seu, quality == "high")
seu.high    <- seu.high %>% NormalizeData() %>% ScaleData(features = rownames(seu.high)) %>% RunPCA(features = HVGs) %>% 
  RunHarmony(group.by.vars = "sample", reduction.use = "pca") %>% FindNeighbors(reduction = "harmony",dims = 1:30) %>% FindClusters() %>% RunUMAP(reduction = "harmony",dims = 1:30)

save("seu.high", file = "./output/section4/subclone/CFMF markers harmony high quality seu.RData")



#############subclone#########
seu.high@meta.data <- seu.high@meta.data %>%
  mutate(
    subclone = case_when(
      seurat_clusters %in% c(2,6,12,13)       ~ "subclone1",
      seurat_clusters %in% c(11)          ~ "subclone2",
      seurat_clusters %in% c(10)         ~ "subclone3",
      seurat_clusters %in% c(0,18,19)        ~ "subclone4",
      seurat_clusters %in% c(4,7,9)        ~ "subclone5",
      seurat_clusters %in% c(5,15,17)       ~ "subclone6",
      seurat_clusters %in% c(1,14)        ~ "subclone7",
      seurat_clusters %in% c(3,8,16)        ~ "subclone8"
    )
  )
seu.high$subclone <- factor(seu.high$subclone, levels = c("subclone1","subclone2","subclone3","subclone4","subclone5","subclone6","subclone7","subclone8"))


seu.high@meta.data <- seu.high@meta.data %>%
  mutate(
    subclone.final = case_when(
      seurat_clusters %in% c(2,6,12,13)       ~ "subclone1",
      seurat_clusters %in% c(11)          ~ "subclone2",
      seurat_clusters %in% c(10)         ~ "subclone3",
      seurat_clusters %in% c(0,18,19)        ~ "subclone4",
      seurat_clusters %in% c(4,7,9)        ~ "subclone5",
      seurat_clusters %in% c(5,15,17)       ~ "subclone6",
      seurat_clusters %in% c(1,14,3,8,16)        ~ "subclone7"
   
    )
  )
seu.high$subclone.final <- factor(seu.high$subclone.final, levels = c("subclone1","subclone2","subclone3","subclone4","subclone5","subclone6","subclone7"))


save("seu.high", file = "./output/section4/subclone/CFMF markers harmony high quality seu.RData")


out         <- scran::scoreMarkers(x = seu.high[["RNA"]]@data, groups = seu.high$subclone.final)
markers.ls  <- lapply(out, \(x) x |> as.data.frame() |> filter(min.AUC > 0.5) |> slice_max(order_by = min.AUC, n = 50) |> rownames())
save(list = c("out","markers.ls"), file = "./output/section4/subclone/scoreMarkers of seu.high subclones.RData")





