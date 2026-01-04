.libPaths("/home/fuchunyang/R/x86_64-pc-linux-gnu-library/4.2/")  ## necessary for finding the packages used in the script.
library(spacexr)
library(foreach)
library(Seurat)
library(tidyverse)
library(readxl)


setwd("/home/yinshumin/pro_cfmf/")


#################RCTD of CRC samples########################################
sample.meta <- readxl::read_excel("./input/CRC_spatial/41586_2024_8087_MOESM3_ESM_modified.xlsx")
samples     <- sample.meta %>% filter(cancer_type == "CRC", sample_type_tumor == "Metastasis") %>% pull(LibraryName)


load("./output/section4/subclone/subclone/CFMF markers harmony high quality seu.RData")
cells        <- colnames(seu.high)
set.seed(1234)
cells.sample <- sample(cells, size = 5000, replace = F)
ref          <- subset(seu.high, cells = cells.sample)


# extract information to pass to the RCTD Reference function
counts         <- ref@assays$RNA@counts %>% as.matrix()
counts         <- round(counts)
cluster        <- as.factor(ref$subclone.final)
names(cluster) <- colnames(ref)
nUMI           <- round(ref$nCount_RNA)
names(nUMI)    <- colnames(ref)
reference      <- Reference(counts, cluster, nUMI)


prefix <- "./output/section4/spatial/RCTD/"
foreach(i = samples) %do% {
  tryCatch(expr = {
    
    load(paste0("./input/CRC_spatial/", i, ".RData"))
    counts                          <- seu[["Spatial"]]@counts
    coords                          <- GetTissueCoordinates(seu)
    colnames(coords)                <- c("x", "y")
    coords[is.na(colnames(coords))] <- NULL
    query                           <- SpatialRNA(coords, counts, colSums(counts))
    
    RCTD           <- create.RCTD(query, reference, max_cores = 50)
    RCTD           <- run.RCTD(RCTD, doublet_mode = "full")
    cell.propotion <- normalize_weights(RCTD@results$weights)
    
    dir.create(file.path(prefix, i))
    saveRDS(RCTD, file = file.path(prefix, i, "RCTD.rds"))
    saveRDS(cell.propotion, file = file.path(prefix, i, "cell.propotion.rds"))
    
  })
}


## k-means clustering
prefix      <- "./input/CRC_spatial/"
sample.meta <- readxl::read_excel("./input/CRC_spatial/41586_2024_8087_MOESM3_ESM_modified.xlsx")
samples     <- sample.meta %>% filter(cancer_type == "CRC", sample_type_tumor == "Metastasis") %>% pull(LibraryName)
files       <- paste0(samples, ".RData")

purity.kmeans.cl.ls        <- foreach(i = files) %do% {
  tryCatch(
    expr = {
      load(file.path(prefix, i))
      cl   <- kmeans(t(as.matrix(seu[["SCT"]]@data[VariableFeatures(seu), ])), centers = 5)
      cl   <- data.frame(kmeans = cl$cluster) |> rownames_to_column("barcode")
      meta <- seu[[]] |> rownames_to_column("barcode")
      
      meta |> select(barcode, purity = Tumor_purity_estimate) |> full_join(cl, by = "barcode")
    },
    
    error = function(e) e$message
  )
}
names(purity.kmeans.cl.ls) <- samples


prefix         <- "./output/section4/spatial/RCTD/"
samples        <- dir(prefix)
rctd.ls        <- foreach(i = samples) %do% {readRDS(file.path(prefix, i, "cell.propotion.rds")) |> as.data.frame() |> rownames_to_column("barcode") }
names(rctd.ls) <- samples

save(list = c("purity.kmeans.cl.ls","rctd.ls"),
     file = "./output/section4/spatial/RCTD.plot.RData")

