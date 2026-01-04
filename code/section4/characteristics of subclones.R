library(magrittr)
library(ComplexHeatmap)
library(Seurat)
library(ggplot2)
library(tidyverse)

setwd("/home/yinshumin/pro_cfmf/")

load("./out/5. cancer subclone/subclone/CFMF markers harmony high quality seu.RData")




########1. hallmark heatmap#########
hallmarks        <- msigdbr::msigdbr(category = "H") %$% split(gene_symbol, gs_name)
names(hallmarks) <- names(hallmarks) |> str_sub(start = 10) |> str_replace_all("_", " ")


AUC.score <- AUCell::AUCell_run(exprMat    = seu.high[["RNA"]]@counts,
                                geneSets   = hallmarks,
                                aucMaxRank = ceiling(0.05 * nrow(seu.high)),
                                BPPARAM    = BiocParallel::MulticoreParam(20))
AUC.score <- AUCell::getAUC(AUC.score)


ls       <- scran::scoreMarkers(x = AUC.score, groups = seu.high$subclone)
top      <- lapply(ls, \(x) x |> as.data.frame() |> filter(median.AUC > 0.5) |> slice_max(order_by = median.AUC, n = 5) |> rownames())
features <- unique(unlist(top))


set.seed(123)
cells   <- split(colnames(seu.high), seu.high$subclone) |> lapply(sample, 1000) |> unlist()
coldata <- seu.high[[]] |> dplyr::select(cluster = subclone)
coldata <- coldata[cells, , drop = FALSE]
mat     <- AUC.score[features, cells] |> t() |> scale() |> t()

save(list = c("coldata","mat"), file = "./output/section4/subclone/hallmark.RData")



#######2. reactome#########
reactome.path <- read.table("/home/yinshumin/Proj_Marker/input/reactome/Reactome_Pathways_2024.txt", header = F, sep = "\t", stringsAsFactors = F, fill = T)
ls            <- GSEABase::getGmt("/home/yinshumin/Proj_Marker/input/reactome/c2.cp.reactome.v2024.1.Hs.symbols.gmt") |> GSEABase::geneIds()
reactome.ls   <- Filter(function(x) length(x) > 10, ls)


AUC.score <- AUCell::AUCell_run(exprMat    = seu.high[["RNA"]]@counts,
                                geneSets   = reactome.ls,
                                aucMaxRank = ceiling(0.05 * nrow(seu.high)),
                                BPPARAM    = BiocParallel::MulticoreParam(40))
AUC.score <- AUCell::getAUC(AUC.score)


ls       <- scran::scoreMarkers(x = AUC.score, groups = seu.high$subclone)
top      <- lapply(ls, \(x) x |> as.data.frame() |> filter(median.AUC > 0.5) |> slice_max(order_by = median.AUC, n = 5) |> rownames())
features <- unique(unlist(top))


set.seed(123)
cells   <- split(colnames(seu.high), seu.high$subclone) |> lapply(sample, 1000) |> unlist()
coldata <- seu.high[[]] |> dplyr::select(cluster = subclone)
coldata <- coldata[cells, , drop = FALSE]


mat           <- AUC.score[features, cells] |> t() |> scale() |> t()
rownames(mat) <- gsub("REACTOME_", "", rownames(mat)) %>% str_replace_all("_", " ")
rownames(mat) <- rownames(mat) |> str_wrap(width = 45)

save(list = c("coldata","mat"), file = "./output/section4/subclone/reactiome.RData")
