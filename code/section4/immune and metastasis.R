library(Seurat)
library(harmony)
library(tidyverse)
library(ggplot2)
library(foreach)

setwd("/home/yinshumin/pro_cfmf/")

load("./output/section4/subclone/CFMF markers harmony high quality seu.RData")
immune.anno     <- readRDS("./output/section4/subclone/cell.anno/immune.anno.ls.rds")
cell.types      <- readRDS("./output/section4/subclone/cell.anno/coarse.cell.type.ls.rds")
cell.types      <- cell.types[c("Data_Colorectal_Lee2020", "Data_Colorectal_Pelka2021", "Data_Colorectal_Qian2020_normalized", "Data_Colorectal_Uhlitz2021_normalized")]
cell.types      <- lapply(cell.types, function(x) {df <- x %>% dplyr::select("sample", "my_coarse_cell_type")})
immune.anno     <- immune.anno[c("Data_Colorectal_Lee2020", "Data_Colorectal_Pelka2021", "Data_Colorectal_Qian2020_normalized", "Data_Colorectal_Uhlitz2021_normalized")]

immune.pro.ls   <- foreach(i=cell.types, j=immune.anno) %do% {
  i      <- i %>% rownames_to_column("cell_id")  
  j      <- j %>% rownames_to_column("cell_id")
  merged <- left_join(i, j, by = "cell_id", suffix = c("", "_b")) %>%
    mutate(immune.cell.type = coalesce(singleR_labels, my_coarse_cell_type)) %>% dplyr::select(-ends_with("_b"))  
}

immune.pro          <- Reduce(rbind, immune.pro.ls) %>% dplyr::select("sample", "immune.cell.type")


##meta##
seu <- readRDS("./input/CRC/Data_Colorectal_Che2021.rds")
meta.anno           <- seu@meta.data %>% dplyr::select(sample, singleR_labels)
meta.anno$singleR_labels[meta.anno$singleR_labels == "Epithelial cells"] <- "Malignant"
colnames(meta.anno) <- c("sample", "immune.cell.type")
rownames(meta.anno) <- NULL
immune.pro          <- rbind(immune.pro, meta.anno)

immune.pro.filtered <- immune.pro %>% filter(!immune.cell.type %in% c("Malignant", "Unknown"))
result              <- immune.pro.filtered %>% group_by(sample) %>% summarise(count = n(), 
                                                                  CD8T = sum(immune.cell.type == "CD8+ T cells"), 
                                                              CD8T.pro = CD8T/count ) %>% replace_na(list(cd8_prop = 0))



subclone.pro.df     <- seu.high@meta.data %>% group_by(sample, subclone.final) %>% dplyr::summarise(count = n()) %>% mutate(pro = count / sum(count)) 
subclone.pro.df     <- subclone.pro.df %>% ungroup() %>% dplyr::select(sample, subclone.final, pro)

samples             <- intersect(unique(result$sample), unique(subclone.pro.df$sample))
subclone.pro.df     <- subclone.pro.df %>% filter(sample %in% samples)
immune.pro.df       <- result %>% filter(sample %in% samples)



merged_df   <- subclone.pro.df %>% left_join(immune.pro.df %>% dplyr::select(sample, CD8T.pro), by = "sample")

cor.results <- merged_df %>% group_by(subclone.final) %>% summarise(n_samples = n(), 
                                                                    cor = cor(pro, CD8T.pro, method = "pearson", use = "complete.obs"),
                                                                    p_value = cor.test(pro, CD8T.pro, method = "pearson")$p.value)


save(list = c("subclone.pro.df", "immune.pro.df", "merged_df", "cor.results"), file = "./output/section4/subclone/subclone.pro.df and immune.pro.df.RData")




########metastasis#######
Roe      <- Startrac::calTissueDist(dat.tb = seu.high[[]], colname.cluster = "subclone.final", colname.tissue = "sample_origin") 
Roe      <- Roe[, c("primary", "metastasis")] |> as.matrix()
## +++, Ro/e > 1; ++, 0.8 < Ro/e ≤ 1; +, 0.2 ≤ Ro/e ≤ 0.8; +/−, 0 < Ro/e < 0.2; −, Ro/e = 0
cell_fun <- function(j, i, x, y, w, h, fill) {
  if (Roe[i, j] > 1) {
    grid.text("+++", x, y)
  } else if (Roe[i, j] > 0.8) {
    grid.text("++", x, y)
  } else if (Roe[i, j] >= 0.2) {
    grid.text("+", x, y)
  } else if (Roe[i, j] > 0) {
    grid.text("+/-", x, y)
  } else if (Roe[i, j] == 0) {
    grid.text("-", x, y)
  }
}

save("Roe", file = "./output/section4/subclone/Roe.RData")




