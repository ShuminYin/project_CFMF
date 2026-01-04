library(Seurat)
library(foreach)
library(tidyverse)
library(ggplot2)
library(patchwork)



setwd("/home/yinshumin/pro_cfmf/")

load("./input/rare_cell/hlca.seu.RData")


filtered_samples <- seu@meta.data %>%
  dplyr::select(sample, ann_level_4) %>%
  filter(ann_level_4 %in% c("Neuroendocrine", "Ionocyte")) %>%
  group_by(sample, ann_level_4) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(sample) %>%
  filter(any(count > 3)) %>%
  distinct(sample) %>% pull(sample) %>% as.character() %>% unique()

seu.rare <- subset(seu, sample %in% filtered_samples)

save("seu.rare", file = "./output/section3/seu.rare.all.RData")

df.neuro <- seu.rare@meta.data %>%
  group_by(sample) %>%
  mutate(total = n()) %>%
  ungroup() %>%
  filter(ann_level_4 == "Neuroendocrine") %>%
  group_by(sample, total) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(proportion = count / total)

df.Iono <- seu.rare@meta.data %>%
  group_by(sample) %>%
  mutate(total = n()) %>%
  ungroup() %>%
  filter(ann_level_4 == "Ionocyte") %>%
  group_by(sample, total) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(proportion = count / total)

save(list = c("df.neuro","df.Iono"), file = "./output/section3/prop.df of neuro and Iono.RData")

#sample:GRO-03_biopsy
#rare.pct:0.02%

######################CFMF##########################
load("./input/rare_cell/hlca.seu.RData")

seu.ls      <- SplitObject(seu, split.by = "sample")
seu.name    <- names(seu.ls)

#####0. FDR and silhouette####
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
  purity     <- calculate.purity(i, gene.pct = j, neighbor = 50, ncores = 60)
}
names(purity.ls) <- names(seu.ls)


source("./code/section1/Silhouette.function.R")
silhouette.ls <- foreach(i = seu.ls, j = purity.ls) %do% {
  silhouette <- calculate.silhouette(i, purity = j, ncores = 60)
}
names(silhouette.ls) <- names(seu.ls)

save(list = c("seu.ls","silhouette.ls"),
     file = "./output/section3/seu and silhouette of samples.RData")


#################marker.pct#######################
df.ls <- lapply(silhouette.ls, function(df) {
  df <- df %>% filter(-log10(FDR) > 200, silhouette > 0)
})

targets <- c("SCG3", "BSND", "ASCL3", 'CLEC3B', 'CLCNKB', 'ADGRF5', 'CLNK',
             'USH1C',  'TAGLN3',   'PCSK1N',  'FGF14',  'PROX1', 'SCG2',  'CHGA')

res <- sum(sapply(df.ls, function(x) any(targets %in% x[[1]])))

res/length(names(seu.ls))



