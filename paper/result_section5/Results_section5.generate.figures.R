.libPaths("/home/yinshumin/R/x86_64-pc-linux-gnu-library/4.2/")  ## necessary for finding the packages used in the script.
library(Seurat)
library(foreach)
library(tibble)
library(magrittr)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(purrr)

setwd("/home/yinshumin/pro_cfmf/")

##############################main.5.A##########################################
load("./output/section5/pbmc3k:AUC.ls.RData")

AUC.max.ls <- lapply(AUC.ls, function(x) {
  rownames(x) <- NULL
  max.va      <- apply(x, 1, max)
})

AUC.max.df      <- do.call(cbind, AUC.max.ls) %>% as.data.frame()

save("AUC.max.df", file = "./paper/result_section5/source_data/pbmc3k.auc.max.df.RData")

plots <- foreach(i = c(30,50,100,200,300,500,1000,2000)) %do% {
  AUC.df      <- head(AUC.max.df, i) |> rownames_to_column("id")
  AUC.df.long <- pivot_longer(AUC.df, cols = -id, names_to = "method", values_to = "AUC") |> mutate(method = gsub(".genes", "", method))
  
  
  m <- AUC.df.long |> group_by(method) |> summarise(x = median(AUC, na.rm = TRUE)) |> arrange(desc(x)) |> pull(method) %>% .[2]
  p <- ggplot(AUC.df.long, aes(x = method, y = AUC, fill = method)) +
    geom_boxplot(outlier.shape = NA) +
    labs(title = paste0("Top ",i), x = NULL, y = "AUC", fill = "Method") +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(comparisons = list(c("CFMF", m)), paired = TRUE, group = "id",
                               label = "p.signif", tip.length = 0) +
    scale_y_continuous(limits = c(0.45, 1.05), breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 1))
  p
  
  return(p)
}

pdf("./paper/result_section5/plot/main.5.A.pdf", height = 7, width = 14)
plots %>% wrap_plots(ncol = 4, guides = "collect") & 
  theme(legend.position = "bottom") &
  guides(fill = guide_legend(nrow = 1))
dev.off()



##############################main.5.C##########################################
load("./output/section5/seu.ls_of_7_methods_with_silhouette.RData")

sil.df <- map2_dfr(seu.feature.gene.ls, names(seu.feature.gene.ls), ~{
  .x@meta.data %>%
    dplyr::select(cell_type, silhouette) %>%
    mutate(method = .y)
})


sil.df$method <- gsub(".genes","",sil.df$method)

save("sil.df", file = "./paper/result_section5/source_data/pbmc3k.silhouette.df.RData")


pdf("./paper/result_section5/plot/main.5.C.pdf", height = 8, width = 10)
ggplot(sil.df, aes(x = method, y = silhouette, fill = method)) +
  geom_boxplot(outlier.shape = NA, color = "black") +
  facet_wrap(~cell_type, scales = "free_y") +
  theme_bw(base_size = 14) +
  labs(x = "method", y = "Silhouette") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


##############################main.5.D##########################################
load("./output/section5/accuracy_of_50pbmc_compared_7methods.RData")

accuracy.df <- do.call(rbind, accuracy.ls)
df          <- accuracy.df |> rownames_to_column(var = "id") |> pivot_longer(-id, names_to = "Method")
ls          <- lapply(c("SE"), c, "CFMF")

save(list = c("df","ls"), file = "./paper/result_section5/source_data/singler:plot.df.RData")

pdf("./paper/result_section5/plot/main.5.D.pdf", height = 5, width = 4)
df |>
  ggplot(aes(x = Method, y = value)) + 
  geom_point(aes(color = Method), shape = 21, fill = "white", size = 1.5, alpha = 1,      
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) +
  geom_boxplot(aes(fill = Method), width = 0.6, alpha = 0.2, linewidth = 0.3, outlier.shape = NA, notch = TRUE) +
  theme_bw() +
  ggpubr::stat_compare_means(comparisons = ls, paired = TRUE, tip.length = 0, step.increase = 0.05,  group = "id") +
  labs(x = NULL, y = "Accuracy") +
  theme(legend.position = "none",  axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()



##############################supp.6.C##########################################
hotspot.results <- read_csv("./output/section5/hotspot_all_genes_results.csv")
load("./output/section5/pbmc3k: gene of 7 mothed.RData")

Get_Method_Specificity <- function(hs.df, gene.ls) {
  
  res.ls <- list()
  
  for (method in names(gene.ls)) {
    genes.select <- gene.ls[[method]]
    
    # Filter for genes present in both the list and the Hotspot results
    valid.data <- hs.df[hs.df$Gene %in% genes.select, ]
    
    if (nrow(valid.data) > 0) {
      res.ls[[method]] <- data.frame(
        Gene = valid.data$Gene,
        Method = method,
        Z.score = valid.data$Z 
      )
    }
  }
  
  final.df <- do.call(rbind, res.ls)
  rownames(final.df) <- NULL
  return(final.df)
}

plot.df <- Get_Method_Specificity(hs.df = hotspot.results, gene.ls = gene.methods)

save("plot.df", file = "./paper/result_section5/source_data/hotspot.df.RData")

plot.df$Method <- gsub(".genes","",plot.df$Method)

ls             <- list(c("CFMF","SE"))

pdf("./paper/result_section5/plot/supp.6.C.pdf", height = 5, width = 5)
plot.df |>
  ggplot(aes(x = Method, y = Z.score, fill = Method)) +
  geom_boxplot(outlier.shape = NA, outlier.alpha = 0.5, width = 0.6) +
  scale_y_log10(breaks = c(0.1, 1, 10, 100), labels = c("0.1", "1", "10", "100")) +
  theme_bw() +
  ggpubr::stat_compare_means(comparisons = ls, paired = TRUE, tip.length = 0, step.increase = 0.05,  group = "id") +
  labs(x = NULL, y = "log10 (Z score)") +
  theme(legend.position = "none",  axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()





##############################main.5.B##########################################
load("./output/section5/Darmanis2017_log2TPM_AUC_of_7methods.RData")
## GBM: BT_S1 ------------------------------------------------------------------
sam          <- "BT_S1"
df           <- auc.sample.ls[[sam]] %>% as.data.frame()
colnames(df) <- names(auc.sample.ls[[sam]])
colnames(df) <- gsub(".genes", "", colnames(df))


save("df", file = "./paper/result_section5/source_data/main.5.B.a_sample_of_GBM.RData")

topN         <- c(30, 50, 100, 200, 300, 500, 1000, 2000)

auc.long     <- data.frame()
for (n in topN) {
  tmp <- df[1:n, ]  |> rownames_to_column("id")
  tmp$TopN <- n           
  auc_long <- pivot_longer(tmp, cols = -c(id, TopN), names_to = "Method", values_to = "AUC")
  auc.long <- bind_rows(auc.long, auc_long)
}


ls        <- split(auc.long, auc.long$TopN)
names(ls) <- paste("Top", names(ls))
plots <- foreach(i = names(ls)) %do% {
  df <- ls[[i]]
  m  <- df|> group_by(Method) |> summarise(x = median(AUC, na.rm = TRUE)) |> arrange(desc(x)) |> pull(Method) %>% .[2]
  
  ggplot(df, aes(x = Method, y = AUC, fill = Method)) +
    geom_boxplot(outlier.shape  = NA) +
    labs(title = i, x = NULL, y = "AUC") +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(comparisons = list(c("CFMF", m)), paired = TRUE, group = "id",
                               label = "p.signif", tip.length = 0) +
    scale_y_continuous(limits = c(0.45, 1.05), breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 1))
}


pdf("./paper/result_section5/plot/main.5.B.pdf", height = 7, width = 14)
plots %>% wrap_plots(ncol = 4, guides = "collect") & 
  theme(legend.position = "bottom") &
  guides(fill = guide_legend(nrow = 1))
dev.off()


##############################supp.5.A##########################################
load("./output/section5/Abdelfattah2022_GBM_AUC_of_7methods.RData")

save("auc.sample.ls", file = "./paper/result_section5/source_data/supp.5.A.RData")

plots <- lapply(names(auc.sample.ls), function(s) {
  
  
  df <- auc.sample.ls[[s]] %>% as.data.frame()
  colnames(df) <- names(auc.sample.ls[[s]])
  colnames(df) <- gsub(".genes", "", colnames(df))
  df <- df[,-8]
  
  topN          <- c(30, 50, 100, 200, 300, 500, 1000, 2000)
  auc.long      <- data.frame()
  
  for (n in topN) {
    tmp        <- df[1:n, ]
    tmp$TopN   <- n
    auc_long   <- pivot_longer(tmp, cols = -TopN, names_to = "Method", values_to = "AUC")
    auc.long   <- bind_rows(auc.long, auc_long)
  }
  
  p <- ggplot(auc.long, aes(x = factor(TopN), y = AUC, fill = Method)) +
    geom_boxplot(linewidth = 0.1, outlier.shape = NA, width = 0.7) +
    labs(x = "Number of genes", y = "AUC score", title = paste0("sample: ", s)) +
    theme_bw() +
    theme(panel.grid = element_blank())
  
  return(p)
})


pdf("./paper/result_section5/plot/supp.5.A.pdf",
    height = 12, width = 20)
plots %>%
  wrap_plots(ncol = 2, guides = "collect") &
  theme(legend.position = "right") &
  plot_annotation(
    title = "Dataset: Abdelfattah2022",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))
  )
dev.off()





##############################supp.5.B##########################################
load("./output/section5/Wang2019_AUC_of_7methods.RData")

save("auc.sample.ls", file = "./paper/result_section5/source_data/supp.5.B.RData")

plots <- lapply(names(auc.sample.ls), function(s) {
  
  
  df <- auc.sample.ls[[s]] %>% as.data.frame()
  colnames(df) <- names(auc.sample.ls[[s]])
  colnames(df) <- gsub(".genes", "", colnames(df))
  df <- df[,-8]
  
  topN          <- c(30, 50, 100, 200, 300, 500, 1000, 2000)
  auc.long      <- data.frame()
  
  for (n in topN) {
    tmp        <- df[1:n, ]
    tmp$TopN   <- n
    auc_long   <- pivot_longer(tmp, cols = -TopN, names_to = "Method", values_to = "AUC")
    auc.long   <- bind_rows(auc.long, auc_long)
  }
  
  p <- ggplot(auc.long, aes(x = factor(TopN), y = AUC, fill = Method)) +
    geom_boxplot(linewidth = 0.1, outlier.shape = NA, width = 0.7) +
    labs(x = "Number of genes", y = "AUC score", title = paste0("sample: ", s)) +
    theme_bw() +
    theme(panel.grid = element_blank())
  
  return(p)
})


pdf("./paper/result_section5/plot/supp.5.B.pdf",
    height = 12, width = 20)
plots %>%
  wrap_plots(ncol = 2, guides = "collect") &
  theme(legend.position = "right") &
  plot_annotation(
    title = "Dataset: Wang2019",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))
  )
dev.off()




##############################supp.6.A##########################################
load("./output/section5/Neftel2019_10x_AUC_of_7methods.RData")

save("auc.sample.ls", file = "./paper/result_section5/source_data/supp.6.A.RData")

plots <- lapply(names(auc.sample.ls), function(s) {
  
  
  df <- auc.sample.ls[[s]] %>% as.data.frame()
  colnames(df) <- names(auc.sample.ls[[s]])
  colnames(df) <- gsub(".genes", "", colnames(df))
  df <- df[,-8]
  
  topN          <- c(30, 50, 100, 200, 300, 500, 1000, 2000)
  auc.long      <- data.frame()
  
  for (n in topN) {
    tmp        <- df[1:n, ]
    tmp$TopN   <- n
    auc_long   <- pivot_longer(tmp, cols = -TopN, names_to = "Method", values_to = "AUC")
    auc.long   <- bind_rows(auc.long, auc_long)
  }
  
  p <- ggplot(auc.long, aes(x = factor(TopN), y = AUC, fill = Method)) +
    geom_boxplot(linewidth = 0.1, outlier.shape = NA, width = 0.7) +
    labs(x = "Number of genes", y = "AUC score", title = paste0("sample: ", s)) +
    theme_bw() +
    theme(panel.grid = element_blank())
  
  return(p)
})


pdf("./paper/result_section5/plot/supp.6.A.pdf",
    height = 12, width = 20)
plots %>%
  wrap_plots(ncol = 2, guides = "collect") &
  theme(legend.position = "right") &
  plot_annotation(
    title = "Dataset: Neftel2019",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))
  )
dev.off()


##############################supp.6.B##########################################
load("./output/section5/Darmanis2017_log2TPM_AUC_of_7methods.RData")

auc.sample.ls <- auc.sample.ls[c(1,3,4)]
save("auc.sample.ls", file = "./paper/result_section5/source_data/supp.6.B.RData")

plots <- lapply(names(auc.sample.ls), function(s) {
  
  
  df <- auc.sample.ls[[s]] %>% as.data.frame()
  colnames(df) <- names(auc.sample.ls[[s]])
  colnames(df) <- gsub(".genes", "", colnames(df))
  df <- df[,-8]
  
  topN          <- c(30, 50, 100, 200, 300, 500, 1000, 2000)
  auc.long      <- data.frame()
  
  for (n in topN) {
    tmp        <- df[1:n, ]
    tmp$TopN   <- n
    auc_long   <- pivot_longer(tmp, cols = -TopN, names_to = "Method", values_to = "AUC")
    auc.long   <- bind_rows(auc.long, auc_long)
  }
  
  p <- ggplot(auc.long, aes(x = factor(TopN), y = AUC, fill = Method)) +
    geom_boxplot(linewidth = 0.1, outlier.shape = NA, width = 0.7) +
    labs(x = "Number of genes", y = "AUC score", title = paste0("sample: ", s)) +
    theme_bw() +
    theme(panel.grid = element_blank())
  
  return(p)
})


pdf("./paper/result_section5/plot/supp.6.B.pdf",
    height = 9, width = 10)
plots %>%
  wrap_plots(ncol = 1, guides = "collect") &
  theme(legend.position = "right") &
  plot_annotation(
    title = "Dataset: Darmanis2017",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))
  )
dev.off()
