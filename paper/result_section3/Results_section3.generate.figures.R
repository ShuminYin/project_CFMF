library(Seurat)
library(foreach)
library(tidyverse)
library(ggplot2)
library(patchwork)


setwd("/home/yinshumin/pro_cfmf/")


load("./output/section3/seu and silhouette of samples.RData")

######################main.3.A#####################################
df  <- silhouette.ls[[13]] %>% filter(pct < 0.05)
seu <- seu.ls[[13]]
df  <- df |> mutate(y = -log10(FDR))

genes <- "SCG3"

p1 <- ggplot(df, aes(x = silhouette, y = y)) +
  geom_point(alpha = 0.5, size = 0.8, color = "grey") +
  geom_point(data = df |> filter(silhouette > 0, y > 200), 
             alpha = 0.7, size = 1, color = "#2166AC") +
  geom_vline(xintercept = 0, color = "grey", linetype = "dashed", linewidth = 0.2) +
  geom_hline(yintercept = 200, color = "grey", linetype = "dashed", linewidth = 0.2) +
  theme_classic() +
  labs(x = "Silhouette Coefficient",
       y = expression(-log[10](meta~italic(p)))) +
  geom_point(data = filter(df, gene %in% genes),
             color = "black", fill = "gold", size = 2, shape = 21) +
  ggrepel::geom_label_repel(
    data = df |> filter(gene %in% genes),
    aes(label = gene),
    color = "black",
    fill = "white",
    segment.color = "grey50",
    size = 5,
    box.padding = 0.3,
    min.segment.length = 0.1,
    max.overlaps = 20,
    force = 3,
    nudge_x = 0,
    fontface = "italic",
    label.size = 0.2
  ) +
  coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1, line = element_line(linewidth = 0.3))

seu$neuro <- ifelse(seu$ann_level_4 == "Neuroendocrine","Neuroendocrine","Other")
Idents(seu) <- "neuro"
p2 <- DimPlot(seu) + 
  scale_color_manual(values = c("Other" = "lightgrey", "Neuroendocrine" = "#377EB8")) +
  theme(legend.position = c(0.05, 0.12), 
        legend.box = "horizontal", 
        legend.box.background = element_rect(color = "black"),
        legend.key.size = unit(0.9, "lines"),
        legend.text = element_text(size = 8))


p3 <- scCustomize::FeaturePlot_scCustom(seu, features = "SCG3", order = TRUE) + scale_color_gradient(low = "lightgrey", high = "blue") + theme(plot.title = element_text(face = "italic"))

colors <- c(
  "#6BAED6", "#2171B5", 
  "#74C476", "#238B45", 
  "#FB6A4A", "#A50F15", 
  "#FDBF6F", "#E6550D", 
  "#9E9AC8", "#54278F", 
  "#F0E442", "#8C510A"
)

p4 <- DimPlot(seu,
              reduction = "umap",
              group.by = "seurat_clusters",
              repel = TRUE) + scale_color_manual(values = colors)


pdf("./paper/result_section3/plot/mian.3.A.pdf", height = 4, width = 16)
p1|p2|p3|p4
dev.off()

save(list = c("seu","df"), file = "./paper/result_section3/source_data/mian:seu_and_sil_of_neuro.RData")


######################main.3.B#####################################
df  <- silhouette.ls[[38]] %>% filter(pct < 0.05)
seu <- seu.ls[[38]]
df  <- df |> mutate(y = -log10(FDR))


genes <- c("ASCL3","BSND")

p1 <- ggplot(df, aes(x = silhouette, y = y)) +
  geom_point(alpha = 0.5, size = 0.8, color = "grey") +
  geom_point(data = df |> filter(silhouette > 0, y > 200), 
             alpha = 0.7, size = 0.8, color = "#2166AC") +
  geom_vline(xintercept = 0, color = "grey", linetype = "dashed", linewidth = 0.2) +
  geom_hline(yintercept = 200, color = "grey", linetype = "dashed", linewidth = 0.2) +
  theme_classic() +
  labs(x = "Silhouette Coefficient",
       y = expression(-log[10](meta~italic(p)))) +
  geom_point(data = filter(df, gene %in% genes),
             color = "black", fill = "gold", size = 2, shape = 21) +
  ggrepel::geom_label_repel(
    data = df |> filter(gene %in% genes),
    aes(label = gene),
    color = "black",
    fill = "white",
    segment.color = "grey50",
    size = 5,
    box.padding = 0.3,
    min.segment.length = 0.1,
    max.overlaps = 20,
    force = 3,
    nudge_x = 0,
    fontface = "italic",
    label.size = 0.2
  ) +
  coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1, line = element_line(linewidth = 0.3))

seu$Iono <- ifelse(seu$ann_level_4 == "Ionocyte","Ionocyte","Other")
Idents(seu) <- "Iono"
p2 <- DimPlot(seu) + 
  scale_color_manual(values = c("Other" = "lightgrey", "Ionocyte" = "#377EB8")) +
  theme(legend.position = c(0.05, 0.12), 
        legend.box = "horizontal", 
        legend.box.background = element_rect(color = "black"),
        legend.key.size = unit(0.9, "lines"),
        legend.text = element_text(size = 8))

p3 <- scCustomize::FeaturePlot_scCustom(seu, features = "BSND", order = TRUE) + scale_color_gradient(low = "lightgrey", high = "blue") + theme(plot.title = element_text(face = "italic"))

colors <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(20)

p4 <- DimPlot(seu,
              reduction = "umap",
              group.by = "seurat_clusters",
              repel = TRUE) + scale_color_manual(values = colors)


pdf("./paper/result_section3/plot/mian.3.B.pdf", height = 6, width = 16)
p1|p2|p3|p4
dev.off()

save(list = c("seu","df"), file = "./paper/result_section3/source_data/mian:seu_and_sil_of_iono.RData")


##################supp######################

sample.sup <- c('D326_Biop_Pro1', 'T120',  'T85',"T154")
seu.supp.ls <- seu.ls[sample.sup]

plots <- lapply(sample.sup, function(sample) {
  seu <- seu.supp.ls[[sample]]
  Idents(seu) <- "rare"
  p1  <- DimPlot(seu) + 
    scale_color_manual(values = c("other" = "lightgrey", "Neuroendocrine" = "#E41A1C", "Ionocyte" = "#377EB8"))
  
  Idents(seu) <- "seurat_clusters"
  p2  <- DimPlot(seu)
  
  p3  <- FeaturePlot(seu, "SCG3") + NoLegend()
  p4  <- FeaturePlot(seu, "ASCL3") + NoLegend()
  
  p   <- p1 | p2 | p3 | p4
  return(p)
})

pdf("./paper/result_section3/plot/sup.A-D.pdf",height = 8, width = 20)
plots[1:2] %>% wrap_plots(ncol = 1)
plots[3:4] %>% wrap_plots(ncol = 1)
dev.off()


save(list = c("seu.supp.ls"), file = "./paper/result_section3/source_data/supp:seu.RData")




