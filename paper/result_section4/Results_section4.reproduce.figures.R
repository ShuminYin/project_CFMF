library(Seurat)
library(harmony)
library(tidyverse)
library(foreach)
library(ggplot2)
library(patchwork)
library(magrittr)
library(ComplexHeatmap)


setwd("/home/yinshumin/pro_cfmf/")

###############################mian.4.B########################################
load("./paper/result_section4/source_data/mian.4.B.RData")
p2 <- DimPlot(seu, raster = T, group.by = "sample") +
  labs(title = "CFMF genes") +
  tidydr::theme_dr() +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5), 
        line = element_line(linewidth = 0.1),
        panel.grid = element_blank())

p1 <- DimPlot(seu.malignant, raster = T, group.by = "sample") +
  labs(title = "HVGs") +
  tidydr::theme_dr() +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5), 
        line = element_line(linewidth = 0.1),
        panel.grid = element_blank())

pdf("./paper/result_section4/plot/mian.4.B.pdf", height = 6, width = 12)
p1 | p2
dev.off()



###############################mian.4.C########################################
load( "./paper/result_section4/source_data/seu.high.RData")
subclone.colors <- c(
  subclone1 = "#7FC8A9",
  subclone2 = "#8BB4D4",
  subclone3 = "#D8A47F",
  subclone4 = "#B2A3D4",
  subclone5 = "#E3A1B5",
  subclone6 = "#DFC380",
  subclone7 = "#6F9EAF",
  subclone8 = "#CC8B79"
)

p2 <- DimPlot(seu.high, group.by = "subclone.final", raster = FALSE) + 
  scale_color_manual(values = subclone.colors) +
  labs(title = "Refined Subclones")

pdf("./paper/result_section4/plot/mian.4.C.pdf", height = 6, width = 7)
p2
dev.off()


#######################mian.4.D#################################################
load("./paper/result_section4/source_data/mian.4.D.RData")

subclone.colors <- c(
  subclone1 = "#7FC8A9",
  subclone2 = "#8BB4D4",
  subclone3 = "#D8A47F",
  subclone4 = "#B2A3D4",
  subclone5 = "#E3A1B5",
  subclone6 = "#DFC380",
  subclone7 = "#6F9EAF",
  subclone8 = "#CC8B79"
)


col.ha <- HeatmapAnnotation(subclone = anno_block(gp = gpar(fill = subclone.colors), 
                                                  labels = names(subclone.colors), 
                                                  labels_gp = gpar(col = "white", fontsize = 10))
)


col.fun <- circlize::colorRamp2(seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out = 100), colorRampPalette(c("skyblue", "white", "firebrick"))(100))
col.ha  <- HeatmapAnnotation(subclone = anno_block(gp = gpar(fill = subclone.colors), 
                                                   labels = names(subclone.colors), 
                                                   labels_gp = gpar(col = "white", fontsize = 10, fontface = "bold"))
)

pdf("./paper/result_section4/plot/mian.4.D.pdf", width = 10, height = 5)
Heatmap(mat, col = col.fun, name = "z-score",
        cluster_rows = TRUE, cluster_columns = FALSE,
        show_row_names = TRUE, show_column_names = FALSE, row_names_gp = gpar(fontsize = 8),
        column_split = coldata$cluster, top_annotation = col.ha, use_raster = F,
        column_title = "MSigDB Hallmark Pathways", show_row_dend = FALSE,
        border = TRUE, border_gp = gpar(col = "black", lwd = 0.01),
        row_split = 4, row_title = " ")
dev.off()

#######################mian.4.E#################################################
load("./paper/result_section4/source_data/seu.high.RData")
load("./paper/result_section4/source_data/mian.4.E.RData")

library(ComplexHeatmap)
subclone.colors <- c(
  subclone1 = "#7FC8A9",
  subclone2 = "#8BB4D4",
  subclone3 = "#D8A47F",
  subclone4 = "#B2A3D4",
  subclone5 = "#E3A1B5",
  subclone6 = "#DFC380",
  subclone7 = "#6F9EAF",
  subclone8 = "#CC8B79"
)


col.fun        <- circlize::colorRamp2(seq(0, 2, length.out = 100), colorRampPalette(c("#FFF0D4", "#FFDAB9", "#FF8C00", "#993300"))(100))
pdf("./paper/result_section4/plot/mian.4.E.pdf", width = 4, height = 6)
Heatmap(Roe, col = col.fun, cell_fun = cell_fun, name = "Ro/e",
        cluster_rows = FALSE, cluster_columns = FALSE,
        column_names_rot = 0, column_names_centered = TRUE, column_names_side = "top",
        border = TRUE, 
        column_title = "", row_names_gp = gpar(col = "black"))
dev.off()


#######################mian.4.F / supp.4.C#################################################
load("./paper/result_section4/source_data/mian.4.F.RData")

plots <- foreach(i = paste0("subclone", 1:7)) %do% {
  merged_df |> 
    filter(subclone.final == i) |>
    ggpubr::ggscatter(x = "pro", y = "CD8T.pro", color = "#2E5984", size = 1, alpha = 0.6,
                      add = "reg.line", add.params = list(color = "#C85200", fill = "gray"),
                      conf.int = TRUE, cor.coef = TRUE, cor.coeff.args = list(method = "spearman", size = 3.5)) + 
    labs(x     = paste(i, " (% malignant cells)"), 
         y     = expression(CD8^"+"~"T cells (% microenvironment)"), 
         title = i) +
    theme_classic(base_size = 12) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5)) +
    coord_fixed(ratio = 1) +
    theme(aspect.ratio = 1)
}

pdf("./paper/result_section4/plot/mian.4.F.pdf", width = 4, height = 4)
plots[c(4)] %>% wrap_plots(ncol = 1)
dev.off()

pdf("./paper/result_section4/plot/supp.4.C.pdf", width = 12, height = 8)
plots[c(1,2,3,5,6,7)] %>% wrap_plots(ncol = 3)
dev.off()



#######################mian.4.G#################################################
load("./paper/result_section4/source_data/mian.4.G.RData")

subclone.colors <- c(
  subclone1 = "#7FC8A9",
  subclone2 = "#8BB4D4",
  subclone3 = "#D8A47F",
  subclone4 = "#B2A3D4",
  subclone5 = "#E3A1B5",
  subclone6 = "#DFC380",
  subclone7 = "#6F9EAF"
)

col.fun <- circlize::colorRamp2(seq(0, 1, length.out = 100), colorRampPalette(c("#F7F7F7", "#D9D9D9", "#F7CACA", "#F1898F", "#D94A6E", "#A62F5E"))(100))
col.ha  <- HeatmapAnnotation(subclone = anno_block(gp = gpar(fill = subclone.colors), 
                                                   labels = names(subclone.colors), 
                                                   labels_gp = gpar(col = "white", fontsize = 10, fontface = "bold"))
)
p <- Heatmap(mat, col = col.fun, name = "proportion (%)",
             cluster_rows = TRUE, cluster_columns = FALSE,
             show_row_names = FALSE, show_column_names = FALSE,
             column_split = colnames(mat), top_annotation = col.ha,
             show_row_dend = FALSE, row_title = "Spots", column_title = " ",
             border = TRUE, border_gp = gpar(col = "black", lwd = 0.01),
             heatmap_legend_param = list(direction = "horizontal", title_position = "leftcenter"))


pdf("./paper/result_section4/plot/mian.4.G.pdf", width = 7, height = 5)
draw(p, heatmap_legend_side = "bottom")
dev.off()



###############################supp.3.A########################################
load("./paper/result_section4/source_data/mian.4.B.RData")

seu$log10_nCount   <- log10(seu$nCount_RNA)
seu$log10_nFeature <- log10(seu$nFeature_RNA)
clusters           <- seu$seurat_clusters
avg.nCount         <- tapply(seu$log10_nCount, clusters, mean)
avg.nFeature       <- tapply(seu$log10_nFeature, clusters, mean)
avg.nCount.df      <- data.frame(Cluster      = names(avg.nCount),
                                 Avg.nCount   = as.numeric(avg.nCount),
                                 Avg.nFeature = as.numeric(avg.nFeature))

Idents(seu) <- "seurat_clusters"
pdf("./paper/result_section4/plot/supp.3.A.pdf", height = 10, width = 14)
p1 <- FeaturePlot(seu, c("nCount_RNA", "nFeature_RNA"),ncol = 1,label = T, raster = F) & scale_color_gradient(trans = "log10", low = "lightgrey", high = "blue")
p2 <- VlnPlot(seu, features = c("nCount_RNA", "nFeature_RNA"),ncol = 1, pt.size = 0, raster = F)
p1|p2
dev.off()


#######################supp.3.B#################################################
load("./paper/result_section4/source_data/supp.3.B.RData")

subclone.colors <- c(
  subclone1 = "#7FC8A9",
  subclone2 = "#8BB4D4",
  subclone3 = "#D8A47F",
  subclone4 = "#B2A3D4",
  subclone5 = "#E3A1B5",
  subclone6 = "#DFC380",
  subclone7 = "#6F9EAF",
  subclone8 = "#CC8B79"
)


col.fun <- circlize::colorRamp2(seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out = 100), colorRampPalette(c("skyblue", "white", "firebrick"))(100))
col.ha  <- HeatmapAnnotation(subclone = anno_block(gp = gpar(fill = subclone.colors), 
                                                   labels = names(subclone.colors), 
                                                   labels_gp = gpar(col = "white", fontsize = 15, fontface = "bold"))
)
pdf("./paper/result_section4/plot/supp.3.B.pdf", width = 18, height = 10)
Heatmap(mat, col = col.fun, name = "z-score",
        cluster_rows = TRUE, cluster_columns = FALSE,
        show_row_names = TRUE, show_column_names = FALSE, row_names_gp = gpar(fontsize = 7),
        column_split = coldata$cluster, top_annotation = col.ha, use_raster = F,
        column_title = "Reactome Pathways", show_row_dend = FALSE,
        border = TRUE, border_gp = gpar(col = "black", lwd = 0.01),
        row_split = 5, row_title = " ")
dev.off()




#######################supp.4.A#################################################
load("./paper/result_section4/source_data/seu.high.RData")
Idents(seu.high) <- "seurat_clusters"
genes <- c("SLC2A1","ADH1C","CEACAM7","KIAA1324","LGR5","HMGB2")
p  <- FeaturePlot(seu.high, features = genes,label = T, raster = F, ncol = 3) 
p1 <- p & theme(plot.title = element_text(face = "italic")) 

pdf("./paper/result_section4/plot/supp.4.A.pdf",height = 8, width = 12)
p1
dev.off()


#######################supp.4.B#################################################
load("./paper/result_section4/source_data/seu.high.RData")
load("./paper/result_section4/source_data/supp.4.B.RData")

pdf("./paper/result_section4/plot/supp.4.B.pdf",height = 8, width = 10)
Idents(seu.high) <- "subclone.final"
DotPlot(seu.high, features = marker.2) +
  scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(7, "RdYlBu"))) +
  RotatedAxis() +
  labs(x = NULL, y = NULL) +
  geom_vline(xintercept = c(3.5, 6.5, 9.5, 12.5, 15.5, 18.5),
             linetype = 2, color = "grey") +
  scale_x_discrete(labels = setNames(
    lapply(marker.2, function(g) bquote(italic(.(g)))), marker.2
  ))
dev.off()


#######################supp.7.A#################################################
load("./paper/result_section4/source_data/mian.4.B.RData")
p1 <- FeaturePlot(seu, features = c("LGR5","TOP2A"), raster=FALSE) & NoLegend()
p2 <- FeaturePlot(seu, features = c("LGR5","TOP2A"), raster=FALSE) & NoLegend()


pdf("./paper/result_section4/plot/supp.7.A.pdf", height = 10, width = 10)
list(p1,p2) %>% wrap_plots(ncol = 1)
dev.off()

