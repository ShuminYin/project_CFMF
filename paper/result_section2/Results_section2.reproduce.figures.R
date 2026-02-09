
library(Seurat)
library(tidyverse)
library(foreach)
library(patchwork)
library(scales)

setwd("/home/yinshumin/pro_cfmf/")

############main.2.A#####################
load("./paper/result_section2/source_data/seu_of_pbmc3k.RData")

df    <- silhouette |> mutate(y = -log10(FDR))
genes <- c("S100A8", "CD3E", "CD14", "CCR7", "LEF1", "PRSS23", "AKR1C3", "FOLR3", "RBP7", "HMOX1")

pdf("./paper/result_section1/plots/main.2.A.pdf", width = 6, height = 6)
ggplot(df, aes(x = silhouette, y = y)) +
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
dev.off()


############main.2.B#####################
load("./paper/result_section2/source_data/seu_of_pbmc3k.RData")

Idents(seu) <- "cell_type"
p1 <- scCustomize::DimPlot_scCustom(seu, label = TRUE, repel = TRUE, label.size = 5) + NoLegend() + scale_color_brewer(palette = "Set2")
p2 <- scCustomize::FeaturePlot_scCustom(seu, "AKR1C3") + theme(plot.title = element_text(face = "italic")) 
p3 <- scCustomize::FeaturePlot_scCustom(seu, "FOLR3") + theme(plot.title = element_text(face = "italic")) 
p4 <- scCustomize::FeaturePlot_scCustom(seu, "HMOX1") + theme(plot.title = element_text(face = "italic")) 
p.list <- list(p1,p2,p3,p4)

pdf("./paper/result_section1/plots/main.2.B.pdf", height = 12, width = 12)
p.list %>% wrap_plots(ncol = 2)
dev.off()


#######mian.2.C#######
load("./paper/result_section1/source_data/sil.df_of_GBM_sample.RData")
genes <- c("CST3", "S100B", "BCAN", "GPR17", "LGALS3", "MLC1", "TNR", "DCX", "TNFRSF1A", "COL20A1", "ELAVL4")

pdf("./paper/result_section1/plots/mian.2.C.pdf", width = 6, height = 6)
ggplot(df, aes(x = silhouette, y = y)) +
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
    size = 3,
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

dev.off()



#######mian.2.D#######
load("./paper/result_section1/source_data/harmony_of_GBM.RData")
x      <- c("MES", "AC", "OPC", "NPC")
umaps1 <- foreach(i = x) %do% {scCustomize::FeaturePlot_scCustom(seu, features = paste0("MP_", i)) + labs(title = i)} |> wrap_plots(ncol = 1)
x      <- c("LGALS3", "MLC1", "TNR", "DCX")
umaps2 <- foreach(i = x) %do% {scCustomize::FeaturePlot_scCustom(seu, features = i) + theme(plot.title = element_text(face = "italic"))} |> wrap_plots(ncol = 1)


pdf("./paper/result_section1/plots/mian.2.D.pdf", width = 6, height = 9.5)
umaps1 | umaps2
dev.off()





############main.2.E#####################
load("./paper/result_section1/source_data/fisher_test_of_pbmc.RData")

p.value        <- fisher.result$p.value
odds.ratio     <- fisher.result$estimate
p_label        <- paste0("p = ", "1.65e-65")
or_label       <- paste0("OR = ", signif(odds.ratio, 3))
label_all      <- paste0(or_label, ", ", p_label)

pdf("./paper/result_section1/plots/main.2.E.pdf", height = 6, width = 6)
df |> mutate(Group = ifelse(Group == "In DEG", "DEGs", "Not DEGs")) |>
  ggplot(aes(x = Category, y = Count, fill = Group)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7, alpha = 0.5) +
  scale_fill_manual(values = c("DEGs" = "#D53E4F", "Not DEGs" = "#4575B4"), name = "") +
  geom_text(aes(label = label, y = prop), 
            position = position_fill(vjust = 0.5),
            color = "black", 
            size = 4) +
  annotate("text", x = 1.5, y = 1.1, label = label_all, size = 4, color = "black") +
  labs(x = "", y = "Proportion", title = "PBMC3K") +
  scale_y_continuous(breaks = seq(0, 1, 0.2), labels = scales::percent_format(), expand =  expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Proportion of DEGs", title = NULL) + 
  theme_classic(base_size = 14)
dev.off()




############main.2.F#####################
load("./paper/result_section1/source_data/fisher_test_of_GBM.RData")
p.value        <- fisher.result$p.value
odds.ratio     <- fisher.result$estimate
p_label  <- paste0("p = ", signif(p.value, 3))
or_label <- paste0("OR = ", signif(odds.ratio, 3))
label_all <- paste0(or_label, ", ", p_label)

pdf("./paper/result_section1/plots/mian.2.F.pdf", height = 6, width = 5)
df |> mutate(Group = ifelse(Group == "In GBM signature genesets", "GBM genes", "Not GBM genes")) |>
  ggplot(aes(x = Category, y = Count, fill = Group)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7, alpha = 0.5) +
  scale_fill_manual(values = c("GBM genes" = "#D53E4F", "Not GBM genes" = "#4575B4"), name = "") +
  geom_text(aes(label = label, y = prop), 
            position = position_fill(vjust = 0.5),
            color = "black", 
            size = 4) +
  annotate("text", x = 1.5, y = 1.1, label = label_all, size = 4, color = "black") +
  labs(x = "", y = "Proportion", title = "PBMC3K") +
  scale_y_continuous(breaks = seq(0, 1, 0.2), labels = scales::percent_format(), expand =  expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Proportion of GMB signature genes", title = NULL) + 
  theme_classic(base_size = 14)
dev.off()




############supp.1.A#####################
load("./paper/result_section2/source_data/seu_of_pbmc3k.RData")

Idents(seu) <- "cell_type"
p1 <- scCustomize::DimPlot_scCustom(seu, label = TRUE, repel = TRUE, label.size = 5) + NoLegend() + scale_color_brewer(palette = "Set2")
p5 <- scCustomize::FeaturePlot_scCustom(seu, "S100A8") + theme(plot.title = element_text(face = "italic")) 
p6 <- scCustomize::FeaturePlot_scCustom(seu, "CD3E") + theme(plot.title = element_text(face = "italic")) 
p7 <- scCustomize::FeaturePlot_scCustom(seu, "CD14") + theme(plot.title = element_text(face = "italic")) 
p8 <- scCustomize::FeaturePlot_scCustom(seu, "CCR7") + theme(plot.title = element_text(face = "italic")) 
p9 <- scCustomize::FeaturePlot_scCustom(seu, "LEF1") + theme(plot.title = element_text(face = "italic")) 
p.list <- list(p1,p5,p6,p7,p8,p9)

pdf("./paper/result_section1/plots/supp.1.A.pdf", height = 12, width = 18)
p.list %>% wrap_plots(ncol = 3)
dev.off()

############supp.1.B#####################
load("./paper/result_section2/source_data/seu_of_pbmc3k.RData")
Idents(seu) <- "cell_type"
p10 <- scCustomize::FeaturePlot_scCustom(seu, "TRAT1") + theme(plot.title = element_text(face = "italic"))
p11 <- scCustomize::FeaturePlot_scCustom(seu, "RBP7") + theme(plot.title = element_text(face = "italic")) 
p12 <- scCustomize::FeaturePlot_scCustom(seu, "VMO1") + theme(plot.title = element_text(face = "italic")) 
p13 <- scCustomize::FeaturePlot_scCustom(seu, "LGALS2") + theme(plot.title = element_text(face = "italic")) 
p14 <- scCustomize::FeaturePlot_scCustom(seu, "FCER2") + theme(plot.title = element_text(face = "italic")) 
p15 <- scCustomize::FeaturePlot_scCustom(seu, "PRSS23") + theme(plot.title = element_text(face = "italic")) 
p.list <- list(p10,p11,p12,p13,p14,p15)

pdf("./paper/result_section1/plots/supp.1.B.pdf", height = 12, width = 18)
p.list %>% wrap_plots(ncol = 3)
dev.off()


#######supp.1.C#######
load("./paper/result_section1/source_data/harmony_of_GBM.RData")
x      <- c("MES", "AC", "OPC", "NPC")
y      <- c("TNFRSF1A", "AQP4", "COL20A1", "ELAVL4")

plots <- foreach(i = x, j = y) %do% {
  p1 <- scCustomize::FeaturePlot_scCustom(seu, features = paste0("MP_", i)) + labs(title = i)
  p2 <- scCustomize::FeaturePlot_scCustom(seu, features = j) + theme(plot.title = element_text(face = "italic"))
  p <- p1|p2
  return(p)
}

pdf("./paper/result_section1/plots/supp.1.C.pdf", width = 12, height = 5)
plots %>% wrap_plots(ncol = 2)
dev.off()


