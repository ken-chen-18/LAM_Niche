library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(ComplexHeatmap)
library(tidyr)
library(patchwork)
library(monocle3)
library(RColorBrewer)
library(viridis)
library(reshape2)
library(cowplot)


setwd("/global/scratch/users/kenchen/LAM/figures")
input_path <- "/global/scratch/users/kenchen/LAM/input/"
cellref_genes <- readRDS("/global/scratch/users/kenchen/LAM/input/files/CellRef.genes 1.rds")
af1_obj <- readRDS("../5_regulon/af1_obj.rds")
lam_ctrl <- readRDS("../5_regulon/lam_ctrl.rds")


#Fig 3I
lam <- subset(lam_ctrl, Diagnosis == "LAM" & celltype_011625 != "LAMCORE-3")
lam$lamcore_laf <- ifelse(lam$celltype_011625 %in% c("LAMCORE-1","LAMCORE-2"), "LAMCORE",
                          ifelse(lam$celltype_011625 %in% c("AF1", "AF2"), "LAF", lam$celltype_011625))
lam$lamcore_laf <- factor(lam$lamcore_laf, levels=c("LAMCORE", "LAF", setdiff(unique(lam$lamcore_laf), c("LAMCORE", "LAF"))))
sym <- rev(c("MEOX2", "MAFB", "MFAP5", "IGFBP7", "VIM", "CXCR4", "CXCL12", "COL5A1", "COL6A1", "DDR2", "PDGFRB", "SPARC",
             "SPARCL1", "S100A4", "FAP", "OSR1", "MITF", "TCF21", "ACTA2", "PMEL", "HOXD11", "HOXA11", "HOXA10", "EMX2"))
genes <- cellref_genes[cellref_genes$symbol %in% sym, ]
genes <- genes[match(sym, genes$symbol), ]
Idents(lam) <- "lamcore_laf"
DotPlot(lam, 
                  assay = "RNA", 
                  features=rownames(genes),
                  group.by = "lamcore_laf",
                  col.min = 0,
                  col.max = 1,
                  scale.min = 0,
                  scale.max = 100,
                  scale.by = "size",
                  #idents = c("LAMCORE", #LAM
                  #           "LAF"),
                  cluster.idents=FALSE) +
  scale_x_discrete(labels=genes$symbol) +
  scale_size_area(max_size = 6)+ 
  cowplot::theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10, family="TT Times New Roman"),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size = 10, family="TT Times New Roman"),
        legend.text = element_text(size=8),
        legend.title = element_text(size = 9)) +
  scale_color_gradientn(colours = colorRampPalette(c("gray", "blue"))(100), limits = c(0,1), oob = scales::squish, name = 'log2 (count + 1)') +
  coord_flip()



#Fig 7K
sym <- c("PMEL", "ACTA2", "DES", "UNC5D", "PCP4", "TAGLN", "HOXA11","DDR2", "MLANA","FAP")
genes <- cellref_genes[match(sym, cellref_genes$symbol), ]


ct <- sort(unique(lam_ctrl$celltype_combined))

DotPlot(lam_ctrl,
        features = rownames(genes),
        group.by = "celltype_combined",
        dot.scale = 3, dot.min = 0.01) +
  scale_x_discrete(labels = genes$symbol) +
  #scale_y_discrete(limits=rev) +
  #scale_color_gradient2(low = "blue", mid="gray", high ="red") +
  scale_radius(breaks = c(25, 50, 75, 100), limits = c(0, 100), range = c(0, 7)) +
  guides(size = guide_legend(title = "Percent\nExpressed"),
         color = guide_colorbar(title = "Average\nExpression")) +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.title = element_text(size = 13, family="Arial"))






