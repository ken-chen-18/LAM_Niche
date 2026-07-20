library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(ComplexHeatmap)
library(tidyr)
library(EnhancedVolcano)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

setwd("/global/scratch/users/kenchen/LAM/figures")
input_path <- "/global/scratch/users/kenchen/LAM/input/"
cellref_genes <- readRDS("/global/scratch/users/kenchen/LAM/input/files/CellRef.genes 1.rds")
ctrl_clean <- readRDS("/global/scratch/users/kenchen/LAM/3_lam_ctrl/ctrl_clean.rds")
lam_ctrl_final <- readRDS("../3_lam_ctrl/lam_ctrl_final.rds")
lam_lung_filt <- readRDS("../1_lam_int/lam_lung_filt.rds")
lam_lung_nonim <- readRDS("../1_lam_int/lam_lung_nonim.rds")
mesen_filt <- readRDS("../1_lam_int/mesen_filt.rds")
load("/global/scratch/users/kenchen/LAM/2_uterus/mesen_ut_lamcore.RData")
load("../4_deg/gse.RData")
samp_meta <- read.table(paste0(input_path, "files/sample_meta.csv"),sep=",",
                        header=T)
at2_obj <- readRDS("../5_regulon/at2_obj.rds")


colors <-  c("AF1"="#FFD700", "AF2"="#68228B", "LAMCORE-1"="#E31A1C","LAMCORE-2"="#1F78B4", "Pericyte"="#771122", "SMC"="#777711", "Mesothelial"="#BEAED4", "MyoFB"="#117744",
             "AEC"="#AA4488", "CAP1"="#6A3D9A", "CAP2"="#B15928", "LEC"="#771155", "SVEC"="#AA7744", "VEC"="#60CC52",
             "AT1"="#386CB0", "AT1/AT2"="#FBB4AE", "AT2"="#66A61E", "Basal"="#FF7F00", "Ciliated"="#FFFF99", "Goblet"="#AA4455", "RAS"="#CAB2D6", "Secretory"="#DDDD77","PNEC"="#F4CAE4",
             "Uterine LAMCORE" = "#E5D8BD", "Low-Quality" = "#000000", "LAMCORE-3"= "#CCCCCC",
             "AM" = "#E7298A", "B"="#FFD92F", "cDC2"="#FB9A99", "Plasma"= "#00CED1", "iMON" = "#FB8072", "CD8 T" = "#666666",
             "NK"="#377EB8", "CD4 T"="#00AAAD", "Treg"= "#A6D854", "Mast/Basophil"="#44AA77", "pMON"="#A65628", "cDC1"="#80B1D3",
             "NKT"="#774411","IM"="#FDBF6F", "Neutrophil"="#CCEBC5", "maDC"="#8DA0CB", "pDC"="#E6AB02", "Rare" = "#9900FF")





#CCI
lam_lung_nonim$LAM1 <- ifelse(lam_lung_nonim$LAM != "Other", lam_lung_nonim$LAM,
                              ifelse(lam_lung_nonim$celltype_011625 == "AF1", "AF1", 
                                     ifelse(lam_lung_nonim$lineage_level1=="Mesenchymal", 
                                            "Other Mesenchymal", "Other")))
lam_lung_nonim$AF <- ifelse(lam_lung_nonim$celltype_011625 == "AF1", "AF1", 
                            ifelse(lam_lung_nonim$lineage_level1=="Mesenchymal", 
                                   "Other Mesenchymal", "Other"))


ctrl_imm <- subset(ctrl_clean, lineage_level1 == "Immune")
lam_imm <- subset(lam_lung_filt, lineage_level1 == "Immune")

ctrl_imm@assays$RNA@counts@Dimnames[[1]] <- rownames(cellref_genes)
ctrl_imm@assays$RNA@data@Dimnames[[1]] <- rownames(cellref_genes)
ctrl_imm@assays$RNA@meta.features <- cellref_genes

imm <- merge(lam_imm, ctrl_imm)
imm@assays[c("SCT", "SoupX", "SoupX_SCT")] <- NULL
imm$celltype_combined <- paste0(imm$Diagnosis, "_", imm$celltype_011625)

lam2_ligs <- c("COL1A1",
               "COL1A2",
               "COL4A1",
               "COL6A1",
               "COL6A2",
               "COL6A3",
               "FN1",
               "TNC",
               "TNXB",
               "LAMA2",
               "LAMB1",
               "TNC",
               "THBS2")

af1_rec <- c("ITGA1", "ITGB1", "CD44", "SDC4", "CD47")

genes <- cellref_genes[match(af1_rec, cellref_genes$symbol),]

plots <- VlnPlot(
  lam_lung_nonim,
  features = genes$id,
  group.by = "AF",
  pt.size = 0,
  combine = FALSE
)
for (i in seq_along(plots)) {
  plots[[i]] <- plots[[i]] + ggtitle(genes[i, "symbol"]) + NoLegend() + xlab(NULL) + ylab(NULL)
}
wrap_plots(plots, ncol=5)



DotPlot(lam_lung_nonim,
        features = rownames(genes),
        group.by="AF",
        dot.scale = 3, dot.min=0.01) +
  scale_x_discrete(labels=genes$symbol) +
  scale_color_gradient2(low = "blue", mid="gray", high ="red") +
  theme(axis.text = element_text(size = 10)) +
  scale_y_discrete(limits=rev)+
  scale_radius(breaks=c(25,50,75,100), limits=c(0,100), range=c(0,6))+
  theme(axis.title = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5),)+
  labs(title=NULL) + coord_flip()


lam_ctrl_final$celltype_combined <- factor(lam_ctrl_final$celltype_combined, levels=c(
  "LAM_LAMCORE-1", "LAM_LAMCORE-2", "LAM_LAMCORE-3",
  "LAM_AEC", "Control_AEC",
  "LAM_AF1", "Control_AF1",
  "LAM_AF2", "Control_AF2",
  "LAM_AT1", "Control_AT1",
  "LAM_AT2", "Control_AT2",
  "LAM_Basal", "Control_Basal",
  "LAM_CAP1", "Control_CAP1",
  "LAM_CAP2", "Control_CAP2",
  "LAM_Ciliated", "Control_Ciliated",
  "LAM_Goblet", "Control_Goblet",
  "LAM_LEC", "Control_LEC",
  "LAM_Mesothelial", "Control_Mesothelial",
  "Control_MyoFB",
  "LAM_Pericyte", "Control_Pericyte",
  "Control_PNEC",
  "Control_Rare",
  "LAM_RAS", "Control_RAS",
  "LAM_Secretory", "Control_Secretory",
  "LAM_SMC", "Control_SMC",
  "LAM_SVEC", "Control_SVEC",
  "LAM_VEC", "Control_VEC"
))

lam_ctrl_final$Diagnosis <- factor(lam_ctrl_final$Diagnosis, levels=c("LAM", "Control"))
#lam_ctrl_final$lineage_combined <- paste0(lam_ctrl_final$Diagnosis, "_", lam_ctrl_final$lineage_level1)
lam_ctrl_final$lamcore <- ifelse(lam_ctrl_final$celltype_011625 %in% c("LAMCORE-1", "LAMCORE-2", "LAMCORE-3"), lam_ctrl_final$celltype_011625,
                                 ifelse(lam_ctrl_final$Diagnosis == "LAM", "Other LAM", "Control"))
lam_ctrl_final$lamcore <- factor(lam_ctrl_final$lamcore, levels=c("LAMCORE-1", "LAMCORE-2", "LAMCORE-3", "Other LAM", "Control"))



imm$celltype_combined <- factor(imm$celltype_combined, levels=c(
  c(rbind(paste0("LAM_", unique(imm$celltype_011625)), paste0("Control_", unique(imm$celltype_011625))))
))
#
#
#
lig <- subset(cellref_genes, symbol %in% c("ITGB2"))
VlnPlot(subset(lam_ctrl_final, lineage_level1 == "Epithelial"), features=rownames(lig), group.by="celltype_combined", pt.size=0) + ggtitle(lig$symbol) +
  scale_fill_manual(values = c(
    "LAM_AEC" = "#AA4488",
    "Control_AEC" = "pink",
    "LAM_CAP1" = "#68228B",
    "Control_CAP1" = "#CEA2FD",
    "LAM_CAP2" = "#4B3621",
    "Control_CAP2" = "#C3B091",
    "LAM_LEC" = "darkorange",
    "Control_LEC" = "#FFBD31",
    "LAM_SVEC" = "darkblue",
    "Control_SVEC" = "lightblue",
    "LAM_VEC" = "darkgreen",
    "Control_VEC" = "lightgreen"
    
  ))


VlnPlot(subset(imm, celltype_011625 %in% c("cDC1", "cDC2","AM", "iMON", "pMON")), features=rownames(lig), group.by="celltype_combined", pt.size=0) + ggtitle(lig$symbol) + 
  NoLegend() + xlab(NULL) +
  scale_fill_manual(values = c(
    "LAM_AEC" = "#AA4488",
    "Control_AEC" = "pink",
    "LAM_CAP1" = "#68228B",
    "Control_CAP1" = "#CEA2FD",
    "LAM_CAP2" = "#4B3621",
    "Control_CAP2" = "#C3B091",
    "LAM_LEC" = "darkorange",
    "Control_LEC" = "#FFBD31",
    "LAM_SVEC" = "darkblue",
    "Control_SVEC" = "lightblue",
    "LAM_VEC" = "darkgreen",
    "Control_VEC" = "lightgreen",
    "LAM_AF1" = "#B8860B",
    "Control_AF1"="#FFD700",
    "LAM_AF2"= "#68228B",
    "Control_AF2"="#CEA2FD",
    "LAM_AT2" = "darkgreen",
    "Control_AT2"="lightgreen",
    "LAM_AT1"= "darkblue",
    "Control_AT1"="lightblue",
    "Control" = "#CCEBC5",
    
    "LAM_AM" = "#AA4488",
    "Control_AM" = "pink",
    "LAM_IM" = "#68228B",
    "Control_IM" = "#CEA2FD",
    "LAM_iMON" = "#4B3621",
    "Control_iMON" = "#C3B091",
    "LAM_pMON" = "darkorange",
    "Control_pMON" = "#FFBD31",
    "LAM_cDC1"= "darkblue",
    "Control_cDC1"="lightblue",
    "LAM_cDC2" = "darkgreen",
    "Control_cDC2"="lightgreen"
    
  ))

# "Other LAM"="#A65628", 
# "Control" = "#CCEBC5"





VlnPlot(lam_ctrl_final, features=rownames(lig), group.by="lamcore", pt.size=0) + ggtitle(lig$symbol) +
  scale_fill_manual(values = c(
    "LAMCORE-1" = "#E31A1C","LAMCORE-2"="#1F78B4", "LAMCORE-3"= "#CCCCCC", "Other LAM"="#A65628", "Control" = "#CCEBC5"
  )) + NoLegend() + xlab(NULL)



colors <-  c("AF1"="#FFD700", "AF2"="#68228B", "LAMCORE-1"="#E31A1C","LAMCORE-2"="#1F78B4", "Pericyte"="#771122", "SMC"="#777711", "Mesothelial"="#BEAED4", "MyoFB"="#117744",
             "AEC"="#AA4488", "CAP1"="#6A3D9A", "CAP2"="#B15928", "LEC"="#771155", "SVEC"="#AA7744", "VEC"="#60CC52",
             "AT1"="#386CB0", "AT1/AT2"="#FBB4AE", "AT2"="#66A61E", "Basal"="#FF7F00", "Ciliated"="#FFFF99", "Goblet"="#AA4455", "RAS"="#CAB2D6", "Secretory"="#DDDD77","PNEC"="#F4CAE4",
             "Uterine LAMCORE" = "#E5D8BD", "Low-Quality" = "#000000", "LAMCORE-3"= "#CCCCCC",
             "AM" = "#E7298A", "B"="#FFD92F", "cDC2"="#FB9A99", "Plasma"= "#00CED1", "iMON" = "#FB8072", "CD8 T" = "#666666",
             "NK"="#377EB8", "CD4 T"="#00AAAD", "Treg"= "#A6D854", "Mast/Basophil"="#44AA77", "pMON"="#A65628", "cDC1"="#80B1D3",
             "NKT"="#774411","IM"="#FDBF6F", "Neutrophil"="#CCEBC5", "maDC"="#8DA0CB", "pDC"="#E6AB02", "Rare" = "#9900FF")




VlnPlot(subset(lam_ctrl_final, lineage_level1 == "Epithelial"), features=rownames(lig), group.by="Diagnosis", pt.size=0) + ggtitle(lig$symbol) +
  scale_x_discrete(labels = c(
    "LAM" = "LAM Epithelial",
    "Control" = "Control Epithelial"
  ))
VlnPlot(subset(lam_ctrl_final, celltype_011625 %in% c("AF1", "AF2")), features=rownames(lig), group.by="celltype_combined", pt.size=0) + ggtitle(lig$symbol)

VlnPlot(lam_ctrl_final, features=rownames(lig), group.by="lamcore", pt.size=0) + ggtitle(lig$symbol)


rec <- subset(cellref_genes, symbol %in% c("TGFBR2"))
VlnPlot(subset(lam_ctrl_final, lineage_level1 == "Endothelial"), features=rownames(rec), group.by="Diagnosis", pt.size=0) + ggtitle(rec$symbol) +
  scale_x_discrete(labels = c(
    "LAM" = "LAM Endothelial",
    "Control" = "Control Endothelial"
  ))
VlnPlot(subset(lam_ctrl_final, celltype_011625 %in% c("VEC", "SVEC","LEC")), features=rownames(rec), group.by="celltype_combined", pt.size=0) + ggtitle(rec$symbol)

VlnPlot(lam_ctrl_final, features=rownames(rec), group.by="lamcore", pt.size=0) + ggtitle(lig$symbol)

DotPlot(lam_ctrl_final,
        features = rownames(lig),
        group.by="a",
        dot.scale = 3, dot.min=0.01) +
  scale_x_discrete(labels=lig$symbol) +
  scale_color_gradient2(low = "blue", mid="gray", high ="red") +
  theme(axis.text = element_text(size = 10)) +
  scale_y_discrete(limits=rev)+
  scale_radius(breaks=c(25,50,75,100), limits=c(0,100), range=c(0,10))+
  theme(axis.title = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5),)+
  labs(title=NULL) +
  coord_flip()




