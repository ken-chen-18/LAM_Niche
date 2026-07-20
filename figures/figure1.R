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
lam_ctrl_final <- readRDS("../3_lam_ctrl/lam_ctrl_final.rds")
lam_lung_filt <- readRDS("../1_lam_int/lam_lung_filt.rds")
mesen_filt <- readRDS("../1_lam_int/mesen_filt.rds")
load("../4_deg/gse.RData")
samp_meta <- read.table(paste0(input_path, "files/sample_meta.csv"),sep=",",
                        header=T)


#######################################
#Fig 1A
colors <-  c("AF1"="#FFD700", "AF2"="#68228B", "LAMCORE-1"="#E31A1C","LAMCORE-2"="#1F78B4", "Pericyte"="#771122", "SMC"="#777711", "Mesothelial"="#BEAED4", "MyoFB"="#117744",
             "AEC"="#AA4488", "CAP1"="#6A3D9A", "CAP2"="#B15928", "LEC"="#771155", "SVEC"="#AA7744", "VEC"="#60CC52",
             "AT1"="#386CB0", "AT1/AT2"="#FBB4AE", "AT2"="#66A61E", "Basal"="#FF7F00", "Ciliated"="#FFFF99", "Goblet"="#AA4455", "RAS"="#CAB2D6", "Secretory"="#DDDD77","PNEC"="#F4CAE4",
             "Uterine LAMCORE" = "#E5D8BD", "Low-Quality" = "#000000", "LAMCORE-3"= "#CCCCCC",
             "AM" = "#E7298A", "B"="#FFD92F", "cDC2"="#FB9A99", "Plasma"= "#00CED1", "iMON" = "#FB8072", "CD8 T" = "#666666",
             "NK"="#377EB8", "CD4 T"="#00AAAD", "Treg"= "#A6D854", "Mast/Basophil"="#44AA77", "pMON"="#A65628", "cDC1"="#80B1D3",
             "NKT"="#774411","IM"="#FDBF6F", "Neutrophil"="#CCEBC5", "maDC"="#8DA0CB", "pDC"="#E6AB02", "Rare" = "#9900FF")


lam_from_int <- subset(lam_ctrl_final, Diagnosis=="LAM")
ctrl_from_int <- subset(lam_ctrl_final, Diagnosis=="Control")

p1 <- DimPlot(lam_from_int, group.by="celltype_011625") + scale_color_manual(values=colors)+ NoLegend()
p2 <- DimPlot(ctrl_from_int, group.by="celltype_011625") + scale_color_manual(values=colors)+ NoLegend()

p1+p2

###########################
#Fig 1B
p3 <- DimPlot(mesen_filt, group.by = "celltype_011625", order=c("LAMCORE-1", "LAMCORE-2", "AF2"))+
  scale_colour_manual(values = colors)+
  theme(axis.text = element_blank(), axis.title = element_blank())+ggtitle("")
p3



###########################
#Fig 1C
lam_lung_filt$lam_lineage <- ifelse(lam_lung_filt$celltype_011625 %in% c("LAMCORE-1", "LAMCORE-2", "LAMCORE-3"),
                                    lam_lung_filt$celltype_011625, lam_lung_filt$lineage_level1)
lam_lung_filt$lam_lineage <- factor(lam_lung_filt$lam_lineage, levels=rev(c("LAMCORE-1", "LAMCORE-2", "LAMCORE-3",
                                                                        "Mesenchymal", "Endothelial", "Epithelial",
                                                                        "Immune")))
p4_genes <- c("CTSK","PMEL", 'PCP4',"ESR1", 'EMX2','PI15', "HOXA10",'LRRC7',
              'PLIN4', "HOXD11",'GAD1','PTCHD4','SPINK13','MLANA','TDO2','POSTN','CTHRC1','MMP11')
p4_genes <- cellref_genes[match(p4_genes, cellref_genes$symbol),]

p4 <- DotPlot(lam_lung_filt,
               features = rownames(p4_genes),
               group.by="lam_lineage",
               dot.scale = 3, dot.min=0.01) +
  scale_x_discrete(labels=p4_genes$symbol) +
  scale_color_gradient2(low = "blue", mid="gray", high ="red") +
  theme(axis.text = element_text(size = 10)) +
  scale_y_discrete(limits=rev)+
  scale_radius(breaks=c(25,50,75,100), limits=c(0,100), range=c(0,10))+
  theme(axis.title = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5),)+
  labs(title=NULL)

p4


#############################
#Fig 1D

lam_ctrl_final <- NormalizeData(lam_ctrl_final)
de_helper <- function(obj, ident.1, ident.2=NULL, only.pos=T, p=0.01, min.pct=0.1, 
                      logfc.threshold=log2(1.5), group.by="celltype_011625") {
  mark <- FindMarkers(obj, ident.1=ident.1, ident.2=ident.2, min.pct=min.pct,
                      logfc.threshold=logfc.threshold, only.pos=only.pos,
                      group.by=group.by)
  de <- subset(mark, p_val < p)
  de$symbol <- cellref_genes[rownames(de), "symbol"]
  de
}
l12 <- de_helper(lam_ctrl_final, ident.1="LAMCORE-1", ident.2="LAMCORE-2", only.pos=F, 
                 logfc.threshold=0, p=1)

keyvals <- ifelse(
  l12$avg_log2FC < -log2(1.5) & l12$p_val<0.01, "#1F78B4",
  ifelse(l12$avg_log2FC > log2(1.5) & l12$p_val<0.01, "#E31A1C",
         'grey'))
names(keyvals)[keyvals == '#E31A1C'] <- 'LAMCORE-1'
names(keyvals)[keyvals == 'grey'] <- 'None'
names(keyvals)[keyvals == '#1F78B4'] <- 'LAMCORE-2'

p5 <- EnhancedVolcano(l12, lab=l12$symbol, x="avg_log2FC", y="p_val", pCutoff = 0.01,xlim = c(-3, 3.1),
                      FCcutoff = log2(1.5), pointSize=2, labSize=0,
                      selectLab= c("CTHRC1", "FAP"),
                      drawConnectors = F, widthConnectors =0.5, typeConnectors="open",colCustom = keyvals,colAlpha=1)+
  theme(axis.text = element_blank(), axis.title = element_blank())+ggtitle("")

p5


################################
#Fig 1E

sel_fun <- read.table("txt/selected_fun.txt", sep="\t", header=T)
sel_fun$name_short <- c("EMT up (Jechlinger)", "TGFB EMT up (Foroutan)", "Uterine SMC (Zhang)",
                        'Uterine leiomyoma (GWAS)', 'Regulation of apoptotic signaling',
                        'Positive regulation of cell-cell adhesion',
                        'Canonical Wnt signaling', 'Collagen biosynthetic process', 'Collagen metabolic process',
                        'Muscle cell development', 'Muscle contraction', 'Myofibril assembly',
                        'TGFB1 interactions', 'MMP2 interactions', 'PI3K/Akt/mTOR signaling','Regulation of IGF/IGFBP signaling',
                        'ECM degradation','ECM organization', 'Vascular SMC contraction', ' Muscle contraction',
                        'Rho GTPases activate PAKs', 'V$SRF_Q4')
sel_fun[c("ID",'Source','LAM1_Hits', 'LAM2_Hits')] <- NULL

sel_fun <- sel_fun[c(3,2,1,4,11,10,12,5,7,6,8,9,13,14,20,21,19,15,16,17,18,22),]
sel_fun$lam1_logp <- -log10(sel_fun$LAM1_FDR)
sel_fun$lam2_logp <- -log10(sel_fun$LAM2_FDR)

sel_long <- sel_fun %>%
  pivot_longer(cols = c(lam1_logp, lam2_logp), names_to = "lam_type", values_to = "value")
sel_long$name_short <- factor(sel_long$name_short, levels=rev(sel_fun$name_short))

p6 <- ggplot(sel_long, aes(x = value, y = name_short, fill = lam_type, width=0.7)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("lam1_logp" = "#E31A1C", "lam2_logp" = "#1F78B4")) + 
  theme_classic() +
  labs(x = NULL, y = NULL, fill = NULL) +
  geom_vline(xintercept = 1, linetype = "dotted", color = "black") +
  theme(axis.text = element_text(size = 14, color = "black"))

p6

####################
#Fig 1F

p7 <- gseaplot2(gsego_lam, color=c("blue", "purple"),geneSetID = c("GO:0010718", #EMT
                                                                   #"GO:0006936", #Musc
                                                                   "GO:0030198" #ECM
), subplots=1) + theme_classic()

p8 <- gseaplot2(gsekegg_lam,color=c("red", "orange"), 
                geneSetID = c("hsa05213", #endometrial cancer
                              "hsa04270"#, #vsmc contraction
                              #"hsa04512" #ECM
                              
                ), subplots=1) + theme_classic()

p8/p7


