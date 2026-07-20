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
library(SLICE)
library(monocle3)
library(viridis)

setwd("/global/scratch/users/kenchen/LAM/figures")
input_path <- "/global/scratch/users/kenchen/LAM/input/"
cellref_genes <- readRDS("/global/scratch/users/kenchen/LAM/input/files/CellRef.genes 1.rds")
lam_ctrl_final <- readRDS("../5_regulon/lam_ctrl.rds")
lam_lung_filt <- readRDS("../1_lam_int/lam_lung_filt.rds")
lam_lung_nonim <- readRDS("../1_lam_int/lam_lung_nonim.rds")
mesen_filt <- readRDS("../1_lam_int/mesen_filt.rds")
samp_meta <- read.table(paste0(input_path, "files/sample_meta.csv"),sep=",",
                        header=T)

load("../2_uterus/mesen_ut_lamcore.RData")
load("../2_uterus/lamcore_only.RData")
load("../2_uterus/cds.RData")
load("/global/scratch/users/kenchen/LAM/2_uterus/pseudotime.RData")


af1_obj <- readRDS("../5_regulon/af1_obj.rds")
at2_obj <- readRDS("../5_regulon/at2_obj.rds")

ctrl_clean <- readRDS("/global/scratch/users/kenchen/LAM/3_lam_ctrl/ctrl_clean.rds")
#######################################

colors <-  c("AF1"="#FFD700", "AF2"="#68228B", "LAMCORE-1"="#E31A1C","LAMCORE-2"="#1F78B4", "Pericyte"="#771122", "SMC"="#777711", "Mesothelial"="#BEAED4", "MyoFB"="#117744",
             "AEC"="#AA4488", "CAP1"="#6A3D9A", "CAP2"="#B15928", "LEC"="#771155", "SVEC"="#AA7744", "VEC"="#60CC52",
             "AT1"="#386CB0", "AT1/AT2"="#FBB4AE", "AT2"="#66A61E", "Basal"="#FF7F00", "Ciliated"="#FFFF99", "Goblet"="#AA4455", "RAS"="#CAB2D6", "Secretory"="#DDDD77","PNEC"="#F4CAE4",
             "Uterine LAMCORE" = "#E5D8BD", "Low-Quality" = "#000000", "LAMCORE-3"= "#CCCCCC",
             "AM" = "#E7298A", "B"="#FFD92F", "cDC2"="#FB9A99", "Plasma"= "#00CED1", "iMON" = "#FB8072", "CD8 T" = "#666666",
             "NK"="#377EB8", "CD4 T"="#00AAAD", "Treg"= "#A6D854", "Mast/Basophil"="#44AA77", "pMON"="#A65628", "cDC1"="#80B1D3",
             "NKT"="#774411","IM"="#FDBF6F", "Neutrophil"="#CCEBC5", "maDC"="#8DA0CB", "pDC"="#E6AB02", "Rare" = "#9900FF")


lam_from_int <- subset(lam_ctrl_final, Diagnosis=="LAM")
ctrl_from_int <- subset(lam_ctrl_final, Diagnosis=="Control")



######################
#Fig S1

lam_lung_filt$DataID2 <- factor(lam_lung_filt$DataID2, levels=samp_meta$DataID2[1:13])
p1 <- DimPlot(lam_lung_filt, group.by="DataID2") + ggtitle("Sample ID")
p2 <- DimPlot(lam_lung_filt, group.by="Technology")
p3 <- DimPlot(lam_lung_filt, group.by="seurat_clusters") + ggtitle("Clusters")
p4 <- DimPlot(lam_lung_filt, group.by="celltype_011625") + scale_color_manual(values=colors) + ggtitle("Cell Type")

(p1 | p2)
(p3 | p4)


ct_dataid <- data.frame(table(lam_lung_filt[[c("celltype_011625", "DataID2")]]))
p5 <- ggplot(ct_dataid, aes(x=Freq, y=celltype_011625, fill=DataID2)) + 
  geom_bar(stat="identity", position="fill") + theme_classic() + 
  labs(x=NULL, y=NULL) + theme(legend.title = element_blank())


ct_tech <- data.frame(table(lam_lung_filt[[c("celltype_011625", "Technology")]]))
p6 <- ggplot(ct_tech, aes(x=Freq, y=celltype_011625, fill=Technology)) + 
  geom_bar(stat="identity", position="fill") + theme_classic() + 
  labs(x=NULL, y=NULL) + theme(legend.title = element_blank())

p5+p6


####################################################################
#Fig S2

p7 <- DimPlot(lam_ctrl_final, group.by="DataID2")
p7

p8 <- DimPlot(lam_ctrl_final, group.by="Technology")
p8

################################
#Fig S3

all_genes <- read.table("txt/dp_genes.txt")$V1
all_ens <- subset(cellref_genes, symbol %in% all_genes)
all_ens <- all_ens[rev(match(all_genes, all_ens$symbol)),]

cluster_order <- read.table("txt/cluster_order.txt")
lam_lung_filt$seurat_clusters1 <- factor(lam_lung_filt$seurat_clusters,
                                         levels=rev(cluster_order$V1))

p9 <- DotPlot(lam_lung_filt,
               features = rownames(all_ens),
               group.by="seurat_clusters1",
               dot.scale = 3, dot.min=0.01) +
  scale_x_discrete(labels=all_ens$symbol) +
  scale_color_gradient2(low = "blue", mid="gray", high ="red") +
  theme(axis.text = element_text(size = 10)) +
  scale_y_discrete(limits=rev)+
  scale_radius(breaks=c(25,50,75,100), limits=c(0,100), range=c(0,4))+
  theme(axis.title = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5),)+
  coord_flip()+
  ggtitle("")


p9

################################
#Fig S4
ct_order <- read.table("txt/ct_order.txt")
lam_lung_filt$celltype <- factor(lam_lung_filt$celltype_011625,
                                 levels=rev(ct_order$V1))
p10 <- DotPlot(lam_lung_filt,
               features = rownames(all_ens),
               group.by="celltype",
               dot.scale = 3, dot.min=0.01) +
  scale_x_discrete(labels=all_ens$symbol) +
  scale_color_gradient2(low = "blue", mid="gray", high ="red") +
  theme(axis.text = element_text(size = 10)) + 
  scale_y_discrete(limits=rev)+
  scale_radius(breaks=c(25,50,75,100), limits=c(0,100), range=c(0,4))+
  theme(axis.title = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5),)+
  coord_flip()+
  ggtitle("")

p10



##################################
#Fig S5

p11 <- VlnPlot(lam_ctrl_final, group.by="celltype_011625", features="nFeature_RNA", pt.size=0) +
  geom_hline(yintercept=c(500,7500), color="red") + ylim(0, 8000) + xlab(NULL) +
  theme(legend.position = "none")

p12 <- VlnPlot(lam_ctrl_final, group.by="celltype_011625", features="nCount_RNA", pt.size=0) +
  geom_hline(yintercept=c(40000), color="red") + ylim(0, 40000) + xlab(NULL)

p13 <- VlnPlot(lam_ctrl_final, group.by="celltype_011625", features="pMT", pt.size=0)+
  geom_hline(yintercept=c(20), color="red") + ylim(0, 20) + xlab(NULL) +
  theme(legend.position = "none")

p11+p12+p13


##################################
#Fig S6A-C
lam_pct_df <- read.table("../4_deg/lam_pct_df.csv", sep=",", header=T,
                         row.names = 1)
lam_lung_filt$lam_lineage <- ifelse(lam_lung_filt$celltype_011625 %in% c("LAMCORE-1", "LAMCORE-2", "LAMCORE-3"),
                                    lam_lung_filt$celltype_011625, lam_lung_filt$lineage_level1)
lam_lung_filt$lam_lineage <- factor(lam_lung_filt$lam_lineage, levels=rev(c("LAMCORE-1", "LAMCORE-2", "LAMCORE-3",
                                                                            "Mesenchymal", "Endothelial", "Epithelial",
                                                                            "Immune")))
l1 <- subset(lam_pct_df, lam1.pct>20 & l1_v_lc>3 & l1_v_nl>4)
l1 <- l1[order(l1$lam1.pct, decreasing = T),]
p14 <- DotPlot(lam_lung_filt,
               features = rownames(l1),
               group.by="lam_lineage",
               dot.scale = 3, dot.min=0.01) +
  scale_x_discrete(labels=l1$symbol) +
  scale_color_gradient2(low = "blue", mid="gray", high ="red") +
  theme(axis.text = element_text(size = 10)) +
  scale_y_discrete(limits=rev)+
  scale_radius(breaks=c(25,50,75,100), limits=c(0,100), range=c(0,10))+
  theme(axis.title = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5, face="italic"),
        axis.text=element_text(size=12))+
  labs(title=NULL)


l2 <- subset(lam_pct_df, lam2.pct>10 & l2_v_lc>1.7 & l2_v_nl>1.7)
l2 <- l2[order(l2$lam2.pct, decreasing = T),]
l2 <- subset(l2, !(symbol %in% c("CD68", "THEMIS2", "EVI2A", "C1orf162",
                                 "MNDA", "EPSTI1", "DDX60L", "SLC31A2",
                                 'FCGR3A', "SAMD3", "AIF1", "CAMK4")))
p15 <- DotPlot(lam_lung_filt,
               features = rownames(l2),
               group.by="lam_lineage",
               dot.scale = 3, dot.min=0.01) +
  scale_x_discrete(labels=l2$symbol) +
  scale_color_gradient2(low = "blue", mid="gray", high ="red") +
  theme(axis.text = element_text(size = 10)) +
  scale_y_discrete(limits=rev)+
  scale_radius(breaks=c(25,50,75,100), limits=c(0,100), range=c(0,10))+
  theme(axis.title = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5, face="italic"),
        axis.text=element_text(size=12))+
  labs(title=NULL)


l3 <- subset(lam_pct_df, lam3.pct>10 & l3_v_nl>1.4)
l3 <- l3[order(l3$lam3.pct, decreasing = T),]
p16 <- DotPlot(lam_lung_filt,
               features = rownames(l3),
               group.by="lam_lineage",
               dot.scale = 3, dot.min=0.01) +
  scale_x_discrete(labels=l3$symbol) +
  scale_color_gradient2(low = "blue", mid="gray", high ="red") +
  theme(axis.text = element_text(size = 10)) +
  scale_y_discrete(limits=rev)+
  scale_radius(breaks=c(25,50,75,100), limits=c(0,100), range=c(0,10))+
  theme(axis.title = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5, face="italic"),
        axis.text=element_text(size=12))+
  labs(title=NULL)


p14/p15/p16





##################################
#Fig S7
sel_fun <- read.table("txt/lam3vsnonlam_fun.txt", sep="\t", header=T)
sel_fun$name_short <- c("cytoplasmic translation", "muscle contraction", 
                        "muscle structure development", "oxidative phosphorylation",
                        "cell migration", "fibroblast migration", "collagen fibril organization", "G1/S transition",
                        "ACTB interactions", "FN1 interactions", "TGFB1 interactions", "LAMTOR1 interactions",
                        "SRF", "PBX1", "MEF2", "Ovary SMC (Fan)", "Lung SMC (Travaglini)",
                        "Ovary Stromal (Jones)", "TGFB EMT up (Foroutan)")
sel_fun$name_short <- factor(sel_fun$name_short, levels=sel_fun$name_short)
sel_fun[c(2:6, 8, 10:12)] <- NULL
names(sel_fun) <- c("category", "FDR", "hc", "name")

sel_fun$logp <- -log10(sel_fun$FDR)

#Fig S7B
ggplot(sel_fun, aes(x = logp, y = name, width=0.7)) +
  geom_bar(stat = "identity", fill="darkgray") +
  theme_classic() +
  labs(x = "-Log10(FDR)", y = NULL, fill = NULL) +
  geom_vline(xintercept = -log10(0.05), linetype = "dotted", color = "black") +
  theme(axis.text = element_text(size = 10, color = "black"))


sel_fun <- read.table("txt/lam3vslam1_fun.txt", sep="\t", header=T)
sel_fun$name_short <- c("Ovarian cancer (Liu)", "Uterine leiomyoma (GWAS)", 
                        "Cell adhesion", "Muscle structure development",
                        "Cell migration", "Neuron differentiation", "ECM organization",
                        "Blood vessel development", "Muscle system process", "Muscle contraction",
                        "Oxidative Phosphorylation", "Translation", "Ribosome biogenesis", 
                        "SRF", "ECM organization (Reactome)",
                        "TGFB signaling (Reactome)", "HOXA2 target genes", "Translation (Reactome)")
sel_fun$name_short <- factor(sel_fun$name_short, levels=sel_fun$name_short)
sel_fun[c(2:5, 8:12)] <- NULL
names(sel_fun) <- c("category", "LAM1_FDR", "LAM3_FDR","name")

sel_fun$lam1_logp <- -log10(sel_fun$LAM1_FDR)
sel_fun$lam3_logp <- -log10(sel_fun$LAM3_FDR)

sel_long <- sel_fun %>%
  pivot_longer(cols = c(lam1_logp, lam3_logp), names_to = "lam_type", values_to = "value")

#Fig S7C
ggplot(sel_long, aes(x = value, y = name, fill = lam_type, width=0.7)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("lam1_logp" = "#E31A1C", "lam3_logp" = "darkgray")) + 
  theme_classic() +
  labs(x = "-Log10(FDR)", y = NULL, fill = NULL) +
  geom_vline(xintercept = 1, linetype = "dotted", color = "black") +
  theme(axis.text = element_text(size = 10, color = "black"))

sel_fun <- read.table("txt/lam3vslam2_fun.txt", sep="\t", header=T)
sel_fun$name_short <- c("Lung SMC (Travaglini)", "Ovary SMC (Fan)", "Ovarian cancer (Liu)", 
                        "Lung adventital fibroblast (Travaglini)",
                        "Cell adhesion", "Cell migration","ECM organization",
                        "Blood vessel development", "Neuron differentiation", 
                        "Collagen fibril organization", "MAPK cascade", "Wound healing",
                        "Muscle structure development", "Oxidative Phosphorylation",
                        "Muscle contraction",
                        "Translation", "Ribosome biogenesis", 
                        "Ribosome (KEGG)", "Translation (Reactome)", "SRF", "LEF1",
                        "ECM organization (Reactome)",
                        "RTK signaling (Reactome)")
sel_fun$name_short <- factor(sel_fun$name_short, levels=sel_fun$name_short)
sel_fun[c(2:5, 8:12)] <- NULL
names(sel_fun) <- c("category", "LAM2_FDR", "LAM3_FDR","name")

sel_fun$lam2_logp <- -log10(sel_fun$LAM2_FDR)
sel_fun$lam3_logp <- -log10(sel_fun$LAM3_FDR)

sel_long <- sel_fun %>%
  pivot_longer(cols = c(lam2_logp, lam3_logp), names_to = "lam_type", values_to = "value")

#Fig S7D
ggplot(sel_long, aes(x = value, y = name, fill = lam_type, width=0.7)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("lam2_logp" = "#1F78B4", "lam3_logp" = "darkgray")) + 
  theme_classic() +
  labs(x = "-Log10(FDR)", y = NULL, fill = NULL) +
  geom_vline(xintercept = 1, linetype = "dotted", color = "black") +
  theme(axis.text = element_text(size = 10, color = "black"))


lam3_reg <- read.table("txt/lam3_reg.txt")


top <- lam3_reg[order(lam3_reg$rra, decreasing = FALSE), ][1:15, ]
top$TF <- factor(rownames(top), levels=rev(rownames(top)))
top$mean_delta <-(top$dc_nonlam_delta+top$sc_nonlam_delta)/2
top$mean_p <- (-log10(top$dc_nonlam_p_val)+-log10(top$sc_nonlam_p_val))/2
top$pct.fc <- ifelse(top$pct.fc>10, 10, top$pct.fc)

#Fig S7E
ggplot(top, aes(x = mean_p, y = TF, size = pct.fc, color = mean_delta)) +
  geom_segment(aes(x = mean_p, xend = 0, y = TF, yend = TF), color = "gray", size = 0.5) +  
  geom_point() +  
  geom_vline(xintercept=-log10(0.05), linetype="dashed") +
  scale_color_gradient(low = "lightgray", high = "blue") +  
  scale_size_continuous(range = c(1, 10), limits=c(0.1,10)) + 
  labs(x = "Mean -Log(p-value)", y = "Transcription Factor", title = "Regulon Activity",
       color = " Mean Cliff's\nDelta", size = "Expression\nFrequency Ratio") +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black", size = 12),
    axis.title=element_text(size=14))


#treatment
treat <- read.table("txt/treatment_meta.txt", header=T, sep="\t")
df <- lam_ctrl_final[[c("celltype_011625", "DataID")]]

prop_df <- df %>%
  dplyr::count(DataID, celltype_011625, name = "Freq") %>%
  group_by(DataID) %>%
  mutate(Freq = Freq / sum(Freq)) %>%
  ungroup()

prop_df <- subset(prop_df, celltype_011625 %in% c("LAMCORE-1", "LAMCORE-2", "LAMCORE-3"))
prop_df <- prop_df %>%
  complete(
    DataID,
    celltype_011625,
    fill = list(Freq = 0)
  )
prop_df$treated <- treat[match(prop_df$DataID, treat$DataID), "Treated"]
prop_df$treated <- ifelse(prop_df$treated =="Yes", "Treated", "Untreated")
prop_df$age <- treat[match(prop_df$DataID, treat$DataID), "Age"]

#Fig S7F
ggplot(prop_df, aes(x = celltype_011625, y = Freq, fill = treated)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    outlier.shape = NA
  ) +
  theme_classic() +
  ylab("Proportion of cells") +
  xlab(NULL) +
  theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5, color="black"),
        axis.text.y=element_text(color="black"))




#S7GH
g0_genes <- c("CFLAR", "CALCOCO1", "YPEL3", "CST3", "SERINC1", "CLIP4", "PCYOX1", 
              "TMEM59", "RGS2", "YPEL5", "CD63", "KIAA1109", "CDH13", "GSN", "MR1", 
              "CYB5R1", "AZGP1", "ZFYVE1", "DMXL1", "EPS8L2", "PTTG1IP", "MIR22HG", 
              "PSAP", "GOLGA8B", "NEAT1", "TXNIP", "MTRNR2L12")
g0_genes <- subset(cellref_genes, symbol %in% g0_genes)


lamcore <- subset(lam_lung_filt, celltype_011625 %in% c("LAMCORE-1", "LAMCORE-2", "LAMCORE-3"))
lamcore$g0 <- PercentageFeatureSet(lamcore, features=g0_genes$id)

g1s <- read.table("txt/g1_s_trans.txt")
g1s <- subset(cellref_genes, symbol %in% g1s$V1)
lamcore$g1s <- PercentageFeatureSet(lamcore, features=g1s$id)


lamcore$Treated <- treat[match(lamcore$DataID2, treat$SampleID), "Treated"]


lamcore$Treated <- ifelse(lamcore$Treated == "No", "Untreated", "Treated")
lamcore$Treated1 <-  factor(lamcore$Treated, levels=c("Treated", "Untreated"))

ggplot(lamcore[[]], aes(y=g0, x=celltype_011625)) + 
  geom_violin(aes(fill=Treated1)) +
  theme_classic() + theme(axis.text.x=element_text(color="black", size=10),
                          axis.text.y=element_text(color="black", size=10)) +
  ylab("G0 enrichment") + xlab(NULL) + NoLegend()


ggplot(lamcore[[]], aes(y=g1s, x=celltype_011625)) + 
  geom_violin(aes(fill=Treated1)) +
  theme_classic() + theme(axis.text.x=element_text(color="black", size=10),
                          axis.text.y=element_text(color="black", size=10)) +
  ylab("G1S transition enrichment") + xlab(NULL) + NoLegend() + ylim(0, 7)


##################################
#Fig S8
#Fig S8A and S9A
DimPlot(mesen_ut_lamcore, group.by="celltype_011625", reduction="lam_mesh_umap") + 
  scale_color_manual(values=colors)

p3_genes <- rev(c("PMEL", "PI15", "PBX1", "MYH9", "MYH11", "UNC5D", "ESR1", "PCP4", "PTGER3", "PGR", 
                  "HAND2-AS1", "EMX2", "RAMP1", "CTHRC1", "MMP11", "FAP", "LUM", "COL1A1"))
p3_genes <- cellref_genes[match(p3_genes, cellref_genes$symbol),]

mesen_ut_lamcore$celltype_011625 <- factor(mesen_ut_lamcore$celltype_011625, levels=rev(c('Uterine LAMCORE',"LAMCORE-1", 
                                                            "LAMCORE-2", "LAMCORE-3", "AF1", "AF2", "Mesothelial", 
                                                            "SMC", "Pericyte")))

p3 <- DotPlot(mesen_ut_lamcore,
              features = rownames(p3_genes),
              group.by="celltype_011625",
              dot.scale = 3, dot.min=0.01) +
  scale_x_discrete(labels=p3_genes$symbol) +
  scale_color_gradient2(low = "blue", mid="gray", high ="red") +
  theme(axis.text = element_text(size = 10)) +
  scale_y_discrete(limits=rev)+
  scale_radius(breaks=c(25,50,75,100), limits=c(0,100), range=c(0,10))+
  theme(axis.title = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5),
        axis.text=element_text(size=12),
        axis.text.y = element_text(face = "italic"))+
  labs(title=NULL) +
  coord_flip()

#Fig S8B
p3

#Fig S8C
slice_obj <- readRDS("/global/scratch/users/kenchen/LAM/0_cihan/slice_obj.rds")
entropy = slice_obj@data@phenoData@data

ggplot(entropy, aes(x=entropy, fill=state)) + geom_density()

#Fig S8D
DimPlot(pseudotime, group.by="celltype_011625")
FeaturePlot(pseudotime, features="monocle3_pseudotime")

#Fig S8E
ut_deg <- read_tsv("txt/ut_lamcore_deg.txt", col_names =T)
pan_ut <- read_tsv("txt/pan_ut_deg.txt", col_names =T)
ut_leio <- read_tsv("txt/ut_leio_deg.txt", col_names =T)

p1 <- FeaturePlot(mesen_ut_lamcore, features = ut_deg$id, order=F)+
  scale_color_gradient2(low = "darkblue", mid="darkgreen", high ="yellow",
                        midpoint = 0.1, limits = c(0.05, 0.15),
                        oob = scales::squish)
p2 <- FeaturePlot(subset(mesen_ut_lamcore, celltype_090523 != "Uterine LAM"), features = ut_deg$id, order=F)+
  scale_color_gradient2(low = "darkblue", mid="darkgreen", high ="yellow",
                        midpoint = 0.1, limits = c(0.05, 0.15),
                        oob = scales::squish)
p1+p2


p3 <- FeaturePlot(mesen_ut_lamcore, features = pan_ut$id, order=F)+
  scale_color_gradient2(low = "darkblue", mid="darkgreen", high ="yellow",
                        midpoint = 0.1, limits = c(0.05, 0.15),
                        oob = scales::squish)
p4 <- FeaturePlot(subset(mesen_ut_lamcore, celltype_090523 != "Uterine LAM"), features = pan_ut$id, order=F)+
  scale_color_gradient2(low = "darkblue", mid="darkgreen", high ="yellow",
                        midpoint = 0.1, limits = c(0.05, 0.15),
                        oob = scales::squish)
p3+p4

p5 <- FeaturePlot(mesen_ut_lamcore, features = ut_leio$id, order=F)+
  scale_color_gradient2(low = "darkblue", mid="darkgreen", high ="yellow",
                        midpoint = 0.1, limits = c(0.05, 0.15),
                        oob = scales::squish)
p6 <- FeaturePlot(subset(mesen_ut_lamcore, celltype_090523 != "Uterine LAM"), features = ut_leio$id, order=F)+
  scale_color_gradient2(low = "darkblue", mid="darkgreen", high ="yellow",
                        midpoint = 0.1, limits = c(0.05, 0.15),
                        oob = scales::squish)
p5+p6


#Fig S8F
DefaultAssay(mesen_ut_lamcore) <- "RNA"
mesen_ut_lamcore <- mesen_ut_lamcore %>% NormalizeData() %>% FindVariableFeatures()
mesen_ut_lamcore <- ScaleData(mesen_ut_lamcore)
a <- subset(mesen_ut_lamcore, celltype_011625 != "LAMCORE-3")
av <- AverageExpression(a, assays="RNA",slot="scale.data",
                        group.by="celltype_011625", 
                        features=VariableFeatures(mesen_ut_lamcore))$RNA
cor.exp <- cor(av)
ord <- rev(c('Uterine LAMCORE',"LAMCORE-1", "LAMCORE-2", "AF1", "AF2", "Mesothelial", "SMC", "Pericyte"))
cor.exp <- cor.exp[ord, ord]

p2 <- Heatmap(cor.exp, cluster_rows = FALSE, cluster_columns = FALSE,
              heatmap_legend_param = list(
                title = NULL,        # remove legend title
                legend_height = unit(4, "cm"),  # increase legend size vertically
                labels_gp = gpar(fontsize = 12)  # increase font size of legend labels
              ))

p2


##################################
#Fig S9
#Fig S9B
DimPlot(mesen_ut_lamcore,reduction="lam_mesh_umap", 
        cells.highlight = WhichCells(mesen_ut_lamcore, idents="Uterine LAMCORE"),
        pt.size = 0.25, sizes.highlight = 0.5,
        cols.highlight = "#FFC8AD") + NoLegend() + ggtitle("Uterine LAMCORE")

#Fig S9C
DimPlot(mesen_ut_lamcore,reduction="lam_mesh_umap", 
        cells.highlight = WhichCells(mesen_ut_lamcore, idents="LAMCORE-1"),
        pt.size = 0.25, sizes.highlight = 0.5,
        cols.highlight = "#E31A1C") + NoLegend() + ggtitle("LAMCORE-1")

#Fig S9D
DimPlot(mesen_ut_lamcore,reduction="lam_mesh_umap", 
        cells.highlight = WhichCells(mesen_ut_lamcore, idents="LAMCORE-2"),
        pt.size = 0.25, sizes.highlight = 0.5,
        cols.highlight = "#1F78B4") + NoLegend() + ggtitle("LAMCORE-2")

#Fig S9E
DimPlot(mesen_ut_lamcore,reduction="lam_mesh_umap", 
        cells.highlight = WhichCells(mesen_ut_lamcore, idents="LAMCORE-3"),
        pt.size = 0.25, sizes.highlight = 0.5,
        cols.highlight = "black") + NoLegend() + ggtitle("LAMCORE-3")


#Partition
mesen_ut_lamcore <- AddMetaData(mesen_ut_lamcore, as.data.frame(Embeddings(mesen_ut_lamcore, reduction = "lam_mesh_umap")))
meta <- mesen_ut_lamcore@meta.data
meta <- subset(meta, celltype_011625 %in% c("LAMCORE-1", "LAMCORE-2", "LAMCORE-3", "Uterine LAMCORE"))
meta$partition <- ifelse(meta$lammeshumap_1< 3, "R2", 
                         ifelse(meta$lammeshumap_1<6.8 & meta$lammeshumap_2<2.4, "R1", 
                                ifelse(meta$lammeshumap_2>2.4, "R3", "R4")))

pcolors <-  c("R1"="#E31A1C","R2"="#1F78B4", "R4"="#777711", "R3"= "#CCCCCC")

#Fig S9F
ggplot(meta, aes(x=celltype_011625, fill=partition)) + geom_bar(position="fill", width=0.5) +
  scale_fill_manual(values = pcolors) + theme_classic() + 
  theme(axis.text.y=element_text(size=15), axis.text.x=element_text(size=0)) +
  xlab(NULL) + ylab(NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))





#Fig S11

pseudotime <- NormalizeData(pseudotime)

lam_pseudotime <- subset(pseudotime, celltype_011625 %in% c("LAMCORE-1", "LAMCORE-2"))
lc_sub <- subset(lam_ctrl_final, cells = Cells(lam_pseudotime))

lam_pseudotime[["pyscenic"]] <- lc_sub[["pyscenic"]]


pseudo_exp <- GetAssayData(lam_pseudotime, slot = "data", assay="RNA")
pseudo_reg <- GetAssayData(lam_pseudotime, slot = "data", assay="pyscenic")

#Fig S11A
lam1_tfs <- c("MYOCD", "SRF", "MEF2A", "HOXA10", "PITX2", "PURA")
lam1_tf_df <- cbind(lam_pseudotime[[c("monocle3_pseudotime", "celltype_011625")]], 
                    t(pseudo_reg[lam1_tfs,]))

plots <- lapply(lam1_tfs, function(x) {
  return(ggplot(lam1_tf_df, aes_string(x="monocle3_pseudotime", y=x, col="celltype_011625")) + geom_point() +
           geom_smooth(method = "loess", color = "black", se=F) + theme_classic()+ ylab(x) +
           xlab("Pseudotime") +  NoLegend())
})
wrap_plots(plots)


lam1_tgs <- c("MYLK", "CNN1", "TPM2", "MYL9", "DES", "TPM1", "RAMP1", "HAND2-AS1", 
              "HOXA10", "LEFTY2", "APBB1", "MBP")
lam1_tgs <- subset(cellref_genes, symbol %in% lam1_tgs)

lam1_tg_df <- cbind(lam_pseudotime[[c("monocle3_pseudotime", "celltype_011625")]], 
                    t(pseudo_exp[rownames(lam1_tgs),]))

plots <- apply(lam1_tgs, 1, function(x) {
  return(ggplot(lam1_tg_df, aes_string(x="monocle3_pseudotime", y=x["id"], col="celltype_011625")) + geom_point() +
           geom_smooth(method = "loess", color = "black", se=F) + theme_classic()+ ylab(x["symbol"]) +
           xlab("Pseudotime") +  NoLegend())
})
wrap_plots(plots)

#Fig S11B
lam2_tfs <- c("TWIST1", "TWIST2", "SNAI2", "RUNX2", "ZEB2", "PRRX2")
lam2_tf_df <- cbind(lam_pseudotime[[c("monocle3_pseudotime", "celltype_011625")]], t(pseudo_reg[lam2_tfs,]))

plots <- lapply(lam2_tfs, function(x) {
  return(ggplot(lam2_tf_df, aes_string(x="monocle3_pseudotime", y=x, col="celltype_011625")) + geom_point() +
           geom_smooth(method = "loess", color = "black", se=F) + theme_classic()+ ylab(x) +
           xlab("Pseudotime") +  NoLegend())
})
wrap_plots(plots)

lam2_tgs <- c("FN1", "DCN", "MMP9", "MMP2", "VIM", "MMP13", "FAP", "VCAN", "COL1A1", 
              "POSTN", "CXCL12", "CDH11")
lam2_tgs <- subset(cellref_genes, symbol %in% lam2_tgs)

lam2_tg_df <- cbind(lam_pseudotime[[c("monocle3_pseudotime", "celltype_011625")]], 
                    t(pseudo_exp[rownames(lam2_tgs),]))

plots <- apply(lam2_tgs, 1, function(x) {
  return(ggplot(lam2_tg_df, aes_string(x="monocle3_pseudotime", y=x["id"], col="celltype_011625")) + geom_point() +
           geom_smooth(method = "loess", color = "black", se=F) + theme_classic()+ ylab(x["symbol"]) +
           xlab("Pseudotime") +  NoLegend())
})
wrap_plots(plots)



#################################################
#S16C-E

RenameFeatures <- function(obj, old, new) {
  obj@meta.data[new] <- obj@meta.data[old]
  obj@meta.data[old] <- NULL
  return(obj)
}

genes <- c("MMP2", "MYL9", "THBS1", "IGFBP7", "FGF7", "EGR1", "DDR2", "STAT1", "STAT3", "COL14A1", "CXCL12", "CNN3")
genes <- subset(cellref_genes, symbol %in% genes)

p1 <- VlnPlot(af1_obj, features=rownames(genes), group.by="Diagnosis", pt.size=0)
p2 <- wrap_plots(lapply(p1, function(x) {
  ens <- x$labels$title
  x$labels$title <- cellref_genes[ens, "symbol"]
  x$labels$x <- NULL
  x$labels$y <- NULL
  return(x)
}), ncol=2)

#Fig S16C
p2

af1.fun <- read.table("txt/af1.fun.full.txt", sep="\t", header=T, check.names = F)

af1.fun.list <- apply(af1.fun, 2, function(tf) {
  tf <- tf[tf != ""]
  sym <- rownames(subset(cellref_genes, symbol %in% tf))
  return(sym)
})

af1_obj <- AddModuleScore(af1_obj, features=af1.fun.list)
af1_obj <- RenameFeatures(af1_obj, paste0("Cluster", 1:6), 
                          names(af1.fun.list))


lam_af1 <- subset(af1_obj, Diagnosis=="LAM")
ctrl_af1 <- subset(af1_obj, Diagnosis=="Control")

fun_names <- names(af1.fun.list)
limits <- list(c(0.05, 0.2), c(0.05, 0.2), c(0.1, 0.25), c(0.1, 0.5), c(0.1,0.25),
               c(0,0.1))

p11 <- Map(function(fun, lim) {
  title_plot <- ggplot() +
    theme_void() +
    annotate("text", x = 0.5, y = 0.5, label = fun, size = 5, hjust = 0.5, vjust = 0.5)
  
  p1 <- FeaturePlot(lam_af1, features = fun, order = FALSE) +
    scale_color_viridis(option = "D", direction = 1,
                        limits = lim, oob = scales::squish) +
    ggtitle(NULL) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none"
    )
  
  p2 <- FeaturePlot(ctrl_af1, features = fun, order = FALSE) +
    scale_color_viridis(option = "D", direction = 1,
                        limits = lim, oob = scales::squish) +
    ggtitle(NULL) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
  
  title_plot / (p1 | p2) + plot_layout(heights = c(0.1, 1))
}, fun = fun_names, lim = limits)

p12 <- wrap_plots(p11, ncol=1)

#Fig S16D
p12


af1_reg <- read.table("txt/af1_reg_rank.txt", sep="\t", header=T, row.names=1)
af1_reg$TF <- factor(rownames(af1_reg), levels=rev(rownames(af1_reg)))
af1_reg$logScore <- -log10(af1_reg$Score)
af1_reg$logDC <- -log10(af1_reg$dc_p_val)
af1_reg$logSC <- -log10(af1_reg$sc_p_val)
af1_reg$mean_p <- rowMeans(af1_reg[c("logDC", "logSC")])
af1_reg$mean_delta <- rowMeans(af1_reg[c("dc_delta", "sc_delta")])

#Fig S16E
ggplot(af1_reg, aes(x = mean_p, y = TF, size = pct.fc, color = mean_delta)) +
  geom_segment(aes(x = mean_p, xend = 0, y = TF, yend = TF), color = "gray", size = 0.5) +  # Thinner line and reversed direction
  geom_point() +  # Circles on top of the line
  geom_vline(xintercept=-log10(0.05), linetype="dashed") +
  scale_color_gradient(low = "lightgray", high = "red") +  # Color scale from gray to red
  scale_size_continuous(range = c(3, 12), limits=c(1,4)) +  # Change the range of point sizes
  labs(x = "Mean -Log(p-value)", y = "Transcription Factor", title = "Regulon Activity",
       color = " Mean Cliff's\nDelta", size = "Expression\nFrequency Ratio") +  # Custom legend titles
  theme_classic() +
  theme(
    axis.text = element_text(color = "black", size = 12),
    axis.title=element_text(size=14), 
    legend.position = "bottom",legend.direction = "vertical")


##################
#Fig S17
#Fig S17A
DimPlot(at2_obj, group.by="Diagnosis") +
  ggtitle(NULL) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.text = element_text(size=25)
  )

at2.fun <- read.table("txt/at2.fun.full.txt", sep="\t", header=T, check.names = F)

at2.fun.list <- apply(at2.fun, 2, function(tf) {
  tf <- tf[tf != ""]
  sym <- rownames(subset(cellref_genes, symbol %in% tf))
  return(sym)
})

at2_obj <- AddModuleScore(at2_obj, features=at2.fun.list)
at2_obj <- RenameFeatures(at2_obj, paste0("Cluster", 1:length(at2.fun.list)), 
                          names(at2.fun.list))


lam_at2 <- subset(at2_obj, Diagnosis=="LAM")
ctrl_at2 <- subset(at2_obj, Diagnosis=="Control")

fun_names <- names(at2.fun.list)
limits <- list(c(0.05, 0.15), c(0, 0.1), c(0.05, 0.1), c(0,0.25), c(0.05, 0.1))

p11 <- Map(function(fun, lim) {
  title_plot <- ggplot() +
    theme_void() +
    annotate("text", x = 0.5, y = 0.5, label = fun, size = 5, hjust = 0.5, vjust = 0.5)
  
  p1 <- FeaturePlot(lam_at2, features = fun, order = F) +
    scale_color_viridis(option = "D", direction = 1,
                        limits = lim, oob = scales::squish
    ) +
    ggtitle(NULL) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
  
  p2 <- FeaturePlot(ctrl_at2, features = fun, order = F) +
    scale_color_viridis(option = "D", direction = 1,
                        limits = lim, oob = scales::squish
    ) +
    ggtitle(NULL) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none"
    )
  
  title_plot / (p2 | p1) + plot_layout(heights = c(0.1, 1))
}, fun = fun_names, lim = limits)

p12 <- wrap_plots(p11, ncol=1)

#Fig S17C
p12


at2_reg <- read.table("txt/at2_reg.txt", sep="\t", header=T, row.names=1)
at2_reg$TF <- factor(rownames(at2_reg), levels=rev(rownames(at2_reg)))
at2_reg$logScore <- -log10(at2_reg$Score)
at2_reg$logDC <- -log10(at2_reg$dc_p_val)
at2_reg$logSC <- -log10(at2_reg$sc_p_val)
at2_reg$mean_p <- rowMeans(at2_reg[c("logDC", "logSC")])
at2_reg$mean_delta <- rowMeans(at2_reg[c("dc_delta", "sc_delta")])
at2_reg$mean_p <- ifelse(at2_reg$mean_p>250, 250, at2_reg$mean_p)

#Fig S17D
ggplot(at2_reg, aes(x = mean_p, y = TF, size = pct.fc, color = mean_delta)) +
  geom_segment(aes(x = mean_p, xend = 0, y = TF, yend = TF), color = "gray", size = 0.5) +  # Thinner line and reversed direction
  geom_point() +  # Circles on top of the line
  geom_vline(xintercept=-log10(0.05), linetype="dashed") +
  scale_color_gradient(low = "lightgray", high = "red") +  # Color scale from gray to red
  scale_size_continuous(range = c(4, 12), limits=c(1,3)) +  # Change the range of point sizes
  labs(x = "Mean -Log(p-value)", y = "Transcription Factor", title = "Regulon Activity",
       color = " Mean Cliff's\nDelta", size = "Expression\nFrequency Ratio") +  # Custom legend titles
  theme_classic() +
  theme(
    axis.text = element_text(color = "black", size = 12),
    axis.title=element_text(size=14))


##################
#Fig S18
#AT2 senescence
at2_obj$Diagnosis <- factor(at2_obj$Diagnosis, levels=c("LAM", "Control"))
reactome <- read.table("txt/cellular_sen.csv", sep="\t", header=T)
reactome$id <- cellref_genes[match(reactome$Symbol, cellref_genes$symbol), "id"]
reactome <- na.omit(reactome)
reactome_fc <- FindMarkers(at2_obj, features=reactome$id, min.pct=0, logfc.threshold=0, ident.1="LAM")
reactome_fc$symbol <- cellref_genes[rownames(reactome_fc), "symbol"]

at2_obj <- AddModuleScore(at2_obj, features=list(reactome$id), name="reactome")

#Fig S18A
VlnPlot(at2_obj, features="reactome1", group.by="Diagnosis", pt.size=0) + 
  ggtitle("Cellular Senescence")

reactome_fc <- reactome_fc[order(reactome_fc$avg_log2FC),]
reactome_fc$pct.fc <- reactome_fc$pct.1/reactome_fc$pct.2
reactome_fc_filt <- subset(reactome_fc, (pct.1>0.1 | pct.2>0.1) & (pct.fc>1.2 | pct.fc<1/1.2))

at2_obj$Diagnosis <- factor(at2_obj$Diagnosis, levels=c("Control", "LAM"))

#Fig S18B
DotPlot(at2_obj,
        features = rownames(reactome_fc_filt),
        group.by="Diagnosis",
        dot.scale = 3, dot.min=0.01) +
  scale_x_discrete(labels=reactome_fc_filt$symbol) +
  scale_color_gradient2(low = "blue", mid="gray", high ="red") +
  theme(axis.text = element_text(size = 10)) +
  scale_y_discrete(limits=rev)+
  scale_radius(breaks=c(25,50,75,100), limits=c(0,100), range=c(0,10))+
  theme(axis.title = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5),)+
  labs(title=NULL) +
  coord_flip()


##################
#Fig S19C

lam_ctrl_final$LAM_ct <- ifelse(lam_ctrl_final$celltype_011625 %in% c("LAMCORE-1", "LAMCORE-2", "LAMCORE-3"), 
                               lam_ctrl_final$celltype_011625,
                              ifelse(lam_ctrl_final$Diagnosis == "LAM", "Other LAM", 
                                     "Control"))
VlnPlot(lam_ctrl_final, features="ENSG00000110492", group.by="LAM_ct", pt.size=0) + ggtitle("MDK")

af_obj <- subset(lam_ctrl_final, celltype_011625 %in% c("AF1", "AF2"))
VlnPlot(af_obj, features="ENSG00000137801", group.by="celltype_combined", pt.size=0) + ggtitle("THBS1")

at_obj <- subset(lam_ctrl_final, celltype_011625 %in% c("AT1", "AT2"))
at_genes <- subset(cellref_genes, symbol %in% c("SDC4", "APP", "ICAM1"))

at_plot <- apply(at_genes, 1, function(x) {
  return(VlnPlot(at_obj, features=x[["id"]], group.by="celltype_combined", pt.size=0) + 
           ggtitle(x[["symbol"]]))
})
wrap_plots(at_plot)


ctrl_imm <- subset(ctrl_clean, lineage_level1 == "Immune")
lam_imm <- subset(lam_lung_filt, lineage_level1 == "Immune")

ctrl_imm@assays$RNA@counts@Dimnames[[1]] <- rownames(cellref_genes)
ctrl_imm@assays$RNA@data@Dimnames[[1]] <- rownames(cellref_genes)
ctrl_imm@assays$RNA@meta.features <- cellref_genes

imm <- merge(lam_imm, ctrl_imm)
imm@assays[c("SCT", "SoupX", "SoupX_SCT")] <- NULL
imm$celltype_combined <- paste0(imm$Diagnosis, "_", imm$celltype_011625)

imac <- subset(imm, celltype_011625 == "IM")
VlnPlot(imac, features="ENSG00000019582", group.by="celltype_combined", pt.size=0) + ggtitle("CD74")

amac <- subset(imm, celltype_011625 == "AM")
VlnPlot(amac, features="ENSG00000160255", group.by="celltype_combined", pt.size=0) + ggtitle("ITGB2")

