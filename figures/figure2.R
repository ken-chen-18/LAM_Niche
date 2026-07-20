library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(ComplexHeatmap)
library(tidyr)
library(patchwork)
library(monocle3)
library(colorRamp2)
library(RColorBrewer)


setwd("/global/scratch/users/kenchen/LAM/figures")
input_path <- "/global/scratch/users/kenchen/LAM/input/"
cellref_genes <- readRDS("/global/scratch/users/kenchen/LAM/input/files/CellRef.genes 1.rds")
lam_ctrl <- readRDS("../5_regulon/lam_ctrl.rds")
load("../2_uterus/lam_ut.RData")
load("/global/scratch/users/kenchen/LAM/2_uterus/pseudotime.RData")
load("/global/scratch/users/kenchen/LAM/2_uterus/cds.RData")


#######################################
#Fig 2A
lam <- subset(lam_ctrl, Diagnosis=="LAM")
DefaultAssay(lam) <- "pyscenic"
lam <- ScaleData(lam)
av.sc <- AverageExpression(lam, assays = 'pyscenic', slot="scale.data",
                           group.by='celltype_011625')$pyscenic
cor.sc <- cor(av.sc)

ord <- c("Goblet", "Secretory", "Basal", "RAS", "AT2", "AT1", "Ciliated", "CAP1", 
         "CAP2", "AEC", "VEC", "SVEC", "LEC", "Mesothelial", "AF2", "AF1", "Pericyte", 
         "SMC", "LAMCORE-2", "LAMCORE-1")
cor.sc <- cor.sc[ord,ord]

p1 <- Heatmap(cor.sc,cluster_rows=F, cluster_columns = F,
              heatmap_legend_param = list(
                title = "Regulon activity\ncorrelation", 
                legend_height = unit(4, "cm"),
                title_gp = gpar(fontsize = 14),
                labels_gp = gpar(fontsize = 12))) 
p1

#######################################
#Fig 2B
cds_l <- cds[,colData(cds) %>%
               subset(
                 !(celltype_011625 %in% c("Uterine LAMCORE", "LAMCORE-3"))
               ) %>%
               row.names
]

sym <- c("MYOCD", "TCF21", "SRF", "HOXA10", "MEF2A", "PITX2", "SMAD2", "MITF","SMAD4", 
         "KLF12", "FGF2","PRRX1", "TWIST2", "RUNX2", "TWIST1", "SKIL", "SNAI2", "ZEB2", "PRRX2")#,"CDH1", "CDH2") #EMT

genes <- subset(cellref_genes, symbol %in% sym)
genes <- genes[match(sym,genes$symbol),]

exp_mat_t <- exprs(cds_l)[match(rownames(genes),rownames(rowData(cds_l))),order(pseudotime(cds_l))]

exp_mat <- t(apply(exp_mat_t,1,function(x){smooth.spline(x,df=3)$y}))
exp_mat <- t(apply(exp_mat,1,function(x){(x-mean(x))/sd(x)}))
rownames(exp_mat) <- genes$symbol
colnames(exp_mat) <- colnames(exp_mat_t)

lc <- subset(lam_ctrl, cells= Cells(cds_l))
DefaultAssay(lc) <- "pyscenic"

reg_mat <- lc@assays$pyscenic@data

reg_mat <- reg_mat[,order(pseudotime(cds_l))]
reg_mat_t <- reg_mat[sym,]
reg_mat <- t(apply(reg_mat_t,1,function(x){smooth.spline(x,df=3)$y}))
reg_mat <- t(apply(reg_mat,1,function(x){(x-mean(x))/sd(x)}))
colnames(reg_mat) <- colnames(reg_mat_t)

lc <- colData(cds_l)[order(pseudotime(cds_l)), "celltype_011625"]
lc_col <- HeatmapAnnotation(celltype=lc, col=list(celltype=c("LAMCORE-1"= "red", "LAMCORE-2"= "blue")),annotation_name_gp = gpar(col = NA))

exp_hm <- Heatmap(
  exp_mat,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 10),
  clustering_method_rows = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = F,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  bottom_annotation = lc_col)

reg_hm <- Heatmap(
  reg_mat,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 10),
  clustering_method_rows = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = F,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  bottom_annotation = lc_col)


exp_hm+reg_hm

#########
#Fig 2C

lam1_reg <- read.table("txt/lam1_reg.txt")
lam2_reg <- read.table("txt/lam2_reg.txt")



lam1 <- c("MYOCD", "TCF21", "SRF", "PURA", "HOXA10", "MEF2A", "PITX2", "SMAD2", "EGR2","FOXF1")
lam2 <- c("PRRX2","ZEB2", "SATB1", "SNAI2", "SKIL", "TWIST1", "RARG", "RUNX2", "HDGF", "NR2F1")

lam1_reg$TF <- factor(rownames(lam1_reg), levels=rev(rownames(lam1_reg)))
lam1_reg <- subset(lam1_reg, TF %in% lam1)
lam1_reg$logScore <- -log10(lam1_reg$Score)
lam1_reg$logDC <- -log10(lam1_reg$dc_lam12_p_val)
lam1_reg$logSC <- -log10(lam1_reg$sc_lam12_p_val)
lam1_reg$mean_p <- rowMeans(lam1_reg[c("logDC", "logSC")])
lam1_reg$mean_delta <- rowMeans(lam1_reg[c("dc_lam12_delta", "sc_lam12_delta")])
lam1_reg$pct.fc.12 <- ifelse(lam1_reg$pct.fc.12>4, 4, lam1_reg$pct.fc.12)

ggplot(lam1_reg, aes(x = logScore, y = TF, size = LAM1.pct, color = mean_delta)) +
  geom_segment(aes(x = logScore, xend = 0, y = TF, yend = TF), color = "gray", size = 0.5) + 
  geom_point() + 
  scale_color_gradient(low = "lightgray", high = "red") + 
  scale_size_continuous(range = c(5, 10), limits=c(0.19, 0.8)) + 
  labs(x = "-Log(RRA)", y = "Transcription Factor", title = "Regulon Activity",
       color = " Mean Cliff's\nDelta", size = "Expression\nFrequency") +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black", size = 12),
    axis.title=element_text(size=14))

lam2_reg$TF <- factor(rownames(lam2_reg), levels=rev(rownames(lam2_reg)))
lam2_reg <- subset(lam2_reg, TF %in% lam2)
lam2_reg$logScore <- -log10(lam2_reg$Score)
lam2_reg$logDC <- -log10(lam2_reg$dc_lam12_p_val)
lam2_reg$logSC <- -log10(lam2_reg$sc_lam12_p_val)
lam2_reg$mean_p <- rowMeans(lam2_reg[c("logDC", "logSC")])
lam2_reg$mean_delta <- rowMeans(lam2_reg[c("dc_lam12_delta", "sc_lam12_delta")])
lam2_reg$pct.fc.12 <- ifelse(lam2_reg$pct.fc.12>4, 4, lam2_reg$pct.fc.12)

ggplot(lam2_reg, aes(x = logScore, y = TF, size = LAM2.pct, color = mean_delta)) +
  geom_segment(aes(x = logScore, xend = 0, y = TF, yend = TF), color = "gray", size = 0.5) + 
  geom_point() + 
  scale_color_gradient(low = "lightgray", high = "blue") + 
  scale_size_continuous(range = c(3, 10), limits=c(0.1, 0.8)) + 
  labs(x = "-Log(RRA)", y = "Transcription Factor", title = "Regulon Activity",
       color = " Mean Cliff's\nDelta", size = "Expression\nFrequency") +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black", size = 12),
    axis.title=element_text(size=14))

###########
#Fig 2D
pseudotime <- NormalizeData(pseudotime)

lam_pseudotime <- subset(pseudotime, celltype_011625 %in% c("LAMCORE-1", "LAMCORE-2"))
lc_sub <- subset(lam_ctrl, cells = Cells(lam_pseudotime))

lam_pseudotime[["pyscenic"]] <- lc_sub[["pyscenic"]]


pseudo_exp <- GetAssayData(lam_pseudotime, slot = "data", assay="RNA")
pseudo_reg <- GetAssayData(lam_pseudotime, slot = "data", assay="pyscenic")

lam1_tfs <- c("SRF", "MYOCD", "HOXA10")
lam1_tf_df <- cbind(lam_pseudotime[[c("monocle3_pseudotime", "celltype_011625")]], 
                    t(pseudo_reg[lam1_tfs,]))

plots <- lapply(lam1_tfs, function(x) {
  return(ggplot(lam1_tf_df, aes_string(x="monocle3_pseudotime", y=x, col="celltype_011625")) + geom_point() +
    geom_smooth(method = "loess", color = "black", se=F) + theme_classic()+ ylab(x) +
    xlab("Pseudotime") +  NoLegend())
})
wrap_plots(plots)


lam1_tgs <- c("ACTA2", "MYH11", "ACTG2", "TAGLN", "PTGER3", "CDKN1A")
lam1_tgs <- subset(cellref_genes, symbol %in% lam1_tgs)

lam1_tg_df <- cbind(lam_pseudotime[[c("monocle3_pseudotime", "celltype_011625")]], 
                    t(pseudo_exp[rownames(lam1_tgs),]))

plots <- apply(lam1_tgs, 1, function(x) {
  return(ggplot(lam1_tg_df, aes_string(x="monocle3_pseudotime", y=x["id"], col="celltype_011625")) + geom_point() +
           geom_smooth(method = "loess", color = "black", se=F) + theme_classic()+ ylab(x["symbol"]) +
           xlab("Pseudotime") +  NoLegend())
})
wrap_plots(plots)

#Fig 2E
lam2_tfs <- c("SNAI2", "TWIST1", "PRRX2")
lam2_tf_df <- cbind(lam_pseudotime[[c("monocle3_pseudotime", "celltype_011625")]], t(pseudo_reg[lam2_tfs,]))

plots <- lapply(lam2_tfs, function(x) {
  return(ggplot(lam2_tf_df, aes_string(x="monocle3_pseudotime", y=x, col="celltype_011625")) + geom_point() +
           geom_smooth(method = "loess", color = "black", se=F) + theme_classic()+ ylab(x) +
           xlab("Pseudotime") +  NoLegend())
})
wrap_plots(plots)

lam2_tgs <- c("CXCL12", "VIM", "FAP", "FN1", "DCN", "LOXL1")
lam2_tgs <- subset(cellref_genes, symbol %in% lam2_tgs)

lam2_tg_df <- cbind(lam_pseudotime[[c("monocle3_pseudotime", "celltype_011625")]], 
                    t(pseudo_exp[rownames(lam2_tgs),]))

plots <- apply(lam2_tgs, 1, function(x) {
  return(ggplot(lam2_tg_df, aes_string(x="monocle3_pseudotime", y=x["id"], col="celltype_011625")) + geom_point() +
           geom_smooth(method = "loess", color = "black", se=F) + theme_classic()+ ylab(x["symbol"]) +
           xlab("Pseudotime") +  NoLegend())
})
wrap_plots(plots)



