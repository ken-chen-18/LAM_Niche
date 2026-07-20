library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(ComplexHeatmap)
library(colorRamp2)
library(RColorBrewer)
library(effsize)
library(RobustRankAggreg)


setwd("/global/scratch/users/kenchen/LAM/5_regulon")
input_path <- "/global/scratch/users/kenchen/LAM/input/"
cellref_genes <- readRDS("/global/scratch/users/kenchen/LAM/input/files/CellRef.genes 1.rds")

lam_ctrl <- readRDS("../3_lam_ctrl/lam_ctrl_final.rds")
load("/global/scratch/users/kenchen/LAM/2_uterus/lam_ut_sym.RData")


lam_ctrl_acts <- read.table("reg_results/combined/comb_acts.txt", check.names = FALSE)
lam_ctrl_auc <- read.table("reg_results/combined/auc_mtx_agg.tsv", header=T, 
                           row.names=1, check.names = FALSE)

lam_ctrl[["decoupler"]] <- CreateAssayObject(data =  t(lam_ctrl_acts))
lam_ctrl[["pyscenic"]] <- CreateAssayObject(data =  t(lam_ctrl_auc))

lam_ut_acts <- read.table("reg_results/combined/comb_acts_ut.txt", check.names = FALSE)
lam_ut_auc <- read.table("reg_results/combined/auc_mtx_agg_ut.tsv", header=T, 
                         row.names=1, check.names = FALSE)
lam_ut_auc <- lam_ut_auc[names(lam_ut_acts)]

lam_ut_sym[["decoupler"]] <- CreateAssayObject(data =  t(lam_ut_acts))
lam_ut_sym[["pyscenic"]] <- CreateAssayObject(data =  t(lam_ut_auc))


#wilcoxon
lam_ctrl_meta <- lam_ctrl@meta.data
lam_ut_meta <- lam_ut@meta.data

tfs <- names(lam_ctrl_acts)
names(tfs) <- tfs

nonlam <- grep("LAMCORE.*", unique(lam_ctrl$celltype_011625), value=T, 
               invert = T)


lam1_acts <- lam_ctrl_acts[rownames(subset(lam_ctrl_meta, celltype_011625=="LAMCORE-1")),]
lam2_acts <- lam_ctrl_acts[rownames(subset(lam_ctrl_meta, celltype_011625=="LAMCORE-2")),]
lam_acts <- lam_ctrl_acts[rownames(subset(lam_ctrl_meta, celltype_011625 %in% 
                                            c("LAMCORE-1", "LAMCORE-2", "LAMCORE-3"))),]
nonlam_acts <- lam_ctrl_acts[rownames(subset(lam_ctrl_meta, Diagnosis=="LAM" & 
                                               celltype_011625 %in% nonlam)),]

lam1_auc <- lam_ctrl_auc[rownames(subset(lam_ctrl_meta, celltype_011625=="LAMCORE-1")),]
lam2_auc <- lam_ctrl_auc[rownames(subset(lam_ctrl_meta, celltype_011625=="LAMCORE-2")),]
lam_auc <- lam_ctrl_auc[rownames(subset(lam_ctrl_meta, celltype_011625 %in% 
                                          c("LAMCORE-1", "LAMCORE-2", "LAMCORE-3"))),]
nonlam_auc <- lam_ctrl_auc[rownames(subset(lam_ctrl_meta, Diagnosis=="LAM" & 
                                             celltype_011625 %in% nonlam)),]

lam1_reg <- data.frame(lapply(tfs, function(tf) {
  lam1_act_val <- lam1_acts[[tf]]
  lam2_act_val <- lam2_acts[[tf]]
  nonlam_act_val <- nonlam_acts[[tf]]
  
  w1 <- wilcox.test(lam1_act_val, lam2_act_val, alternative="greater")$p.value
  d1 <- cliff.delta(lam1_act_val, lam2_act_val)$estimate
  
  w2 <- wilcox.test(lam1_act_val, nonlam_act_val, alternative="greater")$p.value
  d2 <- cliff.delta(lam1_act_val, nonlam_act_val)$estimate
  
  lam1_auc_val <- lam1_auc[[tf]]
  lam2_auc_val <- lam2_auc[[tf]]
  nonlam_auc_val <- nonlam_auc[[tf]]
  
  w3 <- wilcox.test(lam1_auc_val, lam2_auc_val, alternative="greater")$p.value
  d3 <- cliff.delta(lam1_auc_val, lam2_auc_val)$estimate
  
  w4 <- wilcox.test(lam1_auc_val, nonlam_auc_val, alternative="greater")$p.value
  d4 <- cliff.delta(lam1_auc_val, nonlam_auc_val)$estimate
  
  return(c(dc_lam12_p_val = w1, dc_lam12_delta = d1,
           sc_lam12_p_val = w3, sc_lam12_delta = d3,
           dc_nonlam_p_val = w2, dc_nonlam_delta = d2,
           sc_nonlam_p_val = w4, sc_nonlam_delta = d4))
}), check.names = F)
  

lam2_reg <- data.frame(lapply(tfs, function(tf) {
  lam1_act_val <- lam1_acts[[tf]]
  lam2_act_val <- lam2_acts[[tf]]
  nonlam_act_val <- nonlam_acts[[tf]]
  
  w1 <- wilcox.test(lam2_act_val, lam1_act_val, alternative="greater")$p.value
  d1 <- cliff.delta(lam2_act_val, lam1_act_val)$estimate
  
  w2 <- wilcox.test(lam2_act_val, nonlam_act_val, alternative="greater")$p.value
  d2 <- cliff.delta(lam2_act_val, nonlam_act_val)$estimate
  
  lam1_auc_val <- lam1_auc[[tf]]
  lam2_auc_val <- lam2_auc[[tf]]
  nonlam_auc_val <- nonlam_auc[[tf]]
  
  w3 <- wilcox.test(lam2_auc_val, lam1_auc_val, alternative="greater")$p.value
  d3 <- cliff.delta(lam2_auc_val, lam1_auc_val)$estimate
  
  w4 <- wilcox.test(lam2_auc_val, nonlam_auc_val, alternative="greater")$p.value
  d4 <- cliff.delta(lam2_auc_val, nonlam_auc_val)$estimate
  
  return(c(dc_lam12_p_val = w1, dc_lam12_delta = d1,
           sc_lam12_p_val = w3, sc_lam12_delta = d3,
           dc_nonlam_p_val = w2, dc_nonlam_delta = d2,
           sc_nonlam_p_val = w4, sc_nonlam_delta = d4))
}), check.names = F)

lam_reg <- data.frame(lapply(tfs, function(tf) {
  lam_acts_val <- lam_acts[[tf]]
  nonlam_acts_val <- nonlam_acts[[tf]]
  w1 <- wilcox.test(lam_acts_val, nonlam_acts_val, alternative="greater")$p.value
  d1 <- cliff.delta(lam_acts_val, nonlam_acts_val)$estimate
  
  lam_auc_val <- lam_auc[[tf]]
  nonlam_auc_val <-nonlam_auc[[tf]]
  w2 <- wilcox.test(lam_auc_val, nonlam_auc_val, alternative="greater")$p.value
  d2 <- cliff.delta(lam_auc_val, nonlam_auc_val)$estimate
  
  return(c(dc_p_val = w1, dc_delta = d1,
           sc_p_val = w2, sc_delta = d2))
}), check.names = F)


lam_at2_acts <- lam_ctrl_acts[rownames(subset(lam_ctrl_meta, Diagnosis=="LAM" & celltype_011625 == "AT2")),]
ctrl_at2_acts <- lam_ctrl_acts[rownames(subset(lam_ctrl_meta, Diagnosis=="Control" & celltype_011625 == "AT2")),]
lam_at2_auc <- lam_ctrl_auc[rownames(subset(lam_ctrl_meta, Diagnosis=="LAM" & celltype_011625 == "AT2")),]
ctrl_at2_auc <- lam_ctrl_auc[rownames(subset(lam_ctrl_meta, Diagnosis=="Control" & celltype_011625 == "AT2")),]

at2_reg <- data.frame(lapply(tfs, function(tf) {
  lam_acts_val <- lam_at2_acts[[tf]]
  ctrl_acts_val <- ctrl_at2_acts[[tf]]
  w1 <- wilcox.test(lam_acts_val, ctrl_acts_val, alternative="two.sided")$p.value
  d1 <- cliff.delta(lam_acts_val, ctrl_acts_val)$estimate
  
  lam_auc_val <- lam_at2_auc[[tf]]
  ctrl_auc_val <- ctrl_at2_auc[[tf]]
  w2 <- wilcox.test(lam_auc_val, ctrl_auc_val, alternative="two.sided")$p.value
  d2 <- cliff.delta(lam_auc_val, ctrl_auc_val)$estimate
  
  return(c(dc_p_val = w1, dc_delta = d1,
           sc_p_val = w2, sc_delta = d2))
}), check.names = F)

lam_af1_acts <- lam_ctrl_acts[rownames(subset(lam_ctrl_meta, Diagnosis=="LAM" & celltype_011625 == "AF1")),]
ctrl_af1_acts <- lam_ctrl_acts[rownames(subset(lam_ctrl_meta, Diagnosis=="Control" & celltype_011625 == "AF1")),]
lam_af1_auc <- lam_ctrl_auc[rownames(subset(lam_ctrl_meta, Diagnosis=="LAM" & celltype_011625 == "AF1")),]
ctrl_af1_auc <- lam_ctrl_auc[rownames(subset(lam_ctrl_meta, Diagnosis=="Control" & celltype_011625 == "AF1")),]

af1_reg <- data.frame(lapply(tfs, function(tf) {
  lam_acts_val <- lam_af1_acts[[tf]]
  ctrl_acts_val <- ctrl_af1_acts[[tf]]
  w1 <- wilcox.test(lam_acts_val, ctrl_acts_val, alternative="two.sided")$p.value
  d1 <- cliff.delta(lam_acts_val, ctrl_acts_val)$estimate
  
  lam_auc_val <- lam_af1_auc[[tf]]
  ctrl_auc_val <- ctrl_af1_auc[[tf]]
  w2 <- wilcox.test(lam_auc_val, ctrl_auc_val, alternative="two.sided")$p.value
  d2 <- cliff.delta(lam_auc_val, ctrl_auc_val)$estimate
  
  return(c(dc_p_val = w1, dc_delta = d1,
           sc_p_val = w2, sc_delta = d2))
}), check.names = F)

ut_tfs <- names(lam_ut_acts)
names(ut_tfs) <- ut_tfs
lam_acts <- lam_ut_acts[rownames(subset(lam_ut_meta, celltype_4=="LAM")),]
nonlam_acts <- lam_ut_acts[rownames(subset(lam_ut_meta, celltype_4 != "LAM")),]
lam_auc <- lam_ut_auc[rownames(subset(lam_ut_meta, celltype_4=="LAM")),]
nonlam_auc <- lam_ut_auc[rownames(subset(lam_ut_meta, celltype_4 != "LAM")),]

lam_ut_reg <- data.frame(lapply(ut_tfs, function(tf) {

  lam_act_val <- lam_acts[[tf]]
  nonlam_act_val <- nonlam_acts[[tf]]

  w1 <- wilcox.test(lam_act_val, nonlam_act_val, alternative="greater")$p.value
  d1 <- cliff.delta(lam_act_val, nonlam_act_val)$estimate
  
  lam_auc_val <- lam_auc[[tf]]
  nonlam_auc_val <- nonlam_auc[[tf]]
  
  w2 <- wilcox.test(lam_auc_val, nonlam_auc_val, alternative="greater")$p.value
  d2 <- cliff.delta(lam_auc_val, nonlam_auc_val)$estimate
  
  return(c(dc_p_val = w1, dc_delta = d1,
           sc_p_val = w2, sc_delta = d2))
}), check.names = F)

at2_reg <- data.frame(t(at2_reg), check.names = F)
af1_reg <- data.frame(t(af1_reg), check.names = F)
lam1_reg <- data.frame(t(lam1_reg), check.names = F)
lam2_reg <- data.frame(t(lam2_reg), check.names = F)
lam_reg <- data.frame(t(lam_reg), check.names = F)
lam_ut_reg <- data.frame(t(lam_ut_reg), check.names = F)



lam_pct <- FoldChange(subset(lam_ctrl, Diagnosis=="LAM"), 
                      ident.1=c("LAMCORE-1", "LAMCORE-2", "LAMCORE-3"),
                      group.by="celltype_011625")
lam12_pct <- FoldChange(lam_ctrl, ident.1="LAMCORE-1", ident.2="LAMCORE-2",
                      group.by="celltype_011625")
at2_pct <- FoldChange(lam_ctrl, ident.1="LAM_AT2", ident.2="Control_AT2",
                      group.by="celltype_combined")
af1_pct <- FoldChange(lam_ctrl, ident.1="LAM_AF1", ident.2="Control_AF1",
                      group.by="celltype_combined")
lam_ut_pct <- FoldChange(lam_ut_sym, ident.1="LAM",
                      group.by="celltype_4")

lam_pct$symbol <- cellref_genes[rownames(lam_pct), "symbol"]
lam12_pct$symbol <- cellref_genes[rownames(lam12_pct), "symbol"]
at2_pct$symbol <- cellref_genes[rownames(at2_pct), "symbol"]
af1_pct$symbol <- cellref_genes[rownames(af1_pct), "symbol"]


lam_reg[c("LAM.pct", "nonLAM.pct")] <- lam_pct[match(rownames(lam_reg), lam_pct$symbol), 
                                               c("pct.1", "pct.2")]
at2_reg[c("LAM.pct", "Control.pct")] <- at2_pct[match(rownames(at2_reg), at2_pct$symbol), 
                                               c("pct.1", "pct.2")]
af1_reg[c("LAM.pct", "Control.pct")] <- af1_pct[match(rownames(af1_reg), af1_pct$symbol), 
                                               c("pct.1", "pct.2")]
lam_ut_reg[c("LAM.pct", "nonLAM.pct")] <- lam_ut_pct[match(rownames(lam_ut_reg), rownames(lam_ut_pct)), 
                                               c("pct.1", "pct.2")]

lam1_reg[c("LAM1.pct", "LAM2.pct", "nonLAM.pct")] <- cbind(lam12_pct[match(rownames(lam1_reg), lam12_pct$symbol), 
                                                c("pct.1", "pct.2")],
                                                lam_pct[match(rownames(lam1_reg), lam_pct$symbol), 
                                                        "pct.2"])
lam2_reg[c("LAM2.pct", "LAM1.pct", "nonLAM.pct")] <- lam1_reg[c("LAM2.pct", "LAM1.pct", "nonLAM.pct")]


lam_reg <- subset(lam_reg, dc_delta>0 & sc_delta>0 & LAM.pct>0.1 & 
                    dc_p_val<0.05 & sc_p_val<0.05)
lam_ut_reg <- subset(lam_ut_reg, dc_delta>0 & sc_delta>0 & LAM.pct>0.1 & 
                       dc_p_val<0.05 & sc_p_val<0.05)

lam1_reg <- subset(lam1_reg, dc_lam12_delta>0 & sc_lam12_delta>0 &
                     dc_nonlam_delta>0 & sc_nonlam_delta>0 & LAM1.pct>0.1 & 
                     dc_lam12_p_val<0.05 & sc_lam12_p_val<0.05 &
                     dc_nonlam_p_val<0.05 & sc_nonlam_p_val<0.05)
lam2_reg <- subset(lam2_reg, dc_lam12_delta>0 & sc_lam12_delta>0 &
                     dc_nonlam_delta>0 & sc_nonlam_delta>0 & LAM2.pct>0.1 & 
                     dc_lam12_p_val<0.05 & sc_lam12_p_val<0.05 &
                     dc_nonlam_p_val<0.05 & sc_nonlam_p_val<0.05)



lam_reg$pct.fc <- lam_reg$LAM.pct/lam_reg$nonLAM.pct
lam_ut_reg$pct.fc <- lam_ut_reg$LAM.pct/lam_ut_reg$nonLAM.pct
at2_reg$pct.fc <- at2_reg$LAM.pct/at2_reg$Control.pct
af1_reg$pct.fc <- af1_reg$LAM.pct/af1_reg$Control.pct
lam1_reg$pct.fc.12 <- lam1_reg$LAM1.pct/lam1_reg$LAM2.pct
lam1_reg$pct.fc.nl <- lam1_reg$LAM1.pct/lam1_reg$nonLAM.pct
lam2_reg$pct.fc.12 <- lam2_reg$LAM2.pct/lam2_reg$LAM1.pct
lam2_reg$pct.fc.nl <- lam2_reg$LAM2.pct/lam2_reg$nonLAM.pct

at2_lam_reg <- subset(at2_reg, dc_delta>0 & sc_delta>0 & LAM.pct>0.1 & 
                        dc_p_val<0.025 & sc_p_val<0.025)
at2_ctrl_reg <- subset(at2_reg, dc_delta<0 & sc_delta<0 & Control.pct>0.1 & 
                         dc_p_val<0.025 & sc_p_val<0.025)
af1_lam_reg <- subset(af1_reg, dc_delta>0 & sc_delta>0 & LAM.pct>0.1 & 
                        dc_p_val<0.025 & sc_p_val<0.025)
af1_ctrl_reg <- subset(af1_reg, dc_delta<0 & sc_delta<0 & Control.pct>0.1 & 
                         dc_p_val<0.025 & sc_p_val<0.025)



rank_rra <- function(df, increasing_cols, decreasing_cols) {
  gene_names <- rownames(df)
  
  # Helper to create ranked name vector
  rank_column <- function(values, decreasing = FALSE) {
    gene_names[order(values, decreasing = decreasing)]
  }
  
  # Build list of rankings
  ranked_lists <- lapply(increasing_cols, function(col) rank_column(df[, col], decreasing = FALSE))
  ranked_lists <- c(ranked_lists, lapply(decreasing_cols, function(col) rank_column(df[, col], decreasing = TRUE)))
  
  # Run RRA
  aggregateRanks(ranked_lists)
}

lam_rra <- rank_rra(lam_reg, c(1,3), c(2,4,5,7))
lam_ut_rra <- rank_rra(lam_ut_reg, c(1,3), c(2,4,5,7))
at2_lam_rra <- rank_rra(at2_lam_reg, c(1,3), c(2,4,5,7))
at2_ctrl_rra <- rank_rra(at2_ctrl_reg, c(1,2,3,4,7), 6)
af1_lam_rra <- rank_rra(af1_lam_reg, c(1,3), c(2,4,5,7))
af1_ctrl_rra <- rank_rra(af1_ctrl_reg, c(1,2,3,4,7), 6)
lam1_rra <- rank_rra(lam1_reg, c(1,3,5,7), c(2,4,6,8,9,12,13))
lam2_rra <- rank_rra(lam2_reg, c(1,3,5,7), c(2,4,6,8,9,12,13))

lam_reg <- cbind(lam_reg[rownames(lam_rra),], lam_rra["Score"])
lam1_reg <- cbind(lam1_reg[rownames(lam1_rra),], lam1_rra["Score"])
lam2_reg <- cbind(lam2_reg[rownames(lam2_rra),], lam2_rra["Score"])
lam_ut_reg <- cbind(lam_ut_reg[rownames(lam_ut_rra),], lam_ut_rra["Score"])
at2_lam_reg <- cbind(at2_lam_reg[rownames(at2_lam_rra),], at2_lam_rra["Score"])
at2_ctrl_reg <- cbind(at2_ctrl_reg[rownames(at2_ctrl_rra),], at2_ctrl_rra["Score"])
af1_lam_reg <- cbind(af1_lam_reg[rownames(af1_lam_rra),], af1_lam_rra["Score"])
af1_ctrl_reg <- cbind(af1_ctrl_reg[rownames(af1_ctrl_rra),], af1_ctrl_rra["Score"])

save(lam_reg, af1_lam_reg, af1_ctrl_reg, at2_lam_reg, at2_ctrl_reg, lam1_reg, lam2_reg, lam_ut_reg, file="reg.RData")


write.table(lam_ut_reg, file = "reg_rank/lam_ut_reg.csv", row.names=T, col.names=NA,
            sep=",", quote=F)




DefaultAssay(lam_ctrl) <- "RNA"

at2_obj <- subset(lam_ctrl, celltype_011625 == "AT2")
af1_obj <- subset(lam_ctrl, celltype_011625 =="AF1")
Idents(at2_obj) <- "Diagnosis"
Idents(af1_obj) <- "Diagnosis"



af1_obj <- RunUMAP(af1_obj, dims=1:50)
at2_obj <- RunUMAP(at2_obj, dims=1:50)
DimPlot(af1_obj)
DimPlot(at2_obj)

saveRDS(at2_obj, file="at2_obj.rds")
saveRDS(af1_obj, file="af1_obj.rds")
saveRDS(lam_ctrl, file="lam_ctrl.rds")
saveRDS(lam_ut, file="lam_ut.rds")



lam_ctrl_acts <- read.table("../5_regulon/reg_results/combined/comb_acts.txt", check.names = FALSE)
lam_ctrl_auc <- read.table("../5_regulon/reg_results/combined/auc_mtx_agg.tsv", header=T, 
                           row.names=1, check.names = FALSE)

lam_ctrl_meta <- lam_ctrl@meta.data
lam1_acts <- lam_ctrl_acts[rownames(subset(lam_ctrl_meta, celltype_011625=="LAMCORE-1")),]
lam2_acts <- lam_ctrl_acts[rownames(subset(lam_ctrl_meta, celltype_011625=="LAMCORE-2")),]
lam_acts <- lam_ctrl_acts[rownames(subset(lam_ctrl_meta, celltype_011625 %in% 
                                            c("LAMCORE-1", "LAMCORE-2", "LAMCORE-3"))),]
nonlam_acts <- lam_ctrl_acts[rownames(subset(lam_ctrl_meta, Diagnosis=="LAM" & 
                                               !(celltype_011625 %in% 
                                                   c("LAMCORE-1", "LAMCORE-2", "LAMCORE-3")))),]
lam3_acts <- lam_ctrl_acts[rownames(subset(lam_ctrl_meta, celltype_011625=="LAMCORE-3")),]

lam1_auc <- lam_ctrl_auc[rownames(subset(lam_ctrl_meta, celltype_011625=="LAMCORE-1")),]
lam2_auc <- lam_ctrl_auc[rownames(subset(lam_ctrl_meta, celltype_011625=="LAMCORE-2")),]
lam_auc <- lam_ctrl_auc[rownames(subset(lam_ctrl_meta, celltype_011625 %in% 
                                          c("LAMCORE-1", "LAMCORE-2", "LAMCORE-3"))),]
nonlam_auc <- lam_ctrl_auc[rownames(subset(lam_ctrl_meta, Diagnosis=="LAM" & 
                                             !(celltype_011625 %in% 
                                                 c("LAMCORE-1", "LAMCORE-2", "LAMCORE-3")))),]
lam3_auc <- lam_ctrl_auc[rownames(subset(lam_ctrl_meta, celltype_011625=="LAMCORE-3")),]

tfs <- names(lam_ctrl_acts)
names(tfs) <- tfs
lam3_reg <- data.frame(lapply(tfs, function(tf) {
  lam1_act_val <- lam1_acts[[tf]]
  lam2_act_val <- lam2_acts[[tf]]
  lam3_act_val <- lam3_acts[[tf]]
  nonlam_act_val <- nonlam_acts[[tf]]
  
  w1 <- wilcox.test(lam3_act_val, lam1_act_val, alternative="greater")$p.value
  d1 <- cliff.delta(lam3_act_val, lam1_act_val)$estimate
  w2 <- wilcox.test(lam3_act_val, lam2_act_val, alternative="greater")$p.value
  d2 <- cliff.delta(lam3_act_val, lam2_act_val)$estimate
  
  w3 <- wilcox.test(lam3_act_val, nonlam_act_val, alternative="greater")$p.value
  d3 <- cliff.delta(lam3_act_val, nonlam_act_val)$estimate
  
  lam1_auc_val <- lam1_auc[[tf]]
  lam2_auc_val <- lam2_auc[[tf]]
  lam3_auc_val <- lam3_auc[[tf]]
  nonlam_auc_val <- nonlam_auc[[tf]]
  
  w4 <- wilcox.test(lam3_auc_val, lam1_auc_val, alternative="greater")$p.value
  d4 <- cliff.delta(lam3_auc_val, lam1_auc_val)$estimate
  w5 <- wilcox.test(lam3_auc_val, lam2_auc_val, alternative="greater")$p.value
  d5 <- cliff.delta(lam3_auc_val, lam2_auc_val)$estimate
  
  w6 <- wilcox.test(lam3_auc_val, nonlam_auc_val, alternative="greater")$p.value
  d6 <- cliff.delta(lam3_auc_val, nonlam_auc_val)$estimate
  
  return(c(dc_lam31_p_val = w1, dc_lam31_delta = d1,
           sc_lam31_p_val = w4, sc_lam31_delta = d4,
           dc_lam32_p_val = w2, dc_lam32_delta = d2,
           sc_lam32_p_val = w5, sc_lam32_delta = d5,
           dc_nonlam_p_val = w3, dc_nonlam_delta = d3,
           sc_nonlam_p_val = w6, sc_nonlam_delta = d6))
}), check.names = F)

lam3_reg <- lam3_reg[9:12]
lam3_reg$id <- cellref_genes[match(rownames(lam3_reg), cellref_genes$symbol), "id"]
pct <- FoldChange(lam_lung_nonim, features=lam3_reg$id, group.by="LAM", ident.1="LAMCORE-3",
                  ident.2="Other")
lam3_reg$LAM3.pct <- pct[match(lam3_reg$id, rownames(pct)), "pct.1"]
lam3_reg$Other.pct <- pct[match(lam3_reg$id, rownames(pct)), "pct.2"]
lam3_reg$pct.fc <- lam3_reg$LAM3.pct/lam3_reg$Other.pct
lam3_reg <- subset(lam3_reg, LAM3.pct>0.1)

genes <- rownames(lam3_reg)

rank_lists <- list(
  dc_p = genes[order(lam3_reg$dc_nonlam_p_val, decreasing = FALSE)],   # lower better
  sc_p = genes[order(lam3_reg$sc_nonlam_p_val, decreasing = FALSE)],
  
  dc_delta = genes[order(lam3_reg$dc_nonlam_delta, decreasing = TRUE)], # higher better
  sc_delta = genes[order(lam3_reg$sc_nonlam_delta, decreasing = TRUE)],
  
  lam3_pct = genes[order(lam3_reg$LAM3.pct, decreasing = TRUE)],
  pct_fc   = genes[order(lam3_reg$pct.fc,   decreasing = TRUE)]
)
rra_res <- aggregateRanks(rank_lists)

lam3_reg$rra <- rra_res[rownames(lam3_reg), "Score"]

lam3_reg <- data.frame(t(lam3_reg), check.names = F)
save(lam3_reg, file="reg_rank/lam3_reg.RData")
write.table(lam3_reg, file="reg_rank/lam3_reg.txt", sep="\t", quote=F)
