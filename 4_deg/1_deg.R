options(encoding = "UTF-8")
.libPaths( c('/global/scratch/users/kenchen/R_4.2', .libPaths()) )
library(Seurat)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

setwd("/global/scratch/users/kenchen/LAM/4_deg")
input_path <- "/global/scratch/users/kenchen/LAM/input/"
cellref_genes <- readRDS("/global/scratch/users/kenchen/LAM/input/files/CellRef.genes 1.rds")

lam_ctrl <- readRDS("../3_lam_ctrl/lam_ctrl_final.rds")
lam_lung_nonim <- readRDS("../1_lam_int/lam_lung_nonim.rds")
load("../2_uterus/lam_ut.RData")
#load("lam_ut.RData")

de_helper <- function(obj, ident.1, ident.2=NULL, only.pos=T, p=0.01, min.pct=0.1, 
                      logfc.threshold=log2(1.5), group.by="celltype_011625") {
  mark <- FindMarkers(obj, ident.1=ident.1, ident.2=ident.2, min.pct=min.pct,
                      logfc.threshold=logfc.threshold, only.pos=only.pos,
                      group.by=group.by)
  de <- subset(mark, p_val < p)
  de$symbol <- cellref_genes[rownames(de), "symbol"]
  de
}


lam_lung_nonim <- NormalizeData(lam_lung_nonim)
Idents(lam_lung_nonim) <- "celltype_011625"

nonlam <- grep("LAMCORE.*", unique(lam_lung_nonim$celltype_011625), value=T, 
               invert = T)

lam1 <- de_helper(lam_lung_nonim, ident.1="LAMCORE-1", ident.2=nonlam)
lam2 <- de_helper(lam_lung_nonim, ident.1="LAMCORE-2", ident.2=nonlam)
lam3 <- de_helper(lam_lung_nonim, ident.1="LAMCORE-3", ident.2=nonlam)
lam_all <- de_helper(lam_lung_nonim, ident.1=c("LAMCORE-1", "LAMCORE-2", "LAMCORE-3"), 
                     ident.2=nonlam)
lam12 <- de_helper(lam_lung_nonim, ident.1="LAMCORE-1", ident.2="LAMCORE-2",
                   only.pos=F)

Idents(lam_ut) <- "celltype_011625"
lam_ut <- NormalizeData(lam_ut)
lam_ut_deg <- de_helper(subset(lam_ut, celltype_011625 != "Uterine Immune"), 
                        ident.1="Uterine LAMCORE")


l1 <- FoldChange(lam_lung_nonim, ident.1="LAMCORE-1")
l1 <- subset(l1, pct.1>0.1)
l2 <- FoldChange(lam_lung_nonim, ident.1="LAMCORE-2")
l2 <- subset(l2, pct.1>0.1)
l3 <- FoldChange(lam_lung_nonim, ident.1="LAMCORE-3")
l3 <- subset(l3, pct.1>0.1)

ueg_cand <- unique(c(rownames(l1), rownames(l2), rownames(l3)))

lam_plot <- DotPlot(lam_lung_nonim, features =ueg_cand, group.by="celltype_011625")$data

lam_plot$ens_id <- as.character(lam_plot$features.plot)
lam_plot$symbol <- lam_genes[lam_plot$ens_id, "symbol"]
lam_plot$features.plot <- NULL
lam_data <- transform(subset(lam_plot, id == "LAMCORE-1"), lam1.pct=pct.exp, pct.exp=NULL)
lam_data$lam2.pct <- subset(lam_plot, id == "LAMCORE-2")$pct.exp
lam_data$lam3.pct <- subset(lam_plot, id == "LAMCORE-3")$pct.exp

max_nonlam <- subset(lam_plot, id %in% nonlam) %>% group_by(ens_id) %>% summarize(max = max(pct.exp))
lam_data$max_nonlam <- max_nonlam[match(lam_data$ens_id, max_nonlam$ens_id),]$max
rownames(lam_data) <- lam_data$ens_id
lam_data$symbol <- cellref_genes[lam_data$ens_id, "symbol"]
lam_data[c("avg.exp", "id", "avg.exp.scaled", "ens_id")] <- NULL

lam_data$l1_v_lc <- lam_data$lam1.pct/pmax(lam_data$lam2.pct, lam_data$lam3.pct)
lam_data$l1_v_nl <- lam_data$lam1.pct/lam_data$max_nonlam
lam_data$l2_v_lc <- lam_data$lam2.pct/pmax(lam_data$lam1.pct, lam_data$lam3.pct)
lam_data$l2_v_nl <- lam_data$lam2.pct/lam_data$max_nonlam
lam_data$l3_v_lc <- lam_data$lam3.pct/pmax(lam_data$lam1.pct, lam_data$lam2.pct)
lam_data$l3_v_nl <- lam_data$lam3.pct/lam_data$max_nonlam

save(lam1, lam2, lam12, lam_data, lam_plot, file="lam_deg.RData")


lam_ctrl <- NormalizeData(lam_ctrl)

cts <- setdiff(nonlam, c("Mesothelial", "Goblet"))
names(cts) <- cts

lam_ctrl_deg <- lapply(cts, function(x) {
  obj <- subset(lam_ctrl, celltype_011625 == x)
  deg <- de_helper(obj, ident.1="LAM", ident.2="Control", group.by="Diagnosis",
                   only.pos = F)
  return(deg)
})


save(lam_ctrl_deg, file="lam_ctrl_deg.RData")

lapply(names(lam_ctrl_deg), function(ct) {
  deg <- lam_ctrl_deg[[ct]]
  write.table(deg, file = paste0("lam_ctrl_", ct, "_deg.csv"), row.names=T, col.names=NA,
              sep=",", quote=F)
})

write.table(lam_data, file = "lam_pct_df.csv", row.names=T, col.names=NA,
            sep=",", quote=F)




###########
#GSEA

lam_genelist <- lam12$avg_log2FC
names(lam_genelist) <- as.character(lam12$symbol)
lam_genelist <- sort(lam_genelist, decreasing=T)

lam_ent <- bitr(names(lam_genelist), fromType = "SYMBOL", toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)
lam_genelist <- lam_genelist[names(lam_genelist) %in% lam_ent$SYMBOL]
names(lam_genelist) <- lam_ent$ENTREZID



gsego_lam <- gseGO(geneList = lam_genelist,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pvalueCutoff=0.5)



gsekegg_lam <- gseKEGG(geneList = lam_genelist,
                       keyType = "ncbi-geneid",
                       organism = "hsa",
                       pvalueCutoff=0.5,)


gsecomb <- gsego_lam
gsecomb@result <- rbind(gsego_lam@result, gsekegg_lam@result)
gsecomb@geneSets <- c(gsego_lam@geneSets, gsekegg_lam@geneSets)

save(gsego_lam, gsekegg_lam, gsecomb, file="gse.RData")




#LAMCORE-3
lam_lung_nonim <- NormalizeData(lam_lung_nonim)
lam_lung_nonim$LAM <- ifelse(lam_lung_nonim$celltype_011625 %in% paste0("LAMCORE-", 1:3),
                             lam_lung_nonim$celltype_011625, "Other")
lam3 <- FindMarkers(lam_lung_nonim, group.by="LAM", ident.1="LAMCORE-3", ident.2="Other", min.pct=0.1,
                    logfc.threshold = log2(1.5), only.pos=T)
lam3$symbol <- cellref_genes[rownames(lam3), "symbol"]
lam3 <- subset(lam3, p_val<0.01)

lam3v2 <- FindMarkers(lam_lung_nonim, group.by="LAM", ident.1="LAMCORE-3", ident.2="LAMCORE-2", min.pct=0.1,
                      logfc.threshold = log2(1.5))
lam3v2$symbol <- cellref_genes[rownames(lam3v2), "symbol"]
lam3v2$pct.fc <- lam3v2$pct.1/lam3v2$pct.2
lam3v2 <- subset(lam3v2, p_val<0.01)

lam3v1 <- FindMarkers(lam_lung_nonim, group.by="LAM", ident.1="LAMCORE-3", ident.2="LAMCORE-1", min.pct=0.1,
                      logfc.threshold = log2(1.5))
lam3v1$symbol <- cellref_genes[rownames(lam3v1), "symbol"]
lam3v1$pct.fc <- lam3v1$pct.1/lam3v1$pct.2
lam3v1 <- subset(lam3v1, p_val<0.01)

lam_lung_nonim$lam_lineage <- ifelse(lam_lung_nonim$celltype_011625 %in% c("LAMCORE-1", "LAMCORE-2", "LAMCORE-3"),
                                     lam_lung_nonim$celltype_011625, lam_lung_nonim$lineage_level1)
lam_lung_nonim$lam_lineage <- factor(lam_lung_nonim$lam_lineage, levels=rev(c("LAMCORE-1", "LAMCORE-2", "LAMCORE-3",
                                                                              "Mesenchymal", "Endothelial", "Epithelial",
                                                                              "Immune")))
top3v2 <- slice_max(lam3v2, avg_log2FC, n=25)
top3v1 <- slice_max(lam3v1, avg_log2FC, n=25)
top3 <- slice_max(lam3, avg_log2FC, n=25)

DotPlot(lam_lung_nonim,
        features = rownames(top3v2),
        group.by="lam_lineage",
        dot.scale = 3, dot.min=0.01) +
  scale_x_discrete(labels=top3v2$symbol) +
  scale_color_gradient2(low = "blue", mid="gray", high ="red") +
  theme(axis.text = element_text(size = 10)) +
  scale_y_discrete(limits=rev)+
  scale_radius(breaks=c(25,50,75, 100), limits=c(0,100), range=c(0,6))+
  theme(axis.title = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5),)+
  labs(title=NULL)

DotPlot(lam_lung_nonim,
        features = rownames(top3v1),
        group.by="lam_lineage",
        dot.scale = 3, dot.min=0.01) +
  scale_x_discrete(labels=top3v1$symbol) +
  scale_color_gradient2(low = "blue", mid="gray", high ="red") +
  theme(axis.text = element_text(size = 10)) +
  scale_y_discrete(limits=rev)+
  scale_radius(breaks=c(25,50,75, 100), limits=c(0,100), range=c(0,6))+
  theme(axis.title = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5),)+
  labs(title=NULL)

DotPlot(lam_lung_nonim,
        features = rownames(top3),
        group.by="lam_lineage",
        dot.scale = 3, dot.min=0.01) +
  scale_x_discrete(labels=top3$symbol) +
  scale_color_gradient2(low = "blue", mid="gray", high ="red") +
  theme(axis.text = element_text(size = 10)) +
  scale_y_discrete(limits=rev)+
  scale_radius(breaks=c(25,50,75, 100), limits=c(0,100), range=c(0,6))+
  theme(axis.title = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5),)+
  labs(title=NULL)

write.table(lam3v2, file="lam3v2_deg.txt", sep="\t", quote=F)
write.table(lam3v1, file="lam3v1_deg.txt", sep="\t", quote=F)
