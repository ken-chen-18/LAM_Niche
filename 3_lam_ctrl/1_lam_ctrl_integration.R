options(encoding = "UTF-8")
.libPaths( c('/global/scratch/users/kenchen/R_4.2', .libPaths()) )
library(Matrix, lib.loc = "/global/scratch/users/kenchen/R_4.2")
library(ggplot2)
library(Seurat)
library(dplyr)
library(tidyverse)
library(harmony)

options(future.globals.maxSize = 1.4 * 1024^4)

setwd("/global/scratch/users/kenchen/LAM/3_lam_ctrl")
input_path <- "/global/scratch/users/kenchen/LAM/input/"

ctrl_genes <-  readRDS(paste0(input_path, "files/CellRef.genes 1.rds"))
ctrl <- readRDS("../0_cihan/ctrl.rds")
lam_lung_nonim <- readRDS("../1_lam_int/lam_lung_nonim.rds")
samp_meta <- read.table(paste0(input_path, "files/sample_meta.csv"),sep=",",
                        header=T)


lam_lung_nonim$Age <- lam_lung_nonim$age
lam_lung_nonim$Sex <- lam_lung_nonim$sex
lam_lung_nonim$Tissue <- "Lung"

lam_lung_nonim@meta.data[c("nCount_SoupX", "nFeature_SoupX", "nCount_SoupX_SCT",
                           "nFeature_SoupX_SCT", "age", "sex", "cellref_seed_celltype_level3",
                           "cellref_seed_lineage_level1", "cellref_celltype_level3",
                           "cellref_lineage_level1", "celltype_072923", "seurat_clusters",
                           "integrated_snn_res.4.2", "prediction.score.max")] <- NULL



ctrl$Barcode <- paste0(gsub(".*?([ACGT]{5,}).*", "\\1", colnames(ctrl)), "-1")
ctrl <- RenameCells(ctrl, new.names = paste(ctrl$DataID, ctrl$Barcode, sep = "_"))
ctrl$Tissue <- "Lung"
ctrl$Modality <- "RNA"
ctrl$Diagnosis <- "Control"
ctrl$Technology <- samp_meta[match(ctrl$DataID, samp_meta$DataID), "Technology"]
ctrl$DataID2 <- factor(samp_meta[match(ctrl$DataID, samp_meta$DataID), "DataID2"], 
                       levels=paste0("Ctrl", 1:16))
ctrl$Sex <- "Female"
ctrl$celltype_011625 <- ctrl$celltype_level3
ctrl$celltype_011625[ctrl$celltype_011625 == "SCMF"] <- "MyoFB"
ctrl$celltype_011625 <- ifelse(ctrl$celltype_011625 %in% c("Deuterosomal", "Ionocyte",
                                                           "Suprabasal", "Tuft"),
                               "Rare", ctrl$celltype_011625)

ctrl@meta.data[c("lineage_level2", "condition", "celltype_level1", "celltype_level2",
                 "celltype_level3", "celltype_level3_fullname")] <- NULL

DefaultAssay(ctrl) <- "RNA"

ctrl_nonim <- subset(ctrl, lineage_level1 != "Immune")
ctrl_nonim <- subset(ctrl_nonim, nFeature_RNA>500 & nFeature_RNA<7500 & nCount_RNA<40000 & 
                       pMT<20)

ctrl_counts <- GetAssayData(ctrl_nonim, slot="counts")
rownames(ctrl_counts) <- rownames(ctrl_genes)
ctrl_nonim_ens <- CreateSeuratObject(ctrl_counts, meta.data=ctrl_nonim@meta.data)
ctrl_nonim_ens@meta.data <- ctrl_nonim@meta.data

lam_list <- c(SplitObject(lam_lung_nonim, split.by="DataID"), SplitObject(ctrl_nonim_ens, split.by="DataID"))

set.seed(2023)
lam_list <- lapply(X = lam_list, FUN = function(o) {
  o = CreateSeuratObject(counts = o@assays[["RNA"]]@counts, meta.data = o@meta.data)
  o = NormalizeData(o) %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    SCTransform(vars.to.regress = c("pMT", "S.Score", "G2M.Score"))
  o
})


features = SelectIntegrationFeatures(object.list = lam_list, nfeatures = 2000)
lam_list = PrepSCTIntegration(object.list = lam_list, anchor.features=features)
lam_list <- lapply(X = lam_list, FUN = RunPCA, features = features)
anchors = FindIntegrationAnchors(lam_list, anchor.features = features,
                                 reduction = "rpca",
                                 k.anchor=15,
                                 normalization.method = "SCT")

lam_ctrl_int = IntegrateData(anchors, normalization.method = "SCT",
                             k.weight=50)
DefaultAssay(lam_ctrl_int) = "integrated"

lam_ctrl_int = RunPCA(lam_ctrl_int)
lam_ctrl_int = RunUMAP(lam_ctrl_int, dims = 1:50)

DefaultAssay(lam_ctrl_int) <- "RNA"

lam_ctrl_int$celltype_combined <- paste0(lam_ctrl_int$Diagnosis, "_", lam_ctrl_int$celltype_011625)
Idents(lam_ctrl_int) <- "celltype_011625"
lam_ctrl_int <- NormalizeData(lam_ctrl_int)

saveRDS(ctrl, file="ctrl_clean.rds")
saveRDS(ctrl_nonim_ens, file = "ctrl_nonim.rds")
saveRDS(lam_ctrl_int, file = "lam_ctrl_final.rds")

