setwd("X:/cakar/2023_LAM")
rm(list=ls())
library(Seurat)
library(dplyr)
library(tidyverse)
library(reticulate)

options(reticulate.conda_binary = "C:/Users/caksf7/AppData/Local/anaconda3/condabin/conda.bit")
use_condaenv("py8")
leidenalg = import("leidenalg")

set.seed(2023)

sc_obj = readRDS("saved_objects/lam_integrate_rpca.rds")

#import LAM1_scRNA
raw_counts = Read10X_h5(filename = "imported_data/old_lam_data/LAM1_scRNA_raw_feature_bc_matrix.h5", use.names = FALSE)
LAM1_scRNA = CreateSeuratObject(counts = raw_counts)
LAM1_scRNA = LAM1_scRNA[,sc_obj[,WhichCells(sc_obj, expression=DataID=="LAM1_scRNA")]$Barcode]
all(subset(sc_obj, DataID=="LAM1_scRNA")$Barcode %in% colnames(LAM1_scRNA))
all(colnames(LAM1_scRNA) %in% subset(sc_obj, DataID=="LAM1_scRNA")$Barcode)

#import LAM3_scRNA
raw_counts = Read10X_h5(filename = "imported_data/old_lam_data/LAM3_scRNA_raw_feature_bc_matrix.h5", use.names = FALSE)
LAM3_scRNA = CreateSeuratObject(counts = raw_counts)
LAM3_scRNA = LAM3_scRNA[,sc_obj[,WhichCells(sc_obj, expression=DataID=="LAM3_scRNA")]$Barcode]
all(subset(sc_obj, DataID=="LAM3_scRNA")$Barcode %in% colnames(LAM3_scRNA))
all(colnames(LAM3_scRNA) %in% subset(sc_obj, DataID=="LAM3_scRNA")$Barcode)

#import LAM32_Multiome
raw_counts = Read10X_h5(filename = "imported_data/old_lam_data/LAM32_Multiome_raw_feature_bc_matrix.h5", use.names = FALSE)
raw_counts = raw_counts[["Gene Expression"]]
LAM32_Multiome = CreateSeuratObject(counts = raw_counts)
LAM32_Multiome = LAM32_Multiome[,sc_obj[,WhichCells(sc_obj, expression=DataID=="LAM32_Multiome")]$Barcode]
all(subset(sc_obj, DataID=="LAM32_Multiome")$Barcode %in% colnames(LAM32_Multiome))
all(colnames(LAM32_Multiome) %in% subset(sc_obj, DataID=="LAM32_Multiome")$Barcode)

#import LAM44_Multiome
raw_counts = Read10X_h5(filename = "imported_data/old_lam_data/LAM44_Multiome_raw_feature_bc_matrix.h5", use.names = FALSE)
raw_counts = raw_counts[["Gene Expression"]]
LAM44_Multiome = CreateSeuratObject(counts = raw_counts)
LAM44_Multiome = LAM44_Multiome[,sc_obj[,WhichCells(sc_obj, expression=DataID=="LAM44_Multiome")]$Barcode]
all(subset(sc_obj, DataID=="LAM44_Multiome")$Barcode %in% colnames(LAM44_Multiome))
all(colnames(LAM44_Multiome) %in% subset(sc_obj, DataID=="LAM44_Multiome")$Barcode)

#import LAM1158_scRNA
raw_counts = Read10X(data.dir = "imported_data/old_lam_data/LAM1158_scRNA", gene.column = 1)
LAM1158_scRNA = CreateSeuratObject(counts = raw_counts)
LAM1158_scRNA = LAM1158_scRNA[,sc_obj[,WhichCells(sc_obj, expression=DataID=="LAM1158_scRNA")]$Barcode]
all(subset(sc_obj, DataID=="LAM1158_scRNA")$Barcode %in% colnames(LAM1158_scRNA))
all(colnames(LAM1158_scRNA) %in% subset(sc_obj, DataID=="LAM1158_scRNA")$Barcode)

#import LAM1163_scRNA
raw_counts = Read10X(data.dir = "imported_data/old_lam_data/LAM1163_scRNA", gene.column = 1)
LAM1163_scRNA = CreateSeuratObject(counts = raw_counts)
LAM1163_scRNA = LAM1163_scRNA[,sc_obj[,WhichCells(sc_obj, expression=DataID=="LAM1163_scRNA")]$Barcode]
all(subset(sc_obj, DataID=="LAM1163_scRNA")$Barcode %in% colnames(LAM1163_scRNA))
all(colnames(LAM1163_scRNA) %in% subset(sc_obj, DataID=="LAM1163_scRNA")$Barcode)

#import LAM1164_scRNA
raw_counts = Read10X(data.dir = "imported_data/old_lam_data/LAM1164_scRNA", gene.column = 1)
LAM1164_scRNA = CreateSeuratObject(counts = raw_counts)
LAM1164_scRNA = LAM1164_scRNA[,sc_obj[,WhichCells(sc_obj, expression=DataID=="LAM1164_scRNA")]$Barcode]
all(subset(sc_obj, DataID=="LAM1164_scRNA")$Barcode %in% colnames(LAM1164_scRNA))
all(colnames(LAM1164_scRNA) %in% subset(sc_obj, DataID=="LAM1164_scRNA")$Barcode)
                                                                  
#import LAM_UPenn_Rep1_scRNA
raw_counts = Read10X_h5(filename = "imported_data/old_lam_data/LAM_UPenn_Rep1_scRNA_raw_feature_bc_matrix.h5", use.names = FALSE)
LAM_UPenn_Rep1_scRNA = CreateSeuratObject(counts = raw_counts)
LAM_UPenn_Rep1_scRNA = LAM_UPenn_Rep1_scRNA[,sc_obj[,WhichCells(sc_obj, expression=DataID=="LAM_UPenn_Rep1_scRNA")]$Barcode]
all(subset(sc_obj, DataID=="LAM_UPenn_Rep1_scRNA")$Barcode %in% colnames(LAM_UPenn_Rep1_scRNA))
all(colnames(LAM_UPenn_Rep1_scRNA) %in% subset(sc_obj, DataID=="LAM_UPenn_Rep1_scRNA")$Barcode)

#import LAM_UPenn_Rep2_scRNA
raw_counts = Read10X_h5(filename = "imported_data/old_lam_data/LAM_UPenn_Rep2_scRNA_raw_feature_bc_matrix.h5", use.names = FALSE)
LAM_UPenn_Rep2_scRNA = CreateSeuratObject(counts = raw_counts)
LAM_UPenn_Rep2_scRNA = LAM_UPenn_Rep2_scRNA[,sc_obj[,WhichCells(sc_obj, expression=DataID=="LAM_UPenn_Rep2_scRNA")]$Barcode]
all(subset(sc_obj, DataID=="LAM_UPenn_Rep2_scRNA")$Barcode %in% colnames(LAM_UPenn_Rep2_scRNA))
all(colnames(LAM_UPenn_Rep2_scRNA) %in% subset(sc_obj, DataID=="LAM_UPenn_Rep2_scRNA")$Barcode)

#import LAM18_scRNA
raw_counts = Read10X_h5(filename = "imported_data/old_lam_data/LAM18_scRNA_raw_feature_bc_matrix.h5", use.names = FALSE)
LAM18_scRNA = CreateSeuratObject(counts = raw_counts)
LAM18_scRNA = LAM18_scRNA[,sc_obj[,WhichCells(sc_obj, expression=DataID=="LAM18_scRNA")]$Barcode]
all(subset(sc_obj, DataID=="LAM18_scRNA")$Barcode %in% colnames(LAM18_scRNA))
all(colnames(LAM18_scRNA) %in% subset(sc_obj, DataID=="LAM18_scRNA")$Barcode)

#import LAM32_scRNA
raw_counts = Read10X_h5(filename = "imported_data/old_lam_data/LAM32_scRNA_raw_feature_bc_matrix.h5", use.names = FALSE)
LAM32_scRNA = CreateSeuratObject(counts = raw_counts)
LAM32_scRNA = LAM32_scRNA[,sc_obj[,WhichCells(sc_obj, expression=DataID=="LAM32_scRNA")]$Barcode]
all(subset(sc_obj, DataID=="LAM32_scRNA")$Barcode %in% colnames(LAM32_scRNA))
all(colnames(LAM32_scRNA) %in% subset(sc_obj, DataID=="LAM32_scRNA")$Barcode)

#import LAM3_Multiome
raw_counts = Read10X_h5(filename = "raw_data/LAM3_Multiome_raw_feature_bc_matrix.h5", use.names = FALSE)
raw_counts = raw_counts[["Gene Expression"]]
LAM3_Multiome = CreateSeuratObject(counts = raw_counts)
LAM3_Multiome = LAM3_Multiome[,sc_obj[,WhichCells(sc_obj, expression=DataID=="LAM3_Multiome")]$Barcode]
all(subset(sc_obj, DataID=="LAM3_Multiome")$Barcode %in% colnames(LAM3_Multiome))
all(colnames(LAM3_Multiome) %in% subset(sc_obj, DataID=="LAM3_Multiome")$Barcode)

#import LAM50_Multiome
raw_counts = Read10X_h5(filename = "raw_data/LAM50_Multiome_raw_feature_bc_matrix.h5", use.names = FALSE)
raw_counts = raw_counts[["Gene Expression"]]
LAM50_Multiome = CreateSeuratObject(counts = raw_counts)
LAM50_Multiome = LAM50_Multiome[,sc_obj[,WhichCells(sc_obj, expression=DataID=="LAM50_Multiome")]$Barcode]
all(subset(sc_obj, DataID=="LAM50_Multiome")$Barcode %in% colnames(LAM50_Multiome))
all(colnames(LAM50_Multiome) %in% subset(sc_obj, DataID=="LAM50_Multiome")$Barcode)

lam_list = c("LAM1_scRNA"=LAM1_scRNA,
             "LAM3_scRNA"=LAM3_scRNA,
             "LAM32_Multiome"=LAM32_Multiome,
             "LAM44_Multiome"=LAM44_Multiome,
             "LAM1158_scRNA"=LAM1158_scRNA,
             "LAM1163_scRNA"=LAM1163_scRNA,
             "LAM1164_scRNA"=LAM1164_scRNA,
             "LAM_UPenn_Rep1_scRNA"=LAM_UPenn_Rep1_scRNA,
             "LAM_UPenn_Rep2_scRNA"=LAM_UPenn_Rep2_scRNA,
             "LAM18_scRNA"=LAM18_scRNA,
             "LAM32_scRNA"=LAM32_scRNA,
             "LAM3_Multiome"=LAM3_Multiome,
             "LAM50_Multiome"=LAM50_Multiome)

lam_list = lapply(names(lam_list), function(x){
  cat(x)
  cat("\n")
  o = lam_list[x][[1]]
  o = RenameCells(o, new.names=paste(x, colnames(o), sep="_"))
  DefaultAssay(o) = "RNA"
  o@meta.data$pMT = sc_obj@meta.data[rownames(o@meta.data), "pMT"]
  o@meta.data$S.Score = sc_obj@meta.data[rownames(o@meta.data), "S.Score"]
  o@meta.data$G2M.Score = sc_obj@meta.data[rownames(o@meta.data), "G2M.Score"]
  
  o = NormalizeData(o)%>%
    FindVariableFeatures()%>%
    ScaleData()%>%
    SCTransform(vars.to.regress = c("pMT", "S.Score", "G2M.Score"))
  o
})

#integrate
features = SelectIntegrationFeatures(object.list = lam_list, nfeatures = 2000)
lam_list = PrepSCTIntegration(object.list = lam_list, anchor.features=features)
lam_list <- lapply(X = lam_list, FUN = RunPCA, features = features)

anchors = FindIntegrationAnchors(lam_list, anchor.features = features,
                                 reduction = "rpca",
                                 k.anchor = 15,
                                 normalization.method = "SCT")
lam_int = IntegrateData(anchors, normalization.method = "SCT")

sc_obj = RenameAssays(sc_obj, RNA="SoupX")
sc_obj = RenameAssays(sc_obj, SCT="SoupX_SCT")
sc_obj@assays[["SCT"]] = lam_int@assays[["SCT"]]
sc_obj@assays[["RNA"]] = lam_int@assays[["RNA"]]

DefaultAssay(sc_obj) = "RNA"
sc_obj = NormalizeData(sc_obj)
sc_obj = FindVariableFeatures(sc_obj)
sc_obj = ScaleData(sc_obj)

saveRDS(sc_obj, "saved_objects/lam_combined.rds")
writeLines(capture.output(sessionInfo()), "session_info/lam_combined_session_info.txt")
