options(encoding = "UTF-8")
.libPaths( c('/global/scratch/users/kenchen/R_4.2', .libPaths()) )
library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)


setwd("/global/scratch/users/kenchen/LAM/2_uterus")
input_path <- "/global/scratch/users/kenchen/LAM/input/"

load("lam_ut.RData")
lam_lung <- readRDS("../1_lam_int/lam_lung.rds")
lam_lung_filt <- readRDS("../1_lam_int/lam_lung_filt.rds")
cellref_genes <- readRDS(paste0(input_path, "files/CellRef.genes 1.rds"))

lam_ut <- subset(lam_ut, features = rownames(lam_lung))

#########
#INTEGRATING LUNG AND UTERUS
#########

lam_list <- c(SplitObject(lam_lung, split.by="DataID"), c("LAMUT" = lam_ut))
samples <- names(lam_list)

set.seed(2023)
lam_list <- lapply(lam_list, FUN = function(o) {
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
                                 k.anchor = 15,
                                 normalization.method = "SCT")

lam_rna_int = IntegrateData(anchors, normalization.method = "SCT")
DefaultAssay(lam_rna_int) = "integrated"

lam_rna_int = RunPCA(lam_rna_int)
lam_rna_int = RunUMAP(lam_rna_int, dims = 1:50)

DefaultAssay(lam_rna_int) <- "RNA"
Idents(lam_rna_int) <- "celltype_011625"


lam_rna_int_filt <- subset(lam_rna_int, cells=c(Cells(lam_lung_filt), Cells(lam_ut)))
lam_rna_int_filt <- subset(lam_rna_int_filt, features = rownames(cellref_genes))

saveRDS(lam_rna_int, file = "lam_rna_int.rds")
saveRDS(lam_rna_int_filt, file = "lam_rna_int_filt.rds")
