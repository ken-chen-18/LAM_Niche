library(Seurat)
library(SeuratDisk)


setwd("/global/scratch/users/kenchen/LAM/5_regulon")
input_path <- "/global/scratch/users/kenchen/LAM/input/"
cellref_genes <- readRDS("/global/scratch/users/kenchen/LAM/input/files/CellRef.genes 1.rds")

lam_ctrl_final <- readRDS("../3_lam_ctrl/lam_ctrl_final.rds")
load("/global/scratch/users/kenchen/LAM/2_uterus/lam_ut_sym.RData")

lam_ctrl_diet = DietSeurat(object = lam_ctrl_final, counts = T, data = T, scale.data = F, assays = "RNA")

lam_ctrl_diet@assays$RNA@counts@Dimnames[[1]] <- cellref_genes$symbol
lam_ctrl_diet@assays$RNA@data@Dimnames[[1]] <- cellref_genes$symbol


rownames(lam_ctrl_diet@assays$RNA@meta.features) <- cellref_genes$symbol
SaveH5Seurat(lam_ctrl_diet, filename = "converted/lam_ctrl.h5Seurat")
Convert("converted/lam_ctrl.h5Seurat", dest = "h5ad")


lam_ut_diet = DietSeurat(object = lam_ut_sym, counts = T, data = T, scale.data = F, assays = "RNA")
SaveH5Seurat(lam_ut_diet, filename = "converted/lam_ut.h5Seurat")
Convert("converted/lam_ut.h5Seurat", dest = "h5ad")
