options(encoding = "UTF-8")
.libPaths( c('/global/scratch/users/kenchen/R_4.2', .libPaths()) )
library(Seurat)
library(dplyr)

setwd("/global/scratch/users/kenchen/LAM/1_lam_int")
input_path <- "/global/scratch/users/kenchen/LAM/input/"
lam_genes <- readRDS(paste0(input_path, "files/LAM.genes.common.rds"))
cellref_genes <- readRDS("/global/scratch/users/kenchen/LAM/input/files/CellRef.genes 1.rds")
samp_meta <- read.table(paste0(input_path, "files/sample_meta.csv"),sep=",", 
                        header=T)

lam_annotated <- readRDS("../cihan/lam_annotated.rds")

lam_annotated$celltype_011625 <- ifelse(lam_annotated$seurat_clusters==45, "AT2",
                                        ifelse(lam_annotated$seurat_clusters==59, "AT1",
                                               lam_annotated$celltype_072923))
lam_annotated@meta.data <- lam_annotated@meta.data %>%
  mutate(celltype_011625 = case_when(
    celltype_011625 == "SM" ~ "SMC",
    celltype_011625 == "CD4_T" ~ "CD4 T",
    celltype_011625 == "CD8_T" ~ "CD8 T",
    celltype_011625 == "Mast_Basophil" ~ "Mast/Basophil",
    celltype_011625 == "LAMCORE-Type1" ~ "LAMCORE-2",
    celltype_011625 == "LAMCORE-Type2" ~ "LAMCORE-1",
    TRUE ~ celltype_011625
  ))


lam_annotated$DonorID <- gsub("_.*", "", rownames(lam_annotated[[]]))
lam_annotated$Tissue <- "lung"
lam_annotated$Modality <- "RNA"
lam_annotated$Diagnosis <- "LAM"
lam_annotated$Technology <- gsub(".*(scRNA|Multiome).*", "\\1", rownames(lam_annotated[[]]))
lam_annotated$DataID2 <- samp_meta[match(lam_annotated$DataID, samp_meta$DataID), "DataID2"]

lam_lung <- lam_annotated
DefaultAssay(lam_lung) <- "RNA"

mesh_reclustering <- readRDS("../0_cihan/mesh_reclustering.rds")
set.seed(2023)
mesen = FindClusters(mesh_reclustering, resolution=1.6)
DimPlot(mesen, label=T)
DefaultAssay(mesen) <- "RNA"

lam_lung$celltype_011625 <- ifelse(lam_lung$celltype_011625 %in% c("LAMCORE-1", "LAMCORE-2"),
                                   "old LAMCORE", lam_lung$celltype_011625)
lam_lung$celltype_011625[WhichCells(mesen, ident=14)] <- "LAMCORE-1"
lam_lung$celltype_011625[WhichCells(mesen, ident=15)] <- "LAMCORE-2"
lam_lung$celltype_011625[WhichCells(mesen, ident=10)] <- "LAMCORE-3"
lam_lung$celltype_011625[lam_lung$celltype_011625=="old LAMCORE" & lam_lung$seurat_clusters==47] <- "SMC"
lam_lung$celltype_011625[lam_lung$celltype_011625=="old LAMCORE" & lam_lung$seurat_clusters==49] <- "Pericyte"
lam_lung$celltype_011625[lam_lung$celltype_011625=="old LAMCORE" & lam_lung$seurat_clusters==62] <- "Mesothelial"
lam_lung$Technology <- ifelse(lam_lung$Technology == "Multiome", "snRNA", "scRNA")

Idents(lam_lung) <- "celltype_011625"

mesen$celltype_011625 <- lam_lung$celltype_011625
Idents(mesen) <- "celltype_011625"

lam_lung_filt <- subset(lam_lung, features = rownames(cellref_genes))
lam_lung_filt$pMT <- PercentageFeatureSet(lam_lung_filt, pattern = "^MT-")
lam_lung_filt <- subset(lam_lung_filt, nFeature_RNA>500 & nFeature_RNA<7500 & nCount_RNA<40000 & 
                          pMT<20 & celltype_011625 != "Low-Quality")

mesen_filt <- subset(mesen, features = rownames(cellref_genes))
mesen_filt <- subset(mesen_filt, cells = Cells(lam_lung_filt))

lam_lung_nonim <- subset(lam_lung_filt, lineage_level1 != "Immune")


saveRDS(lam_lung, file="lam_lung.rds")
saveRDS(lam_lung_filt, file="lam_lung_filt.rds")
saveRDS(mesen, file="mesen.rds")
saveRDS(mesen_filt, file="mesen_filt.rds")
saveRDS(lam_lung_nonim, file="lam_lung_nonim.rds")



