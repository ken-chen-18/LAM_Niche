options(encoding = "UTF-8")
.libPaths( c('/global/scratch/users/kenchen/R_4.2', .libPaths()) )
library(Matrix)
library(Seurat)
library(dplyr)
library(readr)



setwd("/global/scratch/users/kenchen/LAM/2_uterus")
input_path <- "/global/scratch/users/kenchen/LAM/input/"

###########
#CREATING UTERINE SEURAT OBJECT
############
uterine_genes <- read.table(paste0(input_path, "files/uterine.genes.tsv"), row.names = 1)[1]
names(uterine_genes) <- "symbol"

disease_meta <- read.table(paste0(input_path, "data/uterus/disease/meta.txt"), header=T, row.names = 1)
disease_meta$id <- paste0(substring(rownames(disease_meta), 12), "-1")
disease_meta$cell_name <- paste0("LAMUT_", disease_meta$id) 

lam_d_data <- Read10X_h5(paste0(input_path, "data/h38_uterus/disease/raw_feature_bc_matrix.h5"))

data_sub <- lam_d_data[, disease_meta$id]
lam_ut_sym <- CreateSeuratObject(counts = data_sub, project = "LAMUT")
lam_ut_sym <- subset(lam_ut_sym, cells = disease_meta$id)

lam_ut_sym <- RenameCells(lam_ut_sym, "LAMUT")

################
#CELL CYCLE + MT QC
###############
lam_ut_sym[["pMT"]] <- PercentageFeatureSet(lam_ut_sym, pattern = "^MT-")

lam_ut_sym = subset(lam_ut_sym, subset = 
                      nCount_RNA > 0 & nCount_RNA < 40000 &
                      nFeature_RNA > 500 & nFeature_RNA < 7500 &
                      pMT > -1 & pMT < 20)

set.seed(2023)
lam_ut_sym = NormalizeData(lam_ut_sym)%>%
  FindVariableFeatures() %>%
  ScaleData()%>%
  CellCycleScoring(set.ident = FALSE,
                   s.features = cc.genes.updated.2019$s.genes,
                   g2m.features = cc.genes.updated.2019$g2m.genes) %>%
  SCTransform(vars.to.regress = c("pMT", "S.Score", "G2M.Score")) %>%
  RunPCA() %>%
  RunUMAP(dims=1:50)

d_meta_filt <- subset(disease_meta, cell_name %in% rownames(lam_ut_sym@meta.data))
lam_ut_sym@meta.data <- cbind(lam_ut_sym@meta.data, d_meta_filt)
lam_ut_sym$DataID <- "LAMUT"
lam_ut_sym$DataID2 <- "LAM17"
lam_ut_sym$Technology <- "snRNA"
lam_ut_sym$Tissue <- "uterus"
lam_ut_sym$Modality <- "RNA"
lam_ut_sym$Diagnosis <- "LAM"

lam_ut_meta <- lam_ut_sym@meta.data

###########
#Ensemble IDs
#########
rownames(data_sub) <- rownames(uterine_genes)
lam_ut <- CreateSeuratObject(counts = data_sub, project = "LAMUT")
lam_ut <- subset(lam_ut, cells = d_meta_filt$id)

lam_ut <- RenameCells(lam_ut, "LAMUT")
lam_ut@meta.data <- lam_ut_meta

set.seed(2023)
lam_ut = NormalizeData(lam_ut)%>%
  FindVariableFeatures() %>%
  ScaleData()%>%
  SCTransform(vars.to.regress = c("pMT", "S.Score", "G2M.Score")) %>%
  RunPCA() %>%
  RunUMAP(dims=1:50)

DimPlot(lam_ut_sym, group.by="celltype_4")+DimPlot(lam_ut, group.by="celltype_4")

DefaultAssay(lam_ut) <- "RNA"
DefaultAssay(lam_ut_sym) <- "RNA"


lam_ut$celltype_011625 <- ifelse(lam_ut$celltype_4=="LAM", "Uterine LAMCORE",
                                 ifelse(lam_ut$celltype_4=="LymEndo", "Uterine LEC",
                                        ifelse(lam_ut$celltype_4=="VasEndo", "Uterine VEC",
                                               ifelse(lam_ut$celltype_4=="Myocytes", "Uterine SMC",
                                                      paste0("Uterine ", lam_ut$celltype_4)))))

old_umap <- lam_ut@meta.data[c("UMAP1", "UMAP2")]
lam_ut[["old_umap"]] <- CreateDimReducObject(as.matrix(old_umap), key="old_umap")

save(lam_ut, file = "lam_ut.RData")
save(lam_ut_sym, file = "lam_ut_sym.RData")
