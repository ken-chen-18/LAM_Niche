options(encoding = "UTF-8")
.libPaths( c('/global/scratch/users/kenchen/R_4.2', .libPaths()) )
library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(harmony)
library(monocle3)
library(SeuratWrappers)


setwd("/global/scratch/users/kenchen/LAM/2_uterus")
input_path <- "/global/scratch/users/kenchen/LAM/input/"
cellref_genes <- readRDS("/global/scratch/users/kenchen/LAM/input/files/CellRef.genes 1.rds")

lam_rna_int_filt <- readRDS("lam_rna_int_filt.rds")



#########
#SUBCLUSTERING
#########

lamcore_only <- subset(lam_rna_int_filt, celltype_011625 %in% c("Uterine LAMCORE", "LAMCORE-1", "LAMCORE-2", "LAMCORE-3"))

Idents(lamcore_only) <- "celltype_011625"
DefaultAssay(lamcore_only)= "SCT"
set.seed(2023)
lamcore_only <- RunPCA(lamcore_only, features=rownames(lamcore_only), reduction.name="lam_only_pca")
lamcore_only = RunHarmony(lamcore_only, group.by.vars = c("DataID","Technology"), 
                      project.dim = FALSE, assay.use = "SCT", reduction="lam_only_pca")
set.seed(2023)
lamcore_only = RunUMAP(lamcore_only, reduction = "harmony", dims=1:30, reduction.name="lam_only_umap", n.neighbors = 30)

DefaultAssay(lamcore_only)= "RNA"
DimPlot(lamcore_only, reduction = "lam_only_umap", group.by = "celltype_090523", label = F)


save(lamcore_only, file="lamcore_only.RData")




#Mesen + Uterine LAMCORE
set.seed(2023)
a <- subset(lung_ut_int, subset = DataID == "uterus" | Cells(lung_ut_int) %in% Cells(lam_rna_int_filt))
a <- subset(a, features=rownames(lam_rna_int_filt))
mesen_ut_lamcore <- subset(a, celltype_090523 == "Uterine LAMCORE" | lineage_level1=="Mesenchymal")

Idents(mesen_ut_lamcore) <- "celltype_011625"
DefaultAssay(mesen_ut_lamcore)= "SCT"
mesen_ut_lamcore <- RunPCA(mesen_ut_lamcore, features=rownames(mesen_ut_lamcore), reduction.name="lam_mesh_pca")
mesen_ut_lamcore = RunHarmony(mesen_ut_lamcore, group.by.vars = c("DataID","Technology"), 
                      project.dim = FALSE, assay.use = "SCT", reduction="lam_mesh_pca")
mesen_ut_lamcore = RunUMAP(mesen_ut_lamcore, reduction = "harmony", dims=1:30, reduction.name="lam_mesh_umap", n.neighbors = 15)
DimPlot(mesen_ut_lamcore, reduction = "lam_mesh_umap", group.by = "celltype_090523", label = T, repel = T, pt.size=0.5)



DefaultAssay(mesen_ut_lamcore) <- "RNA"

save(mesen_ut_lamcore, file = "mesen_ut_lamcore.RData")


#MONOCLE3
mn <- lamcore_only
mn@reductions$umap <- mn@reductions$lam_only_umap
cds = as.cell_data_set(mn, assay = "RNA")
set.seed(2023)
cds <- cluster_cells(cds, k=100, cluster_method = "louvain")
plot_cells(cds, color_cells_by = "cluster", cell_size=1)
cds <- learn_graph(cds)
a <- cds@clusters$UMAP$clusters
cluster_cells <- which(a %in% c(2, 3, 5, 6))
cds <- order_cells(cds)

plot_cells(cds,
           color_cells_by = "cluster",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)
plot_cells(cds, color_cells_by = "pseudotime", graph_label_size = 0, cell_size = 1,
           trajectory_graph_color = "black",
           trajectory_graph_segment_size = 1.5,
           label_branch_points = F, label_leaves=T, label_roots=T)
pseudotime <- as.Seurat(cds, assay = NULL)
FeaturePlot(pseudotime, "monocle3_pseudotime")|DimPlot(pseudotime, group.by="celltype_011625")
save(pseudotime, file="pseudotime.RData")
save(cds, file="cds.RData")
