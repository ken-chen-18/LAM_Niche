setwd("X:/cakar/2023_LAM")
rm(list=ls())

library(Seurat)
library(tidyverse)
library(ShinyCell)
library(stringi)

#import output object from 5_lam_combine.R
sc_obj = readRDS("saved_objects/lam_combined.rds")

#import cluster celltype mapping for global clusters
celltypes = read.csv("saved_info/lam_cluster_celltypes.csv")

#import cluster celltype mapping for cluster 58 - 61 reclustering
celltypes_endo = read.csv("saved_info/cluster_58_61_celltypes.csv")

#import cluster celltype mapping for DC cells
celltypes_dc = read.csv("saved_info/dc_celltypes.csv")

#import celltype lineage mapping
celltypes_lineages = read.csv("aux_data/celltype_lineages.csv", row.names = 1)
celltypes_lineages$Celltype = gsub("/", "_", celltypes_lineages$Celltype)
rownames(celltypes_lineages) = celltypes_lineages$Celltype

#add lineage column to cluster celltype mapping and complete missing lineages
celltypes = mutate(celltypes, Lineage=celltypes_lineages[CellType, "Lineage"])
celltypes[celltypes$CellType=="AT1_AT2", "Lineage"] = "Epithelial"
celltypes[celltypes$CellType=="NKT", "Lineage"] = "Immune"
celltypes[celltypes$CellType=="DC", "Lineage"] = "Immune"
celltypes[celltypes$CellType=="LAMCORE", "Lineage"] = "Mesenchymal"
celltypes[celltypes$CellType=="Low-Quality", "Lineage"] = "Unknown"
rownames(celltypes) = celltypes$Cluster

#add celltypes and lineages to the meta.data
sc_obj@meta.data = mutate(sc_obj@meta.data, celltype_level3 = celltypes[seurat_clusters, "CellType"])
sc_obj@meta.data = mutate(sc_obj@meta.data, lineage_level1 = celltypes[seurat_clusters, "Lineage"])

#import mesenchymal reclustering object
mesh = readRDS("saved_objects/mesh_reclustering.rds")

#relabel LAMCORE cells (Cluster 47) using mesenchymal reclustering object 
sc_obj@meta.data[WhichCells(mesh, expression = seurat_clusters==6),"celltype_level3"] = "LAMCORE-Type2"
sc_obj@meta.data[WhichCells(mesh, expression = seurat_clusters==10),"celltype_level3"] = "LAMCORE-Type1"
sc_obj@meta.data[sc_obj@meta.data$seurat_clusters==47 & !(sc_obj@meta.data$celltype_level3 %in% c("LAMCORE-Type1", "LAMCORE-Type2")),"celltype_level3"] = "SM"

#import cluster58_61 object to relabel Clusters 58 - 61
clusters58_61 = readRDS("saved_objects/cluster58_61_reclustering.rds")
rownames(celltypes_endo) = celltypes_endo$Cluster
clusters58_61@meta.data = mutate(clusters58_61@meta.data, celltype_level3 = celltypes_endo[seurat_clusters, "CellType"])
clusters58_61 = subset(clusters58_61, glb_cls %in% c(58,61))

# update cell labels of 58 and 61 using cluster58_61 reclustering object
sc_obj@meta.data[rownames(clusters58_61@meta.data), "celltype_level3"] = clusters58_61$celltype_level3

#import dc clusters object to relabel dc cells
dc_clusters = readRDS("saved_objects/dc_reclustering.rds")
rownames(celltypes_dc) = celltypes_dc$Cluster
dc_clusters@meta.data = mutate(dc_clusters@meta.data, celltype_level3 = celltypes_dc[seurat_clusters, "CellType"])

# update cell labels of dc cells using dc_clusters reclustering object
sc_obj@meta.data[rownames(dc_clusters@meta.data), "celltype_level3"] = dc_clusters$celltype_level3

sc_obj$celltype_072923 = sc_obj$celltype_level3

#save the meta which contains legacy columns
write.csv(sc_obj@meta.data,"saved_info/old_meta.csv", quote=FALSE)

#remove extra meta data columns
sc_obj@meta.data = sc_obj@meta.data[c("orig.ident","Phase","pMT", "S.Score", "G2M.Score",
                                      "nFeature_RNA", "nCount_RNA", 
                                      "cellref_seed_celltype_level3","cellref_seed_lineage_level1",
                                      "cellref_celltype_level3","cellref_lineage_level1",
                                      "celltype_072923","lineage_level1", "seurat_clusters", "integrated_snn_res.4.2",
                                      "DataID", "DonorID", "Barcode", "prediction.score.max")]
sc_obj@meta.data$sex = "Female"
sample_info = read.csv("saved_info/sample_meta.csv", row.names = 1)
sc_obj@meta.data = mutate(sc_obj@meta.data, age = sample_info[DataID, "Age"])
write.csv(sc_obj@meta.data,"saved_info/lam_meta.csv", quote=FALSE)
saveRDS(sc_obj, "saved_objects/lam_annotated.rds")
writeLines(capture.output(sessionInfo()), "session_info/lam_annotations_session_info.txt")
