setwd("X:/cakar/lam_2023/lam_analysis_2023")
rm(list=ls())

library(Seurat)
library(tidyverse)
library(harmony)
library(reticulate)
library(ShinyCell)
options(reticulate.conda_binary = "C:/Users/caksf7/AppData/Local/anaconda3/condabin/conda.bit")
use_condaenv("py8")
leidenalg = import("leidenalg")

set.seed(2023)

label_transfer_Seurat = function(sc_obj, ref, meta_to_transfer,
                                 gene_conv_list=NULL, gene_conv_col=NA,
                                 pred_scores=NA,
                                 meta_prefix=""){
  if(!is.null(gene_conv_list)){
    counts = sc_obj@assays$RNA@counts
    rownames(counts) = gene_conv_list[rownames(counts),gene_conv_col]
    tmp = CreateSeuratObject(counts = counts)
    tmp = NormalizeData(tmp) %>%
      FindVariableFeatures() %>%
      ScaleData()%>%
      SCTransform() %>%
      RunPCA()
  }else{
    counts = sc_obj@assays$RNA@counts
    tmp = CreateSeuratObject(counts = counts)
    tmp = NormalizeData(tmp) %>%
      FindVariableFeatures() %>%
      ScaleData()%>%
      SCTransform() %>%
      RunPCA()
  }
  
  for(i in 1:length(meta_to_transfer)){
    transfer_anchors = FindTransferAnchors(
      reference = ref,
      query = tmp,
      reference.reduction = "pca",
      normalization.method = "SCT",
      dims = 1:50
    )
    predictions = TransferData(
      anchorset = transfer_anchors,
      reference = ref,
      refdata = ref@meta.data[,meta_to_transfer[[i]]],
      dims = 1:50
    )
    if(!is.na(pred_scores) & pred_scores == meta_to_transfer[[i]]){
      sc_obj = AddMetaData(sc_obj, predictions)
      sc_obj@meta.data[,paste0(meta_prefix,"_",meta_to_transfer[[i]])] = predictions$predicted.id
    }else{
      sc_obj@meta.data[,paste0(meta_prefix,"_",meta_to_transfer[[i]])] = predictions$predicted.id
    }
  }
  sc_obj
}
lam_genes = readRDS("aux_data/LAM.genes.rds")
cellref = readRDS("imported_data/cellref.rds")
cellref_seed = readRDS("imported_data/cellref_seed.rds")

sc_obj = readRDS("saved_objects/lam_integrate_rpca.rds")


################################################################################
################################################################################
######################### Mesenchymal Cells Reclustering #######################
################################################################################
################################################################################

mesh = subset(sc_obj, seurat_clusters %in% c(62,22,37,39,47,18,49))
mesh$glb_cls = mesh$seurat_clusters
DefaultAssay(mesh) = "SCT"
mesh = RunPCA(mesh, features = rownames(mesh))
mesh = RunHarmony(mesh, group.by.vars = c("DataID","Technology"), project.dim = FALSE, assay.use = "SCT")
mesh = FindNeighbors(mesh, dims = 1:30, reduction = "harmony")
mesh = FindClusters(mesh, method = "igraph", algorithm = 4, resolution=0.8)
mesh = RunUMAP(mesh, reduction = "harmony", dims=1:30)
DimPlot(mesh, label = T)
scConf = createConfig(mesh, maxLevels = 100)
makeShinyApp(mesh, scConf, gene.mapping = TRUE, shiny.title = "LAM Label Diagnostics",
             gex.assay = "SCT", shiny.dir = "shiny_outs/mesh_reclustering")


saveRDS(mesh, "saved_objects/mesh_reclustering.rds")
mesh = readRDS("saved_objects/mesh_reclustering.rds")

################################################################################
################################################################################
################### Cluster 58/61 Reclustering #################################
################################################################################
################################################################################

clusters58_61 = subset(sc_obj, seurat_clusters %in% c(28,32, 61,4,58))
clusters58_61$glb_cls = clusters58_61$seurat_clusters
clusters58_61 = label_transfer_Seurat(clusters58_61, subset(cellref_seed, celltype_level3 %in% c("CAP1","SVEC", "VEC")), meta_prefix="cellref_seed",
                          meta_to_transfer = c("celltype_level3", "lineage_level1"),
                          pred_scores = "celltype_level3")
clusters58_61 = label_transfer_Seurat(clusters58_61, subset(cellref, celltype_level3 %in% c("CAP1","SVEC", "VEC")), meta_prefix="cellref",
                          meta_to_transfer = c("celltype_level3", "lineage_level1"),
                          pred_scores = "celltype_level3", gene_conv_list = lam_genes, gene_conv_col = "symbol")

DefaultAssay(clusters58_61) = "SCT"
clusters58_61 = RunPCA(clusters58_61, features = rownames(clusters58_61))
clusters58_61 = RunHarmony(clusters58_61, group.by.vars = c("DataID","Technology"), project.dim = FALSE, assay.use = "SCT")
clusters58_61 = FindNeighbors(clusters58_61, dims = 1:30, reduction = "harmony")
clusters58_61 = FindClusters(clusters58_61, method = "igraph", algorithm = 4, resolution=0.8)
clusters58_61 = RunUMAP(clusters58_61, reduction = "harmony", dims=1:30)
DimPlot(clusters58_61, label = T)
scConf = createConfig(clusters58_61, maxLevels = 100)
makeShinyApp(clusters58_61, scConf, gene.mapping = TRUE, shiny.title = "LAM Label Diagnostics",
             gex.assay = "SCT", shiny.dir = "shiny_outs/clusters58_61_reclustering")


saveRDS(clusters58_61, "saved_objects/cluster58_61_reclustering.rds")

################################################################################
################################################################################
################### Dentric Cells Reclustering #################################
################################################################################
################################################################################

dc_clusters = subset(sc_obj, seurat_clusters %in% c(30, 38, 74))
dc_clusters$glb_cls = dc_clusters$seurat_clusters
dc_clusters = label_transfer_Seurat(dc_clusters, subset(cellref_seed, celltype_level3 %in% c("cDC1", "cDC2", "maDC", "pDC")), meta_prefix="cellref_seed",
                          meta_to_transfer = c("celltype_level3", "lineage_level1"),
                          pred_scores = "celltype_level3")
dc_clusters = label_transfer_Seurat(dc_clusters, subset(cellref, celltype_level3 %in% c("cDC1", "cDC2", "maDC", "pDC")), meta_prefix="cellref",
                          meta_to_transfer = c("celltype_level3", "lineage_level1"),
                          pred_scores = "celltype_level3", gene_conv_list = lam_genes, gene_conv_col = "symbol")

DefaultAssay(dc_clusters) = "SCT"
dc_clusters = RunPCA(dc_clusters, features = rownames(dc_clusters))
dc_clusters = RunHarmony(dc_clusters, group.by.vars = c("DataID","Technology"), project.dim = FALSE, assay.use = "SCT")
dc_clusters = FindNeighbors(dc_clusters, dims = 1:30, reduction = "harmony")
dc_clusters = FindClusters(dc_clusters, method = "igraph", algorithm = 4, resolution=0.8)
dc_clusters = RunUMAP(dc_clusters, reduction = "harmony", dims=1:30)
DimPlot(dc_clusters, label = T)
scConf = createConfig(dc_clusters, maxLevels = 100)
makeShinyApp(dc_clusters, scConf, gene.mapping = TRUE, shiny.title = "LAM Label Diagnostics",
             gex.assay = "SCT", shiny.dir = "shiny_outs/dc_reclustering")

saveRDS(dc_clusters, "saved_objects/dc_reclustering.rds")

writeLines(capture.output(sessionInfo()), "session_info/reclustering_session_info.txt")
