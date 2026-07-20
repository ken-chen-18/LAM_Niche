setwd("Y:/cakar/2023_LAM")
rm(list = ls())
library(Seurat)
library(ggplot2)
library(tidyverse)
library(ShinyCell)
library(reticulate)
options(reticulate.conda_binary = "C:/Users/caksf7/AppData/Local/anaconda3/condabin/conda.bit")

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

use_condaenv("py8")
leidenalg = import("leidenalg")

set.seed(2023)

lam_genes = readRDS("aux_data/LAM.genes.rds")
cellref = readRDS("imported_data/cellref.rds")
cellref_seed = readRDS("imported_data/cellref_seed.rds")

lam3 = readRDS("saved_objects/LAM3_Multiome_qc_pipe.rds")
lam3@meta.data$Barcode = str_extract(rownames(lam3@meta.data),"[:upper:]{16}-[:digit:]+")

lam50 = readRDS("saved_objects/LAM50_Multiome_qc_pipe.rds")
lam50@meta.data$Barcode = str_extract(rownames(lam50@meta.data),"[:upper:]{16}-[:digit:]+")

lam_combined = readRDS("imported_data/LAM.combined.rds")
lam_combined$DataID = gsub("LAM44_multiome", "LAM44_Multiome", lam_combined$DataID)
lam_combined@meta.data[,"Barcode"] = str_extract(rownames(lam_combined@meta.data),"[:upper:]{16}-[:digit:]+")

lam_list = c(SplitObject(lam_combined, split.by = "DataID"),
             c("LAM3_Multiome" = lam3, "LAM50_Multiome" = lam50))


lam_list <- lapply(X = names(lam_list), FUN = function(x) {
  o = lam_list[x][[1]]
  o = RenameCells(o, new.names = paste(o$DataID, o$Barcode, sep = "_"))
  cat(paste("Processing ", x,":\n"))
  if(x=="LAM3_Multiome" | x=="LAM50_Multiome"){
    o = CreateSeuratObject(counts = o@assays[["SoupX"]]@counts, meta.data = o@meta.data)
  }else{
    cat("Using RNA assay\n")
    o = CreateSeuratObject(counts = o@assays[["RNA"]]@counts, meta.data = o@meta.data)
  }
  o = NormalizeData(o) %>%
    FindVariableFeatures() %>%
    ScaleData()%>%
    SCTransform(vars.to.regress = c("pMT", "S.Score", "G2M.Score"))
  
  o = label_transfer_Seurat(o, cellref_seed, meta_prefix="cellref_seed",
                               meta_to_transfer = c("celltype_level3", "lineage_level1"),
                               pred_scores = "celltype_level3")
  o = label_transfer_Seurat(o, cellref, meta_prefix="cellref",
                               meta_to_transfer = c("celltype_level3", "lineage_level1"),
                               pred_scores = "celltype_level3", gene_conv_list = lam_genes, gene_conv_col = "symbol")
  
  o
})
features = SelectIntegrationFeatures(object.list = lam_list, nfeatures = 2000)
lam_list = PrepSCTIntegration(object.list = lam_list, anchor.features=features)
lam_list <- lapply(X = lam_list, FUN = RunPCA, features = features)

anchors = FindIntegrationAnchors(lam_list, anchor.features = features,
                                 reduction = "rpca",
                                 k.anchor = 15,
                                 normalization.method = "SCT")
lam_int = IntegrateData(anchors, normalization.method = "SCT")
lam_int = RunPCA(lam_int)
lam_int = RunUMAP(lam_int, dims = 1:50)
lam_int = FindNeighbors(lam_int, dims = 1:50)
lam_int = FindClusters(lam_int, method = "igraph", algorithm = 4, resolution=4.2)

lam_int$Technology = str_extract(lam_int$DataID, "Multiome")
lam_int$Technology = replace_na(lam_int$Technology, "scRNA")

saveRDS(lam_int, "saved_objects/lam_integrate_rpca.rds")
writeLines(capture.output(sessionInfo()), "session_info/integrate_rpca_session_info.txt")
