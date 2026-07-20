library(Seurat)
library(dplyr)
library(ggplot2)
library(SoupX)

setwd("X:/cakar/lam_2023/lam_analysis_2023")
set.seed(2023)

data = read.csv("aux_data/run_parameters.csv")
for(i in 1:dim(data)[1]){
  data_id = data[i, "data_id"]
  raw_counts = Read10X_h5(filename = paste0('raw_data/', data_id, "_raw_feature_bc_matrix.h5"), use.names = FALSE)
  sc_counts = Read10X_h5(filename = paste0('raw_data/', data_id, "_filtered_feature_bc_matrix.h5"))
  raw_counts = raw_counts[["Gene Expression"]]
  sc_counts = sc_counts[["Gene Expression"]]
  
  sc_obj = CreateSeuratObject(counts = sc_counts, project = "lam_integrated_analysis", min.cells = 3, min.features = 200)
  sc_obj[['pMT']] = PercentageFeatureSet(sc_obj, pattern = "^MT-")
  
  sc_obj = ScaleData(sc_obj)%>%
    NormalizeData() %>%
    FindVariableFeatures()%>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims=1:50, verbose = FALSE)%>%
    FindNeighbors(reduction = "pca", dims = 1:50, verbose = FALSE) %>%
    FindClusters(verbose = FALSE)
  sc_obj = CellCycleScoring(sc_obj, set.ident = FALSE,
                          s.features = cc.genes.updated.2019$s.genes,
                          g2m.features = cc.genes.updated.2019$g2m.genes)
  
  
  sc_counts_ens = Read10X_h5(filename = paste0('raw_data/', data_id, "_filtered_feature_bc_matrix.h5"), use.names = FALSE)
  if(data[i,"Technology"]=="Multiome"){
    sc_counts_ens = sc_counts_ens[["Gene Expression"]]
  }
  sc_obj_ens = CreateSeuratObject(counts = sc_counts_ens, project = "lam_integrated_analysis", min.cells = 3, min.features = 200)
  sc_obj_ens@meta.data = sc_obj@meta.data
  sc_obj = sc_obj_ens
  sc_obj = ScaleData(sc_obj)%>%
    NormalizeData() %>%
    FindVariableFeatures()%>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims=1:50, verbose = FALSE)%>%
    FindNeighbors(reduction = "pca", dims = 1:50, verbose = FALSE) %>%
    FindClusters(verbose = FALSE)
  
  raw_counts = raw_counts[, which(Matrix::colSums(raw_counts)>0)]
  selected_counts = raw_counts[, colnames(sc_obj)]
  
  sc = SoupChannel(tod=raw_counts, toc=selected_counts)
  sc = setClusters(sc, setNames(sc_obj@meta.data[, "RNA_snn_res.0.8"], rownames(sc_obj@meta.data)))
  sc = setDR(sc, DR=sc_obj@reductions$umap@cell.embeddings[colnames(sc$toc), ], reductName="umap")
  
  if(data[i,"soupX_cont"]=="auto"){
    sc = autoEstCont(sc)
  } else{
    sc = setContaminationFraction(sc, data[i,"soupX_cont"])
  }
  out = adjustCounts(sc)
  out = round(out)
  
  tmp = CreateSeuratObject(counts=out, meta.data=sc_obj@meta.data)
  
  sc_obj@assays$SoupX = tmp@assays$RNA
  
  sc_obj$orig.ident = data_id
  sc_obj$DataID = data_id
  sc_obj$condition = data[i, "Condition"]
  sc_obj$Technology = data[i,"Technology"]
  
  sc_obj = subset(sc_obj, subset = 
                    nCount_RNA > data[i,"nCount_RNA_lower"] & nCount_RNA < data[i,"nCount_RNA_upper"] &
                    nFeature_RNA > data[i,"nFeature_RNA_lower"] & nFeature_RNA < data[i,"nFeature_RNA_upper"] &
                    pMT > data[i,"pMT_lower"] & pMT < data[i,"pMT_upper"])
  
  doublet_pred = read.csv(paste0("saved_info/",data_id,'_doublet_scores_predictions.csv'))
  rownames(doublet_pred) = doublet_pred$X
  doublet_pred = subset(doublet_pred, select = -c(X))
  if (data[i,"doublet_tresh"] == "auto"){
    sc_obj$doub_tres = ifelse(doublet_pred[colnames(sc_obj),]$predicted_doublet=="True", TRUE, FALSE)
  }else{
    sc_obj$doub_tres = ifelse(doublet_pred[colnames(sc_obj),]$doublet_score>=data[i,"doublet_tresh"], TRUE, FALSE)
  }
  sc_obj = subset(sc_obj, doub_tres == FALSE)
  sc_obj = RenameCells(sc_obj, add.cell.id = data_id)
  saveRDS(sc_obj, paste0("pits/",data_id, "_qc_pipe_last.rds"))
}
writeLines(capture.output(sessionInfo()), "session_info/qc_pipe_session_info.txt")