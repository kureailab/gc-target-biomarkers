# author: hongc2@ccf.org
# input: cellranger output directory
# objective: to compare groups of scRNA
# ref:
# https://satijalab.org/seurat/v3.0/immune_alignment.html
# https://www.biorxiv.org/content/biorxiv/early/2018/11/02/460147.full.pdf

library(Seurat)
library(data.table)
library(openxlsx)
library(reshape2)
library(dplyr)
library(reticulate)
library(patchwork)
library(cowplot)
library(grid)
library(gridExtra)
library(parallel)
library(DoubletFinder)
library(R.utils)
library(ggplot2)
#library(DEsingle)
library(BiocParallel)
library(matrixStats)
library(pheatmap)
library(SingleCellExperiment)
# library(SingleR)
library(future)
library(harmony) #install_github("immunogenomics/harmony")
library(SeuratWrappers) #remotes::install_github('satijalab/seurat-wrappers')
# library(monocle3) #devtools::install_github('cole-trapnell-lab/leidenbase');devtools::install_github('cole-trapnell-lab/monocle3');https://cole-trapnell-lab.github.io/monocle3/docs/installation/
library(scCustomize)
library(scSorter)
# library(UCell)
library(gtools)
library(CellChat)
library(copykat)

source(file.path(Sys.getenv('R_UTIL'),'lib_apps.R'))
# use_condaenv(condaenv="scell",conda="/media/sammy/apps/miniconda3/bin/conda")

#' Convert 10x cellranger count output to seurat3 object
#' more detail
#'
#'    Changjin Hong, hongc2@ccf.org

#'    ref:
#'
#' @param sample_dt A data frame or data.table in the column names (sample,cellranger_count_outd,sgroup). For example,
#'  D1_1, /home/blahblah/1907UNHS-0026_cellranger_count/D1_1, pre_infusion
#' @param min_cell The minimum number of cells to include.
#' @param min_features The min number of features to include.
#' @param ncpu The number of cpus to utilize
#' @param debug 0 (no debug), 1 (debug; test with the first sample only)
#' @return seurat3 object in a list.
#'
sc10x_to_seurat3 <- function(sample_dt,min_cell=3,min_features=200,ncpu=2,debug=0,mt_pct=10.) {
  
  #make sure that sample_dt has 3 columns
  stopifnot(dim(sample_dt)[2]==3)
  if (debug==1){browser()}
  S <- dim(sample_dt)[1]
  message(sprintf("total number of samples [%d]",S))
  
  if (debug==1) {S <- 1}
  
  seurat_raws <- mclapply(1:S,function(i) {
    if (debug==1){browser()}
    sample <- sample_dt[i,sample]
    crcount_out_fpath <- sample_dt[i,cellranger_count_outd]
    message(sample)
    
    bc_matrix_dir <- sprintf("%s/outs/filtered_feature_bc_matrix",crcount_out_fpath)
    
    message(sprintf('start to read 10x count at [%s]',bc_matrix_dir))
    
    sc.data <- Read10X(data.dir = bc_matrix_dir)
    
    sc_raw <- CreateSeuratObject(counts = sc.data,
                                 project = sample,
                                 min.cells = min_cell,
                                 min.features = min_features)
    
    sc_raw[["percent.mt"]] <- PercentageFeatureSet(sc_raw, pattern = "^MT-")
    sc_raw$sgroup <- sample_dt[i,sgroup]
    if (debug==1){browser()}
    sc_raw <- subset(sc_raw,subset = percent.mt < 10.)
    
    message(sprintf("Done[%s]",sample))
    return(sc_raw)
  },mc.cores = ncpu) #,mc.cores = ncpu
  
  names(seurat_raws) <- sample_dt$sample
  if (debug==1){browser()}
  return(seurat_raws)
}



simple_qc_seu <- function(seu,
                          max_mt_pct=15,
                          max_hb_pct=0.,
                          min_rb_pct=0.,
                          min.features=200,
                          max.features=0,
                          min.counts=0,
                          max.counts=0,
                          min.cells=3,
                          debug2=0) {
  
  if (debug2==1) {browser()}
  mtx.dim = dim(seu)
  
  stopifnot('RNA' %in% names(seu@assays) | 'Spatial' %in% names(seu@assays))
  nCount_col=get_meta_colname(seu,col_pat="nCount_")
  nFeature_col=get_meta_colname(seu,col_pat="nFeature_")
  
  message(sprintf("before QC: %d x %d", mtx.dim[1],mtx.dim[2]))
  
  cell_ids <- colnames(seu)
  if (max.features==0){max.features=Inf}
  if (max.counts==0){max.counts=Inf}
  
  celloi <- cell_ids[(seu@meta.data[[nCount_col]] >= min.counts & 
                        seu@meta.data[[nCount_col]] < max.counts) &
                       (seu@meta.data[[nFeature_col]] >= min.features & 
                          seu@meta.data[[nFeature_col]] < max.features)]
  
  genes = rownames(seu)
  genes_oi <- genes[Matrix::rowSums(GetAssayData(seu, slot="counts")) >= min.cells]
  
  seu <- subset(x=seu, cells=celloi, features=genes_oi)
  
  meta_cols = colnames(seu@meta.data)
  
  if (max_mt_pct>0. & ('percent.mt' %in% meta_cols)) {
    cell_ids <- colnames(seu)
    message(sprintf("applying max_mt_pct[%g]",max_mt_pct))
    celloi <- cell_ids[seu@meta.data$percent.mt<=max_mt_pct]
    seu <- subset(x=seu, cells=celloi)
  }
  
  if (max_hb_pct>0. & ('percent.hb' %in% meta_cols)) {
    cell_ids <- colnames(seu)
    message(sprintf("applying max_hb_pct[%g]",max_hb_pct))
    celloi <- cell_ids[seu@meta.data$percent.hb<=max_hb_pct]
    seu <- subset(x=seu, cells=celloi)
  }
  
  if (min_rb_pct>0. & ('percent.ribo' %in% meta_cols)) {
    cell_ids <- colnames(seu)
    message(sprintf("applying min_rb_pct[%g]",min_rb_pct))
    celloi <- cell_ids[seu@meta.data$percent.ribo>=min_rb_pct]
    seu <- subset(x=seu, cells=celloi)
  }
  mtx.dim = dim(seu)
  message(sprintf("after QC: %d x %d", mtx.dim[1],mtx.dim[2]))
  seu
}

get_k_most_distinct_genes_per_cluster <- function(sc3.markers.dt,
                                                  top_k_pos=10,
                                                  max_p_val_adj=0.001) {
  
  sc3.markers.dt2 <- split(sc3.markers.dt,by="cluster")
  
  top_most_vars <- lapply(sc3.markers.dt2,function(sc3m) {
    sc3m <- sc3m[p_val_adj<max_p_val_adj,]
    if (dim(sc3m)[1]<top_k_pos){
      sc3m[order(-avg_logFC)]
    } else {
      sc3m[order(-avg_logFC)][1:top_k_pos,]
    }
  })
  
  top_vars <- rbindlist(top_most_vars)
  
  #
  # if (top_k_neg>0) {
  # 	top_most_vars <- lapply(sc3.markers.dt2,function(sc3m) {
  # 		sc3m <- sc3m[p_val_adj<max_p_val_adj,]
  # 		if (dim(sc3m)[1]<top_k_pos){
  # 			sc3m[order(avg_logFC)]
  # 		} else {
  # 			sc3m[order(avg_logFC)][1:top_k_pos,]
  # 		}
  # 	})
  # 	top_neg_vars <- rbindlist(top_most_vars)
  # 	top_vars <- distinct(rbindlist(list(top_vars,top_neg_vars)))
  # }
  
  
  return(top_vars)
}

find_clusters <- function(seu,genes=NA){
  seu <- NormalizeData(seu)
  if (invalid(genes)){
    seu <- FindVariableFeatures(seu, selection.method = "vst")
  } else {
    VariableFeatures(seu) <- genes
  }
  seu <- ScaleData(seu)
  seu <- RunPCA(seu, 
                npcs = find_max_npcs(seu),
                features = VariableFeatures(object = seu))
  seu <- FindNeighbors(seu, dims = 1:20)
  seu <- FindClusters(seu, resolution = 0.5)
  seu <- RunUMAP(seu, dims = 1:20)
  seu
}

solo_analysis <- function(wkd,sc3,sname,max_cluster=12,vst_feature=2000) {
  
  fig_prefix <- "01_feature_readCnt_mt"
  pdf_file <- file.path(wkd,sprintf('%s_%s.pdf',sname,fig_prefix))
  message(pdf_file)
  pdf(pdf_file)
  p<-VlnPlot(sc3, pt.size=0.1, features = c(get_, get_meta_colname(sc3,col_pat = "nCount_"), "percent.mt"), ncol = 3)
  plot(p)
  dev.off()
  
  # ============
  
  fig_prefix <- "02_readCnt_vs_mt_feature"
  pdf_file <- file.path(wkd,sprintf('%s_%s.pdf',sname,fig_prefix))
  message(pdf_file)
  pdf(pdf_file)
  plot1 <- FeatureScatter(sc3, feature1 = get_meta_colname(sc3,col_pat = "nCount_"), feature2 = "percent.mt")
  plot2 <- FeatureScatter(sc3, feature1 = get_meta_colname(sc3,col_pat = "nCount_"), feature2 = "nFeature_RNA")
  p<-CombinePlots(plots = list(plot1, plot2))
  plot(p)
  dev.off()
  
  fig_prefix <- "03_readCnt_features_hist"
  pdf_file <- file.path(wkd,sprintf('%s_%s.pdf',sname,fig_prefix))
  message(pdf_file)
  pdf(pdf_file)
  plot1 <- ggplot(sc3@meta.data,aes(x=nCount_RNA))+geom_histogram(color="black",fill="white")
  plot2 <- ggplot(sc3@meta.data,aes(x=nFeature_RNA))+geom_histogram(color="black",fill="white")
  p <- CombinePlots(plots = list(plot1,plot2))
  plot(p)
  dev.off()
  
  # ..............................................
  # find variable features if necessary. Otherwise, use default
  sc3 <- FindVariableFeatures(sc3, selection.method = "vst", nfeatures = vst_feature)
  
  
  # Scaling the data across cells before dim reduction
  all.genes <- rownames(sc3)
  sc3 <- ScaleData(sc3, features = all.genes)
  
  # Perform PCA
  num_cells <- dim(sc3@assays$RNA@data)[2]
  if (num_cells < 200) {
    npcs <- floor(num_cells/4)
  } else {
    npcs <- 50
  }
  
  sc3 <- RunPCA(sc3, 
                features = VariableFeatures(object = sc3), 
                npcs = min(npcs,find_max_npcs(sc3)))
  
  # ..............................................
  VizDimLoadings(sc3, dims = 1:2, reduction = "pca")
  
  fig_prefix <- "04_pca"
  pdf_file <- file.path(wkd,sprintf('%s_%s.pdf',sname,fig_prefix))
  message(pdf_file)
  pdf(pdf_file)
  p <- DimPlot(sc3, reduction = "pca")
  plot(p)
  dev.off()
  
  fig_prefix <- "05_pca_heatmap"
  pdf_file <- file.path(wkd,sprintf('%s_%s.pdf',sname,fig_prefix))
  message(pdf_file)
  pdf(pdf_file)
  DimHeatmap(sc3, dims = 1:6, cells = 500, balanced = TRUE)
  dev.off()
  
  # ..............................................
  
  if (max_cluster==0) {
    # Determine the ‘dimensionality’ of the dataset to cluster cells
    sc3 <- JackStraw(sc3, num.replicate = 100)
    sc3 <- ScoreJackStraw(sc3, dims = 1:20)
    L <- length(sc3@reductions$pca@stdev)
    D <- which(sc3@reductions$pca@stdev < sc3@reductions$pca@stdev[L]*1.2)[1]
  } else {
    # ..............................................
    D <- max_cluster
  }
  
  if (D > npcs) {
    D <- npcs
  }
  
  message(sprintf("the number of max clusters [%d]",D))
  #Cluster the cells
  sc3 <- FindNeighbors(sc3, dims = 1:D) #KNN-graph
  sc3 <- FindClusters(sc3, resolution = 0.5) #Louvain algorithm (granularity: 0.4 ~ 1.2 for 3K cells)
  
  # Non-linear dim reduction for visualization
  # for umap installation
  # reticulate::py_install(packages ='umap-learn')
  
  sc3 <- RunUMAP(sc3, dims = 1:D, umap.method = "umap-learn", metric = "correlation")
  
  fig_prefix <- "06_umap"
  pdf_file <- file.path(wkd,sprintf('%s_%s.pdf',sname,fig_prefix))
  message(pdf_file)
  pdf(pdf_file)
  p <- DimPlot(sc3,reduction="umap")
  plot(p)
  dev.off()
  
  # find markers for every cluster compared to all remaining cells, report only the positive ones
  sc3.markers <- FindAllMarkers(sc3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  # sc3.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  
  sc3.markers.dt <- as.data.table(sc3.markers)
  
  top_most_vars <- get_k_most_distinct_genes_per_cluster(sc3.markers.dt,top_k_pos=10)
  
  fig_prefix <- "07_most_vgenes_heatmap"
  pdf_file <- file.path(wkd,sprintf('%s_%s.pdf',sname,fig_prefix))
  message(pdf_file)
  pdf(pdf_file)
  p <- DoHeatmap(sc3, features = top_most_vars$gene) + NoLegend()
  plot(p)
  dev.off()
  
  tab_prefix <- "08_distinct_expr_genes_per_cluster"
  wb <- createWorkbook()
  sheet_name <- sprintf('%s_top%d_pos_expr_genes',sname,10)
  message(sprintf("sheet_name:%s",sheet_name))
  addWorksheet(wb,sheetName = sheet_name)
  writeData(wb,sheet_name,top_most_vars)
  
  top_most_vars <- get_k_most_distinct_genes_per_cluster(sc3.markers.dt,top_k_pos=1000)
  sheet_name <- sprintf('%s_top%d_pos_expr_genes',sname,1000)
  message(sprintf("sheet_name:%s",sheet_name))
  addWorksheet(wb,sheetName = sheet_name)
  writeData(wb,sheet_name,top_most_vars)
  
  saveWorkbook(wb,file=file.path(wkd,sprintf('%s_%s.xlsx',sname,tab_prefix)),overwrite=T)
  
  # library(devtools)
  # devtools::install_github(repo = "hhoeflin/hdf5r")
  # devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")
  # /media/sammy/projects/apps/share/lib/hdf5-1.10.1/bin
  # pfile <- Convert(from = sc3,
  #                  to = "loom",
  #                  filename = "pbmc_small.loom",
  #                  display.progress = FALSE)
  return(sc3)
}

sc3_normalization <- function(sc_raws,norm_method="default",nfeature_rna=0,mt_pct=0,vst_feature=2000,npcs=50,norm_scale_factor=10000,scale_to_regress=c("percent.mt","S.Score","G2M.Score"),ncpu=8) {
  
  sc_norms <- lapply(names(sc_raws),function(name) {
    
    sc <- sc_raws[[name]]
    if (FALSE) {
      sinfo <- samples[sample == name,]
      sc$sample <- sinfo$sample
      sc$patient <- sinfo$Patient
      sc$timepoint <- sinfo$Timepoint
      sc$type <- sinfo$type
    }
    message(sprintf('computing mitochondria pct [%s] ...',name))
    sc <- PercentageFeatureSet(sc,pattern="^MT-",col.name = "percent.mt")
    
    if (FALSE) {
      fig_tag <- "mitochondria"
      pdf_file <- file.path(wkd,sprintf('%s_%s.pdf',fig_tag,name))
      message(pdf_file)
      pdf(pdf_file)
      hist(sc$percent.mt,xlab = "mitochondria[%]",main = name,breaks = 50)
      dev.off()
    }
    
    message("filtering cells containing high mt ...")
    sc <- subset(sc, cells = which(sc$nFeature_RNA > nfeature_rna & sc$percent.mt < mt_pct))
    
    if (norm_method == "sctf") {
      message(sprintf("norm_method[%s]",norm_method))
      sc <- SCTransform(sc, vars.to.regress = scale_to_regress, verbose = FALSE)
    } else {
      message(sprintf("norm_method[%s]",norm_method))
      message(sprintf('perform a normalization [%s] ...',name))
      sc <- NormalizeData(sc, scale.factor = norm_scale_factor, verbose = FALSE)
    }
    sc <- RunPCA(sc, 
                 npcs = min(npcs,find_max_npcs(sc)),
                 verbose = FALSE)
    
    message(sprintf('find variable features [%s] ...',name))
    if (vst_feature>0) {
      sc <- FindVariableFeatures(sc,
                                 selection.method = "vst",
                                 nfeatures = vst_feature,
                                 verbose = FALSE)
    }
    message(sprintf('Done[%s]',name))
    return(sc)
  }
  # ,mc.cores = ncpu
  ) #,mc.cores = ncpu
  
  names(sc_norms) <- names(sc_raws)
  return(sc_norms)
}


runFIA_with_tryCatch <- function(sc.list,sc.features,D,k_filter=200) {
  sc.anchors <- tryCatch(
    {
      sc.anchors <- FindIntegrationAnchors(object.list = sc.list,
                                           normalization.method = "SCT",
                                           anchor.features = sc.features,
                                           dims = 1:D,
                                           k.filter = k_filter)
    },
    error=function(cond){
      message(sprintf('k_filter[%d]',k_filter))
      return(NA)
    },
    finally={
      message(sprintf('Done[%d]',k_filter))
    }
  )
  return(sc.anchors)
}

integration_analysis <- function(seus,
                                 min.cell=50,
                                 num_features=3000,
                                 D=50,
                                 u.reduction = "rpca",
                                 fc.resolution = 0.8,
                                 scale_to_regress=c("percent.mt","S.Score","G2M.Score"),
                                 debug2=0) {
  
  if (debug2==1) {browser()}
  
  stopifnot(all(scale_to_regress %in% colnames(seus[[1]]@meta.data)))
  
  feat_cnts <- sapply(seus,function(seu){return(dim(seu)[1])})
  mean_feat <- mean(feat_cnts)
  
  D <- 30
  pca_approx <- TRUE
  if (D > mean_feat) {
    D <- round(mean_feat * 0.7)
    pca_approx <- FALSE
  }
  if (num_features > mean_feat) {num_features <- mean_feat}
  
  if (num_features>0) {
    message(sprintf('Find [%d] variable features ...',num_features))
    seus <- lapply(seus,function(seu){
      FindVariableFeatures(seu,selection.method = "vst", nfeatures = num_features)
    })
  }
  
  features <- SelectIntegrationFeatures(object.list = seus)
  seus <- lapply(seus,function(seu) {
    ScaleData(seu, features = features, verbose = FALSE) %>%
      RunPCA(features = features, 
             npcs = min(50,find_max_npcs(.)),
             verbose = FALSE)
  })
  
  message(sprintf('integrate the anchors found to each sample ([%d] anchor features)...',num_features))
  sc.anchors <- FindIntegrationAnchors(object.list = seus, dims = 1:D, reduction = u.reduction)
  rm(seus)
  seui <- IntegrateData(anchorset = sc.anchors, dims = 1:D)
  rm(sc.anchors)
  # -----------------
  DefaultAssay(seui) <- "integrated"
  
  message('cell counts per sample')
  summary(seui@active.ident)
  
  message('data matrix dimension')
  dim(seui@assays$integrated@data)
  
  # Run the standard workflow for visualization and clustering
  seui <- regress_out_scaledata(seui,
                                scale_to_regress=scale_to_regress)
  
  D = min(D,find_max_npcs(seui))
  
  seui <- RunPCA(seui, npcs = D, verbose = TRUE, approx=pca_approx)
  
  # UMAP and Clustering
  seui <- RunUMAP(seui, dims = 1:D)
  
  seui <- FindNeighbors(seui, dims = 1:D)
  seui <- FindClusters(seui, resolution = fc.resolution)
  
  return(seui)
}

rm_assay <- function(seus,assay2default="RNA",assay2del="SCT") {
  
  message(sprintf("removing assay [%s] ...",assay2del))
  snames<-names(seus)
  seus <- lapply(seus, function(seu){
    avail_assays <- names(seu@assays)
    stopifnot(assay2default %in% avail_assays)
    DefaultAssay(seu) <- assay2default
    if (assay2del %in% avail_assays) {
      seu[[assay2del]] <- NULL
    }
    return(seu)
  })
  names(seus) <- snames
  return(seus)
}

merge_seurats <- function(seus,add_sample_idx=TRUE,enable_diet=1,debug2=0) {
  if (debug2==1){browser()}
  S <- length(seus)
  stopifnot(S>1)
  
  message("appending an index to distinguish cbc ...")
  snames<-names(seus)
  if (add_sample_idx) {
    for (s in 1:S) {
      seus[[s]] <- RenameCells(seus[[s]],add.cell.id = snames[[s]])
      
      if (enable_diet==1) {
        seus[[s]] = DietSeurat(seus[[s]])
      }
    }
  }
  names(seus) <- snames
  
  message(sprintf("merging [%d] seurat objects ...",S))
  seui <- seus[[1]]
  for (s in 2:S) {
    message(s)
    seui <- merge(seui,seus[[s]])
  }
  seui
}

#####################

aligned_scts_by_harmony <- function(iseu,
                                    D=30,
                                    nfeatures=3000,
                                    # min.cells=50,
                                    fc.resolution=0.8,
                                    harmony_to_regress=c("orig.ident"),
                                    uvar_genes = NA,
                                    excl_pat="^IG[HJKL]",
                                    do_cluster=1,
                                    harmony=1,
                                    do_tsne=0,
                                    debug2=0) {
  if (debug2==1){browser()}
  
  cmeta = colnames(iseu[[1]]@meta.data)
   
  sapply(iseu,function(seu){
    'SCT' %in% names(seu@assays)
  }) %>% stopifnot()
  
  if (harmony==1) {
    stopifnot(all(harmony_to_regress %in% cmeta))
  }
  
  # Find most variable features across samples to integrate
  integ_features <- SelectIntegrationFeatures(object.list = iseu,
                                              nfeatures = nfeatures) 
  
  if (!invalid(uvar_genes)){
    uvar_genes = uvar_genes[uvar_genes %in% rownames(iseu)]
    integ_features = unique(c(integ_features,uvar_genes))
  }
  
  # Merge normalized samples
  iseu <- merge(x = iseu[[1]],
                y = iseu[2:length(iseu)],
                merge.data = TRUE)
  
  DefaultAssay(iseu) <- "SCT"
  
  # Manually set variable features of merged Seurat object
  VariableFeatures(iseu) <- integ_features
  
  if (do_cluster==0 & harmony==1){
    do_cluster=1
  }
  
  D = min(D,find_max_npcs(iseu))
  
  if (do_cluster==1) {
    message("--RunPCA ...")
    # Calculate PCs using manually set variable features
    iseu <- RunPCA(iseu,
                   assay = 'SCT',
                   npcs = D,
                   features = VariableFeatures(object = iseu))
  }
  
  red_method="pca"
  if (harmony==1) {
    message("running harmony ...")
    iseu <- RunHarmony(iseu, assay="SCT", group.by.vars = harmony_to_regress)
    red_method="harmony"
  }
  
  if (do_cluster==1) {
    iseu <- perform_seu_clustering(iseu,D,fc.resolution,red_method=red_method,do_tsne = do_tsne)
  }
  
  message("Done[harmony].")
  if (debug2==1){browser()}
  return(iseu)
}

#####################
#checked 05/14/2021
aligned_by_harmony <- function(seus,
                               D=30,
                               min.cells=50,
                               fc.resolution=0.8,
                               scale_to_regress=c("percent.mt","S.Score","G2M.Score"), #"nCount_RNA"
                               harmony_to_regress=c("orig.ident"),
                               uvar_genes = NA,
                               excl_pat="^IG[HJKL]",
                               do_cluster=1,
                               harmony=1,
                               do_tsne=0,
                               debug2=0) {
  if (debug2==1){browser()}
  
  cmeta = colnames(seus[[1]]@meta.data)
  if (harmony==1) {
    stopifnot(all(c(scale_to_regress,harmony_to_regress) %in% cmeta))
  } else {
    stopifnot(all(scale_to_regress %in% cmeta))
  }
  
  cols_oi = c(scale_to_regress,harmony_to_regress)
  cols_oi = cols_oi[!is.na(cols_oi) & cols_oi!=""]
  
  stopifnot(all(cols_oi %in% cmeta))
  
  message("deleting SCT assay ...")
  seui <- merge_seus_norm_fvar(seus,uvar_genes,excl_pat)
  message("--ScaleData ...")
  seui <- regress_out_scaledata(seui,
                                scale_to_regress=scale_to_regress)
  
  seui = FindVariableFeatures(seui)
  
  if (do_cluster==0 & harmony==1){
    do_cluster=1
  }
  D = min(D,find_max_npcs(seui))
  if (do_cluster==1) {
    message("--RunPCA ...")
    seui <- RunPCA(seui,
                   npcs = find_max_npcs(seui,D=D),
                   features = VariableFeatures(object = seui))
  }
  
  red_method="pca"
  if (harmony==1) {
    message("running harmony ...")
    seui <- RunHarmony(seui, group.by.vars = harmony_to_regress)
    red_method="harmony"
  }
  
  if (do_cluster==1) {
    seui <- perform_seu_clustering(seui,D,fc.resolution,red_method=red_method)
  }
  
  message("Done[harmony].")
  if (debug2==1){browser()}
  return(seui)
}

perform_seu_clustering <- function(seui,D=30,fc.resolution=0.8,do_tsne=0,red_method="pca") {
  
  # Dimensional reduction and plotting
  message(sprintf("--RunUMAP [%s] ...",red_method))
  seui <- RunUMAP(seui, dims = 1:D, reduction = red_method)
  
  if (do_tsne==1) {
    message("--RunTSNE ...")
    seui <- RunTSNE(seui, dims = 1:D, reduction = red_method)
  }
  
  message("--FindNeighbors ...")
  seui <- FindNeighbors(seui, reduction = red_method, dims = 1:D)
  
  message("--FindClusters ...")
  seui <- FindClusters(seui, resolution = fc.resolution)
  seui
}

cellcycle_genes_avail <- function(seu,debug2=0) {
  if (debug2==1){browser()}
  feats.avail <- rownames(GetAssayData(seu))
  ccg <- cc.genes.updated.2019
  ccg$s.genes <- ccg$s.genes[ccg$s.genes %in% feats.avail]
  ccg$g2m.genes <- ccg$g2m.genes[ccg$g2m.genes %in% feats.avail]
  return(ccg)
}

cellcycle_scoring <- function(seu,debug2=0) {
  seu <- tryCatch(
    {
      assays = names(seu@assays)
      if (debug2==1){browser()}
      
      default.assay.bkp = DefaultAssay(seu)
      if ('RNA' %in% assays) {
        temp.assay="RNA"
      } else if ('Spatial' %in% assays) {
        temp.assay="Spatial"
      } else {
        stop('RNA or Spatial assay is required!')
      }
      
      message(sprintf("set assay to [%s] temporaly",temp.assay))
      DefaultAssay(seu)=temp.assay
      
      #check if the assay was normalized prevly
      normalize_done = 1
      if (nrow(GetAssayData(seu,slot="data"))==0) {
        normalize_done = 0
        seu=NormalizeData(seu)
      } else if (identical(GetAssayData(seu,slot="counts")@x,
                           GetAssayData(seu,slot="data")@x)) {
        normalize_done = 0
        seu=NormalizeData(seu)
      }
      
      ccg <- cellcycle_genes_avail(seu,debug2=debug2)
      if (debug2==1){browser()}
      feats.avail <- rownames(GetAssayData(seu))
      
      if (any(ccg$s.genes %in% feats.avail) & any(ccg$g2m.genes %in% feats.avail)) {
        message("assigning cell cycle scores ...")
        if (debug2==1){browser()}
        seu <- CellCycleScoring(seu,
                                s.features = ccg$s.genes[ccg$s.genes %in% feats.avail],
                                g2m.features = ccg$g2m.gene[ccg$g2m.gene %in% feats.avail],
                                set.ident = TRUE)
        
        if (normalize_done==0){
          seu[[temp.assay]]@data@x = rlang::duplicate(seu[[temp.assay]]@counts@x)
        }
        DefaultAssay(seu) = default.assay.bkp
        return(seu)
      }
    },
    error=function(cond){
      message(sprintf("cellcycle_scoring failed![%s]",cond))
      return(NA)
    },
    warning=function(cond){
      message(sprintf("cellcycle_scoring caused warning![%s]",cond))
      return(seu)
    },
    finally={
      message("Done[cellcycle_score].")
    }
  )
  return(seu)
}

regress_out_scaledata <- function(seu,
                                  scale_to_regress=c("percent.mt","nCount_RNA"),
                                  do_center=TRUE,
                                  split_by=NA,
                                  debug2=0) {
  if (debug2==1){browser()}
  
  stopifnot(all(scale_to_regress %in% colnames(seu@meta.data)))
  
  message(sprintf("scaling/regressing out [%s]...",paste0(scale_to_regress,collapse = ",")))
  
  if (length(scale_to_regress)>0) {
    
    if (invalid(split_by)) {
      seu <- ScaleData(seu,
                       do.center=do_center,
                       vars.to.regress = scale_to_regress)
    } else {
      seu <- ScaleData(seu,
                       do.center=do_center,
                       split.by=split_by,
                       vars.to.regress = scale_to_regress)
    }
  } else {
    if (invalid(split_by)) {
      seu <- ScaleData(seu,do.center=do_center)
    } else {
      seu <- ScaleData(seu,do.center=do_center,split.by=split_by)
    }
  }
  return(seu)
}

find_max_npcs <- function(seu,D=50) {
  npca = min(D,round(dim(seu)[2] * 0.65))
  npca
}

find_optimal_npca <- function(seuj) {
  # Find significant PCs
  stopifnot('pca' %in% names(seuj))
  
  stdv <- seuj[["pca"]]@stdev
  sum.stdv <- sum(seuj[["pca"]]@stdev)
  percent.stdv <- (stdv / sum.stdv) * 100
  cumulative <- cumsum(percent.stdv)
  co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
  co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] -
                       percent.stdv[2:length(percent.stdv)]) > 0.1),
              decreasing = T)[1] + 1
  min.pc <- min(co1, co2)
  message(sprintf("found optimal # of PC [%d] ...",min.pc))
  min.pc
}

IntegrateData_with_tryCatch <- function(iseu,k_weight=100,debug2=0) {
  if (debug2==1){browser()}
  iseu <- tryCatch(
    {
      message(sprintf('running integrateData(k.weight=%d)',k_weight))
      iseu <- IntegrateData(anchorset = iseu, 
                            k.weight = k_weight,
                            normalization.method = 'SCT')
    },
    error=function(cond){
      message(sprintf('integrateData failed w/ k.weight[%d]',k_weight))
      return(NULL)
    },
    finally={
      message('<integratedata_with_trycatch>')
    }
  )
  return(iseu)
}

perform_cca <- function(iseu, si_features=3000, k_weight=100,k_filter=200,feats_to_incl=NA,red_method='cca',npca=50,max_iter=10,debug2=0) {
  if (debug2==1){browser()}
  features <- SelectIntegrationFeatures(object.list = iseu,
                                        nfeatures = si_features)
  
  if (!invalid(feats_to_incl)) {
    features = c(features,feats_to_incl) %>% unique()
  }
  
  total_n_samples = length(iseu)
  features = imap(iseu,function(seu,sname) {
    feats = GetAssayData(seu,assay = 'SCT',slot = 'scale.data') %>% rownames()
    data.table(feat=features[features %in% feats],
               sample=sname)
  }) %>% rbindlist() %>%
    dplyr::group_by(feat) %>%
    dplyr::summarise(n_cnt = n()) %>%
    dplyr::filter(n_cnt==total_n_samples) %>%
    dplyr::pull(feat)
  
  iseu <- PrepSCTIntegration(object.list = iseu,
                             anchor.features = features)
  
  iseu <- FindIntegrationAnchors(object.list = iseu,
                                 dims = 1:npca,
                                 reduction = red_method,
                                 k.filter=k_filter,
                                 normalization.method = "SCT",
                                 anchor.features = features)
  
  iter=0
  was_successful=0
  while(iter<max_iter & was_successful==0) {
    ret <- IntegrateData_with_tryCatch(iseu, k_weight = k_weight)
    if (!invalid(ret)) {
      iseu = ret
      was_successful = 1
    }
    iter = iter + 1
    k_weight = k_weight - 1
  }
  
  iseu <- RunPCA(iseu, 
                 verbose = FALSE)
  
  npca.nb = find_optimal_npca(iseu)
  
  message(sprintf("--FindNeighbors ..."))
  iseu <- FindNeighbors(iseu, 
                        dims = 1:npca.nb, 
                        verbose = TRUE)
  
  message(sprintf("--FindClusters ..."))
  iseu <- FindClusters(iseu, 
                       verbose = TRUE)
  
  RunUMAP(iseu, reduction = "pca", dims = 1:npca.nb)
}

# iseu : a list of seurats objects
integration_analysis_sctf <- function(iseu,
                                      si_features=3000,
                                      feats_to_incl=NA,
                                      k_weight=100,
                                      k_filter=200,
                                      max_iter=50,
                                      vtrs=c("nCount_RNA","percent.mt"),
                                      npca=0,
                                      reuse_sct=1, 
                                      red_method='cca',
                                      input_assay="RNA",
                                      ncpu=4,
                                      debug2=0) {
  if (debug2==1){browser()}
  message(sprintf("performing data integration (SCT,%s) ...",red_method))
  # ---------
  if (reuse_sct==0 | !('SCT' %in% names(iseu[[1]]@assays))) {
    # plan2(ncpu=ncpu)
    iseu = imap(iseu,function(seu,sname) {
      stopifnot((input_assay %in% names(seu@assays)))
      stopifnot(all(vtrs %in% colnames(seu@meta.data)))
      
      message(sname)
      seu = DietSeurat(seu) %>%
        SCTransform(object = ., 
                    vars.to.regress = vtrs, 
                    verbose = FALSE, 
                    method="glmGamPoi", 
                    assay=input_assay,
                    return.only.var.genes = FALSE)
    })
    # plan2()
  }
  if (debug2==1){browser()}
  
  if (npca==0) {
    npca = sapply(iseu, function(iseu) {
      round(ncol(iseu)*0.65)
    }) %>%
      min(c(.,50))
    message('max npcs = [%d]',npca)
  }
  
  if (red_method=="rpca") {
    iseu = imap(iseu,function(seu,sname) {
      RunPCA(seu, 
             npcs = npca)
    })
  }
  
  iseu = perform_cca(iseu=iseu,
                     si_features=si_features,
                     feats_to_incl=feats_to_incl,
                     k_weight=k_weight,
                     k_filter=k_filter,
                     max_iter = max_iter,
                     red_method=red_method,
                     npca=npca,
                     debug2 = debug2)
  
  if (debug2==1){browser()}
  
  iseu
}

print_integ_result <- function(sci,
                               pdf_prefix,
                               split_by,
                               group_by="seurat_clusters",
                               ncol = 3,
                               pdfwinch=0,
                               pdfhinch=0,
                               debug=0) {
  
  pdf_file <- sprintf("%s_dimPlot.pdf",pdf_prefix)
  orig_sample_cnt <- length(unique(sci@meta.data$orig.ident))
  
  if (debug==1) {browser()}
  # else {pdf(pdf_file,height=5.5,width=11)}
  else {
    if (pdfwinch>0){
      pdf(pdf_file,width = pdfwinch,height = pdfhinch)
    } else {
      pdf(pdf_file)
    }
  }
  
  p1 <- DimPlot(sci,
                reduction = "umap",
                group.by = split_by,
                label = FALSE,
                combine=TRUE)
  
  plot(p1,asp=1.)
  
  split_by_ids <- unique(sci[[split_by]])
  
  N <- dim(split_by_ids)[1]
  M <- ceiling(N/4)
  for (m in 1:M) {
    st1 <- (m-1)*4 + 1
    ed2 <- m*4
    if (ed2>N) {ed2 <- N}
    
    seu <- subset(sci,cells = Cells(sci)[sci@meta.data[,split_by] %in% split_by_ids[st1:ed2,]])
    
    p2 <- DimPlot(seu,
                  reduction = "umap",
                  split.by = split_by,
                  group.by = "seurat_clusters",
                  ncol=2,
                  label = TRUE,
                  combine=TRUE)
    plot(p2)
  }
  
  if (group_by != "seurat_clusters") {
    if (group_by %in% colnames(sci@meta.data)) {
      for (m in 1:M) {
        st1 <- (m-1)*4 + 1
        ed2 <- m*4
        if (ed2>N) {ed2 <- N}
        
        seu <- subset(sci,cells = Cells(sci)[sci@meta.data[,split_by] %in% split_by_ids[st1:ed2,]])
        
        p3 <- DimPlot(seu,
                      reduction = "umap",
                      split.by = split_by,
                      group.by = group_by,
                      ncol=2,
                      label = FALSE,
                      combine=TRUE)
        plot(p3)
      }
      
    } else {
      message(sprintf("group_by[%s] does not exist in sci@meta.data",group_by))
    }
  }
  
  if (debug==0) {
    dev.off()
  }
}

print_integ_featurePlot <- function(sc.integrated,pdf_prefix,group_by,goi,debug2=0) {
  
  library(randomcoloR)
  if (debug2==1){browser()}
  DefaultAssay(sc.integrated) <- "SCT"
  
  features <- rownames(sc.integrated@assays$SCT@data)
  goi_eff <- sort(goi[goi %in% features])
  
  L <- length(goi_eff)
  if (L == 0) {
    return(FALSE)
  }
  
  S <- dim(unique(sc.integrated[[group_by]]))[1]
  
  if (L > 6) {
    goi_effs <- split(goi_eff, ceiling(seq_along(goi_eff)/(L/ceiling(L/6))))
  } else {
    goi_effs <- list(goi_eff)
  }
  
  for (i in 1:length(goi_effs)) {
    if (debug2==0){
      pdf(sprintf("%s_featurePlot_%d.pdf",pdf_prefix,i),width=5.5,height=11)
    }
    
    if (debug2==1){browser()}
    my_pal <- distinctColorPalette(S)
    
    p <- FeaturePlot(sc.integrated,
                     features = goi_effs[[i]],
                     split.by = group_by)
    plot(p)
    
    plots <- VlnPlot(sc.integrated,
                     features = goi_effs[[i]],
                     split.by = group_by,
                     pt.size = 0,
                     combine = FALSE)
    
    p <- CombinePlots(plots = plots, ncol = 1)
    plot(p)
    
    
    p <- RidgePlot(sc.integrated,
                   features = goi_effs[[i]],
                   ncol = 1,
                   group.by = group_by)
    plot(p)
    
    p <- DoHeatmap(sc.integrated,
                   features = goi_effs[[i]],
                   size = 3,
                   group.by = group_by,
                   assay ="SCT",
                   slot = "scale.data")
    plot(p)
    if (debug2==1){browser()}
    if (debug2==0){dev.off()}
  }
  return(TRUE)
}

#' https://divingintogeneticsandgenomics.rbind.io/post/customize-featureplot-in-seurat-for-multi-condition-comparisons-using-patchwork/
#' p_list<- FeaturePlotSingle(pbmc, feature= "MS4A1", metadata_column = "samples", pt.size = 0.05, order =TRUE)
#'
FeaturePlotSingle<- function(obj, feature, metadata_column, assay2="RNA", red2="umap", ...) {
  
  all_cells<- colnames(obj)
  # groups<- levels(obj@meta.data[, metadata_column])
  groups<- unique(obj@meta.data[[metadata_column]])
  
  # the minimal and maximal of the value to make the legend scale the same.
  minimal<- min(obj[[assay2]]@data[feature, ])
  maximal<- max(obj[[assay2]]@data[feature, ])
  ps<- list()
  for (group in groups) {
    subset_indx<- obj@meta.data[, metadata_column] == group
    subset_cells<- all_cells[subset_indx]
    p<- FeaturePlot(obj, features = feature, cells= subset_cells, reduction=red2, raster = FALSE, ...) +
      scale_color_viridis_c(limits=c(minimal, maximal), direction = 1) +
      ggtitle(group) +
      theme(plot.title = element_text(size = 10, face = "bold"))
    ps[[group]]<- p
  }
  return(ps)
}


markedFeaturePlot_solo <- function(goi,sc3,pdf_prefix,plotName) {
  
  DefaultAssay(sc3) <- "RNA"
  features <- rownames(sc3@assays$RNA@data)
  goi_eff <- goi[goi %in% features]
  
  pdf_file <- sprintf("%s.pdf",pdf_prefix)
  if (plotName == "FeaturePlot") {
    plot_list <- FeaturePlot(sc3,features=goi_eff,pt.size=0.2,coord.fixed=TRUE,combine = FALSE)
  } else if (plotName == "DoHeatmap") {
    plot_list <- DoHeatmap(sc3, features=goi_eff,combine = FALSE)
  }
  
  multi.page <- ggarrange(plotlist=plot_list,
                          ncol = 2,
                          nrow = 2)
  message(pdf_file)
  ggexport(multi.page, filename = pdf_file)
  message("Done.")
}

depreciated_markedFeaturePlot_solo <- function(goi,sc3,pdf_prefix,plotName) {
  source(file.path(Sys.getenv('R_UTIL'),'lib_util.R'))
  
  # browser()
  DefaultAssay(sc3) <- "RNA"
  features <- rownames(sc3@assays$RNA@data)
  goi_eff <- goi[goi %in% features]
  
  L <- length(goi_eff)
  if (L == 0) {
    return(FALSE)
  }
  goi_eff <- sort(goi_eff)
  
  if (plotName=="FeaturePlot") {
    M <- 4
  } else if (plotName == "DoHeatmap") {
    M <- 100
  }
  
  if (L > M) {
    goi_effs <- chunk2(goi_eff,M)
  } else {
    goi_effs <- list(goi_eff)
  }
  
  pdf_file <- sprintf("%s.pdf",pdf_prefix)
  message(pdf_file)
  pdf(pdf_file)
  
  for (i in 1:length(goi_effs)) {
    if (length(goi_effs[[i]]) > 0) {
      if (plotName == "FeaturePlot") {
        p <- FeaturePlot(sc3,features=goi_effs[[i]],ncol=2,coord.fixed=TRUE)
      } else if (plotName == "DoHeatmap") {
        p <- DoHeatmap(sc3, features = goi_effs[[i]]) + NoLegend()
      }
      plot(p)
    }
  }
  dev.off()
}

diff_analy_prep <- function(sci) {
  sci$cluster_sgroup <- paste0(Idents(sci),'_',sci$sgroup)
  sci$cluster <- Idents(sci)
  Idents(sci) <- "cluster_sgroup"
  return(sci)
}


diff_analy_anchors <- function(sci,comps,gsample='group_sample',min_cell_cnts=3) {
  
  #comps: a pair of sgroups to compare with
  #comps <- data.table(ctrl=c('pre-infusion','day14','pre-infusion'),
  #	expr=c('day14','day30','day30'))
  # browser()
  DefaultAssay(sci) <- "RNA"
  
  cluster_labels <- sort(unique(sci@meta.data$seurat_clusters))
  
  sci2 <- diff_analy_prep(sci)
  diffMarkers <- list()
  for (clab in cluster_labels) {
    
    message(sprintf("cluster_label:%s",clab))
    for (r in 1:dim(comps)[1]) {
      message(sprintf("sgroup:%s vs. %s",comps$ctrl[r],comps$expr[r]))
      ctrl_lab <- sprintf("%s_%s",clab,comps$ctrl[r])
      expr_lab <- sprintf("%s_%s",clab,comps$expr[r])
      
      if ((length(which(sci2$cluster_sgroup==ctrl_lab)) >= min_cell_cnts) & (length(which(sci2$cluster_sgroup==expr_lab)) >= min_cell_cnts)) {
        message('enough samples for diff analysis.')
        comp_label <- sprintf("%s_%s_%d",gsample,clab,r)
        fm_result <- FindMarkers(sci2,
                                 ident.1=expr_lab,
                                 ident.2=ctrl_lab,
                                 verbose=TRUE)
        diffMarkers[[comp_label]] <- as.data.table(fm_result)
        diffMarkers[[comp_label]]$gsample <- gsample
        diffMarkers[[comp_label]]$cluster <- clab
        diffMarkers[[comp_label]]$gene <- rownames(fm_result)
        diffMarkers[[comp_label]]$sgroups <- sprintf("%s_%s",comps$ctrl[r],comps$expr[r])
      }
    }
  }
  
  return(rbindlist(diffMarkers))
  
}


tmp_split_multilanes_crcount <- function(gene_names,seurat3_raws,total_sidx) {
  
  rcnts <- lapply(names(seurat3_raws),function(sname) {
    
    message(sname)
    srd1 <- seurat3_raws[[sname]]
    features <- rownames(srd1@assays$RNA@counts)
    cell_bc <- colnames(srd1@assays$RNA@counts)
    
    bc_start_idx <- sapply(1:total_sidx,function(i){which(str_detect(cell_bc,as.character(i)))[1]})
    bc_start_idx <- c(bc_start_idx,length(cell_bc)+1)
    rcnts2 <- list()
    for (i in 1:(length(bc_start_idx)-1)) {
      message(i)
      from_bc_idx <- bc_start_idx[i]
      to_bc_idx <- bc_start_idx[i+1]-1
      
      rcnts2[[i]] <- as.data.table(Matrix::rowSums(srd1@assays$RNA@counts[,(from_bc_idx:to_bc_idx)]))
    }
    
    # browser()
    rcnt <- do.call(cbind,rcnts2)
    colnames(rcnt) <- sprintf("%d",1:total_sidx)
    rcnt2 <- rcnt[match(gene_names,features),]
    rcnt2$gene <- gene_names
    return(rcnt2)
  })
  
  names(rcnts) <- names(seurat3_raws)
  
  return(rcnts)
}


report_sc3_dim <- function(seus,tsv_fpath=NA,qc_comment="N",append2=FALSE,debug2=0) {
  if (debug2==1){browser()}
  message(tsv_fpath)
  rc_dims <- sapply(names(seus),function(sname) {
    return(dim(seus[[sname]]@assays$RNA@counts))
  }, simplify = TRUE)
  
  rc_dims <- as.data.table(t(rc_dims),keep.rownames=TRUE)
  colnames(rc_dims) <- c("sample","features","cells")
  
  if (!invalid(tsv_fpath)) {
    rc_dims$qc_comment <- qc_comment
    fwrite(rc_dims,file=tsv_fpath,append=append2,sep="\t")
  }
  if (debug2==1){browser()}
  return(rc_dims)
}

cell_typing <- function(cell_marker,dge_per_cluster) {
  
  marker_genes <- rownames(cell_marker)
  ctypes <- colnames(cell_marker)
  
  dge_ctypes <- list()
  
  for (j in 1:dim(cell_marker)[2]) {
    # j <- 2
    ctype <- ctypes[j]
    mg_indicator <- cell_marker[,j]>0
    
    mprof <- get_mprof(marker_genes[mg_indicator])
    
    sample_tags <- names(dge_per_cluster)
    dgeis <- list()
    for (sample_tag in sample_tags) {
      message(sprintf("j:%d/sample_tag:%s",j,sample_tag))
      
      dgei <- as.data.table(dge_per_cluster[[sample_tag]])
      dgei <- dgei[order(avg_logFC,-p_val)]
      dgei$sample <- sample_tag
      dgei$ctype <- ctype
      
      dgei$rscore <- 1:dim(dgei)[1]/dim(dgei)[1]
      dgei_up <- dgei[gene %in% mprof$up & avg_logFC>=0.,]
      
      if (length(mprof$down)>0) {
        dgei$rscore <- (1-1:dim(dgei)[1]/dim(dgei)[1])
        dgei_down <- dgei[gene %in% mprof$down & avg_logFC<=0.,]
        dgei <- rbind(dgei_up,dgei_down)
      } else {
        dgei <- dgei_up
      }
      dgeis[[sample_tag]] <- dgei
    }
    dge_ctypes[[j]] <- rbindlist(dgeis)
  }
  dge_dt <- rbindlist(dge_ctypes)
  
  dge_by_samples <- split(dge_dt,by=c("sample"))
  
  sample_cluster_ctype_profs <- list()
  
  for (sname in names(dge_by_samples)) {
    dge_by_clusters <- split(dge_by_samples[[sname]],by='cluster')
    dge_by_clusters <- dge_by_clusters[!isEmpty(dge_by_clusters)]
    
    for (cid in names(dge_by_clusters)) {
      
      name2 <- sprintf("%s_%s",sname,cid)
      message(name2)
      
      dge_cluster <- dge_by_clusters[[cid]]
      
      dge_df <- dcast(dge_cluster,ctype ~ gene, value.var = 'rscore')
      
      ctypes <- dge_df$ctype
      dge_df <- dge_df[,2:dim(dge_df)[2]]
      
      dge_df[is.na(dge_df)] <- 0.
      
      topK <- rowSums(as.data.table(dge_df) > 0)
      topK[topK>1] <- 2
      topK[topK==1] <- 1
      
      dge_df <- cbind(dge_df,topK)
      
      rownames(dge_df) <- ctypes
      
      D <- dim(dge_df)[2]
      
      topKscore <- apply(dge_df,1,function(cvals) {
        rscores <- cvals[1:(D-1)]
        topK <- cvals[D]
        meanScore <- mean(sort(rscores,decreasing = T)[1:topK])
        return(meanScore)
      })
      
      ctype2 <- names(topKscore)
      cscore <- as.data.table(topKscore)
      cscore$ctype <- ctype2
      cscore$cluster <- cid
      cscore$sample <- sname
      
      cscore <- cscore[order(-topKscore)]
      sample_cluster_ctype_profs[[name2]] <- cscore
    }
  }
  
  sample_cluster_ctype_profs <- rbindlist(sample_cluster_ctype_profs)
  
  #dge_by_sample_cids <- dge_by_sample_cids[!isEmpty(dge_by_sample_cids)]
  screp <- as.data.table(dcast(sample_cluster_ctype_profs,sample + cluster ~ ctype, fun.aggreate = max, value.var = "topKscore"))
  
  screps <- split(screp,by="sample")
  
  for (sname in names(screps)) {
    screpj <- screps[[sname]]
    M <- dim(screpj)[2]
    cluster_id <- screpj$cluster
    screpj <- screpj[,3:M]
    ctypes <- colnames(screpj)
    
    rownames(screpj) <- paste0('cluster_',cluster_id)
    colnames(screpj) <- tcell_to_fullname(ctypes)
    pdf_file <- file.path(wkd,sprintf('%s_tcellSubtype_x_cluster.pdf',sname))
    message(pdf_file)
    screpj[is.na(screpj)] <- 0.
    pheatmap(screpj,main = sprintf("sctf_%s[T-cell subtypes vs. cluster]",sname),file=pdf_file)
    
  }
  
  return(screps)
}


integs_to_dgeTable <- function(sc.integs,comps,wkd,min_cell_cnts=3) {
  
  # comps <- data.table(ctrl=c('pre-infusion','day14','pre-infusion'),
  # 										expr=c('day14','day30','day30'))
  
  message("identify differentially expressed genes across conditions ...")
  
  group_by <- "sgroup"
  
  group_samples <- names(sc.integs)
  
  diffMarkers <- list()
  
  for (gsample in group_samples) {
    
    message(sprintf("gsample:%s",gsample))
    
    sci <- sc.integs[[gsample]]
    if (!invalid(sci)) {
      DefaultAssay(sci) <- "RNA"
      sci <- NormalizeData(sci, verbose = FALSE)
      
      cluster_labels <- sort(unique(sci@meta.data$seurat_clusters))
      
      sci2 <- seurat3_diff_analy_prep(sci)
      
      for (clab in cluster_labels) {
        message(sprintf("cluster_label:%s",clab))
        for (r in 1:dim(comps)[1]) {
          message(sprintf("sgroup:%s vs. %s",comps$ctrl[r],comps$expr[r]))
          ctrl_lab <- sprintf("%s_%s",clab,comps$ctrl[r])
          expr_lab <- sprintf("%s_%s",clab,comps$expr[r])
          
          if ((length(which(sci2$cluster_sgroup==ctrl_lab)) >= min_cell_cnts) & (length(which(sci2$cluster_sgroup==expr_lab)) >= min_cell_cnts)) {
            message('enough samples for diff analysis.')
            comp_label <- sprintf("%s_%s_%d",gsample,clab,r)
            fm_result <- FindMarkers(sci2,
                                     ident.1=expr_lab,
                                     ident.2=ctrl_lab,
                                     verbose=TRUE)
            diffMarkers[[comp_label]] <- as.data.table(fm_result)
            diffMarkers[[comp_label]]$gsample <- gsample
            diffMarkers[[comp_label]]$cluster <- clab
            diffMarkers[[comp_label]]$gene <- rownames(fm_result)
            diffMarkers[[comp_label]]$sgroups <- sprintf("%s_%s",comps$ctrl[r],comps$expr[r])
          }
        }
      }
    }
  }
  
  diffMarker_dt <- rbindlist(diffMarkers)
  warnings()
  save(diffMarker_dt,file=file.path(wkd,'seurat3_diff.rd'),compress=TRUE)
  
  wb <- createWorkbook("seurat3_diff_analysis")
  dts <- split(diffMarker_dt,by=c("gsample","sgroups"))
  
  for (i in 1:length(dts)) {
    sheet_name <- names(dts)[i]
    addWorksheet(wb,sheetName = sheet_name)
    writeData(wb,sheet_name,dts[[i]])
  }
  saveWorkbook(wb,file=file.path(wkd,"seurat3_diff.xlsx"),overwrite=T)
}


export_sc3_to_text <- function(seu,wkd,exp_tag="sample",slot="norm",debug=0) {
  
  if (debug==1){browser()}
  
  if (!file.exists(wkd)) {dir.create(wkd,showWarnings = FALSE, recursive=TRUE)}
  
  if (slot == "norm") {
    sc3_dt <- as.data.table(seu@assays$SCT@data)
    genes <- rownames(seu@assays$SCT@data)
  } else {
    sc3_dt <- as.data.table(seu@assays$RNA@counts)
    genes <- rownames(seu@assays$RNA@counts)
  }
  
  cluster_avail <- 0
  if ("seurat_clusters" %in% colnames(seu@meta.data)) {
    cluster_avail <- 1
  }
  
  cb_mat <- colnames(sc3_dt)
  
  meta_info <- data.table(cb=rownames(seu@meta.data),
                          orig.ident=seu@meta.data$orig.ident)
  
  if (cluster_avail==1){
    meta_info$cluster_id <- seu@meta.data$seurat_clusters
  }
  sample_cb <- paste0(meta_info[match(cb_mat,meta_info$cb),orig.ident],'.',meta_info$cb)
  
  colnames(sc3_dt) <- sample_cb
  
  sc3_dt$gene <- genes
  
  N <- dim(sc3_dt)[2]
  sc3_dt <- cbind(sc3_dt[,N,with=F],sc3_dt[,1:(N-1),with=F])
  
  tsv_fpath <- file.path(wkd,sprintf("%s_%s_data_matrix.tsv",exp_tag,slot))
  fwrite(sc3_dt,file=tsv_fpath,sep="\t",row.names=FALSE,col.names = TRUE)
  tsv_gz <- file.path(wkd,sprintf("%s_%s_matrix.tsv.gz",exp_tag,slot))
  if (file.exists(tsv_gz)) {unlink(tsv_gz)}
  gzip(tsv_fpath,destname=tsv_gz)
  
  if (cluster_avail==1) {
    tsv_fpath <- file.path(wkd,sprintf("%s_%s_cluster_info.tsv",exp_tag,slot))
    meta_info$sample_cb <- paste0(meta_info$orig.ident,'.',meta_info$cb)
    fwrite(meta_info,file=tsv_fpath,sep="\t",row.names=FALSE,col.names = TRUE)
    tsv_gz <- file.path(wkd,sprintf("%s_%s_cluster_info.tsv.gz",exp_tag,slot))
    if (file.exists(tsv_gz)) {unlink(tsv_gz)}
    gzip(tsv_fpath,destname=tsv_gz)
  }
}

export_sc3_to_loom <- function(seu,wkd,exp_tag,slot="raw") {
  
  # library(devtools)
  # devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")
  # library(hdf5r)
  # library(loomR)
  
  if (slot=="norm") {
    seu_mtrx <- seu@assays$RNA@data
  } else {
    seu_mtrx <- seu@assays$RNA@counts
  }
  
  loom_fpath <- file.path(wkd,sprintf("%s_data_matrix.loom",exp_tag))
  
  genes<- list(rownames(seu_mtrx))
  names(genes) <- genes
  
  cids<- list(colnames(seu_mtrx))
  names(cids) <- cids
  
  loomR::create(filename=loom_fpath,
                data=seu_mtrx,
                gene.attrs=genes,
                cell.attrs=cids,
                do.transpose=TRUE,
                overwrite=TRUE,
                verbose=FALSE)
  
  return(loom_fpath)
}

export_to_tsv <- function(seu,out_pref,assay2="RNA",slot2="data") {
  
  mtrx <- as.matrix(GetAssayData(object = seu,
                                 assay=assay2,
                                 slot = slot2))
  
  mtrx <- as.data.table(mtrx,keep.rownames=TRUE)
  setnames(mtrx,'rn','gene')
  fwrite(mtrx,
         file=sprintf("%s_mtrx.tsv.gz",out_pref),
         sep="\t",
         col.names=TRUE,
         quote=FALSE)
  
  smeta <- as.data.table(seu@meta.data,keep.rownames=TRUE)
  setnames(smeta,'rn','cell')
  fwrite(smeta,
         file=sprintf("%s_meta.tsv.gz",out_pref),
         sep="\t",
         col.names=TRUE,
         quote=FALSE)
}

seurat3_rename_to_sample_w_cluster <- function(sci) {
  sci$cluster_sgroup <- paste0(Idents(sci),'_',sci$sgroup)
  sci$cluster <- Idents(sci)
  Idents(sci) <- "cluster_sgroup"
  return(sci)
}

diff_from_integration <- function(sci,
                                  cluster_col="seurat_clusters",
                                  assay2="RNA",
                                  min_cell_cnts=3,
                                  ufeatures=NULL,
                                  min.pct=0.1,
                                  logfc.threshold = 0.25,
                                  method="LR",
                                  ncpu=3,
                                  debug=0) {
  
  diffMarker_dt <- NA
  
  if (debug==1){browser()}
  stopifnot(!invalid(method) & !invalid(method))
  
  if (!invalid(sci)) {
    
    DefaultAssay(sci) <- assay2
    
    cluster_labels <- sort(unique(sci@meta.data[[cluster_col]]))
    sci$cluster_sample <- paste0(sci@meta.data[[cluster_col]],'_',sci@meta.data$comp)
    
    Idents(sci) <- "cluster_sample"
    diffMarker2 <- lapply(cluster_labels,function(clab) {
      comp_label <- sprintf("cluster_label [%s]",clab)
      
      ctrl_lab <- sprintf("%s_ctrl",clab)
      expr_lab <- sprintf("%s_expr",clab)
      
      if (debug==1){browser()}
      
      C1=sum(sci$cluster_sample==ctrl_lab)
      E1=sum(sci$cluster_sample==expr_lab)
      
      message(sprintf("Performing DEG on %s ...",clab))
      diffMarker= NA
      if ((C1 >= min_cell_cnts) & (E1 >= min_cell_cnts)) {
        
        # message(sprintf("assay[%s];test.method[%s]",assay2,method))
        fm_result <- FindMarkers(sci,
                                 test.use=method,
                                 assay=assay2,
                                 features=ufeatures,
                                 ident.1=expr_lab,
                                 min.pct = min.pct,
                                 logfc.threshold = logfc.threshold,
                                 ident.2=ctrl_lab,
                                 only.pos=FALSE,
                                 verbose=FALSE) #NOTE that ident.1 is of our interest
        if (nrow(fm_result)>0) {
          diffMarker <- as.data.table(fm_result)
          diffMarker$cluster <- clab
          diffMarker$gexpr = NA
          diffMarker$gene <- rownames(fm_result)
        }
      } else if (C1>3 | E1>3) {
        message(sprintf('not enough samples for diff analysis [%s].',comp_label))
        
        if (C1 > E1) {
          sci2 = subset(sci,cluster_sample == ctrl_lab)
        } else {
          sci2 = subset(sci,cluster_sample == expr_lab)
        }
        avg_expr.dt = as.data.table(AverageExpression(sci2)[[assay2]],keep.rownames = T)
        avg_expr.dt = avg_expr.dt[all>1.,]
        avg_expr.dt = avg_expr.dt[order(-all),][1:round(dim(avg_expr.dt)[1] * 0.2),]
        avg_expr.dt = avg_expr.dt[!is.na(rn),]
        
        if (nrow(avg_expr.dt)>0) {
          data.m = GetAssayData(object = sci2, slot = "data")
          data.m = as.matrix(data.m[avg_expr.dt$rn,])
          data.m[data.m>0]=1.
          
          N = dim(data.m)[2]
          pct.1 = rowSums(data.m)/N
          
          diffMarker <- data.table(p_val=NA,
                                   avg_log2FC = NA,
                                   pct.1=rowSums(data.m)/N,
                                   pct.2=0,
                                   p_val_adj=NA,
                                   cluster=clab,
                                   gexpr=avg_expr.dt$all,
                                   gene=rownames(data.m))
        }
      }
      
      return(diffMarker)
    }
    # ,mc.cores = ncpu
    )
    
    if (debug==1){browser()}
    diffMarker2 <- diffMarker2[!is.na(diffMarker2)]
    diffMarker_dt <- rbindlist(diffMarker2)
  }
  
  if (debug==1){browser()}
  warnings()
  return(diffMarker_dt)
}

#' Perform DEsingle on Seurat Integrated Data
#' NOTE: this is under development. DGE postprocess (table/visualization) will be affected.
#'
#'    Changjin Hong, hongc2@ccf.org

#'    ref: DEsingle for detecting three types of differential expression in single-cell RNA-seq data
#'
#' @param sci seurat integrated object
#' @param cluster_col meta.data column to use for clustering/group
#' @param min_cell_cnts the min number of cells to avoid DGE comp error due to insufficient samples
#' @param ncpu The number of cpus to utilize
#' @param debug 0 (no debug), 1 (debug; test with the first sample only)
#' @return DGE report dt.table
#'
#'
DEsingle_from_integration <- function(sci,
                                      cluster_col="seurat_clusters",
                                      min_cell_cnts=3,
                                      ncpu=3,
                                      debug=0) {
  
  diffMarker_dt <- NA
  
  if (debug==1){browser()}
  
  if (!invalid(sci)) {
    param <- MulticoreParam(workers = ncpu, progressbar = TRUE)
    register(param)
    
    assay2 <- "SCT"
    DefaultAssay(sci) <- assay2
    cluster_labels <- sort(unlist(unique(sci[[cluster_col]]),use.names=F))
    
    if (invalid(cluster_col)) {
      stopifnot(length(unique(sci$comp)),2)
      sci$cluster_sample <- sci$comp
    } else {
      sci$cluster_sample <- paste0(unlist(sci[[cluster_col]],use.names=F),'_',sci$comp)
    }
    
    Idents(sci) <- "cluster_sample"
    
    diffMarker2 <- lapply(cluster_labels,function(clab) {
      comp_label <- sprintf("cluster_label [%s]",clab)
      
      ctrl_lab <- sprintf("%s_ctrl",clab)
      expr_lab <- sprintf("%s_expr",clab)
      
      if ((length(which(sci$cluster_sample==ctrl_lab)) >= min_cell_cnts) &
          (length(which(sci$cluster_sample==expr_lab)) >= min_cell_cnts)) {
        
        message(sprintf('enough samples for diff analysis [%s].',comp_label))
        
        message(sprintf("assay[%s]",assay2))
        cbc.ctrl <- Cells(sci)[sci$cluster_sample==ctrl_lab]
        cbc.expr <- Cells(sci)[sci$cluster_sample==expr_lab]
        
        group2 <- factor(c(rep(ctrl_lab,length(cbc.ctrl)),
                           rep(expr_lab,length(cbc.expr))))
        
        if (debug==1){browser()}
        
        sci.mtx <- as.matrix(GetAssayData(sci,slot = "counts"))[,c(cbc.ctrl,cbc.expr)]
        
        # Detecting the DE genes in parallelization with ncpu cores
        fm_result <- DEsingle(counts = sci.mtx, group = group2, parallel = TRUE, BPPARAM = param)
        rm(sci.mtx)
        fm_result.classified <- DEtype(results = fm_result, threshold = 0.05)
        
        diffMarker <- as.data.table(fm_result.classified)
        diffMarker$cluster_id <- clab
        diffMarker$gene <- rownames(fm_result.classified)
      } else {
        diffMarker <- NA
      }
      return(diffMarker)
    }
    )
    
    if (debug==1){browser()}
    
    diffMarker2 <- diffMarker2[!is.na(diffMarker2)]
    
    diffMarker_dt <- rbindlist(diffMarker2)
  }
  
  if (debug==1){browser()}
  
  warnings()
  return(diffMarker_dt)
}

ctrl_vs_expr_for_comp <- function(snames,comp_sheet,sheet_idx=1,ctrl_by_or=NA,expr_by_or=NA,debug2=0) {
  if (debug2==1) {browser()}
  
  if (!invalid(comp_sheet) & file.exists(comp_sheet)) {
    if (endsWith(comp_sheet,"xlsx")|endsWith(comp_sheet,"xls")) {
      message(sprintf("reading xlsx[%s],sheet_idx[%d]",comp_sheet,sheet_idx))
      comp2_dt <- as.data.table(read.xlsx(comp_sheet,sheet=sheet_idx))
    } else {
      comp2_dt <- fread(comp_sheet)
    }
    
    message(sprintf("column headers avail at %s: %s",comp_sheet,paste0(colnames(comp2_dt),collapse = ',')))
    if (any(snames %in% comp2_dt$sample)) {
      comp2_dt <- comp2_dt[(sample %in% snames),]
      sfield <- "sample"
    } else if (any(snames %in% comp2_dt$sample_tag)) {
      comp2_dt <- comp2_dt[(sample_tag %in% snames),]
      sfield <- "sample_tag"
    } else {
      stop("verify comp_sheet if either sample or sample_tag does match with seurat object ...")
    }
    
    stopifnot(nrow(comp2_dt[comp=="ctrl",])>0 & nrow(comp2_dt[comp=="expr",])>0)
    
    comp_dt <- data.table(ctrl=comp2_dt[comp=="ctrl",paste0(get(sfield),collapse=",")],
                          expr=comp2_dt[comp=="expr",paste0(get(sfield),collapse=",")])
  } else {
    comp_dt <- as.data.table(t(combn(snames,2,simplify = TRUE)))
    colnames(comp_dt) <- c("ctrl","expr")
  }
  if (debug2==1) {browser()}
  return(comp_dt)
}

batch_dge_from_seu_integs_objs <- function(sc.integs,
                                           comp_sheet="",
                                           cluster_col="seurat_clusters",
                                           dge_method="MAST",
                                           assay2 = "RNA",
                                           ncpu=3,
                                           sheet_idx=1,
                                           debug=0) {
  if (debug==1){browser()}
  
  diffMarker_dts <- list()
  sci_names <- names(sc.integs)
  N <- length(sci_names)
  for (gsi in 1:N) {
    gsample <- sci_names[[gsi]]
    # gsample <- names(sc.integs)[[2]]
    message(gsample)
    sci <- sc.integs[[gsample]]
    snames <- sort(unique(sci@meta.data$orig.ident))
    if (length(snames)>1) {
      comp_dt <- ctrl_vs_expr_for_comp(snames,comp_sheet,sheet_idx=sheet_idx,debug2 = debug)
      
      diffMarker_dt2 <- list()
      R <- dim(comp_dt)[1]
      for (r in 1:R) {
        # r <- 1
        job_notice <- sprintf("%s[%d/%d],comp[%d/%d]",gsample,gsi,N,r,R)
        message(job_notice)
        samples_to_comp <- unlist(apply(comp_dt[r,],2,comma_string_to_list),use.names = F)
        
        seu <- subset(sci,cells = Cells(sci)[sci@meta.data$orig.ident %in% samples_to_comp])
        seu@meta.data$comp <- "comp"
        seu@meta.data$comp[which(seu@meta.data$orig.ident %in% comma_string_to_list(comp_dt[r,'ctrl']))] <- "ctrl"
        seu@meta.data$comp[which(seu@meta.data$orig.ident %in% comma_string_to_list(comp_dt[r,'expr']))] <- "expr"
        
        if (dge_method=="desingle") {
          diffMarker_dt <- DEsingle_from_integration(seu,
                                                     cluster_col=cluster_col,
                                                     min_cell_cnts=3,
                                                     ncpu=ncpu,
                                                     debug=debug)
        } else {
          diffMarker_dt <- diff_from_integration(seu,
                                                 cluster_col=cluster_col,
                                                 min_cell_cnts=3,
                                                 method=dge_method,
                                                 assay2 = assay2,
                                                 ncpu=ncpu,
                                                 debug=debug)
        }
        if (!invalid(diffMarker_dt)) {
          diffMarker_dt$ctrl <- comp_dt[r,'ctrl']
          diffMarker_dt$expr <- comp_dt[r,'expr']
          diffMarker_dt$comp <- sprintf("comp_%d",r)
          diffMarker_dt$gsample <- gsample
          diffMarker_dt$group1 <- diffMarker_dt$gsample
          diffMarker_dt$group2 <- NA
          if (!isEmpty(grep("\\.",diffMarker_dt$gsample))) {
            diffMarker_dt$group1 <- tstrsplit(diffMarker_dt$gsample,"\\.")[[1]]
            diffMarker_dt$group2 <- tstrsplit(diffMarker_dt$gsample,"\\.")[[2]]
          }
          
          cmpname <- diffMarker_dt[1,paste0(ctrl,"_",expr)]
          diffMarker_dt2[[cmpname]] <- diffMarker_dt
        }
        message(sprintf("Done[%s]",job_notice))
      }
      diffMarker_dts[[gsample]] <- diffMarker_dt2
    }
  }
  return(diffMarker_dts)
}


findConservedMarkers_with_tryCatch <- function(sci,expr_ident,group_by,method="MAST",debug2=0,exp_tag="fcm") {
  ge_conserved_df <- tryCatch(
    {
      DefaultAssay(sci) <- "RNA"
      if (debug2==1){browser()}
      
      if (method=="desingle") {
        cbc.ctrl <- names(Idents(sci))[Idents(sci)!=expr_ident]
        cbc.expr <- names(Idents(sci))[Idents(sci)==expr_ident]
        
        group2 <- factor(c(rep(sprintf("no_%s",expr_ident),length(cbc.ctrl)),
                           rep(expr_ident,length(cbc.expr))))
        
        sci.mtx <- as.matrix(GetAssayData(sci,slot = "counts"))[,c(cbc.ctrl,cbc.expr)]
        
        # Detecting the DE genes in parallelization with ncpu cores
        fm_result <- DEsingle(counts = sci.mtx, group = group2)
        rm(sci.mtx)
        ge_conserved_df <- DEtype(results = fm_result, threshold = 0.05)
      } else {
        ge_conserved_df <-FindConservedMarkers(sci,
                                               test.use=method,
                                               ident.1 = expr_ident,
                                               grouping.var = group_by,
                                               verbose = FALSE)
      }
    },
    error=function(cond) {
      # browser() #debug
      message(sprintf('FindConservedMarkers(fcm=%s,expr_ident=%s,group_by=%s) is failed!',exp_tag,expr_ident,group_by))
      return(NA)
    },
    finally={
      # message('Done.')
    }
  )
  return(ge_conserved_df)
}

conserved_from_integration <- function(sc.integs,cluster_col="seurat_clusters",group_by="orig.ident",dge_method="MAST",ncpu=4,debug2=0) {
  
  # cons_rd_file <- sprintf("%s.rd",out_pref)
  
  conservedMarkers <- list()
  # for (scomp in c("CART.P12")) {
  for (scomp in names(sc.integs)) { #debug
    
    message(sprintf("scomp[%s]",scomp))
    if (debug2==1) {browser()}
    
    sci <- sc.integs[[scomp]]
    DefaultAssay(sci) <- "RNA"
    Idents(sci) <- cluster_col
    
    if (!invalid(sci)) {
      ident_cluster_df <- table(sci$orig.ident,unlist(sci[[cluster_col]]))
      if (dge_method=="desingle") {
        comm_clusters <- colnames(ident_cluster_df)[colSums(ident_cluster_df)>=10]
      } else {
        comm_clusters <- colnames(ident_cluster_df)[which(colSums(ident_cluster_df>2)>=2)]
      }
      
      cluster_markers <- mclapply(comm_clusters, function(identj) {
        message(identj)
        ret <-findConservedMarkers_with_tryCatch(sci,identj,group_by,method=dge_method,debug2=debug2,exp_tag=scomp)
        return(ret)}
        ,mc.cores=ncpu
      )
      
      names(cluster_markers) <- comm_clusters
      cluster_markers <- cluster_markers[!is.na(cluster_markers)]
      
      conMarkers <- lapply(names(cluster_markers),function(identj) {
        message(identj)
        cmarker <- cluster_markers[[identj]]
        if (dim(cmarker)[1]==0) {
          cmarker <- NA
        } else {
          cmarker$cluster_id<-identj
          cmarker$gene<-rownames(cmarker)
          cmarker <- as.data.table(cmarker)
        }
        return(cmarker)
      })
      names(conMarkers) <- names(cluster_markers)
      conMarkers <- conMarkers[!is.na(conMarkers)]
      
      if (debug2==1) {browser()}
      
      consMarker_dt <- rbindlist2(conMarkers)
      
      conservedMarkers[[scomp]] <- list(shared=consMarker_dt)
      
    }
  }
  if (debug2==1) {browser()}
  # save(conservedMarkers,file=cons_rd_file,compress=TRUE)
  return(conservedMarkers)
}
#' To get DGE tables from seurat3 SCT integration where DGE is defined from one cluster vs. all the other clusters in the same sample
#'
#' To obtain DGE tables from seurat3 SCT integration
#' more detail
#'
#'    Changjin Hong, hongc2@ccf.org

#'    ref:
#'
#' @param sc.integs a list of Seurat3 SCT integration objects

#' @return list of data.table of DGE.
#'
dge_within_sample_from_integration <- function(sc.integs) {
  
  dgei <- list()
  for (cmp_name in names(sc.integs)) {
    # browser()
    sci <- sc.integs[[cmp_name]]
    
    DefaultAssay(sci) <- "RNA"
    
    sci <- NormalizeData(sci, verbose = FALSE)
    
    if (!invalid(sci)) {
      Idents(object = sci) <- "orig.ident"
      
      dgei[[cmp_name]] <- mclapply(unique(Idents(sci)),function(sample) {
        
        message(sprintf('DGE analy from sc3 integ [%s,%s]',cmp_name,sample))
        
        sci_j <- subset(x = sci, idents = sample)
        Idents(sci_j) <- sci_j$seurat_clusters
        
        dge_j <- FindAllMarkers(object = sci_j,
                                only.pos = FALSE)
        return(dge_j)
      },mc.cores = length(unique(Idents(sci))))
      
      names(dgei[[cmp_name]]) <- unique(Idents(sci))
    }
  }
  
  return(dgei)
}

# --------------------
#' To annotate sc3 object integ clusters by cellassign
#'
#' To obtain ca_fits and ca_fits_matrix_info, refer to cellassign_pbmc.ipynb
#'
#'    Changjin Hong, hongc2@ccf.org

#'    ref:
#'
#' @param sc.integs a list of Seurat3 SCT integration objects
#' @param ca_fits a list of cellassign() output variables
#' @param ca_fits_matrix_info: a list of matrix row(genes)/col(cell_ids) information used in ca_fits
#'
#' @return list of Seurat3 integration objects where celltype is added into $meta.data
#'
#'
annotate_integ_cluster_by_cellassign <- function(sc.integs,annot_cell,ncpu=3,debug2=0) {
  #ca_fits,ca_fits_matrix_info
  #
  message("annotating a cell type predicted by cellassign ...")
  if (debug2==1) {browser()}
  sc.integs.annot <- lapply(sc.integs, function(sc.integ) {
    if (debug2==1) {browser()}
    sm <- sc.integ@meta.data
    sm$uuid <- paste0(sm$orig.ident,':',tstrsplit(rownames(sm),"_")[[1]])
    
    anns <- list()
    for (sname in unique(sm$orig.ident)) {
      anns[[sname]] <- data.table(uuid=paste0(sname,':',annot_cell$mtx_meta[[sname]]$cellid),
                                  cell_type=annot_cell$fits[[sname]]$cell_type)
    }
    ann <- rbindlist(anns)
    
    sm$cellassign <- ann[match(sm$uuid,ann$uuid),cell_type]
    sm$uuid <- NULL
    sc.integ@meta.data <- sm
    return(sc.integ)
  }) #,mc.cores = ncpu
  
  return(sc.integs.annot)
}


# --------------------
#' To annotate sc3 object integ clusters by cellassign
#'
#' To obtain ca_fits and ca_fits_matrix_info, refer to cellassign_pbmc.ipynb
#'
#'    Changjin Hong, hongc2@ccf.org

#'    ref:
#'    NOTE: using table(), the implementation can be simpler!
#'
#' @param sc.integ a Seurat3 SCT integration object
#' @param group_by2 column variable of sc.integ meta.data
#'
#' @return a data.frame of count matrix in the format of (group_by2,sample)
#'
cell_counts_by <- function(sc.integ,
                           out_prefix,
                           cell_id="seurat_clusters",
                           plot_title="cell_count",
                           debug=0) {
  
  # meta <- data.table(sample = sc.integ$orig.ident,
  #                    cid = sc.integ$seurat_clusters)
  if (debug==1) {browser()}
  meta <- data.table(sample = sc.integ$orig.ident,
                     ucell_id = sc.integ@meta.data[[cell_id]],
                     cid = unlist(sc.integ[["seurat_clusters"]]))
  
  meta_by_sample <- split(meta,by="sample")
  samples <- names(meta_by_sample)
  
  nn_sizes <- list()
  nn_size.pcts <- list()
  cellid_sizes <- list()
  cellid_size.pcts <- list()
  
  wb <- createWorkbook(plot_title)
  
  pdf_file <- sprintf('%s.pdf',out_prefix)
  message(sprintf("generating pdf_file[%s]",pdf_file))
  pdf(pdf_file)
  
  for (sname in samples) {
    a.ldf <- meta_by_sample[[sname]][,.N,by=c("cid","ucell_id")]
    a.wdf <- dcast(a.ldf,cid ~ ucell_id, value.var = 'N')
    a.mat <- as.matrix(a.wdf[,2:dim(a.wdf)[2]])
    rownames(a.mat) <- a.wdf[,'cid']
    a.mat[is.na(a.mat)] <- 0
    
    sheet_name <- sname
    addWorksheet(wb,sheetName = sheet_name)
    writeData(wb,sheet_name,a.mat,colNames=TRUE,rowNames=TRUE)
    
    nn_size <- data.table(sample=sname,
                          cid=rownames(a.mat),
                          val=rowSums2(a.mat))
    
    nn_size.pct <- data.table(sample=sname,
                              cid=rownames(a.mat),
                              val=100.*nn_size$val/sum(nn_size$val))
    
    cellid_size <- data.table(sample=sname,
                              cid=colnames(a.mat),
                              val=colSums2(a.mat))
    
    cellid_size.pct <- data.table(sample=sname,
                                  cid=colnames(a.mat),
                                  val=100.*cellid_size$val/sum(cellid_size$val))
    
    nn_sizes[[sname]] <- nn_size
    nn_size.pcts[[sname]] <- nn_size.pct
    cellid_sizes[[sname]] <- cellid_size
    cellid_size.pcts[[sname]] <- cellid_size.pct
    
    if (debug==1) {browser()}
    r1c2 <- 1
    a.norm_mat <- apply(a.mat,r1c2,function(vec) {
      S <- sum(vec)
      if (S>0) {
        norm_vec <- vec/S
      } else {
        norm_vec <- vec
      }
      return(norm_vec)
    })
    if (r1c2==1) {
      a.norm_mat <- t(a.norm_mat)
    }
    pheatmap(a.norm_mat,angle_col=45,main=sname,cluster_rows=FALSE,cluster_cols=FALSE)
  }
  if (debug==1) {browser()}
  
  p2 <- list()
  nn_sizes <- rbindlist(nn_sizes)
  p2[['nn_sizes']] <- ggplot(data=nn_sizes, aes(x=sample, y=val, fill=sample)) +
    geom_bar(stat="identity", position=position_dodge()) +
    ggtitle(sprintf("%s: total cells",plot_title)) +
    theme(axis.text.x=element_text(angle=45, hjust=1),legend.position="none") +
    xlab("sample") +
    ylab("cell counts")
  
  nn_size.wdf <- dcast(nn_sizes,cid ~ sample, value.var = 'val')
  sheet_name <- "nn_sizes"
  addWorksheet(wb,sheetName = sheet_name)
  writeData(wb,sheet_name,nn_size.wdf)
  
  # ----------
  nn_size.pcts <- rbindlist(nn_size.pcts)
  
  p2[['nn_size.pcts']] <- ggplot(data=nn_size.pcts, aes(x=cid, y=val, fill=sample)) +
    geom_bar(stat="identity", position=position_dodge()) +
    ggtitle(sprintf("%s: cell counts[%%]",plot_title)) +
    theme(axis.text.x=element_text(angle=45, hjust=1),legend.position="none") +
    xlab("sample") +
    ylab("cell counts [%]")
  
  nn_size.pct.wdf <- dcast(nn_size.pcts,cid ~ sample, value.var = 'val')
  sheet_name <- "nn_size.pcts"
  addWorksheet(wb,sheetName = sheet_name)
  writeData(wb,sheet_name,nn_size.pct.wdf)
  
  p <- CombinePlots(p2,ncol=1)
  plot(p)
  # -----------
  p2 <- list()
  cellid_sizes <- rbindlist(cellid_sizes)
  
  p2[['cellid_sizes']] <- ggplot(data=cellid_sizes, aes(x=cid, y=val, fill=sample)) +
    geom_bar(stat="identity", position=position_dodge()) +
    ggtitle(sprintf("%s: cell ids",plot_title)) +
    theme(axis.text.x=element_text(angle=45, hjust=1),legend.position="none") +
    xlab("cell counts") +
    ylab("cell ids") # + coord_flip()
  
  cellid_size.wdf <- dcast(cellid_sizes,cid ~ sample, value.var = 'val')
  sheet_name <- "cellid_sizes"
  addWorksheet(wb,sheetName = sheet_name)
  writeData(wb,sheet_name,cellid_size.wdf)
  # ------------
  cellid_size.pcts <- rbindlist(cellid_size.pcts)
  
  p2[['cellid_size.pcts']] <- ggplot(data=cellid_size.pcts, aes(x=cid, y=val, fill=sample)) +
    geom_bar(stat="identity", position=position_dodge()) +
    ggtitle(sprintf("%s: cell ids[%%]",plot_title)) +
    theme(axis.text.x=element_text(angle=45, hjust=1),legend.position="none") +
    xlab("cell counts") +
    ylab("cell ids [%]")# + coord_flip()
  
  cellid_size.pct.wdf <- dcast(cellid_size.pcts,cid ~ sample, value.var = 'val')
  sheet_name <- "cellid_size.pcts"
  addWorksheet(wb,sheetName = sheet_name)
  writeData(wb,sheet_name,cellid_size.pct.wdf)
  
  xlsx_file <- sprintf("%s.xlsx",out_prefix)
  message(sprintf("creating an excel file [%s]",xlsx_file))
  saveWorkbook(wb, xlsx_file, overwrite = TRUE)
  
  p <- CombinePlots(p2,ncol = 1)
  plot(p)
  dev.off()
  if (debug==1) {browser()}
  message(sprintf("Done [%s].",plot_title))
}

doublet_finder <- function(seuj,ncpu=1,debug2=0) {
  # browser()
  # Pre-process seurat object with standard seurat workflow
  if (debug2==1) {browser()}
  
  use_sct = FALSE
  if (DefaultAssay(seuj)=="SCT") {
    use_sct = TRUE
  }
  # pK identification (no ground-truth)
  min.pc = find_optimal_npca(seuj)
  
  # sweep.list <- paramSweep_v3(seuj, PCs = 1:min.pc, num.cores = detectCores() - 1)
  message('do [paramSweep_v3]')
  if (debug2==1){browser()}
  
  sweep.list <- paramSweep_v3(seuj, PCs = 1:min.pc, num.cores = ncpu, sct = use_sct)
  
  message('do [summarizeSweep]')
  sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
  
  message('do [find.pK]')
  bcmvn <- find.pK(sweep.stats)
  
  # Optimal pK is the max of the bimodality coefficent (BCmvn) distribution
  bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
  optimal.pk <- bcmvn.max$pK
  optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
  
  d_rate = get_doublet_rate_in_10x2(u_ncell=dim(seuj)[2])*0.01
  message(sprintf("d_rate[%g]",d_rate))
  
  ## Homotypic doublet proportion estimate
  annotations <- seuj@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp.poi <- round(d_rate * nrow(seuj@meta.data))
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
  
  # run DoubletFinder
  message(sprintf("do doubletFinder_v3 [min.pc=%d, pK=%g, nExp=%g]",min.pc, optimal.pk, nExp.poi.adj))
  seuj <- doubletFinder_v3(seu = seuj,
                           PCs = 1:min.pc,
                           pK = optimal.pk,
                           reuse.pANN = FALSE,
                           nExp = nExp.poi.adj,
                           sct = use_sct)
  metadata <- seuj@meta.data
  j=grep('DF.class',colnames(metadata))
  colnames(metadata)[j] <- "doublet_finder"
  seuj@meta.data <- metadata
  
  message(sprintf("before taking singlet [%d]",dim(seuj)[2]))
  seuj=subset(seuj, doublet_finder == "Singlet")
  message(sprintf("after taking singlet [%d]",dim(seuj)[2]))
  
  j=grep('pANN_',colnames(seuj@meta.data))
  seuj@meta.data[[j]] = NULL
  seuj[["doublet_finder"]] = NULL
  DietSeurat(seuj)
}


# replaced by doublet_finder
scrublet_annot_batch <- function(bc_matrix_dir,debug2=0) {
  if (debug2==1){browser()}
  mtx_file <- file.path(bc_matrix_dir,'matrix.mtx')
  if (file.exists(sprintf("%s.gz",mtx_file))) {
    cmd <- sprintf("gunzip -fc %s.gz > %s",mtx_file,mtx_file)
    system(cmd)
  } else {
    stop(sprintf("check if %s.gz exists",mtx_file))
  }
  
  fea_file <- file.path(bc_matrix_dir,'features.tsv')
  if (file.exists(sprintf("%s.gz",fea_file))) {
    cmd <- sprintf("gunzip -fc %s.gz > %s",fea_file,fea_file)
    system(cmd)
  }
  
  scrublet_annot <- annotate_doublets(mtx_fpath=normalizePath(mtx_file),
                                      feature_fpath=normalizePath(fea_file))
  
  if (file.exists(mtx_file)) {
    unlink(mtx_file)
  }
  
  if (file.exists(fea_file)) {
    unlink(fea_file)
  }
  
  return(as.numeric(scrublet_annot[[1]]))
}


get_doublet_rate_in_10x <- function(ncells) {
  # ref: https://www.biotech.wisc.edu/services/gec/services/rochelightcyclerservices
  
  # ~0.8% per 1,000 cells; ~1.6% per 2,000 cells; ~2.3% per 3,000 cells; ~3.1% per 4,000 cells; ~3.9% per 5,000 cells; etc.
  
  library(splines)
  x <- c(1000,  2000, 3000, 4000, 5000, 8000)
  y <- c(0.008,0.016,0.023,0.031,0.039,0.070)
  fit2 <- lm( y~ns(x, 3) )
  # plot(x,y, xlim=c(1000,8001), ylim=c(0,0.1))
  # xx <- seq(1000,8001, length.out=250)
  # lines(xx, predict(fit2, data.frame(x=xx)), col='orange')
  return(predict(fit2,data.frame(x=ncells)))
}

get_doublet_rate_in_10x2 <- function(u_ncell=6000) {
  # ref: https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_01_qc.html#Filter_genes
  d_rate.dt = data.table(multiple_rate=c(0.4,0.8,1.6,2.3,3.1,3.9,4.6,5.4,6.1,6.9,7.6),
                         n_cells=c(500,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000))
  
  a=d_rate.dt[n_cells<u_ncell,]
  if (nrow(a)>0){
    m_rate=a[nrow(a),multiple_rate]
  } else {
    m_rate=d_rate.dt[1,multiple_rate]
  }
  m_rate
}

rbindlist2 <- function(dt_list) {
  
  ucname <- sort(unique(unlist(lapply(dt_list,colnames))))
  head_col<-c("gene","minimump_p_val","max_pval","cluster_id")
  ucname <- c(head_col,setdiff(ucname,head_col))
  
  dts <- lapply(dt_list,function(dt){
    delta_col <- setdiff(ucname,colnames(dt))
    if (!isEmpty(delta_col)) {
      dt[,c(delta_col):=NA]
    }
    setcolorder(dt,ucname)
    return(dt)
  })
  return(rbindlist(dts))
}

solo_FindAllMarkers <- function(seu,ident2=NA,assay2="RNA",de_test="LR",only_pos=FALSE,min_pct=0.1,ncpu=4,rds_fpath=NA,debug2=0) {
  if (debug2==1){browser()}
  
  message(sprintf("performing findallmarkers on the cell group [%s;%s;%s]...",ident2,assay2,de_test))
  
  plan2(ncpu = ncpu)
  if (assay2 %in% Seurat::Assays(seu)) {
    my_assay <- assay2
  } else {
    my_assay <- "RNA"
  }
  
  if (!invalid(ident2)) {
    Idents(seu) <- ident2
  }
  seu.markers.full <- FindAllMarkers(object = seu,
                                     assay = my_assay,
                                     test.use = de_test,
                                     min.pct = min_pct,
                                     only.pos = only_pos)
  plan2()
  
  return(seu.markers.full)
}

batch_solo_FindAllMarkers <- function(seus,ident2=NA, assay2="SCT",only_pos=FALSE,min_pct=0.1,de_test="LR",ncpu=4,debug2=0) {
  if (debug2==1){browser()}
  snames <- names(seus)
  
  cluster_marker_dges <- lapply(snames, function(sname) {
    message(sprintf("running FindAllMarkers[%s];Idents[%s]",sname,ident2))
    
    assays_avail <- names(seus[[sname]]@assays)
    if (!(assay2 %in% assays_avail)) {
      if (("SCT" %in% assays_avail)) {
        assay2 <- "SCT"
      } else {
        if ("RNA" %in% assays_avail) {
          assay2 <- "RNA"
        } else {
          stop(sprintf("no %s,SCT,RNA is available from the seurat object",assay2))
        }
      }
    }
    if (debug2==1){browser()}
    ret <- solo_FindAllMarkers(seus[[sname]],ident2=ident2,assay2=assay2,only_pos=only_pos,min_pct=min_pct,de_test=de_test,ncpu=ncpu,debug2=debug2)
    message(sprintf("Done[%s]",sname))
    return(ret)
  }
  # ,mc.cores = ncpu
  )
  names(cluster_marker_dges) <- snames
  if (debug2==1){browser()}
  return(cluster_marker_dges)
}


cell_id_by_singler <- function(seu,ref.sce,pmid,label2id="label.fine",qc_label="labels",debug2=0) {
  if (debug2==1){browser()}
  DefaultAssay(seu) <- "RNA"
  query.sce <- as.SingleCellExperiment(seu)
  shared_feat <- intersect(rownames(ref.sce),rownames(query.sce))
  
  message(sprintf("ref feat# [%d]",length(rownames(ref.sce))))
  message(sprintf("test feat# [%d]",length(rownames(query.sce))))
  message(sprintf("shared feat# [%d]",length(shared_feat)))
  message("cell ID annotation in progress ...")
  
  
  pred.id <- SingleR(test=query.sce,
                     ref=ref.sce,
                     # method="single", #depreciated
                     labels=ref.sce[[label2id]])
  
  message("Done.")
  seu[[pmid]] <- unlist(pred.id[qc_label])
  if (debug2==1){browser()}
  return(list(seu=seu,pred.id=pred.id))
}


single_seurat_to_sce <- function(seu,assay="SCT") {
  
  DefaultAssay(seu) <- assay
  sce <- SingleCellExperiment(assays = list(counts = GetAssayData(seu,slot = "counts"),
                                            logcounts=GetAssayData(seu,slot = "data")),
                              colData = as.data.frame(seu@meta.data))
  
  if (!invalid(Reductions(seu))) {
    reducedDims(sce) <- SimpleList(PCA=seu@reductions$pca@cell.embeddings[,1:2],
                                   UMAP=seu@reductions$umap@cell.embeddings)
  }
  
  return(sce)
}

seurat3_to_monocle3 <- function(seurat,clid="seurat_clusters",debug2=0) {
  if (debug2==1){browser()}
  # part one, gene annotations
  
  # gene_annotation <- as.data.frame(rownames(seurat@reductions[["pca"]]@feature.loadings),
  #                                  row.names = rownames(seurat@reductions[["pca"]]@feature.loadings))
  gene_annotation <- as.data.frame(rownames(seurat),
                                   row.names = rownames(seurat))
  colnames(gene_annotation) <- "gene_short_name"
  
  # part two, cell information
  
  cell_metadata <- as.data.frame(seurat@assays[["RNA"]]@counts@Dimnames[[2]],
                                 row.names = seurat@assays[["RNA"]]@counts@Dimnames[[2]])
  
  colnames(cell_metadata) <- "barcode"
  
  cell_metadata <- cbind(cell_metadata,seurat@meta.data)
  
  # part three, counts sparse matrix
  New_matrix <- seurat@assays[["RNA"]]@counts
  #New_matrix <- New_matrix[rownames(seurat@reductions[["pca"]]@feature.loadings), ]
  New_matrix <- New_matrix[rownames(seurat), ]
  expression_matrix <- New_matrix
  
  
  ### Construct the basic cds object
  cds_from_seurat <- new_cell_data_set(expression_matrix,
                                       cell_metadata = cell_metadata,
                                       gene_metadata = gene_annotation)
  
  
  ### Construct and assign the made up partition
  
  recreate.partition <- c(rep(1, length(cds_from_seurat@colData@rownames)))
  names(recreate.partition) <- cds_from_seurat@colData@rownames
  recreate.partition <- as.factor(recreate.partition)
  
  cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition
  
  ### Assign the cluster info
  
  list_cluster <- seurat@meta.data[[clid]]
  names(list_cluster) <- seurat@assays[["RNA"]]@data@Dimnames[[2]]
  
  cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
  
  
  ### Could be a space-holder, but essentially fills out louvain parameters
  
  cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"
  
  
  ### Assign UMAP coordinate
  cds_from_seurat@int_colData@listData$reducedDims@listData[["UMAP"]] <-seurat@reductions[["umap"]]@cell.embeddings
  
  ### Assign feature loading for downstream module analysis
  
  cds_from_seurat@preprocess_aux$gene_loadings <- seurat@reductions[["pca"]]@feature.loadings
  
  
  ### Learn graph, this step usually takes a significant period of time for larger samples
  
  print("Learning graph, which can take a while depends on the sample")
  cds_from_seurat <- learn_graph(cds_from_seurat, use_partition = T)
  return(cds_from_seurat)
  
  
}
get_cluster_by_expr <- function(seu2,assay2,slot2="data",debug2=0) {
  if(debug2==1) {browser()}
  avgExpr <- AverageExpression(seu2,assay=assay2,slot = slot2)
  avgExpr[[assay2]]
}


evaluate_metrics_by_cellranger_recommend<-function() {
  lower_bnd=data.table(`Valid Barcodes`=75,
                       `Q30 Bases in RNA Read`=65,
                       `Estimated Number of Cells`=500,
                       `Fraction Reads in Cells`=70,
                       `Mean Reads per Cell`=20000,
                       `Reads Mapped Confidently to Transcriptome`=30,
                       `Reads Mapped Antisense to Gene`=0)
  
  upper_bnd=data.table(`Valid Barcodes`=100,
                       `Q30 Bases in RNA Read`=100,
                       `Estimated Number of Cells`=10000,
                       `Fraction Reads in Cells`=100,
                       `Mean Reads per Cell`=Inf,
                       `Reads Mapped Confidently to Transcriptome`=100,
                       `Reads Mapped Antisense to Gene`=10)
  
  a = rbind(lower_bnd,upper_bnd)
  a$type = c("lower","upper")
  a
}

cell_count_w_meta <- function(sc_raws,sample_meta=NA,qc_comment="N",debug2=0) {
  if (debug2==1) {browser()}
  cell_cnts <- as.data.table(sapply(sc_raws,function(seu) {dim(seu)[2]}),keep.rownames = T)
  colnames(cell_cnts) <- c("sample","cnt")
  
  if (!invalid(sample_meta)) {
    cell_cnts<-merge(cell_cnts,sample_meta,by="sample")
  }
  cell_cnts$qc <- qc_comment
  cell_cnts
}


scina_with_tryCatch <- function(seu,markers,allow_unknown,outd,debug2=0) {
  if (debug2==1){browser()}
  pred_cell_ids <- tryCatch(
    {
      pred_cell_ids <- SCINA(GetAssayData(seu,slot = "data"),
                             markers,
                             max_iter = 100,
                             convergence_n = 10,
                             convergence_rate = 0.999,
                             sensitivity_cutoff = 0.9,
                             rm_overlap=FALSE,
                             allow_unknown=allow_unknown,
                             log_file=file.path(outd,"SCINA.log"))
    },
    error=function(cond){
      message('scina failed in handling input!')
      return(NULL)
    },
    finally={
      message('scina processes input!')
    }
  )
  return(pred_cell_ids)
}

# ------------
split_seurat_by_sample <- function(seu,debug2=0) {
  if (debug2==1){browser()}
  items <- unique(seu@meta.data$orig.ident)
  seus <- lapply(items, function(item){
    subset(seu,orig.ident==`item`)
  })
  names(seus) <- items
  seus
}

# ----------------
get_pseudobulk <- function(seu,uq3=1) {
  counts.m <- GetAssayData(seu,slot = "counts")
  if (uq3==1) {
    uq_mean = get_limit_from_quantile(colSums(counts.m),target_r = 0.75)
    counts.m = round(uq_mean*sweep(counts.m,2,colSums(counts.m),FUN="/"))
  }
  bulk.m = t(Matrix.utils::aggregate.Matrix(t(counts.m),
                                            groupings = seu@meta.data$orig.ident, fun = "sum"))
  bulk.m
}

cluster_n_lowdim <- function(seu, red_method="pca", n_pc=30,fc.resolution=0.8,do_tsne=0) {
  message("--FindNeighbors ...")
  seu <- FindNeighbors(seu, reduction = red_method, dims = 1:n_pc)
  
  message(sprintf("--FindClusters [%g] ...",fc.resolution))
  seu <- FindClusters(seu, resolution = fc.resolution)
  
  # Dimensional reduction and plotting
  message(sprintf("--RunUMAP [%s] ...",red_method))
  seu <- RunUMAP(seu, dims = 1:n_pc, reduction = red_method)
  
  if (do_tsne==1) {
    message("--RunTSNE ...")
    seu <- RunTSNE(seu, dims = 1:n_pc, reduction = red_method)
  }
  seu
}

# ----------------
# take a single seu object and perform a clustering based on RNA assay
sc_cluster = function(seu,
                      uvar_genes=NA,
                      excl_pat=NA,
                      reuse_sct=0,
                      scale_to_regress=c("percent.mt","nCount_RNA"),
                      n_features=3000,
                      fc.resolution=0.8,
                      do_tsne=0,
                      default_assay="SCT",
                      D=0, #find an optimal PC #
                      debug2=0) {
  
  if (debug2==1){browser()}
  
  set.seed(101)
  
  if (DefaultAssay(seu) != "Protein") {
    if ('S.Score' %in% scale_to_regress){
      if (!('Phase' %in% colnames(seu@meta.data))) {
        assay_bkp = DefaultAssay(seu)
        
        if ('RNA' %in% names(seu@assays)){
          DefaultAssay(seu) = 'RNA'
        } else if ('Spatial' %in% names(seu@assays)) {
          DefaultAssay(seu) = 'Spatial'
        }
        seu <- cellcycle_scoring(seu,debug2=debug2)
        DefaultAssay(seu) = assay_bkp
      }
    }
  }
  
  if (default_assay=="SCT") {

    if (reuse_sct==0){
      seu <- SCTransform(seu, 
                         variable.features.n = n_features,
                         vars.to.regress = scale_to_regress, 
                         assay = DefaultAssay(seu), 
                         verbose = FALSE, 
                         method = "glmGamPoi")
    }
  } else {
    stopifnot(default_assay %in% names(seu@assays))
    DefaultAssay(seu) <- default_assay
    seu <- NormalizeData(seu)
    seu <- FindVariableFeatures(seu,selection.method = "vst", nfeatures = n_features)
    
    if (!invalid(uvar_genes)) {
      message("use genes assigned by user for variableFeatures ...")
      VariableFeatures(seu) = c(VariableFeatures(seu),uvar_genes[uvar_genes %in% rownames(seu)]) %>%
        unique()
    }
    
    if (!invalid(excl_pat)) {
      var_genes = VariableFeatures(seu)
      VariableFeatures(seu) = var_genes[grep(excl_pat,var_genes,invert = T)]
    }
    message("--ScaleData ...")
    seu <- regress_out_scaledata(seu,
                                 scale_to_regress=scale_to_regress)
  }
  
  do_cluster = 1
  if (do_cluster==1) {
    message("--RunPCA ...")
    seu <- RunPCA(seu,
                  npcs = find_max_npcs(seu),
                  features = VariableFeatures(object = seu))

  }
  
  # Find significant PCs
  if (D==0) {
    D = find_optimal_npca(seu)
  }
  message(sprintf("D=%d",D))
  red_method="pca"
  if (do_cluster==1) {
  seu = cluster_n_lowdim(seu = seu, 
                         red_method = red_method, 
                         do_tsne=do_tsne, 
                         n_pc=D,
                         fc.resolution = fc.resolution)
  }
  message("Done[sc_cluster].")
  if (debug2==1){browser()}
  return(seu)
}


# split seurat object by a meta.data column
split_seurats = function(seu,split_by="orig.ident",debug2=0) {
  if (debug2==1){browser()}
  
  stopifnot(split_by %in% colnames(seu@meta.data))
  if ('split_by' %in% colnames(seu@meta.data)) {
    seu@meta.data$split_by=NULL
  }
  
  seu@meta.data$split_by=seu@meta.data[[split_by]]
  items <- unique(seu@meta.data[[split_by]])
  
  seus <- lapply(items, function(item){
    seuj=subset(seu,split_by==`item`) %>%
      DietSeurat()
    seuj$split_by=NULL
    seuj
  })
  names(seus) <- items
  seus
}

gexpr_prof_by_cgroup <- function(seu, 
                                  min_expr_co = 0.,
                                  slot='data',
                                 features = NULL,
                                  group.by = 'orig.ident',
                                 debug2=0) {
  if (debug2==1) {browser()}
  message(sprintf('gene expression profiling on [%s]',group.by))
  cmeta = get_cmeta_from_seu(seu)
  stopifnot(group.by %in% colnames(cmeta))
  setnames(cmeta, group.by, 'cellgroup1234')
  seu[['cellgroup1234']] = cmeta$cellgroup1234
  group_cnt = cmeta[,.N,by=c("cellgroup1234")]
  group_cnt = group_cnt[N>0,]
  cgroups = group_cnt$cellgroup1234
  names(cgroups) = cgroups
  
  imap(cgroups, function(cgroup,dmy) {
    if (debug2==1) {browser()}
    seu.2=subset(seu,cellgroup1234==cgroup)
    
    if (invalid(features)) {
      data.m = as.matrix(GetAssayData(object = seu.2, slot = slot))
    } else {
      data.m = as.matrix(GetAssayData(object = seu.2, slot = slot)[features,])
    }
    
    # compute expression fraction
    expr_cell_frac = rowSums(data.m>min_expr_co)/dim(data.m)[2]
    
    avg_expr = apply(data.m,1,function(r2) {
      log1p(mean(expm1(r2)))
    })
    avg_expr[is.nan(avg_expr)] = 0
    
    avg_expr_nz = apply(data.m,1,function(r2) {
      log1p(mean(expm1(r2[r2>0])))
    })
    
    avg_expr_nz[is.nan(avg_expr_nz)] = 0
    
    expr_summary = data.table(gene = names(expr_cell_frac),
                              pct.1 = expr_cell_frac,
                              cell_count = group_cnt[cellgroup1234==cgroup,N],
                              avg_expr = avg_expr,
                              avg_expr_nz = avg_expr_nz,
                              cell_group = cgroup)
  }) %>% rbindlist()
}

# -----
# this function to be replaced by gexpr_prof_by_cgroup()
average_expression_by_cellgroup_origident <- function(seu0,
                                                      group.by,
                                                      slot="data",
                                                      ncpu=4,
                                                      min_expr_co=0.,
                                                      debug2=0){
  
  if (debug2==1){browser()}
  seu0@meta.data$cellgroup1234 <- seu0@meta.data[[group.by]]
  
  group_cnt = as.data.table(seu0@meta.data)[,.N,by=c("cellgroup1234","orig.ident")]
  group_cnt = group_cnt[N>0,]
  total_cell_groups = nrow(group_cnt)
  
  row_vectors = 1:total_cell_groups
  names(row_vectors) = row_vectors
  cmeta = get_cmeta_from_seu(seu0)
  
  # plan2(ncpu = ncpu)
  expr_summary.dt = imap(row_vectors,function(r,dmy){
    message(sprintf("%d/%d in progress... ",r,total_cell_groups))
    cgroupj = group_cnt[r,cellgroup1234]
    samplej = group_cnt[r,orig.ident]
    # if (r==7){browser()}
    
    expr_summary = NA
    check_if_any = cmeta %>%
      dplyr::filter(cellgroup1234==cgroupj & orig.ident==samplej)
    if (nrow(check_if_any)>0) {
    seu0j=subset(seu0,cellgroup1234==cgroupj & orig.ident==samplej)
    
    data.m = as.matrix(GetAssayData(object = seu0j, slot = slot))
    
    # compute expression fraction
    expr_cell_frac = rowSums(data.m>min_expr_co)/dim(data.m)[2]
    
    # mean(expm1(Transcript_exp))
    
    avg_expr = apply(data.m,1,function(r2) {
      # mean(r2[r2>0])
      log1p(mean(expm1(r2)))
    })
    avg_expr[is.nan(avg_expr)] = 0
    
    avg_expr_nz = apply(data.m,1,function(r2) {
      # mean(r2[r2>0])
      log1p(mean(expm1(r2[r2>0])))
    })
    
    avg_expr_nz[is.nan(avg_expr_nz)] = 0

    expr_summary = data.table(gene = names(expr_cell_frac),
                              pct.1 = expr_cell_frac,
                              cell_count = group_cnt[r,N],
                              avg_expr = avg_expr,
                              avg_expr_nz = avg_expr_nz,
                              cell_group = cgroupj,
                              orig.ident = samplej)
    }
    expr_summary
    
  }
  # ,mc.cores = ncpu
  ) %>% rbindlist2()
  # plan2()
  # seu0@meta.data$cellgroup1234 = NULL
  expr_summary.dt
}

merge_seus_norm_fvar <- function(seus,default_assay="RNA",uvar_genes=NA,excl_pat=NA,debug2=0) {
  seus <- lapply(seus,function(seu) {
    if ('SCT' %in% names(seu@assays)) {
      seu@assays[['SCT']] <- NULL
      seu@reductions <- list()
      seu@graphs <- list()
    }
    DietSeurat(seu)
  })
  
  message("running a cell alignment for data integration ...")
  seui <- merge_seurats(seus)
  rm(seus)
  
  # DefaultAssay(seui) <- "RNA"
  seui <- NormalizeData(seui)
  if (invalid(uvar_genes)) {
    seui <- FindVariableFeatures(seui,selection.method = "vst")
  } else {
    message("use genes assigned by user for variableFeatures ...")
    VariableFeatures(seui) <- uvar_genes[uvar_genes %in% rownames(seui)]
  }
  
  if (!invalid(excl_pat)) {
    var_genes = VariableFeatures(seui)
    VariableFeatures(seui) = var_genes[grep(excl_pat,var_genes,invert = T)]
  }
  seui
}

perform_seu_clustering.2 <- function(seu,
                                     norm.method="LogNormalize",
                                     fv.sel_method = "vst",
                                     fv.mean_cutoff= c(0.1, 8),
                                     fv.dispersion_cutoff = c(1,Inf),
                                     genes = NA,
                                     scale.vars_to_regress = NULL,
                                     pca.npcs = 30,
                                     fneigh.npcs = 30,
                                     rtsne.dim = 5,
                                     fc.resolution=0.8,
                                     do_tsne=0,
                                     seed=12345,
                                     debug2=0,
                                     red_method="pca") {
  
  set.seed(seed)
  
  seu <- NormalizeData(seu,normalization.method = norm.method)
  
  if (invalid(genes)){
    seu <- FindVariableFeatures(seu,
                                selection.method = fv.sel_method,
                                mean.cutoff = fv.mean_cutoff,
                                dispersion.cutoff = fv.dispersion_cutoff)
  } else {
    VariableFeatures(seu) <- genes
  }
  
  meta.cols = colnames(seu@meta.data)
  if ('G2M.Score' %in% scale.vars_to_regress) {
    if (!('G2M.Score' %in% meta.cols)) {
      seu <- cellcycle_scoring(seu)
    }
  }
  
  seu <- ScaleData(seu,vars.to.regress = scale.vars_to_regress)
  
  seu <- RunPCA(seu,npcs=find_max_npcs(seu))
  
  # Dimensional reduction and plotting
  message(sprintf("--RunUMAP [%s] ...",red_method))
  
  if (debug2==1){browser()}
  seu <- RunUMAP(seu, dims = 1:pca.npcs, reduction = red_method)
  
  if (do_tsne==1) {
    message("--RunTSNE ...")
    seu <- RunTSNE(seu,
                   dims = 1:rtsne.dim,
                   seed.use=seed,
                   reduction = red_method)
  }
  
  message("--FindNeighbors ...")
  seu <- FindNeighbors(seu, reduction = red_method, dims = 1:fneigh.npcs)
  
  message("--FindClusters ...")
  seu <- FindClusters(seu,
                      random.seed = seed,
                      resolution = fc.resolution)
  seu
}

# the main issue of this method: Error in asMethod(object) :
# Cholmod error 'problem too large' at file ../Core/cholmod_dense.c, line 102
######################
combat_on_seurats <- function(seus,
                              covars,
                              batch.col="orig.ident",
                              D=30,
                              min.cells=50,
                              fc.resolution=0.8,
                              scale_to_regress=c("percent.mt","S.Score","G2M.Score"), #"nCount_RNA"
                              uvar_genes = NA,
                              excl_pat="^IG[HJKL]",
                              do_cluster=1,
                              do_tsne=0,
                              debug2=0) {
  if (debug2==1){browser()}
  
  cmeta = colnames(seus[[1]]@meta.data)
  stopifnot(all(c(scale_to_regress,batch.col) %in% cmeta))
  
  
  message("dieting seurats ...")
  seui <- merge_seus_norm_fvar(seus,uvar_genes,excl_pat,debug2)
  rm(seus)
  
  #-----------
  # run combat
  smeta = seui@meta.data
  
  for(covar in covars) {
    smeta[[covar]] = as.factor(smeta[[covar]])
  }
  
  formula_str <- sprintf("~%s",paste0(covars,collapse=" + "))
  message(sprintf("formula=%s",formula_str))
  covar.m = model.matrix(as.formula(formula_str), data=smeta)
  if (debug2==1){browser()}
  
  #back up data slot
  data.before_combat <- rlang::duplicate(GetAssayData(seui,slot="data"))
  
  SetAssayData(seui,slot="data") = GetAssayData(seui,slot="data") %>%
    as.matrix() %>%
    as.data.frame() %>%
    sva::ComBat(dat=.,
                batch=seui@meta.data[[batch.col]],
                mod=covar.m,
                par.prior=TRUE,
                prior.plots=FALSE) %>%
    as.matrix() %>%
    Matrix()
  
  message("--ScaleData ...")
  seui <- regress_out_scaledata(seui,
                                scale_to_regress=scale_to_regress)
  
  if (do_cluster==1) {
    message("--RunPCA ...")
    seui <- RunPCA(seui,
                   npcs = find_max_npcs(seui),
                   features = VariableFeatures(object = seui))
  }
  
  red_method="pca"
  if (do_cluster==1) {
    seui <- perform_seu_clustering(seui,D,fc.resolution)
  }
  message("Done[combat].")
  if (debug2==1){browser()}
  
  SetAssayData(seui,slot="data") <- data.before_combat
  return(seui)
}


subsample_cells_same_c <- function(seu,comp_group="orig.ident",n_sample=0,debug2=0){
  if (debug2==1) {browser()}
  set.seed(678)
  smeta = get_cmeta_from_seu(seu)
  
  cgroup_cnt = smeta %>%
    dplyr::count(comp_v=get(comp_group))
  
  min_c = cgroup_cnt[,min(n)]
  if (n_sample==0 | (n_sample>=min_c)) {
    n_sample = min_c
  }
  
  sampled_cbc = imap(split(smeta,by=comp_group),function(by_cg,sname) {
    message(sprintf("sampling cells in %s",sname))
    N = dim(by_cg)[1]
    if (N>n_sample){
      cbc = sample(by_cg$cbc, size=n_sample, replace=FALSE)
      a=data.table(comp_v=sname,cbc=cbc)
    } else if (N==0) {
      a = NA
    } else {
      cbc = by_cg$cbc
      a=data.table(comp_v=sname, cbc=cbc)
    }
    a
  }) %>% `[`(!is.na(.)) %>% rbindlist()
  
  sampled_cbc = sampled_cbc[order(comp_v,cbc),cbc]
  seu=subset(seu,cells = sampled_cbc)
  seu
}

########
subsample_by_cellgroup <- function(seu,sample.by,comp_group="orig.ident",n_sample=500,debug2=0){
  if (debug2==1) {browser()}
  set.seed(678)
  
  smeta = as.data.table(seu@meta.data,keep.rownames=T)
  setnames(smeta,'rn','cbc')
  
  sampled_cbc = imap(split(smeta,by=comp_group),function(smetaj,sname){
    imap(split(smetaj,by=sample.by),function(by_cg,cg){
      message(sprintf("sampling cells in %s|%s",sname,cg))
      N = dim(by_cg)[1]
      if (N>n_sample){
        cbc = sample(by_cg$cbc,size=n_sample,replace=FALSE)
        a=data.table(orig.ident=sname,cbc=cbc)
      } else if (N==0) {
        a = NA
      } else {
        cbc = by_cg$cbc
        a=data.table(orig.ident=sname,cbc=cbc)
      }
      a
    }) %>% `[`(!is.na(.)) %>% rbindlist()
  }) %>% `[`(!is.na(.)) %>% rbindlist()
  
  sampled_cbc = sampled_cbc[order(orig.ident,cbc),cbc]
  seu=subset(seu,cells = sampled_cbc)
  seu
}

iseu_to_pseudobulk <- function(iseu,assay2="RNA") {
  
  stopifnot(assay2 %in% Seurat::Assays(iseu))
  stopifnot(nrow(GetAssayData(seu,assay=assay2,slot="counts"))>0)
  
  split_seurats(iseu) %>%
    imap(.,function(seu,sname) {
      GetAssayData(seu,assay=assay2,slot="counts") %>%
        rowSums()
    }) %>% do.call('cbind',.)
}

seus_to_pseudobulk <- function(seus,assay2="RNA",debug2=0) {
  if (debug2==1) {browser()}
  imap(seus,function(seu,sname){
    message(sname)
    stopifnot(assay2 %in% Seurat::Assays(seu))
    stopifnot(nrow(GetAssayData(seu,assay=assay2,slot="counts"))>0)
    GetAssayData(seu,assay=assay2,slot="counts") %>%
      rowSums()
  }) %>% do.call('cbind',.)
}

pseudobulk_with_replicate2 <- function(by_sample,sname,uassay="RNA",uslot="counts",agg_fun="sum",rep_prefix=".",R=3,uq3=0,debug2=0) {
  if (debug2==1) {browser()}
  counts.m <- as.matrix(GetAssayData(by_sample,assay=uassay,slot = uslot))
  counts.m[is.na(counts.m)] = 0
  cbc_on = colnames(counts.m)
  
  if (length(cbc_on) >= 2*R) {
    pseudo_repi <- (sample(1:R, length(cbc_on), replace=T))
  } else {
    pseudo_repi = rep(1:R,length(cbc_on))
    pseudo_repi = pseudo_repi[1:length(cbc_on)]
  }
  
  pbk = NA
  if (length(cbc_on)>0) {
    if (uq3==1) {
      uq_mean = get_limit_from_quantile(colSums(counts.m),target_r = 0.75)
      counts.m = round(uq_mean*sweep(counts.m,2,colSums(counts.m),FUN="/"))
      counts.m[is.nan(counts.m)] = 0
    }
    
    pbk = t(Matrix.utils::aggregate.Matrix(t(counts.m),
                                           groupings = pseudo_repi, fun = agg_fun))
    
    colnames(pbk) = sprintf("%s%sR%s",sname,rep_prefix,colnames(pbk))
  }
  pbk
}

# this is temp since hd5r was not installed and refer to Load10X_Spatial() instead later
Load10X_bcdirs_Spatial <- function(data.dir,
                                   bc.dir="filtered_feature_bc_matrix",
                                   assay = "Spatial",
                                   slice = "slice1",
                                   filter.matrix = TRUE,
                                   orig.ident = "sample2",
                                   debug2= 0,
                                   to.upper = FALSE, image = NULL, ...)
{
  if (length(x = data.dir) > 1) {
    warning("'Load10X_Spatial' accepts only one 'data.dir'",
            immediate. = TRUE)
    data.dir <- data.dir[1]
  }
  
  bc_dir = file.path(data.dir,bc.dir)
  message(bc_dir)
  stopifnot(file.exists(bc_dir))
  data = Read10X(bc_dir)
  
  has_fb <- 0
  data_fb = NA
  if (invalid(names(data))) {
    data_gex = data
  } else {
    has_fb = 1
    data_gex = data$`Gene Expression`
    data_fb = data$`Antibody Capture`
  }
  
  if (to.upper) {
    rownames(x = data_gex) <- toupper(x = rownames(x = data_gex))
  }
  
  object = CreateSeuratObject(counts = data_gex,
                              project = orig.ident,
                              assay = "RNA",
                              min.cells = 0,
                              min.features = 0)
  
  if (has_fb==1) {
    message("Storing Protein Feature Barcodes ...")
    object[['Protein']] = CreateAssayObject(counts = data_fb)
  }
  
  if (debug2==1){browser()}
  
  if (invalid(x = image)) {
    spatial_d = file.path(data.dir,"spatial")
    
    if (file.exists(file.path(spatial_d,'scalefactors_json.json'))) {
      image <- Read10X_Image(image.dir = spatial_d,
                             filter.matrix = filter.matrix)
      image <- image[Cells(x = object)]
      DefaultAssay(object = image) <- assay
      object[[slice]] <- image
    } else {
      image=NULL
    }
    
  } else {
    if (!inherits(x = image, what = "VisiumV1")) {
      stop("Image must be an object of class 'VisiumV1'.")
    }
    image <- image[Cells(x = object)]
    DefaultAssay(object = image) <- assay
    object[[slice]] <- image
  }
  
  return(object)
}

Read10X_Image_hd <- function() {
  
  image <- readPNG(source = file.path(image.dir, "tissue_lowres_image.png"))
  
  scale.factors <- fromJSON(txt = file.path(image.dir, "scalefactors_json.json"))
  
  tissue.positions.path <- Sys.glob(paths = file.path(image.dir,
                                                      "tissue_positions*"))
  
  tissue.positions <- read.csv(file = tissue.positions.path, 
                               col.names = c("barcodes", "tissue", "row", "col", "imagerow", 
                                             "imagecol"), header = ifelse(test = basename(tissue.positions.path) == "tissue_positions.csv", yes = TRUE, no = FALSE), 
                               as.is = TRUE, row.names = 1)
  if (filter.matrix) {
    tissue.positions <- tissue.positions[which(x = tissue.positions$tissue == 
                                                 1), , drop = FALSE]
  }
  unnormalized.radius <- scale.factors$fiducial_diameter_fullres * 
    scale.factors$tissue_lowres_scalef
  spot.radius <- unnormalized.radius/max(dim(x = image))
  return(new(Class = "VisiumV1", 
             image = image, 
             scale.factors = scalefactors(spot = scale.factors$spot_diameter_fullres,
                                          fiducial = scale.factors$fiducial_diameter_fullres, 
                                          hires = scale.factors$tissue_hires_scalef, scale.factors$tissue_lowres_scalef), 
             coordinates = tissue.positions, spot.radius = spot.radius))
}

get_cmeta_from_seu <- function(seu) {
  seu@meta.data %>%
    as.data.table(keep.rownames = T) %>%
    rename(cbc="rn")
}

# compute gene enrichment score for seurat (JASMIN)
# ##############
get_jasmine_on_seurat <- function(seu,gset_lst,uassay="RNA",uslot="data",jasmin_method="likelihood") {
  gxs.mtx = as.matrix(GetAssayData(seu,assay=uassay,slot=uslot))
  compute_jasmin_on_matrix(gxs.mtx,gset_lst,jasmin_method)
}

# compute avg expression and pct expressed from seurat object
# ##############
DotPlots_dt <- function(seuj,cgroup="seurat_clusters",uassay="RNA",expr_co=0.,debug2=0){
  
  if (debug2==1){browser()}
  
  data.m = GetAssayData(object = seuj, assay=uassay, slot = "count") %>%
    as.matrix()
  # data.m[is.na(data.m)]=0
  cgroup_ids = seuj@meta.data[[cgroup]] %>% as.character()
  names(cgroup_ids) = rownames(seuj@meta.data)
  cgroup_ids[is.na(cgroup_ids)]="not_determined"
  
  avg_expr.dt = t(Matrix.utils::aggregate.Matrix(t(data.m),
                                                 groupings = cgroup_ids, fun = "sum")) %>%
    get_seurat_norm() %>%
    melt(varnames = c("gene","cell_group"),value.name = "mean") %>%
    as.data.table()
  
  # ----
  data.m[data.m>expr_co]=1.
  cnt_by_cgroup = table(cgroup_ids)
  expr_cnt.m = t(Matrix.utils::aggregate.Matrix(t(data.m),
                                                groupings = cgroup_ids, fun = "sum"))
  
  cnt_by_cgroup = cnt_by_cgroup[match(colnames(expr_cnt.m),names(cnt_by_cgroup))]
  
  pct.1 = sweep(expr_cnt.m, 2, cnt_by_cgroup, FUN = "/") %>%
    as.matrix()%>%
    melt(varnames = c("gene","cell_group"),value.name = "pct.1") %>%
    as.data.table()
  
  avg_expr.dt = merge(avg_expr.dt,pct.1,by=c("gene","cell_group"))
  avg_expr.dt
}

avg_expr_w_repl.1 <- function(seuj,sname,uassay="RNA",uslot="data",group.by="seurat_clusters",rep_prefix="#",debug2=0) {
  if (debug2==1) {browser()}
  if (uslot=="data"){
    # message("processing AE ...")
    seuj = AverageExpression(object=seuj,assays = "RNA",slot = "data",group.by = group.by, return.seurat = T)
    
    # message("processing GAD ...")
    pbk <- GetAssayData(seuj,assay="RNA",slot="data")
  } else {
    pbk <- GetAssayData(seuj,assay=uassay,slot = "counts")
    if (uq3==1) {
      uq_mean = get_limit_from_quantile(colSums(counts.m),target_r = 0.75)
      counts.m = round(uq_mean*sweep(counts.m,2,colSums(counts.m),FUN="/"))
    }
    pbk = t(Matrix.utils::aggregate.Matrix(t(pbk),
                                           groupings = seuj@meta.data[[group.by]], fun = "sum")) %>%
      as.matrix()
  }
  colnames(pbk) = sprintf("%s%s%s",sname,rep_prefix,colnames(pbk))
  pbk[sort(rownames(pbk)),]
}


# TODO: to be replaced by avg_expr_w_repl.1()
avg_expr_w_repl <- function(seuj,sname,uassay="RNA",uslot="data",rep_prefix=".",R=3,uq3=0,min_cells=2,debug2=0) {
  if (debug2==1) {browser()}
  cbc_on = Cells(seuj)
  if (length(cbc_on) >= 2*R) {
    set.seed(101)
    pseudo_repi <- (sample(1:R, length(cbc_on), replace=T))
  } else {
    pseudo_repi = rep(1:R,length(cbc_on))
    pseudo_repi = pseudo_repi[1:length(cbc_on)]
  }
  
  pbk = NA
  if (length(cbc_on)>=min_cells) {
    message(sprintf("reducing the single-cell gene expression matrix[%s] to #[%d] replicate in slot[%s] ...",sname,R,uslot))
    seuj=AddMetaData(seuj,metadata=pseudo_repi,col.name = "replicate")
    if (uslot=="data"){
      seuj = AverageExpression(object=seuj,assays = "RNA",slot = "data",group.by = "replicate",return.seurat = T)
      pbk <- GetAssayData(seuj,assay="RNA",slot="data")
    } else {
      pbk <- GetAssayData(seuj,assay=uassay,slot = "counts")
      if (uq3==1) {
        uq_mean = get_limit_from_quantile(colSums(counts.m),target_r = 0.75)
        counts.m = round(uq_mean*sweep(counts.m,2,colSums(counts.m),FUN="/"))
      }
      pbk = t(Matrix.utils::aggregate.Matrix(t(pbk),
                                             groupings = seuj@meta.data$replicate, fun = "sum")) %>%
        as.matrix()
    }
    
    if (R>1) {
      colnames(pbk) = sprintf("%s%sR%s",sname,rep_prefix,colnames(pbk))
    } else {
      colnames(pbk) = sname
    }
    pbk = pbk[sort(rownames(pbk)),]
  }
  pbk
}

get_seurat_meta_dt <- function(seu) {
  seu@meta.data %>%
    as.data.table(keep.rownames = T) %>%
    rename(cbc="rn")
}

Cluster_Highlight_Plot.2 <- function(seurat_object, cluster_name, highlight_color = "navy",
                                     background_color = "lightgray", pt.size = NULL, raster = NULL, ...) {
  
  # Is_Seurat(seurat_object = seurat_object)
  
  raster <- raster %||% (length(x = colnames(x = seurat_object)) >
                           2e+05)
  
  cells_to_highlight = NA
  if (any(cluster_name %in% Idents(seurat_object))) {
    cells_to_highlight <- CellsByIdentities(seurat_object, idents = cluster_name)
    if (invalid(x = pt.size)) {
      pt.size <- AutoPointSize_scCustom.2(data = sum(lengths(cells_to_highlight)),
                                          raster = raster)
    }
  }
  
  plot <- DimPlot(object = seurat_object, cells.highlight = cells_to_highlight,
                  cols.highlight = highlight_color, cols = background_color,
                  sizes.highlight = pt.size, pt.size = pt.size, order = TRUE,
                  raster = raster, ...)
  return(plot)
}

AutoPointSize_scCustom.2 = function(data, raster = NULL) {
  if (invalid(x = nrow(x = data)) && length(x = data) == 1 &&
      is.numeric(x = data)) {
    return(ifelse(test = isTRUE(x = raster), yes = 1, no = min(1583/data,
                                                               1)))
  }
  else {
    return(ifelse(test = isTRUE(x = raster), yes = 1, no = min(1583/nrow(x = data),
                                                               1)))
  }
}


# task: Want to viz two gene expression status where gene1 is only activated, gene2 is only activated, and both are activated
FeaturePlot_2genes <- function(seu,two_genes,co=0.5,debug2=0) {
  if (debug2==1) {browser()}
  gene1=two_genes[1]
  gene2=two_genes[2]
  
  cmeta = get_cmeta_from_seu(seu)
  gxs.m = GetAssayData(seu,slot="data")
  stopifnot(all(two_genes %in% rownames(gxs.m)))
  cmeta[,gene1.tmp:=gxs.m[gene1,]]
  cmeta[,gene2.tmp:=gxs.m[gene2,]]
  cmeta[,activation:="not_expressed"]
  
  cmeta[gene1.tmp>=co,activation:=gene1]
  cmeta[gene2.tmp>=co,activation:=gene2]
  cmeta[gene1.tmp>=co & gene2.tmp>=co, activation:="Both"]
  
  seu=AddMetaData(seu,metadata=cmeta$activation,col.name = "activation")
  
  Idents(seu) = "seurat_clusters"
  p = Meta_Highlight_Plot(seurat_object = seu, meta_data_column = "activation", pt.size = 0.01,
                          meta_data_highlight = c(two_genes,"Both"),
                          highlight_color = c("orange", "skyblue", "black"),
                          background_color = "lightgray")
  
  p = p + theme(legend.position = c(0.1, 0.9))
  
  p
}

get_density_of_feat <- function(seu.2, u.gene, assay="RNA") {
  umap.m = Embeddings(seu.2,reduction = "umap")
  DefaultAssay(seu.2) = assay
  v=Nebulosa:::calculate_density(w = GetAssayData(seu.2,slot="data")[u.gene,], x = umap.m, method = "wkde")
  names(v) = colnames(seu.2)
  v
}


scsorter_core <- function(seu, 
                          cellid.dt, 
                          ref_info="pmid000000",
                          default_assay=NA,
                          var_genes=NULL,
                          ret_fmt="seu",
                          debug2=0) {
  if (debug2==1){browser()}
  
  stopifnot(!invalid(ref_info))
  stopifnot(all(c('Type','Marker') %in% colnames(cellid.dt)))
  
  if (invalid(default_assay)) {
    default_assay = DefaultAssay(seu) 
  } else {
    DefaultAssay(seu) = default_assay
  }
  
  sname = seu@meta.data$orig.ident[1]
  cellid_avail.dt = cellid.dt[cellid.dt$Marker %in% rownames(seu),]
  genes_avail = unique(cellid_avail.dt$Marker)
  
  if(invalid(var_genes)) {
    if (invalid(VariableFeatures(seu))) {
      seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000, verbose = F)
    }
    var_genes=VariableFeatures(seu)
  }
  
  topgenes <- c(var_genes,genes_avail) %>%
    unique()
  
  sname = seu@meta.data$orig.ident[1]
  
  expr = GetAssayData(seu,slot = "data")
  topgene_filter = rowSums(as.matrix(expr)[topgenes,]!=0) > ncol(expr)*.1
  topgenes = topgenes[topgene_filter]
  
  picked_genes = unique(c(cellid_avail.dt$Marker, topgenes))
  expr = expr[rownames(expr) %in% picked_genes,]
  
  message(sprintf('running scSorter [%s] ...',sname))
  rts <- scSorter(expr, cellid_avail.dt)
  rm(expr)
  
  # #############
  if (debug2==1){browser()}
  seu[[ref_info]] = rts$Pred_Type
  seu
}

get_cellgroup_colors <- function(cell_groups) {
  
  # stopifnot(cell_group %in% colnames(seu@meta.data))
  # cell_groups = seu@meta.data[[cell_group]] %>%
  #   unique()
  
  S = length(cell_groups)
  
  cg_colors=scCustomize_Palette(num_groups = S, ggplot_default_colors = FALSE, color_seed = 123)
  
  cg_colors = cg_colors[1:S]
  names(cg_colors) = cell_groups
  cg_colors
}


assign_cellgroup_colors <- function(seu,cg_colors,group_by="seurat_clusters") {
  cg_names = names(cg_colors)
  cg_names = cg_names[cg_names %in% seu@meta.data[[group_by]]]
  cg_colors.2 = cg_colors[cg_names]
  seu@meta.data[[group_by]] = factor(seu@meta.data[[group_by]],levels=cg_names)
  Idents(seu) = group_by
  list(seu=seu,cg_colors=cg_colors.2)
}

AddModuleScore_with_tryCatch <- function(seu, gset_lst, n_bin=24) {
  seu.2 <- tryCatch(
    {
      seu.2 <- AddModuleScore(seu, 
                              features = gset_lst, 
                              name = "gs_prefix", 
                              nbin = n_bin)
    },
    error=function(cond){
      message(sprintf('n_bin[%d]',n_bin))
      return(NULL)
    },
    finally={
      message(sprintf('Done[%d]',n_bin))
    }
  )
  return(seu.2)
}

addmodulescore_in_dt <- function(seu, gset_lst,group_by="seurat_clusters",max_iter=10,debug2=0) {
  if (debug2==1){browser()}
  cmeta2.dt = get_cmeta_from_seu(seu)
  Idents(seu) = group_by
  
  n_bin=24
  iter=0
  was_successful=0

  gset_lst.cvd = imap(gset_lst,function(genes,gs_name){
    genes[genes %in% rownames(seu)]
  })
  while(iter<max_iter & was_successful==0) {
    seu.t2 <- AddModuleScore_with_tryCatch(seu, 
                                           gset_lst = gset_lst.cvd, 
                                           n_bin = n_bin)
    if (!invalid(seu.t2)) {
      was_successful = 1
    }
    iter = iter + 1
    n_bin = n_bin - 1
  }
  
  if (debug2==1){browser()}
  # dplyr::select(cbc,orig.ident,)
  cmeta.2 = NULL
  if (!invalid(seu.t2)){
    colnames(seu.t2@meta.data)[grep("gs_prefix", colnames(seu.t2@meta.data))] <- names(gset_lst.cvd)
    cmeta.2=get_cmeta_from_seu(seu.t2) %>%
      dplyr::select(cbc, names(gset_lst.cvd)) %>%
      melt(id.vars = 'cbc',variable.name = 'module_name', value.name = 'escore') %>%
      merge(x=.,y=cmeta2.dt,by="cbc") %>%
      data.table()
  }
  cmeta.2
}

# compute addmodulescore and extract the score data.table
AddModuleScore_Assay <- function(seu,
                                 gset_lst,
                                 max_iter=20,
                                 cell_group_by="seurat_clusters",
                                 new_assay_name="user_gset",
                                 debug2=0) {
  if (debug2==1){browser()}
  cvd_lst_fmt = imap(gset_lst,function(genes,gs_name) {
    genes[genes %in% rownames(seu)]
  })
  
  escore.m=addmodulescore_in_dt(seu,
                                gset_lst=cvd_lst_fmt,
                                debug2=debug2,
                                max_iter=max_iter,
                                group_by = cell_group_by) %>%
    dplyr::select(cbc,module_name,escore)
  
  if (invalid(escore.m)) {
    seu = NULL
  } else {
    escore.m = escore.m %>%
      dplyr::select(cbc,module_name,escore) %>%
      acast(module_name ~ cbc, value.var = 'escore')
    
    escore.m=escore.m[,colnames(seu)]
    rownames(escore.m) = gsub('_','-',rownames(escore.m))
    seu[[new_assay_name]] = CreateAssayObject(data=escore.m)
    seu = ScaleData(seu, assay=new_assay_name)
  }
  seu
}

compute_cell_frac <- function(cmeta,
                              sample_col="orig.ident",
                              cell_group_cols=c("seurat_clusters"),
                              debug2=0) {
  if (debug2==1){browser()}
  
  stopifnot(all(c(cell_group_cols,sample_col) %in% colnames(cmeta)))
  
  total_n = cmeta[,.(total_n=.N),by=sample_col]
  cell_n = cmeta[,.(cnt=.N),by=c(sample_col,cell_group_cols)]
  cell_n = merge(total_n,cell_n,by=sample_col)
  cell_n[,frac:=cnt/total_n] %>%
    as.data.table()
}

RegionNeighbors_out_nhop <- function(se, cell_group_column, query_label,nb_column=NA,max_hop=3,debug2=0) {
  
  if (debug2==1){browser()}
  if (invalid(nb_column)) {
    nb_column = sprintf("outer.%d_from_%s",max_hop,query_label)
  }
  
  se.2=rlang::duplicate(se)
  
  cmeta.2 <- get_cmeta_from_seu(se.2)
  # check sanity
  stopifnot(cell_group_column %in% colnames(cmeta.2))
  stopifnot(any(query_label %in% cmeta.2[[cell_group_column]]))
  
  cmeta.2$from_spot="else"
  cmeta.2[get(cell_group_column)==query_label,from_spot:="nb"]
  se.2[['from_spot']] = cmeta.2$from_spot
  cmeta.2 <- get_cmeta_from_seu(se.2)
  nb_lst = list()
  
  for(n_hop in 1:max_hop) {
    
    message(n_hop)
    
    se.2 <- RegionNeighbors(se.2, column_name ="from_spot", column_labels="nb",
                            column_key="aio_",mode = "all_inner_outer")
    
    cmeta.1 = get_cmeta_from_seu(se.2)
    
    cmeta.1[is.na(aio_nb),aio_nb:="else"]
    
    cmeta.1[,n_hop:=0]
    cmeta.1[aio_nb=="aio_nb",n_hop:=1]
    
    cmeta.1[,.N,by="aio_nb"]
    
    nb_lst[[n_hop]] = cmeta.1 %>%
      dplyr::select(cbc,aio_nb,n_hop)
    
    cmeta.1[,from_spot:="else"]
    cmeta.1[aio_nb!="else",from_spot:="nb"]
    se.2[["aio_nb"]]=NULL
    se.2[["from_spot"]]=cmeta.1$from_spot
  }
  se.2[["from_spot"]]=NULL
  
  nhop.m = imap(nb_lst,function(nb.dt,n_hop) {
    nb.dt %>%
      pull(n_hop)
  }) %>% do.call('cbind',.)
  rownames(nhop.m) = nb_lst[[1]]$cbc
  
  nhop.m2[is.na(nhop.m)]=0
  boundary_flag=rowSums(nhop.m)>0
  
  cmeta.2 <- get_cmeta_from_seu(se.2)
  cmeta.2[,max_hop_n:=0]
  cmeta.2[boundary_flag,max_hop_n:=max_hop]
  
  if (debug2==1){browser()}
  cmeta.2[get(cell_group_column)==query_label,max_hop_n:=0]
  
  se.2[[nb_column]] = cmeta.2$max_hop_n
  list(se=se.2, nb_lst=nb_lst)
}

seus_to_qced <- function(seus, outd=".", max_mt_pct=15, max_hb_pct=0, min_rb_pct=0, min.features=200, max.features=Inf, min.cells=3, assay="SCT", doublet_finder=1,exptag="initial",ncpu=4,debug2=0) {
  
  create_dir_if_not_exist(outd)
  plan2(ncpu=ncpu)
  seus = furrr::future_imap(seus,function(seu,sname) {
    message(sprintf("--- computing fraction of mt,ribo,hb [%s]",sname))
    seu <- PercentageFeatureSet(seu, pattern = "^MT-",col.name = "percent.mt") %>%
      PercentageFeatureSet(pattern="^RP[SL]", col.name = "percent.ribo") %>%
      PercentageFeatureSet(pattern="^HB[^(P)]", col.name = "percent.hb") %>%
      cellcycle_scoring()
    seu
  })
  
  # plot/table before QC
  # #########
  qc_comment="before.QC"
  pdf_fn = file.path(outd,"before_qc.pdf")
  print_qc_figures(seus,pdf_file = pdf_fn,qc_comment=qc_comment,debug2=0)
  
  qc_report_tsv <- file.path(outd,"cell_count_beforeQC.tsv")
  report_sc3_dim(seus,tsv_fpath=qc_report_tsv,qc_comment=qc_comment,append2=FALSE,debug2=0)
  
  # perform QC on mt/hb/rb/min.features/min.cells
  # #########
  seus = imap(seus,function(seu,sname) {
    message(sprintf("--- QC [%s]",sname))
    simple_qc_seu(seu,
                  max_mt_pct=max_mt_pct,
                  max_hb_pct=max_hb_pct,
                  min_rb_pct=min_rb_pct,
                  min.features=min.features,
                  max.features=max.features,
                  min.cells=min.cells,
                  debug2=0)
  })
  
  seus = imap(seus,function(seu,sname) {
    message(sprintf("--- normalizing [%s@%s]",sname,assay))
    if (assay=="SCT") {
      seu <- SCTransform(seu, vars.to.regress = "percent.mt", verbose = FALSE)
    } else {
      seu <- NormalizeData(seu) %>%
        FindVariableFeatures(selection.method="vst") %>%
        regress_out_scaledata()
    }
    seu
  })
  plan2()
  
  # # perform QC on doublets
  # #########
  if (doublet_finder==1) {
    seus = imap(seus,function(seu,sname) {
      message(sprintf("--- doublet finders [%s]",sname))
      seu=doublet_finder(seu,ncpu=1,debug2=0)
      seu
    })
  }
  plan2()
  
  message(sprintf("saving RDS file for qced seurat objects ..."))
  rds_fn = file.path(outd,"seus.rds")
  saveRDS(seus,file=rds_fn)
  
  if (F) {
    seus=imap(seus,function(seu,sname) {rename_orig_ident_seu(seu)})
    rename.map = rename_orig_ident()
    names(seus)=rename.map[match(names(seus),old_sname),new_sname]
    
    #restore slide image name
    seus=imap(seus,function(seu,sname) {
      names(seu@images) = sname
      seu
    })
    
    rds_file <- file.path(outd,sprintf("seus%s.rds",exptag))
    
    message(sprintf("saving RDS file [%s] ...",rds_file))
    saveRDS(sc_raws,file=rds_file)
  }
  
  # plot/table before QC
  # #########
  qc_comment=sprintf("mc%d_mf%d_mtr%g_hbr%g_rbr%g_dblt%d",min.cells,min.features,max_mt_r,max_hb_r,min_rb_r,doublet_filter)
  pdf_fn = file.path(outd,"after_qc.pdf")
  print_qc_figures(seus,pdf_fn,qc_comment=qc_comment,debug2=0)
  
  qc_report_tsv <- file.path(outd,"cell_count_afterQC.tsv")
  report_sc3_dim(seus,tsv_fpath=qc_report_tsv,qc_comment=qc_comment,append2=FALSE,debug2=0)
  seus
}

print_qc_figures <- function(cmeta,
                             pdf_fn = "./seus_qced_plots.pdf",
                             numerical_cols = c("nCount_RNA","nFeature_RNA","percent.mt","percent.ribo","percent.hb"),
                             categorical_cols = c('Phase'),
                             num_pws_cols = list(c('nCount_RNA','nFeature_RNA')),
                             qc_comment="before.QC",
                             debug2=0) {
  
  # ---------------------
  if (debug2==1){browser()}
  message("gathering single cell read/gene metrics ...")
  
  #sanity check
  stopifnot('orig.ident' %in% colnames(cmeta))
  stopifnot(all(numerical_cols %in% colnames(cmeta)))
  stopifnot(all(categorical_cols %in% colnames(cmeta)))
  lapply(num_pws_cols,function(num_pws_col){
    stopifnot(all(num_pws_col %in% colnames(cmeta)))
  })
  
  
  cell_n.df = cmeta %>%
    dplyr::group_by(orig.ident) %>%
    dplyr::summarise(n_cnt = n())
  
  names(numerical_cols) = numerical_cols
  stat_numeric.lst = imap(numerical_cols,function(numerical_col,dmy) {
    cmeta %>%
      dplyr::group_by(orig.ident) %>%
      dplyr::summarise(meta_col=numerical_col, 
                       mean_v=mean(get(numerical_col)),
                       sd_v=sd(get(numerical_col)))
  })
  
  stat_categorical.lst=NULL
  if (!invalid(categorical_cols)) {
    names(categorical_cols) = categorical_cols
    stat_categorical.lst=imap(categorical_cols,function(cat_col,dmy) {
      a=cmeta %>%
        dplyr::group_by(orig.ident,get(cat_col)) %>%
        dplyr::summarise(n_cnt=n())
      colnames(a) = c('orig.ident','categorical','n_cnt')
      a
    })
  }
  
  # ---------------------
  message("violin plot (read count) ...")
  # pdf_fn = sprintf("%s_plots.pdf",outd_pref)
  pdf(file=pdf_fn, width=8, height=8)
  imap(numerical_cols,function(num_col,dmy) {
    mean_v=stat_numeric.lst[[num_col]] %>%
      dplyr::summarize(mean_v=mean(mean_v)) %>%
      dplyr::pull(mean_v)
    
    if (mean_v>100) {
      y_scale = "log10(y)"
      p <- ggplot(cmeta, aes(x = orig.ident,
                             y = log10(get(num_col)))) +
        geom_violin(trim=TRUE)
    } else {
      y_scale = 'y'
      p <- ggplot(cmeta, aes(x = orig.ident,
                             y = get(num_col))) +
        geom_violin(trim=TRUE)
    }
    
    p <- p + stat_summary(fun.data=mean_sdl,
                          geom="pointrange", 
                          color="red") +
      labs(title=sprintf("%s [QC:%s]",num_col,qc_comment),
           x="samples",
           y=y_scale)+
      coord_flip()
    
    plot(p)
    NA
  })
  
  # pairwise
  message("scatter plots between read count vs. feature counts ...")
  cmetas = split(cmeta, by="orig.ident")
  orig.idents=names(cmetas)
  
  batches=split_vec(orig.idents,chunk_len=4)
  names(batches)=1:length(batches)
  
  lapply(num_pws_cols,function(num_pws_col) {
    num_col.x=num_pws_col[[1]]
    num_col.y=num_pws_col[[2]]
    imap(batches,function(snames,pgn){
      plst=imap(cmetas[snames],function(cmeta2,sname) {
        p=ggplot(cmeta2, aes(x=get(num_col.x),y=get(num_col.y))) +
          geom_point()+
          coord_trans(x="log10", y="log10")+
          labs(x=num_col.x,y=num_col.y,title=sname)
        p
      })
      p=wrap_plots(plst,ncol=2,nrow=2) + plot_annotation(paste0(num_pws_col,collapse='_'))
      plot(p)
      NA
    })
  })
  
  if (!invalid(categorical_cols)) {
    imap(stat_categorical.lst,function(stat_dt,cat_col) {
      p=ggplot(stat_dt, aes(fill=categorical, y=n_cnt, x=orig.ident)) + 
        geom_bar(position="fill", stat="identity") +
        coord_flip() +
        labs(title=sprintf("%s [QC:%s]",cat_col,qc_comment))
      
      plot(p)
      NA
    })
  }
  dev.off()
  
  s.lst = list()
  s.lst[['cell_n']] = cell_n.df
  s.lst =c(s.lst, stat_numeric.lst)
  if (!invalid(stat_categorical.lst)) {
    s.lst =c(s.lst, stat_categorical.lst)
  }
  s.lst
}

transfer_masked_slot.ams_to_cmeta <- function(seu,assay_name,min_ams=0.15,qco=0.95) {
  DefaultAssay(seu)=assay_name
  data.m = GetAssayData(seu,slot="data")
  N = dim(data.m)[2]
  data.masked = lapply(1:nrow(data.m),function(r2) {
    v=rep('else',N)
    ams_co=quantile(data.m[r2,],qco,na.rm=TRUE)
    if (min_ams>ams_co){ams_co = min_ams}
    v[data.m[r2,]>=ams_co]='higher'
    v
  })
  data.masked = do.call('rbind',data.masked) %>%
    t() %>%
    as.data.table()
  colnames(data.masked) = rownames(data.m)
  data.masked[,cbc:=colnames(data.m)]
  
  cmeta=merge(get_cmeta_from_seu(seu),data.masked,by="cbc")
  cbc.0=colnames(seu)
  seu@meta.data = data.frame(cmeta[match(cbc.0,cbc),-'cbc'],row.names = cbc.0)
  seu
}

SpatiallyVariableFeatures_workaround <- function(object, assay="SCT", selection.method = "moransi") {
  #' This is work around function to replace SeuratObject::SpatiallyVariableFeatures function.
  #' return ranked list of Spatially Variable Features
  
  # Check if object is a Seurat object
  if (!inherits(object, "Seurat")) {
    stop("object must be a Seurat object")
  }
  
  # Check if assay is a valid assay
  if (!assay %in% names(object@assays)) {
    stop("assay must be a valid assay")
  }
  
  # Extract meta.features from the specified object and assay
  data <- object@assays[[assay]]@meta.features
  
  # Select columns starting with the provided col_prefix
  moransi_cols <- grep(paste0("^", selection.method), colnames(data), value = TRUE)
  
  # Filter rows where "moransi.spatially.variable" is TRUE
  filtered_data <- data[data[[paste0(selection.method, ".spatially.variable")]], moransi_cols]
  
  # Sort filtered data by "moransi.spatially.variable.rank" column in ascending order
  sorted_data <- filtered_data[order(filtered_data[[paste0(selection.method, ".spatially.variable.rank")]]), ]
  
  # Return row names of the sorted data frame
  rownames(sorted_data)
}
cellchat_with_tryCatch <- function(seu, scale.factors=NULL, assay="SCT",trim_r=0.1,debug2=0) {
  cellchat <- tryCatch(
    {
      cellchat=cellchat_on_visium_seu(seu, scale.factors=scale.factors, assay=assay,trim_r=trim_r,debug2=debug2)
    },
    error=function(cond){
      message('cellcat_was_not_successful')
      return(NA)
    },
    finally={
      message('Done')
    }
  )
  return(cellchat)
}


cellchat_on_visium_seu <- function(seu, 
                                   scale.factors=NULL,
                                   assay="SCT",
                                   impute_on_ppi=1,
                                   trim_r=0.1,
                                   nboot=100,
                                   population_size=FALSE,
                                   debug2=0) {
  
  if (debug2==1){browser()}
  
  cmeta=get_cmeta_from_seu(seu)
  cnames=colnames(cmeta)
  # >=======
  
  # Idents(seu) = args$cellgroup
  # message(sprintf("Idents(seu)=[%s]",Idents(seu)))
  
  data.input = GetAssayData(seu, 
                            slot = "data", 
                            assay = assay) # normalized data matrix
  
  meta = data.frame(labels = Idents(seu), 
                    row.names = names(Idents(seu))) # manually create a dataframe consisting of the cell labels
  
  # unique(meta$labels) # check the cell labels
  
  # load spatial imaging information
  # Spatial locations of spots from full (NOT high/low) resolution images are required
  if (invalid(scale.factors)) {
    cellchat <- createCellChat(object = data.input, 
                               meta = meta, 
                               group.by = "labels")
  } else {
    spatial.locs = GetTissueCoordinates(seu, 
                                        scale = NULL, 
                                        cols = c("imagerow", "imagecol")) 
    
    # Scale factors and spot diameters of the full resolution images 
    
    
    scale.factors = list(spot.diameter = 65, 
                         spot = scale.factors$spot_diameter_fullres, # these two information are required
                         fiducial = scale.factors$fiducial_diameter_fullres, 
                         hires = scale.factors$tissue_hires_scalef, 
                         lowres = scale.factors$tissue_lowres_scalef # these three information are not required
    )
    cellchat <- createCellChat(object = data.input,
                               meta = meta, 
                               group.by = "labels",
                               datatype = "spatial", 
                               coordinates = spatial.locs,
                               scale.factors = scale.factors)
  }
  
  cellchat@DB <- CellChatDB.human
  # <-------
  
  
  # >=======
  # subset the expression data of signaling genes for saving computation cost
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  
  #> parallelly::supportsMulticore
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
  raw_use=TRUE
  if (impute_on_ppi==1) {
    cellchat <- projectData(cellchat, PPI.human)
    raw_use=FALSE
  }
  
  if (invalid(scale.factors)) {
    cellchat <- computeCommunProb(cellchat,
                                  type = "triMean", 
                                  raw.use = raw_use)
  } else {
    cellchat <- computeCommunProb(cellchat,
                                  type = "truncatedMean", 
                                  trim = trim_r, 
                                  raw.use = raw_use,
                                  population.size = population_size,
                                  distance.use = TRUE,
                                  interaction.length = 200,
                                  nboot=nboot,
                                  scale.distance = 0.01)
  }
  
  
  #> truncatedMean is used for calculating the average gene expression per cell group. 
  #> [1] ">>> Run CellChat on spatial imaging data using distances as constraints <<< [2022-11-12 07:49:23]"
  #> The suggested minimum value of scaled distances is in [1,2], and the calculated value here is  1.30553 
  #> [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2022-11-12 08:10:42]"
  
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 5)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  # Compute the network centrality scores
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
  # <-------
  cellchat
}

# use cellranger_to_qcseurats
load_read_matrix <- function(input_type, input_loc, sname, excl_feat_pat="^DEPRECATED",debug2=0) {
  if (debug2==1){browser()}
  seu <- NULL
  only_ge <- 0
  is_matrices <- 1
  if (input_type=="bc.matrix.dir") {
    seu <- Read10X(data.dir = input_loc)
    if (invalid(names(seu))) {only_ge <- 1}
  } else if (input_type=="bc.matrix.dir_spatial") {
    stop('TODO: to support bc.matrix.dir and spatial')
  } else if (input_type=="h5_spatial") {
    seu = Load10X_Spatial(data.dir = dirname(input_loc),
                          filename = basename(input_loc))
    
    if (!invalid(excl_feat_pat)) {
      seu <- subset(x=seu,
                    features=rownames(seu)[!grepl(excl_feat_pat,rownames(seu))])
    }
    is_matrices=0
  } else {
    seu <- fread_as_matrix(input_loc)
    if (invalid(names(seu))) {only_ge <- 1}
  }
  
  
  
  if (is_matrices==1){
    if (only_ge==1) {
      seu <- CreateSeuratObject(counts = seu,
                                project = sname,
                                min.cells = 0,
                                min.features = 0)
      if (!invalid(excl_feat_pat)) {
        seu <- subset(x=seu,
                      features=rownames(seu)[!grepl(excl_feat_pat,rownames(seu))])
      }
    } else {
      seu <- CreateSeuratObject(counts = seu$`Gene Expression`,
                                project = sname,
                                min.cells = 0,
                                min.features = 0)
      # browser()
      if ('Antibody Capture' %in% names(seu)) {
        message("Storing Protein Feature Barcodes ...")
        seu[['Protein']] = CreateAssayObject(counts = seu$`Antibody Capture`)
      }
      
      if ('Pathoscope2' %in% names(seu)) {
        message("Storing Pathoscope2 ...")
        seu[['Pathoscope2']] = CreateAssayObject(counts = seu$`Pathoscope2`)
      }
      seu <- NormalizeData(seu)
      
      assay_names = names(seu)
      if ("Custom" %in% assay_names) {
        assay_names[assay_names=="Custom"] <- "HTO"
        names(seu) <- assay_names
        seu[['HTO']] <- CreateAssayObject(counts = seu$HTO)
        
        seu <- NormalizeData(seu, assay = "HTO", normalization.method = "CLR")
        # seu <- HTODemux_with_tryCatch(seu,sample,pos_quantile=args$pos_quantile,debug2=0)
        # seu <- MULTIseqDemux(seu,assay="HTO",autoThresh=TRUE,maxiter = 10)
      }
    }
  }
  seu
}

#NOTE: this func is supported in only Seurat v5+
downsample_by_sketch <- function(seu, dns_cell_n=50000) {
  if (dns_cell_n < dim(seu)[2]) {
    def_assay.bkp = DefaultAssay(seu)
    
    seu=NormalizeData(seu) %>%
      FindVariableFeatures() %>%
      SketchData(ncells = dns_cell_n,
                 method = "LeverageScore",
                 sketched.assay = "sketch")
    
    DefaultAssay(seu) = def_assay.bkp
    seu <- subset(x=seu, cells=colnames(seu@assays$sketch$counts))
    seu@assays$sketch <- list()
    seu@assays$sketch <- NULL
  }
  seu
}


get_meta_colname<-function(seu,col_pat="nCount_") {
  cmeta_cols=colnames(seu@meta.data)
  column_name=cmeta_cols[grepl(sprintf('^%s',col_pat),cmeta_cols)]
  stopifnot(!invalid(column_name))
  column_name
}

subset_seu <- function(seu,cgroup,cgroup_col="seurat_clusters") {
  cmeta=get_cmeta_from_seu(seu)
  stopifnot(cgroup_col %in% colnames(cmeta))
  j= cmeta[[cgroup_col]] %in% cgroup
  subset(seu, cells = cmeta$cbc[j])
}

ikarus_normal_cells <- function(seu,top_r=0.25) {
  DefaultAssay(seu) = 'ikarus.ref'
  # FeaturePlot_scCustom(iseu, features = rownames(iseu))
  
  top_r = 0.25
  normal_topk = data.table(normal_margin=GetAssayData(seu,slot = "scale.data")['Normal',] 
                           - GetAssayData(seu,slot = "scale.data")['Tumor',],
                           cbc=colnames(seu)) %>%
    dplyr::arrange(-normal_margin) %>%
    dplyr::slice(1: round(dim(seu)[2]*top_r))
  
  # normal_topk
  
  cmeta=get_cmeta_from_seu(seu) %>%
    merge(normal_topk,by="cbc",all.x=T) %>%
    dplyr::mutate(normal_margin=if_else(is.na(normal_margin),0,normal_margin))
  
  seu[['normal_margin']]=cmeta[match(colnames(seu),cbc),normal_margin]
  # FeaturePlot(iseu, features = 'normal_margin')
  
  norm.cbc = normal_topk %>%
    dplyr::pull(cbc)
}

knn_recover_unknown_cid <- function(seu, cellid="pmid33398161", unk_cellid="unknown", knn=5, lowdim='umap', debug2=0) {
  
  if (debug2==1){browser()}
  
  cmeta = get_cmeta_from_seu(seu)
  setnames(cmeta,cellid,'cellid2')
  
  umap.dist =seu[[lowdim]]@cell.embeddings %>%
    dist() %>%
    as.matrix() %>%
    reshape2::melt(varnames=c('cbc','cbc2')) %>%
    dplyr::filter(cbc!=cbc2) %>%
    merge(cmeta[cellid2==`unk_cellid`,list(cbc,cellid2)],by="cbc") %>%
    dplyr::group_by(cbc) %>%
    dplyr::top_n(n=knn, wt = -value)
  
  umap.dist2 = merge(umap.dist,
                     cmeta[cellid2!=`unk_cellid`,list(cbc,cellid2)],
                     by.x="cbc2",by.y="cbc") %>%
    dplyr::group_by(cbc,cellid2.y) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::arrange(cbc, cellid2.y) %>%
    as.data.table()
  
  umap.dist2_n = umap.dist2[,.(total_n=sum(n)),by="cbc"] %>%
    as.data.table() %>%
    merge(umap.dist2,by="cbc")
  
  umap.dist2_n[,frac:= (n/total_n)]
  
  umap.dist2_n_max = umap.dist2_n %>%
    dplyr::group_by(cbc) %>%
    dplyr::top_n(n=1,wt=frac) %>%
    dplyr::filter(frac>=0.5) %>%
    as.data.table()
  
  cmeta = umap.dist2_n_max[,.N,by="cbc"] %>%
    dplyr::filter(N==1) %>%
    merge(umap.dist2_n_max,by="cbc") %>%
    dplyr::select(cbc, cellid2.y) %>%
    merge(cmeta,by="cbc",all.y=T)
  
  cmeta[!is.na(cellid2.y) & cellid2==`unk_cellid`,cellid2:=cellid2.y]
  
  seu[[cellid]] = cmeta[match(colnames(seu),cbc),cellid2]
  seu
}


annotate_cellgroup_w_compv <- function(iseu,
                                       cellgroup="pmid33398161",
                                       comp_var="Response",
                                       comp_map=list(ctrl='Responder',
                                                     expr='NonResponder'),
                                       comp_short=list(ctrl='R',
                                                       expr='NR')) {
  
  # cellgroup="pmid33398161"
  # comp_var="Response"
  # 
  # comp_map=list(ctrl='Responder',
  #               expr='NonResponder')
  # 
  # comp_short=list(ctrl='R',
  #                 expr='NR')
  
  comp_map.dt = imap(comp_map,function(orig1,comp1) {
    data.table(comp_group=comp1,
               orig=orig1,
               brief=comp_short[[comp1]])
  }) %>% rbindlist()
  
  cmeta = get_cmeta_from_seu(iseu)
  
  stopifnot(cellgroup %in% colnames(cmeta))
  stopifnot(any(comp_map %in% names(table(cmeta[[comp_var]]))))
  
  cmeta <- cmeta %>%
    merge(comp_map.dt[,list(orig,brief)],by.x=`comp_var`,by.y='orig',all.x=T) %>%
    dplyr::mutate(cellgroup.comp = sprintf('%s.%s',get(cellgroup),brief))
  
  cr_sorted = sort(unique(cmeta$cellgroup.comp))
  
  cmeta[,cellgroup.comp:=factor(cellgroup.comp, levels = cr_sorted)]
  
  iseu[['cellgroup.comp']] = cmeta[match(colnames(iseu),cbc),cellgroup.comp]
  
  # DefaultAssay(iseu) = 'SCT'
  Idents(iseu) = 'cellgroup.comp'
  iseu
}



copykat_core <- function(seu, 
                         assay_to_use='SCT',
                         ikarus_gset.lst=NULL, 
                         out_pref='./copykat_out',
                         top_r=0.25,
                         ncpu=8,
                         debug2=0) {
  
  if (debug2==1){browser()}
  default_assay.bkp = DefaultAssay(seu)
  stopifnot(assay_to_use %in% names(seu@assays))
  stopifnot('seurat_clusters' %in% colnames(seu@meta.data))
  # create_dir_if_not_exist(outd)
  # stopifnot(dir.exists(outd))
  
  DefaultAssay(seu) = assay_to_use
  
  if (invalid(ikarus_gset.lst)) {
    ikarus_gset.lst = load_ikarus_gset_list()
  }
  
  
  seu.2=AddModuleScore_Assay(seu = seu,
                             max_iter = max_iter,
                             gset_lst = ikarus_gset.lst,
                             cell_group_by = "seurat_clusters",
                             new_assay_name = "ikarus.ref")
  copykat_pred.v = NULL
  if (invalid(seu.2)) {
    copykat_pred.v = rep('not.defined',dim(seu)[2])
  } else {
    
    normal.cbc=ikarus_normal_cells(seu.2,top_r = top_r)
    
    copykat.test <- copykat(rawmat=GetAssayData(seu.2,
                                                assay = assay_to_use,
                                                slot='count'),
                            id.type="S",
                            ngene.chr=5,
                            win.size=25,
                            KS.cut=0.1, 
                            sam.name=out_pref,
                            distance="euclidean", 
                            norm.cell.names=normal.cbc,
                            output.seg="FLASE", 
                            plot.genes="TRUE", 
                            genome="hg20",
                            n.cores=ncpu)
    
    ck.pred=as.data.table(copykat.test$prediction)
    copykat_pred.v = ck.pred[match(colnames(seu.2),cell.names),copykat.pred]
  }
  copykat_pred.v
}

copykat_w_tryCatch <- function(seu,
                               top_r=0.25,
                               out_pref='./copykat_sname',
                               assay_to_use='SCT',
                               ikarus_gset.lst=NULL,
                               ncpu=4,
                               debug2=0) {
  if (debug2==1){browser()}
  copykat_pred.v <- tryCatch(
    {
      copykat_pred.v = copykat_core(seu,
                                    top_r=top_r,
                                    out_pref=out_pref,
                                    assay_to_use=assay_to_use,
                                    ncpu=ncpu,
                                    ikarus_gset.lst = ikarus_gset.lst)
    },
    error=function(cond) {
      message(sprintf("[%s] not successful",out_pref))
      return(NULL)
    },
    finally={
      message(sprintf('Done[%s]',out_pref))
    }
  )
  
  if (invalid(copykat_pred.v)) {
    seu[["copykat"]] = rep('not.defined',dim(seu)[2])
  } else {
    seu[["copykat"]] = copykat_pred.v
  }
  seu
}


