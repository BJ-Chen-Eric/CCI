library(RcppParallel)
# library(FastGGM)
library(matrixcalc)
library(MASS)
library(data.table)
library(dplyr)
library(Matrix)
library(scTenifoldNet)
library(Seurat)
library(cutpointr)
library(stringr)
library(ggplot2)
library(igraph)
library(qgraph)
library(ggrepel)
library(RcisTarget)
library(ggpubr)
# library(SeuratData)
library(patchwork)
library(progress)

seu_preprocess <- function(seurat, specy='human', remove_cell=T, dim=1:20, res=0.5, scale=T
                           ,pre_process=T , decomposition=F, max='', mt_p=10)  {
  t1 <- Sys.time()
  library(Seurat)
  library(dplyr)
  if(isTRUE(pre_process))  {
    if(identical(specy, 'mouse'))  {mito <- '^mt-'}
    if(identical(specy, 'human'))  {mito <- '^MT-'}
    seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = mito)
    # A <- as.data.frame(seurat@assays$RNA@counts@Dimnames[[1]])
    # grep(pattern = '^IGHA', x = seurat@assays$RNA@counts@Dimnames[[1]], value = T)
    print(VlnPlot(seurat, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3))
    # plot1 <- FeatureScatter(seurat, feature1 = "nCount_Spatial", feature2 = "percent.mt")
    # plot2 <- FeatureScatter(seurat, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")
    # plot1 + plot2
    if(isTRUE(remove_cell))  {
      if(identical('', max))  {
        max <- min(boxplot.stats(seurat$nCount_Spatial)$out)
      }
      seurat <- subset(seurat, subset = nFeature_Spatial> 200 & nCount_Spatial < max & percent.mt < 10)
    }
    dim(seurat)
    seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)
    # seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
    
    # Identify the 10 most highly variable genes
    # top10 <- head(VariableFeatures(seurat), 10)
    
    # plot variable features with and without labels
    # plot1 <- VariableFeaturePlot(seurat)
    # plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
    # plot1 + plot2
    if(isTRUE(scale))  {
      seurat <- ScaleData(seurat, features = rownames(seurat))
    }
    t2 <- Sys.time()
    print(paste('preprocess time is ', t2-t1, sep = ''))
  }
  
  A <- as.data.frame(rownames(seurat))
  row <- ''
  for(i in c('^AC[0-9]+','^AF[0-9]+','^AL[0-9]+','^AP[0-9]+','^AC[0-9]+'
             ,'^BX[0-9]+','^CU[0-9]+','^FP[0-9]+','^LINC[0-9]+'
             ,'^MT-','^Z[0-9]+'))  { #, '^IG[HKL]'
    row <- c(row, grep(rownames(seurat), pattern = i))
  }
  row <- row[-c(1)]
  seurat <- seurat[-c(as.numeric(row)),]
  
  
  # VizDimLoadings(seurat, dims = 1:2, reduction = "pca")
  # DimPlot(seurat, reduction = "pca")
  # 
  # seurat <- JackStraw(seurat, num.replicate = 100)
  # seurat <- ScoreJackStraw(seurat, dims = 1:20)
  # ElbowPlot(seurat)
  if(decomposition==T)  {
    t1 <- Sys.time()
    # seurat <- seurat[VariableFeatures(all), ]
    seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 3000)
    seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))
    t2 <- Sys.time()
    print(paste('PCA time is ', t2-t1, sep=' '))
    
    t1 <- Sys.time()
    seurat <- FindNeighbors(seurat, dims = dim)
    seurat <- FindClusters(seurat, resolution = res)
    
    seurat <- RunTSNE(seurat, dims = dim, check_duplicates = FALSE)
    seurat <- RunUMAP(seurat, dims = dim)
    t2 <- Sys.time()
    print(paste('de-dimension time is ', t2-t1, sep = ''))
  }
  
  return(seurat)
}

prepare_data <- function(seu, HVG=T, ori_ident, cut_off=4, nHVG=1000
                         , detectborder=T)  {
  A <- as.data.frame(rownames(seu))
  row <- ''
  for(i in c('^AC[0-9]+','^AF[0-9]+','^AL[0-9]+','^AP[0-9]+','^AC[0-9]+'
             ,'^BX[0-9]+','^CU[0-9]+','^FP[0-9]+','^LINC[0-9]+'
             ,'^MT-','^Z[0-9]+', '^IG[HKL]'))  {
    row <- c(row, grep(rownames(seu), pattern = i))
  }
  row <- row[-c(1)]
  seu <- seu[-c(as.numeric(row)),]
  # seu <- seu[-c(grep(rownames(seu), pattern = '^AF[0-9]+')),]
  # seu <- seu[-c(grep(rownames(seu), pattern = '^AL[0-9]+')),]
  # seu <- seu[-c(grep(rownames(seu), pattern = '^AP[0-9]+')),]
  # seu <- seu[-c(grep(rownames(seu), pattern = '^AC[0-9]+')),]
  # seu <- seu[-c(grep(rownames(seu), pattern = '^BX[0-9]+')),]
  # seu <- seu[-c(grep(rownames(seu), pattern = '^CU[0-9]+')),]
  # seu <- seu[-c(grep(rownames(seu), pattern = '^FP[0-9]+')),]
  # seu <- seu[-c(grep(rownames(seu), pattern = '^LINC[0-9]+')),]
  # seu <- seu[-c(grep(rownames(seu), pattern = '^MT-')),]
  # seu <- seu[-c(grep(rownames(seu), pattern = '^Z[0-9]+')),]
  if(HVG==T)  {
    seu <- FindVariableFeatures(seu, nfeature=nHVG, assay='Spatial')
    seu <- seu[rownames(seu) %in% VariableFeatures(seu)[1:nHVG],]
  }
  
  if(isTRUE(detectborder))  {
    meta <- seu@images[[1]]@coordinates[,2:3]
    A <- as.data.frame(seu[[ori_ident]])
    meta$seu <- A[rownames(A) %in% rownames(meta),1]
    meta$xct_ident <- ''
    meta$index <- seq(1, nrow(meta))
    
    for(i in unique(meta[,3]))  {
      cell_i <- rownames(meta[meta[,3] %in% i,])
      A <- as.data.frame(as.character(unique(meta[,3])))
      nest <- A[!(A[,1] %in% i),]
      for(j in nest)  {
        cell_j <- rownames(meta[meta[,3] %in% j,])
        meta_subset_i <- meta[rownames(meta) %in% cell_i, ]
        meta_subset_j <- meta[rownames(meta) %in% cell_j, ]
        
        row_dis <- abs(matrix(meta_subset_i$row, nrow = nrow(meta_subset_j), ncol = nrow(meta_subset_i), byrow = T)-
                         matrix(meta_subset_j$row, nrow = nrow(meta_subset_j), ncol = nrow(meta_subset_i), byrow = F))
        
        col_dis <- abs(matrix(meta_subset_i$col, nrow = nrow(meta_subset_j), ncol = nrow(meta_subset_i), byrow = T)-
                         matrix(meta_subset_j$col, nrow = nrow(meta_subset_j), ncol = nrow(meta_subset_i), byrow = F))
        
        dis <- t(row_dis+col_dis)
        adj <- ifelse(test = dis<cut_off, yes = 1, no = 0)
        rownames(adj) <- rownames(meta_subset_i)
        colnames(adj) <- rownames(meta_subset_j)
        
        dim(adj)
        A_sum <- as.data.frame(colSums(adj)) 
        A_sum[!(A_sum[,1] %in% 0), 'class'] <- paste(j, '_close_to_', i, sep = '')
        A_sum <- A_sum[!(is.na(A_sum$class)), ]
        
        B_sum <- as.data.frame(rowSums(adj))
        B_sum[!(B_sum[,1] %in% 0), 'class'] <- paste(i, '_close_to_', j, sep = '')
        B_sum <- B_sum[!(is.na(B_sum$class)), ]
        
        meta[rownames(meta) %in% rownames(A_sum), 'xct_ident'] <- A_sum$class
        meta[rownames(meta) %in% rownames(B_sum), 'xct_ident'] <- B_sum$class
      }
    }
    
    result <- as.data.frame(matrix(data = '', nrow = 0, 0))
    for(i in unique(meta$seu))  {
      cell_i <- meta[meta$seu %in% i, ]
      cell_i[cell_i$xct_ident %in% '', 4] <- paste(i, '_inside', sep = '')
      result <- rbind(result, cell_i)
    }
    
    result <- result[order(result$index, decreasing = F), ]
    # table(rownames(result) == colnames(seu))
    A <- as.data.frame(seu$seurat_clusters)
    A$xct <- 'discard'
    A$index <- seq(1, nrow(A))
    A$sp <- rownames(A) %in% rownames(result)
    row <- A$index[A$sp %in% TRUE]
    A[row,'xct'] <- result[,'xct_ident']
    seu$xct_ident <- A$xct
  }
  return(seu)
}

# LR production
# spatialXct <- function(seu, sender, receiver, cluster_label
#                        , weight=T, adjacent= T, specy='human')  {
#   db <- read.csv('~/Analysis/spatial/data_baseLR.csv')
#   if(identical(specy, 'mouse'))  {
#     db[,2:3] <- apply(db[,2:3], MARGIN = 2, tolower)
#     db[,2:3] <- apply(db[,2:3], MARGIN = 2, stringr::str_to_title)
#   }
#   
#   t1 <- Sys.time()
#   assay <- DefaultAssay(seu)
#   sender_cell <- subset(seu, subset= xct_ident== sender)
#   sender_cell <- sender_cell@assays[[assay]]@counts
#   
#   receiver_cell <- subset(seu, subset= xct_ident== receiver)
#   receiver_cell <- receiver_cell@assays[[assay]]@counts
#   
#   
#   meta <- as.data.frame(seu@images$slice1@coordinates$row)
#   colnames(meta) <- 'row'
#   meta$col <- seu@images$slice1@coordinates$col
#   A <- seu[[cluster_label]]
#   meta$xct_ident <- A[,1]
#   meta$index <- seq(1, nrow(meta))
#   rownames(meta)<- colnames(seu)
#   
#   
#   cell_i <- rownames(meta[meta$xct_ident %in% sender,])
#   cell_j <- rownames(meta[meta$xct_ident %in% receiver,])
#   meta_subset_i <- meta[rownames(meta) %in% cell_i, ]
#   meta_subset_j <- meta[rownames(meta) %in% cell_j, ]
#   
#   row_dis <- abs(matrix(meta_subset_i$row, nrow = nrow(meta_subset_j), ncol = nrow(meta_subset_i), byrow = T)-
#                    matrix(meta_subset_j$row, nrow = nrow(meta_subset_j), ncol = nrow(meta_subset_i), byrow = F))
#   
#   col_dis <- abs(matrix(meta_subset_i$col, nrow = nrow(meta_subset_j), ncol = nrow(meta_subset_i), byrow = T)-
#                    matrix(meta_subset_j$col, nrow = nrow(meta_subset_j), ncol = nrow(meta_subset_i), byrow = F))
#   
#   dis <- t(row_dis+col_dis)
#   
#   if(weight==F)  {
#     adj <- ifelse(test = dis<=4, yes = 1, no = 0)
#   }
#   
#   if(weight==T)  {
#     adj <- ifelse(test = dis<=30, yes = dis, no = 0)
#     adj <- ifelse(test = adj==0, yes = 0, no = 1/adj)
#   }
#   
#   dim(adj)
#   adj <- as.vector(t(adj))
#   pair <- length(adj)- table(adj)['0']
#   
#   B <- as.data.frame(rowSums(receiver_cell) !=0) %>%
#     tibble::rownames_to_column()
#   B <- B[B[,2] %in% T, ]
#   receiver_cell <- receiver_cell[rownames(receiver_cell) %in% B$rowname, ]
#   B <- intersect(rownames(receiver_cell), db$receptor)
#   R_count <- db[db$receptor %in% B,c(1,3)]
#   
#   
#   B <- as.data.frame(rowSums(sender_cell) !=0) %>% tibble::rownames_to_column()
#   B <- B[B[,2] %in% T, ]
#   sender_cell <- sender_cell[rownames(sender_cell) %in% B$rowname, ]
#   B <- intersect(rownames(sender_cell), db$ligand)
#   L_count <- db[db$ligand %in% B,1:2]
#   
#   B <- intersect(L_count$X, R_count$X)
#   L_count <- L_count[L_count$X %in% B, ]
#   R_count <- R_count[R_count$X %in% B, ]
#   
#   library(Matrix)
#   # result <- as(matrix(data = 0, nrow = nrow(L_count), ncol = ncol(receiver_cell)*ncol(sender_cell)), "dgCMatrix")
#   
#   result <- matrix(data = 0, nrow = (nrow(L_count)+1), ncol = ncol(receiver_cell)*ncol(sender_cell))
#   rownames(result) <- c(paste(L_count$ligand, R_count$receptor, sep = '_'), 'cell_pairs')
#   
#   for(j in seq(1,ncol(sender_cell)))  {
#     L <- as.data.frame(sender_cell[rownames(sender_cell) %in% unique(L_count$ligand),j])
#     for(i in unique(L_count$ligand))  {
#       L_count[L_count$ligand %in% i,j+2] <- L[rownames(L) %in% i, ]
#     }
#   }
#   
#   for(k in seq(1,ncol(receiver_cell)))  {
#     R <- as.data.frame(receiver_cell[rownames(receiver_cell) %in% unique(R_count$receptor),k])
#     for(i in unique(R_count$receptor))  {
#       R_count[R_count$receptor %in% i,k+2] <- R[rownames(R) %in% i, ]
#     }
#   }
#   
#   x=0
#   for(i in seq(1, ncol(sender_cell)))  {
#     for(j in seq(1, ncol(receiver_cell)))  {
#       x=x+1
#       result[1:(nrow(result)-1),x] <- L_count[,(i+2)]*R_count[,(j+2)]
#     }
#   }
#   
#   
#   if(adjacent==T)  {
#     for(i in seq(1, (nrow(result)-1)))  {
#       result[i,] <- result[i,]*adj
#     }
#   }
#   t2 <- Sys.time()
#   print(t2-t1)
#   
#   # filter 0 LR production 
#   A <- as.data.frame(rowSums(result))
#   A$index <- seq(1, nrow(A))
#   result <- list(LR=as(result[A$index[!(A[,1] %in% 0)],], 'dgCMatrix')
#                  , pair=pair)
#   return(result)
#   gc()
# }

merge_run.test <- function(close_in, far_in, test='ks')  {
  close_pair <- close_in[[2]]; far_pair <- far_in[[2]]
  close <- close_in[[1]]; far <- far_in[[1]]
  if(is.na(close_pair))  {close_pair <- ncol(close)}
  if(is.na(far_pair))  {far_pair <- ncol(far)}
  inter <- intersect(rownames(far), rownames(close))
  close <- close[rownames(close) %in% inter, ]
  far <- far[rownames(far) %in% inter, ]
  
  out <- as.data.frame(rownames(far))
  if(test=='ks')  {
    for(i in seq(nrow(out)))  {
      stest <- ks.test(close[i,], far[i,])
      out[i,paste(test, 'p', sep = '.')] <- stest$p.value
      out[i,'D'] <- stest$statistic
    }
  }
  
  if(test=='t')  {
    for(i in seq(nrow(far)))  {
      stest <- t.test(close[i,], far[i,])
      out[i,paste(test, 'p', sep = '.')] <- stest$p.value
      out[i,'D'] <- stest$statistic
    }
  }
  
  if(test=='wc')  {
    for(i in seq(nrow(far)))  {
      stest <- wilcox.test(close[i,], far[i,])
      out[i,paste(test, 'p', sep = '.')] <- stest$p.value
      out[i,'D'] <- stest$statistic
    }
  }
  
  # out$close <- apply(X = close_filter, MARGIN = 1, FUN = mean)
  # out$far <- apply(X = far_filter, MARGIN = 1, FUN = mean)
  
  out$close <- rowSums(close)/close_pair
  out$far <- rowSums(far)/far_pair
  
  out <- out[order(out[,2], decreasing = F), ]
  out[,'adj.p'] <- p.adjust(out[,2], method='BH')
  
  out[,'direction'] <- ifelse(out[,'close']>out[,'far'], yes = 'close', no = 'far')
  
  A <- out[out$adj.p < 0.05, ]
  print(nrow(A)/nrow(out))
  print(table(out$direction[out$adj.p < 0.05]))
  A <- as.list(A %>% group_split(direction))
  for(i in seq_along(A))  {
    B <- A[[i]]
    B <- B[order(B$D, decreasing = T), ]
    A[[i]] <- as.data.frame(B)
    names(A)[i] <- B[1,"direction"]
  }
  return(A)
}



# Xct_output <- function(seu, cluster1, cluster2, specy='human'
#                        , specific=F, border1, border2, inner1, inner2)  {
#   
#   if(isTRUE(specific))  {
#     close <- spatialXct(seu = seu, sender = border1, cluster_label = 'xct_ident'
#                         , receiver = border2, weight = T, adjacent = T, specy=specy)
#     far <- spatialXct(seu = seu, sender = inner1, cluster_label = 'xct_ident'
#                       , receiver = inner2, weight = T, adjacent = T, specy=specy)
#     
#   } else{
#     border1 <- paste(cluster1, '_close_to_', cluster2, sep = '')
#     border2 <- paste(cluster2, '_close_to_', cluster1, sep = '')
#     
#     inner1 <- paste(cluster1, '_inside', sep = '')
#     inner2 <- paste(cluster2, '_inside', sep = '')
#     
#     
#     close <- spatialXct(seu = seu, sender = border1, cluster_label = 'xct_ident'
#                         , receiver = border2, weight = T, adjacent = T, specy=specy)
#     far <- spatialXct(seu = seu, sender = inner1, cluster_label = 'xct_ident'
#                       , receiver = inner2, weight = T, adjacent = T, specy=specy)
#   }
#   
#   
#   
#   result <- merge_run.test(close_in = close, far_in = far)
#   c12 <- result[['close']]
#   
#   if(isTRUE(specific))  {
#     close <- spatialXct(seu = seu, sender = border2, cluster_label = 'xct_ident'
#                         , receiver = border1, weight = T, adjacent = T, specy=specy)
#     far <- spatialXct(seu = seu, sender = inner2, cluster_label = 'xct_ident'
#                       , receiver = inner1, weight = T, adjacent = T, specy=specy)
#     
#   } else{
#     border1 <- paste(cluster1, '_close_to_', cluster2, sep = '')
#     border2 <- paste(cluster2, '_close_to_', cluster1, sep = '')
#     
#     inner1 <- paste(cluster1, '_inside', sep = '')
#     inner2 <- paste(cluster2, '_inside', sep = '')
#     
#     
#     close <- spatialXct(seu = seu, sender = border2, cluster_label = 'xct_ident'
#                         , receiver = border1, weight = T, adjacent = T, specy=specy)
#     far <- spatialXct(seu = seu, sender = inner2, cluster_label = 'xct_ident'
#                       , receiver = inner1, weight = T, adjacent = T, specy=specy)
#   }
#   
#   result <- merge_run.test(close_in = close, far_in = far)
#   
#   c21 <- result[['close']]
#   
#   c12_out <- setdiff(c12[,1], c21[,1]) # get the L_R in one of them
#   c12_out <- c12[c12[,1] %in% c12_out, ]
#   
#   c21_out <- setdiff(c21[,1], c12[,1])
#   c21_out <- c21[c21[,1] %in% c21_out, ]
#   
#   
#   out_list <- list(c12_out, c21_out)
#   names(out_list) <- c(border2, border1)
#   for(i in seq_along(out_list))  {
#     A <- out_list[[i]]
#     B <- do.call(rbind, strsplit(A[,1], split = '_'))
#     colnames(B) <- c('Ligand', 'Receptor')
#     out_list[[i]] <- cbind(B, A[,2:ncol(A)])
#   }
#   return(out_list)
# }


gene_subset <- function(seu, drug_path, drug_b=T, tf_path, tf_b=T
                       , nHVG=100, L_R, specy='human')  {
  library(Seurat)
  library(dplyr)
  if(isTRUE(drug_b))  {
    drug <- fread(drug_path)
    drug <- drug[grep(drug$interaction_claim_source, pattern = 'Cancer'), ]
    drug <- drug$gene_name %>% unique()
  }else{drug <- ''}
  
  if(isTRUE(tf_b))  {
    tf <- fread(tf_path)
    tf <- tf[grep(tf$tissue, pattern = 'lung'), ]
    tf <- unique(tf$TF)
  }else{tf <- ''}
  
  gene <- c(drug, tf, VariableFeatures(seu)[1:nHVG]) %>% unique() #, unique(L_R)
  gene <- gene[!(gene %in% '')]
  if(identical(specy, 'mouse'))  {
    gene <- tolower(gene)
    gene <- stringr::str_to_title(gene)
  }
  return(gene)
  # seu <- prepare_data(seu = seu, HVG = F, ori_ident = 'seurat_clusters'
  #              , cut_off = 4 ,nHVG = nHVG)
  # 
  # 
  # gex <- subset(seu, subset = xct_ident==xct_cluster)
  # gex <- gex[gene, ]
  # gex <- gex[['Spatial']]@data
}


seu_preprocess <- function(seurat, specy='human', dim=1:20, res=0.5, scale=T
                           ,pre_process=T , decomposition=F, max='', mt_p=10)  {
  t1 <- Sys.time()
  library(Seurat)
  library(dplyr)
  if(isTRUE(pre_process))  {
    mito <- '^MT-'
    if(specy=='mouse')  {mito <- '^mt-'}
    seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern =mito)
    # A <- as.data.frame(seurat@assays$RNA@counts@Dimnames[[1]])
    # grep(pattern = '^IGHA', x = seurat@assays$RNA@counts@Dimnames[[1]], value = T)
    print(VlnPlot(seurat, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3))
    # plot1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
    # plot2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    # plot1 + plot2
    # if(identical('', max))  {
    #   max <- min(boxplot.stats(seurat$nCount_Spatial)$out)
    # }
    # seurat <- subset(seurat, subset = nFeature_Spatial > 200 & nCount_Spatial < max & percent.mt < 10)
    # dim(seurat)
    seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)
    seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
    
    # Identify the 10 most highly variable genes
    # top10 <- head(VariableFeatures(seurat), 10)
    
    # plot variable features with and without labels
    # plot1 <- VariableFeaturePlot(seurat)
    # plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
    # plot1 + plot2
    if(isTRUE(scale))  {
      seurat <- ScaleData(seurat, features = rownames(seurat))
    }
    t2 <- Sys.time()
    print(paste('preprocess time is ', t2-t1, sep = ''))
  }
  
  A <- as.data.frame(rownames(seurat))
  row <- ''
  for(i in c('^AC[0-9]+','^AF[0-9]+','^AL[0-9]+','^AP[0-9]+','^AC[0-9]+'
             ,'^BX[0-9]+','^CU[0-9]+','^FP[0-9]+','^LINC[0-9]+'
             ,'^MT-','^Z[0-9]+'))  { #, '^IG[HKL]'
    row <- c(row, grep(rownames(seurat), pattern = i))
  }
  row <- row[-c(1)]
  seurat <- seurat[-c(as.numeric(row)),]
  
  
  # VizDimLoadings(seurat, dims = 1:2, reduction = "pca")
  # DimPlot(seurat, reduction = "pca")
  # 
  # seurat <- JackStraw(seurat, num.replicate = 100)
  # seurat <- ScoreJackStraw(seurat, dims = 1:20)
  # ElbowPlot(seurat)
  if(decomposition==T)  {
    t1 <- Sys.time()
    # seurat <- seurat[VariableFeatures(all), ]
    seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 3000)
    seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))
    t2 <- Sys.time()
    print('PCA time is ')
    t2-t1
    
    t1 <- Sys.time()
    seurat <- FindNeighbors(seurat, dims = dim)
    seurat <- FindClusters(seurat, resolution = res)
    
    seurat <- RunTSNE(seurat, dims = dim, check_duplicates = FALSE)
    seurat <- RunUMAP(seurat, dims = dim)
    t2 <- Sys.time()
    print(paste('de-dimension time is ', t2-t1, sep = ''))
  }
  
  return(seurat)
}


as.tensor <- function(x,drop=FALSE){
  stopifnot(is.array(x)||is.vector(x))
  tnsr <- list()
  if (is.vector(x)){
    modes <- c(length(x))
    num_modes <- 1L
    tnsr$modes <- modes
    tnsr$num_modes <- num_modes
    tnsr$data <- x
    #new("Tensor", num_modes, modes, data = x)
  }
  else {
    modes <- dim(x)
    num_modes <- length(modes)
    dim1s <- which(modes==1)
    if (drop && (length(dim1s)>0)){
      modes <- modes[-dim1s]
      num_modes <- num_modes-length(dim1s)
      tnsr$modes <- modes
      tnsr$num_modes <- num_modes
      tnsr$data <- array(x,dim=modes)
      #new("list",num_modes,modes,data=array(x,dim=modes))
    } else {
      tnsr$modes <- modes
      tnsr$num_modes <- num_modes
      tnsr$data <- x
    }
  }
  return(tnsr)
}



cpDecomposition <- function(tnsr, num_components=NULL,max_iter=25, tol=1e-5){
  kronecker_list <- function(L){
    isvecORmat <- function(x){is.matrix(x) || is.vector(x)}
    stopifnot(all(unlist(lapply(L,isvecORmat))))
    retmat <- L[[1]]
    for(i in 2:length(L)){
      retmat <- kronecker(retmat,L[[i]])
    }
    retmat
  }
  
  
  fnorm <- function(tnsr){
    arr<-tnsr$data
    sqrt(sum(arr*arr))
  }
  
  rs_unfold <- function(tnsr,m=NULL){
    if(is.null(m)) stop("mode m must be specified")
    num_modes <- tnsr$num_modes
    rs <- m
    cs <- (1:num_modes)[-m]
    unfold(tnsr,row_idx=rs,col_idx=cs)
  }
  
  unfold <- function(tnsr,row_idx=NULL,col_idx=NULL){
    #checks
    rs <- row_idx
    cs <- col_idx
    if(is.null(rs)||is.null(cs)) stop("row and column indices must be specified")
    num_modes <- tnsr$num_modes
    if (length(rs) + length(cs) != num_modes) stop("incorrect number of indices")
    if(any(rs<1) || any(rs>num_modes) || any(cs < 1) || any(cs>num_modes)) stop("illegal indices specified")
    perm <- c(rs,cs)
    if (any(sort(perm,decreasing=TRUE) != num_modes:1)) stop("missing and/or repeated indices")
    modes <- tnsr$modes
    mat <- tnsr$data
    new_modes <- c(prod(modes[rs]),prod(modes[cs]))
    #rearranges into a matrix
    mat <- aperm(mat,perm)
    dim(mat) <- new_modes
    as.tensor(mat)
  }
  
  hadamard_list <- function(L){
    isvecORmat <- function(x){is.matrix(x) || is.vector(x)}
    stopifnot(all(unlist(lapply(L,isvecORmat))))
    retmat <- L[[1]]
    for (i in 2:length(L)){
      retmat <- retmat*L[[i]]
    }
    retmat
  }
  
  khatri_rao_list <- function(L,reverse=FALSE){
    stopifnot(all(unlist(lapply(L,is.matrix))))
    ncols <- unlist(lapply(L,ncol))
    stopifnot(length(unique(ncols))==1)
    ncols <- ncols[1]
    nrows <- unlist(lapply(L,nrow))
    retmat <- matrix(0,nrow=prod(nrows),ncol=ncols)
    if (reverse) L <- rev(L)
    for(j in 1:ncols){
      Lj <- lapply(L,function(x) x[,j])
      retmat[,j] <- kronecker_list(Lj)
    }
    retmat
  }
  
  superdiagonal_tensor <- function(num_modes,len,elements=1L){
    modes <- rep(len,num_modes)
    arr <- array(0, dim = modes)
    if(length(elements)==1) elements <- rep(elements,len)
    for (i in 1:len){
      txt <- paste("arr[",paste(rep("i", num_modes),collapse=","),"] <- ", elements[i],sep="")
      eval(parse(text=txt))
    }
    as.tensor(arr)
  }
  
  ttl<-function(tnsr,list_mat,ms=NULL){
    if(is.null(ms)||!is.vector(ms)) stop ("m modes must be specified as a vector")
    if(length(ms)!=length(list_mat)) stop("m modes length does not match list_mat length")
    num_mats <- length(list_mat)
    if(length(unique(ms))!=num_mats) warning("consider pre-multiplying matrices for the same m for speed")
    mat_nrows <- vector("list", num_mats)
    mat_ncols <- vector("list", num_mats)
    for(i in 1:num_mats){
      mat <- list_mat[[i]]
      m <- ms[i]
      mat_dims <- dim(mat)
      modes_in <- tnsr$modes
      stopifnot(modes_in[m]==mat_dims[2])
      modes_out <- modes_in
      modes_out[m] <- mat_dims[1]
      tnsr_m <- rs_unfold(tnsr,m=m)$data
      retarr_m <- mat%*%tnsr_m
      tnsr <- rs_fold(retarr_m,m=m,modes=modes_out)
    }	
    tnsr	
  }
  
  rs_fold <- function(mat,m=NULL,modes=NULL){
    if(is.null(m)) stop("mode m must be specified")
    if(is.null(modes)) stop("Tensor modes must be specified")
    num_modes <- length(modes)
    rs <- m
    cs <- (1:num_modes)[-m]
    fold(mat,row_idx=rs,col_idx=cs,modes=modes)
  }
  
  fold <- function(mat, row_idx = NULL, col_idx = NULL, modes=NULL){
    #checks
    rs <- row_idx
    cs <- col_idx
    if(is.null(rs)||is.null(cs)) stop("row space and col space indices must be specified")
    if(is.null(modes)) stop("Tensor modes must be specified")
    if(!is(mat,"list")){
      if(!is.matrix(mat))  stop("mat must be of class 'matrix'")
    }else{
      stopifnot(mat$num_modes==2)
      mat <- mat$data			
    }
    num_modes <- length(modes)
    stopifnot(num_modes==length(rs)+length(cs))
    mat_modes <- dim(mat)
    if((mat_modes[1]!=prod(modes[rs])) || (mat_modes[2]!=prod(modes[cs]))) stop("matrix nrow/ncol does not match Tensor modes")
    #rearranges into array
    iperm <- match(1:num_modes,c(rs,cs))
    as.tensor(aperm(array(mat,dim=c(modes[rs],modes[cs])),iperm))
  }
  
  
  if(is.null(num_components)) stop("num_components must be specified")
  stopifnot(is(tnsr,"list"))
  #if (.is_zero_tensor(tnsr)) stop("Zero tensor detected")
  
  #initialization via truncated hosvd
  num_modes <- tnsr$num_modes
  modes <- tnsr$modes
  U_list <- vector("list",num_modes)
  unfolded_mat <- vector("list",num_modes)
  tnsr_norm <- fnorm(tnsr)
  for(m in 1:num_modes){
    unfolded_mat[[m]] <- rs_unfold(tnsr,m=m)$data
    U_list[[m]] <- matrix(rnorm(modes[m]*num_components), nrow=modes[m], ncol=num_components)
  }
  est <- tnsr
  curr_iter <- 1
  converged <- FALSE
  #set up convergence check
  fnorm_resid <- rep(0, max_iter)
  CHECK_CONV <- function(est){
    curr_resid <- fnorm(as.tensor(est$data - tnsr$data))
    fnorm_resid[curr_iter] <<- curr_resid
    if (curr_iter==1) return(FALSE)
    if (abs(curr_resid-fnorm_resid[curr_iter-1])/tnsr_norm < tol) return(TRUE)
    else{ return(FALSE)}
  }	
  #progress bar
  pb <- txtProgressBar(min=0,max=max_iter,style=3)
  #main loop (until convergence or max_iter)
  norm_vec <- function(vec){
    norm(as.matrix(vec))
  }
  while((curr_iter < max_iter) && (!converged)){
    setTxtProgressBar(pb,curr_iter)
    for(m in 1:num_modes){
      V <- hadamard_list(lapply(U_list[-m],function(x) {t(x)%*%x}))
      V_inv <- solve(V)			
      tmp <- unfolded_mat[[m]]%*%khatri_rao_list(U_list[-m],reverse=TRUE)%*%V_inv
      lambdas <- apply(tmp,2,norm_vec)
      U_list[[m]] <- sweep(tmp,2,lambdas,"/")	
      Z <- superdiagonal_tensor(num_modes=num_modes,len=num_components,elements=lambdas)
      est <- ttl(Z,U_list,ms=1:num_modes)
    }
    #checks convergence
    if(CHECK_CONV(est)){
      converged <- TRUE
      setTxtProgressBar(pb,max_iter)
    }else{
      curr_iter <- curr_iter + 1
    }
  }
  if(!converged){setTxtProgressBar(pb,max_iter)}
  close(pb)
  #end of main loop
  #put together return list, and returns
  fnorm_resid <- fnorm_resid[fnorm_resid!=0]
  norm_percent<- (1-(tail(fnorm_resid,1)/tnsr_norm))*100
  invisible(list(lambdas=lambdas, U=U_list, conv=converged, est=est, norm_percent=norm_percent, fnorm_resid = tail(fnorm_resid,1),all_resids=fnorm_resid))
}


makeL_R_list <- function(seu, specy='human', times=100, sender, receiver, ncell=100, max_dis=30
                         , weight = T, adjacent = T, cluster_label = 'xct_ident', sampling=T)  {
  db <- read.csv('~/Analysis/spatial/data_baseLR.csv')
  if(identical(specy, 'mouse'))  {
    db[,2:3] <- apply(db[,2:3], MARGIN = 2, tolower)
    db[,2:3] <- apply(db[,2:3], MARGIN = 2, stringr::str_to_title)
  }
  assay <- DefaultAssay(seu)
  # gex matrix
  sender_cell_raw <- subset(seu, subset= xct_ident== sender)
  sender_cell_raw <- sender_cell_raw@assays[[assay]]@counts
  ncell_s <- ncol(sender_cell_raw)
  # set.seed(713)
  sender_ind_list <- list()
  if(isTRUE(sampling))  {
    if(ncol(sender_cell_raw) < ncell)  {ncell <- ncol(sender_cell_raw)}
    for(i in 1:times)  {
      sender_ind_list[[i]] <- sample(colnames(sender_cell_raw), size = ncell, replace = F)
    }
  }else{sender_ind_list[[1]] <- colnames(sender_cell_raw)}
  
  
  
  receiver_cell_raw <- subset(seu, subset= xct_ident== receiver)
  receiver_cell_raw <- receiver_cell_raw@assays[[assay]]@counts
  ncell_r <- ncol(receiver_cell_raw)
  # set.seed(713)
  receiver_ind_list <- list()
  if(isTRUE(sampling))  {
    if(ncol(receiver_cell_raw) < ncell)  {ncell <- ncol(receiver_cell_raw)}
    for(i in 1:times)  {
      receiver_ind_list[[i]] <- sample(colnames(receiver_cell_raw), size = ncell, replace = F)
    }
  }else{receiver_ind_list[[1]] <- colnames(receiver_cell_raw)}
  
  
  if(isFALSE(sampling))  {times=1}
  out_list <- list()
  for(time in 1:times)  {
    sender_cell <- sender_cell_raw[, colnames(sender_cell_raw) %in% sender_ind_list[[time]]]
    receiver_cell <- receiver_cell_raw[, colnames(receiver_cell_raw) %in% receiver_ind_list[[time]]]
    
    meta <- as.data.frame(seu@images$slice1@coordinates$row)
    colnames(meta) <- 'row'
    meta$col <- seu@images$slice1@coordinates$col
    A <- seu[[cluster_label]]
    meta$xct_ident <- A[,1]
    meta$index <- seq(1, nrow(meta))
    rownames(meta)<- colnames(seu)
    
    
    # cell_i <- rownames(meta[meta$xct_ident %in% sender,])
    # cell_j <- rownames(meta[meta$xct_ident %in% receiver,])
    meta_subset_i <- meta[rownames(meta) %in% colnames(sender_cell), ]
    meta_subset_j <- meta[rownames(meta) %in% colnames(receiver_cell), ]
    
    row_dis <- abs(matrix(meta_subset_i$row, nrow = nrow(meta_subset_j), ncol = nrow(meta_subset_i), byrow = T)-
                     matrix(meta_subset_j$row, nrow = nrow(meta_subset_j), ncol = nrow(meta_subset_i), byrow = F))
    
    col_dis <- abs(matrix(meta_subset_i$col, nrow = nrow(meta_subset_j), ncol = nrow(meta_subset_i), byrow = T)-
                     matrix(meta_subset_j$col, nrow = nrow(meta_subset_j), ncol = nrow(meta_subset_i), byrow = F))
    
    # dis=1st norm
    dis <- t(row_dis+col_dis)
    rownames(dis) <- colnames(sender_cell)
    colnames(dis) <- colnames(receiver_cell)
    
    # col_ind <- sample(1:ncol(dis), size = 200, replace = F) %>% as.data.frame()
    # col_ind <- col_ind[order(col_ind$., decreasing = F), ]
    # row_ind <- sample(1:nrow(dis), size = 200, replace = F) %>% as.data.frame()
    # row_ind <- row_ind[order(row_ind$., decreasing = F), ]
    # 
    # dis <- dis[row_ind, col_ind]
    
    if(weight==F)  {
      adj <- ifelse(test = dis<=4, yes = 1, no = 0)
    }
    
    # dis to adj (inverse distance)
    if(weight==T)  {
      adj <- ifelse(test = dis<=max_dis, yes = dis, no = 0)
      adj <- ifelse(test = adj==0, yes = 0, no = 1/adj)
    }
    
    
    dim(adj)
    # adj <- as.vector(t(adj))
    # pair <- length(adj)- table(adj)['0']
    
    # filter 0 expression genes and grep from L_R meta_data
    B <- ifelse(as.matrix(receiver_cell) ==0, yes = 0, no = 1)
    B <- as.data.frame(rowSums(B) >= ncell*0.05) %>% tibble::rownames_to_column()
    
    B <- B[B[,2] %in% T, ]
    receiver_cell <- receiver_cell[rownames(receiver_cell) %in% B$rowname, ]
    B <- intersect(rownames(receiver_cell), db$receptor)
    R_count <- db[db$receptor %in% B,c(1,3)]
    
    B <- ifelse(as.matrix(sender_cell) ==0, yes = 0, no = 1)
    B <- as.data.frame(rowSums(B) >= 10) %>% tibble::rownames_to_column()
    B <- B[B[,2] %in% T, ]
    sender_cell <- sender_cell[rownames(sender_cell) %in% B$rowname, ]
    B <- intersect(rownames(sender_cell), db$ligand)
    L_count <- db[db$ligand %in% B,1:2]
    
    
    B <- intersect(L_count$X, R_count$X)
    L_count <- L_count[L_count$X %in% B, ]
    R_count <- R_count[R_count$X %in% B, ]
    
    library(Matrix)
    # result <- as(matrix(data = 0, nrow = nrow(L_count), ncol = ncol(receiver_cell)*ncol(sender_cell)), "dgCMatrix")
    
    # result <- matrix(data = 0, nrow = (nrow(L_count)+1), ncol = ncol(receiver_cell)*ncol(sender_cell))
    # rownames(result) <- c(paste(L_count$ligand, R_count$receptor, sep = '_'), 'cell_pairs')
    
    for(j in seq(1,ncol(sender_cell)))  {
      L <- as.data.frame(sender_cell[rownames(sender_cell) %in% unique(L_count$ligand),j])
      for(i in unique(L_count$ligand))  {
        L_count[L_count$ligand %in% i,j+2] <- L[rownames(L) %in% i, ]
      }
    }
    
    for(j in seq(1,ncol(receiver_cell)))  {
      R <- as.data.frame(receiver_cell[rownames(receiver_cell) %in% unique(R_count$receptor),j])
      for(i in unique(R_count$receptor))  {
        R_count[R_count$receptor %in% i,j+2] <- R[rownames(R) %in% i, ]
      }
    }
    
    
    if(adjacent==T)  {
      subLR_list <- list()
      for(k in seq_len(nrow(L_count)))  {
        A <- L_count[k, c(3:ncol(L_count))] %>% as.matrix()
        B <- R_count[k, c(3:ncol(R_count))] %>% as.matrix()
        
        subLR_list[[k]] <- t(A) %*% B * adj %>% as('dgCMatrix')
      }
    }
    if(adjacent==F)  {
      subLR_list <- list()
      for(k in seq_len(nrow(L_count)))  {
        A <- L_count[k, c(3:ncol(L_count))] %>% as.matrix()
        B <- R_count[k, c(3:ncol(R_count))] %>% as.matrix()
        
        subLR_list[[k]] <- t(A) %*% B %>% as('dgCMatrix')
      }
    }
    names(subLR_list) <- paste(L_count$ligand, R_count$receptor, sep = "_")
    out_list[[time]] <- subLR_list
  }
  return(out_list)
}



# Tensor_ks_test
# hbc_xct <- spatialXct(seu = hbc, specy = 'human', cluster_1 = 'Immune', cluster_2 = 'Invasive', Tensor_dec = T)
spatialXct <- function(seu, specy='human', cluster_1, cluster_2, Tensor_dec=T, nHVGs=2000, sampling)  {
  makeL_R_list_meta <- function(seu, specy='human', sender_i, sender, receiver_i, receiver
                                , Tensor_dec=T, HVGs=nHVGs)  {
    db <- read.csv('~/Analysis/spatial/data_baseLR.csv')
    db$L_R <- paste(db$ligand, db$receptor)
    db <- distinct(db, L_R, .keep_all = T)[, 1:3]
    if(identical(specy, 'mouse'))  {
      db[,2:3] <- apply(db[,2:3], MARGIN = 2, tolower)
      db[,2:3] <- apply(db[,2:3], MARGIN = 2, stringr::str_to_title)
    }
    assay <- DefaultAssay(seu)
    
    cluster_s <- subset(seu, subset= xct_ident %in% c(sender_i, sender))
    cluster_s <- FindVariableFeatures(cluster_s, nfeature=HVGs)
    cluster_s <- VariableFeatures(cluster_s)
    # gex matrix
    sender_cell_i <- subset(seu, subset= xct_ident== sender_i)
    sender_cell_i <- sender_cell_i@assays[[assay]]@data[cluster_s,]
    sender_cell <- subset(seu, subset= xct_ident== sender)
    sender_cell <- sender_cell@assays[[assay]]@data[cluster_s,]
    ncell <- min(c(ncol(sender_cell), ncol(sender_cell_i)))


    test <- sample(1:ncol(sender_cell_i), size = ncell, replace = F) %>% as.data.frame()
    sender_cell_i <- sender_cell_i[, test[order(test$., decreasing = F), ]]
    test <- sample(1:ncol(sender_cell), size = ncell, replace = F) %>% as.data.frame()
    sender_cell <- sender_cell[, test[order(test$., decreasing = F), ]]
    # if(ncell != 'all')  {
    #   test <- sample(1:ncol(sender_cell_i), size = ncell, replace = F) %>% as.data.frame()
    #   sender_cell_i <- sender_cell_i[, test[order(test$., decreasing = F), ]]
    #   test <- sample(1:ncol(sender_cell), size = ncell, replace = F) %>% as.data.frame()
    #   sender_cell <- sender_cell[, test[order(test$., decreasing = F), ]]
    # }

    cluster_r <- subset(seu, subset= xct_ident %in% c(receiver_i, receiver))
    cluster_r <- FindVariableFeatures(cluster_r, nfeature=HVGs)
    cluster_r <- VariableFeatures(cluster_r)
    receiver_cell_i <- subset(seu, subset= xct_ident== receiver_i)
    receiver_cell_i <- receiver_cell_i@assays[[assay]]@data[cluster_r, ]
    receiver_cell <- subset(seu, subset= xct_ident== receiver)
    receiver_cell <- receiver_cell@assays[[assay]]@data[cluster_r, ]
    ncell <- c(ncol(receiver_cell), ncol(receiver_cell_i)) %>% min()


    test <- sample(1:ncol(receiver_cell_i), size = ncell, replace = F) %>% as.data.frame()
    receiver_cell_i <- receiver_cell_i[, test[order(test$., decreasing = F), ]]
    test <- sample(1:ncol(receiver_cell), size = ncell, replace = F) %>% as.data.frame()
    receiver_cell <- receiver_cell[, test[order(test$., decreasing = F), ]]
    # if(ncell != 'all')  {
    #   test <- sample(1:ncol(receiver_cell_i), size = ncell, replace = F) %>% as.data.frame()
    #   receiver_cell_i <- receiver_cell_i[, test[order(test$., decreasing = F), ]]
    #   test <- sample(1:ncol(receiver_cell), size = ncell, replace = F) %>% as.data.frame()
    #   receiver_cell <- receiver_cell[, test[order(test$., decreasing = F), ]]
    # }


    sender_cell <- cbind(sender_cell, sender_cell_i)
    receiver_cell <- cbind(receiver_cell, receiver_cell_i)
    cat(paste('sender_cell: ', ncol(sender_cell), ' receiver_cell: ', ncol(receiver_cell), '\n', sep = ''))

    meta <- as.data.frame(seu@images$slice1@coordinates$row)
    colnames(meta) <- 'row'
    meta$col <- seu@images$slice1@coordinates$col
    A <- seu[['xct_ident']]
    meta$xct_ident <- A[,1]
    meta$index <- seq(1, nrow(meta))
    rownames(meta) <- colnames(seu)


    meta_subset_i <- meta[rownames(meta) %in% colnames(sender_cell), ]
    meta_subset_j <- meta[rownames(meta) %in% colnames(receiver_cell), ]

    row_dis <- abs(matrix(meta_subset_i$row, nrow = nrow(meta_subset_j), ncol = nrow(meta_subset_i), byrow = T)-
                     matrix(meta_subset_j$row, nrow = nrow(meta_subset_j), ncol = nrow(meta_subset_i), byrow = F))

    col_dis <- abs(matrix(meta_subset_i$col, nrow = nrow(meta_subset_j), ncol = nrow(meta_subset_i), byrow = T)-
                     matrix(meta_subset_j$col, nrow = nrow(meta_subset_j), ncol = nrow(meta_subset_i), byrow = F))

    # dis=1st norm
    dis <- t(row_dis+col_dis) %>% as.data.frame()
    rownames(dis) <- colnames(sender_cell)
    colnames(dis) <- colnames(receiver_cell)

    adj <- ifelse(test = dis<=4, yes = 1, no = 0)
    inn <- ifelse(test = dis>4, yes = 1, no = 0)


    dim(adj)
    # adj <- as.vector(t(adj))
    # pair <- length(adj)- table(adj)['0']

    ## filter 0 expression genes and grep from L_R meta_data
    # B <- ifelse(as.matrix(receiver_cell) ==0, yes = 0, no = 1)
    # B <- as.data.frame(rowSums(B) >= ncol(receiver_cell)*0.05) %>% tibble::rownames_to_column() # express in at least 5 % cells

    # B <- B[B[,2] %in% T, ]
    # receiver_cell <- receiver_cell[rownames(receiver_cell) %in% B$rowname, ]
    B <- intersect(rownames(receiver_cell), db$receptor)
    R_count <- db[db$receptor %in% B,c(1,3)]

    # B <- ifelse(as.matrix(sender_cell) ==0, yes = 0, no = 1)
    # # B <- as.data.frame(rowSums(B) >= ncol(sender_cell)*0.05) %>% tibble::rownames_to_column()
    # B <- B[B[,2] %in% T, ]
    # sender_cell <- sender_cell[rownames(sender_cell) %in% B$rowname, ]
    B <- intersect(rownames(sender_cell), db$ligand)
    L_count <- db[db$ligand %in% B,1:2]


    B <- intersect(L_count$X, R_count$X)
    cat(paste('Number of LR pairs filtered: ', length(B), ', HVGs: ', HVGs, '\n', sep = ''))
    L_count <- L_count[L_count$X %in% B, ]
    R_count <- R_count[R_count$X %in% B, ]

    library(Matrix)
    # result <- as(matrix(data = 0, nrow = nrow(L_count), ncol = ncol(receiver_cell)*ncol(sender_cell)), "dgCMatrix")

    # result <- matrix(data = 0, nrow = (nrow(L_count)+1), ncol = ncol(receiver_cell)*ncol(sender_cell))
    # rownames(result) <- c(paste(L_count$ligand, R_count$receptor, sep = '_'), 'cell_pairs')

    for(j in seq(1,ncol(sender_cell)))  {
      L <- as.data.frame(sender_cell[rownames(sender_cell) %in% unique(L_count$ligand),j])
      for(i in unique(L_count$ligand))  {
        L_count[L_count$ligand %in% i,j+2] <- L[rownames(L) %in% i, ]
      }
    }

    for(j in seq(1,ncol(receiver_cell)))  {
      R <- as.data.frame(receiver_cell[rownames(receiver_cell) %in% unique(R_count$receptor),j])
      for(i in unique(R_count$receptor))  {
        R_count[R_count$receptor %in% i,j+2] <- R[rownames(R) %in% i, ]
      }
    }


    subLR_list <- list()
    for(k in seq_len(nrow(L_count)))  {
      A <- L_count[k, c(3:ncol(L_count))] %>% as.matrix()
      B <- R_count[k, c(3:ncol(R_count))] %>% as.matrix()

      subLR_list[[k]] <- t(A) %*% B %>% as('dgCMatrix')
      # subLR_list_n[[k]] <- LR_prod[which(adj==1)]
      # subLR_list_i[[k]] <- LR_prod[which(adj==0)]
    }
    names(subLR_list) <- paste(L_count$ligand, R_count$receptor, sep = "_")

    if(isTRUE(Tensor_dec))  {
      nL <- nrow(subLR_list[[1]])
      nR <- ncol(subLR_list[[1]])
      tensor <- array(data = 0, dim = c(nL,nR,1, length(subLR_list)))
      for(i in seq_along(subLR_list)){
        tensor[,,,i] <- subLR_list[[i]] %>% as.matrix()
      }
      tensor <- as.tensor(tensor)
      gc()
      tensor <- cpDecomposition(tnsr = tensor, num_components = 5, max_iter = 500, tol = 1e-5) #, max_iter = 1e1, tol = 1e-5
      gc()
      tensor_data <- lapply(seq(dim(tensor$est$data)[4]), function(x) tensor$est$data[ , , ,x])
      names(tensor_data) <- names(subLR_list)
      subLR_list <- tensor_data
    }

    subLR_list_i <- list()
    subLR_list_n <- list()
    for(i in seq_along(subLR_list))  {
      subLR_list_n[[i]] <- subLR_list[[i]][which(adj==1)]
      subLR_list_i[[i]] <- subLR_list[[i]][which(adj==0)]
    }
    names(subLR_list_n) <- names(subLR_list_i) <- names(subLR_list)

    out_list <- list(neigh_propuct=subLR_list_n, inner_propuct=subLR_list_i)
    return(out_list)
  }
  
  sender_i = paste(cluster_1, '_inside', sep = '')
  sender = paste(cluster_1, '_close_to_', cluster_2, sep = '')
  receiver_i = paste(cluster_2, '_inside', sep = '')
  receiver = paste(cluster_2, '_close_to_', cluster_1, sep = '')
  cat(paste('Computing ', sender_i, ' communicating to ', receiver_i, '\n', sep = ''))
  one2two <- makeL_R_list_meta(seu = seu, specy = specy, HVGs = nHVGs
                               , sender_i=sender_i, sender = sender
                               , receiver_i = receiver_i, receiver = receiver, Tensor_dec = Tensor_dec)
  gc()
  one2two_test <- meta_test(meta_LR_list = one2two, sampling=sampling)
  
  cat(paste('Computing ', receiver_i, ' communicating to ', sender_i, '\n', sep = ''))
  two2one <- makeL_R_list_meta(seu = seu, specy = specy, HVGs = nHVGs
                               , sender_i=receiver_i, sender = receiver
                               , receiver_i = sender_i, receiver = sender, Tensor_dec = Tensor_dec)
  gc()
  two2one_test <- meta_test(meta_LR_list = two2one, sampling=sampling)
  out_list <- list(one2two, one2two_test, two2one, two2one_test)
  names(out_list) <- c(paste(cluster_1, '_close_to_', cluster_2, sep = ''), 
                       paste(cluster_1, '_close_to_', cluster_2, '_test', sep = ''), 
                       paste(cluster_2, '_close_to_', cluster_1, sep = ''), 
                       paste(cluster_2, '_close_to_', cluster_1, '_test', sep = ''))
  return(out_list)
}


diff_genes <- function(seu, drug_path, tf_path, nHVG, specy='human'
                       , spatialxct, ori_ident = 'seurat_clusters'
                       , xct_clu='xct_ident', cluster_1, cluster_2 
                       , alpha=1, tf_b=T, drug_b=T)  {
  ggm_comparison <- function(seu, drug_path, tf_path, nHVG, specy='human'
                             , L_R, ori_ident = 'seurat_clusters'
                             , xct_clu='xct_ident', R_border, R_inside 
                             , alpha=1, tf_b, drug_b)  {
    A <- seu[[xct_clu]] %>% as.data.frame() %>% tibble::rownames_to_column('cell')
    seu <-prepare_data(seu = seu, HVG = F, detectborder = F)
    inside <- R_inside
    subcell <- A[A$xct_ident %in% inside, 'cell']
    seu_inside <- subset(seu, cells = subcell)
    seu_inside <- FindVariableFeatures(seu_inside, nfeatures = nHVG)
    gene_inside <- gene_subset(seu = seu_inside, drug_path = drug_path , tf_path = tf_path, specy=specy
                               , nHVG = nHVG, L_R = L_R, drug_b = drug_b, tf_b = tf_b)
    
    
    border <- R_border
    subcell <- A[A$xct_ident %in% border, 'cell']
    seu_border <- subset(seu, cells = subcell)
    gene_border <- gene_subset(seu = seu_border, drug_path = drug_path , tf_path = tf_path, specy=specy
                               , nHVG = nHVG, L_R = L_R, drug_b = drug_b, tf_b = tf_b)
    
    gene <- union(gene_inside, gene_border) # inside and neighbor cell
    
    devide <- sample(1:ncol(seu_inside), size = ncol(seu_inside)/2, replace = F) # sampling inner cell, building inner cell control
    c1 <- seu_inside[gene, -devide]
    c2 <- seu_inside[gene, devide]
    print(paste('Control GRN gene number=', nrow(c1), sep=''))
    # c1 <- c2 <- seu_inside[gene,]
    c1_cell <- ncol(c1)
    c2_cell <- ncol(c2)
    c1 <- c1[['Spatial']]@data %>% as.matrix()
    c2 <- c2[['Spatial']]@data %>% as.matrix()
    output <- scTenifoldNet(X = c1, Y = c2, qc = F,
                            nc_nNet = 10, nc_nCells = 500,
                            td_K = 3, qc_minLibSize = 30, 
                            nc_symmetric = F, td_nDecimal = 10)
    control = output
    
    
    border_cell <- ncol(seu_border)
    seu_border <- seu_border[gene, ]
    seu_inside <- seu_inside[gene, ]
    seu_border <- seu_border[['Spatial']]@data %>% as.matrix()
    seu_inside <- seu_inside[['Spatial']]@data %>% as.matrix()
    print(paste('B_I GRN gene number=', nrow(seu_inside), sep=''))
    output <- scTenifoldNet(X = seu_border, Y = seu_inside, qc = F,
                            nc_nNet = 10, nc_nCells = 500,
                            td_K = 3, qc_minLibSize = 30,
                            nc_symmetric = F, td_nDecimal = 10)
    
    B <- output[["diffRegulation"]]
    # B <- sct3000[["SCT_GRN_B_I"]][["diffRegulation"]]
    B <- B[B$p.adj <0.05, ]
    
    C <- control$diffRegulation
    # C <- sct3000$SCT_GRN_control$diffRegulation
    C <- C[C$p.adj <0.05, ]
    
    
    inner_g <- setdiff(C$gene, B$gene) # control specific degs
    neighbor_g <- setdiff(B$gene, C$gene) # neighbor specific degs
    
    
    results = list(SCT_GRN_B_I = output,
                   SCT_GRN_control = control,
                   genes= gene,
                   inner_g =inner_g, 
                   neighbor_g=neighbor_g
                   # GRN_p = pvalue_SCT_Diff,
                   # GRN_FDR = GRN_FDR,
                   # FDR_genepair=results
    )
    return(results)
  }
  t1 <- Sys.time()
  # out <- spatialxct[[2]][[1]]
  # out <- replace(out, is.na(out), 0)
  # nei_r <- out[out$ks_pos_times >10 & out$fc >50 , ] #& out$fc >50 
  # nei_r <- nei_r %>% slice_max(ks_pos_times, n = 10) 
  # nei_r <- strsplit(nei_r$L_R, split = '_') %>% unlist()
  # nei_r <- sub(nei_r[,1], pattern = '.*_', replacement = '') %>% unique()
  R_border <- paste(cluster_1, '_close_to_', cluster_2, sep = '')
  R_inside <- paste(cluster_1, '_inside', sep = '')
  one2two <- ggm_comparison(seu = seu ,drug_path = drug_path
                            , tf_path = tf_path
                            , nHVG = nHVG, L_R = nei_r
                            , R_border = R_border
                            , R_inside = R_inside
                            , alpha = 1, drug_b=T, tf_b = T, specy = specy)
  gc()
  
  # out <- spatialxct[[4]][[1]]
  # out <- replace(out, is.na(out), 0)
  # nei_r <- out[out$ks_pos_times >10 & out$fc >50 , ] #& out$fc >50 
  # nei_r <- nei_r %>% slice_max(ks_pos_times, n = 10) 
  # nei_r <- sub(nei_r[,1], pattern = '.*_', replacement = '') %>% unique()
  R_border <- paste(cluster_2, '_close_to_', cluster_1, sep = '')
  R_inside <- paste(cluster_2, '_inside', sep = '')
  two2one <- ggm_comparison(seu = seu ,drug_path = drug_path
                            , tf_path = tf_path
                            , nHVG = nHVG, L_R = nei_r
                            , R_border = R_border
                            , R_inside = R_inside
                            , alpha = 1, drug_b=T, tf_b = T, specy = specy)
  out <- list(c1_neighbor_degs=one2two, c2_neighbor_degs=two2one)
  names(out) <- c(paste(cluster_1, '_close_to_', cluster_2, sep = '')
                  , paste(cluster_2, '_close_to_', cluster_1, sep = ''))
  t2 <- Sys.time()
  print(t2-t1)
  return(out)
}


neighbor_tf <- function(neighbor_deg, inner_deg, Rcis_feather_dir, nTfs=10)  {
    tf_prediction <- function(degs, nTfs=nTfs)  {
    geneLists <- list(geneListName=degs)
    # Motif enrichment analysis:
    motifEnrichmentTable_wGenes <- cisTarget(geneSets = geneLists
                                             , motifRankings = motifRankings
                                             , motifAnnot=motifAnnotations)
    
    
    A <- motifEnrichmentTable_wGenes[!(motifEnrichmentTable_wGenes$TF_highConf %in% ''), ]
    if(nTfs>nrow(A))  {
      print('Ask too many Tfs, set to maximun')
      nTfs <- nrow(A)
    }
    A <- A %>% top_n(NES, n = nTfs)
    B <- sub(A$TF_highConf, pattern = ' \\(.*', replacement = '')
    tf <- strsplit(B, split = '; ') %>% unlist() %>% unique()
    return(tf)
  }
  
  data(motifAnnotations_hgnc)
  motifAnnotations <- motifAnnotations
  motifRankings <- importRankings(Rcis_feather_dir) # "~/Analysis/spatial/tf_df/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather"
  
  if(length(neighbor_deg) !=0)  {
    print("performing the TFs prediction of neighbor_deg")
    neighbor_tf <- tryCatch(
      {
        neighbor_tf <- tf_prediction(degs=neighbor_deg, nTfs=nTfs)
      },
      error = function(error_message) {
        message("Fewer than 80% of the genes. No TF will be predicted")
        return('')
      }
    )
  } 
  else {
    print("Zero neighbor_deg, pass TFs identification")
    neighbor_tf <- ''
  }
  if(length(inner_deg) !=0)  {
    print("Performing the TFs prediction of inner_degs")
    inner_tf <- tryCatch(
      {
        inner_tf <- tf_prediction(degs=inner_deg, nTfs=nTfs)
      },
      error = function(error_message) {
        message("Fewer than 80% of the genes. No TF will be predicted")
        return('')
      }
    )
  }
  else {
    print("Zero neighbor_deg, pass TFs identification")
    inner_tf <- ''
  }
  
  neighbor_tf <- setdiff(neighbor_tf, inner_tf)
  
  return(neighbor_tf)
}


write_output <-  function(xct, n_receptor, n_ligand, deg, c1, c2
                          , c1_nei_tf, c2_nei_tf)  {
  
  c1_c2 <- paste(c1, '_close_to_', c2, sep = '')
  c2_c1 <- paste(c2, '_close_to_', c1, sep = '')
  
  out <- xct[[paste(c1_c2, '_test', sep = '')]]#[[1]] # Immune to invasive receptor
  out <- replace(out, is.na(out), 0)
  out$dir <- ifelse(out$fc >=50, yes = 'H', no = 'L')
  out <- out %>% group_split(dir) %>% as.list() %>% lapply(as.data.frame)
  nei_r <- c()
  nei_l <- c()
  for(i in seq_along(out))  {
    A <- out[[i]]
    A <- A[A$ks_pos_times >50 & A$fc >70 | A$fc < 20 , ] #& out$fc >50
    A <- A %>% slice_max(ks_pos_times, n = n_receptor)
    nei_r <- c(nei_r, paste(sub(A[,1], pattern = '.*_', replacement = '') %>% unique(), ' ', 'R_', A$dir, sep = ''))
    nei_l <- c(nei_l, paste(sub(A[,1], pattern = '_.*', replacement = '') %>% unique(), ' ', 'L_', A$dir, sep = ''))
  }
  
  
  # out <- xct[[paste(c2_c1, '_test', sep = '')]][[1]]# Immune to invasive ligand
  # out <- replace(out, is.na(out), 0)
  # out$dir <- ifelse(out$fc >=50, yes = 'H', no = 'L')
  # out <- out %>% group_split(dir) %>% as.list() %>% lapply(as.data.frame)
  # nei_l <- c()
  # for(i in seq_along(out))  {
  #   A <- out[[i]]
  #   A <- A[A$ks_pos_times >50 & A$fc >70 | A$fc < 20 , ] #& out$fc >50
  #   A <- A %>% slice_max(ks_pos_times, n = n_receptor)
  #   nei_l <- c(nei_l, paste(sub(A[,1], pattern = '.*_', replacement = '') %>% unique(), A$dir, sep = '\t'))
  # }
  
  
  dg <- deg[[paste(c1, c2, sep = '_')]]$neighbor_g # Immune to invasive degs
  # dg <- setdiff(deg[[1]]$neighbor_g, deg[[2]]$neighbor_g) # remove the degs in the neighbor cell against it
  
  
  meta_c1 <- data.frame(gene=c(nei_r, c1_nei_tf, nei_l, dg)
                        , type=c(rep('', length(nei_r))
                                 , rep('TFs', length(c1_nei_tf))
                                 , rep('', length(nei_l))
                                 , rep('DEGs', length(dg))))
  meta_c1$clean <- sub(meta_c1$gene, replacement = '', patter=' .*')
  
  
  write.table(meta_c1$clean, file = paste('~/', names(xct[1]), '_r_deg_tf_l.txt', sep = '')
              , quote = F, col.names = F, row.names = F)
  
  write.table(meta_c1[, 1:2], file = paste('~/', names(xct[1]), '_r_deg_tf_l_label.txt', sep = '')
              , quote = F, col.names = F, row.names = F)
  
  
  out <- xct[[paste(c2_c1, '_test', sep = '')]][[1]] # Immune to invasive receptor
  out <- replace(out, is.na(out), 0)
  out$dir <- ifelse(out$fc >=50, yes = 'H', no = 'L')
  out <- out %>% group_split(dir) %>% as.list() %>% lapply(as.data.frame)
  nei_r <- c()
  nei_l <- c()
  for(i in seq_along(out))  {
    A <- out[[i]]
    A <- A[A$ks_pos_times >50 & A$fc >70 | A$fc < 20 , ] #& out$fc >50
    A <- A %>% slice_max(ks_pos_times, n = n_receptor)
    if(nrow(A) >0)  {
      nei_r <- c(nei_r, paste(sub(A[,1], pattern = '.*_', replacement = '') %>% unique(), ' ', 'R_', A$dir, sep = ''))
      nei_l <- c(nei_l, paste(sub(A[,1], pattern = '_.*', replacement = '') %>% unique(), ' ', 'L_', A$dir, sep = ''))
    }
  }
  
  
  dg <- deg[[paste(c2, c1, sep = '_')]]$neighbor_g
  # dg <- setdiff(deg[[2]]$neighbor_g, deg[[1]]$neighbor_g)
  
  
  meta_c2 <- data.frame(gene=c(nei_r, c2_nei_tf, nei_l, dg)
                        , type=c(rep('', length(nei_r))
                                 , rep('TFs', length(c2_nei_tf))
                                 , rep('', length(nei_l))
                                 , rep('DEGs', length(dg))))
  meta_c2$clean <- sub(meta_c2$gene, replacement = '', patter=' .*')
  
  write.table(meta_c2$clean, file = paste('~/', names(xct[3]), '_r_deg_tf_l.txt', sep = '')
              , quote = F, col.names = F, row.names = F)
  
  write.table(meta_c2[, 1:2], file = paste('~/', names(xct[3]), '_r_deg_tf_l_label.txt', sep = '')
              , quote = F, col.names = F, row.names = F)
  out <- list(meta_c1=meta_c1, meta_c2=meta_c2)
  names(out) <- c(names(xct[1]), names(xct[3]))
  return(out)
}


feature_vio_vol <- function(seu, clu1, clu2, clu1_genes_lable, clu2_genes_lable, out_path)  {
  draw_output <- function(c1, c2, genes_lable, out_path)  {
    c1_nei <- paste(c1, '_close_to_', c2, sep = '')
    c1_inn <- paste(c1, '_inside', sep = '')
    A <- subset(seu, subset = xct_ident %in% c1_nei)
    B <- subset(seu, subset = xct_ident %in% c1_inn)
    
    c2 <- scCustomize::Merge_Seurat_List(list(A,B))
    c2 <- seu[, colnames(seu) %in% colnames(c2)]
    p3 <- SpatialDimPlot(c2, group.by = 'xct_ident')+gg_theme
    
    plot_list <- list()
    x=0
    in_gene <- rownames(c2)[rownames(c2) %in% genes_lable[,1]]
    location <- c2$xct_ident %>% as.data.frame() %>% tibble::rownames_to_column('cell') 
    colnames(location)[2] <- 'xct'
    location <- group_split(location, xct) %>% as.list()
    clu <- c(location[[1]][1,2], location[[2]][1,2]) %>% unlist()
    location <- location %>% lapply(function(x){x[,1] %>% unlist()})
    names(location) <- clu
    p_table <- matrix(0, nrow = length(in_gene), ncol = ) %>% as.data.frame()
    p_table$gene <- in_gene
    for(i in in_gene)  {
      x=x+1
      nei <- c2@assays$Spatial@data[i, location[[c1_nei]]]
      inn <- c2@assays$Spatial@data[i, location[[c1_inn]]] # 0_close_to_2
      p <- wilcox.test(nei, inn, alternative = 'g')$p.value
      p_table[x, 'p'] <- p
      p_table[x, 'logp'] <- -log10(p)
      p_table[x, 'fc'] <- (mean(nei)/mean(inn))
      p_table[x, 'logfc'] <- log2(mean(nei)/mean(inn))
      p2 <- SpatialFeaturePlot(c2, features = i)+gg_theme
      
      p4 <- VlnPlot(c2, features = i, group.by = 'xct_ident', )+NoLegend()+
        labs(subtitle = paste(i, ' ', genes_lable[genes_lable[, 1] %in% i, 2], ', wilcoxon test (greater): ', p, sep = ''))+
        gg_theme #theme(axis.text.x = element_text(size = 18, angle=90,hjust=0.95,vjust=0.2))+
      plot_list[[x]] <- ggarrange(plotlist = list(p2, p4), ncol = 2, nrow = 1)
      names(plot_list)[x] <- i
    }
    colnames(genes_lable)[1] <- 'gene'
    p_table <- merge(p_table[, -c(1)], genes_lable , by.x = 'gene')
    # colnames(p_table)[1] <- c1_nei
    
    plot_list[[paste(c1_nei, '_Volcano', sep = '')]] <- 
      p_table %>% distinct(gene, .keep_all = T) %>% as.data.frame() %>% 
      ggplot(aes(x=logfc, y=logp, label=gene, color=type))+ #=V1
      geom_point(size=3)+
      geom_vline(xintercept = c(1,-1), color= '#F6A707', linetype='dashed', size=1.5)+
      geom_hline(yintercept = -log(0.05), color= '#F6A707', linetype='dashed', size=1.5)+
      scale_color_manual(values = c('#5124CD','#24CD50', '#CD5124'))+
      labs(title = paste('Volcano plot of detected genes (receptors, DEGs, TFs)', sep = '')
           , subtitle = paste('FC=Neighor/Inner, P: wilcoxon test')
           , x='log2FC', y='-log(p-value)', color='')+
      geom_text_repel(nudge_x = 0.1, direction = "y", hjust = "left", size=6)+
      guides(color = guide_legend(title = "", override.aes = aes(label = "")))+
      theme_bw()+gg_theme
    
    
    dir.create(file.path(out_path), showWarnings = FALSE)
    out_path <- paste(out_path, c1_nei, '/', sep = '')
    dir.create(file.path(out_path), showWarnings = FALSE)
    for(i in seq_along(plot_list))  {
      png(paste(out_path, names(plot_list)[i], '.png', sep = ''), width = 1600, height =  900)
      print(plot_list[[i]])
      dev.off()
    }
    return(p_table)
  }
  
  p1 <- draw_output(c1 = clu1, c2 = clu2
              , genes_lable = clu1_genes_lable
              , out_path=out_path)
  
  p2 <- draw_output(c1 = clu2, c2 = clu1
              , genes_lable = clu2_genes_lable
              , out_path=out_path)
  return(list(clu1=p1, clu2=p2))
}





# test_list ---------------------------------------------------------------
meta_test <- function(meta_LR_list, sampling)  {
  pb <- progress_bar$new(total = length(meta_LR_list[[1]]))
  if(isTRUE(sampling))  {
    cat(c('Performing test by subsampling'))
    out_df_tensor <- matrix('', nrow = length(meta_LR_list[[1]])) %>% as.data.frame()
    out_df_tensor$L_R <- names(meta_LR_list[[1]])
    out_df_tensor_list <- list()
    
    for(j in seq_along(meta_LR_list[[1]]))  {
      nei <- meta_LR_list[[1]][[j]]
      inn <- meta_LR_list[[2]][[j]]
      
      out_ks <- c()
      out_wc <- c()
      fc_vec <- c()
      for(i in 1:100)  {
        sub_inn <- inn[sample(1:length(inn), size = length(nei), replace = F)]
        fc_vec <- c(fc_vec, sum(nei)/sum(sub_inn))
        out_ks <- c(out_ks, ks.test(nei, sub_inn)$p.value)
        out_wc <- c(out_wc, wilcox.test(nei, sub_inn)$p.value)
      }
      out_df_tensor[j, 'ks_pos_times'] <- table(out_ks < 0.05)['TRUE']
      out_df_tensor[j, 'ks_neg_times'] <- table(out_ks < 0.05)['FALSE']
      out_df_tensor[j, 'wc_pos_times'] <- table(out_wc < 0.05)['TRUE']
      out_df_tensor[j, 'wc_neg_times'] <- table(out_wc < 0.05)['FALSE']
      out_df_tensor[j, 'fc'] <- table(fc_vec > 1)['TRUE']
      
      pb$tick()
    }
    out_df_tensor <- replace(out_df_tensor, is.na(out_df_tensor), 0)
    out_df_tensor_list[['ks']] <- out_df_tensor[, c(2,3,4,7)]
    out_df_tensor_list[['wc']] <- out_df_tensor[, c(2,5,6,7)]
    return(out_df_tensor_list)
  }
  if(isFALSE(sampling))  {
    cat(c('Performing test by whole data'))
    out <- matrix('', nrow = 0, ncol = 0) %>% as.data.frame()
    for(j in seq_along(meta_LR_list[[1]]))  {
      nei <- meta_LR_list[[1]][[j]]
      inn <- meta_LR_list[[2]][[j]]
      
      out[j, 'fc'] <- mean(nei)/mean(inn)
      out[j, 'ks'] <- ks.test(nei, inn)$p.value
      out[j, 'wc'] <- wilcox.test(nei, inn)$p.value
      pb$tick()
    }
    rownames(out) <- names(meta_LR_list[[1]])
    return(out)
  }
  
}

runtest_sample_cell <- function(neighbor_list, inner_list, test, times)  {
  # if(is.null(method)) stop("Specify the method")
  if(identical('', times)|times>length(neighbor_list))  {times <- length(neighbor)}
  loopL_R <- list()
  for(j in 1:times)  {
    sub_neighbor <- neighbor_list[[j]]
    sub_inner <- inner_list[[j]]
    share_LR <- intersect(names(sub_neighbor), names(sub_inner))
    
    sub_neighbor <- sub_neighbor[share_LR]
    sub_inner <- sub_inner[share_LR]
    
    
    out <- as.data.frame(names(sub_neighbor))
    if(identical('wc', test))  {
      for(i in seq_along(sub_neighbor))  {
        stest <- wilcox.test(as.matrix(sub_neighbor[[i]]), as.matrix(sub_inner[[i]]))
        out[i, 2] <- stest$p.value
      }
    }
    if(identical('ks', test))  {
      for(i in seq_along(sub_neighbor))  {
        stest <- ks.test(as.matrix(sub_neighbor[[i]]), as.matrix(sub_inner[[i]]))
        out[i, 2] <- stest$p.value
      }
      
    }
    out$adj.p <- p.adjust(out[,2], method='BH')
    out$fc <- (lapply(sub_neighbor, sum) %>% as.numeric())/(lapply(sub_inner, sum) %>% as.numeric())
    loopL_R[[j]] <- out
  }
  return(loopL_R)
}



runtest_sample_cell_tensor <- function(neighbor_list, inner_list, test, times
                                       , num_components = 5, max_iter = 1e3, tol = 1e-5)  {
  # if(is.null(method)) stop("Specify the method")
  if(identical('', times)|times>length(neighbor_list))  {times <- length(neighbor)}
  loopL_R <- list()
  for(j in 1:times)  {
    sub_neighbor <- neighbor_list[[j]]
    sub_inner <- inner_list[[j]]
    share_LR <- intersect(names(sub_neighbor), names(sub_inner))
    
    sub_neighbor <- sub_neighbor[share_LR]
    sub_inner <- sub_inner[share_LR]
    
    
    ncell <- nrow(sub_neighbor[[1]])
    tensorI <- tensorN <- array(data = 0, dim = c(ncell, ncell, 1, length(sub_neighbor)))
    for(i in seq_len(length(sub_neighbor)))  {
      tensorN[,,,i] <- sub_neighbor[[i]] %>% as.matrix()
      tensorI[,,,i] <- sub_inner[[i]] %>% as.matrix()
    }
    tensorN <- as.tensor(tensorN)
    tensorI <- as.tensor(tensorI)
    tensorN <- cpDecomposition(tnsr = tensorN, num_components = 5, max_iter = 1e3, tol = 1e-5)
    tN <- lapply(seq(dim(tensorN$est$data)[4]), function(x) tensorN$est$data[ , , ,x])
    tensorI <- cpDecomposition(tnsr = tensorI, num_components = 5, max_iter = 1e3, tol = 1e-5)
    tI <- lapply(seq(dim(tensorI$est$data)[4]), function(x) tensorI$est$data[ , , ,x])
    names(tI) <- names(tN) <- names(sub_neighbor)
    
    out <- as.data.frame(names(tI))
    if(identical('wc', test))  {
      for(i in seq_along(sub_neighbor))  {
        stest <- wilcox.test(tN[[i]], tI[[i]])
        out[i, 2] <- stest$p.value
      }
    }
    if(identical('ks', test))  {
      for(i in seq_along(sub_neighbor))  {
        stest <- ks.test(tN[[i]], tI[[i]])
        out[i, 2] <- stest$p.value
      }
    }
    out$adj.p <- p.adjust(out[,2], method='BH')
    out$fc <- (lapply(sub_neighbor, sum) %>% as.numeric())/(lapply(sub_inner, sum) %>% as.numeric())
    loopL_R[[j]] <- out
    
  }
  return(loopL_R)
}


keep_sgi_LR <- function(list, p=0.1)  {
  out_list <- list()
  for(i in seq_along(list))  {
    rep <- list[[i]]
    out_list[[i]] <- rep[rep$adj.p <p & rep$fc<1, 1]
  }
  return(out_list)
}



# test2
draw_hist <- function(data1, data2, l1='g1', l2='g2', title='')  {
  data1 <- keep_sgi_LR(data1) %>% unlist() %>% table() %>% as.data.frame(stringsAsFactors=F) 
  data1$group <- l1
  freq1 <- data1[data1$Freq >90, ] %>% nrow()
  sig1 <- round(freq1/nrow(data1)*100, digits = 1)
  data2 <- keep_sgi_LR(data2) %>% unlist() %>% table() %>% as.data.frame(stringsAsFactors=F)
  data2$group <- l2
  freq2 <- data2[data2$Freq >90, ] %>% nrow()
  sig2 <- round(freq2/nrow(data2)*100, digits = 1)
  
  rbind(data1, data2) %>% 
    ggplot()+
    geom_histogram(aes(x=Freq, fill=group), color="#e9ecef", alpha=0.5, binwidth = 1, position = 'identity')+
    labs(title = title
         , subtitle = paste(l1, ' identified total LR pair=', nrow(data1), ', > 90% iteration=', freq1, ', Significant rate is ', sig1, '%', '\n',
                            l2, ' identified total LR pair=', nrow(data2), ', > 90% iteration=', freq2, ', Significant rate is ', sig2, '%', sep = ''))+
    scale_fill_manual(values=c("#404080", "#69b3a2"))+
    theme_bw()+gg_theme
}

autocurve.edges2 <-function (graph, start = 0.5)  {
  cm <- count.multiple(graph)
  mut <-is.mutual(graph)  #are connections mutual?
  el <- apply(get.edgelist(graph, names = FALSE), 1, paste,
              collapse = ":")
  ord <- order(el)
  res <- numeric(length(ord))
  p <- 1
  while (p <= length(res)) {
    m <- cm[ord[p]]
    mut.obs <-mut[ord[p]] #are the connections mutual for this point?
    idx <- p:(p + m - 1)
    if (m == 1 & mut.obs==FALSE) { #no mutual conn = no curve
      r <- 0
    }
    else {
      r <- seq(-start, start, length = m)
    }
    res[ord[idx]] <- r
    p <- p + m
  }
  res
}


library(ggplot2)
library(dplyr)
gg_theme= theme(
  title = element_text(size = 30),
  axis.title.x = element_text(size = 28),
  axis.text.x = element_text(size = 18-6, angle=90,hjust=0.95,vjust=0.2), # -6, angle=90,hjust=0.95,vjust=0.2
  axis.title.y = element_text(size = 28),
  axis.text.y = element_text(size = 18),
  plot.subtitle = element_text(size = 24),
  plot.caption = element_text(size = 30), 
  legend.text = element_text(size = 16), 
  legend.key.size = unit(3, 'lines'), 
  # legend.key.height = unit(5, "cm"),
  strip.text = element_text(size = 20), 
  strip.background = element_blank())
