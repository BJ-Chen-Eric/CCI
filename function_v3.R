library(RcppParallel)
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
library(patchwork)
library(progress)
library(ReactomePA)
library(STRINGdb)
string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = 200, input_directory="", protocol="http")

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



detectborder_seucluster <- function(seu, cell_list=cell_list) {
  ann <- data.frame(type=seu$seurat_clusters) %>% mutate(anno='unc')
  ann[ann$type %in% cell_list[[1]], 'anno'] <- names(cell_list)[1]
  # ann[ann$type %in% c(0,3,6,7), 'anno'] <- 'IC'
  ann[ann$type %in% cell_list[[2]], 'anno'] <- names(cell_list)[2]
  # ann[ann$type %in% c(5), 'anno'] <- 'CIS'
  
  ann <- ann[colnames(seu), ]
  seu$anno <- ann$anno
  A <- seu[, colnames(seu) %in% rownames(ann[!ann$anno %in% 'unc', ])]
  
  SpatialDimPlot(A, label = F, group.by = 'anno')
  # DotPlot(hbc, features = c('IL2RG', 'CEACAM6'), group.by = 'anno')
  # DotPlot(hbc, features = c('IL2RG', 'CDH2'), group.by = 'anno')
  # SpatialFeaturePlot(hbc, features = c('IL2RG', 'CDH2'))
  return(prepare_data(seu = A, HVG = F, ori_ident = 'anno', cut_off = 4, detectborder = T))
}


prepare_data <- function(seu, HVG=T, ori_ident, cut_off=4, nHVG=1000, detectborder=T)  {
  A <- as.data.frame(rownames(seu))
  row <- ''
  for(i in c('^AC[0-9]+','^AF[0-9]+','^AL[0-9]+','^AP[0-9]+','^AC[0-9]+'
             ,'^BX[0-9]+','^CU[0-9]+','^FP[0-9]+','^LINC[0-9]+'
             ,'^MT-','^Z[0-9]+', '^IG[HKL]'))  {
    row <- c(row, grep(rownames(seu), pattern = i))
  }
  row <- row[-c(1)]
  seu <- seu[-c(as.numeric(row)),]
  
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
        A_sum[!(A_sum[,1] %in% 0), 'class'] <- paste(j, '_PC_', i, sep = '')
        A_sum <- A_sum[!(is.na(A_sum$class)), ]
        
        B_sum <- as.data.frame(rowSums(adj))
        B_sum[!(B_sum[,1] %in% 0), 'class'] <- paste(i, '_PC_', j, sep = '')
        B_sum <- B_sum[!(is.na(B_sum$class)), ]
        
        meta[rownames(meta) %in% rownames(A_sum), 'xct_ident'] <- A_sum$class
        meta[rownames(meta) %in% rownames(B_sum), 'xct_ident'] <- B_sum$class
      }
    }
    
    result <- as.data.frame(matrix(data = '', nrow = 0, 0))
    for(i in unique(meta$seu))  {
      cell_i <- meta[meta$seu %in% i, ]
      cell_i[cell_i$xct_ident %in% '', 4] <- paste(i, '_CI', sep = '')
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
    rownames(dis) <- rownames(meta_subset_i) # colnames(sender_cell)rownames(meta_subset_i)
    colnames(dis) <- rownames(meta_subset_j) # colnames(receiver_cell)rownames(meta_subset_j)
    # rownames(dis) <- colnames(sender_cell)
    # colnames(dis) <- colnames(receiver_cell)
    
    adj <- ifelse(test = dis<=4, yes = 1, no = 0)
    inn <- ifelse(test = dis>4, yes = 1, no = 0)
    
    B <- data.frame(colnames(seu))
    B$exam <- NA
    for(k in 1:ncol(adj))  {
      B[B[,1] %in% c(adj[adj[, k] ==1, 3] %>% names()), 'exam'] <- 'sender_cell'
    }
    B[B[,1] %in% (colSums(adj) %>% as.data.frame() %>% filter(. != 0) %>% rownames()), 'exam'] <- 'receiver_cell' #
    seu$exam <- B$exam
    exam_plot <- SpatialDimPlot(seu, group.by = 'exam')+gg_theme
    print(exam_plot)
    # return(exam_plot)
    
    
    dim(adj)
    # adj <- as.vector(t(adj))
    # pair <- length(adj)- table(adj)['0']
    
    ## filter 0 expression genes and grep from L_R meta_data
    # B <- ifelse(as.matrix(receiver_cell) ==0, yes = 0, no = 1)
    # B <- as.data.frame(rowSums(B) >= ncol(receiver_cell)*0.05) %>% tibble::rownames_to_column() # express in at least 5 % cells
    
    # B <- B[B[,2] %in% T, ]
    # receiver_cell <- receiver_cell[rownames(receiver_cell) %in% B$rowname, ]
    B <- intersect(rownames(receiver_cell), db$receptor)
    R_count <- db[db$receptor %in% B, c(1,3)]
    
    # B <- ifelse(as.matrix(sender_cell) ==0, yes = 0, no = 1)
    # # B <- as.data.frame(rowSums(B) >= ncol(sender_cell)*0.05) %>% tibble::rownames_to_column()
    # B <- B[B[,2] %in% T, ]
    # sender_cell <- sender_cell[rownames(sender_cell) %in% B$rowname, ]
    B <- intersect(rownames(sender_cell), db$ligand)
    L_count <- db[db$ligand %in% B, 1:2]
    
    
    B <- intersect(L_count$X, R_count$X)
    cat(paste('Number of LR pairs filtered: ', length(B), ', HVGs: ', HVGs, '\n', sep = ''))
    L_count <- L_count[L_count$X %in% B, ]
    R_count <- R_count[R_count$X %in% B, ]
    
    library(Matrix)
    # result <- as(matrix(data = 0, nrow = nrow(L_count), ncol = ncol(receiver_cell)*ncol(sender_cell)), "dgCMatrix")
    
    # result <- matrix(data = 0, nrow = (nrow(L_count)+1), ncol = ncol(receiver_cell)*ncol(sender_cell))
    # rownames(result) <- c(paste(L_count$ligand, R_count$receptor, sep = '_'), 'cell_pairs')
    
    L_count_mat <- matrix(0,ncol = 1, nrow = unique(L_count$ligand) %>% length()) %>% as.data.frame()
    rownames(L_count_mat) <- unique(L_count$ligand)
    for(j in rownames(adj))  {
      L <- data.frame(sender_cell[rownames(sender_cell) %in% unique(L_count$ligand),j], index=1)
      L <- L[rownames(L_count_mat), ]
      L_count_mat[, j] <- L[, 1]
    }
    L_count_mat <- L_count_mat[, -1]
    
    
    R_count_mat <- matrix(0,ncol = 1, nrow = unique(R_count$receptor) %>% length()) %>% as.data.frame()
    rownames(R_count_mat) <- unique(R_count$receptor)
    for(j in  colnames(adj))  {
      R <- data.frame(receiver_cell[rownames(receiver_cell) %in% unique(R_count$receptor),j], index=1)
      R <- R[rownames(R_count_mat), ]
      R_count_mat[, j] <- R[, 1]
    }
    R_count_mat <- R_count_mat[, -1]
    
    subLR_list <- list()
    for(k in seq_len(nrow(L_count)))  {
      A <- L_count_mat[L_count[k, 2], ] %>% as.matrix()
      B <- R_count_mat[R_count[k, 2], ] %>% as.matrix()
      
      subLR_list[[k]] <- t(A) %*% B %>% as('dgCMatrix')
      names(subLR_list)[k] <- paste(L_count[k, 2], R_count[k, 2], sep = "_")
    }
    
    if(isTRUE(Tensor_dec))  {
      cat(paste("Performing Tensor decomposition.\n"))
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
      # subLR_list_n[[i]] <- subLR_list[[i]] %*% adj
      subLR_list_n[[i]] <- subLR_list[[i]][which(adj==1)]
      subLR_list_i[[i]] <- subLR_list[[i]][which(adj==0)]
    }
    names(subLR_list_n) <- names(subLR_list_i) <- names(subLR_list)
    
    out_list <- list(out=list(neigh_propuct=subLR_list_n, inner_propuct=subLR_list_i), p=exam_plot)
    return(out_list)
  }
  
  meta_test <- function(meta_LR_list, sampling=FALSE)  {
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
  
  sender_i = paste(cluster_1, '_inside', sep = '')
  sender = paste(cluster_1, '_close_to_', cluster_2, sep = '')
  receiver_i = paste(cluster_2, '_inside', sep = '')
  receiver = paste(cluster_2, '_close_to_', cluster_1, sep = '')
  cat(paste('Computing ', sender_i, ' communicating to ', receiver_i, '\n', sep = ''))
  one2two <- makeL_R_list_meta(seu = seu, specy = specy, HVGs = nHVGs
                               , sender_i=sender_i, sender = sender
                               , receiver_i = receiver_i, receiver = receiver, Tensor_dec = Tensor_dec)
  gc()
  one2two_test <- meta_test(meta_LR_list = one2two[[1]], sampling=sampling)
  
  cat(paste('Computing ', receiver_i, ' communicating to ', sender_i, '\n', sep = ''))
  two2one <- makeL_R_list_meta(seu = seu, specy = specy, HVGs = nHVGs
                               , sender_i=receiver_i, sender = receiver
                               , receiver_i = sender_i, receiver = sender, Tensor_dec = Tensor_dec)
  # gc()
  two2one_test <- meta_test(meta_LR_list = two2one[[1]], sampling=sampling)
  out_list <- list(one2two, one2two_test, two2one, two2one_test)
  # out_list <- list(one2two, two2one)
  names(out_list) <- c(paste(cluster_1, '_PC', cluster_2, sep = ''),
                       paste(cluster_1, '_PC', cluster_2, '_test', sep = ''),
                       paste(cluster_2, '_PC', cluster_1, sep = ''),
                       paste(cluster_2, '_PC', cluster_1, '_test', sep = ''))
  return(out_list)
}




diff_genes <- function(seu, drug_path, tf_path, nHVG, specy='human'
                       , spatialxct, ori_ident = 'seurat_clusters'
                       , xct_clu='xct_ident', cluster_1, cluster_2 
                       , alpha=1, tf_b=T, drug_b=T)  {
  
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
  
  ggm_comparison <- function(seu, drug_path, tf_path, nHVG, specy='human'
                             , L_R, ori_ident = 'seurat_clusters'
                             , xct_clu='xct_ident', R_border, R_inside 
                             , alpha=1, tf_b, drug_b)  {
    A <- seu[[xct_clu]] %>% as.data.frame() %>% tibble::rownames_to_column('cell')
    seu <- prepare_data(seu = seu, HVG = F, detectborder = F)
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
    
    gene <- union(gene_inside, gene_border) # inside and neighbor cell union geneset
    seu_border <- seu_border[gene, ]
    seu_inside <- seu_inside[gene, ]
    
    devide <- sample(1:ncol(seu_inside), size = ncol(seu_inside)/2, replace = F) # sampling inner cell, building inner cell control
    c1 <- seu_inside[gene, -devide]
    c2 <- seu_inside[gene, devide]
    # print(paste('Control GRN gene number=', nrow(c1), sep=''))
    cat(paste(sub(pattern = '_.*', replacement = '', R_inside), ' cluster control GRN inferece'
              ,'\nGRN gene number=', nrow(seu_inside)
              , '\nControl cell number=', ncol(seu_inside)
              , '\nBorder cell number=', ncol(seu_border), '\n', sep=''))
    # c1 <- c2 <- seu_inside[gene,]
    c1_cell <- ncol(c1)
    c2_cell <- ncol(c2)
    c1 <- c1[['Spatial']]@data %>% as.matrix()
    c2 <- c2[['Spatial']]@data %>% as.matrix()
    # A <- colnames(seu_inside) %>% as.data.frame() 
    # A[A[, 1] %in% colnames(c1), 'test'] <- 'c1'
    # A[A[, 1] %in% colnames(c2), 'test'] <- 'c2'
    # seu_inside$test <- A$test
    control <- scTenifoldNet(X = c1, Y = c2, qc = F,
                             nc_nNet = 10, nc_nCells = 500,
                             td_K = 3, qc_minLibSize = 30, 
                             nc_symmetric = F, td_nDecimal = 10)
    
    seu_border <- seu_border[['Spatial']]@data %>% as.matrix()
    seu_inside <- seu_inside[['Spatial']]@data %>% as.matrix()
    cat(paste(sub(pattern = '_.*', replacement = '', R_inside), ' cluster border/inner GRN inferece'
              ,'\nGRN gene number=', nrow(seu_inside)
              , '\nControl cell number=', ncol(seu_inside)
              , '\nBorder cell number=', ncol(seu_border), '\n', sep=''))
    output <- scTenifoldNet(X = seu_border, Y = seu_inside, qc = F,
                            nc_nNet = 10, nc_nCells = 500,
                            td_K = 3, qc_minLibSize = 30,
                            nc_symmetric = F, td_nDecimal = 10)
    
    B <- output[["diffRegulation"]] %>% filter(p.adj <0.05)
    # B <- sct3000[["SCT_GRN_B_I"]][["diffRegulation"]]
    # B <- B[B$p.adj <0.05, ]
    
    C <- control$diffRegulation %>% filter(p.adj <0.05)
    # C <- sct3000$SCT_GRN_control$diffRegulation
    # C <- C[C$p.adj <0.05, ]
    
    
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
                            , nHVG = nHVG
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
  names(out) <- c(paste(cluster_1, '_PC', sep = '')
                  , paste(cluster_2, '_PC', sep = ''))
  t2 <- Sys.time()
  print(t2-t1)
  return(out)
}


# neighbor_tf <- function(neighbor_deg, inner_deg, Rcis_feather_dir, nTfs=10)  {
#   tf_prediction <- function(degs, nTfs=nTfs)  {
#     geneLists <- list(geneListName=degs)
#     # Motif enrichment analysis:
#     motifEnrichmentTable_wGenes <- cisTarget(geneSets = geneLists
#                                              , motifRankings = motifRankings
#                                              , motifAnnot=motifAnnotations)
#     
#     
#     A <- motifEnrichmentTable_wGenes[!(motifEnrichmentTable_wGenes$TF_highConf %in% ''), ]
#     if(nTfs>nrow(A))  {
#       print('Ask too many Tfs, set to maximun')
#       nTfs <- nrow(A)
#     }
#     A <- A %>% top_n(NES, n = nTfs)
#     B <- sub(A$TF_highConf, pattern = ' \\(.*', replacement = '')
#     tf <- strsplit(B, split = '; ') %>% unlist() %>% unique()
#     return(tf) 
#   }
#   
#   data(motifAnnotations_hgnc)
#   motifAnnotations <- motifAnnotations
#   motifRankings <- importRankings(Rcis_feather_dir) # "~/Analysis/spatial/tf_df/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather"
#   
#   if(length(neighbor_deg) !=0)  {
#     print("performing the TFs prediction of neighbor_deg")
#     neighbor_tf <- tryCatch(
#       {
#         neighbor_tf <- tf_prediction(degs=neighbor_deg, nTfs=nTfs)
#       },
#       error = function(error_message) {
#         message("Fewer than 80% of the genes. No TF will be predicted")
#         return('')
#       }
#     )
#   } 
#   else {
#     print("Zero neighbor_deg, pass TFs identification")
#     neighbor_tf <- ''
#   }
#   if(length(inner_deg) !=0)  {
#     print("Performing the TFs prediction of inner_degs")
#     inner_tf <- tryCatch(
#       {
#         inner_tf <- tf_prediction(degs=inner_deg, nTfs=nTfs)
#       },
#       error = function(error_message) {
#         message("Fewer than 80% of the genes. No TF will be predicted")
#         return('')
#       }
#     )
#   }
#   else {
#     print("Zero neighbor_deg, pass TFs identification")
#     inner_tf <- ''
#   }
#   
#   neighbor_tf <- setdiff(neighbor_tf, inner_tf)
#   
#   return(neighbor_tf)
# }


write_output <- function(xct, n_receptor=10, n_ligand=10, deg, write_path=paste(getwd(), '/', sep = ''))  {
  n <- names(deg)[1]
  out <- xct[[paste(n, '_test', sep = '')]] # Immune to invasive receptor
  out <- out[complete.cases(out), ]
  out$dir <- ifelse(out$fc >=1, yes = 'H', no = 'L')
  out$wcadj <- p.adjust(out$wc)
  out <- out %>% tibble::rownames_to_column('L_R') %>% group_split(dir) %>% as.list() %>% lapply(as.data.frame)
  nei_lr <- c()
  # nei_l <- c()
  for(i in seq_along(out))  {
    A <- out[[i]]
    A <- A[A$wcadj < 0.1  , ] #& out$fc >50
    # row <- nrow(A)
    # fc_cut <- 1.3
    # while(row< 5)  {
    #   A <- out[[i]]
    #   A <- A[A$wcadj < 0.1 & A$fc >fc_cut | A$fc < 2-fc_cut , ]
    #   row <- nrow(A)
    #   fc_cut <- fc_cut-0.05
    #   if(fc_cut==1.05) {break}
    # }
    A <- A %>% slice_max(wcadj, n = n_receptor)
    if(nrow(A) > 0)  {
      # nei_r <- c(nei_r, paste(sub(A[,1], pattern = '.*_', replacement = '') %>% unique(), ' ', 'R_', A$dir, sep = ''))
      # nei_l <- c(nei_l, paste(sub(A[,1], pattern = '_.*', replacement = '') %>% unique(), ' ', 'L_', A$dir, sep = ''))
      nei_lr <- c(nei_lr, paste(A[,1] %>% unique(), '_', A$dir, sep = ''))
    }
  }
  
  if(length(nei_lr) > 0)  {
    nei_lr <- strsplit(nei_lr, split = '_') %>% do.call(what=rbind) %>% as.data.frame()
    colnames(nei_lr) <- c(paste(names(deg)[2] %>% str_remove(pattern = '_PC'), '_L', sep = '')
                          , paste(names(deg)[1] %>% str_remove(pattern = '_PC'), '_R', sep = ''), 'direction')
  }
  if(length(nei_lr) == 0)  {
    nei_lr <- data.frame('', '', '')
    colnames(nei_lr) <-c(paste(names(deg)[2] %>% str_remove(pattern = '_PC'), '_L', sep = '')
                         , paste(names(deg)[1] %>% str_remove(pattern = '_PC'), '_R', sep = ''), 'direction')
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
  
  dg <- deg[[n]]$neighbor_g # Immune to invasive degs
  # dg <- setdiff(deg[[1]]$neighbor_g, deg[[2]]$neighbor_g) # remove the degs in the neighbor cell against it
  
  
  meta_c1 <- data.frame(gene=c(nei_lr[, 2], dg)
                        , type=c(rep('R', length(nei_lr[, 2]))
                                 , rep('DEGs', length(dg))))
  colnames(meta_c1)[1] <- paste(names(deg)[1] %>% str_remove(pattern = '_PC'), 'receptor', sep = '_')
  meta_c1[1:nrow(nei_lr), c(paste(names(deg)[2] %>% str_remove(pattern = '_PC'), '_ligand', sep=''), 'dir')] <- nei_lr[, c(1,3)]
  meta_c1 <- replace(meta_c1, is.na(meta_c1), '')
  
  write.csv(meta_c1[,1], file = paste(write_path, n, '_r_deg_tf_l.csv', sep = '')
              , quote = F, row.names = F)
  
  write.csv(meta_c1, file = paste(write_path, n, '_r_deg_tf_l_label.csv', sep = '')
              , quote = F, row.names = F)
  
  n <- names(deg)[2]
  out <- xct[[paste(n, '_test', sep = '')]] # Immune to invasive receptor
  out <- out[complete.cases(out), ]
  out$dir <- ifelse(out$fc >=1, yes = 'H', no = 'L')
  out$wcadj <- p.adjust(out$wc)
  out <- out %>% tibble::rownames_to_column('L_R') %>% group_split(dir) %>% as.list() %>% lapply(as.data.frame)
  nei_lr <- c()
  # nei_l <- c()
  for(i in seq_along(out))  {
    A <- out[[i]]
    A <- A[A$wcadj < 0.1 , ] #& out$fc >50
    # row <- nrow(A)
    # fc_cut <- 1.3
    # while(row < 5)  {
    #   A <- out[[i]]
    #   A <- A[A$wcadj < 0.1 & A$fc >fc_cut | A$fc < 2-fc_cut , ]
    #   row <- nrow(A)
    #   fc_cut <- fc_cut-0.05
    #   if(fc_cut==1.05) {break}
    # }
    A <- A %>% slice_max(wcadj, n = n_receptor)
    if(nrow(A) > 0)  {
      # nei_r <- c(nei_r, paste(sub(A[,1], pattern = '.*_', replacement = '') %>% unique(), ' ', 'R_', A$dir, sep = ''))
      # nei_l <- c(nei_l, paste(sub(A[,1], pattern = '_.*', replacement = '') %>% unique(), ' ', 'L_', A$dir, sep = ''))
      nei_lr <- c(nei_lr, paste(A[,1] %>% unique(), '_', A$dir, sep = ''))
    }
  }
  if(length(nei_lr) > 0)  {
    nei_lr <- strsplit(nei_lr, split = '_') %>% do.call(what=rbind) %>% as.data.frame()
    colnames(nei_lr) <- c(paste(names(deg)[1] %>% str_remove(pattern = '_PC'), '_L', sep = '')
                          , paste(names(deg)[2] %>% str_remove(pattern = '_PC'), '_R', sep = ''), 'direction')
  }
  if(length(nei_lr) == 0)  {
    nei_lr <- data.frame('', '', '')
    colnames(nei_lr) <- c(paste(names(deg)[1] %>% str_remove(pattern = '_PC'), '_L', sep = '')
                          , paste(names(deg)[2] %>% str_remove(pattern = '_PC'), '_R', sep = ''), 'direction')
  }
  
  
  
  
  dg <- deg[[n]]$neighbor_g
  # dg <- setdiff(deg[[2]]$neighbor_g, deg[[1]]$neighbor_g)
  
  meta_c2 <- data.frame(gene=c(nei_lr[, 2],  dg)
                        , type=c(rep('R', length(nei_lr[, 2]))
                                 , rep('DEGs', length(dg))))
  colnames(meta_c2)[1] <- paste(names(deg)[2] %>% str_remove(pattern = '_PC'), 'receptor', sep = '_')
  meta_c2[1:nrow(nei_lr), c(paste(names(deg)[1] %>% str_remove(pattern = '_PC'), '_ligand', sep = ''), 'dir')] <- nei_lr[, c(1,3)]
  meta_c2 <- replace(meta_c2, is.na(meta_c2), '')
  
  write.csv(meta_c2[,1], file = paste(write_path, n, '_r_deg_tf_l.csv', sep = '')
              , quote = F, row.names = F)
  
  write.csv(meta_c2, file = paste(write_path, n, '_r_deg_tf_l_label.csv', sep = '')
              , quote = F, row.names = F)
  
  out <- list(meta_c1=meta_c1, meta_c2=meta_c2)
  names(out) <- names(deg)
  return(out)
}


infer_pathway <- function(seu_raw, xct_raw, out_path=paste(getwd(), '/', sep = ''), write=T) {
  inferTF <- function(gene_list, TF_number=8, prefix='') {
    geneLists <- list(geneListName=gene_list)
    # Motif enrichment analysis:
    data(motifAnnotations_hgnc)
    motifAnnotations <- motifAnnotations
    motifRankings <- importRankings("~/Analysis/spatial/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather") # "~/Analysis/spatial/tf_df/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather"
    
    motifEnrichmentTable_wGenes <- cisTarget(geneSets = geneLists
                                             , motifRankings = motifRankings
                                             , motifAnnot=motifAnnotations)
    
    A <- motifEnrichmentTable_wGenes[!(motifEnrichmentTable_wGenes$TF_highConf %in% ''), ]
    A <- A %>% top_n(NES, n = TF_number)
    
    TF <- str_remove(A$TF_highConf, pattern = ' \\(.*') %>% strsplit(split = '; ') %>% unlist() %>% unique()
    result <- data.frame(gene=c(gene_list, TF), type=c(rep('Receptor', length(gene_list)), rep('TF', length(TF))))
    # write.csv(result, file = paste('~/', prefix, 'R-TF.csv', sep = ''), quote = F, row.names = F)
    return(result)
  }
  
  individual_pathway <- function(seu, xct, receptor_cluster, ligand_cluster
                                 , out_path=paste(getwd(), '/', sep = ''), write=write) {
    A <- xct[, 1] %>% unlist() %>% clusterProfiler::bitr(fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    x <- enrichPathway(gene=A[, 2], pvalueCutoff = 0.05, readable=TRUE)
    receptor_pathway <- x@result
    receptor_pathway <- receptor_pathway[grep(pattern = paste(xct[xct$type %in% 'DEGs', 1], collapse = '|')
                                              , x = receptor_pathway$geneID), ]
    
    if(xct[!xct[, 3] %in% '', 3] %>% length >0 ) {
      A <- xct[!xct[, 3] %in% '', 3] %>% unlist() %>% clusterProfiler::bitr(fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
      x2 <- enrichPathway(gene=A[, 2], pvalueCutoff = 0.05, readable=TRUE)
      x2@result <- x2@result[x2@result$ID %in% receptor_pathway$ID, ]
      ligand_pathway <- x2@result
      
      receptor_pathway <- receptor_pathway[intersect(receptor_pathway$ID, x2@result$ID), ]
      x@result <- receptor_pathway
    }else{
      cat(paste('No ligand identifed\n'))
      ligand_pathway <- matrix('', nrow = 0, ncol = 9) %>% as.data.frame()
      colnames(ligand_pathway) <- c("ID", "Description", "GeneRatio", "BgRatio","pvalue","p.adjust","qvalue","geneID","Count")
      x2 <- ''
    }
    
    L <- ligand_pathway[, c('ID', 'Description', 'GeneRatio', 'p.adjust', 'geneID')]
    colnames(L) <- c('ID', 'Description', 'Ligand_GeneRatio', 'Ligand_adjust', 'Ligand_geneID')
    R <- receptor_pathway[, c('ID', 'GeneRatio', 'p.adjust', 'geneID')]
    colnames(R) <- c('ID', 'Receptor_GeneRatio', 'Receptor_p.adjust', 'Receptor_geneID')
    A <- merge(L, R, by='ID') 
    if(nrow(A) > 0) {
      A <- A %>% mutate(Ligand_geneID=str_replace_all(Ligand_geneID, pattern = '\\/', replacement = ', ')
                 , Receptor_geneID=str_replace_all(Receptor_geneID, pattern = '\\/', replacement = ', ')
                 , Ligand_adjust=round(Ligand_adjust, digits = 3)
                 , Receptor_p.adjust=round(Receptor_p.adjust, digits = 3))
    }else(A <- A)
    
    pathway_result <- list(ligand=ligand_pathway, receptor=receptor_pathway, merged=A)
    
    
    if(all(receptor_pathway$p.adjust>0.05)) {p1 <- p2 <- p3<- '' } else {
      p1 <- barplot(x, showCategory=20, x="count")+
        ggtitle("Receptor/DRGs significant, shared pathways")+
        theme(plot.title = element_text(size = 36, face = "bold"))
      p2 <- enrichplot::pairwise_termsim(x) %>%
        emapplot(showCategory=20, cex_label_category=0.7)+
        ggtitle("Receptor/DRGs significant, shared\npathways interaction")+
        theme(plot.title = element_text(size = 36, face = "bold"))
      
      R <- receptor_pathway[receptor_pathway$p.adjust < 0.05, 'geneID'] %>% unlist() %>% 
        str_split(pattern = '/') %>% unlist() %>% unique
      
      A <- data.frame(ident=seu$xct_ident) %>% tibble::rownames_to_column('cell')
      A <- merge(A[grep(A$ident, pattern=receptor_cluster), ]
                 , seu[R, A$cell]@assays$Spatial@data %>% t() %>% as.data.frame() %>% tibble::rownames_to_column('cell')
                 , by='cell') %>% group_split(ident) %>% as.list() %>% lapply(as.data.frame)
      
      result <- data.frame(colMeans(A[[1]][, -c(1, 2)]), colMeans(A[[2]][, -c(1, 2)])) %>% tibble::rownames_to_column('gene')
      colnames(result)[-1] <- c(A[[1]][1, 2], A[[2]][1, 2])
      result$fc <- result[, 2]/result[, 3]
      
      A <- merge(result
                 , R %>% unlist() %>% clusterProfiler::bitr(fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
                 , by.x='gene', by.y='SYMBOL')
      gene <- A$fc; names(gene) <- A$ENTREZID
      
      p3 <- cnetplot(x, showCategory = 10, foldChange = gene,
                     categorySize="pvalue", color_category = "#A4ACB8",)+
        ggtitle("Interaction of Receptor/DRGs and pathways")+
        theme(plot.title = element_text(size = 36, face = "bold"))+
        scale_color_gradientn(name = "fold change",colours = c("#52ABFF", "#FFED9A", "#FF66C7"),
                              limits= c(0, 2), breaks=c(0 , 1, 2))
      
    }
    
    
    
    if(all(ligand_pathway$p.adjust>0.05)) {p1_1 <- p2_1 <- p3_1 <- '' } else {
      p1_1 <- barplot(x2, showCategory=20, x="count")+
        ggtitle("Ligands significant, shared pathways")+
        theme(plot.title = element_text(size = 36, face = "bold"))
      p2_1 <- enrichplot::pairwise_termsim(x2) %>%
        emapplot(showCategory=20, cex_label_category=0.7)+
        ggtitle("Ligand significant, shared\npathways interaction")+
        theme(plot.title = element_text(size = 36, face = "bold"))
      
      L <- ligand_pathway[ligand_pathway$p.adjust < 0.05, 'geneID'] %>% unlist() %>%
        str_split(pattern = '/') %>% unlist() %>% unique
      
      A <- data.frame(ident=seu$xct_ident) %>% tibble::rownames_to_column('cell')
      A <- merge(A[grep(A$ident, pattern=ligand_cluster), ]
                 , seu[L, A$cell]@assays$Spatial@data %>% t() %>% as.data.frame() %>% tibble::rownames_to_column('cell')
                 , by='cell') %>% group_split(ident) %>% as.list() %>% lapply(as.data.frame)
      
      result <- data.frame(colMeans(A[[1]][, -c(1, 2)]), colMeans(A[[2]][, -c(1, 2)])) %>% tibble::rownames_to_column('gene')
      colnames(result)[-1] <- c(A[[1]][1, 2], A[[2]][1, 2])
      result$fc <- result[, 2]/result[ , 3]
      
      A <- merge(result
                 , L %>% unlist() %>% clusterProfiler::bitr(fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
                 , by.x='gene', by.y='SYMBOL')
      gene <- A$fc; names(gene) <- A$ENTREZID
      
      p3_1 <- cnetplot(x2, showCategory = 10, foldChange = gene,
                       categorySize="pvalue", color_category = "#A4ACB8",)+
        ggtitle("Interaction of ligands and pathways")+
        theme(plot.title = element_text(size = 36, face = "bold"))+
        scale_color_gradientn(name = "fold change",
                              colours = c("#52ABFF", "#FFED9A", "#FF66C7"),
                              limits= c(0, 2),
                              breaks=c(0 , 1, 2))
      
    }
    
    
    p <- list(receptor_bar=p1, receptor_ema=p2, receptor_cne=p3, ligand_bar=p1_1, ligand_ema=p2_1, ligand_cne=p3_1)
    
    
    if(identical(character(0), receptor_pathway$geneID)) {
      R <- data.frame(gene='', type='')
    }else{
      R <- receptor_pathway$geneID %>% unlist() %>% str_split(pattern = '/') %>% unlist() %>% unique %>%
        inferTF(TF_number=8, prefix='')
      R[R$gene %in% xct[xct$type %in% 'DEGs', 1], 'type'] <- 'DEG'
    }
    
    
    if(identical(character(0), ligand_pathway$geneID)) {
      L <- data.frame(gene='', type='')
    }else{
      L <- ligand_pathway$geneID %>% unlist() %>% str_split(pattern = '/') %>% unlist() %>% unique %>%
        inferTF(TF_number=8, prefix='')
      L[L$type %in% 'Receptor', 'type'] <- 'Ligand'
    }
    
    
    if(isTRUE(write)) {
      for(j in 1:3) {
        png(filename = paste(out_path, colnames(xct)[1], '_', names(p)[j], '.png', sep = ''), res = 600, width = 9600, height = 5400)
        print(p[[j]])
        dev.off()
      }
      
      for(j in 4:6) {
        png(filename = paste(out_path, colnames(xct)[3], '_', names(p)[j], '.png', sep = ''), res = 600, width = 9600, height = 5400)
        print(p[[j]])
        dev.off()
      }
      
      write.csv(R %>% arrange(type), paste(out_path, colnames(xct)[1], '_Receptor_DRG.csv', sep = ''), row.names = F)
      write.csv(L %>% arrange(type), paste(out_path, colnames(xct)[3], '.csv', sep = ''), row.names = F)
      write.csv(pathway_result[['merged']] 
                , paste(out_path, colnames(xct)[1], '_', colnames(xct)[3], '_merged_pathway', '.csv', sep = '')
                , row.names = F)
    }
    
    
    gene <- list(R=R, L=L)
    names(gene) <- c(colnames(xct)[c(1, 3)])
    
    out <- list(pathway=pathway_result, plot=p, gene_TF=gene)
    return(out)
  }
  
  g1 <- paste('^', names(xct_raw)[1] %>% str_sub(1, 3), sep = '')
  g2 <- paste('^', names(xct_raw)[2] %>% str_sub(1, 3), sep = '')
  
  result <- list(individual_pathway(seu = seu_raw, xct = xct_raw[[1]], receptor_cluster = g1, ligand_cluster=g2
                                    , out_path=out_path, write = write)
                 , individual_pathway(seu = seu_raw, xct = xct_raw[[2]], receptor_cluster = g2, ligand_cluster=g1
                                      , out_path=out_path, write = write))
  
  names(result) <- names(xct_raw)
  
  return(result)
}


string_network <- function(xct_result, out_path=paste(getwd(), '/', sep = '')) {
  for(i in seq_along(xct_result)) {
    n <- names(xct_result)[i]
    
    test <- xct_result[[i]][['gene_TF']]
    for(j in seq_along(test)) {
      gene <- data.frame((test[[j]]), pvalue=0, logFC=0)
      if(gene$gene %>% length !=1) {
        example1_mapped <- string_db$map(gene, "gene", removeUnmappedRows = TRUE)
        
        example1_mapped_sig <- string_db$add_diff_exp_color(subset(example1_mapped, pvalue >= 0, logFcColStr="logFC" ))
        example1_mapped_sig[example1_mapped_sig$gene %in% gene[gene$type %in% 'DEG', 'gene'], 'color'] <- '#D870AD'
        example1_mapped_sig[example1_mapped_sig$gene %in% gene[gene$type %in% 'Receptor', 'gene'], 'color'] <- '#6DBFE1'
        example1_mapped_sig[example1_mapped_sig$gene %in% gene[gene$type %in% 'TF', 'gene'], 'color'] <- '#93E191'
        example1_mapped_sig[example1_mapped_sig$gene %in% gene[gene$type %in% 'Ligand', 'gene'], 'color'] <- '#AC92ED'
        
        payload_id <- string_db$post_payload(example1_mapped_sig$STRING_id,
                                             colors=example1_mapped_sig$color, )
        # display a STRING network png with the "halo"
        png(filename = paste(out_path, names(test)[j], '_STRING.png', sep = ''), res = 600, width = 9600, height = 5400)
        string_db$plot_network(example1_mapped$STRING_id, payload_id=payload_id, network_flavor='confidence')
        dev.off()
      }else{cat(paste('No gene to infer PPI\n'));next}
      
    }
  }
}

gg_theme= theme(
  title = element_text(size = 30),
  axis.title.x = element_text(size = 28),
  axis.text.x = element_text(size = 28), # -6, angle=90,hjust=0.95,vjust=0.2
  axis.title.y = element_text(size = 28),
  axis.text.y = element_text(size = 28),
  plot.subtitle = element_text(size = 24),
  plot.caption = element_text(size = 30), 
  legend.text = element_text(size = 26), 
  legend.key.size = unit(2.5, 'lines'), 
  legend.key.height = unit(1.25, "cm"),
  strip.text = element_text(size = 20), 
  strip.background = element_blank())
