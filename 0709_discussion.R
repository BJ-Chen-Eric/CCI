DC_xct <- readRDS('~/Analysis/spatial/DC_xct.RData')
hbc_xct_old <- readRDS('~/Analysis/spatial/hbc_xct.RData')

A <- hbc_xct$Immune_to_Invasive_test$ks
A <- hbc_xct$Immune_to_Invasive_test$wc
A <- DC_xct$`2_close_to_0_test`
A <- DC_xct$`0_close_to_2_test`


ggplot(A)+
  geom_point(aes(x=fc, y=-log10(ks)))+
  geom_hline(yintercept=-log10(0.05), color='#E9573E', linetype='dashed', size=1.5)+
  geom_vline(xintercept=c(1), color='#E9573E', linetype='dashed', size=1.5)+
  labs(title = paste('Immune neighbor cells LR pairs')
       , subtitle = paste(nrow(A), ' LR pairs are examined (ks test)', sep = ''))
  scale_y_continuous(labels = c(seq(0,100,10)), breaks = c(seq(0,100,10)), limits = c(0,100))+
  theme_bw()+gg_theme


A <- hbc_xct$Invasive_to_Immune_test$ks
A <- hbc_xct$Invasive_to_Immune_test$wc
A <- dc_xct$`2_close_to_0_test`$ks
A <- dc_xct$`2_close_to_0_test`$wc

colnames(A)[2] <- 'postime'


ggplot(A)+
  geom_point(aes(x=fc, y=postime))+
  geom_hline(yintercept=c(50), color='#E9573E', linetype='dashed', size=1.5)+
  geom_vline(xintercept=c(50), color='#E9573E', linetype='dashed', size=1.5)+
  labs(title = paste('Invasive neighbor cells LR pairs')
       , subtitle = paste(nrow(A), ' LR pairs are examined (ks test)', sep = ''))+
  scale_y_continuous(labels = c(seq(0,100,10)), breaks = c(seq(0,100,10)), limits = c(0,100))+
  theme_bw()+gg_theme




# ks/wc convergence LR pairs ----------------------------------------------
DC_xct <- readRDS('~/Analysis/spatial/DC_xct.RData')
A <- DC_xct$`2_close_to_0_test` %>% tibble::rownames_to_column('L_R')
B <- DC_xct$`0_close_to_2_test` %>% tibble::rownames_to_column('L_R')

out <- matrix(0,0,0) %>% as.data.frame()
x=0
for(i in c(seq(10, nrow(A)-5, 10), nrow(A)))  {
  sub_ks <- A %>% top_n(ks, n = i)
  sub_wc <- B %>% top_n(wc, n = i)
  x=x+1
  out[x, 'shared'] <- intersect(sub_ks$L_R, sub_wc$L_R) %>% length()
  out[x, 'ks_only'] <- setdiff(sub_ks$L_R, sub_wc$L_R) %>% length()
  out[x, 'wc_only'] <- setdiff(sub_wc$L_R, sub_ks$L_R) %>% length()
}
out$x <- c(seq(10, nrow(A)-5, 10), nrow(A))

out %>% tidyr::gather(key, value, c(1:3)) %>% 
  ggplot()+
  geom_line(aes(x=x, y=value, color=key), linewidth=2)+
  geom_point(aes(x=x, y=value, color=key), size=5)+
  labs(title = 'ks test and wc test shared LR pairs', x='Top LR pairs in ks and wc test'
       , y='LR_pairs')+
  scale_color_discrete(name=NULL)+
  # labs(title = paste('sender: ', cluster_1, ', receiver: ', cluster_2, sep = ''))+
  scale_x_continuous(breaks = c(seq(10, nrow(A)-5, 10), nrow(A)))+
  theme_bw()+gg_theme

# LR overlap --------------------------------------------------------------
# DC_xct <- spatialXct(seu = DC_subset, specy = 'human', cluster_1 = '0', cluster_2 = '2', Tensor_dec = T)


db <- read.csv('~/Analysis/spatial/data_baseLR.csv')

### sender =0, receiver =2
LR_overlap <- function(cluster_1, cluster_2, seu, nHVGs=3000)  {
  sender_i = paste(cluster_1, '_inside', sep = '')
  sender = paste(cluster_1, '_close_to_', cluster_2, sep = '')
  receiver_i = paste(cluster_2, '_inside', sep = '')
  receiver = paste(cluster_2, '_close_to_', cluster_1, sep = '')
  db <- read.csv('~/Analysis/spatial/data_baseLR.csv')
  
  db <- distinct(db, ligand, .keep_all = T)[, 1:3]
  assay <- DefaultAssay(seu)
  cluster_s <- subset(seu, subset= xct_ident %in% c(sender_i, sender))
  cluster_s_all <- FindVariableFeatures(cluster_s, nfeature=nHVGs) %>% VariableFeatures()
  
  cluster_r <- subset(seu, subset= xct_ident %in% c(receiver_i, receiver))
  cluster_r_all <- FindVariableFeatures(cluster_r, nfeature=nHVGs) %>% VariableFeatures()
  
  
  out <- matrix(0,nrow = 0,0) %>% as.data.frame()
  x=0
  for(i in seq(100,nHVGs,100))  {
    cluster_s <- cluster_s_all[1:i]
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
    
    cluster_r <- cluster_r_all[1:i]
    
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
    
    
    R_count <- intersect(rownames(receiver_cell), db$receptor)
    
    L_count <- intersect(rownames(sender_cell), db$ligand)
    
    
    x=x+1
    out[x, c('R_count', 'L_count')] <- c(length(R_count), length(L_count))
  }
  out$HVGs <- seq(100,nHVGs,100)
  
  plot1 <- ggplot(out)+
    geom_line(aes(x=HVGs, y=R_count))+
    geom_point(aes(x=HVGs, y=R_count))+
    labs(title = 'LR pairs overlap among\nnumbers of HVGs'
         , subtitle = paste('sender: ', cluster_1, ', receiver: ', cluster_2, sep = '')
         , x='numbers of HVGs'
         , y='R_count')+
    scale_x_continuous(breaks = seq(100,nHVGs,100))+
    theme_bw()+gg_theme+
    theme(axis.text.x = element_text(size = 18-6, angle=90,hjust=0.95,vjust=0.2))
  
  plot2 <- ggplot(out)+
    geom_line(aes(x=HVGs, y=L_count))+
    geom_point(aes(x=HVGs, y=L_count))+
    labs(title = 'LR pairs overlap among\nnumbers of HVGs'
         , subtitle = paste('sender: ', cluster_1, ', receiver: ', cluster_2, sep = '')
         , x='numbers of HVGs'
         , y='L_count')+
    scale_x_continuous(breaks = seq(100,nHVGs,100))+
    theme_bw()+gg_theme+
    theme(axis.text.x = element_text(size = 18-6, angle=90,hjust=0.95,vjust=0.2))
  return(cowplot::plot_grid(plot1, plot2, ncol = 2))
}

LR_overlap(cluster_1 = '0', cluster_2 = '2', seu = DC_subset, nHVGs = 5000)
LR_overlap(cluster_1 = '2', cluster_2 = '0', seu = DC_subset, nHVGs = 5000)
LR_overlap(cluster_1 = 'Immune', cluster_2 = 'Invasive', seu = hbc, nHVGs = 5000)
LR_overlap(cluster_1 = 'Invasive', cluster_2 = 'Immune', seu = hbc, nHVGs = 5000)




shareLR_ST <- function(cluster_1_1, cluster_1_2, cluster_2_1, cluster_2_2, seu1_1, seu1_2, title)  {
  LR_genes <- function(cluster_1, cluster_2, seu, nHVGs=3000)  {
    sender_i = paste(cluster_1, '_inside', sep = '')
    sender = paste(cluster_1, '_close_to_', cluster_2, sep = '')
    receiver_i = paste(cluster_2, '_inside', sep = '')
    receiver = paste(cluster_2, '_close_to_', cluster_1, sep = '')
    db <- read.csv('~/Analysis/spatial/data_baseLR.csv')
    assay <- DefaultAssay(seu)
    cluster_s <- subset(seu, subset= xct_ident %in% c(sender_i, sender))
    cluster_s <- FindVariableFeatures(cluster_s, nfeature=nHVGs) %>% VariableFeatures()
    
    cluster_r <- subset(seu, subset= xct_ident %in% c(receiver_i, receiver))
    cluster_r <- FindVariableFeatures(cluster_r, nfeature=nHVGs) %>% VariableFeatures()
    
    
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
    
    
    R_count <- intersect(rownames(receiver_cell), db$receptor)
    
    L_count <- intersect(rownames(sender_cell), db$ligand)
    
    B <- list(R=R_count, L=L_count)
    
    return(B)
  }
  out <- matrix(0,nrow = 0,0) %>% as.data.frame()
  x=0
  for(i in seq(100,5000,100))  {
    imm_inva1 <- LR_genes(cluster_1 = cluster_1_1, cluster_2 = cluster_2_1, seu = seu1_1, nHVGs = i)
    imm_inva2 <- LR_genes(cluster_1 = cluster_1_2, cluster_2 = cluster_2_2, seu = seu1_2, nHVGs = i)
    
    x=x+1
    out[x, 'Shared_R_count'] <- intersect(imm_inva1$R, imm_inva2$R) %>% length()
    out[x, 'data1_R_count'] <- imm_inva1$R %>% length()
    out[x, 'data2_R_count'] <- imm_inva2$R %>% length()
    out[x, 'Shared_L_count'] <- intersect(imm_inva1$L, imm_inva2$L) %>% length()
    out[x, 'data1_L_count'] <- imm_inva1$L %>% length()
    out[x, 'data2_L_count'] <- imm_inva2$L %>% length()
  }
  out$HVGs <- seq(100,5000,100)
  
  plot1 <- out[,c(1:3, 7)] %>% tidyr::gather(key, value, c(1:3)) %>% 
    ggplot()+
    geom_line(aes(x=HVGs, y=value, color=key), linewidth=2)+
    geom_point(aes(x=HVGs, y=value, color=key), size=5)+
    labs(title = 'Overlapped L_R pairs among\ndifferent clusters'
         , subtitle = title
         , x='numers of HVGs'
         , y='LR_pairs')+
    scale_color_manual(values = c(data1_R_count='#4FC0E8', data2_R_count='#EC87BF', Shared_R_count='#86DC7B'))+
    theme_bw()+gg_theme
  
  plot2 <- out[, c(4,5,6,7)] %>% tidyr::gather(key, value, c(1:3)) %>% 
    ggplot()+
    geom_line(aes(x=HVGs, y=value, color=key), linewidth=2)+
    geom_point(aes(x=HVGs, y=value, color=key), size=5)+
    labs(title = 'Overlapped L_R pairs among\ndifferent clusters'
         , subtitle = title
         , x='numers of HVGs'
         , y='LR_pairs')+
    scale_color_manual(values = c(data1_L_count='#4FC0E8', data2_L_count='#EC87BF', Shared_L_count='#86DC7B'))+
    theme_bw()+gg_theme
  return(list(p=cowplot::plot_grid(plot1, plot2), out=out))
}

A <- shareLR_ST(cluster_1_1=0, cluster_1_2=0, cluster_2_1 = 1, cluster_2_2 = 2
                , seu1_1 = DC_subset, seu1_2 = DC_subset
                , title = paste('Data1: ', 'sender-', 0, ', receiver-', 1
                                , '\n', 'Data2: ', 'sender-', 0, ', receiver-', 2, sep = ''))

A[['out']][, c(1:3, 7)] %>% tidyr::gather(key, value, c(1:3)) %>% 
  ggplot()+
  geom_line(aes(x=HVGs, y=value, color=key), linewidth=2)+
  geom_point(aes(x=HVGs, y=value, color=key), size=5)+
  labs(title = 'Overlapped L_R pairs among different clusters'
       , subtitle = paste('Data1: ', 'sender-', 0, ', receiver-', 1
                          , '\n', 'Data2: ', 'sender-', 0, ', receiver-', 2, sep = '')
       , x='numers of HVGs'
       , y='LR_pairs')+
  scale_color_manual(values = c(data1_R_count='#4FC0E8', data2_R_count='#EC87BF', Shared_R_count='#86DC7B'))+
  theme_bw()+gg_theme

B <- shareLR_ST(cluster_1_1=0, cluster_1_2='Immune', cluster_2_1 = 1, cluster_2_2 = 'Invasive'
                , seu1_1 = DC_subset, seu1_2 = hbc
                , title = paste('Data1: Ductal carcinoma', 'sender-', 0, ', receiver-', 1
                                , '\n', 'Data2: human breast cancer', 'sender-', 'Immune', ', receiver-', 'Invasive', sep = ''))
B[['out']] %>% tidyr::gather(key, value, c(1:3)) %>% 
  ggplot()+
  geom_line(aes(x=HVGs, y=value, color=key), linewidth=2)+
  geom_point(aes(x=HVGs, y=value, color=key), size=5)+
  labs(title = 'Overlapped L_R pairs among different clusters'
       , subtitle = paste('Data1: Ductal carcinoma', 'sender-', 0, ', receiver-', 1
                          , '\n', 'Data2: human breast cancer', 'sender-', 'Immune', ', receiver-', 'Invasive', sep = '')
       , x='numers of HVGs'
       , y='LR_pairs')+
  scale_color_manual(values = c(data1='#4FC0E8', data2='#EC87BF', Shared='#86DC7B'))+
  theme_bw()+gg_theme

C <- A[['out']]
C$data1 <- C$data1-C$Shared
C$data2 <- C$data2-C$Shared
C$all <- C$data2+C$data1+C$Shared
D <- data.frame(C[, c(1,2,3)]/C[,'all']*100, HVGs=C[, 4]) %>% tidyr::gather(key = color, value, c(1:3)) 
D$color <- factor(D$color, levels = unique(D$color) %>% rev())
ggplot(D)+
  geom_bar(aes(x=HVGs, y=value, fill=color), color='black', stat = 'identity')+
  labs(title = 'unique and shared LR pairs in two exam groups'
       , subtitle = paste('Data1: ', 'sender-', 0, ', receiver-', 1
                          , '\n', 'Data2: ', 'sender-', 0, ', receiver-', 2, sep = '')
       , x='HVGs'
       , y='Percentage (%)')+
  scale_y_continuous(breaks = c(seq(0, 100, 10)))+
  scale_fill_manual(values = c('#4FC0E8', '#EC87BF', '#86DC7B'))+
  # labs(title = paste('sender: ', cluster_1, ', receiver: ', cluster_2, sep = ''))+
  theme_linedraw()+gg_theme


C <- B[['out']]
C$data1 <- C$data1-C$Shared
C$data2 <- C$data2-C$Shared
C$all <- C$data2+C$data1+C$Shared
D <- data.frame(C[, c(1,2,3)]/C[,'all']*100, HVGs=C[, 4]) %>% tidyr::gather(key = color, value, c(1:3)) 
D$color <- factor(D$color, levels = unique(D$color) %>% rev())
ggplot(D)+
  geom_bar(aes(x=HVGs, y=value, fill=color), color='black', stat = 'identity')+
  labs(title = 'unique and shared LR pairs in two exam groups'
       , subtitle = paste('Data1: Ductal carcinoma', 'sender-', 0, ', receiver-', 1
                          , '\n', 'Data2: human breast cancer', 'sender-', 'Immune', ', receiver-', 'Invasive', sep = '')
       , x='HVGs'
       , y='Percentage (%)')+
  scale_y_continuous(breaks = c(seq(0, 100, 10)))+
  scale_fill_manual(values = c('#4FC0E8', '#EC87BF', '#86DC7B'))+
  # labs(title = paste('sender: ', cluster_1, ', receiver: ', cluster_2, sep = ''))+
  theme_linedraw()+gg_theme


out <- matrix(0,nrow = 0,0) %>% as.data.frame()
x=0
for(i in seq(100,5000,100))  {
  imm_inva2 <- LR_genes(cluster_1 = '0', cluster_2 = '2', seu = DC_subset, nHVGs = i)
  imm_inva1 <- LR_genes(cluster_1 = 'Immune', cluster_2 = 'Invasive', seu = hbc, nHVGs = i)
  
  x=x+1
  out[x, 'LR_pairs'] <- intersect(imm_inva1$LR, imm_inva2$LR) %>% length()
  out[x, 'data1'] <- imm_inva2 %>% nrow()
  out[x, 'data2'] <- imm_inva1 %>% nrow()
}
out$HVGs <- seq(100,5000,100)

out %>% tidyr::gather(key, value, c(1:3)) %>% 
  ggplot()+
  geom_line(aes(x=HVGs, y=value, color=key))+
  geom_point(aes(x=HVGs, y=value, color=key))+
  labs(title = 'ks test and wc test shared LR pairs\nsender: Immune, receiver: Invasive'
       , subtitle = paste('data1: ', 'Ductal caricnoma'
                          , '\n', 'data2: ', 'Human Breast Cancer', sep = '')
       , x='numers of HVGs'
       , y='LR_pairs')+
  theme_bw()+gg_theme




# tensor efffect ----------------------------------------------------------
### no-tensor subsampling
# cluster_1 = '0'; cluster_2 = '2'; seu = DC_subset
# 
# sender_i = paste(cluster_1, '_inside', sep = '')
# sender = paste(cluster_1, '_close_to_', cluster_2, sep = '')
# receiver_i = paste(cluster_2, '_inside', sep = '')
# receiver = paste(cluster_2, '_close_to_', cluster_1, sep = '')
# 
# one2two <- makeL_R_list_meta(seu = seu, specy = 'human'
#                              , sender_i=sender_i, sender = sender
#                              , receiver_i = receiver_i, receiver = receiver, Tensor_dec = F)
# gc()
# one2two_test <- meta_test(one2two, sampling = F)
# two2one <- makeL_R_list_meta(seu = seu, specy = 'human'
#                              , sender_i=receiver_i, sender = receiver
#                              , receiver_i = sender_i, receiver = sender, Tensor_dec = F)
# gc()
# two2one_test <- meta_test(two2one, sampling = F)
# out_list_noT_whole <- list(one2two, one2two_test, two2one, two2one_test)
# names(out_list_noT_whole) <- c(paste(cluster_1, '_close_to_', cluster_2, sep = ''), 
#                      paste(cluster_1, '_close_to_', cluster_2, '_test', sep = ''), 
#                      paste(cluster_2, '_close_to_', cluster_1, sep = ''), 
#                      paste(cluster_2, '_close_to_', cluster_1, '_test', sep = ''))
# 
# 
# cluster_1 = '0'; cluster_2 = '2'; seu = DC_subset
# 
# sender_i = paste(cluster_1, '_inside', sep = '')
# sender = paste(cluster_1, '_close_to_', cluster_2, sep = '')
# receiver_i = paste(cluster_2, '_inside', sep = '')
# receiver = paste(cluster_2, '_close_to_', cluster_1, sep = '')
# 
# one2two <- makeL_R_list_meta(seu = seu, specy = 'human'
#                              , sender_i=sender_i, sender = sender
#                              , receiver_i = receiver_i, receiver = receiver, Tensor_dec = T)
# gc()
# one2two_test <- meta_test(one2two, sampling = F)
# two2one <- makeL_R_list_meta(seu = seu, specy = 'human'
#                              , sender_i=receiver_i, sender = receiver
#                              , receiver_i = sender_i, receiver = sender, Tensor_dec = T)
# gc()
# two2one_test <- meta_test(two2one, sampling = F)
# out_list_T_whole <- list(one2two, one2two_test, two2one, two2one_test)
# names(out_list_T_whole) <- c(paste(cluster_1, '_close_to_', cluster_2, sep = ''), 
#                                paste(cluster_1, '_close_to_', cluster_2, '_test', sep = ''), 
#                                paste(cluster_2, '_close_to_', cluster_1, sep = ''), 
#                                paste(cluster_2, '_close_to_', cluster_1, '_test', sep = ''))
# 

out_list_T_whole <- spatialXct(seu = DC_subset, specy = 'human', cluster_1 = '0', cluster_2 = '2', Tensor_dec = T, sampling = F)
out_list_noT_whole <- spatialXct(seu = DC_subset, specy = 'human', cluster_1 = '0', cluster_2 = '2', Tensor_dec = F, sampling = F)

A <- out_list_T_whole$`0_close_to_2_test`
A <- out_list_T_whole$`2_close_to_0_test`
A <- hbc_xct$Immune_close_to_Invasive_test
A <- hbc_xct$Invasive_close_to_Immune_test
A <- hbc_xct_raw$Invasive_close_to_Immune_test
A <- DC_xct$Immune_close_to_IC_test
A <- DC_xct$IC_close_to_Immune_test
A <- DC_xct$CIS_close_to_Immune_test


A <- A[!is.na(A$wc), ]
A$logfc <- log(A$fc)
A$ks.adj <- A$ks %>% p.adjust()
A$logks <- -log10(A$ks.adj)
A$wc.adj <- A$wc %>% p.adjust()
A$logwc <- -log10(A$wc.adj)

A[, 'label_wc'] <- '0'
A[A$wc.adj <=0.1 & A$logfc >0, 'label_wc'] <- '1'
A[A$wc.adj <=0.1 & A$logfc <0, 'label_wc'] <- '-1'
A[A$wc.adj <=0.1, 'repel_wc'] <- A[A$wc.adj <=0.1, ] %>% rownames()


A[, 'label_ks'] <- '0'
A[A$ks.adj <=0.05 & A$logfc >0, 'label_ks'] <- '1'
A[A$ks.adj <=0.05 & A$logfc <0, 'label_ks'] <- '-1'
A[A$ks.adj <=0.05, 'repel_ks'] <- A[A$ks.adj <=0.05, ] %>% rownames()

A[A$fc >1.2 , ] %>% top_n(wt = wc.adj, n = 20) %>% rownames() %>% sub(pattern= '_.*', replacement = '') %>% unique()
A[A$fc >1.2 , ] %>% top_n(wt = wc.adj, n = 20) %>% rownames() %>% sub(pattern= '.*_', replacement = '') %>% unique()


g <- A %>% tibble::rownames_to_column('L_R')


ggplot(A, aes(x=logfc, y=logwc, color=label_wc))+
  geom_point(size=5)+
  geom_hline(yintercept=c(-log10(0.1)), color='#E9573E', linetype='dashed', linewidth=1.5)+
  geom_vline(xintercept=c(0), color='#E9573E', linetype='dashed', linewidth=1.5)+
  labs(title = paste('Invasive neighbor cells LR pairs, Tensor denoised')
       , subtitle = paste(nrow(A), ' LR pairs are examined (wc test)', sep = '')
       , y= '-log(p-value)')+
  # scale_y_continuous(labels = c(seq(0,100,10)), breaks = c(seq(0,100,10)), limits = c(0,100))+
  geom_text_repel(mapping = aes(label=repel_wc), size=5, nudge_x = 0.01)+ # , nudge_x = 0.1, direction = "y", hjust = "left"
  theme_bw()+gg_theme


ggplot(A, aes(x=logfc, y=logks, color=label_ks))+
  geom_point(size=5)+
  geom_hline(yintercept=c(-log10(0.1)), color='#E9573E', linetype='dashed', linewidth=1.5)+
  geom_vline(xintercept=c(0), color='#E9573E', linetype='dashed', linewidth=1.5)+
  labs(title = paste('Invasive neighbor cells LR pairs, Tensor denoised')
       , subtitle = paste(nrow(A), ' LR pairs are examined (ks test)', sep = '')
       , y= '-log(p-value)')+
  # scale_y_continuous(labels = c(seq(0,100,10)), breaks = c(seq(0,100,10)), limits = c(0,100))+
  geom_text_repel(mapping = aes(label=repel_ks), size=5, nudge_x = 0.01)+ # , nudge_x = 0.1, direction = "y", hjust = "left"
  theme_bw()+gg_theme


B <- out_list_noT_whole$`0_close_to_2_test` %>% tibble::rownames_to_column('L_R')
B <- out_list_noT_whole$`2_close_to_0_test` %>% tibble::rownames_to_column('L_R')
B$logfc <- log(B$fc)
B$logks <- -log10(B$ks)
B$logwc <- -log10(B$wc)

B <- merge(B, g[, c(1,8:11)], by = 'L_R', all = T)


ggplot(B, aes(x=logfc, y=logwc, color=label_wc))+
  geom_point(size=5)+
  geom_hline(yintercept=c(-log10(0.05)), color='#E9573E', linetype='dashed', linewidth=1.5)+
  geom_vline(xintercept=c(0), color='#E9573E', linetype='dashed', linewidth=1.5)+
  labs(title = paste('Invasive neighbor cells LR pairs, Raw data')
       , subtitle = paste(nrow(B), ' LR pairs are examined (wc test)', sep = '')
       , y= '-log(p-value)')+
  # scale_y_continuous(labels = c(seq(0,100,10)), breaks = c(seq(0,100,10)), limits = c(0,100))+
  geom_text_repel(mapping = aes(label=repel_wc), size=5, nudge_x = 0.01)+ # , nudge_x = 0.1, direction = "y", hjust = "left"
  theme_bw()+gg_theme

### tensor whole data test
whole_sub_vol <- function(whole_test, sub_test, cluster)  {
  A <- whole_test
  
  A$logfc <- log(A$fc)
  A$logks <- -log10(A$ks)
  A$logwc <- -log10(A$wc)
  
  A[, 'label_wc'] <- '0'
  A[A$wc <=0.05 & A$logfc >0, 'label_wc'] <- '1'
  A[A$wc <=0.05 & A$logfc <0, 'label_wc'] <- '-1'
  A[A$wc <=0.05, 'repel_wc'] <- A[A$wc <=0.05, ] %>% rownames()
  
  
  A[, 'label_ks'] <- '0'
  A[A$ks <=0.05 & A$logfc >0, 'label_ks'] <- '1'
  A[A$ks <=0.05 & A$logfc <0, 'label_ks'] <- '-1'
  A[A$ks <=0.05, 'repel_ks'] <- A[A$ks <=0.05, ] %>% rownames()
  
  g <- A %>% tibble::rownames_to_column('L_R')
  
  
  
  ks_w <- ggplot(A, aes(x=logfc, y=logks, color=label_ks))+
    geom_point(size=5)+
    geom_hline(yintercept=c(-log10(0.05)), color='#E9573E', linetype='dashed', linewidth=1.5)+
    geom_vline(xintercept=c(0), color='#E9573E', linetype='dashed', linewidth=1.5)+
    labs(title = paste(cluster, 'neighbor cells LR pairs, perform statistic by whole data')
         , subtitle = paste(nrow(A), ' LR pairs are examined (ks test)', sep = '')
         , y= '-log(p-value)')+
    # scale_y_continuous(labels = c(seq(0,100,10)), breaks = c(seq(0,100,10)), limits = c(0,100))+
    geom_text_repel(mapping = aes(label=repel_ks), size=5, nudge_x = 0.01)+ # , nudge_x = 0.1, direction = "y", hjust = "left"
    theme_bw()+gg_theme
  
  
  
  wc_w <- ggplot(A, aes(x=logfc, y=logwc, color=label_wc))+
    geom_point(size=5)+
    geom_hline(yintercept=c(-log10(0.05)), color='#E9573E', linetype='dashed', linewidth=1.5)+
    geom_vline(xintercept=c(0), color='#E9573E', linetype='dashed', linewidth=1.5)+
    labs(title = paste(cluster, 'neighbor cells LR pairs, perform statistic by whole data')
         , subtitle = paste(nrow(A), ' LR pairs are examined (wc test)', sep = '')
         , y= '-log(p-value)')+
    geom_text_repel(mapping = aes(label=repel_wc), size=5, nudge_x = 0.01)+
    # scale_y_continuous(labels = c(seq(0,100,10)), breaks = c(seq(0,100,10)), limits = c(0,100))+
    theme_bw()+gg_theme
  
  
  
  # A <- sub_test$ks
  # colnames(A)[2] <- 'postime'
  # A <- merge(A, g[, c(1,8:11)], by = 'L_R', all = T)
  # 
  # ks_sub <- ggplot(A, aes(x=fc, y=postime, color=label_ks, label=L_R))+
  #   geom_point(size=5)+
  #   geom_hline(yintercept=c(50), color='#E9573E', linetype='dashed', size=1.5)+
  #   geom_vline(xintercept=c(50), color='#E9573E', linetype='dashed', size=1.5)+
  #   labs(title = paste(cluster, 'neighbor cells LR pairs, perform statistic by subsampling')
  #        , subtitle = paste(nrow(A), ' LR pairs are examined (ks test)', sep = '')
  #        , x= 'Over/under expression direction', y= 'Siginificant times')+
  #   scale_y_continuous(labels = c(seq(0,100,10)), breaks = c(seq(0,100,10)), limits = c(0,100))+
  #   geom_text_repel(mapping = aes(label=repel_ks), size=5, nudge_x = 0.01)+
  #   theme_bw()+gg_theme
  # 
  # A <- sub_test$wc
  # 
  # colnames(A)[2] <- 'postime'
  # A <- merge(A, g[, c(1,8:11)], by = 'L_R', all = T)
  # 
  # wc_sub <- ggplot(A, aes(x=fc, y=postime, color=label_wc, label=L_R))+
  #   geom_point(size=5)+
  #   geom_hline(yintercept=c(50), color='#E9573E', linetype='dashed', size=1.5)+
  #   geom_vline(xintercept=c(50), color='#E9573E', linetype='dashed', size=1.5)+
  #   labs(title = paste(cluster, 'neighbor cells LR pairs, perform statistic by subsampling')
  #        , subtitle = paste(nrow(A), ' LR pairs are examined (wc test)', sep = '')
  #        , x= 'Over/under expression direction', y= 'Siginificant times')+
  #   scale_y_continuous(labels = c(seq(0,100,10)), breaks = c(seq(0,100,10)), limits = c(0,100))+
  #   geom_text_repel(mapping = aes(label=repel_wc), size=5, nudge_x = 0.01)+
  #   theme_bw()+gg_theme
  
  # return(list(ks_w=ks_w, ks_sub=ks_sub, wc_w=wc_w, wc_sub=wc_sub))
  return(list(ks_w=ks_w, wc_w=wc_w))
}


A <- whole_sub_vol(whole_test = out_list_T_whole$`0_close_to_2_test`, sub_test = dc_xct$`0_close_to_2_test`
                   , cluster='Immune')
A <- whole_sub_vol(whole_test = out_list_T_whole$`2_close_to_0_test`, sub_test = dc_xct$`2_close_to_0_test`
                   , cluster='Invasive')

A <- whole_sub_vol(whole_test = hbc_xct$Immune_close_to_Invasive_test, sub_test = hbc_xct$Immune_close_to_Invasive_test
                   , cluster='Immune')

A <- whole_sub_vol(whole_test = hbc_xct$Invasive_close_to_Immune_test, sub_test = hbc_xct$Invasive_close_to_Immune_test
                   , cluster='Invasive')
# ks/wc convergence LR pairs ----------------------------------------------
A <- out_list_T_whole$`2_close_to_0_test`
A <- out_list_T_whole$`0_close_to_2_test`
A <- DC_xct$`0_close_to_2_test`
A <- DC_xct$`2_close_to_0_test`

A$logfc <- log(A$fc)
A$logks <- round(-log10(A$ks), digits = 3)
A$logwc <- round(-log10(A$wc), digits = 3)
A <- A %>% tibble::rownames_to_column('L_R')

out <- matrix(0,0,0) %>% as.data.frame()
x=0
for(i in c(seq(10, nrow(A)-5, 10), nrow(A)))  {
  sub_ks <- A %>% top_n(logks, n = i)
  sub_wc <- A %>% top_n(logwc, n = i)
  x=x+1
  out[x, 'shared'] <- intersect(sub_ks$L_R, sub_wc$L_R) %>% length()
  out[x, 'ks_only'] <- setdiff(sub_ks$L_R, sub_wc$L_R) %>% length()
  out[x, 'wc_only'] <- setdiff(sub_wc$L_R, sub_ks$L_R) %>% length()
}
out$x <- c(seq(10, nrow(A)-5, 10), nrow(A))

out$all <- out$shared+out$ks_only+out$wc_only

out <- data.frame(out[, c(1,2,3)]/out[, 'all']*100, x=out[,4]) %>% 
  tidyr::gather(key, value, c(1:3))

out$key <- factor(out$key, levels = unique(out$key) %>% rev())

ggplot(out)+
  geom_bar(aes(x=x, y=value, fill=key), color='black', stat = 'identity')+
  labs(title = 'ks test and wc test shared LR pairs'
              , subtitle = 'sender: invasive, receiver: immune' #'sender: immune, receiver: invasive'
              , x='Top n siganifacant L_R pairs'
              , y='LR_pairs')+
  scale_y_continuous(breaks = c(seq(0, 100, 10)))+
  scale_fill_manual(values = c('#4FC0E8', '#EC87BF', '#86DC7B'))+
  # labs(title = paste('sender: ', cluster_1, ', receiver: ', cluster_2, sep = ''))+
  theme_linedraw()+gg_theme

# out %>% tidyr::gather(key, value, c(1:3)) %>% 
#   ggplot()+
#   geom_line(aes(x=x, y=value, color=key), linewidth=2)+
#   geom_point(aes(x=x, y=value, color=key), size=5)+
#   labs(title = 'ks test and wc test shared LR pairs'
#        , subtitle = 'sender: immune, receiver: invasive' #'sender: invasive, receiver: immune' 
#        , x='Top n siganifacant L_R pairs'
#        , y='LR_pairs')+
#   scale_color_discrete(name=NULL)+
#   # labs(title = paste('sender: ', cluster_1, ', receiver: ', cluster_2, sep = ''))+
#   scale_x_continuous(breaks = c(seq(10, nrow(A)-5, 10), nrow(A)))+
#   theme_bw()+gg_theme




A <- DC_xct$`0_close_to_2_test`
A <- DC_xct$`2_close_to_0_test`
A <- hbc_xct$Immune_close_to_Invasive_test
A <- hbc_xct$Invasive_close_to_Immune_test
A$logfc <- log(A$fc)
A$logks <- -log10(A$ks)
A$logwc <- -log10(A$wc)
A <- A %>% tibble::rownames_to_column('L_R') %>% arrange(logwc)
A$index <- seq(1, nrow(A),1)

out <- data.frame(x=c(A$logks[order(A$logks)], A$logwc[order(A$logwc)]), y=c(seq(1, nrow(A),1), seq(1, nrow(A),1))
                  , color=c(rep('ks-test', nrow(A)), rep('wc-test', nrow(A))))

out <- data.frame(x=c(A$ks[order(A$ks, decreasing = F)], A$wc[order(A$wc, decreasing = F)]), y=c(seq(1, nrow(A),1), seq(1, nrow(A),1))
                  , color=c(rep('ks-test', nrow(A)), rep('wc-test', nrow(A))))


ggplot(out)+
  geom_line(aes(x=x, y=y, color=color), linewidth=0.5)+
  geom_point(aes(x=x, y=y, color=color), size=1.5)+
  labs(title = 'ks test and wc test siginificant level'
       , subtitle = 'sender: invasive, receiver: immune'  #'sender: immune, receiver: invasive'
       , x='P-value'
       , y='LR_pairs')+
  geom_vline(xintercept=c(0.05), color='#3CE8AC', linetype='dashed', size=1.5)+
  scale_color_manual(name=NULL, values = c('#4FC0E8', '#E9573E'))+
  scale_y_continuous(breaks = c(seq(0, nrow(A),25)))+
  # labs(title = paste('sender: ', cluster_1, ', receiver: ', cluster_2, sep = ''))+
  theme_linedraw()+gg_theme

