# method benchmarking -----------------------------------------------------

neighbor_wei <- makeL_R_list(seu = DC_subset, times = 100, sender = '0_close_to_2', ncell = 100, max_dis = 30
                             , receiver = '2_close_to_0', weight = T, adjacent = T)
inner_wei <- makeL_R_list(seu = DC_subset, times = 100, sender = '0_inside', ncell = 100, max_dis = 30
                          , receiver = '2_inside', weight = T, adjacent = T)


sample_cell_raw_ks <- runtest_sample_cell(neighbor_list = neighbor_wei, inner_list = inner_wei
                                          , test = 'ks', times = 100)
sample_cell_raw_wc <- runtest_sample_cell(neighbor_list = neighbor_wei, inner_list = inner_wei
                                          , test = 'wc', times = 100)

sample_cell_tensor_ks <- runtest_sample_cell_tensor(neighbor_list = neighbor_wei, inner_list = inner_wei
                                                    , test = 'ks', times = 100)
sample_cell_tensor_wc <- runtest_sample_cell_tensor(neighbor_list = neighbor_wei, inner_list = inner_wei
                                                    , test = 'wc', times = 100)

## unwieghted 
neighbor <- makeL_R_list(seu = DC_subset, times = 100, sender = '0_close_to_2', ncell = 100, max_dis = 30
                         , receiver = '2_close_to_0', weight = F, adjacent = F)
inner <- makeL_R_list(seu = DC_subset, times = 100, sender = '0_inside', ncell = 100, max_dis = 30
                      , receiver = '2_inside', weight = F, adjacent = F)

uw_sample_cell_raw_ks <- runtest_sample_cell(neighbor_list = neighbor, inner_list = inner, test = 'ks', times = 100)
uw_sample_cell_raw_wc <- runtest_sample_cell(neighbor_list = neighbor, inner_list = inner, test = 'wc', times = 100)

#error
uw_sample_cell_tensor_ks <- runtest_sample_cell_tensor(neighbor_list = neighbor, inner_list = inner
                                                       , test = 'ks', times = 100)
uw_sample_cell_tensor_wc <- runtest_sample_cell_tensor(neighbor_list = neighbor, inner_list = inner
                                                       , test = 'wc', times = 100)


meta_LR_tensor <- makeL_R_list_meta(seu = DC_subset, sender_i = '0_inside', sender = '0_close_to_2'
                                    , receiver_i = '2_inside', receiver = '2_close_to_0', Tensor_dec = T)
meta_LR <- makeL_R_list_meta(seu = DC_subset, sender_i = '0_inside', sender = '0_close_to_2'
                             , receiver_i = '2_inside', receiver = '2_close_to_0',Tensor_dec =  )


##### compare methods
draw_hist(data1 = sample_cell_raw_ks, data2 = sample_cell_raw_wc, l1 = 'ks', l2 ='wc'
          , title = 'cell sampling, distance weighted, without tensor\nks and wc test comparison')

draw_hist(data1 = sample_cell_tensor_raw_ks, data2 = sample_cell_tensor_wc, l1 = 'ks', l2 ='wc'
          , title = 'cell sampling, distance weighted, with tensor\nks and wc test comparison')

draw_hist(data1 = sample_cell_raw_ks, data2 = sample_cell_tensor_raw_ks, l1='no tensor', l2 = 'Tensor'
          , title='cell sampling, distance weighted, ks test/nTensor and without tensor comparison')

draw_hist(data1 = sample_cell_raw_wc, data2 = sample_cell_tensor_wc, l1='no tensor', l2 = 'Tensor'
          , title='cell sampling, distance weighted, wc test/nTensor and without tensor comparison')

draw_hist(data1 = uw_sample_cell_tensor_ks, data2 = sample_cell_tensor_raw_ks
          , l1='no weight', l2 = 'weight'
          , title='cell sampling, ks test, with tensor\nWeighted and without nWeighted comparison')

draw_hist(data1 = uw_sample_cell_tensor_wc, data2 = sample_cell_tensor_wc
          , l1='no weight', l2 = 'weight'
          , title='cell sampling, wc test, with tensor\nWeighted and without nWeighted comparison')
draw_hist(data1 = uw_sample_cell_tensor_ks, data2 = uw_sample_cell_tensor_wc
          , l1='ks', l2 = 'wc'
          , title='unweighted, cell sampling, with tensor\nks and wc comparison')


# upset -------------------------------------------------------------------
L_R <- list()
for(i in seq_along(methods)[1:4])  {
  L_R[[i]] <- (keep_sgi_LR(methods[[i]]) %>% unlist() %>% table() %>% 
             as.data.frame(stringsAsFactors=F) %>% filter(Freq >90))[, 1]
}
names(L_R) <- names(methods)[1:4]
L_R[[2]] <- (keep_sgi_LR(methods[[2]]) %>% unlist() %>% table() %>% 
               as.data.frame(stringsAsFactors=F) %>% filter(Freq >50))[, 1]
L_R[[5]] <- paste(all_weighted_cell_ks[[1]]$Ligand, '_', all_weighted_cell_ks[[1]]$Receptor, sep = '')
names(L_R)[5] <- 'All_cell'
L_R_upset <- unlist(L_R) %>% unique()
out <- matrix(data = 0, nrow = length(L_R_upset), ncol = 8) %>% as.data.frame()
out$V1 <- L_R_upset
x <- 1
for(i in seq_along(L_R))  {
  x=x+1
  A <- L_R[[i]]
  out[out$V1 %in% A, x] <- 1
}
colnames(out)[2:6] <- names(L_R)

library(ComplexUpset)
upset(out ,colnames(out)[2:6]
      , width_ratio = 0.4 #ratio between horizen and vertical bar
      , name='Conditions' # names of verticle bar
      , base_annotations=list('Intersection size'=
                                intersection_size(text=list(vjust=-0.1,hjust=0.5, size=5))) # text in group bar
      , set_sizes=(upset_set_size()+ theme(axis.text.x=element_text(size=16)
                                           , axis.title = element_text(size=20)))
      , themes=upset_modify_themes(
        list('Intersection size'=theme(axis.text=element_text(size=15)
                                       , axis.title=element_text(size=30))
             , 'overall_sizes'=theme(axis.text.x=element_text(size = 24))
             , 'intersections_matrix'=theme(text=element_text(size=24))))
      ,  wrap=TRUE
      , matrix=intersection_matrix(geom=geom_point(stroke=4))
      , stripes=upset_stripes(geom = geom_segment(size=15)))+
  labs(title = 'Significant L_R union among weighted methon'
       , subtitle = 'sample_cell_raw_ks Significance threshold: > 50% iteration\nother methods Significance threshold: > 90% iteration')+ 
  theme(title = element_text(size=36))




L_R <- list()
x=0
for(i in seq_along(methods)[5:8])  {
  x=x+1
  L_R[[x]] <- (keep_sgi_LR(methods[[i]]) %>% unlist() %>% table() %>% 
                 as.data.frame(stringsAsFactors=F) %>% filter(Freq >90))[, 1]
}
names(L_R) <- names(methods)[5:8]

L_R[[5]] <- paste(all_weighted_cell_ks[[1]]$Ligand, '_', all_weighted_cell_ks[[1]]$Receptor, sep = '')
names(L_R)[5] <- 'All_cell'
L_R_upset <- unlist(L_R) %>% unique()
out <- matrix(data = 0, nrow = length(L_R_upset), ncol = 8) %>% as.data.frame()
out$V1 <- L_R_upset
x <- 1
for(i in seq_along(L_R))  {
  x=x+1
  A <- L_R[[i]]
  out[out$V1 %in% A, x] <- 1
}
colnames(out)[2:6] <- names(L_R)

library(ComplexUpset)
upset(out ,colnames(out)[2:6]
      , width_ratio = 0.4 #ratio between horizen and vertical bar
      , name='Conditions' # names of verticle bar
      , base_annotations=list('Intersection size'=
                                intersection_size(text=list(vjust=-0.1,hjust=0.5, size=5))) # text in group bar
      , set_sizes=(upset_set_size()+ theme(axis.text.x=element_text(size=16)
                                           , axis.title = element_text(size=20)))
      , themes=upset_modify_themes(
        list('Intersection size'=theme(axis.text=element_text(size=15)
                                       , axis.title=element_text(size=30))
             , 'overall_sizes'=theme(axis.text.x=element_text(size = 24))
             , 'intersections_matrix'=theme(text=element_text(size=24))))
      ,  wrap=TRUE
      , matrix=intersection_matrix(geom=geom_point(stroke=4))
      , stripes=upset_stripes(geom = geom_segment(size=15)))+
  labs(title = 'Significant L_R union among unweighted methon'
       , subtitle = 'Significance threshold: > 90% iteration')+ 
  theme(title = element_text(size=36))


## all method
L_R <- list()
for(i in seq_along(methods))  {
  L_R[[i]] <- (keep_sgi_LR(methods[[i]]) %>% unlist() %>% table() %>% 
                 as.data.frame(stringsAsFactors=F) %>% filter(Freq >90))[, 1]
}
names(L_R) <- names(methods)
L_R[[2]] <- (keep_sgi_LR(methods[[2]]) %>% unlist() %>% table() %>% 
               as.data.frame(stringsAsFactors=F) %>% filter(Freq >50))[, 1]
L_R[[9]] <- paste(all_weighted_cell_ks[[1]]$Ligand, '_', all_weighted_cell_ks[[1]]$Receptor, sep = '')
names(L_R)[9] <- 'All_cell'
L_R_upset <- unlist(L_R) %>% unique()
out <- matrix(data = 0, nrow = length(L_R_upset), ncol = 8) %>% as.data.frame()
out$V1 <- L_R_upset
x <- 1
for(i in seq_along(L_R))  {
  x=x+1
  A <- L_R[[i]]
  out[out$V1 %in% A, x] <- 1
}
colnames(out)[2:10] <- names(L_R)
out <- replace(out,is.na(out),0)

library(ComplexUpset)
upset(out ,colnames(out)[2:10]
      , width_ratio = 0.4 #ratio between horizen and vertical bar
      , name='Conditions' # names of verticle bar
      , base_annotations=list('Intersection size'=
                                intersection_size(text=list(vjust=-0.1,hjust=0.5, size=5))) # text in group bar
      , set_sizes=(upset_set_size()+ theme(axis.text.x=element_text(size=16)
                                           , axis.title = element_text(size=20)))
      , themes=upset_modify_themes(
        list('Intersection size'=theme(axis.text=element_text(size=15)
                                       , axis.title=element_text(size=30))
             , 'overall_sizes'=theme(axis.text.x=element_text(size = 24))
             , 'intersections_matrix'=theme(text=element_text(size=24))))
      ,  wrap=TRUE
      , matrix=intersection_matrix(geom=geom_point(stroke=4))
      , stripes=upset_stripes(geom = geom_segment(size=15)))+
  labs(title = 'Significant L_R union among weighted methon'
       , subtitle = 'sample_cell_raw_ks Significance threshold: > 50% iteration\nother methods Significance threshold: > 90% iteration')+ 
  theme(title = element_text(size=36))

# stacked bar -------------------------------------------------------------
out_bar <- matrix() %>% as.data.frame()
for(i in seq_along(methods))  {
  A <- keep_sgi_LR(methods[[i]]) %>% unlist() %>% table() %>% 
         as.data.frame(stringsAsFactors=F) 
  out_bar[i,1] <- A %>% filter(Freq>90) %>% nrow() 
  out_bar[i,2] <- nrow(A)-out_bar[i,1]
  
}
rownames(out_bar) <- sub(names(methods), pattern = 'cell_', replacement = 'cell\n')
out_bar[,3] <- out_bar[,1]/out_bar[i,2]*100
out_bar[,4] <- 100-out_bar[,3]
out_bar <- tibble::rownames_to_column(out_bar, 'method')
out_bar$test <- sub(out_bar$method, pattern = '.*\\n[a-z]+_', replacement = '')
out_bar$wei <- c(rep('Weighted', 4), rep('Unweighted', 4))
out_bar$wei <- factor(out_bar$wei, levels = c('Weighted', 'Unweighted'))
colnames(out_bar)[c(2,3)] <- c('Sig (> 90% iteration)', 'Un_sig')

out_bar[, c(1,2,3)] %>% tidyr::gather('Significance', 'value', colnames(out_bar)[2:3]) %>% 
  ggplot(aes(x=method, y = value, fill=Significance))+
  geom_bar(stat="identity")+
  labs(title = 'L-R pairs among methods', subtitle = 'Significance threshold: > 90% iteration'
       , y='Frequence')+
  scale_fill_manual(values = c('#EC5D62', '#5D9CEC'))+
  theme_bw()+gg_theme#+theme(axis.text.x = element_blank())

ggplot(out_bar[, c(1,4,5,7)], aes(x=wei, y = V3, fill=method))+
  geom_bar(stat = 'identity', position = 'dodge')+
  geom_text(aes(label=method), position=position_dodge(width=0.9), vjust=-0.1)+
  labs(title = 'Significance proportion L-R pairs among methods'
       , subtitle = 'Significance threshold: > 90% iteration'
       , y='Significance proportion (%)')+
  scale_fill_brewer(palette = 'Paired')+
  theme_bw()+gg_theme#+theme(axis.text.x = element_blank())



# tensor 100 ks -----------------------------------------------------------
# meta_LR_tensor <- makeL_R_list_meta(seu = DC_subset, sender_i = '0_inside', sender = '0_close_to_2'
#                                     , receiver_i = '2_inside', receiver = '2_close_to_0', Tensor_dec = T)
meta_LR_tensor <- readRDS('~/R/CCC/meta_tensor.RData')
out_tensor <- meta_test(meta_LR_list = meta_LR_tensor)


out <- cbind(out_tensor[[1]][, 1:2], out_tensor[[2]][, 2:3])  
colnames(out) <- c('L_R', 'ks', 'wc', 'fc')
out <- replace(out,is.na(out),0)

# nei_r <- sub(out[out$ks >50, 1], pattern = '.*_', replacement = '') %>% unique()
nei_r <- sub(out[out$ks >50 & out$fc >50 , 1], pattern = '.*_', replacement = '') %>% unique()

out %>% tidyr::gather('test', 'value', colnames(out)[2:3], na.rm = T) %>% 
  ggplot()+
  geom_histogram(aes(x=value, fill=test), color="#e9ecef", alpha=0.5, binwidth = 1, position = 'identity')+
  labs(title = 'All cell, unweighted, with tensor, sampling denoised product, ks and wc test comparison'
       , subtitle = paste('ks test', ' identified total LR pair=', table(0==out[,2])[['FALSE']], ', > 50% iteration = ', table(out[,2]>50)[['TRUE']], ' pairs', '\n',
                          'wc test', ' identified total LR pair=', table(0==out[,3])[['FALSE']], ', > 50% iteration = ', table(out[,3]>50)[['TRUE']], ' pairs', sep = ''))+
  scale_fill_manual(values=c("#404080", "#69b3a2"))+
  theme_bw()+gg_theme



sct_3000 <- readRDS('~/R/CCC/DC_3000g_sct.RData')
B <- sct_3000[["SCT_GRN_B_I"]][["diffRegulation"]]
# B <- sct3000[["SCT_GRN_B_I"]][["diffRegulation"]]
B <- B[B$p.value <0.05, ]

C <- sct_3000$SCT_GRN_control$diffRegulation
# C <- sct3000$SCT_GRN_control$diffRegulation
C <- C[C$p.value <0.05, ]


inner_g <- C$gene[!(C$gene %in% intersect(B$gene, C$gene))]
close_g <- B$gene[!(B$gene %in% intersect(B$gene, C$gene))]


A <- neighbor_tf(neighbor_deg = close_g, inner_deg = inner_g, nTfs = 10
            , Rcis_feather_dir = "~/Analysis/spatial/tf_df/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather")


write.table(c(close_g, nei_r, close_tf), file = '~/tumor_neighbor_new.txt', quote = F, col.names = F, row.names = F)

meta <- data.frame(gene=c(close_g, nei_r, close_tf), type=c(rep('dg', length(close_g))
                                                     , rep('r', length(nei_r))
                                                     , rep('tf', length(close_tf))))

A <- fread('~/Analysis/spatial/string_color', sep = '\t', header = F) %>% as.data.frame()
string <- fread('~/Analysis/spatial/color_table.tsv', sep = '\t', header = F) %>% as.data.frame() 

for(i in seq_len(nrow(string)))  {
  out <- meta[meta$gene %in% string[i,1], 2]
  
  if(out=='dg')  {
    string[i,2] <- A[1,3]
  }
  if(out=='tf')  {
    string[i,2] <- A[2,3]
  }
  if(out=='r')  {
    string[i,2] <- A[3,3]
  }
}
write.table(string, file = '~/Analysis/spatial/color_table.tsv', sep = '\t', col.names = F, row.names = F, quote = F)



A <- subset(DC, subset = xct_ident %in% c('Immune_close_to_Invasive'))
B <- subset(DC, subset = xct_ident %in% c('Immune_inside'))
C <- subset(DC, subset = xct_ident %in% c('Invasive_close_to_Immune'))
D <- subset(DC, subset = xct_ident %in% c('Invasive_inside'))

c2 <- scCustomize::Merge_Seurat_List(list(A,B,C))
c2 <- DC[, colnames(DC) %in% colnames(c2)]
p3 <- SpatialDimPlot(c2, group.by = 'xct_ident')


plot_list <- list()
x=0
in_gene <- rownames(c2)[rownames(c2) %in% c(nei_r, close_tf, close_g)]
location <- c2$xct_ident %>% as.data.frame() %>% tibble::rownames_to_column('cell') 
colnames(location)[2] <- 'xct'
location <- group_split(location, xct) %>% as.list()
clu <- c(location[[1]][1,2], location[[2]][1,2], location[[3]][1,2]) %>% unlist()
location <- location %>% lapply(function(x){x[,1] %>% unlist()})
names(location) <- clu
for(i in in_gene)  {
  x=x+1
  nei <- c2@assays$SCT@data[i, location[['2_close_to_0']]]
  inn <- c2@assays$SCT@data[i, location[['0_close_to_2']]] # 0_close_to_2
  p <- round(wilcox.test(nei, inn, alternative = 'greater')$p.value, digits = 10)
  p2 <- SpatialFeaturePlot(c2, features = i)
  p4 <- VlnPlot(c2, features = i, group.by = 'xct_ident')+labs(subtitle = paste('wilcoxon test: ', p, sep = ''))
  plot_list[[x]] <- p3+p2+p3+p4
  names(plot_list)[x] <- i
}


dir.create(file.path('~/Analysis/spatial/plot/0410/'), showWarnings = FALSE)
for(i in seq_along(plot_list))  {
  png(paste('~/Analysis/spatial/plot/0410/'
            , names(plot_list)[i], '.png', sep = ''), width = 1600, height =  900)
  print(plot_list[[i]])
  dev.off()
}





out <- cbind(out_tensor[[1]][, 1:2], out_tensor[[2]][, 2:3])  
colnames(out) <- c('L_R', 'ks', 'wc', 'fc')
nei_r <- sub(out[out$ks >40 & out$fc >50 , 1], pattern = '.*_', replacement = '') %>% unique()

in_gene <- rowSums(c2)
in_gene <- in_gene[in_gene != 0]
in_gene <- names(in_gene)[names(in_gene) %in% c(nei_r, close_tf, close_g)]
location <- c2$xct_ident %>% as.data.frame() %>% tibble::rownames_to_column('cell') 
colnames(location)[2] <- 'xct'
location <- group_split(location, xct) %>% as.list() %>% lapply(function(x){x[,1] %>% unlist()})
border <- DC@assays$SCT@data[in_gene, location[[1]]] %>% as.matrix()
inside <- DC@assays$SCT@data[in_gene, location[[2]]] %>% as.matrix()
print(paste('B_I GRN gene number=', nrow(inside), sep=''))
output <- scTenifoldNet(X = border, Y = inside, qc = F,
                        nc_nNet = 10, nc_nCells = 500,
                        td_K = 3, qc_minLibSize = 30,
                        nc_symmetric = F, td_nDecimal = 10)

grn <- output$tensorNetworks$X %>% as.matrix()
grn[1:nrow(grn), ] <- ifelse(grn[1:nrow(grn), ]<0.05, yes = 0, no = grn[1:nrow(grn), ] )
gene <- rownames(grn)[colSums(abs(grn)) !=0]
grn <- grn[gene, gene]
tr <- igraph::graph_from_adjacency_matrix(grn, weighted = TRUE, diag = F)
E(tr)$weight <- (1+abs(E(tr)$weight)*50/12)
E(tr)$weight <- abs(E(tr)$weight)
plot(tr, vertex.size=10, vertex.label=NA) 
e <- get.edgelist(tr,names=FALSE)
l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(tr),
                                       area=8*(vcount(tr)^5),repulse.rad=(vcount(tr)^3))

curves <- autocurve.edges2(tr)
set.seed(713)

plot(tr, layout=l, edge.arrow.size=0.3, vertex.size=12, vertex.frame.width=1
     , vertex.color='#F9D182', vertex.frame.color='#F6BB43', label.color=''
     , edge.arrow.mode=3
     , edge.color=ifelse(border <= 0, no = "blue", yes = "red")
     , edge.width=E(tr)$weight
     , edge.curved=seq(-0.7, 0.7, length = ecount(tr))
     , cex.main=20)



 
