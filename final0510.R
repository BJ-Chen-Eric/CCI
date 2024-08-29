source('/home/eric/R/CCC/function_v3.R')
hbc <- Load10X_Spatial(
  data.dir='/home/eric/Analysis/spatial/10x/human_breast_cancer_X/CytAssist_FFPE_Human_Breast_Cancer_spatial/',
  filename = "CytAssist_FFPE_Human_Breast_Cancer_filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  filter.matrix = TRUE,
  to.upper = FALSE,
  image = NULL
)


SpatialDimPlot(hbc, label = TRUE, label.size = 3)


hbc <- NormalizeData(hbc, assay = "Spatial", verbose = FALSE) %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA(assay = "Spatial", verbose = FALSE) %>% 
  FindNeighbors(reduction = "pca", dims = 1:30) %>% 
  FindClusters(resolution = 0.2) %>% 
  RunUMAP(reduction = "pca", dims = 1:20)

# hbc <- NormalizeData(hbc, assay = "Spatial", verbose = FALSE) %>% 
# hbc <- FindVariableFeatures(hbc)
# hbc <- ScaleData(hbc)
# hbc <- RunPCA(hbc, assay = "Spatial", verbose = FALSE)
# hbc <- FindNeighbors(hbc, reduction = "pca", dims = 1:30)
# hbc <- FindClusters(hbc, resolution = 0.2)
# hbc <- RunUMAP(hbc, reduction = "pca", dims = 1:20)
# DimPlot(open_data, reduction = "umap", label = TRUE)

SpatialDimPlot(hbc)

# new.cluster.ids <- c("Stromal", "Adipocyte", "Invasive", "DCIS1_2", "Immune", "DCIS2")
# names(new.cluster.ids) <- levels(hbc)
# hbc <- RenameIdents(hbc, new.cluster.ids)
# SpatialDimPlot(hbc, label = TRUE)

# mgene <- fread('~/Analysis/spatial/10x/human_breast_cancer_X/markergene', header = F) %>% unlist()
# cluster_10x <- fread('~/Analysis/spatial/10x/human_breast_cancer_X/clusters.csv', header = T) %>% as.data.frame()
# cluster_10x <- cluster_10x[cluster_10x$Barcode %in% intersect(cluster_10x$Barcode, colnames(hbc)), ]
# rownames(cluster_10x) <- cluster_10x[,1]
# cluster_10x <- cluster_10x[colnames(hbc), ]
# 
# hbc$cluster_10x <- cluster_10x[,2]
# SpatialDimPlot(hbc, label = F, label.size = 3, group.by = 'cluster_10x')
# 
# SpatialFeaturePlot(hbc, features = 'SCGB2A2') # hbcIS1
# SpatialFeaturePlot(hbc, features = 'CPB1') # hbcIS2
# SpatialFeaturePlot(hbc, features = 'KRT17') # Myoepithelial
# SpatialFeaturePlot(hbc, features = 'FABP4') # Adipocytes
# SpatialFeaturePlot(hbc, features = 'IL2RG') # Immune
# SpatialFeaturePlot(hbc, features = 'SFRP2') # Stromal
# SpatialFeaturePlot(hbc, features = 'CDH2') # Invasive


# A <- subset(hbc, subset = IL2RG > boxplot.stats(hbc@assays$Spatial@scale.data['IL2RG', ])$stats[4], slot = 'scale.data')
A <- subset(hbc, subset = IL2RG > boxplot.stats(hbc@assays$Spatial@data['IL2RG', ])$stats[4]) #, slot = 'data'

B <- subset(hbc, subset = CDH2 > boxplot.stats(hbc@assays$Spatial@data['CDH2', ])$stats[4]) #, slot = 'data'
# B <- subset(hbc, subset = CEACAM6 > boxplot.stats(hbc@assays$Spatial@data['CEACAM6', ])$stats[4]) #, slot = 'data'


ann <- data.frame(cell= colnames(hbc), anno='undefinded')
ann[ann$cell %in% colnames(A), 'anno'] <- 'Immune'
ann[ann$cell %in% colnames(B), 'anno'] <- 'Invasive'
# ann[ann$cell %in% colnames(B), 'anno'] <- 'DCIS'
# ann[ann$cell %in% colnames(B), 'anno'] <- 'DCIS2'
rownames(ann) <- ann$cell
ann <- ann[colnames(hbc), ]
hbc$anno <- ann$anno
A <- hbc[, colnames(hbc) %in% c(colnames(A), colnames(B))]


png('~/data2_IC/cluster.png', res=600, width = 9600, height = 5400)
SpatialDimPlot(A, label = F, group.by = 'anno')+
  labs(fill='')+
  guides(fill = guide_legend(override.aes = list(size=12)))+
  gg_theme+theme(axis.title.x = element_blank(),axis.title.y = element_blank())
dev.off()
# DotPlot(hbc, features = c('IL2RG', 'CEACAM6'), group.by = 'anno')
png('~/data2_IC/dotplot.png', res=600, width = 9600, height = 5400)
DotPlot(hbc, features = c('IL2RG', 'CDH2'), group.by = 'anno', )+
  scale_size(range = c(3, 25)) +
  guides(fill = guide_legend(override.aes = list(size=12)))+
  gg_theme
dev.off()

png('~/data2_IC/feature.png', res=600, width = 9600, height = 5400)
(SpatialFeaturePlot(hbc, features = c('IL2RG'))+gg_theme+theme(axis.title.x = element_blank(),axis.title.y = element_blank()))+
(SpatialFeaturePlot(hbc, features = c('CDH2'))+gg_theme+theme(axis.title.x = element_blank(),axis.title.y = element_blank()))
dev.off()



hbc_subset <- prepare_data(seu = A, HVG = T, ori_ident = 'anno'
                           , cut_off = 4, nHVG = 3000)
png('~/data2_IC/cluster_neighbor.png', res=600, width = 9600, height = 5400)
png('~/data2_IC/cluster_neighbor.png', width = 1600, height = 900)
SpatialDimPlot(hbc_subset, label = F, group.by = 'xct_ident')+gg_theme+
  labs(fill='')+
  guides(fill = guide_legend(override.aes = list(size=12)))+
  gg_theme+theme(axis.title.x = element_blank(),axis.title.y = element_blank())
dev.off()

hbc$xct_ident <- hbc_subset$xct_ident
hbc <- prepare_data(seu = hbc, HVG = F, detectborder = F)

# sender: immune, receptor: invasive
### build tensor xct 
hbc_xct <- spatialXct(seu = hbc, specy = 'human', nHVGs = 2000
                      , cluster_1 = 'Immune', cluster_2 = 'Invasive', Tensor_dec = T, sampling = F)
saveRDS(hbc_xct, '~/Analysis/spatial/hbc_xct_v3.RData')
# hbc_xct <- A[1:4]
hbc_xct <- spatialXct(seu = hbc, specy = 'human', nHVGs = 2000
                      , cluster_1 = 'Immune', cluster_2 = 'DCIS1', Tensor_dec = T, sampling = F)
saveRDS(hbc_xct, '~/Analysis/spatial/hbc_xct_v3_DCIS1.RData')


hbc_xct <- spatialXct(seu = hbc, specy = 'human', nHVGs = 1000
                      , cluster_1 = 'Immune', cluster_2 = 'DCIS', Tensor_dec = T, sampling = F)
saveRDS(hbc_xct, '~/Analysis/spatial/hbc_xct_v3_DCIS.RData')

gc()

### identify degs
hbc_degs <- diff_genes(seu = hbc, drug_path = '~/Analysis/spatial/interactions.tsv'
                       , tf_path = '~/Analysis/spatial/tf-target-infomation.txt'
                       , nHVG = 1000, spatialxct = hbc_xct
                       , cluster_1 = 'Immune', cluster_2 = 'Invasive')
saveRDS(hbc_degs, '~/Analysis/spatial/hbc_deg_invasive.RData')
hbc_degs <- readRDS('~/Analysis/spatial/hbc_1000grn.RData')


hbc_degs <- diff_genes(seu = hbc, drug_path = '~/Analysis/spatial/interactions.tsv'
                       , tf_path = '~/Analysis/spatial/tf-target-infomation.txt'
                       , nHVG = 1000, spatialxct = hbc_xct
                       , cluster_1 = 'Immune', cluster_2 = 'DCIS1')
saveRDS(hbc_degs, '~/Analysis/spatial/hbc_deg_DCIS1.RData')


hbc_degs <- diff_genes(seu = hbc, drug_path = '~/Analysis/spatial/interactions.tsv'
                       , tf_path = '~/Analysis/spatial/tf-target-infomation.txt'
                       , nHVG = 1000, spatialxct = hbc_xct
                       , cluster_1 = 'Immune', cluster_2 = 'DCIS')
saveRDS(hbc_degs, '~/Analysis/spatial/hbc_deg_DCIS.RData')
# hbc_degs <- A[5:6]
### identify TFs by degs
hbc_xct <- readRDS('~/Analysis/spatial/hbc_xct_v3.RData')
hbc_degs <- readRDS('~/Analysis/spatial/hbc_deg_invasive.RData')

# invasive_nei_tf <- neighbor_tf(neighbor_deg = hbc_degs$Invasive_close_to_Immune$neighbor_g
#                                , inner_deg = hbc_degs$Invasive_close_to_Immune$inner_g
#                                , Rcis_feather_dir = "~/Analysis/spatial/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather"
#                                , nTfs = 8)
# 
# 
# immune_nei_tf <- neighbor_tf(neighbor_deg = hbc_degs$Immune_close_to_Invasive$neighbor_g
#                              , inner_deg = hbc_degs$Immune_close_to_Invasive$inner_g
#                              , Rcis_feather_dir = "~/Analysis/spatial/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather"
#                              , nTfs = 8)


hbc_ic_out <- write_output(xct = hbc_xct, deg = hbc_degs, n_receptor = 50, n_ligand = 50, write_path = '~/d2_IC_test/')


# meta_genes_vol <- feature_vio_vol(seu = hbc, clu1 = 'Immune', clu2 = 'Invasive'
#                                   , clu1_genes_lable = meta_genes[[2]], clu2_genes_lable = meta_genes[[1]]
#                                   , out_path = '~/Analysis/spatial/plot/0820/')

hbc_IC <- hbc_subset %>% subset(subset=xct_ident %in% c('Immune_CI', 'Immune_PC_Invasive', 'Invasive_CI', 'Invasive_PC_Immune'))
d2_IC <- infer_pathway(seu_raw = hbc_IC, xct_raw = hbc_ic_out, out_path = '~/d2_IC/', write = T)
string_network(xct_result = d2_IC, '~/d2_IC/')

### DCIS
hbc_xct <- readRDS('~/Analysis/spatial/hbc_xct_v3_DCIS.RData')
hbc_degs <- readRDS('~/Analysis/spatial/hbc_deg_DCIS.RData')
invasive_nei_tf <- neighbor_tf(neighbor_deg = hbc_degs$DCIS_PC$neighbor_g
                               , inner_deg = hbc_degs$DCIS_PC$inner_g
                               , Rcis_feather_dir = "~/Analysis/spatial/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather"
                               , nTfs = 8)


immune_nei_tf <- neighbor_tf(neighbor_deg = hbc_degs$Immune_PC$neighbor_g
                             , inner_deg = hbc_degs$Immune_PC$inner_g
                             , Rcis_feather_dir = "~/Analysis/spatial/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather"
                             , nTfs = 8)


hbc_meta_genes <- write_output(xct = hbc_xct, n_receptor = 50, n_ligand = 50, c1 = 'Immune', c2 = 'DCIS'
                               , deg = hbc_degs, c1_nei_tf = immune_nei_tf, c2_nei_tf = invasive_nei_tf)


in_gene <- rowSums(c2)
in_gene <- in_gene[in_gene != 0]
in_gene <- names(in_gene)[names(in_gene) %in% meta$gene]
location <- c2$xct_ident %>% as.data.frame() %>% tibble::rownames_to_column('cell') 
colnames(location)[2] <- 'xct'
location <- group_split(location, xct) %>% as.list() %>% lapply(function(x){x[,1] %>% unlist()})
border <- hbc@assays$Spatial@data[in_gene, location[[1]]] %>% as.matrix()
inside <- hbc@assays$Spatial@data[in_gene, location[[2]]] %>% as.matrix()

output <- scTenifoldNet(X = border, Y = inside, qc = F,
                        nc_nNet = 10, nc_nCells = 500,
                        td_K = 3, qc_minLibSize = 30,
                        nc_symmetric = F, td_nDecimal = 10)

grn <- output$tensorNetworks$X %>% as.matrix()
grn[1:nrow(grn), ] <- ifelse(abs(grn[1:nrow(grn), ])<0.1, yes = 0, no = grn[1:nrow(grn), ] )
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




hbc_imm_to_tumor_degs <- ggm_comparison(seu = hbc ,drug_path = '~/Analysis/spatial/interactions.tsv'
                                        , tf_path = '~/Analysis/spatial/tf-target-infomation.txt'
                                        , nHVG = 3000, L_R = nei_r
                                        , R_border = 'Invasive_close_to_Immune'
                                        , R_inside = 'Invasive_inside'
                                        , alpha = 1, method = 'SCT', drug_b=T, tf_b = T, specy = 'human')




# out_tensor <- meta_test(meta_LR_list = hbc_imm_to_tumor)



test <- makeL_R_list_meta(seu = hbc_subset
                          , sender_i = 'Immune_inside', sender = 'Immune_close_to_Invasive'
                          , receiver_i = 'Invasive_inside', receiver = 'Invasive_close_to_Immune'
                          , Tensor_dec = T)






# genes_10x_intersect -----------------------------------------------------
path <- '/home/eric/Analysis/spatial/plot/old_10x/'

gene_upset <- function(dir, file)  {
  gene_list <- list()
  x <- 0
  for(i in c(100, 300, 1000))  {
    g_p <- paste(dir, 'deg', i, '/', file, sep = '')
    x=x+1
    gene_list[[x]] <- fread(g_p, header = F) %>% as.data.frame() %>% filter(!V2 %in% c('TFs', 'R', 'L'))
  }
  names(gene_list) <- paste('deg', c(100, 300, 1000), sep = '')
  
  genes <- lapply(gene_list, function(x){x[,1]}) %>% unlist() %>% unique()
  genes <- data.frame(genes=genes)
  
  x=2
  for(i in names(gene_list))  {
    A <- gene_list[[i]]
    genes[genes$genes %in% A[, 1], i] <- 1
    x=x+1
  }
  
  genes <- replace(genes, is.na(genes), 0)
  
  library(ComplexUpset)
  upset(genes ,colnames(genes)[2:4]
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
    labs(title = 'DVG recombination among conditions'
         , subtitle = 'Recombination more than 1 ')+ 
    theme(title = element_text(size=36))
}

# invasive
gene_upset(dir = '/home/eric/Analysis/spatial/plot/old_10x/'
           , file = '2_close_to_0r_deg_tf_l_label.txt') %>% print()

# immune
gene_upset(dir = '/home/eric/Analysis/spatial/plot/old_10x/'
           , file = '0_close_to_2_r_deg_tf_l_label.txt') %>% print()


### open 10x
# invasive
gene_upset(dir = '/home/eric/Analysis/spatial/plot/hbc/'
           , file = 'Invasive_close_to_Immuner_deg_tf_l_label.txt') %>% print()

# immune
gene_upset(dir = '/home/eric/Analysis/spatial/plot/hbc/'
           , file = 'Immune_close_to_Invasive_r_deg_tf_l_label.txt') %>% print()





### DCIS 1
hbc <- Load10X_Spatial(
  data.dir='/home/eric/Analysis/spatial/10x/human_breast_cancer_X/CytAssist_FFPE_Human_Breast_Cancer_spatial/',
  filename = "CytAssist_FFPE_Human_Breast_Cancer_filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  filter.matrix = TRUE,
  to.upper = FALSE,
  image = NULL
)


SpatialDimPlot(hbc, label = TRUE, label.size = 3)

hbc <- NormalizeData(hbc, assay = "Spatial", verbose = FALSE)
hbc <- FindVariableFeatures(hbc)
hbc <- ScaleData(hbc)
hbc <- RunPCA(hbc, assay = "Spatial", verbose = FALSE)
hbc <- FindNeighbors(hbc, reduction = "pca", dims = 1:30)
hbc <- FindClusters(hbc, resolution = 0.2)
hbc <- RunUMAP(hbc, reduction = "pca", dims = 1:20)
# DimPlot(open_data, reduction = "umap", label = TRUE)

SpatialDimPlot(hbc)

# new.cluster.ids <- c("Stromal", "Adipocyte", "Invasive", "DCIS1_2", "Immune", "DCIS2")
# names(new.cluster.ids) <- levels(hbc)
# hbc <- RenameIdents(hbc, new.cluster.ids)
# SpatialDimPlot(hbc, label = TRUE)

# mgene <- fread('~/Analysis/spatial/10x/human_breast_cancer_X/markergene', header = F) %>% unlist()
# cluster_10x <- fread('~/Analysis/spatial/10x/human_breast_cancer_X/clusters.csv', header = T) %>% as.data.frame()
# cluster_10x <- cluster_10x[cluster_10x$Barcode %in% intersect(cluster_10x$Barcode, colnames(hbc)), ]
# rownames(cluster_10x) <- cluster_10x[,1]
# cluster_10x <- cluster_10x[colnames(hbc), ]
# 
# hbc$cluster_10x <- cluster_10x[,2]
# SpatialDimPlot(hbc, label = F, label.size = 3, group.by = 'cluster_10x')
# 
# SpatialFeaturePlot(hbc, features = 'SCGB2A2') # hbcIS1
# SpatialFeaturePlot(hbc, features = 'CPB1') # hbcIS2
# SpatialFeaturePlot(hbc, features = 'KRT17') # Myoepithelial
# SpatialFeaturePlot(hbc, features = 'FABP4') # Adipocytes
# SpatialFeaturePlot(hbc, features = 'IL2RG') # Immune
# SpatialFeaturePlot(hbc, features = 'SFRP2') # Stromal
# SpatialFeaturePlot(hbc, features = 'CDH2') # Invasive


# A <- subset(hbc, subset = IL2RG > boxplot.stats(hbc@assays$Spatial@scale.data['IL2RG', ])$stats[4], slot = 'scale.data')
A <- subset(hbc, subset = IL2RG > boxplot.stats(hbc@assays$Spatial@data['IL2RG', ])$stats[4]) #, slot = 'data'
# B <- subset(hbc, subset = CDH2 > boxplot.stats(hbc@assays$Spatial@scale.data['CDH2', ])$stats[4])
B <- subset(hbc, subset = SCGB2A2 > boxplot.stats(hbc@assays$Spatial@data['SCGB2A2', ])$stats[4]) #, slot = 'data' 

ann <- data.frame(cell= colnames(hbc), anno='unc')
ann[ann$cell %in% colnames(A), 'anno'] <- 'Immune'
ann[ann$cell %in% colnames(B), 'anno'] <- 'DCIS1'
rownames(ann) <- ann$cell
ann <- ann[colnames(hbc), ]

hbc$anno <- ann$anno
A <- hbc[, colnames(hbc) %in% c(colnames(A), colnames(B))]
SpatialDimPlot(A, label = F, group.by = 'anno')
DotPlot(hbc, features = c('IL2RG', 'SCGB2A2'), group.by = 'anno')
SpatialFeaturePlot(hbc, features = c('IL2RG', 'SCGB2A2'))


hbc_subset <- prepare_data(seu = A, HVG = T, ori_ident = 'anno'
                           , cut_off = 4, nHVG = 3000)
SpatialDimPlot(hbc_subset, label = F, group.by = 'xct_ident')+gg_theme
hbc$xct_ident <- hbc_subset$xct_ident
hbc <- prepare_data(seu = hbc, HVG = F, detectborder = F)

# sender: immune, receptor: invasive
### build tensor xct 
hbc_xct <- spatialXct(seu = hbc, specy = 'human', nHVGs = 2000
                      , cluster_1 = 'Immune', cluster_2 = 'DCIS1', Tensor_dec = T, sampling = F)
saveRDS(hbc_xct, '~/Analysis/spatial/hbc_xct_dcis.RData')
# hbc_xct <- A[1:4]
gc()

hbc_xct <- readRDS('~/Analysis/spatial/hbc_xct_v2.RData')
### identify degs
hbc_degs <- diff_genes(seu = hbc, drug_path = '~/Analysis/spatial/interactions.tsv'
                       , tf_path = '~/Analysis/spatial/tf-target-infomation.txt'
                       , nHVG = 1000, spatialxct = hbc_xct
                       , cluster_1 = 'Immune', cluster_2 = 'Invasive')
saveRDS(hbc_degs, '~/Analysis/spatial/hbc_deg_invasive.RData')
hbc_degs <- readRDS('~/Analysis/spatial/hbc_1000grn.RData')

# hbc_degs <- A[5:6]
### identify TFs by degs
invasive_nei_tf <- neighbor_tf(neighbor_deg = hbc_degs$Invasive_close_to_Immune$neighbor_g
                               , inner_deg = hbc_degs$Invasive_close_to_Immune$inner_g
                               , Rcis_feather_dir = "~/Analysis/spatial/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather"
                               , nTfs = 8)


immune_nei_tf <- neighbor_tf(neighbor_deg = hbc_degs$Immune_close_to_Invasive$neighbor_g
                             , inner_deg = hbc_degs$Immune_close_to_Invasive$inner_g
                             , Rcis_feather_dir = "~/Analysis/spatial/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather"
                             , nTfs = 8)


hbc_meta_genes <- write_output(xct = hbc_xct, n_receptor = 50, n_ligand = 50, c1 = 'Immune', c2 = 'Invasive'
                               , deg = hbc_degs, c1_nei_tf = immune_nei_tf, c2_nei_tf = invasive_nei_tf)



meta_genes_vol <- feature_vio_vol(seu = hbc, clu1 = 'Immune', clu2 = 'Invasive'
                                  , clu1_genes_lable = meta_genes[[2]], clu2_genes_lable = meta_genes[[1]]
                                  , out_path = '~/Analysis/spatial/plot/0820/')
