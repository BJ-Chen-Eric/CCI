feature_plot <- function(pathway_result, seu) {
  plot_list <- list()
  
  for(i in names(pathway_result)) {
    x <- c(2, 1); names(x) <- names(pathway_result[[i]][[3]])
    for(k in names(pathway_result[[i]][[3]]))  {
      p <- paste(str_extract(k, pattern = '.*_'), 'inside', '|'
                 , str_extract(k, pattern = '.*_'), 'PC_'
                 , names(pathway_result[[i]][[3]])[x[k]] %>% str_remove(pattern = '_.*'), sep = '')
      A <- seu %>% subset(cell=grep(pattern = p, seu$xct_ident))
      SpatialDimPlot(A, label = F, group.by = 'xct_ident')
      plot_tf <- plot <- list()
      gene <- pathway_result[[i]][[3]][[k]][, 1]
      gene <- gene[!gene %in% '']
      gene_table <- pathway_result[[i]][[3]][[k]]
      if(identical(gene, '')) {
        plot=''; next
      }else{
        for(j in seq_along(gene)) {
          p1 <- SpatialDimPlot(A, label = F, group.by = 'xct_ident')+guides(fill = guide_legend(override.aes = list(size=12)))+
            gg_theme
          p2 <- SpatialFeaturePlot(A, features = gene_table[j, 1])+
            labs(title = paste(gene_table[j, 2], ': ', gene_table[j, 1], sep = ''))+
            theme(legend.position="right")+gg_theme
          p3 <- VlnPlot(A, features = gene_table[j, 1], group.by = 'xct_ident')+
            labs(title = paste(gene_table[j, 2], ': ', gene_table[j, 1], sep = ''))+
            gg_theme
          if(gene_table[j, 2] == 'TF') {
            plot_tf[[gene_table[j, 1]]] <- ggarrange(plotlist = list(p1, p2, p3), ncol = 2, nrow = 2)
          }
          if(gene_table[j, 2] != 'TF') {
            plot[[gene_table[j, 1]]] <- ggarrange(plotlist = list(p1, p2, p3), ncol = 2, nrow = 2)
          }
        }
      }
      plot_list[[paste(k, 'TF', sep = '_')]] <- plot_tf
      plot_list[[paste(k, sep = '_')]] <- plot
    }
  }
  return(plot_list)
}




[1] "S100P"    "SERPINA3" "CALML5"   "HSP90AB1" "HSP90AA1" "TFF3"     "CFB"      "PBX1"     "BAMBI"    "DDR1"     "IGFBP5"  
[12] "ZNF263"   "SP1"      "MED30"    "TAL1"     "EBF1"     "ZNF674"  
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



cell_list <- list(Immune=c(1, 2), Benign=4)
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



# SpatialDimPlot(DC_subset, label = F, group.by = 'xct_ident')
# 
# DC$xct_ident <- DC_subset$xct_ident
# SpatialDimPlot(DC, label = F, group.by = 'xct_ident')





path <- fread('Downloads/Reactome_2022_table-3.txt')

path <- path[grep(pattern = 'THY1|COL5A2|SH3BGRL|LUM|COL3A1', path$Genes), ]



path1 <- fread('Downloads/enrich_b_only_L.txt') %>% filter(`Adjusted P-value` < 0.1)

path2 <- fread('Downloads/enrich_b_only_RD.txt') %>% filter(`Adjusted P-value` < 0.1)
path2 <- path2[grep(pattern = 'COX6C|SCGB2A2|RPL23A|CYB5A|CSTA|HEBP1|SNHG25|COL1A1|MUC1|DNAH5|COL6A2|H2AFJ', path2$Genes)]

path1 <- path1[path1$Term %in% intersect(path1$Term, path2$Term), ]
path2 <- path2[path2$Term %in% path1$Term, ]




path1 <- fread('Downloads/enrich_cis_only_L.txt') %>% filter(`Adjusted P-value` < 0.1)

path2 <- fread('Downloads/enrich_cis_only_RD.txt') %>% filter(`Adjusted P-value` < 0.1)
path2 <- path2[grep(pattern = 'SNHG25|RPL38|TFF3|BAMBI|TMEM150C|FOXA1|CD24|TOB1|MGP|RPL23|MUC19|RPL23A|PVALB', path2$Genes)]

path1 <- path1[path1$Term %in% intersect(path1$Term, path2$Term), ]
path2 <- path2[path2$Term %in% path1$Term, ]




path1 <- fread('Downloads/enrich_ic_only_L.txt') %>% filter(`Adjusted P-value` < 0.1)

path2 <- fread('Downloads/enrich_ic_only_RD.txt') %>% filter(`Adjusted P-value` < 0.1)
path2 <- path2[grep(pattern = 'THY1|COL5A2|SH3BGRL|LUM|COL3A1', path2$Genes)]

path1 <- path1[path1$Term %in% intersect(path1$Term, path2$Term), ]
path2 <- path2[path2$Term %in% path1$Term, ]




benign <- (fread('Desktop/project/Cai/ccc_result/1000ccc/Benign_close_to_Immune_r_deg_tf_l_label.txt', fill=TRUE) %>% filter(type=='R'))[,1] %>% unlist()
cis <- (fread('Desktop/project/Cai/ccc_result/1000ccc/CIS_close_to_Immune_r_deg_tf_l_label.txt', fill=TRUE) %>% filter(type=='R'))[,1] %>% unlist()
ic <- (fread('Desktop/project/Cai/ccc_result/1000ccc/IC_close_to_Immune_r_deg_tf_l_label.txt', fill=TRUE) %>% filter(type=='R'))[,1] %>% unlist()



benign <- (fread('Desktop/project/Cai/ccc_result/1000ccc/Benign_close_to_Immune_r_deg_tf_l_label.txt', fill=TRUE) %>% filter(type=='R'))[,3] %>% unlist()
cis <- (fread('Desktop/project/Cai/ccc_result/1000ccc/CIS_close_to_Immune_r_deg_tf_l_label.txt', fill=TRUE) %>% filter(type=='R'))[,3] %>% unlist()
ic <- (fread('Desktop/project/Cai/ccc_result/1000ccc/IC_close_to_Immune_r_deg_tf_l_label.txt', fill=TRUE) %>% filter(type=='R'))[,3] %>% unlist()



benign_imm <- (fread('Desktop/project/Cai/ccc_result/1000ccc/Immune_close_to_Benign_r_deg_tf_l_label.txt', fill=TRUE))
cis_imm <- (fread('Desktop/project/Cai/ccc_result/1000ccc/Immune_close_to_CIS_r_deg_tf_l_label.txt', fill=TRUE))
ic_imm <- (fread('Desktop/project/Cai/ccc_result/1000ccc/Immune_close_to_IC_r_deg_tf_l_label.txt', fill=TRUE))


benign <- (fread('Desktop/project/Cai/ccc_result/1000ccc/Benign_close_to_Immune_r_deg_tf_l_label.txt', fill=TRUE))
cis <- (fread('Desktop/project/Cai/ccc_result/1000ccc/CIS_close_to_Immune_r_deg_tf_l_label.txt', fill=TRUE))
ic <- (fread('Desktop/project/Cai/ccc_result/1000ccc/IC_close_to_Immune_r_deg_tf_l_label.txt', fill=TRUE))

write.table(setdiff(ic$Immune_ligand, union(benign$Immune_ligand, cis$Immune_ligand)), 'Downloads/ic_only_L.txt', quote = F, row.names = F, col.names = F)
C <- ic[ic$Immune_ligand %in% c(setdiff(ic$Immune_ligand, union(benign$Immune_ligand, cis$Immune_ligand)), ''), ]
write.table(C$IC_gene, 'Downloads/ic_only_RD.txt', quote = F, row.names = F, col.names = F)



write.table(setdiff(cis$Immune_ligand, union(benign$Immune_ligand, ic$Immune_ligand)), 'Downloads/cis_only_L.txt', quote = F, row.names = F, col.names = F)
C <- cis[cis$Immune_ligand %in% c(setdiff(cis$Immune_ligand, union(benign$Immune_ligand, ic$Immune_ligand)), ''), ]
write.table(C[, 1], 'Downloads/cis_only_RD.txt', quote = F, row.names = F, col.names = F)



write.table(setdiff(benign$Immune_ligand, union(cis$Immune_ligand, ic$Immune_ligand)), 'Downloads/b_only_L.txt', quote = F, row.names = F, col.names = F)
C <- benign[benign$Immune_ligand %in% c(setdiff(benign$Immune_ligand, union(cis$Immune_ligand, ic$Immune_ligand)), ''), ]
write.table(C[, 1], 'Downloads/b_only_RD.txt', quote = F, row.names = F, col.names = F)





cis <- (fread('Desktop/project/Cai/ccc_result/hbc/DCIS_close_to_Immune_r_deg_tf_l_label.txt', fill=TRUE))
ic <- (fread('Desktop/project/Cai/ccc_result/hbc/Invasive_close_to_Immune_r_deg_tf_l_label.txt', fill=TRUE))

write.table(setdiff(ic$Immune_ligand, cis$Immune_ligand), 'Downloads/hbc_ic_only_L.txt', quote = F, row.names = F, col.names = F)
C <- ic[ic$Immune_ligand %in% c(setdiff(ic$Immune_ligand, cis$Immune_ligand), ''), ]
write.table(C[, 1], 'Downloads/hbc_ic_only_RD.txt', quote = F, row.names = F, col.names = F)



write.table(setdiff(cis$Immune_ligand, ic$Immune_ligand), 'Downloads/hbc_cis_only_L.txt', quote = F, row.names = F, col.names = F)
C <- cis[cis$Immune_ligand %in% c(setdiff(cis$Immune_ligand, ic$Immune_ligand), ''), ]
write.table(C[, 1], 'Downloads/hbc_cis_only_RD.txt', quote = F, row.names = F, col.names = F)




# IC
path1 <- fread('Downloads/enrichr_ic_L_only.txt') %>% filter(`Adjusted P-value` < 0.1)

path2 <- fread('Downloads/enrichr_ic_RD_only.txt') %>% filter(`Adjusted P-value` < 0.1)
path2 <- path2[grep(pattern = 'CDH2|FABP4|G0S2|TACO1|HMGCS2|CFD|PVALB|MUC1|GPX3|PLIN4|ADIPOQ|MLPH|PLIN1|GSR|LPL|TFAP2C|COL6A3|SRPX'
                    , path2$Genes)]

path2 <- path2[grep(pattern = ic[ic$type %in% 'TFs', 1] %>% unlist() %>% paste(collapse = '|')
                    , path2$Genes)]

path1 <- path1[path1$Term %in% intersect(path1$Term, path2$Term), ]
path2 <- path2[path2$Term %in% path1$Term, ]





# CIS
path1 <- fread('Downloads/enrichr_cis_L_only.txt') %>% filter(`Adjusted P-value` < 0.1)

path2 <- fread('Downloads/enrichr_cis_RD_only.txt') %>% filter(`Adjusted P-value` < 0.1)
path2 <- path2[grep(pattern = 'PTGDS|BGN|PDGFRB'
                    , path2$Genes)]

path1 <- path1[path1$Term %in% intersect(path1$Term, path2$Term), ]
path2 <- path2[path2$Term %in% path1$Term, ]



## Reactome
# BiocManager::install('ReactomePA')
# BiocManager::install('reactome.db')
library(ReactomePA)
data(geneList, package="DOSE")
de <- names(geneList)[abs(geneList) > 1.5]
head(de)

# A <- str_split(c('RET;NTRK2;FLT3;ITGA3;CD80;FLT4;ITGA2;INSR;MYLK;IGF1R;ERBB3;CCRL2;MYB;KIT;ERBB2;KDR;ITGB8;COL6A3;CCR7;ITGB6;PLIN1;ACKR1;BMPR1B;FGFR2')
#                , pattern = ';') %>% unlist() %>% clusterProfiler::bitr(fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
A <- ic[,1] %>% unlist() %>% clusterProfiler::bitr(fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
x <- enrichPathway(gene=A[,2], pvalueCutoff = 0.05, readable=TRUE)
test <- x@result
test <- test[grep(pattern = 'THY1|COL5A2|SH3BGRL|LUM|COL3A1', x = test$geneID), ]

# EGF/AREG/EREG/FGFR4
# KIT/FGFR2/NRG1/FGFR3/ERBB3/ERBB2/PDGFRB

A <- ic[!ic$Immune_ligand %in% '',3] %>% unlist() %>% clusterProfiler::bitr(fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
x2 <- enrichPathway(gene=A[,2], pvalueCutoff = 0.05, readable=TRUE)
x2@result <- x2@result[x2@result$ID %in% test$ID, ]

test <- test[intersect(test$ID, x2@result$ID), ]
x@result <- test


barplot(x, showCategory=10)
dotplot(x, showCategory=15, x='p.adjust')
enrichplot::pairwise_termsim(x) %>% emapplot()
cnetplot(x, categorySize="pvalue", foldChange=geneList)


barplot(x2, showCategory=10)
dotplot(x2, showCategory = 10, x='p.adjust')
enrichplot::pairwise_termsim(x2) %>% emapplot()
cnetplot(x2, categorySize="pvalue", foldChange=geneList)

viewPathway("Extracellular matrix organization", readable=TRUE, foldChange=B)
viewPathway("E2F mediated regulation of DNA replication", readable=TRUE, foldChange=geneList)

cnetplot(x2, categorySize="pvalue", foldChange=B)

require(clusterProfiler)
data(gcSample)
res <- compareCluster(gcSample, fun="enrichPathway")
dotplot(res)


y <- gsePathway(A[,2], nPerm=10000,
                pvalueCutoff=0.2,
                pAdjustMethod="BH", verbose=FALSE)
res <- as.data.frame(y)
head(res)
emapplot(y, color="pvalue")
gseaplot(y, geneSetID = "R-HSA-69242")


test <- readRDS('Downloads/DC_IC/DC_IC_LR.RData')




B <- (c(ic[!ic$Immune_ligand %in% '',3], ic[,1]) %>% unlist() %>% clusterProfiler::bitr(fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db"))[,2]
x <- enrichPathway(gene=B, pvalueCutoff = 0.05, readable=TRUE)
test <- x@result
test <- test[grep(pattern = 'THY1|COL5A2|SH3BGRL|LUM|COL3A1', x = test$geneID), ]
x@result <- test
x@result <- x@result[x@result$ID %in% x2@result$ID, ]


barplot(x, showCategory=8)
dotplot(x, showCategory=15, x='p.adjust')
enrichplot::pairwise_termsim(x) %>% emapplot()
cnetplot(x, categorySize="pvalue", foldChange=geneList)




B <- (c(cis[!cis$Immune_ligand %in% '',3], cis[,1]) %>% unlist() %>%
        clusterProfiler::bitr(fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db"))[,2]
x2 <- enrichPathway(gene=B, pvalueCutoff = 0.05, readable=TRUE)
test <- x2@result
test <- test[grep(pattern = 'SNHG25|RPL38|TFF3|BAMBI|TMEM150C|FOXA1|CD24|TOB1|MGP|RPL23|MUC19|RPL23A|PVALB'
                  , x = test$geneID), ]
x2@result <- test


barplot(x2, showCategory=8)
dotplot(x, showCategory=15, x='p.adjust')
enrichplot::pairwise_termsim(x) %>% emapplot()
cnetplot(x, categorySize="pvalue", foldChange=geneList)




# random LR ---------------------------------------------------------------
lrdb <- fread('Downloads/data_baseLR.csv', header = T) %>% as.data.frame()
genes <- fread('Downloads/genes.txt', header = F) %>% as.data.frame()
row <- c()
for(i in c('^AC[0-9]+','^AF[0-9]+','^AL[0-9]+','^AP[0-9]+','^AC[0-9]+'
           ,'^BX[0-9]+','^CU[0-9]+','^FP[0-9]+','^LINC[0-9]+'
           ,'^MT-','^Z[0-9]+'))  { #, '^IG[HKL]'
  row <- c(row, grep(genes[, 1], pattern = i))
}
genes <- genes[-row, ] %>% as.data.frame()

out <- c()
for(i in 1:1000) {
  test <- lrdb[sample(1:nrow(lrdb), 50), ]
  test[51:65, 3] <- genes[sample(1:nrow(genes), 15), ] %>% unlist()
  
  A <- test[,3] %>% unlist() %>% clusterProfiler::bitr(fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  x <- enrichPathway(gene=A[,2], pvalueCutoff = 0.05, readable=TRUE)
  test <- x@result
  test <- test[grep(pattern = paste(test[51:65, 3], collapse = '|'), x = test$geneID), ]
  
  out[i] <- nrow(test)
}

# EGF/AREG/EREG/FGFR4
# KIT/FGFR2/NRG1/FGFR3/ERBB3/ERBB2/PDGFRB

A <- test[,2] %>% unlist() %>% clusterProfiler::bitr(fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
x2 <- enrichPathway(gene=A[,2], pvalueCutoff = 0.05, readable=TRUE)

x2@result <- x2@result[x2@result$ID %in% intersect(x2@result$ID, x@result$ID), ]

x@result <- x@result[x@result$ID %in% intersect(x2@result$ID, x@result$ID), ]


barplot(x2, showCategory=8)
barplot(x, showCategory=8)
dotplot(x, showCategory=15, x='p.adjust')
enrichplot::pairwise_termsim(x) %>% emapplot()
cnetplot(x, categorySize="pvalue", foldChange=geneList)





benign_imm <- (fread('Desktop/project/Cai/ccc_result/1000ccc/Immune_close_to_Benign_r_deg_tf_l_label.txt', fill=TRUE)) %>% as.data.frame()
benign_imm <- benign_imm[!benign_imm$type %in% 'TFs', ]
cis_imm <- (fread('Desktop/project/Cai/ccc_result/1000ccc/Immune_close_to_CIS_r_deg_tf_l_label.txt', fill=TRUE)) %>% as.data.frame()
cis_imm <- cis_imm[!cis_imm$type %in% 'TFs', ]
ic_imm <- (fread('Desktop/project/Cai/ccc_result/1000ccc/Immune_close_to_IC_r_deg_tf_l_label.txt', fill=TRUE)) %>% as.data.frame()
ic_imm <- ic_imm[!ic_imm$type %in% 'TFs', ]


benign <- (fread('Desktop/project/Cai/ccc_result/1000ccc/Benign_close_to_Immune_r_deg_tf_l_label.txt', fill=TRUE)) %>% as.data.frame()
benign <- benign[!benign$type %in% 'TFs', ]
cis <- (fread('Desktop/project/Cai/ccc_result/1000ccc/CIS_close_to_Immune_r_deg_tf_l_label.txt', fill=TRUE)) %>% as.data.frame()
cis <- cis[!cis$type %in% 'TFs', ]

# ic <- (fread('IC_close_to_Immune_r_deg_tf_l_label.txt', fill=TRUE)) %>% as.data.frame()
# ic <- ic[!ic$type %in% 'TFs', ]





geneLists <- list(geneListName=L)
# Motif enrichment analysis:
data(motifAnnotations_hgnc)
motifAnnotations <- motifAnnotations
motifRankings <- importRankings("~/Analysis/spatial/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather") # "~/Analysis/spatial/tf_df/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather"

motifEnrichmentTable_wGenes <- cisTarget(geneSets = geneLists
                                         , motifRankings = motifRankings
                                         , motifAnnot=motifAnnotations)

A <- motifEnrichmentTable_wGenes[!(motifEnrichmentTable_wGenes$TF_highConf %in% ''), ]
A <- A %>% top_n(NES, n = 8)

B <- str_remove(A$TF_highConf, pattern = ' \\(.*') %>% strsplit(B, split = '; ') %>% unlist() %>% unique()
write.table(c(L, B), 'L-TF.txt', quote = F, row.names = F,  col.names = F)


motifEnrichmentTable_wGenes_wLogo <- addLogo(A[1:8, ])
resultsSubset <- motifEnrichmentTable_wGenes_wLogo

library(DT)
datatable(resultsSubset[,-c("enrichedGenes", "TF_lowConf"), with=FALSE], 
          escape = FALSE, # To show the logo
          filter="top", options=list(pageLength=5))


A <- motifEnrichmentTable_wGenes[!(motifEnrichmentTable_wGenes$TF_highConf %in% ''), ]
A <- A %>% top_n(NES, n = 10)
B <- sub(A$TF_highConf, pattern = ' \\(.*', replacement = '')
tf <- strsplit(B, split = '; ') %>% unlist() %>% unique()


signifMotifNames <- motifEnrichmentTable$motif[1:5]
motifEnrichmentTable <- addMotifAnnotation(motifs_AUC, nesThreshold=3,
                                           motifAnnot=motifAnnotations_hgnc,
                                           highlightTFs=list(Bcell="RFX5",NaiveCD4T='ZNF274'))
incidenceMatrix <- getSignificantGenes(geneLists$geneListName, 
                                       motifRankings,
                                       signifRankingNames=signifMotifNames,
                                       plotCurve=TRUE, maxRank=5000, 
                                       genesFormat="incidMatrix",
                                       method="aprox")$incidMatrix



library(reshape2)
edges <- melt(incidenceMatrix)
edges <- edges[which(edges[,3]==1),1:2]
colnames(edges) <- c("from","to")


library(visNetwork)
motifs <- unique(as.character(edges[,1]))
genes <- unique(as.character(edges[,2]))
nodes <- data.frame(id=c(motifs, genes),   
                    label=c(motifs, genes),    
                    title=c(motifs, genes), # tooltip 
                    shape=c(rep("diamond", length(motifs)), rep("elypse", length(genes))),
                    color=c(rep("purple", length(motifs)), rep("skyblue", length(genes))))
visNetwork(nodes, edges) %>% visOptions(highlightNearest = TRUE, 
                                        nodesIdSelection = TRUE)




A <- cis[,1] %>% unlist() %>% clusterProfiler::bitr(fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
x <- enrichPathway(gene=A[,2], pvalueCutoff = 0.05, readable=TRUE)
test <- x@result
test <- test[grep(pattern = paste(cis[cis$type %in% 'DEGs', 1], collapse = '|')
                  , x = test$geneID), ]

# EGF/AREG/EREG/FGFR4
# KIT/FGFR2/NRG1/FGFR3/ERBB3/ERBB2/PDGFRB

A <- cis[!cis[, 3] %in% '',3] %>% unlist() %>% clusterProfiler::bitr(fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
x2 <- enrichPathway(gene=A[,2], pvalueCutoff = 0.05, readable=TRUE)
x2@result <- x2@result[x2@result$ID %in% test$ID, ]
test2 <- x2@result

test <- test[intersect(test$ID, x2@result$ID), ]
x@result <- test


barplot(x, showCategory=20, x="p.adjust")
enrichplot::dotplot(x, showCategory=20, x='p.adjust')
# ridgeplot(x, showCategory=10)
p1 <- enrichplot::pairwise_termsim(x) %>%
  emapplot(showCategory=20, cex_label_category=0.7)+
  ggtitle("Receptor/DRGs\n\n")
cnetplot(x, categorySize="pvalue", foldChange=geneList)


enrichplot::barplot(x2, showCategory=10)
enrichplot::dotplot(x2, showCategory = 10, x='p.adjust')
p2 <- enrichplot::pairwise_termsim(x2) %>%
  emapplot(showCategory=20, cex_label_category=0.7)+
  ggtitle("Ligand\n\n")
cnetplot(x2, categorySize="pvalue", foldChange=geneList)





A <- benign[,1] %>% unlist() %>% clusterProfiler::bitr(fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
x <- enrichPathway(gene=A[,2], pvalueCutoff = 0.05, readable=TRUE)
test <- x@result
test <- test[grep(pattern = paste(benign[benign$type %in% 'DEGs', 1], collapse = '|')
                  , x = test$geneID), ]

# EGF/AREG/EREG/FGFR4
# KIT/FGFR2/NRG1/FGFR3/ERBB3/ERBB2/PDGFRB

A <- benign[!benign[, 3] %in% '',3] %>% unlist() %>% clusterProfiler::bitr(fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
x2 <- enrichPathway(gene=A[,2], pvalueCutoff = 0.05, readable=TRUE)
x2@result <- x2@result[x2@result$ID %in% test$ID, ]
test2 <- x2@result

test <- test[intersect(test$ID, x2@result$ID), ]
x@result <- test


barplot(x, showCategory=20, x="p.adjust")
enrichplot::dotplot(x, showCategory=20, x='p.adjust')
# ridgeplot(x, showCategory=10)
p1 <- enrichplot::pairwise_termsim(x) %>%
  emapplot(showCategory=20, cex_label_category=0.7)+
  ggtitle("Receptor/DRGs\n\n")
cnetplot(x, categorySize="pvalue", foldChange=geneList)


enrichplot::barplot(x2, showCategory=10)
enrichplot::dotplot(x2, showCategory = 10, x='p.adjust')
p2 <- enrichplot::pairwise_termsim(x2) %>%
  emapplot(showCategory=20, cex_label_category=0.7)+
  ggtitle("Ligand\n\n")
cnetplot(x2, categorySize="pvalue", foldChange=geneList)



# Importing the list with the positional candiate genes
# remove.packages("STRINGdb")
# detach("package:STRINGdb", unload=TRUE)
# BiocManager::install("STRINGdb")

string_network <- function(xct_result, out_path=paste(getwd(), '/', sep = '')) {
  string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = 200, input_directory="", protocol="http")
  
  for(i in seq_along(xct_result)) {
    n <- names(xct_result)[i]
    
    test <- xct_result[[i]][['gene_TF']]
    for(j in seq_along(test)) {
      gene <- data.frame((test[[j]]), pvalue=0, logFC=0)
      example1_mapped <- string_db$map(gene, "gene", removeUnmappedRows = TRUE, )
      
      example1_mapped_sig <- string_db$add_diff_exp_color(subset(example1_mapped, pvalue >= 0, logFcColStr="logFC" ))
      example1_mapped_sig[example1_mapped_sig$gene %in% gene[gene$type %in% 'DEG', 'gene'], 'color'] <- '#D870AD'
      example1_mapped_sig[example1_mapped_sig$gene %in% gene[gene$type %in% c('R'), 'gene'], 'color'] <- '#6DBFE1'
      example1_mapped_sig[example1_mapped_sig$gene %in% gene[gene$type %in% 'TF', 'gene'], 'color'] <- '#93E191'
      example1_mapped_sig[example1_mapped_sig$gene %in% gene[gene$type %in% 'L', 'gene'], 'color'] <- '#AC92ED'
      
      payload_id <- string_db$post_payload(example1_mapped_sig$STRING_id,
                                           colors=example1_mapped_sig$color, )
      # display a STRING network png with the "halo"
      png(filename = paste(out_path, n, '_', names(test)[j], '_STRING.png', sep = ''), width = 1600, height = 900)
      string_db$plot_network(example1_mapped$STRING_id, payload_id=payload_id, network_flavor='confidence')
      dev.off()
    }
  }
}


hits <- example1_mapped$STRING_id

string_db$plot_network(hits, network_flavor='confidence')


example1_mapped_sig <- string_db$add_diff_exp_color(subset(example1_mapped, pvalue >= 0, logFcColStr="logFC" ))
head(example1_mapped_sig)


table(example1_mapped_sig$color)
dim(example1_mapped_sig)




hits <- example1_mapped$STRING_id
string_db$plot_network(hits, network_flavor='confidence')




A <- cisTarget(geneSets = geneLists
               , motifRankings = motifRankings
               , motifAnnot=motifAnnotations)

A <- A[!(A$TF_highConf %in% ''), ]
A <- A %>% top_n(NES, n = 8)

B <- str_remove(A$TF_highConf, pattern = ' \\(.*') %>% strsplit(B, split = '; ') %>% unlist() %>% unique()
write.table(c(L, B), 'L-TF.txt', quote = F, row.names = F,  col.names = F)


motifEnrichmentTable_wGenes_wLogo <- addLogo(A[1:8, ])
resultsSubset <- motifEnrichmentTable_wGenes_wLogo

library(DT)
datatable(resultsSubset[,-c("enrichedGenes", "TF_lowConf"), with=FALSE], 
          escape = FALSE, # To show the logo
          filter="top", options=list(pageLength=5))

motifs_AUC <- calcAUC(geneLists, motifRankings, nCores=1)
motifs_AUC

auc <- getAUC(motifs_AUC)[1, ]
hist(auc, main=1, xlab="AUC histogram",
     breaks=100, col="#ff000050", border="darkred")
abline(v=(3*sd(auc)) + mean(auc), col="red")


motifEnrichmentTable <- addMotifAnnotation(motifs_AUC, nesThreshold=3,
                                           motifAnnot=motifAnnotations)
head(motifEnrichmentTable[,-"TF_lowConf", with=FALSE])


motifEnrichmentTable_wGenes <- addSignificantGenes(motifEnrichmentTable,
                                                   rankings=motifRankings, 
                                                   geneSets=geneLists)

motifEnrichmentTable_wGenes[1:4,]


selectedMotifs <-sample(motifEnrichmentTable$motif, 3)
getSignificantGenes(geneLists[[1]], 
                    motifRankings,
                    signifRankingNames=selectedMotifs,
                    plotCurve=TRUE, maxRank=5000, genesFormat="none",
                    method="aprox")


motifEnrichmentTable_wGenes_wLogo <- addLogo(motifEnrichmentTable_wGenes)

resultsSubset <- motifEnrichmentTable_wGenes_wLogo[1:10,]

library(DT)
datatable(resultsSubset[,-c("enrichedGenes", "TF_lowConf"), with=FALSE], 
          escape = FALSE, # To show the logo
          filter="top", options=list(pageLength=5))



anotatedTfs  <- lapply(split(motifEnrichmentTable_wGenes$TF_highConf,
                             motifEnrichmentTable$geneSet),
                       function(x) {
                         genes <- gsub(" \\(.*\\). ", "; ", x, fixed=FALSE)
                         genesSplit <- unique(unlist(strsplit(genes, "; ")))
                         return(genesSplit)
                       })
anotatedTfs  


signifMotifNames <- motifEnrichmentTable$motif[1:30]

incidenceMatrix <- getSignificantGenes(geneLists$geneListName, 
                                       motifRankings,
                                       signifRankingNames=signifMotifNames,
                                       plotCurve=TRUE, maxRank=5000, 
                                       genesFormat="incidMatrix",
                                       method="aprox")$incidMatrix
incidenceMatrix


library(igraph)
mugh<- graph_from_incidence_matrix(incidenceMatrix, directed = F)
mugh

library(reshape2)
edges <- melt(incidenceMatrix)
edges <- edges[which(edges[,3]==1),1:2]
colnames(edges) <- c("from","to")


library(visNetwork)
motifs <- unique(as.character(edges[,1]))
genes <- unique(as.character(edges[,2]))
nodes <- data.frame(id=c(motifs, genes),   
                    label=c(motifs, genes),    
                    title=c(motifs, genes), # tooltip 
                    shape=c(rep("diamond", length(motifs)), rep("elypse", length(genes))),
                    color=c(rep("purple", length(motifs)), rep("skyblue", length(genes))))
visNetwork(nodes, edges) %>% visOptions(highlightNearest = TRUE, 
                                        nodesIdSelection = TRUE)

