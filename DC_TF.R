library(RcisTarget)
# sct <- readRDS('~/Analysis/spatial/DC_3000g_sct.RData')
geneLists <- list(geneListName=close_g)
geneLists <- list(geneListName=inner_g)


data(motifAnnotations_hgnc)
motifAnnotations <- motifAnnotations_hgnc
motifRankings <- importRankings("~/Analysis/spatial/tf_df/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather")
# Motif enrichment analysis:
motifEnrichmentTable_wGenes <- cisTarget(geneLists, motifRankings,
                                         motifAnnot=motifAnnotations)

A <- motifEnrichmentTable_wGenes[!(motifEnrichmentTable_wGenes$TF_highConf %in% ''), ]
A <- A %>% top_n(NES, n = 20)
B <- sub(A$TF_highConf, pattern = ' \\(.*', replacement = '')
B <- strsplit(B, split = '; ') %>% unlist()
write.table(c(B, close_g, unique(out[['0_close_to_2']]$Receptor))
            , file = '~/gene_immune.txt', quote = F, col.names = F, row.names = F)


A <- motifEnrichmentTable[!(motifEnrichmentTable$TF_highConf %in% ''), ]
signifMotifNames <- A$motif

incidenceMatrix <- getSignificantGenes(geneLists$geneListName, 
                                       motifRankings,
                                       signifRankingNames=signifMotifNames,
                                       plotCurve=TRUE, maxRank=5000, 
                                       genesFormat="incidMatrix",
                                       method="aprox")$incidMatrix
incidenceMatrix <- as.data.frame(incidenceMatrix)
rownames(incidenceMatrix) <- paste(A$TF_highConf, seq(1:100), sep = '_')

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


############
string_link <- fread('~/Downloads/9606.protein.links.v11.5.txt.gz') %>% as.data.frame()
string_link_sub <- string_link[string_link$combined_score >990, ]

string_gene <- fread('~/Downloads/9606.protein.aliases.v11.5.txt.gz') %>% as.data.frame()
string_gene <- string_gene[string_gene$source %in% 'Ensembl_WikiGene', ]

A <- string_gene[string_gene$alias %in% unique(out[['2_close_to_0']]$Receptor), 1]
A <- string_link_sub$protein2[string_link_sub$protein1 %in% A]

R_related <- string_gene[string_gene[,1] %in% A, 2]

intersect(R_related, motifAnnotations$TF)

A <- out[['2_close_to_0']]
tf <- intersect(R_related, motifAnnotations$TF)
A[(nrow(A)+1):(nrow(A)+length(tf)), 'Receptor'] <- tf
sct_p <- ggm_comparison(seu_raw = DC ,drug_path = '~/Analysis/spatial/interactions.tsv'
                        , tf_path = '~/Analysis/spatial/tf-target-infomation.txt'
                        , nHVG = 3000, L_R = A
                        , R_border = '2_close_to_0', R_inside = '2_inside'
                        , alpha = 1, method = 'SCT', drug_b=F, tf_b = F, specy = 'human')

sct_0 <- ggm_comparison(seu_raw = DC ,drug_path = '~/Analysis/spatial/interactions.tsv'
                        , tf_path = '~/Analysis/spatial/tf-target-infomation.txt'
                        , nHVG = 3000, L_R = out[['0_close_to_2']]
                        , R_border = '0_close_to_1', R_inside = '0_inside'
                        , alpha = 1, method = 'SCT', drug_b=F, tf_b = F, specy = 'human')

B <- sct[["SCT_GRN_B_I"]][["diffRegulation"]]
# B <- sct3000[["SCT_GRN_B_I"]][["diffRegulation"]]
B <- B[B$p.value <0.05, ]

C <- sct$SCT_GRN_control$diffRegulation
# C <- sct3000$SCT_GRN_control$diffRegulation
C <- C[C$p.value <0.05, ]

close_g <- B$gene[!(B$gene %in% intersect(B$gene, C$gene))]
write.table(c(tf, close_g), file = '~/gene.txt', quote = F, col.names = F, row.names = F)


