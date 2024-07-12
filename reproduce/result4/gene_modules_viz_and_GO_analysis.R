library(tidyverse)
library(ggsignif)
library(Seurat)

setwd('~/Desktop/Research/role_singlecell/results/')
set.seed(123)

my_colors <- c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
		"#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB",
		"#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
		"#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

glb_cluster_distance <- read.table('single_cell/glioblastoma/2_1_2_distance_info.csv', sep = '\t', header = TRUE)

plotEmbedding <- function(embedding){
	p1 <- ggplot(data = embedding, mapping = aes(x = nan_ratio, y = mean, color = as.factor(community), size = size)) + 
			 				 geom_point() + 
			 				 scale_color_manual(values = my_colors) + 
			 				 theme_bw() +
			 				 theme(panel.grid.major = element_blank(),
			 				 	   panel.grid.minor = element_blank(),
			 				 	   legend.position = "right"
			 				 	   ) +
			 				 labs(x = 'NA%', y = 'Mean distance',size = "No. of genes", color = "Gene module")
	return(p1)
}

p1 <- plotEmbedding(glb_cluster_distance)

ggsave('result3/glb_gene_modules.pdf',p1,dpi = 300, width = 6, height = 5)


library(clusterProfiler)
library('org.Hs.eg.db')
glb_index <- read.table('data/singleCell/glioblastoma/splitMatrix/index_tracker.tsv', sep = '\t', header = TRUE)
rownames(glb_index) <- glb_index$X0h

glb_lovain <- read.table('data/Figure4Result3Clusterinfo/Louvain_Exp2_1_2_0h_eeisp_7.csv', sep = '\t', header = TRUE)
glb_lovain_0_5 <- glb_lovain[glb_lovain$community %in% c(0,5),]
glb_lovain_0_5$gene_name <- glb_index[glb_lovain_0_5$gene,]$X

glb_lovain_0_5_list <- split(glb_lovain_0_5$gene_name, glb_lovain_0_5$community)

glb_lovain_0_5_BP <- compareCluster(glb_lovain_0_5_list,
                             fun=enrichGO,
                             OrgDb = 'org.Hs.eg.db',
                             keyType = 'SYMBOL',
                             ont = 'BP', 
                             pAdjustMethod = 'BH',
                             qvalueCutoff = 0.05)

p1 <- dotplot(glb_lovain_0_5_BP, showCategory=5) + 
    scale_y_discrete(labels=function(x) str_wrap(x, width=100)) + 
    theme(axis.title.x = element_blank())

ggsave('result3/glb_0_5_BP.pdf', p1, dpi = 200, units = 'cm', width = 15, height = 10)

#------
glb_cluster_distance12 <- read.table('single_cell/glioblastoma/2_1_2_distance_info_12.csv', sep = '\t', header = TRUE)

p1 <- plotEmbedding(glb_cluster_distance12)

ggsave('result3/glb_gene_modules12.pdf',p1,dpi = 300, width = 6, height = 5)

rownames(glb_index) <- as.character(glb_index$t12)
glb_lovain_12 <- read.table('data/Figure4Result3Clusterinfo/Louvain_Exp2_1_2_t12_eeisp.csv', sep = '\t', header = TRUE)
glb_lovain_12 <- glb_lovain_12[glb_lovain_12$community %in% unique(glb_cluster_distance12$community),]

glb_lovain_12_1_9 <- glb_lovain_12[glb_lovain_12$community %in% c(1,9),]
glb_lovain_12_1_9$gene_name <- glb_index[as.character(glb_lovain_12_1_9$gene),]$X

glb_lovain_12_1_9_list <- split(glb_lovain_12_1_9$gene_name, glb_lovain_12_1_9$community)

glb_lovain_1_9_BP <- compareCluster(glb_lovain_12_1_9_list,
                             fun=enrichGO,
                             OrgDb = 'org.Hs.eg.db',
                             keyType = 'SYMBOL',
                             ont = 'BP', 
                             pAdjustMethod = 'BH',
                             qvalueCutoff = 0.05)

p1 <- dotplot(glb_lovain_1_9_BP, showCategory=5) + 
    scale_y_discrete(labels=function(x) str_wrap(x, width=100)) + 
    theme(axis.title.x = element_blank())

ggsave('result3/glb_lovain_1_9_BP.pdf', p1, dpi = 200, units = 'cm', width = 15, height = 10)

###---
library('org.Mm.eg.db')

multiomics_0_3_distance_info <- read.table('result3/multiomics_0_3_distance_info.csv', sep = '\t', header = TRUE)
p1 <- plotEmbedding(multiomics_0_3_distance_info)
ggsave('result3/multiomics_0_3_modules.pdf',p1,dpi = 300, width = 6, height = 5)

multiomics_0_6_distance_info <- read.table('result3/multiomics_0_6_distance_info.csv', sep = '\t', header = TRUE)
p1 <- plotEmbedding(multiomics_0_6_distance_info)
ggsave('result3/multiomics_0_6_modules.pdf',p1,dpi = 300, width = 6, height = 5)

multiomics_0_9_distance_info <- read.table('result3/multiomics_0_9_distance_info.csv', sep = '\t', header = TRUE)
p1 <- plotEmbedding(multiomics_0_9_distance_info)
ggsave('result3/multiomics_0_9_modules.pdf',p1,dpi = 300, width = 6, height = 5)



multiomics_louvain <- read.table('data/Figure4Result3Clusterinfo/Louvain_Ery_0_network_12.csv', sep = '\t', header = TRUE)
multiomics_louvain_9_7 <- multiomics_louvain[multiomics_louvain$community %in% c(9,7),]

multiomics_index <- read.table('data/mult-omics/index_tracker.tsv', sep = '\t', header = TRUE)
multiomics_index <- multiomics_index[!is.na(multiomics_index$Ery_0_network),c('Gene', 'Ery_0_network')]
colnames(multiomics_index) <- c('gene_symbol', 'gene')

class(multiomics_index$gene) <- 'integer'
index_gene <- inner_join(multiomics_index, multiomics_louvain_9_7,by = 'gene')

index_gene_list <- split(index_gene$gene_symbol, index_gene$community)

index_gene_list_BP <- compareCluster(index_gene_list,
                             fun=enrichGO,
                             OrgDb = 'org.Mm.eg.db',
                             keyType = 'SYMBOL',
                             ont = 'BP', 
                             pAdjustMethod = 'BH',
                             qvalueCutoff = 0.05)

p1 <- dotplot(index_gene_list_BP, showCategory=5) + 
    scale_y_discrete(labels=function(x) str_wrap(x, width=40)) + 
    theme(axis.title.x = element_blank())

ggsave('result3/multiomics9_7_BP.pdf', p1, dpi = 200, units = 'cm', width = 15, height = 12)

###---
multiomics_0_7_gmp_distance_info <- read.table('result3/multiomics_index_mep_0_7_distance.csv', sep = '\t', header = TRUE)
p1 <- plotEmbedding(multiomics_0_7_gmp_distance_info)
ggsave('result3/multiomics_0_7_gmp_modules.pdf',p1,dpi = 300, width = 6, height = 5)

multiomics_louvain_gmp <- read.table('data/Figure4Result3Clusterinfo/Louvain_GMP_0_network.csv', sep = '\t', header = TRUE)
multiomics_louvain_gmp_2_6 <- multiomics_louvain_gmp[multiomics_louvain_gmp$community %in% c(2,6),]

multiomics_index <- read.table('data/mult-omics/index_tracker.tsv', sep = '\t', header = TRUE)
multiomics_index <- multiomics_index[!is.na(multiomics_index$GMP_0_network),c('Gene', 'GMP_0_network')]
colnames(multiomics_index) <- c('gene_symbol', 'gene')

class(multiomics_index$gene) <- 'integer'
index_gene <- inner_join(multiomics_index, multiomics_louvain_gmp_2_6,by = 'gene')

index_gene_list <- split(index_gene$gene_symbol, index_gene$community)

index_gene_list_BP <- compareCluster(index_gene_list,
                             fun=enrichGO,
                             OrgDb = 'org.Mm.eg.db',
                             keyType = 'SYMBOL',
                             ont = 'BP', 
                             pAdjustMethod = 'BH',
                             qvalueCutoff = 0.05)

p1 <- dotplot(index_gene_list_BP, showCategory=5) + 
    scale_y_discrete(labels=function(x) str_wrap(x, width=40)) + 
    theme(axis.title.x = element_blank())

ggsave('result3/multiomics_louvain_gmp_2_6_BP.pdf', p1, dpi = 200, units = 'cm', width = 15, height = 12)

