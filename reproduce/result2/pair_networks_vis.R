library(tidyverse)
library(ggsignif)
library(Seurat)

setwd('~/Desktop/Research/role_singlecell/results/')
set.seed(123)




################################################################################################

glb_distance <- read.table('single_cell/glioblastoma/glb_eeisp_distance.csv', sep = '\t', header = TRUE)

threshold <- quantile(glb_distance$distance, 0.90, na.rm = TRUE)
glb_distance <- glb_distance[!is.na(glb_distance$distance),]
glb_distance$color <- ifelse(glb_distance$distance > threshold, "DTGs", "non-DTGs")

p1 <- ggplot(glb_distance,aes(x=distance, fill = color)) + 
		    geom_histogram(binwidth=0.05, alpha=0.9)+
		    scale_fill_manual(values=c("DTGs"="#FF0000", "non-DTGs"="#69b3a2")) +
		    theme_bw() +
		    theme(panel.grid.major = element_blank(),
		          panel.grid.minor = element_blank(),
		          legend.position = c(0.9, 0.9),
		          legend.background = element_rect(fill = NA)) +
		    labs(x = 'Distance', y = 'Count') +
		    guides(fill = guide_legend(title = NULL))

ggsave('single_cell/glioblastoma/hist.pdf',p1,dpi = 300, width = 6, height = 5)

glb_dtgs <- glb_distance[glb_distance$color == 'DTGs',]
colnames(glb_dtgs)[1] <- 'symbol'

non_dtgs <- glb_distance[glb_distance$color == 'non-DTGs',]

#
glb_degs <- read.table('../data/singleCell/glioblastoma/merge_seurat.markers.txt', 
											 sep = ' ', 
											 header = TRUE)

symbol_ensembl <- read.table('../../metadata/formatted_genes.txt')
colnames(symbol_ensembl) <- c('ensembl', 'symbol')
rownames(symbol_ensembl) <- symbol_ensembl$ensembl
glb_degs$symbol <- symbol_ensembl[glb_degs$gene,]$symbol


library(ggVennDiagram)
gene_lists <- list(DTGs = glb_dtgs$symbol,
									 DEGs = na.omit(glb_degs$symbol))


##### Comparing DTGs with DEGs

# Seurat Object that have 
merge_seurat <- readRDS('../data/singleCell/glioblastoma/merge_seurat.rds')
merge_metadata <- merge_seurat@meta.data

merge_metadata$CD164 <- merge_seurat@assays$RNA$data['ENSG00000135535',]
merge_metadata$DKK3 <- merge_seurat@assays$RNA$data['ENSG00000050165',]
merge_metadata$PBK <- merge_seurat@assays$RNA$data['ENSG00000168078',]
merge_metadata$MCM4 <- merge_seurat@assays$RNA$data['ENSG00000104738',]

p1 <- ggplot(merge_metadata, aes(x = orig.ident, y = DKK3)) + 
     			geom_violin(aes(fill = orig.ident)) +
     			geom_jitter(width = 0.15, size = 0.1) + 
    			geom_signif(comparisons = list(c("0h", "t12")), 
                			test = "wilcox.test") + 
    			labs(x = '', y = 'Normalized expression level') + 
    			theme_bw() + 
    			theme(panel.grid.major = element_blank(),
		          panel.grid.minor = element_blank()) + NoLegend()


ggsave('single_cell/glioblastoma/DKK3.pdf',p1,dpi = 300, width = 6, height = 3)

p1 <- ggplot(merge_metadata, aes(x = orig.ident, y = CD164)) + 
     			geom_violin(aes(fill = orig.ident)) +
     			geom_jitter(width = 0.15, size = 0.1) + 
    			geom_signif(comparisons = list(c("0h", "t12")), 
                			test = "wilcox.test") + 
    			labs(x = '', y = 'Normalized expression level') + 
    			theme_bw() + 
    			theme(panel.grid.major = element_blank(),
		          panel.grid.minor = element_blank()) + NoLegend()


ggsave('single_cell/glioblastoma/CD164.pdf',p1,dpi = 300, width = 6, height = 3)

p1 <- ggplot(merge_metadata, aes(x = orig.ident, y = PBK)) + 
     			geom_violin(aes(fill = orig.ident)) +
     			geom_jitter(width = 0.15, size = 0.1) + 
    			geom_signif(comparisons = list(c("0h", "t12")), 
                			test = "wilcox.test") + 
    			labs(x = '', y = 'Normalized expression level') + 
    			theme_bw() + 
    			theme(panel.grid.major = element_blank(),
		          panel.grid.minor = element_blank()) + NoLegend()


ggsave('single_cell/glioblastoma/PBK.pdf',p1,dpi = 300, width = 6, height = 3)

p1 <- ggplot(merge_metadata, aes(x = orig.ident, y = MCM4)) + 
     			geom_violin(aes(fill = orig.ident)) +
     			geom_jitter(width = 0.15, size = 0.1) + 
    			geom_signif(comparisons = list(c("0h", "t12")), 
                			test = "wilcox.test") + 
    			labs(x = '', y = 'Normalized expression level') + 
    			theme_bw() + 
    			theme(panel.grid.major = element_blank(),
		          panel.grid.minor = element_blank()) + NoLegend()


ggsave('single_cell/glioblastoma/MCM4.pdf',p1,dpi = 300, width = 6, height = 3)


################################################################################################

pbmc_distance <- read.table('result2/pbmc_NM_CD4_distance.csv', sep = '\t', header = TRUE)

threshold <- quantile(pbmc_distance$distance, 0.90, na.rm = TRUE)
pbmc_distance <- pbmc_distance[!is.na(pbmc_distance$distance),]
pbmc_distance$color <- ifelse(pbmc_distance$distance > threshold, "DTGs", "non-DTGs")

p1 <- ggplot(pbmc_distance,aes(x=distance, fill = color)) + 
		    geom_histogram(binwidth=0.05, alpha=0.9)+
		    scale_fill_manual(values=c("DTGs"="#FF0000", "non-DTGs"="#69b3a2")) +
		    theme_bw() +
		    theme(panel.grid.major = element_blank(),
		          panel.grid.minor = element_blank(),
		          legend.position = c(0.9, 0.9),
		          legend.background = element_rect(fill = NA)) +
		    labs(x = 'Distance', y = 'Count') +
		    guides(fill = guide_legend(title = NULL))

ggsave('result2/pbmc_pairs/hist.pdf',p1,dpi = 300, width = 6, height = 5)

DTGs <- pbmc_distance[pbmc_distance$color == 'DTGs',]$Unnamed..0


##### Comparing DTGs with DEGs

pbmc <- readRDS('data/singleCell/pbmc_rna/pbmc_processed.rds')

Idents(pbmc) <- pbmc$celltype
NM_CD4_markers <- FindMarkers(pbmc,ident.1 = 'CD4.N', ident.2 = '22_CD4.M')
NM_CD4_markers$gene_name <- rownames(NM_CD4_markers)
NM_CD4_markers <- NM_CD4_markers[NM_CD4_markers$p_val_adj < 0.05,]

#head(intersect(DTGs,unique(rownames(NM_CD4_markers))))
#[1] "CARD16" "SAT1"   "GPRIN3" "UTS2"   "MLEC"   "PHLDA1"

head(setdiff(DTGs,unique(rownames(NM_CD4_markers))))
# [1] "PSTPIP2" "RRM1"    "KLF9"    "REC8"    "PCK2"    "GTPBP3" 

non_dtgs <- pbmc_distance[pbmc_distance$color == 'non-DTGs',]
non_dtgs <- non_dtgs[non_dtgs$Unnamed..0 %in% NM_CD4_markers$gene_name,]

head(non_dtgs[order(non_dtgs$distance),1])
#[1] "KDM5A"  "SP140"  "MT2A"   "MTHFD2" "TXK"    "CD8A"  

merge_metadata <- pbmc@meta.data
merge_metadata$CARD16 <- pbmc@assays$RNA$data['CARD16',]
merge_metadata$KLF6 <- pbmc@assays$RNA$data['KLF6',]
merge_metadata$SP140 <- pbmc@assays$RNA$data['SP140',]
merge_metadata$PSTPIP2 <- pbmc@assays$RNA$data['PSTPIP2',]

merge_metadata <- merge_metadata[merge_metadata$celltype %in% c('CD4.N', '22_CD4.M'),]
merge_metadata$celltype <- factor(merge_metadata$celltype, level = c('CD4.N', '22_CD4.M'))

p1 <- ggplot(merge_metadata, aes(x = celltype, y = KLF6)) + 
     			geom_violin(aes(fill = celltype)) +
     			geom_jitter(width = 0.15, size = 0.1) + 
    			geom_signif(comparisons = list(c("CD4.N", "22_CD4.M")), 
                			test = "wilcox.test") + 
    			labs(x = '', y = 'Normalized expression level') + 
    			theme_bw() + 
    			theme(panel.grid.major = element_blank(),
		          panel.grid.minor = element_blank()) + NoLegend()
ggsave('result2/pbmc_pairs/KLF6.pdf',p1,dpi = 300, width = 6, height = 3)

p1 <- ggplot(merge_metadata, aes(x = celltype, y = CARD16)) + 
     			geom_violin(aes(fill = celltype)) +
     			geom_jitter(width = 0.15, size = 0.1) + 
    			geom_signif(comparisons = list(c("CD4.N", "22_CD4.M")), 
                			test = "wilcox.test") + 
    			labs(x = '', y = 'Normalized expression level') + 
    			theme_bw() + 
    			theme(panel.grid.major = element_blank(),
		          panel.grid.minor = element_blank()) + NoLegend()
ggsave('result2/pbmc_pairs/CARD16.pdf',p1,dpi = 300, width = 6, height = 3)

p1 <- ggplot(merge_metadata, aes(x = celltype, y = SP140)) + 
     			geom_violin(aes(fill = celltype)) +
     			geom_jitter(width = 0.15, size = 0.1) + 
    			geom_signif(comparisons = list(c("CD4.N", "22_CD4.M")), 
                			test = "wilcox.test") + 
    			labs(x = '', y = 'Normalized expression level') + 
    			theme_bw() + 
    			theme(panel.grid.major = element_blank(),
		          panel.grid.minor = element_blank()) + NoLegend()
ggsave('result2/pbmc_pairs/SP140.pdf',p1,dpi = 300, width = 6, height = 3)

p1 <- ggplot(merge_metadata, aes(x = celltype, y = PSTPIP2)) + 
     			geom_violin(aes(fill = celltype)) +
     			geom_jitter(width = 0.15, size = 0.1) + 
    			geom_signif(comparisons = list(c("CD4.N", "22_CD4.M")), 
                			test = "wilcox.test") + 
    			labs(x = '', y = 'Normalized expression level') + 
    			theme_bw() + 
    			theme(panel.grid.major = element_blank(),
		          panel.grid.minor = element_blank()) + NoLegend()
ggsave('result2/pbmc_pairs/PSTPIP2.pdf',p1,dpi = 300, width = 6, height = 3)



################################################################################################

bmmc_distance <- read.table('result2/BMMC_distance.csv', sep = '\t', header = TRUE)

threshold <- quantile(bmmc_distance$distance, 0.90, na.rm = TRUE)
bmmc_distance <- bmmc_distance[!is.na(bmmc_distance$distance),]
bmmc_distance$color <- ifelse(bmmc_distance$distance > threshold, "DTGs", "non-DTGs")

p1 <- ggplot(bmmc_distance,aes(x=distance, fill = color)) + 
		    geom_histogram(binwidth=0.05, alpha=0.9)+
		    scale_fill_manual(values=c("DTGs"="#FF0000", "non-DTGs"="#69b3a2")) +
		    theme_bw() +
		    theme(panel.grid.major = element_blank(),
		          panel.grid.minor = element_blank(),
		          legend.position = c(0.9, 0.9),
		          legend.background = element_rect(fill = NA)) +
		    labs(x = 'Distance', y = 'Count') +
		    guides(fill = guide_legend(title = NULL))

ggsave('result2/bmmc_pairs/hist.pdf',p1,dpi = 300, width = 6, height = 5)

#bmmc  <- readRDS('data/singleCell/bmmc_rna/bmmc_rna.rds')

# bmmc <-  bmmc %>% 
#     SCTransform() %>%
#     RunPCA() %>%
#     FindNeighbors(dims = 1:20) %>%
#     RunUMAP(dims = 1:20)

#saveRDS(bmmc,'data/singleCell/bmmc_rna/bmmc_processed.rds')	

##### Comparing DTGs with DEGs

bmmc <- readRDS('data/singleCell/bmmc_rna/bmmc_processed.rds')
Idents(bmmc) <- bmmc$celltype

bmmc_markers <- FindMarkers(bmmc,ident.1 = 'GMP', ident.2 = 'CD14.Mono')
bmmc_markers <- bmmc_markers[bmmc_markers$p_val_adj < 0.05,]
bmmc_markers$gene_name <- rownames(bmmc_markers)
DTGs <- bmmc_distance[bmmc_distance$color == 'DTGs',]$Unnamed..0


head(intersect(DTGs,unique(rownames(bmmc_markers))))
#[1] "ETS1"  "CENPM" "CKAP2" "CCT7"  "SYNE2" "PA2G4"

head(setdiff(DTGs,unique(rownames(bmmc_markers))))
# [1] "PSMD2"    "HLA-DQA2" "TRAT1"    "TC2N"     "TMEM70"  
# [6] "CD6"  

non_dtgs <- bmmc_distance[bmmc_distance$color == 'non-DTGs',]
non_dtgs <- non_dtgs[non_dtgs$Unnamed..0 %in% rownames(bmmc_markers),]

head(non_dtgs[order(non_dtgs$distance),1])
# [1] "RETN"   "S100A8" "OXCT1"  "VDAC3"  "FANCL"  "SMTN"  

merge_metadata <- bmmc@meta.data
merge_metadata$ETS1 <- bmmc@assays$RNA$data['ETS1',]
merge_metadata$RASGEF1B <- bmmc@assays$RNA$data['RASGEF1B',]
merge_metadata$RETN <- bmmc@assays$RNA$data['RETN',]

merge_metadata <- merge_metadata[merge_metadata$celltype %in% c('GMP', 'CD14.Mono'),]
merge_metadata$celltype <- factor(merge_metadata$celltype, level = c('GMP', 'CD14.Mono'))

p1 <- ggplot(merge_metadata, aes(x = celltype, y = ETS1)) + 
     			geom_violin(aes(fill = celltype)) +
     			geom_jitter(width = 0.15, size = 0.1) + 
    			geom_signif(comparisons = list(c("GMP", "CD14.Mono")), 
                			test = "wilcox.test") + 
    			labs(x = '', y = 'Normalized expression level') + 
    			theme_bw() + 
    			theme(panel.grid.major = element_blank(),
		          panel.grid.minor = element_blank()) + NoLegend()
ggsave('result2/bmmc_pairs/ETS1.pdf',p1,dpi = 300, width = 6, height = 3)

p1 <- ggplot(merge_metadata, aes(x = celltype, y = RASGEF1B)) + 
     			geom_violin(aes(fill = celltype)) +
     			geom_jitter(width = 0.15, size = 0.1) + 
    			geom_signif(comparisons = list(c("GMP", "CD14.Mono")), 
                			test = "wilcox.test") + 
    			labs(x = '', y = 'Normalized expression level') + 
    			theme_bw() + 
    			theme(panel.grid.major = element_blank(),
		          panel.grid.minor = element_blank()) + NoLegend()
ggsave('result2/bmmc_pairs/RASGEF1B.pdf',p1,dpi = 300, width = 6, height = 3)

p1 <- ggplot(merge_metadata, aes(x = celltype, y = RETN)) + 
     			geom_violin(aes(fill = celltype)) +
     			geom_jitter(width = 0.15, size = 0.1) + 
    			geom_signif(comparisons = list(c("GMP", "CD14.Mono")), 
                			test = "wilcox.test") + 
    			labs(x = '', y = 'Normalized expression level') + 
    			theme_bw() + 
    			theme(panel.grid.major = element_blank(),
		          panel.grid.minor = element_blank()) + NoLegend()
ggsave('result2/bmmc_pairs/RETN.pdf',p1,dpi = 300, width = 6, height = 3)


### for test
# merge_metadata <- bmmc@meta.data
# merge_metadata$RASGEF1B <- bmmc@assays$RNA$data['RASGEF1B',]
# merge_metadata <- merge_metadata[merge_metadata$celltype %in% c('GMP', 'CD14.Mono'),]
# merge_metadata$celltype <- factor(merge_metadata$celltype, level = c('GMP', 'CD14.Mono'))

# ggplot(merge_metadata, aes(x = celltype, y = RASGEF1B)) + 
#      			geom_violin(aes(fill = celltype)) +
#      			geom_jitter(width = 0.15, size = 0.1) + 
#     			geom_signif(comparisons = list(c("GMP", "CD14.Mono")), 
#                 			test = "wilcox.test") + 
#     			labs(x = '', y = 'Normalized expression level') + 
#     			theme_bw() + 
#     			theme(panel.grid.major = element_blank(),
# 		          panel.grid.minor = element_blank()) + NoLegend()