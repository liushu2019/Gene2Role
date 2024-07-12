library(Seurat)
library(data.table)

setwd('~/Desktop/Research/role_singlecell/')
set.seed(123)

#-------------------PBMC----------------


pbmc_orig <- readRDS('data/singleCell/pbmc_rna/pbmc_rna.rds')
# write count
pbmc_count <- as.data.frame(as.matrix(pbmc_orig@assays$RNA@counts))
pbmc <- CreateSeuratObject(pbmc_count)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)


# Generate count matrix for multiple GRN comparison.
pbmc_count <- pbmc_count[VariableFeatures(pbmc),]
fwrite(x = pbmc_count, file = "data/singleCell/pbmc_rna/count.csv",row.names = TRUE)

# write metadata
pbmc_metadata <- pbmc_orig@meta.data
write.table(pbmc_metadata, 
	    file = 'data/singleCell/pbmc_rna/metadata.csv', 
	    sep = ',',
	    row.names = TRUE, 
	    quote = FALSE)
 

# generate pair  GRN comparison.
Idents(pbmc_orig) <- pbmc_orig$celltype

pbmc_CD4_orig <- subset(pbmc_orig, celltype %in% c('22_CD4.M', 'CD4.N'))
pbmc_CD4_count <- as.data.frame(as.matrix(pbmc_CD4_orig@assays$RNA@counts))

pbmc_CD4 <- CreateSeuratObject(pbmc_CD4_count)
pbmc_CD4 <- NormalizeData(pbmc_CD4, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc_CD4 <- FindVariableFeatures(pbmc_CD4, selection.method = "vst", nfeatures = 2000)

pbmc_CD4_count <- pbmc_CD4_count[VariableFeatures(pbmc_CD4),]
fwrite(x = pbmc_CD4_count, file = "data/singleCell/pbmc_rna/PBMC_CD4_count.csv",row.names = TRUE)

pbmc_CD4_metadata <- pbmc_CD4_orig@meta.data
write.table(pbmc_CD4_metadata, 
						file = 'data/singleCell/pbmc_rna/PBMC_CD4_metadata.csv', 
						sep = ',',
						row.names = TRUE, 
						quote = FALSE)


#-------------------BMMC----------------
bmmc_orig <- readRDS('data/singleCell/bmmc_rna/bmmc_rna.rds')

bmmc_count <- as.data.frame(as.matrix(bmmc_orig@assays$RNA@counts))
bmmc <- CreateSeuratObject(bmmc_count)
bmmc <- NormalizeData(bmmc, normalization.method = "LogNormalize", scale.factor = 10000)
bmmc <- FindVariableFeatures(bmmc, selection.method = "vst", nfeatures = 2000)

# write count
bmmc_count <- bmmc_count[VariableFeatures(bmmc),]
fwrite(x = bmmc_count, file = "data/singleCell/bmmc_rna/count.csv",row.names = TRUE)

# write metadata
bmmc_metadata <- bmmc_orig@meta.data
write.table(bmmc_metadata, 
						file = 'data/singleCell/bmmc_rna/metadata.csv', 
						sep = ',',
						row.names = TRUE, 
						quote = FALSE)

# pair cell types
Idents(bmmc_orig) <- bmmc_orig$celltype
bmmc_orig_gmp <- subset(bmmc_orig, celltype %in% c('GMP', 'CD14.Mono'))
bmmc_gmp_count <- as.data.frame(as.matrix(bmmc_orig_gmp@assays$RNA@counts))

bmmc_gmp <- CreateSeuratObject(bmmc_gmp_count)
bmmc_gmp <- NormalizeData(bmmc_gmp, normalization.method = "LogNormalize", scale.factor = 10000)
bmmc_gmp <- FindVariableFeatures(bmmc_gmp, selection.method = "vst", nfeatures = 2000)

bmmc_gmp_count <- bmmc_gmp_count[VariableFeatures(bmmc_gmp),]
fwrite(x = bmmc_gmp_count, file = "data/singleCell/bmmc_rna/bmmc_gmp_count.csv",row.names = TRUE)

bmmc_gmp_metadata <- bmmc_orig_gmp@meta.data
write.table(bmmc_gmp_metadata, 
						file = 'data/singleCell/bmmc_rna/bmmc_gmp_metadata.csv', 
						sep = ',',
						row.names = TRUE, 
						quote = FALSE)


#---------------- Glioblastoma --------------------
library(biomaRt)


t0 <- read.csv('data/singleCell/glioblastoma/GSM4292484_filtered_gene_bc_matrices_P10.csv', 
								header = TRUE, row.names = 1)
colnames(t0) <- paste0('0h_', colnames(t0))

t12 <- read.csv('data/singleCell/glioblastoma/GSM4292485_filtered_gene_bc_matrices_12h.csv', 
                header = TRUE, row.names = 1)
colnames(t12) <- paste0('t12_', colnames(t12))
merge_count <- cbind(t0, t12)
merge_seurat <- CreateSeuratObject(merge_count)


merge_seurat <- NormalizeData(merge_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
merge_seurat <- FindVariableFeatures(merge_seurat, selection.method = "vst", nfeatures = 2000)

merge_seurat <- ScaleData(merge_seurat)
merge_seurat <- RunPCA(merge_seurat, features = VariableFeatures(object = merge_seurat))

merge_seurat <- RunUMAP(merge_seurat, dims = 1:20)
saveRDS(merge_seurat,'data/singleCell/glioblastoma/merge_seurat.rds')

merge_seurat.markers <- FindAllMarkers(merge_seurat, only.pos = TRUE)
merge_seurat.markers <- merge_seurat.markers[merge_seurat.markers$p_val_adj < 0.05,]
write.table(merge_seurat.markers,
					 'data/singleCell/glioblastoma/merge_seurat.markers.txt', 
					  quote = FALSE, 
					  row.names = FALSE)

merge_count <- merge_count[VariableFeatures(merge_seurat),]

ensembl_gene <- read.table('../metadata/formatted_genes.txt', sep = '\t')

rownames(ensembl_gene) <- ensembl_gene$ensembl_gene_id
shared_genes <- intersect(ensembl_gene$V1, rownames(merge_count))
merge_count <- merge_count[shared_genes,]

rownames(merge_count) <- ensembl_gene[rownames(merge_count),]$V2

fwrite(x = merge_count, file = "data/singleCell/glioblastoma/merge_count.csv",row.names = TRUE)

merge_seurat$celltype <- merge_seurat$orig.ident
merge_metadata <- merge_seurat@meta.data

write.table(merge_metadata, 
						file = 'data/singleCell/glioblastoma/metadata.csv', 
						sep = ',',
						row.names = TRUE, 
						quote = FALSE)






