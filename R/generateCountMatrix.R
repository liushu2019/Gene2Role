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

pbmc_count <- pbmc_count[VariableFeatures(pbmc),]
fwrite(x = pbmc_count, file = "data/singleCell/pbmc_rna/count.csv",row.names = TRUE)

# write metadata
pbmc_metadata <- pbmc_orig@meta.data
write.table(pbmc_metadata, 
						file = 'data/singleCell/pbmc_rna/metadata.csv', 
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


#------------------MMBRAIN----------------
mmbrain_orig <- readRDS('data/singleCell/mmbrain_rna/mmbrain_rna.rds')

# write count
mmbrain_count <- as.data.frame(as.matrix(mmbrain_orig@assays$RNA@counts))
mmbrain <- CreateSeuratObject(mmbrain_count)
mmbrain <- NormalizeData(mmbrain, normalization.method = "LogNormalize", scale.factor = 10000)
mmbrain <- FindVariableFeatures(mmbrain, selection.method = "vst", nfeatures = 2000)

mmbrain_count <- mmbrain_count[VariableFeatures(mmbrain),]
fwrite(x = mmbrain_count, file = "data/singleCell/mmbrain_rna/count.csv", row.names = TRUE)

# write metadata
mmbrain_metadata <- mmbrain_orig@meta.data
mmbrain_metadata$celltype <- gsub("/", "_", mmbrain_metadata$celltype)
mmbrain_metadata$celltype <- gsub(" ", "", mmbrain_metadata$celltype)

write.table(mmbrain_metadata, 
						file = 'data/singleCell/mmbrain_rna/metadata.csv', 
						sep = ',',
						row.names = TRUE, 
						quote = FALSE)


#----------------MMLUNG--------------------
mmlung_orig <- readRDS('data/singleCell/mmlung_rna/mmlung_rna.rds')
share_genes <- readRDS('data/singleCell/mmlung_rna/share_genes.rds')

# write count
mmlung_count <- as.data.frame(as.matrix(mmlung_orig@assays$RNA@counts))
rownames(mmlung_count) <- share_genes
mmlung <- CreateSeuratObject(mmlung_count)
mmlung <- NormalizeData(mmlung, normalization.method = "LogNormalize", scale.factor = 10000)
mmlung <- FindVariableFeatures(mmlung, selection.method = "vst", nfeatures = 2000)

# write count
mmlung_count <- mmlung_count[VariableFeatures(mmlung),]
fwrite(x = mmlung_count, file = "data/singleCell/mmlung_rna/count.csv",row.names = TRUE)

# write metadata
mmlung_metadata <- mmlung_orig@meta.data
write.table(mmlung_metadata, 
						file = 'data/singleCell/mmlung_rna/metadata.csv', 
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






