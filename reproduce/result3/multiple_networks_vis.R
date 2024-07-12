library(tidyverse)
library(ggsignif)
library(Seurat)

setwd('~/Desktop/Research/role_singlecell/results/')
set.seed(123)


################################################################################################

pbmc_distance <- read.table('result2/pbmc_all_distance.txt', sep = '\t', header = TRUE)
colnames(pbmc_distance)[1] <- 'name'
pbmc_distance$label <- ifelse(pbmc_distance$name  == 'TNFRSF13B', pbmc_distance$name, NA)

##
plotEmbedding <- function(embedding){
	p1 <- ggplot(data = embedding, mapping = aes(x = distance_avg, y = distance_sd, color = color)) + 
			 				 geom_point(alpha = 0.9,size =3) + 
			 				 theme_bw() +
			 				 theme(panel.grid.major = element_blank(),
			 				 	   panel.grid.minor = element_blank(),
			 				 	   legend.position = "none"
			 				 	   ) +
			 				 labs(x = 'Average distance', y = 'Standard deviation of distance')
	return(p1)
}

pbmc_distance$color <- ifelse((pbmc_distance$distance_avg > quantile(pbmc_distance$distance_avg, 
															 					0.95, na.rm = TRUE)  | pbmc_distance$distance_sd > 0.4), 
															'DTGs', 'non-DTGs')

pbmc_DTGs <- pbmc_distance[pbmc_distance$color == 'DTGs',]$name

p1 <- plotEmbedding(pbmc_distance) +
			 				 geom_hline(yintercept = 0.4, linetype = "dashed", color = "red", size = 0.5) +
			 				 scale_color_manual(values=c("DTGs"="#FF0000", "non-DTGs"="#69b3a2"))  +
			 				 geom_vline(xintercept=3.35,linetype = "dashed", color = "red", size = 0.5) +
			 				 # geom_text_repel(aes(label = label), 
			 				 # 				 size = 7,
			 				 # 				 max.overlaps = Inf,
			 				 # 				 #nudge_x = 0.5,
			 				 # 				# nudge_y = 0.5,
			 				 # 				 box.padding = 0.5, 
			 				 # 				 segment.size = 0.1) + 
			 				theme(legend.position = "right")

ggsave('single_cell/pbmc/mean_sd.pdf',p1,dpi = 300, width = 4.8, height = 4)

# pbmc <- readRDS('data/singleCell/pbmc_rna/pbmc_rna.rds')

# pbmc <- pbmc %>% 
#     SCTransform() %>%
#     RunPCA() %>%
#     FindNeighbors(dims = 1:20) %>%
#     RunUMAP(dims = 1:20)

# saveRDS(pbmc,'data/singleCell/pbmc_rna/pbmc_processed.rds')	

Idents(pbmc) <- pbmc$celltype
pbmc.marker <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.marker <- pbmc.marker[pbmc.marker$p_val_adj < 0.05,]

write.table(pbmc.marker,
					 'data/singleCell/pbmc_rna/pbmc.markers.txt', 
					  quote = FALSE, 
					  row.names = FALSE)

#intersect(unique(pbmc.marker$gene), pbmc_DTGs)

pbmc_non_DTGs <- pbmc_distance[pbmc_distance$color == 'non-DTGs',]

pbmc_non_DTGs <- pbmc_non_DTGs[pbmc_non_DTGs$name %in% unique(pbmc.marker$gene),]


################################################################################################

multi_omics_mep <- read.table('result2/multiomics_exp331_mep_pair.txt', sep = '\t', header = TRUE)
multi_omics_mep_sub <- multi_omics_mep[,c(1,13:22)]
rownames(multi_omics_mep_sub) <- multi_omics_mep_sub$Gene
multi_omics_mep_sub <- multi_omics_mep_sub[,-1]

multi_omics_mep_sub$mean <- rowMeans(multi_omics_mep_sub)
multi_omics_mep_sub$sd <- apply(multi_omics_mep_sub[,-11], 1, sd)

#scatter plot
multi_omics_mep_sub_point <- multi_omics_mep_sub[,c(11,12)]
colnames(multi_omics_mep_sub_point) <- c('distance_avg', 'distance_sd')

multi_omics_mep_sub_point$color <- ifelse((multi_omics_mep_sub_point$distance_avg > quantile(multi_omics_mep_sub_point$distance_avg, 
															 					0.95, na.rm = TRUE)  | multi_omics_mep_sub_point$distance_sd > 1.5), 
															'DTGs', 'non-DTGs')

p1 <- plotEmbedding(multi_omics_mep_sub_point) +
			 				 geom_hline(yintercept = 1.5, linetype = "dashed", color = "red", size = 0.5) +
			 				 scale_color_manual(values=c("DTGs"="#FF0000", "non-DTGs"="#69b3a2"))  +
			 				 geom_vline(xintercept=4.96,linetype = "dashed", color = "red", size = 0.5) +
			 				 theme(legend.position = "right")

ggsave('result2/mep_mean_sd.pdf',p1,dpi = 300, width = 4.8, height = 3.5)

high_mean_genes <- rownames(multi_omics_mep_sub[order(multi_omics_mep_sub$mean, decreasing = TRUE),])[5:1]
high_sd_genes <-  rownames(multi_omics_mep_sub[order(multi_omics_mep_sub$sd, decreasing = TRUE),])[1:5]

complex_plot_df  <- multi_omics_mep_sub[c(high_sd_genes,high_mean_genes), 1:10]

library(ComplexHeatmap)

complex_plot_df <- as.matrix(complex_plot_df)

pdf('result2/MEPtopgenes_heatmap.pdf', width = 6.5, height = 2)
Heatmap(complex_plot_df,
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        show_column_names = FALSE,
        row_names_side = "left",
        row_split = c(rep("A",5),rep("B", 5)),
        heatmap_legend_param = list(
            title = "Distance",
            position = "bottom",
            direction = "horizontal"
            )
)
dev.off()


# gmp 
multi_omics_gmp <- read.table('result2/multiomics_exp331_gmp_pair.txt', sep = '\t', header = TRUE)

multi_omics_gmp_sub_point <- multi_omics_gmp[,c(15,16)]
colnames(multi_omics_gmp_sub_point) <- c('distance_avg', 'distance_sd')
multi_omics_gmp_sub_point$color <- ifelse((multi_omics_gmp_sub_point$distance_avg > quantile(multi_omics_gmp_sub_point$distance_avg, 
															 					0.95, na.rm = TRUE)  | multi_omics_gmp_sub_point$distance_sd > 1.5), 
															'DTGs', 'non-DTGs')
p1 <- plotEmbedding(multi_omics_gmp_sub_point) +
			 				 geom_hline(yintercept = 1.5, linetype = "dashed", color = "red", size = 0.5) +
			 				 scale_color_manual(values=c("DTGs"="#FF0000", "non-DTGs"="#69b3a2"))  +
			 				 geom_vline(xintercept=5.00502,linetype = "dashed", color = "red", size = 0.5) +
			 				 theme(legend.position = "right")
ggsave('result2/gmp_mean_sd.pdf',p1,dpi = 300, width = 4.8, height = 3.5)


high_mean_genes <- rownames(multi_omics_mep_sub[order(multi_omics_mep_sub$mean, decreasing = TRUE),])[1:5]
high_sd_genes <-  rownames(multi_omics_mep_sub[order(multi_omics_mep_sub$sd, decreasing = TRUE),])[1:5]

rownames(multi_omics_gmp) <- multi_omics_gmp$Gene
multi_omics_gmp <- multi_omics_gmp[,-1]

high_mean_genes <- rownames(multi_omics_gmp[order(multi_omics_gmp$mean, decreasing = TRUE),])[1:5]
high_sd_genes <-  rownames(multi_omics_gmp[order(multi_omics_gmp$sd, decreasing = TRUE),])[1:5]

complex_plot_df <- multi_omics_gmp[c(high_sd_genes,high_mean_genes),c(8:13)]
complex_plot_df <- as.matrix(complex_plot_df)

pdf('result2/GMPtopgenes_heatmap.pdf', width = 6.5, height = 2)
Heatmap(complex_plot_df,
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        show_column_names = FALSE,
        row_names_side = "left",
        row_split = c(rep("A",5),rep("B", 5)),
        heatmap_legend_param = list(
            title = "Distance",
            position = "bottom",
            direction = "horizontal"
            )
)
dev.off()
