library(tidyverse)

setwd('~/Desktop/Research/role_singlecell/results/')
set.seed(123)
library(stats)
library(igraph)


my_colors = c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB",
					 "#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00")

#-------------------------------------UMAP--------------------------


#--------------------------function

plotEmbedding <- function(embedding){
	p1 <- ggplot(data = embedding, mapping = aes(x = UMAP1, y = UMAP2)) + 
			 				 geom_point(aes(color = Cluster),size =3) + 
			 				 theme_bw() +
			 				 scale_color_manual(values=my_colors) + 
			 				 theme(panel.grid.major = element_blank(),
			 				 	   panel.grid.minor = element_blank(),
			 				 	   axis.text.x = element_blank(),
			 				 	   axis.text.y = element_blank(),
			 				 	   axis.ticks = element_blank(),
			 				 	   legend.position = "right",
			 				 	   legend.text = element_text(size = 13)
			 				 	   ) + 
			 				guides(color = guide_legend(override.aes = list(size=5)))
	return(p1)
}
#--------------------------

B_eeisp_umap <- read.table('data/Figure2 Result1 UMAP _ Heatmap/df_fig2/df_cluster_umap_Exp1_2_2_10_2025.csv', 
												sep = '\t', 
												header = TRACE)
B_eeisp_umap$Cluster <- as.factor(B_eeisp_umap$Cluster)
p1 <- plotEmbedding(B_eeisp_umap)
ggsave('single_cell/B_eeisp_umap.pdf',p1,dpi = 300, width = 6, height = 5)

B_spearman_umap <- read.table('data/Figure2 Result1 UMAP _ Heatmap/df_fig2/df_cluster_umap_Exp1_2_1_10_2025.csv', 
												sep = '\t', 
												header = TRACE)
B_spearman_umap$Cluster <- as.factor(B_spearman_umap$Cluster)
p1 <- plotEmbedding(B_spearman_umap)
ggsave('single_cell/B_spearman_umap.pdf',p1,dpi = 300, width = 6, height = 5)

Ery0_umap <- read.table('data/Figure2 Result1 UMAP _ Heatmap/df_fig2/df_cluster_umap_Exp1_3_1_10_2025.csv', 
												sep = '\t', 
												header = TRACE)
Ery0_umap$Cluster <- as.factor(Ery0_umap$Cluster)
p1 <- plotEmbedding(Ery0_umap)
ggsave('single_cell/Ery0_umap.pdf',p1,dpi = 300, width = 6, height = 5)



#-------------------------------------HEATMAP--------------------------
library(ComplexHeatmap)

#--------------------------function

normalize_minmax <- function(x) {
    # x is a row of the matrix
    min_val <- min(x, na.rm = TRUE)
    max_val <- max(x, na.rm = TRUE)
    # Initialize a result vector of the same length, default to NA
    result <- rep(NA, length(x))
    # Set non-NA values to 0.5
    non_na_indices <- !is.na(x)
    result[non_na_indices] <- 0.5
    
    # If the maximum value is not equal to the minimum value, perform normalization
    if (min_val != max_val) {
        result[non_na_indices] <- (x[non_na_indices] - min_val) / (max_val - min_val)
    }
    
    return(result)
}

plot_complex_heatmap <- function(normalized_matrix, original_matrix){

    min_value <- min(normalized_matrix, na.rm = TRUE)
    max_value <- max(normalized_matrix, na.rm = TRUE)
    Heatmap(normalized_matrix, 
            cluster_rows = FALSE, 
            cluster_columns = TRUE, 
            na_col = 'grey', 
            show_column_dend = FALSE,
            column_names_side = "top",
            column_names_gp = gpar(fontsize = 9),
            column_names_rot = 0,
            heatmap_legend_param = list(
                title = "", # Remove legend title
                at = c(min_value, max_value), # Set legend tick positions
                labels = c("Min", "Max"), # Set legend tick labels
                legend_position = "topright" # Set legend position to top right
            ),
            cell_fun = function(j, i, x, y, width, height, fill) {
                if(is.na(original_matrix[i, j])) {
                    # If cell value is NA, display "NA"
                    grid.text("NA", x, y, gp = gpar(fontsize = 7))
                } else {
                    # If cell value is not NA, display the value and adjust font size to fit the cell
                    grid.text(sprintf("%.3f", original_matrix[i, j]), x, y, gp = gpar(fontsize = 7))
                }
            })
}

#--------------------------

rownames <- c('Degree centrality (+)', 'Degree centrality (-)', 'Betweenness centrality',
							'Eigenvector centrality', 'Degree assortativity (+)', 'Degree assortativity (-)',
							'Cluster coefficient (+)', 'Cluster coefficient (-)')

#122
B_eeisp_heatmap <- read.table('single_cell/df_fig2_r2/df_Heatmap_Exp1_2_2_10_2025_r2.csv', 
													 sep = '\t', 
													 header = TRUE)

colnames(B_eeisp_heatmap) <- c(0:9)
B_eeisp_heatmap <- round(B_eeisp_heatmap,3)
rownames(B_eeisp_heatmap) <- rownames
B_eeisp_heatmap <- as.matrix(B_eeisp_heatmap)

normalized_matrix <- t(apply(B_eeisp_heatmap, 1, normalize_minmax))
colnames(normalized_matrix) <- c(0:9)

pdf('single_cell/B_eeisp_heatmap.pdf', width = 6.5, height = 2)
plot_complex_heatmap(normalized_matrix, B_eeisp_heatmap)
dev.off()

#131
Ery0_heatmap <- read.table('single_cell/df_fig2_r2/df_Heatmap_Exp1_3_1_10_2025_r2.csv', 
													 sep = '\t', 
													 header = TRUE)

colnames(Ery0_heatmap) <- c(0:9)
Ery0_heatmap <- round(Ery0_heatmap,3)
rownames(Ery0_heatmap) <- rownames
Ery0_heatmap <- as.matrix(Ery0_heatmap)

normalized_matrix <- t(apply(Ery0_heatmap, 1, normalize_minmax))
colnames(normalized_matrix) <- c(0:9)

pdf('single_cell/131_heatmap.pdf', width = 6.5, height = 2)
plot_complex_heatmap(normalized_matrix, Ery0_heatmap)
dev.off()

#121
B_spearman_heatmap <- read.table('single_cell/df_fig2_r2/df_Heatmap_Exp1_2_1_10_2025_r2.csv', 
													 sep = '\t', 
													 header = TRUE)

colnames(B_spearman_heatmap) <- c(0:9)
B_spearman_heatmap <- round(B_spearman_heatmap,3)
rownames(B_spearman_heatmap) <- rownames
B_spearman_heatmap <- as.matrix(B_spearman_heatmap)

normalized_matrix <- t(apply(B_spearman_heatmap, 1, normalize_minmax))
colnames(normalized_matrix) <- c(0:9)

pdf('single_cell/121_heatmap.pdf', width = 6.5, height = 2)
plot_complex_heatmap(normalized_matrix, B_spearman_heatmap)
dev.off()


