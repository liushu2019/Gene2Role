library(tidyverse)
# For some of the plots:
library(igraph)
library(ggraph)
library(tidygraph)
library(ggrepel)

setwd('~/Desktop/Research/role_singlecell/')
set.seed(123)
#--------------------------------------------- Function ---------------------------------------------

createNetworkNode <- function(edge_df){

	#Input: edge_df - A data frame containing the edges of a network. 
	#Output: A data frame representing all the nodes in the network.

	node_df <- data.frame(label = unique(c(edge_df$from,edge_df$to)),
												index = 1:length(unique(c(edge_df$from,edge_df$to))))
												# The 'index' column is a sequence of integers 
												# from 1 to the number of unique labels.
	rownames(node_df) <- node_df$label
	return(node_df)

}

createEdgeList <- function(edge_df, node_df){

	edgelist <- edge_df

	# change the sign representation
	edgelist$newsigned <- '1'
	edgelist[edgelist$signed == '-',]$newsigned <- -1

	edgelist <- edgelist[,c(1,2,4)]

	# change the name
	rownames(node_df) <- node_df$label
	edgelist$from <- node_df[edgelist$from,]$index
	edgelist$to <- node_df[edgelist$to,]$index

	return(edgelist)
}

createGraph <- function(edge_df,node_df){

	#Input:
		#edge_df - A data frame representing the edges of a network. 
		#node_df - A data frame representing the nodes of the network.
	# Output: 
	#	A graph object representing the network

	net <- graph_from_data_frame(d=edge_df,
																			 vertices=node_df,
																			 directed=FALSE)

	net <- as_tbl_graph(net) %>% 
  							 mutate(deg = centrality_degree(mode='all'))
  return(net)
}

plotGraph <- function(net){
	# Input: A graph object.
	# Output: A ggplot object representing the plotted network graph.

	p1 <- ggraph(net,layout = 'kk') + 
				geom_edge_link(aes(color = signed), show.legend=FALSE) +
				scale_edge_color_manual(values = c('+' = "red", '-' = "blue")) +
    		geom_node_point(aes(size = deg),fill = 'grey',shape = 21,color = 'grey' ,show.legend=FALSE)  + 
    		geom_node_text(aes(label = name)) + 
    		theme_graph()

  return(p1)
}

plotEmbedding <- function(embedding){
	p1 <- ggplot(data = embedding, mapping = aes(x = Dim1, y = Dim2)) + 
			 				 geom_point() + 
			 				 geom_text_repel(aes(label = name), size = 3) + 
			 				 theme_bw() +
			 				 theme(panel.grid.major = element_blank(),
			 				 			 panel.grid.minor = element_blank())
	return(p1)
}

#------------------------------------------------------------------------------------------------



#-----------------------------------------
gsd <- read.table('data/curated/GSD.csv', sep = ',', header = TRUE)
colnames(gsd) <- c('from', 'to', 'signed')
#preprocessing
gsd <- gsd[!(gsd$from == gsd$to),]
gsd$from_to <- apply(gsd[,1:2],1,function(x) paste(sort(x), collapse = '-'))
gsd <- gsd[!duplicated(gsd[,c(3,4)]),]
gsd <- gsd[,-4]

gsd_node <- createNetworkNode(gsd)

# save to edgelist.
gsd_edgelist <- createEdgeList(gsd, gsd_node)

write.table(gsd_edgelist,
						'data/curated/gsd.edgelist', 
						sep = '\t', 
						row.names = FALSE, 
						col.names = FALSE, 
						quote = FALSE)

write.table(gsd_node,
						'data/curated/gsd_index.txt', 
						sep = '\t', 
						row.names = FALSE, 
						col.names = FALSE, 
						quote = FALSE)

net <- createGraph(gsd,gsd_node)
p1 <- plotGraph(net)
ggsave('results/curated/gsd_network.pdf',p1,dpi = 300, width = 5, height = 5)

# read the embedding results.
gsd_embed <- read.table('tools/SignedS2V/emb/gsd.emb', sep = ' ', header = FALSE, skip = 1)
colnames(gsd_embed) <- c('index', 'Dim1', 'Dim2')

rownames(gsd_node) <- gsd_node$index
gsd_embed$name <- gsd_node[gsd_embed$index,]$label

p1 <- plotEmbedding(gsd_embed)
ggsave('results/curated/gsd_embedding.pdf',p1,dpi = 300, width = 5, height = 5)

#-----------------------------------------
vsc <- read.table('data/curated/VSC.csv', sep = ',', header = TRUE)
colnames(vsc) <- c('from', 'to', 'signed')

#preprocessing
vsc <- vsc[!(vsc$from == vsc$to),]
vsc$from_to <- apply(vsc[,1:2],1,function(x) paste(sort(x), collapse = '-'))
vsc <- vsc[!duplicated(vsc[,c(3,4)]),]
vsc <- vsc[,-4]

vsc_node <- createNetworkNode(vsc)

# save to edgelist.
vsc_edgelist <- createEdgeList(vsc, vsc_node)

write.table(vsc_edgelist,
						'data/curated/vsc.edgelist', 
						sep = '\t', 
						row.names = FALSE, 
						col.names = FALSE, 
						quote = FALSE)

write.table(vsc_node,
						'data/curated/vsc_index.txt', 
						sep = '\t', 
						row.names = FALSE, 
						col.names = FALSE, 
						quote = FALSE)

net <- createGraph(vsc,vsc_node)
p1 <- plotGraph(net)

ggsave('results/curated/vsc_network.pdf',p1,dpi = 300, width = 5, height = 5)


#-----------------------------------------
hsc <- read.table('data/curated/hsc.csv', sep = ',', header = TRUE)
colnames(hsc) <- c('from', 'to', 'signed')

#preprocessing
hsc <- hsc[!(hsc$from == hsc$to),]
hsc$from_to <- apply(hsc[,1:2],1,function(x) paste(sort(x), collapse = '-'))
hsc <- hsc[!duplicated(hsc[,c(3,4)]),]
hsc <- hsc[,-4]

hsc_node <- createNetworkNode(hsc)

hsc_edgelist <- createEdgeList(hsc, hsc_node)
write.table(hsc_edgelist,
						'data/curated/hsc.edgelist', 
						sep = '\t', 
						row.names = FALSE, 
						col.names = FALSE, 
						quote = FALSE)

write.table(hsc_node,
						'data/curated/hsc_index.txt', 
						sep = '\t', 
						row.names = FALSE, 
						col.names = FALSE, 
						quote = FALSE)

net <- createGraph(hsc,hsc_node)
p1 <- plotGraph(net)

ggsave('results/curated/hsc_network.pdf',p1,dpi = 300, width = 5, height = 5)

# read the embedding results.
hsc_embed <- read.table('tools/SignedS2V/emb/hsc.emb', sep = ' ', header = FALSE, skip = 1)
colnames(hsc_embed) <- c('index', 'Dim1', 'Dim2')

rownames(hsc_node) <- hsc_node$index
hsc_embed$name <- hsc_node[hsc_embed$index,]$label

p1 <- plotEmbedding(hsc_embed)
ggsave('results/curated/hsc_embedding.pdf',p1,dpi = 300, width = 5, height = 5)
#-----------------------------------------

mcad <- read.table('data/curated/mCAD.csv', sep = ',', header = TRUE)
colnames(mcad) <- c('from', 'to', 'signed')

#preprocessing
mcad <- mcad[!(mcad$from == mcad$to),]
mcad$from_to <- apply(mcad[,1:2],1,function(x) paste(sort(x), collapse = '-'))
mcad <- mcad[!duplicated(mcad[,c(3,4)]),]
mcad <- mcad[,-4]

mcad_node <- createNetworkNode(mcad)

mcad_edgelist <- createEdgeList(mcad, mcad_node)
write.table(mcad_edgelist,
						'data/curated/mcad.edgelist', 
						sep = '\t', 
						row.names = FALSE, 
						col.names = FALSE, 
						quote = FALSE)

write.table(mcad_node,
						'data/curated/mcad_index.txt', 
						sep = '\t', 
						row.names = FALSE, 
						col.names = FALSE, 
						quote = FALSE)

net <- createGraph(mcad,mcad_node)
p1 <- plotGraph(net)

ggsave('results/curated/mcad_network.pdf',p1,dpi = 300, width = 5, height = 5)

# read the embedding results.
mcad_embed <- read.table('tools/SignedS2V/emb/mcad.emb', sep = ' ', header = FALSE, skip = 1)
colnames(mcad_embed) <- c('index', 'Dim1', 'Dim2')

rownames(mcad_node) <- mcad_node$index
mcad_embed$name <- mcad_node[mcad_embed$index,]$label

p1 <- plotEmbedding(mcad_embed)
ggsave('results/curated/mcad_embedding.pdf',p1,dpi = 300, width = 5, height = 5)







