library(tidyverse)
# For some of the plots:
library(igraph)
library(ggraph)
library(tidygraph)
library(ggrepel)

setwd('~/Desktop/Research/role_singlecell/')
set.seed(123)


#--------------------------------------------- Function ---------------------------------------------
createNetworkNode <- function(edge_df, anchor){

	#Input: edge_df - A data frame containing the edges of a network. 
	#Output: A data frame representing all the nodes in the network.

	node_df <- data.frame(label = unique(c(edge_df$from,edge_df$to)),
												index = 1:length(unique(c(edge_df$from,edge_df$to))))
												# The 'index' column is a sequence of integers 
												# from 1 to the number of unique labels.
	node_df$color <- 'orange'
	node_df$color[node_df$label %in% anchor] <- 'purple'

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

	p1 <- ggraph(net,layout = 'kk', maxiter = 1000) + 
				geom_node_text(aes(label = name), size = 6) +
				geom_edge_link(aes(color = as.factor(signed)), 
											 start_cap = circle(5, 'mm'),
											 end_cap = circle(5, 'mm'),
											 alpha = 0.5,
											 width = 2, 
											 show.legend=FALSE) +
				scale_edge_color_manual(values = c('+' = "black", '-' = "#8A2BE2")) + #BlueViolet
    		geom_node_point(aes(fill = color, color = color),size = 10, shape = 21,alpha = 0.6, show.legend=FALSE)  + 
    		scale_color_manual(values = c('orange' = "#FF4500", 'purple' = "purple")) + 
    		scale_fill_manual(values = c('orange' = "#FF4500", 'purple' = "purple")) + #OrangeRed1
    		theme_graph()

  return(p1)
}

plotEmbedding <- function(embedding){
	p1 <- ggplot(data = embedding, mapping = aes(x = Dim1, y = Dim2)) + 
			 				 geom_point(aes(color = color), alpha = 0.6,size =8) + 
			 				 scale_color_manual(values = c('orange' = "#FF4500", 'purple' = "purple")) + 
			 				 geom_text_repel(aes(label = name), 
			 				 				 size = 8,
			 				 				 max.overlaps = Inf,
			 				 				 #nudge_x = 0.5,
			 				 				# nudge_y = 0.5,
			 				 				 box.padding = 0.5, 
			 				 				 segment.size = 0.1) + 
			 				 theme_bw() +
			 				 theme(panel.grid.major = element_blank(),
			 				 	   panel.grid.minor = element_blank(),
			 				 	   axis.title.x = element_blank(),
			 				 	   axis.title.y = element_blank(),
			 				 	   axis.text.x = element_blank(),
			 				 	   axis.text.y = element_blank(),
			 				 	   axis.ticks = element_blank(),
			 				 	   legend.position = "none"
			 				 	   )
	return(p1)
}

#------------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------------

#--------------------------gsd
gsd <- read.table('data/curated/GSD.csv', sep = ',', header = TRUE)
colnames(gsd) <- c('from', 'to', 'signed')
#preprocessing
gsd <- gsd[!(gsd$from == gsd$to),]
gsd$from_to <- apply(gsd[,1:2],1,function(x) paste(sort(x), collapse = '-'))
gsd <- gsd[!duplicated(gsd[,c(3,4)]),]
gsd <- gsd[,-4]

anchor <- c('DHH', 'PGD2')
gsd_node <- createNetworkNode(gsd,anchor)

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
ggsave('results/curated/gsd_network.pdf',p1,dpi = 300, width = 7, height = 5)

#-------------------------vsc
vsc <- read.table('data/curated/VSC.csv', sep = ',', header = TRUE)
colnames(vsc) <- c('from', 'to', 'signed')

#preprocessing
vsc <- vsc[!(vsc$from == vsc$to),]
vsc$from_to <- apply(vsc[,1:2],1,function(x) paste(sort(x), collapse = '-'))
vsc <- vsc[!duplicated(vsc[,c(3,4)]),]
vsc <- vsc[,-4]

anchor <- c('Nkx62', 'Dbx1')
vsc_node <- createNetworkNode(vsc,anchor)

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

ggsave('results/curated/vsc_network.pdf',p1,dpi = 300, width = 7, height = 5)

#----------------------------------hsc
hsc <- read.table('data/curated/hsc.csv', sep = ',', header = TRUE)
colnames(hsc) <- c('from', 'to', 'signed')

#preprocessing
hsc <- hsc[!(hsc$from == hsc$to),]
hsc$from_to <- apply(hsc[,1:2],1,function(x) paste(sort(x), collapse = '-'))
hsc <- hsc[!duplicated(hsc[,c(3,4)]),]
hsc <- hsc[,-4]

anchor <- c('Fli1', 'Eklf', 'EgrNab', 'cJun')
hsc_node <- createNetworkNode(hsc, anchor)

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
ggsave('results/curated/hsc_network.pdf',p1,dpi = 300, width = 7, height = 5)

#--------------------------mcad
mcad <- read.table('data/curated/mCAD.csv', sep = ',', header = TRUE)
colnames(mcad) <- c('from', 'to', 'signed')

#preprocessing
mcad <- mcad[!(mcad$from == mcad$to),]
mcad$from_to <- apply(mcad[,1:2],1,function(x) paste(sort(x), collapse = '-'))
mcad <- mcad[!duplicated(mcad[,c(3,4)]),]
mcad <- mcad[,-4]

anchor <- c('Sp8')
mcad_node <- createNetworkNode(mcad, anchor)

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

ggsave('results/curated/mcad_network.pdf',p1,dpi = 300, width = 7, height = 5)

#----------------------------------sim_tree
sim_tree <- read.table('data/curated/tree_5_sign.edgelist', sep = ' ', header = FALSE)
colnames(sim_tree) <- c('from', 'to', 'signed')

sim_tree$from <- paste0('S',sim_tree$from)
sim_tree$to <- paste0('S',sim_tree$to)

anchor <- c('S15', 'S9')
sim_tree_node <- createNetworkNode(sim_tree,anchor)

sim_tree$signed <- ifelse(sim_tree$signed == -1, "-", "+")
sim_tree_edgelist <- createEdgeList(sim_tree, sim_tree_node)

write.table(sim_tree_edgelist,
						'data/curated/tree_5.edgelist', 
						sep = '\t', 
						row.names = FALSE, 
						col.names = FALSE, 
						quote = FALSE)

write.table(sim_tree_node,
						'data/curated/tree_5_index.txt', 
						sep = '\t', 
						row.names = FALSE, 
						col.names = FALSE, 
						quote = FALSE)

net <- createGraph(sim_tree,sim_tree_node)
p1 <- plotGraph(net)

ggsave('results/curated/tree_5_network.pdf',p1,dpi = 300, width = 6, height = 5)
x <- 'tree_5'
beside_embed <- read.table(paste0('mailbox/emb/curated/', x,'_BESIDE.emb'), 
														sep = ' ', 
														header = FALSE, 
														skip = 1)[,c(1,2,3)]

sim_tree_node$index <- sim_tree_node$index-1

colnames(beside_embed) <- c('index', 'Dim1', 'Dim2')

#------------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------------

#----------------plot embedding
network_name <- c('gsd', 'vsc', 'hsc', 'mcad', 'tree_5')

node_list <- list()
for (x in network_name){
	node_index <- read.table(paste0('data/curated/', x, '_index.txt'), sep = '\t')
	colnames(node_index) <- c('label', 'index', 'color')
	if (x  == 'tree_5'){
		node_index$index <- node_index$index-1
	}
	rownames(node_index) <- node_index$index
	node_list[[x]] <- node_index
}

# beside
for (x in network_name){
	embedding_path <- paste0('mailbox/emb/curated/', x,'_BESIDE.emb')
	embed <- read.table(embedding_path, sep = ' ', header = FALSE, skip = 1)[,c(1,2,3)]
	colnames(embed) <- c('index', 'Dim1', 'Dim2')

	if(x == 'tree_5'){
		embed$index <- as.character(embed$index)
	}
	embed$name <- node_list[[x]][embed$index,]$label
	embed$color <- node_list[[x]][embed$index,]$color

	p1 <- plotEmbedding(embed)
	save_path <- paste0('results/curated/', x, '_beside.pdf')
	ggsave(save_path,p1,dpi = 300, width = 5, height = 5)
}

# struc2vec
for (x in network_name){
	embedding_path <- paste0('mailbox/emb/curated/', x,'_graph_struc2vec.emb')
	embed <- read.table(embedding_path,  sep = ' ', header = FALSE, skip = 1)
	colnames(embed) <- c('index', 'Dim1', 'Dim2')
	if(x == 'tree_5'){
		embed$index <- as.character(embed$index)
	}
	embed$name <- node_list[[x]][embed$index,]$label
	embed$color <- node_list[[x]][embed$index,]$color

	p1 <- plotEmbedding(embed)
	save_path <- paste0('results/curated/', x, '_s2v.pdf')
	ggsave(save_path,p1,dpi = 300, width = 5, height = 5)
}


# g2r
for (x in network_name){
	embedding_path <- paste0('mailbox/emb/curated/', x,'.emb')
	embed <- read.table(embedding_path, sep = ' ', header = FALSE, skip = 1)
	colnames(embed) <- c('index', 'Dim1', 'Dim2')
	if(x == 'tree_5'){
		embed$index <- as.character(embed$index)
	}
	embed$name <- node_list[[x]][embed$index,]$label
	embed$color <- node_list[[x]][embed$index,]$color

	p1 <- plotEmbedding(embed)
	save_path <- paste0('results/curated/', x, '_g2r.pdf')
	ggsave(save_path,p1,dpi = 300, width = 5, height = 5)
}
