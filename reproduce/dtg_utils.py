import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from scipy.spatial.distance import euclidean
from itertools import combinations


def cal_pair_distance(row, embedding):

	# Extract indices
	index1, index2 = row[1], row[2]

	# Check if indices exist in embedding
	if index1 in embedding.index and index2 in embedding.index:
		return euclidean(embedding.loc[index1], embedding.loc[index2])
	else:
		return np.nan

def calculate_avg_sd_distances(index, embedding):
    
    # Prepare columns for average and standard deviation of distances
	index_file_copy = index.copy()
	index_file_copy['distance_avg'] = np.nan
	index_file_copy['distance_sd'] = np.nan
    
    # Get the names of all cell type columns
	cell_names = list(index.columns)[1:]
	#index_file_copy[cell_names] = index_file_copy[cell_names].astype('float64')
    # Initialize a list to collect indices of rows to drop
	rows_to_drop = []
    
	for idx, row in index_file_copy.iterrows():
		gene_embeddings = []# To store embeddings for the current gene's cell types
		gene_index_list = row[cell_names] # Get cell type index for the current gene
		distances = []
        
		if gene_index_list.isin(embedding.index).all():
			for gene_index in gene_index_list:
				gene_embeddings.append(embedding.loc[gene_index].values)
    
			#calculate the center point
			center_point = np.mean(gene_embeddings, axis=0)
    
            # Calculate the distance from each cell type's embedding to the center point
			for cell_name, gene_index in zip(cell_names, gene_index_list):
                
				gene_embedding = embedding.loc[gene_index].values
				distance = np.linalg.norm(gene_embedding - center_point)
	
				#index_file_copy.at[idx, cell_name] = distance
				distances.append(distance)
			
			index_file_copy.at[idx, 'distance_avg'] = np.mean(distances)
			index_file_copy.at[idx, 'distance_sd'] = np.std(distances)
		
		else:
			rows_to_drop.append(idx)
            
	index_file_copy.drop(rows_to_drop, inplace=True)
            
	return index_file_copy

def cal_node_distances(index_file, embedding):

	num_cell_types = index_file.shape[1]-1
    
    
	if num_cell_types == 2:
        
		index_file_copy = index_file.copy()
		# Calculate distance for pairs of indices
		index_file_copy['distance'] = index_file_copy.apply(cal_pair_distance, axis=1, embedding=embedding)
		index_file_copy = index_file_copy.sort_values(by='distance', ascending=False, na_position='last')
		plot_frequencies(index_file_copy, 'distance')

	elif num_cell_types >= 3:

		index_file_copy = calculate_avg_sd_distances(index_file, embedding)
        
		plot_avg_sd(index_file_copy)
        
	else:
		"The input index file is wrong!"

	return index_file_copy

def plot_frequencies(df, column_names = 'distance'):
	plt.hist(df[column_names], bins=30)
	plt.xlabel('Distance')
	plt.ylabel('Frequency')
	plt.show()

def plot_avg_sd(df,avg = 'distance_avg', sd = 'distance_sd'):
	ax = df.plot(kind='scatter', x=avg, y=sd, color='blue', alpha=0.5)
	ax.set_xlabel('Average Distance')
	ax.set_ylabel('SD of Distance')

	plt.show()





