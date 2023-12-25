import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

import argparse
import os 


print('---- Calculating gene-gene correlation by spearman correlation ----')

parser = argparse.ArgumentParser(description="Calculate and filter gene-gene Spearman correlations.")
parser.add_argument("file_path", 
                    type=str, 
                    help="Path to the CSV file containing gene expression data.")

parser.add_argument("out_file_name", 
                    type=str, 
                    help="out put file name")

parser.add_argument("-correlation_threshold", 
                    type=float, 
                    help="Threshold for filtering out top value*100\% correlations. (default:0.001)",
                    default=0.001)
parser.add_argument('--reindex', action='store_true', 
                    help='Flag for reindex gene names.')

args = parser.parse_args()

df = pd.read_csv(args.file_path, index_col=0).transpose()
# Calculate spearman correlation.
spearman_corr_matrix = df.corr(method='spearman')

# Flatten the correlation matrix and filter out self-correlations
correlations = spearman_corr_matrix.values.flatten()
correlations = correlations[~np.isnan(correlations)]
correlations = correlations[correlations != 1]

# Determine top pairs for positive and negative correlations
positive_threshold = np.percentile(correlations, 100 - int(args.correlation_threshold*100))
negative_threshold = np.percentile(correlations, int(args.correlation_threshold*100))
# print the number of positive and negative edges
num_positive_edges_top_percent = np.sum(correlations >= positive_threshold)
num_negative_edges_top_percent = np.sum(correlations <= negative_threshold)
print(f"Number of edges in the top {int(args.correlation_threshold*100)}% of positive correlations: {num_positive_edges_top_percent}")
print(f"Number of edges in the top {int(args.correlation_threshold*100)}% of negative correlations: {num_negative_edges_top_percent}")
# Identify high positive and negative correlations
high_positive_corr = spearman_corr_matrix >= positive_threshold
high_negative_corr = spearman_corr_matrix <= negative_threshold

high_corr_indices = np.where(high_positive_corr | high_negative_corr)

#Prepare the gene pairs and their correlation signs
gene_pairs = []
for i, j in zip(*high_corr_indices):
   if i < j:  # Ensure we don't include diagonal and duplicate pairs
       sign = "1" if high_positive_corr.iloc[i, j] else "-1"
       gene_pairs.append((df.columns[i], df.columns[j], sign))

gene_pairs_df = pd.DataFrame(gene_pairs)
# rename node to integers and save the mapping dataframe.
node_list = list(set(gene_pairs_df[0]) | set(gene_pairs_df[1]))
output_file = os.path.join(os.path.dirname(args.file_path), args.out_file_name + "_spearman.edgelist")
mapping = dict(zip(node_list, range(len(node_list))))
print (f'Number of genes: {len(mapping)}')
if args.reindex:
    gene_pairs_df[0] = gene_pairs_df[0].map(mapping)
    gene_pairs_df[1] = gene_pairs_df[1].map(mapping)
    pd.DataFrame(mapping, index=[1]).T.reset_index()[[1,'index']].to_csv(os.path.join(os.path.dirname(args.file_path), args.out_file_name + "_spearman_nodeID_mapping.tsv"), sep='\t', index=False,header=False)
#Output the high correlation gene pairs
gene_pairs_df.to_csv(output_file, sep='\t', index=False,header=False)

print('------------------- Finish -------------------')