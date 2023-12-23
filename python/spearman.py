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
                    help="Threshold for filtering high correlations. (default:0.4)",
                    default=0.4)

args = parser.parse_args()

df = pd.read_csv(args.file_path, index_col=0).transpose()
# Calculate spearman correlation.
spearman_corr_matrix = df.corr(method='spearman')

# Flatten the correlation matrix and filter out self-correlations
correlations = spearman_corr_matrix.values.flatten()
correlations = correlations[~np.isnan(correlations)]
correlations = correlations[correlations != 1]

# Calculate histogram of correlations
#hist, bin_edges = np.histogram(correlations, bins='auto', density=True)

# Determine top 5% thresholds for positive and negative correlations

positive_threshold = np.percentile(correlations, 99)
negative_threshold = np.percentile(correlations, 0.1)
print(positive_threshold)
print(negative_threshold)
num_positive_edges_top_1_percent = np.sum(correlations >= positive_threshold)
num_negative_edges_top_1_percent = np.sum(correlations <= negative_threshold)

# Plotting
plt.hist(correlations, bins='auto', density=True, alpha=0.7, color='blue')
plt.axvline(x=positive_threshold, color='red', linestyle='dashed', linewidth=2)
plt.axvline(x=negative_threshold, color='green', linestyle='dashed', linewidth=2)
plt.title('Frequency Distribution of Correlations in ' + args.out_file_name)
plt.xlabel('Spearman Correlation Coefficient')
plt.ylabel('Frequency')

plot_file_name = os.path.join(os.path.dirname(args.file_path), args.out_file_name + "_correlation_histogram.png")
plt.savefig(plot_file_name)

print(f"Number of edges in the top 1% of positive correlations: {num_positive_edges_top_1_percent}")
print(f"Number of edges in the top 0.1% of negative correlations: {num_negative_edges_top_1_percent}")

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
output_file = os.path.join(os.path.dirname(args.file_path), args.out_file_name + "_spearman.edgelist")

#Output the high correlation gene pairs
gene_pairs_df.to_csv(output_file, sep='\t', index=False,header=False)

print('------------------- Finish -------------------')