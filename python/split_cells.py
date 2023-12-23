import pandas as pd
import argparse
import os
import sys

parser = argparse.ArgumentParser(description="Split gene-cell matrix based on cell types.")
parser.add_argument("count_matrix", type=str, help="CSV file for gene-cell matrix.")
parser.add_argument("cell_metadata", type=str, help="CSV file with cell type information.")

args = parser.parse_args()

# Read the gene-cell count matrix and cell type information
gene_cell_df = pd.read_csv(args.count_matrix, index_col=0)
cell_types_df = pd.read_csv(args.cell_metadata)

# Get directory of the input file to save output files in the same directory
output_directory = os.path.dirname(args.count_matrix)

output_directory = os.path.join(output_directory, 'splitMatrix')
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

# Check if 'celltype' column exists in cell_types_df
if 'celltype' not in cell_types_df.columns:
    print("Error: No 'celltype' column found in cell types file.")
    sys.exit(1)
else:
	cell_types_df = cell_types_df.reset_index()[['index', 'celltype']]

# Merge the datasets on the cell name
merged_df = pd.merge(gene_cell_df.transpose(), 
										cell_types_df, 
										left_on=gene_cell_df.columns, 
										right_on='index')
gene_names = gene_cell_df.index

index_tracker = pd.DataFrame(index=gene_names)

# Loop through each cell type and process the data
print(" ")
print("---------- cell type info ----------")
for i, cell_type in enumerate(merged_df['celltype'].unique()):
    
    # Filter data for the specific cell type
    cell_type_df = merged_df[merged_df['celltype'] == cell_type]
    cell_type_df = cell_type_df.drop(columns=['celltype'])

    # Print the cell type and its count
    print(f"{cell_type}, Cell No.: {len(cell_type_df)}")
    # Transpose back to the original format (genes as rows, cells as columns)
    
    cell_type_df = cell_type_df.set_index('index')

    new_gene_indices = range(i * len(gene_names) + 1, (i + 1) * len(gene_names) + 1)
    cell_type_df.columns = new_gene_indices
    cell_type_df = cell_type_df.transpose()

    index_tracker[cell_type] = new_gene_indices

    # Save to a new CSV file
    cell_type = cell_type.replace(" ", "_")
    output_file = os.path.join(output_directory, f'{cell_type}.csv')
    cell_type_df.to_csv(output_file)

index_tracker.to_csv(os.path.join(output_directory, 'index_tracker.tsv'), sep='\t')

print("---------- count matrix split complete. ----------")