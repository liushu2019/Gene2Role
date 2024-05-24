# Gene2role
Hi! This repository contains the code for our paper “Gene2role: A Role-Based Gene Embedding Method for Comparative Analysis of Signed Gene Regulatory Networks.” You can read the paper on BioRxiv: https://www.biorxiv.org/content/10.1101/2024.05.18.594807v1.abstract.

The repository is organized into two main sections:

	1.	Tutorial for Gene2role
	2.	Reproduction of the Paper’s Results

# Tutorial


## folder structure:

- pipeline.py
- tools
  - [SignedS2V](https://github.com/liushu2019/SignedS2V)
- codes
  - split_cells.py
  - spearman.py
  - eeisp.py

## command line
```
python pipeline.py TaskMode CellType EmbeddingMode input [**]

TaskMode.      1: run SignedS2V for an edgelist file. 
               2: run spearman and SignedS2V from gene X cell count matrix. 
               3: run eeisp and SignedS2V from gene X cell count matrix.
CellType.      1: single cell-type. 
               2: multiple cell-type.
EmbeddingMode. 1: single network embedding. 
               2: multiple network embedding, only work if the previous argument is 2.
input          Input file, either a gene X cell matrix for TaskMode 2 and 3, or edgelist for TaskMode 1.

Other arguments, check "python pipeline.py --help"
```
## Exp:
```
python pipeline.py 3 2 1 data/singleCell/test_rna/count.csv --project TEST321 --OPT1 --OPT2 --OPT3 --workers 61 --until_layer 1 --threCDI 0.0001 --threEEI 0.0001 --cell_metadata data/singleCell/test_rna/metadata.csv 
```
## Parameter combination
```
TBD
```

# Reproduction 

All the processed data generated in the paper can be downloaded from https://figshare.com/articles/dataset/data/25852915. 

The basic information about the data we used is as follows:
![TableS1](TableS1.jpg)

### Gene2role hyperparameters
We used Gene2role to generate embedding by the following hyperparameters using the `pipeline.py` mention above.

### Downstream analysis
In **result1**, we primarily analyzed gene embeddings from a single GRN generated using four different data types. The analysis processes for simulated and curated gene regulatory networks can be found in `\reproduce\result1\simulated_and_curated_networks.R`, and for single-cell RNA-seq and single-cell multi-omics data in `\reproduce\result1\single_cell_networks.R`.

In **result2**, we analyzed the differentially topological genes from pair cell types. you can find the codes in '`\reproduce\result2\`

In **result3**, we analyzed the differentially topological genes from multiple cell types. you can find the codes in '`\reproduce\result3\`

In **result4**, we analyszed the gene module stability between two cell states. you can find the codes in `\reproduce\result4\`

# Contact
Feel free to reach Xin Zeng (wstxinzeng@gmail) and Shu Liu (Shu.liu.eq@gmail.com)to request files and more details from the analysis process!
