# g2r_work


# folder structure:

-pipeline.py
-tools
 -[SignedS2V](https://github.com/liushu2019/SignedS2V)
-codes
 -split_cells.py
 -spearman.py
 -eeisp.py

# command line
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
# Exp:
```
python pipeline.py 3 2 1 data/singleCell/test_rna/count.csv --project TEST321 --OPT1 --OPT2 --OPT3 --workers 61 --until_layer 1 --threCDI 0.0001 --threEEI 0.0001 --cell_metadata data/singleCell/test_rna/metadata.csv 
```
# [Memo] Test progress records:
```
1   1  1 O
2   1  1 O
2   2  1 
2   2  2
3   1  1 O
3   2  1 O
3   2  2 i
```
