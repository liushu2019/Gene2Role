# SignedS2V
Code for paper titled as "SignedS2V: structural embedding method for signed networks" in Complex Networks 2022

SignedS2V learns the embeddings for nodes based on the structural features.

# Installation
Packages needed: 

	pip install futures
	pip install fastdtw
	pip install gensim

# Usage
## Input
Edgelist of a signed network, 3 columns (node1 node2 sign). sign={-1,1}

## Parameters

	'--input': Input graph path
	'--output': Output emb path
	'--dimensions': Number of dimensions. Default is 128.
	'--walk-length': Length of walk per source. Default is 80.
	'--num-walks': Number of walks per source. Default is 10.
	'--window-size': Context size for optimization. Default is 10.
	'--until-layer':Calculation until the layer.
	'--iter': Number of epochs in SGD
	'--workers': Number of parallel workers. Default is 8.
	'--OPT1': optimization 1
	'--OPT2': optimization 2
	'--OPT3': optimization 3
	'--scalefree': scale free flag
  
## Command example
  	python src/main.py --input graph/mirrored_karate_sign.edgelist --num-walks 20 --walk-length 80 --window-size 5 --dimensions 2 --until-layer 5 --workers 8 --OPT1 --OPT2 --OPT3 --scalefree

# Miscellaneous
  Feel free to send email to liu@torilab.net for any questions about the code or the paper.
  
  Note that this framework only work for undirected signed network.
  
  We would like to thank Leonardo F. R. Ribeiro et al., authors of struc2vec, for providing their code.
