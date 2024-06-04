# -*- coding: utf-8 -*-

import argparse, logging
from gensim.models import Word2Vec
from gensim.models.word2vec import LineSentence
from time import time
import os
import signeds2v
import graph
import warnings

# This will suppress all warnings
warnings.filterwarnings("ignore")
logging.basicConfig(filename='signeds2v.log',filemode='w',level=logging.DEBUG,format='%(asctime)s %(message)s')

def parse_args():
	'''
	Parses the signeds2v arguments.
	'''
	parser = argparse.ArgumentParser(description="Run signeds2v.")

	parser.add_argument('--input', nargs='?', default='karate-mirrored.edgelist',
	                    help='Input graph path')
	parser.add_argument('--output', nargs='?', default=None,
	                    help='Output emb path, if Not given, follow input file name')

	parser.add_argument('--dimensions', type=int, default=128,
	                    help='Number of dimensions. Default is 128.')

	parser.add_argument('--walk-length', type=int, default=80,
	                    help='Length of walk per source. Default is 80.')

	parser.add_argument('--num-walks', type=int, default=20,
	                    help='Number of walks per source. Default is 10.')

	parser.add_argument('--window-size', type=int, default=10,
                    	help='Context size for optimization. Default is 10.')

	parser.add_argument('--until-layer', type=int, default=6,
                    	help='Calculation until the layer.')

	parser.add_argument('--iter', default=5, type=int,
                      help='Number of epochs in SGD')

	parser.add_argument('--workers', type=int, default=8,
	                    help='Number of parallel workers. Default is 8.')

	parser.add_argument('--weighted', dest='weighted', action='store_true',
	                    help='Boolean specifying (un)weighted. Default is unweighted.')
	parser.add_argument('--unweighted', dest='unweighted', action='store_false')
	parser.set_defaults(weighted=False)

	parser.add_argument('--directed', dest='directed', action='store_true',
	                    help='Graph is (un)directed. Default is undirected.')
	parser.add_argument('--undirected', dest='undirected', action='store_false')
	parser.set_defaults(directed=False)

	parser.add_argument('--OPT1', action='store_true',
                      help='optimization 1')
	parser.add_argument('--OPT2', action='store_true',
                      help='optimization 2')
	parser.add_argument('--OPT3', action='store_true',
                      help='optimization 3')
	parser.add_argument('--scalefree', action='store_true',
                      help='scale free flag')
	return parser.parse_args()


def read_graph_signed():
	'''
	Reads the input signed network.
	'''
	logging.info(" - Loading signed graph...")
	Gp, Gm = graph.load_edgelist_signed(args.input,undirected=True)
	logging.info(" - Signed Graph loaded.")
	return Gp, Gm

def learn_embeddings(basename, output):
	'''
	Learn embeddings by optimizing the Skipgram objective using SGD.
	'''
	logging.info("Initializing creation of the representations...")
	walks = LineSentence('random_walks.txt')
	model = Word2Vec(walks, vector_size=args.dimensions, window=args.window_size, min_count=0, hs=0, sg=1,negative=5, ns_exponent=0.75, workers=args.workers, epochs=args.iter)
	if (output is not None):
		basename = output
		directory_path = os.path.dirname(basename)
		if not os.path.exists(directory_path):
			try:
				os.makedirs(directory_path)
				print(f"Directory '{directory_path}' created.")
			except OSError as e:
				print(f"Error creating directory '{directory_path}': {e}")
				return
	else:
		basename = "emb/{}.emb".format(basename)
	model.wv.save_word2vec_format(basename)
	logging.info("Representations created.")
	return

def exec_signeds2v_complex(args):
	'''
	Pipeline for representational learning for all nodes in a graph.
	'''
	if(args.OPT3):
		until_layer = args.until_layer
	else:
		until_layer = None

	Gp, Gm = read_graph_signed()
	print('read complete')
	G = signeds2v.Graph_complex(Gp, Gm, args.directed, args.workers, untilLayer = until_layer)
	print('signeds2v graph complete')
	if(args.OPT1):
		G.preprocess_neighbors_with_bfs_compact()
	else:
		G.preprocess_neighbors_with_bfs()

	if(args.OPT2):
		G.create_vectors_complex()
		G.calc_distances_complex(compactDegree = args.OPT1, scale_free = args.scalefree)
	else:
		G.calc_distances_all_vertices(compactDegree = args.OPT1, scale_free = args.scalefree)
	G.create_distances_network()
	G.preprocess_parameters_random_walk()
	print('multi-layer network generation complete')
	G.simulate_walks(args.num_walks, args.walk_length)

	print('random walk complete')
	return G

def main(args):
	print('Process start')
	G = exec_signeds2v_complex(args)
	print('complete network generations.')
	learn_embeddings(args.input.split("/")[1].split(".")[0], args.output)

if __name__ == "__main__":
	args = parse_args()
	print (args)
	main(args)