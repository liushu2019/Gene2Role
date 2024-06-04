# -*- coding: utf-8 -*-

import numpy as np
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import Manager
from time import time
from collections import deque, defaultdict

from utils import *
from algorithm import *
from algorithm_distance import *

class Graph_complex():
	def __init__(self, g1, g2, is_directed, workers, untilLayer = None):

		logging.info(" - Converting graph to dict...")
		self.G_p = defaultdict(mylambda)
		self.G_n = defaultdict(mylambda)
		self.G_p.update(g1.gToDict())
		self.G_n.update(g2.gToDict())
		logging.info("Graph converted.")
		# self.p_max_degree = g1.maxDegree()
		# self.n_max_degree = g2.maxDegree()

		self.num_vertices = len(set(list(g1.nodes()) + list(g2.nodes())))
		self.num_edges = g1.number_of_edges() + g2.number_of_edges()
		self.is_directed = is_directed
		self.workers = workers
		self.calcUntilLayer = untilLayer
		logging.info('Graph - Number of vertices: {}'.format(self.num_vertices))
		logging.info('Graph - Number of edges: {}'.format(self.num_edges))

	def preprocess_neighbors_with_bfs(self):
		# exec_bfs_complex(self.G_p,self.G_n,self.workers,self.calcUntilLayer)
		with ProcessPoolExecutor(max_workers=self.workers) as executor:
			job = executor.submit(exec_bfs_complex,self.G_p,self.G_n,self.workers,self.calcUntilLayer)

			job.result()

		return

	def preprocess_neighbors_with_bfs_compact(self): # TODO
		# exec_bfs_compact_complex(self.G_p,self.G_n,self.workers,self.calcUntilLayer)
		with ProcessPoolExecutor(max_workers=self.workers) as executor:
			job = executor.submit(exec_bfs_compact_complex,self.G_p,self.G_n,self.workers,self.calcUntilLayer)

			job.result()

		return

	def create_vectors(self):
		logging.info("Creating degree vectors...")
		degrees = {}
		degrees_sorted = set()
		G = self.G
		for v in G.keys():
			degree = len(G[v])
			degrees_sorted.add(degree)
			if(degree not in degrees):
				degrees[degree] = {}
				degrees[degree]['vertices'] = deque() 
			degrees[degree]['vertices'].append(v)
		degrees_sorted = np.array(list(degrees_sorted),dtype='int')
		degrees_sorted = np.sort(degrees_sorted)

		l = len(degrees_sorted)
		for index, degree in enumerate(degrees_sorted):
			if(index > 0):
				degrees[degree]['before'] = degrees_sorted[index - 1]
			if(index < (l - 1)):
				degrees[degree]['after'] = degrees_sorted[index + 1]
		logging.info("Degree vectors created.")
		logging.info("Saving degree vectors...")
		saveVariableOnDisk(degrees,'degrees_vector')

	def create_vectors_complex(self):
		logging.info("Creating degree vectors...")
		degrees_p = {} # matrix // TODO
		degrees_n = {} # matrix //TODO
		degrees = {} # matrix
		degrees_sorted_p = set()
		degrees_sorted_p.add(0)
		degrees_sorted_n = set()
		degrees_sorted_n.add(0)
		G = self.G_p
		for v in G.keys():
			degree = len(G[v])
			degrees_sorted_p.add(degree)
			if(degree not in degrees_p):
				degrees_p[degree] = {}
				degrees_p[degree]['vertices'] = deque() 
			degrees_p[degree]['vertices'].append(v)
		degrees_sorted_p = np.array(list(degrees_sorted_p),dtype='int')
		degrees_sorted_p = np.sort(degrees_sorted_p)

		G = self.G_n
		for v in G.keys():
			degree = len(G[v])
			degrees_sorted_n.add(degree)
			if(degree not in degrees_n):
				degrees_n[degree] = {}
				degrees_n[degree]['vertices'] = deque() 
			degrees_n[degree]['vertices'].append(v)
		degrees_sorted_n = np.array(list(degrees_sorted_n),dtype='int')
		degrees_sorted_n = np.sort(degrees_sorted_n)

		degrees_p[0] = dict(vertices=deque(set(self.G_n.keys()) - set(self.G_p.keys())) )
		degrees_n[0] = dict(vertices=deque(set(self.G_p.keys()) - set(self.G_n.keys())) )

		l_p = len(degrees_sorted_p)
		l_n = len(degrees_sorted_n)

		for index_p, degree_p in enumerate(degrees_sorted_p):
			for index_n, degree_n in enumerate(degrees_sorted_n):
				degree = complex(degree_p, degree_n)
				if(degree not in degrees):
					degrees[degree] = {}
				degrees[degree]['vertices'] = list(set(degrees_n[degree_n]['vertices']) & set(degrees_p[degree_p]['vertices']))
				# if(index_p > 0):
				# 	degrees[degree]['before_p'] = degrees_sorted_p[index_p - 1]
				# if(index_p < (l_p - 1)):
				# 	degrees[degree]['after_p'] = degrees_sorted_p[index_p + 1]
				# if(index_n > 0):
				# 	degrees[degree]['before_n'] = degrees_sorted_n[index_n - 1]
				# if(index_n < (l_n - 1)):
				# 	degrees[degree]['after_n'] = degrees_sorted_n[index_n + 1]
		degrees_sorted_n = list(degrees_sorted_n)
		degrees_sorted_p = list(degrees_sorted_p)
		logging.info("Degree vectors created.")
		logging.info("Saving degree vectors...")
		saveVariableOnDisk(degrees,'degrees_vector')
		saveVariableOnDisk(degrees_sorted_n,'degrees_vector_negativeList')
		saveVariableOnDisk(degrees_sorted_p,'degrees_vector_positiveList')

	def calc_distances_all_vertices(self,compactDegree = False, scale_free = False):
		# maxA = np.log(np.sqrt(self.n_max_degree**2+self.p_max_degree**2)+1)
		# set_maxA(maxA)
		# logging.info("Using maxA: {}".format(maxA))
		logging.info("Using compactDegree: {}".format(compactDegree))
		if(self.calcUntilLayer):
			logging.info("Calculations until layer: {}".format(self.calcUntilLayer))

		futures = {}

		count_calc = 0

		vertices = list(reversed(sorted(list(set(list(self.G_p.keys()) + list(self.G_n.keys()))))))
		# list(reversed(sorted(list(self.G.keys()))))

		if(compactDegree):
			logging.info("Recovering compactDegreeList from disk...")
			degreeList = restoreVariableFromDisk('compactDegreeList')
		else:
			logging.info("Recovering degreeList from disk...")
			degreeList = restoreVariableFromDisk('degreeList')

		parts = self.workers
		chunks = partition(vertices,parts)

		t0 = time()
#debug
		# part = 1
		# for c in chunks:
		# 	logging.info("Executing part {}...".format(part))
		# 	list_v = []
		# 	for v in c:
		# 		list_v.append([vd for vd in degreeList.keys() if vd > v])
		# 	calc_distances_all_complex( c, list_v, degreeList,part, compactDegree = compactDegree)


		with ProcessPoolExecutor(max_workers = self.workers) as executor:

			part = 1
			for c in chunks:
				logging.info("Executing part {}...".format(part)+str(c))
				list_v = []
				for v in c:
					list_v.append([vd for vd in degreeList.keys() if vd > v])
				job = executor.submit(calc_distances_all_complex, c, list_v, degreeList,part, compactDegree = compactDegree, scale_free=scale_free)
				futures[job] = part
				part += 1


			logging.info("Receiving results...")

			for job in as_completed(futures):
				job.result()
				r = futures[job]
				logging.info("Part {} Completed.".format(r))
# end debug
		logging.info('Distances calculated.')
		t1 = time()
		logging.info('Time : {}m'.format((t1-t0)/60))

		return


	def calc_distances(self, compactDegree = False):

		logging.info("Using compactDegree: {}".format(compactDegree))
		if(self.calcUntilLayer):
			logging.info("Calculations until layer: {}".format(self.calcUntilLayer))

		futures = {}
		#distances = {}

		count_calc = 0

		G = self.G
		vertices = list(G.keys())

		parts = self.workers
		chunks = partition(vertices,parts)

		with ProcessPoolExecutor(max_workers = 1) as executor:

			logging.info("Split degree List...")
			part = 1
			for c in chunks:
				job = executor.submit(splitDegreeList,part,c,G,compactDegree)
				job.result()
				logging.info("degreeList {} completed.".format(part))
				part += 1


		with ProcessPoolExecutor(max_workers = self.workers) as executor:

			part = 1
			for c in chunks:
				logging.info("Executing part {}...".format(part))
				job = executor.submit(calc_distances, part, compactDegree = compactDegree)
				futures[job] = part
				part += 1

			logging.info("Receiving results...")
			for job in as_completed(futures):
				job.result()
				r = futures[job]
				logging.info("Part {} completed.".format(r))


		return

	def calc_distances_complex(self, compactDegree = False, scale_free = False):
    
		logging.info("Using compactDegree: {}".format(compactDegree))
		if(self.calcUntilLayer):
			logging.info("Calculations until layer: {}".format(self.calcUntilLayer))

		futures = {}
		#distances = {}

		count_calc = 0

		# G = self.G
		vertices = list(set(list(self.G_p.keys()) + list(self.G_n.keys())))

		parts = self.workers
		chunks = partition(vertices,parts)

		# part = 1
		# for c in chunks:
		# 	# print ('debug part = '+str(part)+ ' chunks'+str(chunks))
		# 	splitDegreeList_complex(part,c,self.G_p,self.G_n,compactDegree)
		# 	# print ('debug part = '+str(part))
		# 	part += 1
		with ProcessPoolExecutor(max_workers = 1) as executor:

			logging.info("Split degree List...")
			part = 1
			for c in chunks:
				job = executor.submit(splitDegreeList_complex,part,c,self.G_p,self.G_n,compactDegree)
				job.result()
				logging.info("degreeList {} completed.".format(part))
				part += 1

		# part = 1
		
		# for c in chunks:
		# 	logging.info("Executing part {}...".format(part))
		# 	calc_distances_complex( part, compactDegree = compactDegree)
		# 	part += 1

		with ProcessPoolExecutor(max_workers = self.workers) as executor:

			part = 1
			for c in chunks:
				logging.info("Executing part {}...".format(part))
				job = executor.submit(calc_distances_complex, part, compactDegree = compactDegree, scale_free = scale_free)
				futures[job] = part
				part += 1

			logging.info("Receiving results...")
			for job in as_completed(futures):
				job.result()
				r = futures[job]
				logging.info("Part {} completed.".format(r))


		return

	def consolide_distances(self):

		distances = {}

		parts = self.workers
		for part in range(1,parts + 1):
			d = restoreVariableFromDisk('distances-'+str(part))
			preprocess_consolides_distances(distances)
			distances.update(d)


		preprocess_consolides_distances(distances)
		saveVariableOnDisk(distances,'distances')


	def create_distances_network(self):

		with ProcessPoolExecutor(max_workers=1) as executor:
			job = executor.submit(generate_distances_network,self.workers)

			job.result()

		return

	def preprocess_parameters_random_walk(self):

		with ProcessPoolExecutor(max_workers=1) as executor:
			job = executor.submit(generate_parameters_random_walk,self.workers)

			job.result()

		return


	def simulate_walks(self,num_walks,walk_length):

		# for large graphs, it is serially executed, because of memory use.
		if(self.num_vertices > 500000):

			with ProcessPoolExecutor(max_workers=1) as executor:
				job = executor.submit(generate_random_walks_large_graphs,num_walks,walk_length,self.workers,list(set(list(self.G_p.keys()) + list(self.G_n.keys()))))

				job.result()

		else:

			with ProcessPoolExecutor(max_workers=1) as executor:
				job = executor.submit(generate_random_walks,num_walks,walk_length,self.workers,list(set(list(self.G_p.keys()) + list(self.G_n.keys()))))
				job.result()

		return	
