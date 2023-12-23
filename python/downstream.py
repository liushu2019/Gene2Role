import numpy as np
from scipy.spatial import distance

def find_closest_nodes(embeddings, target_node, n):
    """
    Find the closest n nodes to a target node in a graph embedding.

    :param embeddings: Dictionary of node embeddings {node: embedding_vector}
    :param target_node: The target node for which to find closest nodes
    :param n: The number of closest nodes to find
    :return: List of tuples (node, distance) for the n closest nodes
    """
    if target_node not in embeddings:
        raise ValueError("Target node not found in embeddings")

    target_embedding = embeddings[target_node]
    distances = []

    for node, embedding in embeddings.items():
        if node != target_node:
            dist = distance.euclidean(target_embedding, embedding)
            distances.append((node, dist))

    # Sort by distance and return the closest n nodes
    distances.sort(key=lambda x: x[1])
    return distances[:n]