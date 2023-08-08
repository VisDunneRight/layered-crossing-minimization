import src.read_data
import src.vis
from src import graph
import random
import os


def random_layered_graph(k, n, d):
	"""
	:param k: number of layers
	:param n: number of nodes per layer
	:param d: average edge density of the resultant graph
	:return: LayeredGraph object g
	"""

	if d >= 0.5:
		raise Exception("wrong implementation for high density. d should be in [0, 0.5)")

	g = graph.LayeredGraph()
	n_edges_per_layer = round(d * (n ** 2))

	for i in range(k):  # add nodes
		for j in range(n):
			g.add_node(i + 1)

	for i in range(1, k):  # randomly, uniformly select edges
		n_added = 0
		while n_added < n_edges_per_layer:
			n1 = random.randint(0, n-1) + ((i - 1) * n)
			n2 = random.randint(0, n-1) + (i * n)
			if (n1, n2) not in g.edge_names:
				g.add_edge(n1, n2)
				n_added += 1
	return g


def generate_gange_dataset(seed=None):

	if seed is not None:
		random.seed(seed)

	if "random graphs" not in os.listdir(".."):
		os.mkdir("../random graphs")
	if "gange" not in os.listdir("../random graphs"):
		os.mkdir("../random graphs/gange")

	for k in range(3, 11):
		for n in range(7, 11 if k < 5 else (10 if k < 8 else 9)):
			if f"g{k}_{n}" not in os.listdir("../random graphs/gange"):
				os.mkdir(f"../random graphs/gange/g{k}_{n}")
			for i in range(10):
				ng = random_layered_graph(k, n, 0.2)
				ng.write_out(f"../random graphs/gange/g{k}_{n}/graph{i}.lgbin")


if __name__ == '__main__':
	# generate_gange_dataset(seed=22)
	g = src.read_data.read("../random graphs/gange/g3_9/graph2.lgbin")
	src.vis.draw_graph(g, "rand", gravity=True)
