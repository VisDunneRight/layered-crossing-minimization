import time

import src.read_data
import src.vis
from src import graph
import random
import os


def true_random_connected_layered_graph(k, n, d):
	"""
	:param k: number of layers
	:param n: number of nodes per layer
	:param d: average edge density of the resultant graph
	:return: LayeredGraph object g

	Uniformly random connected layered graph. Will retry until the sampled graph is connected.
	"""

	n_edges_per_layer = round(d * (n ** 2))
	assert n_edges_per_layer * (k - 1) >= n * k, "graph will not be connected"

	flip_edges = d > 0.5
	d = 1 - d if flip_edges else d

	while True:
		g = graph.LayeredGraph()
		n_edges_per_layer = round(d * (n ** 2))

		for i in range(k):  # add nodes
			for j in range(n):
				g.add_node(i)

		for i in range(k - 1):  # randomly, uniformly select edges
			n_added = 0
			while n_added < n_edges_per_layer:
				n1 = random.randint(0, n-1) + (i * n)
				n2 = random.randint(0, n-1) + ((i + 1) * n)
				if (n1, n2) not in g.edge_ids:
					g.add_edge(n1, n2)
					n_added += 1
		if flip_edges:
			edges = set(g.edge_ids.keys())
			g.edges = []
			g.edge_ids = {}
			for i in range(1, k):
				for n1 in range(i * n, n + (i * n)):
					for n2 in range((i + 1) * n, (i + 2) * n):
						if (n1, n2) not in edges:
							g.add_edge(n1, n2)

		if g.is_connected():
			return g


def random_layered_graph_connect_help(k, n, d):
	"""
	:param k: number of layers
	:param n: number of nodes per layer
	:param d: average edge density of the resultant graph
	:return: LayeredGraph object g

	Random connected layered graph, but will select from the set of unconnected nodes if necessary to try to ensure graph is connected.
	"""

	n_edges_per_layer = round(d * (n ** 2))
	assert n_edges_per_layer * (k - 1) >= n * k, "graph will not be connected"

	flip_edges = d > 0.5
	d = 1 - d if flip_edges else d

	while True:
		g = graph.LayeredGraph()
		n_edges_per_layer = round(d * (n ** 2))

		for i in range(k):  # add nodes
			for j in range(n):
				g.add_node(i)

		not_seen = set(range(n * k))
		for i in range(k - 1):  # randomly, uniformly select edges
			n_added = 0
			not_seen_l1 = set((x for x in not_seen if i * n <= x < n + (i * n)))
			if i == k - 1:
				not_seen_l2 = set((x for x in not_seen if (i + 1) * n <= x < (i + 2) * n))
			while n_added < n_edges_per_layer:
				n1 = random.randint(0, n-1) + (i * n)
				while len(not_seen_l1) == n_edges_per_layer - n_added and n1 not in not_seen_l1:
					n1 = random.randint(0, n - 1) + (i * n)
				n2 = random.randint(0, n-1) + ((i + 1) * n)
				if i == k - 1:
					while len(not_seen_l2) == n_edges_per_layer - n_added and n2 not in not_seen_l2:
						n2 = random.randint(0, n-1) + ((i + 1) * n)
				if (n1, n2) not in g.edge_ids:
					g.add_edge(n1, n2)
					n_added += 1
					if n1 in not_seen:
						not_seen.remove(n1)
						not_seen_l1.remove(n1)
					if n2 in not_seen:
						not_seen.remove(n2)
						if i == k - 1:
							not_seen_l2.remove(n2)
		if flip_edges:
			edges = set(g.edge_ids.keys())
			g.edges = []
			g.edge_ids = {}
			for i in range(1, k):
				for n1 in range((i - 1) * n, n + ((i - 1) * n)):
					for n2 in range(i * n, n + (i * n)):
						if (n1, n2) not in edges:
							g.add_edge(n1, n2)

		if g.is_connected():
			return g
		else:
			print("fail")


def random_layered_graph_connect_help_edgecount(k, n, n_edges):
	"""
	:param k: number of layers
	:param n: number of nodes per layer
	:param n_edges: number of edges in the resultant graph
	:return: LayeredGraph object g

	Randomly samples edges over the full network instead of keeping the edge count constant across layers
	as in the above methods.
	"""

	assert n_edges >= n * k, "graph will not be connected"

	max_edges = (k - 1) * (n * n)

	assert n_edges <= max_edges, f"Max edges for this graph size is {max_edges}"

	flip_edges = True if n_edges > max_edges // 2 else False
	if flip_edges:
		n_edges = max_edges - n_edges

	while True:
		g = graph.LayeredGraph()

		for i in range(k):  # add nodes
			for j in range(n):
				g.add_node(i)

		not_seen = set(range(n * k))
		seen_edges = set()
		n_added = 0
		while n_added < n_edges:
			n1_l = random.randint(0, g.n_layers - 2)
			n1 = random.randint(0, n - 1) + (n1_l * n)
			# if len(not_seen) == n_edges - n_added and not flip_edges:
			# 	n1 = random.choice(sorted(not_seen))
			# 	n1_l = g[n1].layer
			n2_l = n1_l + 1 if n1_l != g.n_layers - 1 else n1_l - 1
			n2 = random.randint(0, n - 1) + (n2_l * n)
			if (n1, n2) not in seen_edges and (n2, n1) not in seen_edges:
				g.add_edge(n1, n2)
				if n1 in not_seen:
					not_seen.remove(n1)
				if n2 in not_seen:
					not_seen.remove(n2)
				seen_edges.add((n1, n2))
				n_added += 1
		adj = g.get_adj_list()
		for unseen in not_seen:
			unconnected = True
			while unconnected:
				ot_nd_in_layer = random.choice()
		if flip_edges:
			pre_flip = set(g.edge_ids.keys())
			print(pre_flip)
			g.edges = []
			g.edge_ids = {}
			g.invalidate_data()
			for i in range(k - 1):
				for n1 in range(i * n, n + (i * n)):
					for n2 in range((i + 1) * n, ((i + 2) * n)):
						if (n1, n2) not in pre_flip:
							print(n1, n2)
							g.add_edge(n1, n2)

		if g.is_connected():
			return g
		else:
			print("fail")


def random_layered_graph_edgecount_drop_unconnected(k, n, n_edges):
	"""
	:param k: number of layers
	:param n: number of nodes per layer
	:param n_edges: number of edges in the resultant graph
	:return: LayeredGraph object g

	Randomly samples edges over the entire network. Any nodes left unconnected are removed
	"""

	max_edges = (k - 1) * (n * n)

	assert n_edges <= max_edges, f"Max edges for this graph size is {max_edges}"

	flip_edges = True if n_edges > max_edges // 2 else False
	if flip_edges:
		n_edges = max_edges - n_edges

	while True:
		g = graph.LayeredGraph()

		for i in range(k):  # add nodes
			for j in range(n):
				g.add_node(i)

		not_seen = set(range(n * k))
		seen_edges = set()
		n_added = 0
		while n_added < n_edges:
			n1_l = random.randint(0, g.n_layers - 2)
			n1 = random.randint(0, n - 1) + (n1_l * n)
			n2_l = n1_l + 1
			n2 = random.randint(0, n - 1) + (n2_l * n)
			if (n1, n2) not in seen_edges and (n2, n1) not in seen_edges:
				g.add_edge(n1, n2)
				if n1 in not_seen:
					not_seen.remove(n1)
				if n2 in not_seen:
					not_seen.remove(n2)
				seen_edges.add((n1, n2))
				n_added += 1
		if flip_edges:
			pre_flip = set(g.edge_ids.keys())
			print(pre_flip)
			g.edges = []
			g.edge_ids = {}
			g.invalidate_data()
			for i in range(k - 1):
				for n1 in range(i * n, n + (i * n)):
					for n2 in range((i + 1) * n, ((i + 2) * n)):
						if (n1, n2) not in pre_flip:
							g.add_edge(n1, n2)
		else:
			gp = graph.LayeredGraph()
			nd_map = []
			for nd in g.nodes:
				if nd.id not in not_seen:
					x = gp.add_node(nd.layer)
					nd_map.append(x.id)
				else:
					nd_map.append(0)
			for ed in seen_edges:
				gp.add_edge(nd_map[ed[0]], nd_map[ed[1]])
			g = gp

		if g.is_connected():
			return g
		else:
			print("fail")


def random_layered_graph_edgecount_difflayers(n_edges, layercounts: list):
	"""
	:param n_edges: number of edges in the resultant graph
	:param layercounts: list of #nodes for each layer
	:return: LayeredGraph object g

	Randomly samples edges over the entire network. Any nodes left unconnected are removed
	"""

	max_edges = sum((layercounts[i] * layercounts[i + 1] for i in range(len(layercounts) - 1)))
	n = sum(layercounts)

	assert n_edges <= max_edges, f"Max edges for this graph size is {max_edges}"

	flip_edges = True if n_edges > max_edges // 2 else False
	if flip_edges:
		n_edges = max_edges - n_edges

	while True:
		g = graph.LayeredGraph()

		for i, nl in enumerate(layercounts):  # add nodes
			for _ in range(nl):
				g.add_node(i)

		not_seen = set(range(len(g.nodes)))
		seen_edges = set()
		n_added = 0
		while n_added < n_edges:
			n1 = random.randint(0, n - 1)
			n1_l = g[n1].layer
			n2_l = n1_l + random.choice([-1, 1])
			if n2_l < 0 or n2_l >= g.n_layers:
				continue
			n2 = random.choice(g.layers[n2_l]).id
			if (n1, n2) not in seen_edges and (n2, n1) not in seen_edges:
				g.add_edge(n1, n2)
				if n1 in not_seen:
					not_seen.remove(n1)
				if n2 in not_seen:
					not_seen.remove(n2)
				seen_edges.add((n1, n2))
				n_added += 1
		if flip_edges:
			pre_flip = set(g.edge_ids.keys())
			print(pre_flip)
			g.edges = []
			g.edge_ids = {}
			g.invalidate_data()
			for i in range(g.n_layers - 1):
				for n1 in g.layers[i]:
					for n2 in g.layers[i + 1]:
						if (n1.id, n2.id) not in pre_flip:
							g.add_edge(n1.id, n2.id)
		else:
			# rebuild graph without disconnected nodes/2-cliques
			# print([len(g.get_edges_by_layer()[ls]) for ls in range(g.n_layers - 1)])
			gp = graph.LayeredGraph()
			nd_map = []
			adj = g.get_adj_list()
			for nd in g.nodes:
				if nd.id not in not_seen and not (len(adj[nd.id]) == 1 and len(adj[adj[nd.id][0]]) == 1):
					x = gp.add_node(nd.layer)
					nd_map.append(x.id)
				else:
					nd_map.append(-1)
			for ed in seen_edges:
				if nd_map[ed[0]] != -1 and nd_map[ed[1]] != -1:
					gp.add_edge(nd_map[ed[0]], nd_map[ed[1]])
			g = gp
			# print([len(g.get_edges_by_layer()[ls]) for ls in range(g.n_layers - 1)])

		if g.is_connected():
			return g
		else:
			print("fail")


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
				ng = true_random_connected_layered_graph(k, n, 0.2)
				ng.write_out(f"../random graphs/gange/g{k}_{n}/graph{i}.lgbin")


def generate_random_density_set(seed=None):
	if seed is not None:
		random.seed(seed)

	if "random graphs" not in os.listdir(".."):
		os.mkdir("../random graphs")
	if "density_exp" not in os.listdir("../random graphs"):
		os.mkdir("../random graphs/density_exp")

	for d in range(14, 51, 2):
		if f"d{d}" not in os.listdir("../random graphs/density_exp"):
			os.mkdir(f"../random graphs/density_exp/d{d}")
		for i in range(10):
			ng = true_random_connected_layered_graph(5, 10, d / 100)
			print(f"d={d} graph {i}")
			ng.write_out(f"../random graphs/density_exp/d{d}/graph{i}.lgbin")


def generate_random_fixed_density_set(seed=None):
	if seed is not None:
		random.seed(seed)

	if "random graphs" not in os.listdir(".."):
		os.mkdir("../random graphs")
	if "fixed_density_exp" not in os.listdir("../random graphs"):
		os.mkdir("../random graphs/fixed_density_exp")

	for k in range(3, 21):
		if f"k{k}" not in os.listdir("../random graphs/fixed_density_exp"):
			os.mkdir(f"../random graphs/fixed_density_exp/k{k}")
		for i in range(10):
			ng = random_layered_graph_connect_help(k, 10, 0.15)
			print(f"k={k} graph {i}")
			ng.write_out(f"../random graphs/fixed_density_exp/k{k}/graph{i}.lgbin")


def generate_extended_matuszewski_datsets(seed=None):
	if seed is not None:
		random.seed(seed)
	if "random graphs" not in os.listdir(".."):
		os.mkdir("../random graphs")
	if "matuszewski" not in os.listdir("../random graphs"):
		os.mkdir("../random graphs/matuszewski")
		os.mkdir("../random graphs/matuszewski/5_by_n")
		os.mkdir("../random graphs/matuszewski/k_by_10")
		os.mkdir("../random graphs/matuszewski/10_by_10_density")

	for n in range(10, 101, 10):
		if f"n{n}" not in os.listdir("../random graphs/matuszewski/5_by_n"):
			os.mkdir(f"../random graphs/matuszewski/5_by_n/n{n}")
			for i in range(100):
				ng = random_layered_graph_connect_help_edgecount(5, n, 8 * n)
				print(f"n={n} graph {i}")
				ng.write_out(f"../random graphs/matuszewski/5_by_n/n{n}/graph{i}.lgbin")

	for k in range(2, 26):
		if f"k{k}" not in os.listdir("../random graphs/matuszewski/k_by_10"):
			os.mkdir(f"../random graphs/matuszewski/k_by_10/k{k}")
			for i in range(100):
				ng = random_layered_graph_connect_help_edgecount(k, 10, 20 * (k-1))
				print(f"k={k} graph {i}")
				ng.write_out(f"../random graphs/matuszewski/k_by_10/k{k}/graph{i}.lgbin")

	for d in range(15, 96, 5):
		if f"d{d}" not in os.listdir("../random graphs/matuszewski/10_by_10_density"):
			os.mkdir(f"../random graphs/matuszewski/10_by_10_density/d{d}")
			for i in range(100):
				ng = random_layered_graph_connect_help(10, 10, d / 100)
				print(f"d={d} graph {i}")
				ng.write_out(f"../random graphs/matuszewski/10_by_10_density/d{d}/graph{i}.lgbin")


def generate_big_n_by_n_graphs(seed=220):
	random.seed(seed)
	if "n_by_n" not in os.listdir("../random graphs"):
		os.mkdir("../random graphs/n_by_n")

	for n in range(20, 51, 5):
		if f"n{n}" not in os.listdir("../random graphs/n_by_n"):
			os.mkdir(f"../random graphs/n_by_n/n{n}")
			for i in range(10):
				ng = random_layered_graph_connect_help_edgecount(n, n, 2 * n * (n - 1))
				print(f"n={n} graph {i+1}")
				ng.write_out(f"../random graphs/n_by_n/n{n}/graph{i}.lgbin")


def generate_rectangle_graphs(seed=2200):
	random.seed(seed)
	if "rectangles" not in os.listdir("../random graphs"):
		os.mkdir("../random graphs/rectangles")

	for n in range(10, 31):
		k = 100 - (3*n)
		if f"k{k}n{n}" not in os.listdir("../random graphs/rectangles"):
			os.mkdir(f"../random graphs/rectangles/k{k}n{n}")
			for i in range(10):
				ng = random_layered_graph_edgecount_drop_unconnected(k, n, 2 * n * (k - 1))
				print(f"k={k}/n={n} graph {i+1}")
				ng.write_out(f"../random graphs/rectangles/k{k}n{n}/graph{i}.lgbin")


def generate_ratio_graphs(seed=22200):
	random.seed(seed)
	if "ratio" not in os.listdir("../random graphs"):
		os.mkdir(f"../random graphs/ratio")

	for kn in [10, 15, 20, 25, 30, 35]:
		if f"r1k{kn}n{kn}" not in os.listdir("../random graphs/ratio"):
			os.mkdir(f"../random graphs/ratio/r1k{kn}n{kn}")
			for i in range(10):
				ng = random_layered_graph_edgecount_drop_unconnected(kn, kn, 2 * kn * (kn - 1))
				print(f"k={kn}/n={kn} graph {i + 1}")
				ng.write_out(f"../random graphs/ratio/r1k{kn}n{kn}/graph{i}.lgbin")
	for n in [6, 10, 13, 17, 20, 24]:
		k = n * 2
		if f"r2k{k}n{n}" not in os.listdir("../random graphs/ratio"):
			os.mkdir(f"../random graphs/ratio/r2k{k}n{n}")
			for i in range(10):
				ng = random_layered_graph_edgecount_drop_unconnected(k, n, 2 * n * (k - 1))
				print(f"k={k}/n={n} graph {i + 1}")
				ng.write_out(f"../random graphs/ratio/r2k{k}n{n}/graph{i}.lgbin")
	for n in [5, 8, 11, 14, 17, 20]:
		k = n * 3
		if f"r3k{k}n{n}" not in os.listdir("../random graphs/ratio"):
			os.mkdir(f"../random graphs/ratio/r3k{k}n{n}")
			for i in range(10):
				ng = random_layered_graph_edgecount_drop_unconnected(k, n, 2 * n * (k - 1))
				print(f"k={k}/n={n} graph {i + 1}")
				ng.write_out(f"../random graphs/ratio/r3k{k}n{n}/graph{i}.lgbin")


def generate_ratio1dot5_graphs_d3_with_big_layers(seed=222201):
	random.seed(seed)
	if "big_layer" not in os.listdir("../random graphs"):
		os.mkdir("../random graphs/big_layer")

	for n in [8, 12, 16, 20, 24]:
		k = int(n * 1.5)
		if f"k{k}n{n}" not in os.listdir("../random graphs/big_layer"):
			os.mkdir(f"../random graphs/big_layer/k{k}n{n}")
			for i in range(50):
				n_big = round(k / 10) - 1
				lcounts = [n] * k
				blayer = random.choice(range(n_big // 2, k - (n_big // 2)))
				lcounts[blayer] *= 3
				j = 1
				while n_big > 0:
					if blayer + j < len(lcounts):
						lcounts[blayer + j] *= 3
						n_big -= 1
					if n_big > 0 and blayer - j >= 0:
						lcounts[blayer - j] *= 3
						n_big -= 1
					j += 1
				print(lcounts)
				ng = random_layered_graph_edgecount_difflayers(round(1.5 * n * (k - 1)), lcounts)
				# print([len(lay) for lay in ng.layers.values()])
				print(f"k={k}/n={n} graph {i + 1}")
				ng.write_out(f"../random graphs/big_layer/k{k}n{n}/graph{i}.lgbin")


def generate_ratio1dot5_graphs_d3_triangle(seed=222222):
	# random.seed(seed)
	if "triangle" not in os.listdir("../random graphs"):
		os.mkdir("../random graphs/triangle")

	for n in [12, 16, 20, 24]:
		k = int(n * 1.5)
		if f"k{k}n{n}" not in os.listdir("../random graphs/triangle"):
			os.mkdir(f"../random graphs/triangle/k{k}n{n}")
			lcounts = []
			for xv in range(k):
				lcounts.append(round(xv * ((2 * n - 1) / (k - 1)) + 1))
			print(lcounts)
			for i in range(50):
				ng = random_layered_graph_edgecount_difflayers(round(1.5 * n * (k - 1)), lcounts)
				# print([len(lay) for lay in ng.layers.values()])
				print(f"k={k}/n={n} graph {i + 1}")
				ng.write_out(f"../random graphs/triangle/k{k}n{n}/graph{i}.lgbin")


def generate_ratio_graphs_degree_3(seed=22201):
	# random.seed(seed)
	if "ratio_d3" not in os.listdir("../random graphs"):
		os.mkdir(f"../random graphs/ratio_d3")

	# for kn in [10, 15, 20, 25, 30, 35]:
	# 	# if f"r1k{kn}n{kn}" not in os.listdir("../random graphs/ratio_d3"):
	# 	# 	os.mkdir(f"../random graphs/ratio_d3/r1k{kn}n{kn}")
	# 	for i in range(10, 20):
	# 		ng = random_layered_graph_edgecount_drop_unconnected(kn, kn, round(1.5 * kn * (kn - 1)))
	# 		print(f"k={kn}/n={kn} graph {i + 1}")
	# 		ng.write_out(f"../random graphs/ratio_d3/r1k{kn}n{kn}/graph{i}.lgbin")
	for n in [8, 12, 16, 20, 24, 28]:
		k = int(n * 1.5)
		# if f"r1.5k{k}n{n}" not in os.listdir("../random graphs/ratio_d3"):
		# 	os.mkdir(f"../random graphs/ratio_d3/r1.5k{k}n{n}")
		for i in range(20, 50):
			ng = random_layered_graph_edgecount_drop_unconnected(k, n, round(1.5 * n * (k - 1)))
			print(f"k={k}/n={n} graph {i + 1}")
			ng.write_out(f"../random graphs/ratio_d3/r1.5k{k}n{n}/graph{i}.lgbin")
	# for n in [6, 10, 13, 17, 20, 24]:
	# 	k = n * 2
	# 	# if f"r2k{k}n{n}" not in os.listdir("../random graphs/ratio_d3"):
	# 	# 	os.mkdir(f"../random graphs/ratio_d3/r2k{k}n{n}")
	# 	for i in range(10, 20):
	# 		ng = random_layered_graph_edgecount_drop_unconnected(k, n, round(1.5 * n * (k - 1)))
	# 		print(f"k={k}/n={n} graph {i + 1}")
	# 		ng.write_out(f"../random graphs/ratio_d3/r2k{k}n{n}/graph{i}.lgbin")
	# for n in [5, 8, 11, 14, 17, 20]:
	# 	k = n * 3
	# 	if f"r3k{k}n{n}" not in os.listdir("../random graphs/ratio_d3"):
	# 		os.mkdir(f"../random graphs/ratio_d3/r3k{k}n{n}")
	# 		for i in range(10):
	# 			ng = random_layered_graph_edgecount_drop_unconnected(k, n, round(1.5 * n * (k - 1)))
	# 			print(f"k={k}/n={n} graph {i + 1}")
	# 			ng.write_out(f"../random graphs/ratio_d3/r3k{k}n{n}/graph{i}.lgbin")


if __name__ == '__main__':
	# generate_gange_dataset(seed=22)
	# generate_random_density_set(seed=49)
	# generate_random_fixed_density_set(seed=71)
	# generate_extended_matuszewski_datsets()
	# generate_big_n_by_n_graphs()
	# generate_rectangle_graphs()
	# generate_ratio_graphs()
	# generate_ratio_graphs_degree_3()
	# generate_ratio1dot5_graphs_d3_with_big_layers()
	generate_ratio1dot5_graphs_d3_triangle()

	# gr = src.read_data.read("../random graphs/big_layer/k36n24/graph1.lgbin")
	# gr = src.read_data.read("../random graphs/ratio_d3/r1.5k24n16/graph44.lgbin")
	# print(gr.layers)
	# print(gr.nodes)
	# print(gr.node_ids)
	# print(gr.edge_ids)
	# adj = gr.get_adj_list()
	# print(sum((len(adj[v.id]) for v in gr))/gr.n_nodes)
	# print([len(lay) for lay in gr.layers.values()])
	# print([len(gr.get_edges_by_layer()[ls]) for ls in range(gr.n_layers - 1)])

	# gr = src.read_data.read("../random graphs/matuszewski/10_by_10_density/d15/graph0.lgbin")
	# gr = src.read_data.read("../random graphs/rectangles/k70n20/graph0.lgbin")
	# gr = src.read_data.read("../random graphs/triangle/k30n20/graph9.lgbin")
	# gr = src.read_data.read("../random graphs/ratio_d3/r1.5k18n12/graph8.lgbin")
	# print(gr.n_layers)
	# gr = random_layered_graph_connect_help_edgecount(3, 10, 35)
	# src.vis.draw_graph(gr, "rand", gravity=True, nested=True)
