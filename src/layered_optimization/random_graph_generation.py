import time
import networkx as nx
from layered_optimization.read_data import read
from layered_optimization import vis, graph
import random
import os
import networkx
import math
import shutil


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


def random_layered_graph_connect_three_halves(n, m, pct_long_edges=0):
	"""
	:param n: number of nodes
	:param m: number of edges
	:param pct_long_edges: % of edges which will be between any 
	:return: LayeredGraph object g

	Generates a uniformly random, connected layered graph with n nodes and m edges
	where the ratio #layers:#nodes per layer is roughly 3:2
	"""

	assert m >= n, "graph will not be connected"
	assert n > 5, "too small"

	k = round(math.sqrt(3 * n / 2))
	n_min = math.floor(math.sqrt(2 * n / 3))
	n_extra = n - (n_min * k)
	n_fails = 0

	while True:
		g = graph.LayeredGraph()

		n_added = 0
		for i in range(k):  # add nodes
			for j in range(n_min):
				if n_added < n:
					g.add_node(i)
					n_added += 1
		if n_extra > 0:
			for _ in range(n_extra):
				lid = random.randint(0, k)  # can add to an extra layer at the end
				g.add_node(lid)

		g_n_by_l = g.get_ids_by_layer()
		all_edges = [(nd1, nd2) for lid in range(0, g.n_layers - 1) for nd1 in g_n_by_l[lid] for nd2 in g_n_by_l[lid + 1]]
		num_short_edges = round(m * (1 - pct_long_edges))
		edges_to_add = random.sample(all_edges, num_short_edges)
		for edge in edges_to_add:
			g.add_edge(edge[0], edge[1])

		n_edges = num_short_edges
		while n_edges < m:
			nd1, nd2 = random.randint(0, n - 1), random.randint(0, n - 1)
			if k > 3:
				while abs(g[nd1].layer - g[nd2].layer) <= 1:
					nd2 = random.randint(0, n - 1)
			else:
				while g[nd1].layer == g[nd2].layer:
					nd2 = random.randint(0, n - 1)
			if (nd1, nd2) not in g.edge_ids and (nd2, nd1) not in g.edge_ids:
				g.add_edge(nd1, nd2)
				n_edges += 1
		
		if g.is_connected():
			g.add_anchors()
			g.relayer()
			g.y_val_setup()
			if n_fails > 0:
				print("> success")
			return g
		else:
			n_fails += 1
			print(f"failures: {n_fails}", end='\r')


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
			g.y_val_setup()
			return g
		else:
			print("fail")


def random_layered_graph_edgecount_difflayers(n_edges, layercounts: list):
	"""
	:param n_edges: number of edges in the resultant graph
	:param layercounts: list of #nodes for each layer
	:return: LayeredGraph object g

	Randomly samples edges over the entire network. Any nodes left unconnected are removed
	TODO sample from all adjacent nodes
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

	if "random_graphs" not in os.listdir("datasets"):
		os.mkdir("datasets/random_graphs")
	if "gange" not in os.listdir("datasets/random_graphs"):
		os.mkdir("datasets/random_graphs/gange")

	for k in range(3, 11):
		for n in range(7, 11 if k < 5 else (10 if k < 8 else 9)):
			if f"g{k}_{n}" not in os.listdir("datasets/random_graphs/gange"):
				os.mkdir(f"datasets/random_graphs/gange/g{k}_{n}")
			for i in range(10):
				ng = true_random_connected_layered_graph(k, n, 0.2)
				ng.write_out(f"datasets/random_graphs/gange/g{k}_{n}/graph{i}.lgbin")


def generate_random_density_set(seed=None):
	if seed is not None:
		random.seed(seed)

	if "random_graphs" not in os.listdir("datasets"):
		os.mkdir("datasets/random_graphs")
	if "density_exp" not in os.listdir("datasets/random_graphs"):
		os.mkdir("datasets/random_graphs/density_exp")

	for d in range(14, 51, 2):
		if f"d{d}" not in os.listdir("datasets/random_graphs/density_exp"):
			os.mkdir(f"datasets/random_graphs/density_exp/d{d}")
		for i in range(10):
			ng = true_random_connected_layered_graph(5, 10, d / 100)
			print(f"d={d} graph {i}")
			ng.write_out(f"datasets/random_graphs/density_exp/d{d}/graph{i}.lgbin")


def generate_random_fixed_density_set(seed=None):
	if seed is not None:
		random.seed(seed)

	if "random_graphs" not in os.listdir("datasets"):
		os.mkdir("datasets/random_graphs")
	if "fixed_density_exp" not in os.listdir("datasets/random_graphs"):
		os.mkdir("datasets/random_graphs/fixed_density_exp")

	for k in range(3, 21):
		if f"k{k}" not in os.listdir("datasets/random_graphs/fixed_density_exp"):
			os.mkdir(f"datasets/random_graphs/fixed_density_exp/k{k}")
		for i in range(10):
			ng = random_layered_graph_connect_help(k, 10, 0.15)
			print(f"k={k} graph {i}")
			ng.write_out(f"datasets/random_graphs/fixed_density_exp/k{k}/graph{i}.lgbin")


def generate_extended_matuszewski_datsets(seed=None):
	if seed is not None:
		random.seed(seed)
	if "random_graphs" not in os.listdir("datasets"):
		os.mkdir("datasets/random_graphs")
	if "matuszewski" not in os.listdir("datasets/random_graphs"):
		os.mkdir("datasets/random_graphs/matuszewski")
		os.mkdir("datasets/random_graphs/matuszewski/5_by_n")
		os.mkdir("datasets/random_graphs/matuszewski/k_by_10")
		os.mkdir("datasets/random_graphs/matuszewski/10_by_10_density")

	for n in range(10, 101, 10):
		if f"n{n}" not in os.listdir("datasets/random_graphs/matuszewski/5_by_n"):
			os.mkdir(f"datasets/random_graphs/matuszewski/5_by_n/n{n}")
			for i in range(100):
				ng = random_layered_graph_connect_help_edgecount(5, n, 8 * n)
				print(f"n={n} graph {i}")
				ng.write_out(f"datasets/random_graphs/matuszewski/5_by_n/n{n}/graph{i}.lgbin")

	for k in range(2, 26):
		if f"k{k}" not in os.listdir("datasets/random_graphs/matuszewski/k_by_10"):
			os.mkdir(f"datasets/random_graphs/matuszewski/k_by_10/k{k}")
			for i in range(100):
				ng = random_layered_graph_connect_help_edgecount(k, 10, 20 * (k-1))
				print(f"k={k} graph {i}")
				ng.write_out(f"datasets/random_graphs/matuszewski/k_by_10/k{k}/graph{i}.lgbin")

	for d in range(15, 96, 5):
		if f"d{d}" not in os.listdir("datasets/random_graphs/matuszewski/10_by_10_density"):
			os.mkdir(f"datasets/random_graphs/matuszewski/10_by_10_density/d{d}")
			for i in range(100):
				ng = random_layered_graph_connect_help(10, 10, d / 100)
				print(f"d={d} graph {i}")
				ng.write_out(f"datasets/random_graphs/matuszewski/10_by_10_density/d{d}/graph{i}.lgbin")


def generate_big_n_by_n_graphs(seed=220):
	random.seed(seed)
	if "n_by_n" not in os.listdir("datasets/random_graphs"):
		os.mkdir("datasets/random_graphs/n_by_n")

	for n in range(20, 51, 5):
		if f"n{n}" not in os.listdir("datasets/random_graphs/n_by_n"):
			os.mkdir(f"datasets/random_graphs/n_by_n/n{n}")
			for i in range(10):
				ng = random_layered_graph_connect_help_edgecount(n, n, 2 * n * (n - 1))
				print(f"n={n} graph {i+1}")
				ng.write_out(f"datasets/random_graphs/n_by_n/n{n}/graph{i}.lgbin")


def generate_rectangle_graphs(seed=2200):
	random.seed(seed)
	if "rectangles" not in os.listdir("datasets/random_graphs"):
		os.mkdir("datasets/random_graphs/rectangles")

	for n in range(10, 31):
		k = 100 - (3*n)
		if f"k{k}n{n}" not in os.listdir("datasets/random_graphs/rectangles"):
			os.mkdir(f"datasets/random_graphs/rectangles/k{k}n{n}")
			for i in range(10):
				ng = random_layered_graph_edgecount_drop_unconnected(k, n, 2 * n * (k - 1))
				print(f"k={k}/n={n} graph {i+1}")
				ng.write_out(f"datasets/random_graphs/rectangles/k{k}n{n}/graph{i}.lgbin")


def generate_ratio_graphs(seed=22200):
	random.seed(seed)
	if "ratio" not in os.listdir("datasets/random_graphs"):
		os.mkdir(f"datasets/random_graphs/ratio")

	for kn in [10, 15, 20, 25, 30, 35]:
		if f"r1k{kn}n{kn}" not in os.listdir("datasets/random_graphs/ratio"):
			os.mkdir(f"datasets/random_graphs/ratio/r1k{kn}n{kn}")
			for i in range(10):
				ng = random_layered_graph_edgecount_drop_unconnected(kn, kn, 2 * kn * (kn - 1))
				print(f"k={kn}/n={kn} graph {i + 1}")
				ng.write_out(f"datasets/random_graphs/ratio/r1k{kn}n{kn}/graph{i}.lgbin")
	for n in [6, 10, 13, 17, 20, 24]:
		k = n * 2
		if f"r2k{k}n{n}" not in os.listdir("datasets/random_graphs/ratio"):
			os.mkdir(f"datasets/random_graphs/ratio/r2k{k}n{n}")
			for i in range(10):
				ng = random_layered_graph_edgecount_drop_unconnected(k, n, 2 * n * (k - 1))
				print(f"k={k}/n={n} graph {i + 1}")
				ng.write_out(f"datasets/random_graphs/ratio/r2k{k}n{n}/graph{i}.lgbin")
	for n in [5, 8, 11, 14, 17, 20]:
		k = n * 3
		if f"r3k{k}n{n}" not in os.listdir("datasets/random_graphs/ratio"):
			os.mkdir(f"datasets/random_graphs/ratio/r3k{k}n{n}")
			for i in range(10):
				ng = random_layered_graph_edgecount_drop_unconnected(k, n, 2 * n * (k - 1))
				print(f"k={k}/n={n} graph {i + 1}")
				ng.write_out(f"datasets/random_graphs/ratio/r3k{k}n{n}/graph{i}.lgbin")


def generate_ratio1dot5_graphs_d3_with_big_layers(seed=222201):
	random.seed(seed)
	if "big_layer" not in os.listdir("datasets/random_graphs"):
		os.mkdir("datasets/random_graphs/big_layer")

	for n in [8, 12, 16, 20, 24]:
		k = int(n * 1.5)
		if f"k{k}n{n}" not in os.listdir("datasets/random_graphs/big_layer"):
			os.mkdir(f"datasets/random_graphs/big_layer/k{k}n{n}")
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
				ng.write_out(f"datasets/random_graphs/big_layer/r1.5k{k}n{n}/graph{i}.lgbin")


def generate_ratio1dot5_graphs_d3_triangle(seed=222222):
	# random.seed(seed)
	if "triangle" not in os.listdir("datasets/random_graphs"):
		os.mkdir("datasets/random_graphs/triangle")

	for n in [12, 16, 20, 24, 28]:
		k = int(n * 1.5)
		if f"k{k}n{n}" not in os.listdir("datasets/random_graphs/triangle"):
			os.mkdir(f"datasets/random_graphs/triangle/k{k}n{n}")
			lcounts = []
			for xv in range(k):
				lcounts.append(round(xv * ((2 * n - 1) / (k - 1)) + 1))
			print(lcounts)
			for i in range(50):
				ng = random_layered_graph_edgecount_difflayers(round(1.5 * n * (k - 1)), lcounts)
				# print([len(lay) for lay in ng.layers.values()])
				print(f"k={k}/n={n} graph {i + 1}")
				ng.write_out(f"datasets/random_graphs/triangle/r1.5k{k}n{n}/graph{i}.lgbin")


def generate_ratio_graphs_degree_3(seed=22201):
	# random.seed(seed)
	if "ratio_d3" not in os.listdir("datasets/random_graphs"):
		os.mkdir(f"datasets/random_graphs/ratio_d3")

	# for kn in [10, 15, 20, 25, 30, 35]:
	# 	# if f"r1k{kn}n{kn}" not in os.listdir("datasets/random_graphs/ratio_d3"):
	# 	# 	os.mkdir(f"datasets/random_graphs/ratio_d3/r1k{kn}n{kn}")
	# 	for i in range(10, 20):
	# 		ng = random_layered_graph_edgecount_drop_unconnected(kn, kn, round(1.5 * kn * (kn - 1)))
	# 		print(f"k={kn}/n={kn} graph {i + 1}")
	# 		ng.write_out(f"datasets/random_graphs/ratio_d3/r1k{kn}n{kn}/graph{i}.lgbin")
	for n in [8, 12, 16, 20, 24, 28]:
		k = int(n * 1.5)
		# if f"r1.5k{k}n{n}" not in os.listdir("datasets/random_graphs/ratio_d3"):
		# 	os.mkdir(f"datasets/random_graphs/ratio_d3/r1.5k{k}n{n}")
		for i in range(20, 50):
			ng = random_layered_graph_edgecount_drop_unconnected(k, n, round(1.5 * n * (k - 1)))
			print(f"k={k}/n={n} graph {i + 1}")
			ng.write_out(f"datasets/random_graphs/ratio_d3/r1.5k{k}n{n}/graph{i}.lgbin")
	# for n in [6, 10, 13, 17, 20, 24]:
	# 	k = n * 2
	# 	# if f"r2k{k}n{n}" not in os.listdir("datasets/random_graphs/ratio_d3"):
	# 	# 	os.mkdir(f"datasets/random_graphs/ratio_d3/r2k{k}n{n}")
	# 	for i in range(10, 20):
	# 		ng = random_layered_graph_edgecount_drop_unconnected(k, n, round(1.5 * n * (k - 1)))
	# 		print(f"k={k}/n={n} graph {i + 1}")
	# 		ng.write_out(f"datasets/random_graphs/ratio_d3/r2k{k}n{n}/graph{i}.lgbin")
	# for n in [5, 8, 11, 14, 17, 20]:
	# 	k = n * 3
	# 	if f"r3k{k}n{n}" not in os.listdir("datasets/random_graphs/ratio_d3"):
	# 		os.mkdir(f"datasets/random_graphs/ratio_d3/r3k{k}n{n}")
	# 		for i in range(10):
	# 			ng = random_layered_graph_edgecount_drop_unconnected(k, n, round(1.5 * n * (k - 1)))
	# 			print(f"k={k}/n={n} graph {i + 1}")
	# 			ng.write_out(f"datasets/random_graphs/ratio_d3/r3k{k}n{n}/graph{i}.lgbin")


def generate_small_ratio_graphs(seed=22201):
	if "ratio_small_d3" not in os.listdir("datasets/random_graphs"):
		os.mkdir(f"datasets/random_graphs/ratio_small_d3")

	n = 4
	k = int(n * 1.5)
	if f"r1.5k{k}n{n}" not in os.listdir("datasets/random_graphs/ratio_small_d3"):
		os.mkdir(f"datasets/random_graphs/ratio_small_d3/r1.5k{k}n{n}")
	for i in range(10):
		ng = random_layered_graph_edgecount_drop_unconnected(k, n, round(1.5 * n * (k - 1)))
		print(f"k={k}/n={n} graph {i + 1}")
		ng.write_out(f"datasets/random_graphs/ratio_small_d3/r1.5k{k}n{n}/graph{i}.lgbin")


def generate_nx_random_graph_dataset():
	for n in range(5, 101):
		gid = 0
		ngraphs = 100
		foldername = "datasets/random_graphs/networkx2"
		print(f"N = {n}")
		for i in range(ngraphs):
			gr = networkx.gnm_random_graph(n, round(n * (1 + i / ngraphs)))  # uniform, |E| in [|V|, 2|V|)
			with open(f"{foldername}/graph_{n}_{gid}", 'w') as fd:
				for ed in gr.edges:
					fd.write(f"{ed[0]},{ed[1]}\n")
			gid += 1


def generate_three_to_two_dataset():
	if "uniform_layered" in os.listdir("datasets/random_graphs"):
		overwrite = input("Overwrite existing files? [y]/n: ")
		if overwrite.lower().strip() == 'n':
			return
		else:
			shutil.rmtree("datasets/random_graphs/uniform_layered")
	
	os.mkdir("datasets/random_graphs/uniform_layered")

	for n in range(5, 121, 5):
		os.mkdir(f"datasets/random_graphs/uniform_layered/n{n}")
		print(f"\nBucket {n}")
		print('='*75)
		for i in range(300):
			nnds = random.randint(n, n+4)
			gr = random_layered_graph_connect_three_halves(nnds if nnds > 5 else 6, round(nnds * 1.4), pct_long_edges=0.2)
			gr.write_out(f"datasets/random_graphs/uniform_layered/n{n}/graph{i}.lgbin")


def make_up_extra_three_two(bins_to_cur_num):
	for binval, cur_num in bins_to_cur_num.items():
		cur_count = cur_num
		desired_num = 100
		while cur_count < desired_num:
			nnds = random.randint(round(0.35 * binval), round(0.5 * binval))
			nnds_bin = int(nnds / 5) * 5
			cur_num_in_nnds_bin = len(os.listdir(f"datasets/random_graphs/uniform_layered/n{nnds_bin}"))
			gr = random_layered_graph_connect_three_halves(nnds if nnds > 5 else 6, round(nnds * 1.4), pct_long_edges=0.2)
			graph_tnodes_bin = int(gr.n_nodes / 5) * 5
			if graph_tnodes_bin == binval:
				cur_count += 1
				print(f"Success! {cur_count} / {desired_num} for bin {binval}")
				gr.write_out(f"datasets/random_graphs/uniform_layered/n{nnds_bin}/graph{cur_num_in_nnds_bin}.lgbin")
			else:
				print(f"Fail. Target bin was {binval}, graph tnodes={gr.n_nodes}, actual bin={graph_tnodes_bin}")


def generate_scale_free_networkx_dataset():
	if "scale_free" not in os.listdir("datasets/random_graphs"):
		os.mkdir("datasets/random_graphs/scale_free")

	for n in range(5, 151, 5):
		os.mkdir(f"datasets/random_graphs/scale_free/n{n}")
		for i in range(70):
			nnds = random.randint(n, n+4)
			gr = networkx.scale_free_graph(nnds, alpha=0.35, beta=0.6, gamma=0.05)
			networkx.write_adjlist(gr, f"datasets/random_graphs/scale_free/n{n}/graph{i}.adjlist")
