import os
import graph
import read_data
from random import choice, randint
import heuristics
from src.optimization import LayeredOptimizer
from sklearn.cluster import SpectralClustering

import src.vis


def test_all_unique_ids(test_graphs):
	print("testing if all node names are unique:")
	for grf in test_graphs:
		print(f"\t{grf}... ", end='')
		gr = read_data.read(grf)
		seen_ids = set()
		for nd in gr.nodes:
			assert nd.id not in seen_ids, f"failed on {grf} with double node {nd.id}"
			seen_ids.add(nd.id)
		print("pass")


def test_ascending_names_from_zero(test_graphs):
	print("testing if node names are 0, 1, 2,...:")
	for grf in test_graphs:
		print(f"\t{grf}... ", end='')
		gr = read_data.read(grf)
		nds = sorted([nd.id for nd in gr.nodes])
		comp = list(range(len(gr.nodes)))
		assert nds == comp, f"failed on {grf}"
		print("pass")


def test_method_consistency_for_counts(test_graphs):
	print("testing if num nodes and num layers is consistent after method calls:")
	for grf in test_graphs:
		print(f"\t{grf}... ", end='')
		gr = read_data.read(grf)
		start = gr.n_nodes
		lstart = gr.n_layers
		matrix = gr.adjacency_matrix()
		assert gr.n_nodes == start and gr.n_layers == lstart and len(gr.nodes) == start and len(gr.layers) == lstart, f"failed on {grf}"
		adj = gr.create_normal_adj_list()
		assert gr.n_nodes == start and gr.n_layers == lstart and len(gr.nodes) == start and len(gr.layers) == lstart, f"failed on {grf}"
		asdf1 = gr.get_edge_names_by_layer()
		assert gr.n_nodes == start and gr.n_layers == lstart and len(gr.nodes) == start and len(gr.layers) == lstart, f"failed on {grf}"
		asdf2 = gr.get_names_by_layer()
		assert gr.n_nodes == start and gr.n_layers == lstart and len(gr.nodes) == start and len(gr.layers) == lstart, f"failed on {grf}"
		asdf3 = gr.get_edges_by_layer()
		assert gr.n_nodes == start and gr.n_layers == lstart and len(gr.nodes) == start and len(gr.layers) == lstart, f"failed on {grf}"
		# heuristics.barycenter(gr, 5)
		# assert gr.n_nodes == start and gr.n_layers == lstart and len(gr.nodes) == start and len(gr.layers) == lstart, f"failed on {grf}"
		gr.add_anchors()
		assert gr.n_nodes == start and gr.n_layers == lstart and len(gr.nodes) == start and len(gr.layers) == lstart, f"failed on {grf}"
		gr.relayer()
		assert gr.n_nodes == start and gr.n_layers == lstart and len(gr.nodes) == start and len(gr.layers) == lstart, f"failed on {grf}"
		rassgn = list(SpectralClustering(n_clusters=3, assign_labels="discretize", affinity="precomputed").fit(gr.adjacency_matrix()).labels_)
		asdf4 = gr.stacked_graph_from_subgraph_nodes(rassgn)
		assert gr.n_nodes == start and gr.n_layers == lstart and len(gr.nodes) == start and len(gr.layers) == lstart, f"failed on {grf}"
		gr.add_node(999999)
		start += 1
		lstart += 1
		assert gr.n_nodes == start and gr.n_layers == lstart and len(gr.nodes) == start and len(gr.layers) == lstart, f"failed on {grf}"
		print("pass")


def extract_graph_raw_data(layered_g):
	nds = []
	eds = []
	for nd in layered_g:
		nds.append((nd.id, nd.layer, nd.is_anchor_node, nd.stacked))
	for ed in layered_g.edges:
		eds.append((ed.n1.id, ed.n2.id))
	return nds, eds


def test_method_consistency_for_graph_structure(test_graphs):
	print("testing if the node/edge structure of the graph (including ID consistency, order of node/edge lists) is unchanged by method calls:")
	for grf in test_graphs:
		print(f"\t{grf}... ", end='')
		gr = read_data.read(grf)
		n_orig, e_orig = extract_graph_raw_data(gr)
		matrix = gr.adjacency_matrix()
		n_test, e_test = extract_graph_raw_data(gr)
		assert n_orig == n_test and e_orig == e_test, f"failed on {grf}"
		adj = gr.create_normal_adj_list()
		n_test, e_test = extract_graph_raw_data(gr)
		assert n_orig == n_test and e_orig == e_test, f"failed on {grf}"
		asdf1 = gr.get_edge_names_by_layer()
		n_test, e_test = extract_graph_raw_data(gr)
		assert n_orig == n_test and e_orig == e_test, f"failed on {grf}"
		asdf2 = gr.get_names_by_layer()
		n_test, e_test = extract_graph_raw_data(gr)
		assert n_orig == n_test and e_orig == e_test, f"failed on {grf}"
		asdf3 = gr.get_edges_by_layer()
		n_test, e_test = extract_graph_raw_data(gr)
		assert n_orig == n_test and e_orig == e_test, f"failed on {grf}"
		gr.add_anchors()
		n_test, e_test = extract_graph_raw_data(gr)
		assert n_orig == n_test and e_orig == e_test, f"failed on {grf}"
		gr.relayer()
		n_test, e_test = extract_graph_raw_data(gr)
		assert n_orig == n_test and e_orig == e_test, f"failed on {grf}"
		gr.add_node(999999)
		n_test, e_test = extract_graph_raw_data(gr)
		n_orig.append((gr.n_nodes - 1, 999999, False, False))
		assert n_orig == n_test and e_orig == e_test, f"failed on {grf}"
		print("pass")


def test_optimization_on_graphs_with_non_ascending_names(test_graphs):
	print("testing if the result of optimizing a graph is unchanged if the node names are not 0, 1, 2,...:")
	for grf in test_graphs:
		print(f"\t{grf}... ", end='')
		gr = read_data.read(grf)
		optim = LayeredOptimizer(gr)
		optim.symmetry_breaking = True
		optim.optimize_layout()

		nd1 = gr.nodes.pop(2)
		nd2 = gr.nodes.pop(5)
		nd3 = gr.nodes.pop(7)
		nd1.id = 230
		nd2.id = 304
		nd3.id = 455
		gr.nodes.append()


if __name__ == '__main__':
	n_cases = 5
	test_folds = [f"../Rome-Lib/graficon" + str(choice(list(range(10, 65)))) + "nodi/" for i in range(n_cases)]
	test_set = [l1 + choice(list(os.listdir(l1))) for l1 in test_folds]
	test_all_unique_ids(test_set)
	test_ascending_names_from_zero(test_set)
	test_method_consistency_for_counts(test_set)
	test_method_consistency_for_graph_structure(test_set)
	print("all cases pass")
