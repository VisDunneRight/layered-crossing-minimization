from src import graph
import networkx as nx
import igraph as ig


def layered_graph_to_nx_graph(g: graph.LayeredGraph):
	edge_list = [(e.n1.name, e.n2.name) for e in g.edges]
	return nx.Graph(incoming_graph_data=edge_list)


def layered_graph_to_igraph(g: graph.LayeredGraph):
	edges = [(e.n1.name - 1, e.n2.name - 1) for e in g.edges]
	ret_g = ig.Graph(len(g.nodes), edges)
	ret_g.vs["layer"] = [v.layer for v in g.nodes]
	return ret_g
