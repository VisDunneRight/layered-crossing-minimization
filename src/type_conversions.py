from src import graph, layering
import networkx as nx
import igraph as ig


def layered_graph_to_nx_graph(g: graph.LayeredGraph):
	edge_list = [(e.n1.name, e.n2.name) for e in g.edges]
	return nx.Graph(incoming_graph_data=edge_list)


def layered_graph_to_igraph(g: graph.LayeredGraph):
	edges = [(e.n1.name - 1, e.n2.name - 1) for e in g.edges]
	ret_g = ig.Graph(n=len(g.nodes))
	ret_g.add_edges(edges)
	ret_g.vs["layer"] = [v.layer for v in g.nodes]
	return ret_g


def dagmar_nx_to_layered_graph(nxg: nx.Graph):
	g = graph.LayeredGraph()
	lv = nx.get_node_attributes(nxg, "hierarchy.level")
	for v in nxg.nodes:
		g.add_node(int(lv[v])+1, name=int(v[1:])+1)
	for edge in nxg.edges:
		g.add_edge(int(edge[0][1:])+1, int(edge[1][1:])+1)
	g.add_anchors()
	g.relayer()
	g.y_val_setup()
	return g


def north_nx_to_layered_graph(nxg: nx.Graph):
	return layering.create_layered_graph_from_directed_nx_graph(nxg, 4, 2)
