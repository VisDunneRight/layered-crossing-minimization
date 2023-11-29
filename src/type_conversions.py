from src import graph, layering
import networkx as nx
import igraph as ig


def layered_graph_to_nx_graph(g: graph.LayeredGraph):
	edge_list = [(e.n1.id, e.n2.id) for e in g.edges]
	return nx.Graph(incoming_graph_data=edge_list)


def layered_graph_to_igraph(g: graph.LayeredGraph):
	old_id_to_idx = {v.id: i for i, v in enumerate(g.nodes)}
	edges = [(old_id_to_idx[e.n1.id], old_id_to_idx[e.n2.id]) for e in g.edges]
	ret_g = ig.Graph(n=len(g.nodes))
	ret_g.add_edges(edges)
	ret_g.vs["layer"] = [v.layer for v in g.nodes]
	ret_g.vs["id"] = [v.id for v in g.nodes]
	return ret_g


def dagmar_nx_to_layered_graph(nxg: nx.Graph):
	g = graph.LayeredGraph()
	lv = nx.get_node_attributes(nxg, "hierarchy.level")
	for v in nxg.nodes:
		g.add_node(int(lv[v])+1, idx=int(v[1:]))
	for edge in nxg.edges:
		g.add_edge(int(edge[0][1:]), int(edge[1][1:]))
	g.add_anchors()
	g.relayer()
	g.y_val_setup()
	return g


def north_nx_to_layered_graph(nxg: nx.Graph, w, c):
	return layering.create_layered_graph_from_directed_nx_graph(nxg, w, c)


def nx_with_separate_layerings_to_layered_graph(nxg: nx.Graph, layer_assign):
	g = graph.LayeredGraph()
	for v in nxg.nodes:
		g.add_node(layer_assign[v], idx=int(v))
	for edge in nxg.edges:
		g.add_edge(int(edge[0]), int(edge[1]))
	g.add_anchors()
	g.y_val_setup()
	return g
