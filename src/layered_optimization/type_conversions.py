from layered_optimization import graph, layering
import networkx as nx


def layered_graph_to_nx_graph(g: graph.LayeredGraph):
	edge_list = [(e.n1.id, e.n2.id) for e in g.edges]
	return nx.Graph(incoming_graph_data=edge_list)


def dagmar_nx_to_layered_graph(nxg: nx.Graph, remove_sl=True):
	g = graph.LayeredGraph()
	lv = nx.get_node_attributes(nxg, "hierarchy.level")
	for v in nxg.nodes:
		g.add_node(int(lv[v])+1, idx=int(v[1:]))
	for edge in nxg.edges:
		g.add_edge(int(edge[0][1:]), int(edge[1][1:]))
	g.add_anchors()
	g.relayer(remove_sl=remove_sl)
	g.y_val_setup()
	return g


def nx_to_layered_graph(nxg: nx.Graph, w, c, layer_assign=None, remove_sl=True):
	if layer_assign is None:
		return layering.create_layered_graph_from_directed_nx_graph(nxg, w, c, remove_sl=remove_sl)
	else:
		g = graph.LayeredGraph()
		for v in nxg.nodes:
			g.add_node(layer_assign[v], idx=int(v))
		for edge in nxg.edges:
			g.add_edge(int(edge[0]), int(edge[1]))
		g.add_anchors()
		g.y_val_setup()
		return g
