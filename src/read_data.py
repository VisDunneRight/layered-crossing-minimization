from src import layering, type_conversions
import networkx as nx


def read(filepath):
	collection = filepath[:filepath.index('/')]
	if collection == "Rome-Lib":
		g, tv = layering.create_better_layered_graph(filepath[filepath.index('/') + 1:], 4, 2)
	elif collection == "DAGmar":
		g = type_conversions.dagmar_nx_to_layered_graph(nx.read_graphml(filepath, node_type=str))
	elif collection == "north":
		g = type_conversions.north_nx_to_layered_graph(nx.read_graphml(filepath, node_type=str))
	else:
		print("Invalid path")
		return
	return g
