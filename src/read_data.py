import os.path
import pickle
import networkx.drawing.nx_pydot
from src import layering, type_conversions
import networkx as nx
import pydot


def read(filepath, w=4, c=2, layer_assignments=None):
	assert os.path.isfile(filepath), f"invalid file path '{filepath}'"
	collection = ""
	if '/' in filepath:
		if filepath[:2] == "..":
			collection = filepath[filepath.index('/') + 1:filepath.index('/', 3)]
		else:
			collection = filepath[:filepath.index('/')]
	if collection == "Rome-Lib":
		g, tv = layering.create_better_layered_graph(filepath, w, c)
	elif collection == "DAGmar":
		g = type_conversions.dagmar_nx_to_layered_graph(nx.read_graphml(filepath, node_type=str))
	elif collection == "north":
		g = type_conversions.north_nx_to_layered_graph(nx.read_graphml(filepath, node_type=str), w, c)
	elif collection == "control-flow-graphs":
		gp = pydot.graph_from_dot_file(filepath)[0]
		gnx = networkx.drawing.nx_pydot.from_pydot(gp)
		if '\\n' in gnx:
			gnx.remove_node('\\n')
		g = layering.create_layered_graph_from_directed_nx_graph(gnx, w, c)
	else:
		print("Reading graph...")
		f_ext = os.path.splitext(filepath)[1]
		if f_ext == ".graphml":
			if layer_assignments is not None:
				g = type_conversions.nx_with_separate_layerings_to_layered_graph(nx.read_graphml(filepath, node_type=str), layer_assignments)
			else:
				g = type_conversions.north_nx_to_layered_graph(nx.read_graphml(filepath, node_type=str), w, c)
		elif f_ext == ".lgbin":
			with open(filepath, 'rb') as fdb:
				g = pickle.load(fdb)
				g.names_by_layer = {}
				g.edge_names_by_layer = {}
				g.edges_by_layer = {}
		else:
			if layer_assignments is not None:
				g = layering.create_edge_list_layered_graph_given_layering(filepath, layer_assignments)
			else:
				g = layering.create_edge_list_layered_graph(filepath, w, c)
	return g
