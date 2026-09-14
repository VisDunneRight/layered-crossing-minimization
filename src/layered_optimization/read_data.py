import os.path
import pickle
import networkx.drawing.nx_pydot
from layered_optimization import layering, type_conversions
import networkx as nx
import pydot
import inspect
from pathlib import Path


def read(filepath: str, w=4, c=2, layer_assignments=None, remove_sl=True):
	path = Path(filepath)
	if path.is_absolute():
		fp = path
	else:
		# 1. Inspect the call stack to find the file that called this function
		# frame[1] represents the immediate caller of read_file()
		caller_frame = inspect.stack()[1]
		caller_file = caller_frame.filename

		# 2. Get the directory of that caller script
		caller_dir = Path(caller_file).parent.resolve()
		
		# 3. Combine the caller's directory with the relative path input
		fp = (caller_dir / path).resolve()

	assert fp.is_file(), f"invalid file path '{fp}'"
	collection = ""
	if "Rome-Lib" in filepath:
		g, tv = layering.create_better_layered_graph(fp, w, c, remove_sl=remove_sl)
	elif "DAGmar" in filepath:
		g = type_conversions.dagmar_nx_to_layered_graph(nx.read_graphml(fp, node_type=str), remove_sl=remove_sl)
	elif "north" in filepath:
		g = type_conversions.north_nx_to_layered_graph(nx.read_graphml(fp, node_type=str), w, c, remove_sl=remove_sl)
	elif "control-flow-graphs" in filepath:
		gp = pydot.graph_from_dot_file(fp)[0]
		gnx = networkx.drawing.nx_pydot.from_pydot(gp)
		if '\\n' in gnx:
			gnx.remove_node('\\n')
		g = layering.create_layered_graph_from_directed_nx_graph(gnx, w, c, remove_sl=remove_sl)
	else:
		print("Reading graph... ", end="")
		f_ext = os.path.splitext(filepath)[1]
		if f_ext == ".graphml":
			if layer_assignments is not None:
				g = type_conversions.nx_with_separate_layerings_to_layered_graph(nx.read_graphml(fp, node_type=str), layer_assignments)
			else:
				g = type_conversions.north_nx_to_layered_graph(nx.read_graphml(fp, node_type=str), w, c, remove_sl=remove_sl)
		elif f_ext == ".lgbin":
			with open(fp, 'rb') as fdb:
				g = pickle.load(fdb)
		else:
			if layer_assignments is not None:
				g = layering.create_edge_list_layered_graph_given_layering(fp, layer_assignments)
			else:
				g, _ = layering.create_edge_list_layered_graph(fp, w, c, remove_sl=remove_sl, remove_disconnected_nodes="networkx" not in fp)
		if min(g.layers) != 0:
			g.relayer()
		print("done")
	return g
