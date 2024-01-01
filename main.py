from src.optimization import LayeredOptimizer
from src.vis import draw_graph
from src.optimization_open_src import HiGHSLayeredOptimizer


if __name__ == '__main__':
	""" Example use case. Optimizes and draws a 40-node Rome-Lib graph with direct transitivity (default) and a selection of switches. """
	optimizer = LayeredOptimizer("Rome-Lib/graficon40nodi/grafo3216.40")
	optimizer.locally_optimal_heuristic = True
	optimizer.symmetry_breaking = True
	# optimizer.mip_relax = True
	optimizer.bendiness_reduction = True
	optimizer.n_partitions = 2
	# optimizer.collapse_leaves = True
	optimizer.draw_graph = True
	optimizer.cutoff_time = 60
	optimizer.name = "grafo3216.40"
	optimizer.optimize_layout()

	# gr = read_data.read("Rome-Lib/graficon12nodi/grafo193.12")
	# optimizer.optimize_layout_standard(graph_arg=gr)

	""" Example use case with the open source solver HiGHS """
	# optimizer = HiGHSLayeredOptimizer("Rome-Lib/graficon12nodi/grafo193.12")
	# optimizer.vertical_transitivity = True
	# optimizer.symmetry_breaking = True
	# optimizer.cycle_constraints = True
	# optimizer.mirror_vars = True
	# optimizer.collapse_leaves = True
	# optimizer.bendiness_reduction = True
	# optimizer.draw_graph = True
	# optimizer.cutoff_time = 60
	# optimizer.name = "grafo193.12"
	# optimizer.optimize_layout()
