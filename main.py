from src.optimization import LayeredOptimizer
from src.vis import draw_graph


if __name__ == '__main__':
	""" Example use case. Optimizes and draws a 40-node Rome-Lib graph with vertical position transitivity (default) and a selection of switches. """
	# optimizer = LayeredOptimizer("Rome-Lib/graficon40nodi/grafo3216.40")
	# optimizer.symmetry_breaking = True
	# optimizer.aggro_presolve = True
	# optimizer.mip_relax = True
	# optimizer.draw_graph = True
	# optimizer.bendiness_reduction = True
	# optimizer.name = "grafo3216.40"
	# optimizer.optimize_layout()

	optimizer = LayeredOptimizer("Rome-Lib/graficon55nodi/grafo3318.55")
	# optimizer.symmetry_breaking = True
	# optimizer.aggro_presolve = True
	# optimizer.mip_relax = True
	optimizer.draw_graph = True
	optimizer.bendiness_reduction = True
	# optimizer.cycle_constraints = True
	optimizer.mirror_vars = True
	optimizer.collapse_subgraphs = True
	optimizer.name = "grafo3318.55"
	optimizer.optimize_layout()

	# optimizer = LayeredOptimizer("Rome-Lib/graficon19nodi/grafo976.19")
	# optimizer.bendiness_reduction = True
	# optimizer.optimize_layout()
	# draw_graph(optimizer.g, "no_butterfly")
	# optimizer.mirror_vars = True
	# optimizer.optimize_layout()
	# draw_graph(optimizer.g, "with_butterfly")
