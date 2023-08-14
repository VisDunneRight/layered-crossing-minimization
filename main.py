from src.optimization import LayeredOptimizer
from src.vis import draw_graph
from src.optimization_open_src import HiGHSLayeredOptimizer


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

	# optimizer = LayeredOptimizer("Rome-Lib/graficon55nodi/grafo3318.55")

	# optimizer = HiGHSLayeredOptimizer("Rome-Lib/graficon40nodi/grafo3216.40")
	optimizer = HiGHSLayeredOptimizer("Rome-Lib/graficon36nodi/grafo11126.36")
	# optimizer = HiGHSLayeredOptimizer("Rome-Lib/graficon12nodi/grafo193.12")

	# optimizer.direct_transitivity = True

	# optimizer = LayeredOptimizer("random graphs/fixed_density/k11/graph0.lgbin")
	# optimizer = LayeredOptimizer("Rome-Lib/graficon35nodi/grafo6035.35")
	optimizer.symmetry_breaking = True
	# optimizer.aggro_presolve = True
	# optimizer.mip_relax = True
	optimizer.draw_graph = True
	optimizer.bendiness_reduction = True
	optimizer.cycle_constraints = True
	optimizer.butterfly_reduction = True
	optimizer.mirror_vars = True
	optimizer.collapse_subgraphs = True
	# optimizer.heuristic_start = True
	# optimizer.xvar_branch_priority = True
	optimizer.cutoff_time = 60
	# optimizer.name = "grafo3318.55"
	optimizer.optimize_layout()

	# optimizer = LayeredOptimizer("Rome-Lib/graficon19nodi/grafo976.19")
	# optimizer.bendiness_reduction = True
	# optimizer.optimize_layout()
	# draw_graph(optimizer.g, "no_butterfly")
	# optimizer.mirror_vars = True
	# optimizer.optimize_layout()
	# draw_graph(optimizer.g, "with_butterfly")
