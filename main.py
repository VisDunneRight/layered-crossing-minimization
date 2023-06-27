from src.optimization import LayeredOptimizer


if __name__ == '__main__':
	""" Example use case. Optimizes and draws a 40-node Rome-Lib graph with vertical position __transitivity (default) and a selection of switches. """
	optimizer = LayeredOptimizer("Rome-Lib/graficon40nodi/grafo3216.40")
	optimizer.fix_one_var = True
	optimizer.aggro_presolve = True
	optimizer.mip_relax = True
	optimizer.draw_graph = True
	optimizer.bendiness_reduction = True
	optimizer.sequential_bendiness = True
	optimizer.name = "grafo3216.40"
	optimizer.optimize_layout()
