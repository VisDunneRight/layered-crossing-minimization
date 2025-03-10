from src.optimization import LayeredOptimizer
from src.vis import draw_graph
from src.optimization_open_src import HiGHSLayeredOptimizer
from src.neighborhood import *
from src.graph import LayeredGraph


if __name__ == '__main__':
	""" Example use case. Optimizes and draws a 40-node Rome-Lib graph with direct transitivity (default) and a selection of switches. """
	optimizer = LayeredOptimizer("Rome-Lib/graficon40nodi/grafo3216.40")
	# optimizer = LayeredOptimizer("Rome-Lib/graficon96nodi/grafo3510.96")
	# optimizer = LayeredOptimizer("random graphs/ratio_d3/r1.5k12n8/graph5.lgbin")
	optimizer.optimize_layout(crossing_minimization=True)
	# draw_graph(optimizer.g, "testing")

	gr = LayeredGraph()
	gr.add_nodes([(0, 0), (1, 1), (2, 1), (3, 1), (4, 2), (5, 2), (6, 2), (7, 3), (8, 3), (9, 3), (10, 4), (11, 4), (12, 4), (13, 5)])
	gr.add_edges([(0, 1), (0, 3), (1, 4), (2, 5), (3, 6), (4, 9), (5, 7), (6, 8), (7, 10), (8, 11), (9, 12), (10, 13), (12, 13)])
	# gr.add_edges([(0, 1), (0, 3), (1, 4), (2, 5), (3, 6), (4, 8), (5, 7), (5, 9), (6, 8), (7, 10), (8, 11), (9, 12), (10, 13), (12, 13)])
	opt = LayeredOptimizer(gr)
	opt.optimize_layout(crossing_minimization=True, bendiness_reduction=True)
	draw_graph(opt.g, "stussy")


	# x = vertical_neighborhood(optimizer.g, 5, 2500)
	# draw_graph(optimizer.g, "testing", groups=[int(v) for v in x])

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
