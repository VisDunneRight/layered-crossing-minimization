from src.optimization import LayeredOptimizer
from src.vis import draw_graph
from src.optimization_open_src import HiGHSLayeredOptimizer
from src.neighborhood import *
from src.graph import LayeredGraph


if __name__ == '__main__':
	""" All CINDER metrics """


	""" Example use case. Optimizes and draws a 40-node Rome-Lib graph with direct transitivity (default) and a selection of switches. """
	optimizer = HiGHSLayeredOptimizer("Rome-Lib/graficon40nodi/grafo3216.40")
	# optimizer = LayeredOptimizer("Rome-Lib/graficon96nodi/grafo3510.96")
	# optimizer = LayeredOptimizer("random graphs/ratio_d3/r1.5k12n8/graph5.lgbin")
	optimizer.optimize_layout()
	draw_graph(optimizer.g, "testing")
