from layered_optimization import LayeredOptimizer
from layered_optimization.vis import draw_graph
from layered_optimization.read_data import read
from layered_optimization.heuristics import improved_sifting, weighted_median
from layered_optimization.random_graph_generation import generate_three_to_two_dataset, make_up_extra_three_two

if __name__ == '__main__':
	""" Example use case. Optimizes and draws a 40-node Rome-Lib graph. """
	# graph = read("datasets/Rome-Lib/graficon40nodi/grafo1710.40")
	# optimizer = LayeredOptimizer(graph, name="rome-40", store_optimization_results=True)
	# print(optimizer.optimize_layout(crossing_minimization=True, streamline=True, bend_minimization=True, edge_length_minimization=True))
	# for l_e in graph.get_long_edges():
	# 	print(l_e, [graph[nd].y for nd in l_e])
	# draw_graph(graph, "rome-40-layout")

	# lit_graph = read("localfiles/thesis_lit.txt", layer_assignments=[0,0,1,1,2,2,2,3])
	# print([node.is_anchor_node for node in lit_graph])
	# optimizer = LayeredOptimizer(lit_graph, name="thesis-lit", store_optimization_results=False)
	# print(optimizer.optimize_layout(streamline=True, bend_minimization=True,
	# 								edge_length_minimization=True, anchor_proximity=0.5, hybrid_constraints=[("bend_minimization", 0), ("crossings", 0)]))
	# draw_graph(lit_graph, "thesis-lit", label_nodes=False, node_radius=20, dot_anchors=False, node_color=(222, 216, 251), node_outline=True, node_outline_color=(49, 53, 151), node_outline_width=2)

	""" Optimize using a more complex model """
	# graph2 = read("datasets/Rome-Lib/graficon30nodi/grafo1237.30")
	# optimizer = LayeredOptimizer(graph2, name="rome-30")
	# result = optimizer.optimize_layout(crossing_minimization=True)
	# optimizer.optimize_layout(
	# 	crossing_angle=True,  # also make crossing angles as close as possible to 90 degrees
	# 	symmetry_maximization=True,  # maximize node symmetry
	# 	hybrid_constraints=["crossings", result.objval.crossings]  # keep the crossing number optimal
	# )
	# draw_graph(graph2, "rome-30-layout")

	""" Optimize large graph using our large neighborhood search technique for 60 seconds. """
	# graph3 = read("datasets/random_graphs/ratio_d3/r1.5k36n24/graph7.lgbin")
	# optimizer = LayeredOptimizer(graph3, name="large-graph")
	# optimizer.optimize_layout(crossing_minimization=True, use_lns=True, cutoff_time=60, bucket_size=1000)
	# draw_graph(graph3, "large-layout")

	# graph3 = read("datasets/scale-free/n35/graph0.graphml")
	# optimizer = LayeredOptimizer(graph3)
	# optimizer.optimize_layout(crossing_minimization=True)
	# draw_graph(graph3, "scalefree")

	# generate_three_to_two_dataset()
	make_up_extra_three_two({225: 96, 230: 88, 245: 98})

	# graph4 = read("datasets/random_graphs/uniform_layered/n25/graph4.lgbin")
	# print(graph4.n_nodes, len(graph4.edges))
	# draw_graph(graph4, "random_layered")
