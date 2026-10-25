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
	# 	horizontal_node_symmetry=True,  # maximize node symmetry
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
	# make_up_extra_three_two({225: 96, 230: 88, 245: 98})

	graph4 = read("datasets/random_graphs/uniform_layered/n25/graph40.lgbin")
	
	opt3 = LayeredOptimizer(graph4)
	rs1 = opt3.optimize_layout(crossing_minimization=True)
	rs2 = opt3.optimize_layout(fix_x_vars=True, horizontal_edge_symmetry=True, edge_length_minimization=True)
	print("SEQ", rs1.metrics, rs2.metrics)
	print("SEQ TOTAL", rs1.metrics.crossings + rs2.metrics.edge_symmetry + rs2.metrics.length)
	draw_graph(opt3.g, "un1_seq")

	opt2 = LayeredOptimizer(graph4)
	rh1 = opt2.optimize_layout(crossing_minimization=True)
	rh2 = opt2.optimize_layout(horizontal_edge_symmetry=True, fix_x_vars=True)
	rh3 = opt2.optimize_layout(edge_length_minimization=True, fix_x_vars=True, hybrid_constraints=[("horizontal_edge_symmetry", rh2.metrics.horizontal_edge_symmetry * 1)])
	print("HYBRID", rh3.metrics)
	print("HYBRID TOTAL", rh1.metrics.crossings + rh3.metrics.edge_symmetry + rh3.metrics.length)
	draw_graph(opt2.g, "un1_hybrid")

	opt1 = LayeredOptimizer(graph4)
	rc1 = opt1.optimize_layout(crossing_minimization=True, horizontal_edge_symmetry=True, edge_length_minimization=True)
	print("COMBO", rc1.metrics)
	print("COMBO TOTAL", rc1.metrics.crossings + rc1.metrics.edge_symmetry + rc1.metrics.length)
	draw_graph(opt1.g, "un1_combo")


	# opt3 = LayeredOptimizer(graph4)
	# rs1 = opt3.optimize_layout(crossing_minimization=True, planarization=True)
	# rs2 = opt3.optimize_layout(fix_x_vars=True, edge_length_minimization=True)
	# print("SEQ", rs1.metrics, rs2.metrics)
	# print("SEQ TOTAL", rs1.metrics.crossings + rs1.metrics.planar + rs2.metrics.length)
	# draw_graph(opt3.g, "un2_seq")

	# opt2 = LayeredOptimizer(graph4)
	# rh1 = opt2.optimize_layout(crossing_minimization=True)
	# rh2 = opt2.optimize_layout(planarization=True, hybrid_constraints=[("crossings", rh1.metrics.crossings)])
	# rh3 = opt2.optimize_layout(edge_length_minimization=True, hybrid_constraints=[("crossings", rh1.metrics.crossings), ("planarization", rh2.metrics.planar)])
	# print("HYBRID", rh3.metrics)
	# print("HYBRID TOTAL", rh3.metrics.crossings + rh3.metrics.planar + rh3.metrics.length)
	# draw_graph(opt2.g, "un2_hybrid")

	# opt1 = LayeredOptimizer(graph4)
	# rc1 = opt1.optimize_layout(crossing_minimization=True, planarization=True, edge_length_minimization=True)
	# print("COMBO", rc1.metrics)
	# print("COMBO TOTAL", rc1.metrics.crossings + rc1.metrics.planar + rc1.metrics.length)
	# draw_graph(opt1.g, "un2_combo")
