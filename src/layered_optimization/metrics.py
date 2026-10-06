from layered_optimization.graph import LayeredGraph
import numpy as np


def edge_crossings(g: LayeredGraph):
	e_b_l = g.get_edges_by_layer()
	n_ec = 0
	for edge_list in e_b_l.values():
		for e1, e2 in itertools.combinations(edge_list, 2):
			if len({e1.n1.id, e1.n2.id, e2.n1.id, e2.n2.id}) == 4:
				if e1.same_layer_edge and e2.same_layer_edge:
					u1, w1, u2, w2 = e1.n1.y, e1.n2.y, e2.n1.y, e2.n2.y
					if (u1 > u2 > w1 > w2) or (u1 > w2 > w1 > u2) or (w1 > u2 > u1 > w2) or (w1 > w2 > u1 > u2) or (u2 > u1 > w2 > w1) or (u2 > w1 > w2 > u1) or (w2 > u1 > u2 > w1) or (w2 > w1 > u2 > u1):
						n_ec += 1
				elif e1.same_layer_edge:
					if (e1.n1.y > e2.n1.y > e1.n1.y) or (e1.n2.y > e2.n1.y > e1.n1.y):
						n_ec += 1
				elif e2.same_layer_edge:
					if (e2.n1.y > e1.n1.y > e2.n2.y) or (e2.n2.y > e1.n1.y > e2.n1.y):
						n_ec += 1
				elif (e1.n1.y > e2.n1.y and e1.n2.y < e2.n2.y) or (e1.n1.y < e2.n1.y and e1.n2.y > e2.n2.y):
					n_ec += 1
	return n_ec


def edge_length(g: LayeredGraph):
	t_length = 0
	for ed in g.edges:
		t_length += abs(ed.n1.y - ed.n2.y)
	return t_length


def crossing_angle(g: LayeredGraph, layer_dist=1.5, use_approximation=False):
	t_angle_radians = 0
	crossed_edges = g.edge_crossing_edges()
	for ed1, ed2 in crossed_edges:
		e1m = (g[ed1[0]].y - g[ed1[1]].y) / layer_dist
		e2m = (g[ed2[0]].y - g[ed2[1]].y) / layer_dist
		if use_approximation:
			t_angle_radians += abs(1 + e1m * e2m)
			print(e1m, e2m, ed1, ed2)
		else:
			t_angle_radians += np.pi / 2 - np.atan(abs((e1m - e2m) / (1 + e1m * e2m))) if 1 + e1m * e2m != 0 else np.pi / 2
	return t_angle_radians


def minmax_crossings(g: LayeredGraph):
	pass


def node_symmetry(g: LayeredGraph):
	pass


def edge_symmetry(g: LayeredGraph):
	pass


def bend_minimization(g: LayeredGraph, tolerance=10**-6):
	t_edge_bends = 0
	l_es = g.get_long_edges()
	for l_e in l_es:
		for i, v_1 in enumerate(l_e[:-2]):
			v_2 = l_e[i + 1]
			v_3 = l_e[i + 2]
			if abs((g[v_2].y - g[v_1].y) - (g[v_3].y - g[v_2].y)) > tolerance:
				t_edge_bends += 1
	return t_edge_bends


def fairness(g: LayeredGraph, fairness_metric="crossings"):
	fair_group_1 = [v for v in g.node_ids if g.node_data["fairness"][v] == 0]
	if fairness_metric == "crossings":
		g_cr = [0, 0]
		len_g0, len_g1 = len(fair_group_1), g.n_nodes - len(fair_group_1)
		e_b_l = g.get_edge_ids_by_layer()
		for l in range(g.n_layers - 1):
			for i, e1 in enumerate(e_b_l[l]):
				for e2 in e_b_l[l][i+1:]:
					if (g[e1[0]].y > g[e2[0]].y and g[e1[1]].y < g[e2[1]].y) or (g[e1[0]].y < g[e2[0]].y and g[e1[1]].y > g[e2[1]].y):
						fsum = sum(1 for v in [e1[0], e1[1], e2[0], e2[1] if v in fair_group_1 else 0])
						if fsum <= 1:
							g_cr[1] += 1
						elif fsum >= 3:
							g_cr[0] += 1
		return abs(g_cr[0] - (g_cr[1] * len_g0 / len_g1))
	elif fairness_metric == "edge_length":
		g_bends = [0, 0]
		len_g0, len_g1 = len(fair_group_1), g.n_nodes - len(fair_group_1)
		for ed in g.edges:
			if ed.n1.id in fair_group_1 and ed.n2.id in fair_group_1:
				g_bends[0] += abs(ed.n1.y - ed.n2.y)
			elif ed.n1.id not in fair_group_1 and ed.n2.id not in fair_group_1:
				g_bends[1] += abs(ed.n1.y - ed.n2.y)
		return abs((g_bends[0] / len_g0) - (g_bends[1] / len_g1))


def planarization(g: LayeredGraph):
	print("Minimum edge deletion cannot be verified. Run solver algorithm instead.")
