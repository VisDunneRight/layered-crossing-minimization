from layered_optimization.graph import LayeredGraph
import numpy as np
import itertools
from sys import maxsize


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


def crossing_angle(g: LayeredGraph, layer_dist=1.5, use_approximation=False, use_min_approx=False, return_angles=False, return_slopes=False):
	t_angle_radians = 0
	crossed_edges = g.edge_crossing_edges()
	angles = {}
	slopes = {}
	for ed1, ed2 in crossed_edges:
		e1m = (g[ed1[0]].y - g[ed1[1]].y) / layer_dist
		e2m = (g[ed2[0]].y - g[ed2[1]].y) / layer_dist
		if return_slopes:
			slopes[ed1], slopes[ed2] = e1m, e2m
		if use_approximation:
			t_angle_radians += abs(1 + e1m * e2m)
			# print(e1m, e2m, ed1, ed2)
			if return_angles:
				angles[ed1, ed2] = abs(1 + e1m * e2m)
		elif use_min_approx:
			t_angle_radians += min(abs(1 + e1m * e2m), np.pi / 2)
			# print(e1m, e2m, ed1, ed2)
			if return_angles:
				angles[ed1, ed2] = min(abs(1 + e1m * e2m), np.pi / 2)
		else:
			t_angle_radians += np.pi / 2 - np.arctan(abs((e1m - e2m) / (1 + e1m * e2m))) if 1 + e1m * e2m != 0 else 0
			if return_angles:
				angles[ed1, ed2] = np.pi / 2 - np.arctan(abs((e1m - e2m) / (1 + e1m * e2m))) if 1 + e1m * e2m != 0 else 0
	if return_angles:
		return angles
	if return_slopes:
		return slopes
	return t_angle_radians


def min_max_crossings(g: LayeredGraph, return_minmax_edge=False, return_counts_per_edge=False):
	cur_max_cr = 0
	counts_per_edge = {}
	e_with_max_cr = None
	long_edge_map = {}
	e_b_l = g.get_edges_by_layer()
	for l_e in g.get_long_edges():  # map (first, second) node in long edge to list of edge ids
		long_edge_map[l_e[0], l_e[1]] = []
		for i in range(0, len(l_e) - 1):
			long_edge_map[l_e[0], l_e[1]].append((l_e[i], l_e[i+1]))

	for lid in range(g.n_layers - 1):
		for e1 in e_b_l[lid]:
			n_ec = 0
			if not e1.n1.is_anchor_node and not e1.n2.is_anchor_node:
				for e2 in e_b_l[lid]:
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
				counts_per_edge[e1.n1.id, e1.n2.id] = n_ec
			elif (e1.n1.id, e1.n2.id) in long_edge_map:
				for e1long_id in long_edge_map[e1.n1.id, e1.n2.id]:
					e1l = g.get_edge(e1long_id[0], e1long_id[1])
					for e2 in e_b_l[e1l.n1.layer]:
						if len({e1l.n1.id, e1l.n2.id, e2.n1.id, e2.n2.id}) == 4:
							if e1l.same_layer_edge and e2.same_layer_edge:
								u1, w1, u2, w2 = e1l.n1.y, e1l.n2.y, e2.n1.y, e2.n2.y
								if (u1 > u2 > w1 > w2) or (u1 > w2 > w1 > u2) or (w1 > u2 > u1 > w2) or (w1 > w2 > u1 > u2) or (u2 > u1 > w2 > w1) or (u2 > w1 > w2 > u1) or (w2 > u1 > u2 > w1) or (w2 > w1 > u2 > u1):
									n_ec += 1
							elif e1l.same_layer_edge:
								if (e1l.n1.y > e2.n1.y > e1l.n1.y) or (e1l.n2.y > e2.n1.y > e1l.n1.y):
									n_ec += 1
							elif e2.same_layer_edge:
								if (e2.n1.y > e1l.n1.y > e2.n2.y) or (e2.n2.y > e1l.n1.y > e2.n1.y):
									n_ec += 1
							elif (e1l.n1.y > e2.n1.y and e1l.n2.y < e2.n2.y) or (e1l.n1.y < e2.n1.y and e1l.n2.y > e2.n2.y):
								n_ec += 1
				counts_per_edge[e1.n1.id, e1.n2.id] = n_ec

			if n_ec > cur_max_cr:
				cur_max_cr = n_ec
				e_with_max_cr = (e1.n1.id, e1.n2.id)

	if return_counts_per_edge:
		return counts_per_edge
	elif return_minmax_edge:
		return e_with_max_cr
	else:
		return cur_max_cr


def horizontal_node_symmetry(g: LayeredGraph, return_act_sym=False, return_sym_pairs_and_axis=False, h_axis=None, tolerance=0.035):
	best_tsym = maxsize
	best_act_sym = 0
	nodes_by_layer = g.get_ids_by_layer()
	best_v_axis, best_h_axis = -1, -1
	sym_pairs, best_sym_pairs = [], []

	n_pairs = [v for nlist in nodes_by_layer.values() for v in itertools.combinations_with_replacement(nlist, 2)
		if (g[v[0]].is_anchor_node and g[v[1]].is_anchor_node) or (not g[v[0]].is_anchor_node and not g[v[1]].is_anchor_node)]
	axes = [(g[n_pair[0]].y + g[n_pair[1]].y) / 2 for n_pair in n_pairs] if h_axis is None else [h_axis]
	for axis in axes:
		tsym = 0
		act_sym = 0
		sym_pairs.clear()
		for n_p in n_pairs:
			y1, y2 = g[n_p[0]].y, g[n_p[1]].y
			if abs((axis - y1) - (y2 - axis)) > 2 * tolerance:
				tsym += 1  # count num of non-symmetric nodes
			else:  # else pair is symmetric
				if return_sym_pairs_and_axis:
					sym_pairs.append(n_p)
				act_sym += 1
		if act_sym > best_act_sym:
			best_tsym = tsym
			best_act_sym = act_sym
			best_h_axis = axis
			if return_sym_pairs_and_axis:
				best_sym_pairs = sym_pairs.copy()
	if h_axis is None:
		print(f"Best H-axis: Y={best_h_axis}")

	if return_sym_pairs_and_axis:
		return best_sym_pairs, best_h_axis
	else:
		return best_tsym if not return_act_sym else best_act_sym


def vertical_node_symmetry(g: LayeredGraph, return_act_sym=False, return_sym_pairs_and_axis=False, v_axis=None, tolerance=0.035):
	best_tsym = maxsize
	best_act_sym = 0
	nodes_by_layer = g.get_ids_by_layer()
	best_v_axis, best_h_axis = -1, -1
	sym_pairs, best_sym_pairs = [], []

	axes = [v / 2 for v in range(1, 2 * len(g.layers) - 3)] if v_axis is None else [v_axis]
	for axis in axes:
		layer_pair = [round(axis - 0.001), round(axis + 0.001)]
		n_pairs = [n_pair for l_i in range(min(layer_pair[0], g.n_layers - layer_pair[1]) + 1)
			if layer_pair[0] - l_i in g.layers and layer_pair[1] + l_i in g.layers
			for n_pair in itertools.product(nodes_by_layer[layer_pair[0] - l_i], nodes_by_layer[layer_pair[1] + l_i])
			if g[n_pair[0]].is_anchor_node == g[n_pair[1]].is_anchor_node]
		n_pairs = [nv for nv in n_pairs if (nv[1], nv[0]) not in n_pairs or nv[0] == nv[1] or nv[0] < nv[1]]
		tsym = 0
		act_sym = 0
		sym_pairs.clear()
		for n_p in n_pairs:
			y1, y2 = g[n_p[0]].y, g[n_p[1]].y
			if abs(y1 - y2) > tolerance:
				tsym += 1
			else:
				if return_sym_pairs_and_axis:
					sym_pairs.append(n_p)
				act_sym += 1
		# print(axis, tsym, act_sym)
		if act_sym > best_act_sym:
			best_tsym = tsym
			best_act_sym = act_sym
			best_v_axis = axis
			if return_sym_pairs_and_axis:
				best_sym_pairs = sym_pairs.copy()
	if v_axis is None:
		print(f"Best V-axis: X={best_v_axis}")

	if return_sym_pairs_and_axis:
		return best_sym_pairs, best_v_axis
	else:
		return best_tsym if not return_act_sym else best_act_sym


def rotation_node_symmetry(g: LayeredGraph, return_act_sym=False, return_sym_pairs_and_axis=False, h_axis=None, v_axis=None, tolerance=0.035):
	best_tsym = maxsize
	best_act_sym = 0
	nodes_by_layer = g.get_ids_by_layer()
	best_v_axis, best_h_axis = -1, -1
	sym_pairs, best_sym_pairs = [], []

	v_axes = [v / 2 for v in range(1, 2 * len(g.layers) - 3)] if v_axis is None else [v_axis]
	for v_ax in v_axes:
		layer_pair = [round(v_ax - 0.001), round(v_ax + 0.001)]
		n_pairs = [n_pair for l_i in range(min(layer_pair[0], g.n_layers - layer_pair[1]) + 1)
			if layer_pair[0] - l_i in g.layers and layer_pair[1] + l_i in g.layers
			for n_pair in itertools.product(nodes_by_layer[layer_pair[0] - l_i], nodes_by_layer[layer_pair[1] + l_i])
			if g[n_pair[0]].is_anchor_node == g[n_pair[1]].is_anchor_node]
		n_pairs = [nv for nv in n_pairs if (nv[1], nv[0]) not in n_pairs or nv[0] == nv[1] or nv[0] < nv[1]]
		h_axes = [(g[n_pair[0]].y + g[n_pair[1]].y) / 2 for n_pair in n_pairs] if h_axis is None else [h_axis]
		for h_ax in h_axes:
			tsym = 0
			act_sym = 0
			sym_pairs.clear()
			for n_p in n_pairs:
				y1, y2 = g[n_p[0]].y, g[n_p[1]].y
				if abs((h_ax - y1) - (y2 - h_ax)) > 2 * tolerance:
					tsym += 1
				else:
					if return_sym_pairs_and_axis:
						sym_pairs.append(n_p)
					act_sym += 1
			if act_sym > best_act_sym:
				best_tsym = tsym
				best_act_sym = act_sym
				best_v_axis, best_h_axis = v_ax, h_ax
				if return_sym_pairs_and_axis:
					best_sym_pairs = sym_pairs.copy()
	if h_axis is None:
		print(f"Best H-axis: Y={best_h_axis}")
	if v_axis is None:
		print(f"Best V-axis: X={best_v_axis}")

	if return_sym_pairs_and_axis:
		return best_sym_pairs, best_h_axis, best_v_axis
	else:
		return best_tsym if not return_act_sym else best_act_sym


def horizontal_edge_symmetry(g: LayeredGraph, return_sym_pairs_and_axis=False, return_act_sym=False, h_axis=None, tolerance=0.035):
	edges_by_layer = g.get_edge_ids_by_layer()
	nodes_by_layer = g.get_ids_by_layer()
	best_tsym = maxsize
	best_v_axis, best_h_axis = -1, -1
	best_act_sym = 0
	sym_pairs, best_sym_pairs = [], []

	sym_e_pairs = [v for elist in edges_by_layer.values() for v in itertools.combinations(elist, 2)
		if (g[v[0][0]].is_anchor_node == g[v[1][0]].is_anchor_node) and (g[v[0][1]].is_anchor_node == g[v[1][1]].is_anchor_node)]
	n_pairs = [v for nlist in nodes_by_layer.values() for v in itertools.combinations_with_replacement(nlist, 2)
		if (g[v[0]].is_anchor_node and g[v[1]].is_anchor_node) or (not g[v[0]].is_anchor_node and not g[v[1]].is_anchor_node)]
	axes = [(g[n_pair[0]].y + g[n_pair[1]].y) / 2 for n_pair in n_pairs] if h_axis is None else [h_axis]
	for axis in axes:
		n_sym_pairs, _ = horizontal_node_symmetry(g, return_sym_pairs_and_axis=True, h_axis=axis, tolerance=tolerance)
		tsym = 0
		act_sym = 0
		sym_pairs.clear()
		for e_pair in sym_e_pairs:
			if not (((e_pair[0][0], e_pair[1][0]) in n_sym_pairs or (e_pair[1][0], e_pair[0][0]) in n_sym_pairs) and ((e_pair[0][1], e_pair[1][1]) in n_sym_pairs or (e_pair[1][1], e_pair[0][1]) in n_sym_pairs)):
				tsym += 1
			else:
				act_sym += 1
				if return_sym_pairs_and_axis:
					sym_pairs.append(e_pair)
		if act_sym > best_act_sym:
			best_tsym = tsym
			best_act_sym = act_sym
			best_h_axis = axis
			if return_sym_pairs_and_axis:
				best_sym_pairs = sym_pairs.copy()
	if h_axis is None:
		print(f"Best H-axis: Y={best_h_axis}")

	if return_sym_pairs_and_axis:
		return best_sym_pairs, best_h_axis
	return best_tsym if not return_act_sym else best_act_sym


def vertical_edge_symmetry(g: LayeredGraph, return_sym_pairs_and_axis=False, return_act_sym=False, v_axis=None, tolerance=0.035):
	edges_by_layer = g.get_edge_ids_by_layer()
	nodes_by_layer = g.get_ids_by_layer()
	best_tsym = maxsize
	best_v_axis, best_h_axis = -1, -1
	best_act_sym = 0
	sym_pairs, best_sym_pairs = [], []

	axes = [v / 2 for v in range(1, 2 * len(g.layers) - 3)] if v_axis is None else [v_axis]
	for axis in axes:
		layer_pair = [round(axis - 1 + 0.001), round(axis - 0.001)]
		sym_e_pairs = [e_pair for l_i in range(min(layer_pair[0], g.n_layers - layer_pair[1]) + 1)
			if layer_pair[0] - l_i in edges_by_layer and layer_pair[1] + l_i in edges_by_layer
			for e_pair in itertools.product(edges_by_layer[layer_pair[0] - l_i], edges_by_layer[layer_pair[1] + l_i])
			if (g[e_pair[0][0]].is_anchor_node == g[e_pair[1][1]].is_anchor_node) and (g[e_pair[0][1]].is_anchor_node == g[e_pair[1][0]].is_anchor_node)]
		sym_e_pairs = [ev for ev in sym_e_pairs if (ev[1], ev[0]) not in sym_e_pairs or ev[0] == ev[1] or ev[0][0] < ev[1][0] or (ev[0][0] == ev[1][0] and ev[0][1] < ev[1][1])]
		layer_pair = [round(axis - 0.001), round(axis + 0.001)]
		n_pairs = [n_pair for l_i in range(min(layer_pair[0], g.n_layers - layer_pair[1]) + 1)
			if layer_pair[0] - l_i in g.layers and layer_pair[1] + l_i in g.layers
			for n_pair in itertools.product(nodes_by_layer[layer_pair[0] - l_i], nodes_by_layer[layer_pair[1] + l_i])
			if g[n_pair[0]].is_anchor_node == g[n_pair[1]].is_anchor_node]
		n_sym_pairs, _ = vertical_node_symmetry(g, return_sym_pairs_and_axis=True, v_axis=axis, tolerance=tolerance)
		tsym = 0
		act_sym = 0
		sym_pairs.clear()
		for e_pair in sym_e_pairs:
			if not (((e_pair[0][0], e_pair[1][1]) in n_sym_pairs or (e_pair[1][1], e_pair[0][0]) in n_sym_pairs) and ((e_pair[0][1], e_pair[1][0]) in n_sym_pairs or (e_pair[1][0], e_pair[0][1]) in n_sym_pairs)):
				tsym += 1
			else:
				act_sym += 1
				if return_sym_pairs_and_axis:
					sym_pairs.append(e_pair)
		if act_sym > best_act_sym:
			best_tsym = tsym
			best_act_sym = act_sym
			best_v_axis = axis
			if return_sym_pairs_and_axis:
				best_sym_pairs = sym_pairs.copy()
	if v_axis is None:
		print(f"Best V-axis: X={best_v_axis}")

	if return_sym_pairs_and_axis:
		return best_sym_pairs, best_v_axis
	return best_tsym if not return_act_sym else best_act_sym


def rotation_edge_symmetry(g: LayeredGraph, return_sym_pairs_and_axis=False, return_act_sym=False, h_axis=None, v_axis=None, tolerance=0.035):
	edges_by_layer = g.get_edge_ids_by_layer()
	nodes_by_layer = g.get_ids_by_layer()
	best_tsym = maxsize
	best_v_axis, best_h_axis = -1, -1
	best_act_sym = 0
	sym_pairs, best_sym_pairs = [], []

	v_axes = [v / 2 for v in range(1, 2 * len(g.layers) - 3)] if v_axis is None else [v_axis]
	for v_ax in v_axes:
		layer_pair = [round(v_ax - 0.001), round(v_ax + 0.001)]
		n_pairs = [n_pair for l_i in range(min(layer_pair[0], g.n_layers - layer_pair[1]) + 1)
			if layer_pair[0] - l_i in g.layers and layer_pair[1] + l_i in g.layers
			for n_pair in itertools.product(nodes_by_layer[layer_pair[0] - l_i], nodes_by_layer[layer_pair[1] + l_i])
			if g[n_pair[0]].is_anchor_node == g[n_pair[1]].is_anchor_node]
		layer_pair = [round(v_ax - 1 + 0.001), round(v_ax - 0.001)]
		sym_e_pairs = [e_pair for l_i in range(min(layer_pair[0], g.n_layers - layer_pair[1]) + 1)
			if layer_pair[0] - l_i in edges_by_layer and layer_pair[1] + l_i in edges_by_layer
			for e_pair in itertools.product(edges_by_layer[layer_pair[0] - l_i], edges_by_layer[layer_pair[1] + l_i])
			if (g[e_pair[0][0]].is_anchor_node == g[e_pair[1][1]].is_anchor_node) and (g[e_pair[0][1]].is_anchor_node == g[e_pair[1][0]].is_anchor_node)]
		sym_e_pairs = [ev for ev in sym_e_pairs if (ev[1], ev[0]) not in sym_e_pairs or ev[0] == ev[1] or ev[0][0] < ev[1][0] or (ev[0][0] == ev[1][0] and ev[0][1] < ev[1][1])]
		h_axes = [(g[n_pair[0]].y + g[n_pair[1]].y) / 2 for n_pair in n_pairs] if h_axis is None else [h_axis]
		for h_ax in h_axes:
			tsym = 0
			act_sym = 0
			sym_pairs.clear()
			n_sym_pairs, _, _ = rotation_node_symmetry(g, return_sym_pairs_and_axis=True, h_axis=h_ax, v_axis=v_ax, tolerance=tolerance)
			for e_pair in sym_e_pairs:
				if not (((e_pair[0][0], e_pair[1][1]) in n_sym_pairs or (e_pair[1][1], e_pair[0][0]) in n_sym_pairs) and ((e_pair[0][1], e_pair[1][0]) in n_sym_pairs or (e_pair[1][0], e_pair[0][1]) in n_sym_pairs)):
					tsym += 1
				else:
					act_sym += 1
					if return_sym_pairs_and_axis:
						sym_pairs.append(e_pair)
			if act_sym > best_act_sym:
				best_tsym = tsym
				best_act_sym = act_sym
				best_v_axis, best_h_axis = v_ax, h_ax
				if return_sym_pairs_and_axis:
					best_sym_pairs = sym_pairs.copy()
	if h_axis is None:
		print(f"Best H-axis: Y={best_h_axis}")
	if v_axis is None:
		print(f"Best V-axis: X={best_v_axis}")

	if return_sym_pairs_and_axis:
		return best_sym_pairs, best_h_axis, best_v_axis
	return best_tsym if not return_act_sym else best_act_sym


def bend_minimization(g: LayeredGraph, return_bends=False, tolerance=10**-6):
	t_edge_bends = 0
	l_es = g.get_long_edges()
	bends = {}
	for l_e in l_es:
		for i, v_1 in enumerate(l_e[:-2]):
			v_2 = l_e[i + 1]
			v_3 = l_e[i + 2]
			if abs((g[v_2].y - g[v_1].y) - (g[v_3].y - g[v_2].y)) > tolerance:
				t_edge_bends += 1
				if return_bends:
					bends[(v_1, v_2, v_3)] = 1
			elif return_bends:
				bends[(v_1, v_2, v_3)] = 0
	if return_bends:
		return bends
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


def planarization(g: LayeredGraph, return_counts_per_edge=True):
	print("Warning: Planarization is unverifiable. Run optimization algorithm for maximal planar edge set.\nThis function can only be used to return the total number of crossings per edge in the input graph.")
	counts_per_edge = {(ed.n1.id, ed.n2.id): 0 for ed in g.edges if not e.n1.is_anchor_node}
	e_by_layer = g.get_edges_by_layer()
	for ed in counts_per_edge:
		if g[ed[1]].is_anchor_node:  # long edge
			for i in range(0, len(l_e) - 1):
				e1 = g.get_edge(ed[i], ed[i + 1])
				for e2 in e_by_layer[e1.n1.layer]:
					if len({e1.n1.id, e1.n2.id, e2.n1.id, e2.n2.id}) == 4:
						if e1.same_layer_edge and e2.same_layer_edge:
							u1, w1, u2, w2 = e1.n1.y, e1.n2.y, e2.n1.y, e2.n2.y
							if (u1 > u2 > w1 > w2) or (u1 > w2 > w1 > u2) or (w1 > u2 > u1 > w2) or (w1 > w2 > u1 > u2) or (u2 > u1 > w2 > w1) or (u2 > w1 > w2 > u1) or (w2 > u1 > u2 > w1) or (w2 > w1 > u2 > u1):
								counts_per_edge[ed] += 1
						elif e1.same_layer_edge:
							if (e1.n1.y > e2.n1.y > e1.n1.y) or (e1.n2.y > e2.n1.y > e1.n1.y):
								counts_per_edge[ed] += 1
						elif e2.same_layer_edge:
							if (e2.n1.y > e1.n1.y > e2.n2.y) or (e2.n2.y > e1.n1.y > e2.n1.y):
								counts_per_edge[ed] += 1
						elif (e1.n1.y > e2.n1.y and e1.n2.y < e2.n2.y) or (e1.n1.y < e2.n1.y and e1.n2.y > e2.n2.y):
							counts_per_edge[ed] += 1
		else:
			e1 = g.get_edge(ed[0], ed[1])
			for e2 in e_by_layer[e1.n1.layer]:
				if len({e1.n1.id, e1.n2.id, e2.n1.id, e2.n2.id}) == 4:
					if e1.same_layer_edge and e2.same_layer_edge:
						u1, w1, u2, w2 = e1.n1.y, e1.n2.y, e2.n1.y, e2.n2.y
						if (u1 > u2 > w1 > w2) or (u1 > w2 > w1 > u2) or (w1 > u2 > u1 > w2) or (w1 > w2 > u1 > u2) or (u2 > u1 > w2 > w1) or (u2 > w1 > w2 > u1) or (w2 > u1 > u2 > w1) or (w2 > w1 > u2 > u1):
							counts_per_edge[ed] += 1
					elif e1.same_layer_edge:
						if (e1.n1.y > e2.n1.y > e1.n1.y) or (e1.n2.y > e2.n1.y > e1.n1.y):
							counts_per_edge[ed] += 1
					elif e2.same_layer_edge:
						if (e2.n1.y > e1.n1.y > e2.n2.y) or (e2.n2.y > e1.n1.y > e2.n1.y):
							counts_per_edge[ed] += 1
					elif (e1.n1.y > e2.n1.y and e1.n2.y < e2.n2.y) or (e1.n1.y < e2.n1.y and e1.n2.y > e2.n2.y):
						counts_per_edge[ed] += 1
	
	if return_counts_per_edge:
		return counts_per_edge
	return None


def planarization(g: LayeredGraph):
	raise Exception("Maximal planar subgraph is unverifiable. Run optimization algorithm for maximal planar edge set.")
