from sys import maxsize
from src import graph, read_data, vis
import random
import copy


def __layer_sweep(g: graph.LayeredGraph, n_iter, fn):
	cr_best = g.num_edge_crossings()
	order = [0] * g.n_nodes
	for lay in g.layers.values():
		for i, nd in enumerate(lay):
			order[nd.id] = i
	best = copy.deepcopy(order)
	for i in range(n_iter):
		if i % 2 == 0:
			for j in range(2, g.n_layers + 1):
				fn(g, g.layers[j], order, True)
		else:
			for j in range(g.n_layers - 1, 0, -1):
				fn(g, g.layers[j], order, False)
		cr_order = g.num_edge_crossings()
		print(cr_order)
		if cr_order < cr_best:
			cr_best = cr_order
			best = copy.deepcopy(order)
	for nd, yv in enumerate(best):
		g[nd].y = yv
	for lay in g.layers.values():
		lay.sort(key=lambda node: node.y)


def weighted_median(g: graph.LayeredGraph, n_iter):  # Gansner et al. suggest using n_iter=24
	""" The ordering() function defined by Gansner et al., also called the "Weighted Median Heuristic" """

	order = __gansner_init_ordering(g)
	cr_best = g.num_edge_crossings()
	best = copy.deepcopy(order)
	for i in range(n_iter):  # Even i is a forward pass, odd i is a reverse pass
		__gansner_wmedian(g, order, i)
		flat_order = [0] * (len(g.nodes) + 1)
		for lay in order:
			for j, v in enumerate(lay):
				flat_order[v] = j
		__gansner_transpose(g, order, flat_order)
		for lay in order:
			for j, nd in enumerate(lay):
				g[nd].y = j
		cr_order = g.num_edge_crossings()
		if cr_order < cr_best:
			cr_best = cr_order
			best = copy.deepcopy(order)
		print(cr_order)
	for lay in best:
		for i, nd in enumerate(lay):
			g[nd].y = i
	for lay in g.layers.values():
		lay.sort(key=lambda node: node.y)


def __gansner_init_ordering(gr: graph.LayeredGraph):
	""" Uses BFS to assign initial orderings based on order discovered """
	adj = gr.get_adj_list()
	bfsq = [random.choice(gr.layers[1]).id]
	seen = [False] * (len(gr.nodes) + 1)
	seen[bfsq[0]] = True
	layer_seencts = [0] * (len(gr.layers) + 1)
	layer_seencts[1] += 1
	ordering = [[] for _ in range(len(gr.layers) + 1)]
	ordering[1].append(bfsq[0])
	while bfsq:
		next_layer = []
		for cur in bfsq:
			for to_explore in adj[cur]:
				if not seen[to_explore]:
					seen[to_explore] = True
					next_layer.append(to_explore)
					ordering[gr[to_explore].layer].append(to_explore)
					gr[to_explore].y = layer_seencts[gr[to_explore].layer]
					layer_seencts[gr[to_explore].layer] += 1
		bfsq = next_layer.copy()
	return ordering


def __gansner_wmedian(gr: graph.LayeredGraph, order, c_iter):
	""" weighted median calculation proposed by Gansner et al. """
	for lnum, lay in enumerate(order):
		for nd in lay:
			mval = __gansner_median_value(gr, nd, gr[nd].layer + 2 * (c_iter % 2) - 1)  # y = 2x-1 maps 0 to -1 and 1 to 1
			gr[nd].y = mval if mval != -1 else gr[nd].y
		lay.sort(key=lambda x: gr[x].y)


def __gansner_median_value(gr: graph.LayeredGraph, v, adj_rank):
	if adj_rank > gr[v].layer:
		p = [gr[nd].y for nd in gr.get_double_adj_list()[v][1]]
		p.sort()
	elif adj_rank < gr[v].layer:
		p = [gr[nd].y for nd in gr.get_double_adj_list()[v][0]]
		p.sort()
	else:
		raise Exception("incorrect adjacent layer calculation")
	if len(p) == 0:
		return -1
	elif len(p) % 2 == 1:
		return p[len(p)//2]
	elif len(p) == 2:
		return sum(p) / 2
	else:
		left = p[len(p)//2 - 1] - p[0]
		right = p[len(p) - 1] - p[len(p)//2]
		return (p[len(p)//2 - 1] * right + p[len(p)//2] * left) / (left + right)


def __gansner_transpose(gr: graph.LayeredGraph, rank, flat_rank):
	"""
	Additional refinement heuristic, one of 2 improvements on the basic median heuristic
	Essentially just the greedy switching heuristic
	"""
	double_adj = gr.get_double_adj_list()
	improved = True
	while improved:
		improved = False
		for r in range(len(rank)):
			for i in range(len(rank[r]) - 1):
				v = rank[r][i]
				w = rank[r][i+1]
				if __calc_if_swap_improves(flat_rank, double_adj, v, w):
					improved = True
					rank[r][i], rank[r][i+1] = w, v
					flat_rank[v], flat_rank[w] = flat_rank[w], flat_rank[v]
					break


def __calc_if_swap_improves(rank, d_adj, v, w):  # requirement: v is before w in the ranking
	n_cr_noswap = 0
	n_cr_swap = 0
	for i in range(2):
		for v_adj in d_adj[v][i]:
			for w_adj in d_adj[w][i]:
				if v_adj != w_adj:
					if rank[v_adj] < rank[w_adj]:
						n_cr_swap += 1
					else:
						n_cr_noswap += 1
	return n_cr_swap < n_cr_noswap


def __calc_if_swap_improves_1layer(rank, d_adj, v, w, ln):  # requirement: v is before w in the ranking
	n_cr_noswap = 0
	n_cr_swap = 0
	for v_adj in d_adj[v][ln]:
		for w_adj in d_adj[w][ln]:
			if v_adj != w_adj:
				if rank[v_adj] < rank[w_adj]:
					n_cr_swap += 1
				else:
					n_cr_noswap += 1
	return n_cr_swap < n_cr_noswap


def median(g: graph.LayeredGraph, n_iter=20):
	""" Original Median heuristic plus layer-by-layer sweep, P. Eades and N. C. Wormald, 1986 """
	__layer_sweep(g, n_iter, __median_sort_fn)


def __median_sort_fn(g, layer, order, forward):
	adj = g.get_double_adj_list()
	ln = 1 - int(forward)

	for nd in layer:  # first resort adjacencies by their index in the order
		adj[nd.id][ln].sort(key=lambda x: order[x])

	layer.sort(key=lambda x: order[x.id] if len(adj[x.id][ln]) == 0 else (order[adj[x.id][ln][len(adj[x.id][ln]) // 2]] if len(adj[x.id][ln]) % 2 == 1 else order[adj[x.id][ln][len(adj[x.id][ln]) // 2 - 1]]))

	for i, node in enumerate(layer):
		node.y = i
		order[node.id] = i


def barycenter(g: graph.LayeredGraph, n_iter=20):
	""" Barycenter heuristic, K. Sugiyama et al., 1981 """
	__layer_sweep(g, n_iter, __barycenter_sort_fn)


def __barycenter_sort_fn(g, layer, order, forward):
	adj = g.get_double_adj_list()
	ln = 1 - int(forward)

	# The next line is a generator for the barycentric sorted order, ignoring vertices with no adjacent nodes in the previous layer
	lsort = iter(sorted((v for v in layer if len(adj[v.id][ln]) != 0), key=lambda nd: (sum(order[ndk] for ndk in adj[nd.id][ln]) / len(adj[nd.id][ln]))))
	layer = [v if len(adj[v.id][ln]) == 0 else next(lsort) for v in layer]
	for i, node in enumerate(layer):
		node.y = i
		order[node.id] = i


def global_sifting(g: graph.LayeredGraph, maxfails=0):
	""" Global Sifting heuristic, C. Matuszewski, 1999 """

	adj = g.get_adj_list()
	d_adj = g.get_double_adj_list()
	flat_order = [0] * len(g.nodes)  # nodeID -> order (index in layer)
	for lay in g.layers.values():
		for j, v in enumerate(lay):
			flat_order[v.id] = j
	sift_order = sorted([v.id for v in g.nodes], key=lambda x: -len(adj[x]))
	n_cr = g.num_edge_crossings()
	best_cr = n_cr
	fails = 0
	while fails <= maxfails:
		for nd in sift_order:
			n_cr = __sift(nd, g.layers[g.node_names[nd].layer], flat_order, d_adj, n_cr)
		if n_cr < best_cr:
			best_cr = n_cr
		else:
			fails += 1
			sift_order = sift_order[::-1]
		for nd in sift_order:
			n_cr = __sift(nd, g.layers[g.node_names[nd].layer], flat_order, d_adj, n_cr)
		if n_cr < best_cr:
			best_cr = n_cr
		else:
			fails += 1
		sift_order = sift_order[::-1]
	for nd, yv in enumerate(flat_order):
		g[nd].y = yv
	print(best_cr)


def __sift(v, layer, ranks, d_adj, cr_num):
	best_cr = cr_num
	idx = next((i for i, nd in enumerate(layer) if nd.id == v))
	for i in range(idx, 0, -1):
		cr_num += __cr_diff(layer[i-1].id, layer[i].id, ranks, d_adj)
		layer[i-1], layer[i] = layer[i], layer[i-1]
		if cr_num < best_cr:
			best_cr = cr_num
	for i in range(0, len(layer) - 1):
		cr_num += __cr_diff(layer[i].id, layer[i + 1].id, ranks, d_adj)
		layer[i], layer[i + 1] = layer[i + 1], layer[i]
		if cr_num < best_cr:
			best_cr = cr_num
	i = len(layer) - 1
	while cr_num != best_cr:
		cr_num += __cr_diff(layer[i - 1].id, layer[i].id, ranks, d_adj)
		layer[i - 1], layer[i] = layer[i], layer[i - 1]
		i -= 1
	for j, nd in enumerate(layer):
		ranks[nd.id] = j
	return best_cr


def __cr_diff(v, w, rank, d_adj):  # returns change to CRnum after swapping v, w
	diff = 0
	for i in range(2):
		for v_adj in d_adj[v][i]:
			for w_adj in d_adj[w][i]:
				if v_adj != w_adj:
					if rank[v_adj] < rank[w_adj]:
						diff += 1
					else:
						diff -= 1
	return diff


def greedy_insert(g: graph.LayeredGraph, n_iter=20):
	""" Greedy-Insert heuristic, P. Eades and D. Kelly, 1986 """
	__layer_sweep(g, n_iter, __insert_sort_fn)


def __insert_sort_fn(g, layer, order, forward):
	adj = g.get_double_adj_list()
	ln = 1 - int(forward)
	new_layer_order = {}
	to_select = [v.id for v in layer]
	sum_c = {nd: 0 for nd in to_select}
	c_matrix = [[0] * len(layer) for _ in range(len(layer))]
	for i, v in enumerate(to_select):
		for j, w in enumerate(to_select):
			if v != w:
				c_matrix[i][j] = sum((1 for ed1 in adj[v][ln] for ed2 in adj[w][ln] if order[ed1] > order[ed2]))
	for v in to_select:
		sum_c[v] = sum(c_matrix[order[v]])
	idx = 0
	while len(new_layer_order) < len(layer):
		min_v = maxsize
		min_nd = -1
		for v in sum_c:
			if sum_c[v] < min_v:
				min_v = sum_c[v]
				min_nd = v
		new_layer_order[min_nd] = idx
		del sum_c[min_nd]
		for v in sum_c:
			sum_c[v] -= c_matrix[order[v]][order[min_nd]]
		idx += 1
	layer.sort(key=lambda x: new_layer_order[x.id])
	for i, node in enumerate(layer):
		node.y = i
		order[node.id] = i


def greedy_switching(g: graph.LayeredGraph, n_iter=20):
	""" Greedy Switching heuristic, P. Eades and D. Kelly, 1986 """
	__layer_sweep(g, n_iter, __switching_sort_fn)


def __switching_sort_fn(g, layer, order, forward):
	adj = g.get_double_adj_list()
	ln = 1 - int(forward)
	for i in range(len(layer) - 1):
		if __calc_if_swap_improves_1layer(order, adj, layer[i].id, layer[i + 1].id, ln):
			layer[i], layer[i + 1] = layer[i + 1], layer[i]
	for i, node in enumerate(layer):
		node.y = i
		order[node.id] = i


def split(g: graph.LayeredGraph, n_iter=20):
	""" Split heuristic, P. Eades and D. Kelly, 1986 """
	__layer_sweep(g, n_iter, __split_sort_fn)


def __split_sort_fn(g, layer, order, forward):
	adj = g.get_double_adj_list()
	ln = 1 - int(forward)
	c_matrix = [[0] * len(layer) for _ in range(len(layer))]
	for i, v in enumerate(layer):
		for j, w in enumerate(layer):
			if v.id != w.id:
				c_matrix[i][j] = sum((1 for ed1 in adj[v.id][ln] for ed2 in adj[w.id][ln] if order[ed1] > order[ed2]))
	new_order = [v.id for v in layer]
	node_to_idx = {v.id: k for k, v in enumerate(layer)}

	def split_recur(list_to_order, idx_start):
		if len(list_to_order) > 1:
			pivot = 0
			low = 0
			high = len(list_to_order) - 1
			p_low = []
			p_high = []
			for id_w in range(1, len(list_to_order)):
				if c_matrix[node_to_idx[list_to_order[id_w]]][node_to_idx[list_to_order[pivot]]] < c_matrix[node_to_idx[list_to_order[pivot]]][node_to_idx[list_to_order[id_w]]]:
					p_low.append(list_to_order[id_w])
					low += 1
				else:
					p_high.append(list_to_order[id_w])
					high -= 1
			idx = idx_start
			for l1 in p_low:
				new_order[idx] = l1
				idx += 1
			new_order[idx] = list_to_order[pivot]
			idx += 1
			for h1 in p_high:
				new_order[idx] = h1
				idx += 1
			split_recur(p_low, idx_start)
			split_recur(p_high, idx_start + len(p_low) + 1)

	split_recur([vid for vid in new_order], 0)
	n_ord_rev = {v: k for k, v in enumerate(new_order)}
	layer.sort(key=lambda x: n_ord_rev[x.id])
	for i, node in enumerate(layer):
		node.y = i
		order[node.id] = i


def degree_weighted_barycenter(g: graph.LayeredGraph, threshold=0.05):
	"""
	Degree-weighted barycenter, P. Eades, X. Lin, and R. Tamassia, 1996.
	Note: This algorithm has been modified to return the result with best crossing number seen over all iterations
	"""
	adj = g.get_double_adj_list()
	converged = False
	best_cr = g.num_edge_crossings()
	best_yvals = [node.y for node in g.nodes]
	while not converged:
		converged = True
		for i in range(2, g.n_layers):
			for nd in g.layers[i]:
				old_y = nd.y
				if len(adj[nd.id][1]) > 0 and len(adj[nd.id][0]) > 0:
					nd.y = sum(g[nd_adj].y for nd_adj in adj[nd.id][0]) / (2 * len(adj[nd.id][0])) + sum(g[nd_adj].y for nd_adj in adj[nd.id][1]) / (2 * len(adj[nd.id][1]))
				elif len(adj[nd.id][1]) > 0:
					nd.y = sum(g[nd_adj].y for nd_adj in adj[nd.id][1]) / len(adj[nd.id][1])
				elif len(adj[nd.id][0]) > 0:
					nd.y = sum(g[nd_adj].y for nd_adj in adj[nd.id][0]) / len(adj[nd.id][0])
				if abs(nd.y - old_y) > threshold:
					converged = False
		cur_cr = g.num_edge_crossings()
		print(cur_cr)
		if cur_cr < best_cr:
			best_cr = cur_cr
			best_yvals = [node.y for node in g.nodes]
	for k, node in enumerate(g.nodes):
		node.y = best_yvals[k]
	for lay in g.layers.values():
		lay.sort(key=lambda node: node.y)
		for i in range(len(lay) - 1):
			if lay[i].y == lay[i + 1].y:
				if i > 0 and abs(lay[i].y - lay[i - 1].y) > threshold/2:
					lay[i].y -= threshold/2
				elif i < len(lay) - 2 and abs(lay[i + 2].y - lay[i + 1].y) > threshold/2:
					lay[i + 1].y += threshold/2
				else:
					print("Two nodes converged to same position:", lay[i], lay[i+1])
	print(g.num_edge_crossings())


def switching_with_preprocessing(g: graph.LayeredGraph, n_iter=20):
	""" Greedy Switching with BC preprocessing, suggested by E. Makinen, 1990 """
	barycenter(g, n_iter=n_iter)
	print("here", g.num_edge_crossings())
	greedy_switching(g, n_iter=n_iter)


if __name__ == '__main__':
	graph = read_data.read("../Rome-Lib/graficon70nodi/grafo1233.70")
	# barycenter(graph)
	# median(graph)
	# global_sifting(graph)
	# weighted_median(graph)
	# greedy_insert(graph)
	# greedy_switching(graph)
	# split(graph)
	# switching_with_preprocessing(graph)
	degree_weighted_barycenter(graph)

	# op = optimization.LayeredOptimizer(graph)
	# op.draw_graph = True
	# op.optimize_layout()
	vis.draw_graph(graph, "interim", nested=True)
	print(graph.num_edge_crossings())
