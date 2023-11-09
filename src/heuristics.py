import sys

from src import graph, read_data, vis
import random
import copy


def __layer_sweep(g: graph.LayeredGraph, n_iter, fn):
	cr_best = g.num_edge_crossings()
	order = [0] * g.n_nodes
	for lay in g.layers.values():
		for i, nd in enumerate(lay):
			order[nd.name] = i
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


def __gansner_init_ordering(gr: graph.LayeredGraph):
	""" Uses BFS to assign initial orderings based on order discovered """
	adj = gr.get_adj_list()
	bfsq = [random.choice(gr.layers[1]).name]
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
		adj[nd.name][ln].sort(key=lambda x: order[x])
	# The next line is a generator for the median sorted order, ignoring vertices with no adjacent nodes in the previous layer
	lsort = iter(sorted((v for v in layer if len(adj[v.name][ln]) != 0), key=lambda x: adj[x.name][ln][len(adj[x.name][ln])//2] if len(adj[x.name][ln])%2==1 else (adj[x.name][ln][len(adj[x.name][ln])//2] + adj[x.name][ln][len(adj[x.name][ln])//2-1])/2))
	layer = [v if len(adj[v.name][ln]) == 0 else next(lsort) for v in layer]
	for i, node in enumerate(layer):
		node.y = i
		order[node.name] = i


def barycenter(g: graph.LayeredGraph, n_iter=20):
	""" Barycenter heuristic, K. Sugiyama et al., 1981 """
	__layer_sweep(g, n_iter, __barycenter_sort_fn)


def __barycenter_sort_fn(g, layer, order, forward):
	adj = g.get_double_adj_list()
	ln = 1 - int(forward)

	# The next line is a generator for the barycentric sorted order, ignoring vertices with no adjacent nodes in the previous layer
	lsort = iter(sorted((v for v in layer if len(adj[v.name][ln]) != 0), key=lambda nd: (sum(order[ndk] for ndk in adj[nd.name][ln]) / len(adj[nd.name][ln]))))
	layer = [v if len(adj[v.name][ln]) == 0 else next(lsort) for v in layer]
	for i, node in enumerate(layer):
		node.y = i
		order[node.name] = i


def global_sifting(g: graph.LayeredGraph, maxfails=0):
	""" Global Sifting heuristic, C. Matuszewski, 1999 """

	adj = g.get_adj_list()
	d_adj = g.get_double_adj_list()
	flat_order = [0] * len(g.nodes)  # nodeID -> order (index in layer)
	for lay in g.layers.values():
		for j, v in enumerate(lay):
			flat_order[v.name] = j
	sift_order = sorted([v.name for v in g.nodes], key=lambda x: -len(adj[x]))
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


def __sift(v, layer, ranks, d_adj, cr_num):
	best_cr = cr_num
	idx = next((i for i, nd in enumerate(layer) if nd.name == v))
	for i in range(idx, 0, -1):
		cr_num += __cr_diff(layer[i-1].name, layer[i].name, ranks, d_adj)
		layer[i-1], layer[i] = layer[i], layer[i-1]
		if cr_num < best_cr:
			best_cr = cr_num
	for i in range(0, len(layer) - 1):
		cr_num += __cr_diff(layer[i].name, layer[i + 1].name, ranks, d_adj)
		layer[i], layer[i + 1] = layer[i + 1], layer[i]
		if cr_num < best_cr:
			best_cr = cr_num
	i = len(layer) - 1
	while cr_num != best_cr:
		cr_num += __cr_diff(layer[i - 1].name, layer[i].name, ranks, d_adj)
		layer[i - 1], layer[i] = layer[i], layer[i - 1]
		i -= 1
	for j, nd in enumerate(layer):
		ranks[nd.name] = j
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
	to_select = [v.name for v in layer]
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
		min_v = sys.maxsize
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
	layer.sort(key=lambda x: new_layer_order[x.name])
	for i, node in enumerate(layer):
		node.y = i
		order[node.name] = i


def greedy_switching(g: graph.LayeredGraph, n_iter=20):
	""" Greedy Switching heuristic, P. Eades and D. Kelly, 1986 """
	__layer_sweep(g, n_iter, __switching_sort_fn)


def __switching_sort_fn(g, layer, order, forward):
	adj = g.get_double_adj_list()
	ln = 1 - int(forward)
	for i in range(len(layer) - 1):
		if __calc_if_swap_improves_1layer(order, adj, layer[i].name, layer[i + 1].name, ln):
			layer[i], layer[i + 1] = layer[i + 1], layer[i]
	for i, node in enumerate(layer):
		node.y = i
		order[node.name] = i


if __name__ == '__main__':
	graph = read_data.read("../Rome-Lib/graficon70nodi/grafo1233.70")
	# barycenter(graph)
	# median(graph)
	# global_sifting(graph)
	# weighted_median(graph)
	# greedy_insert(graph)
	greedy_switching(graph)

	# op = optimization.LayeredOptimizer(graph)
	# op.optimize_layout()
	vis.draw_graph(graph, "interim", nested=True)
	print(graph.num_edge_crossings())
