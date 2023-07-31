from src import graph
import random
import copy


def gansner_ordering(g: graph.LayeredGraph, n_iter):
	g.create_double_adj_list()
	order = __gansner_init_ordering(g)
	cr_best = g.num_edge_crossings()
	best = copy.deepcopy(order)
	for i in range(n_iter):
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
	for lay in best:
		for i, nd in enumerate(lay):
			g[nd].y = i


def __gansner_init_ordering(gr: graph.LayeredGraph):
	adj = gr.create_normal_adj_list()
	bfsq = [random.choice(gr.layers[1]).name]
	seen = [False] * (len(gr.nodes) + 1)
	seen[bfsq[0]] = True
	layer_seencts = [0] * (len(gr.layers) + 1)
	layer_seencts[1] += 1
	ordering = [[] for i in range(len(gr.layers) + 1)]
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
	for lnum, lay in enumerate(order):
		for nd in lay:
			mval = __gansner_median_value(gr, nd, gr[nd].layer + 2 * (c_iter % 2) - 1)  # y = 2x-1 maps 0 to -1 and 1 to 1
			gr[nd].y = mval if mval != -1 else gr[nd].y
		lay.sort(key=lambda x: gr[x].y)


def __gansner_median_value(gr: graph.LayeredGraph, v, adj_rank):
	if adj_rank > gr[v].layer:
		p = [gr[nd].y for nd in gr.double_adj_list[v][1]]
		p.sort()
	elif adj_rank < gr[v].layer:
		p = [gr[nd].y for nd in gr.double_adj_list[v][0]]
		p.sort()
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
	improved = True
	while improved:
		improved = False
		for r in range(len(rank)):
			for i in range(len(rank[r]) - 1):
				v = rank[r][i]
				w = rank[r][i+1]
				if calc_if_swap_improves(flat_rank, gr, v, w):
					improved = True
					rank[r][i], rank[r][i+1] = w, v
					flat_rank[v], flat_rank[w] = flat_rank[w], flat_rank[v]


def calc_if_swap_improves(rank, gr: graph.LayeredGraph, v, w):  # requirement: v is before w in the ranking
	n_cr_noswap = 0
	n_cr_swap = 0
	# indices = {v_adj: rank[gr[v].layer - 1].index(v_adj) for v_adj in gr.double_adj_list[v][0]}
	# indices.update({v_adj: rank[gr[v].layer + 1].index(v_adj) for v_adj in gr.double_adj_list[v][1]})
	# indices.update({w_adj: rank[gr[v].layer - 1].index(w_adj) for w_adj in gr.double_adj_list[w][0]})
	# indices.update({w_adj: rank[gr[v].layer + 1].index(w_adj) for w_adj in gr.double_adj_list[w][1]})
	for i in range(2):
		for v_adj in gr.double_adj_list[v][i]:
			for w_adj in gr.double_adj_list[w][i]:
				if v_adj != w_adj:
					if rank[v_adj] < rank[w_adj]:
						n_cr_swap += 1
					else:
						n_cr_noswap += 1
	return n_cr_swap < n_cr_noswap


def sugiyama_barycenter(g: graph.LayeredGraph, n_iter):
	g.create_double_adj_list()
	order = __gansner_init_ordering(g)
	cr_best = g.num_edge_crossings()
	best = copy.deepcopy(order)
	for i in range(n_iter):
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
	for lay in best:
		for i, nd in enumerate(lay):
			g[nd].y = i


if __name__ == '__main__':
	pass
