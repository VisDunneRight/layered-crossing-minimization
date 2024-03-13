import random

from src.graph import *
from collections import defaultdict
import heapq
import networkx as nx


def simple_vertical_partition(g: LayeredGraph, cutoff):
    partition = [0] * g.n_nodes
    cur_id = 0
    acc_cv = 0
    adj = g.get_adj_list()
    # layer_edge_cts = defaultdict(int)
    layer_edge_cts = [len(g.get_edges_by_layer()[i]) for i in range(g.n_layers - 1)]
    for i in range(g.n_layers):
        for nd in g.layers[i]:
            for nd_adj in adj[nd.id]:
                if partition[nd_adj] == cur_id:  # C-VAR CALCULATION
                    acc_cv += 2 * layer_edge_cts[nd.layer]
                # if partition[nd_adj] == cur_id:  # C-VAR CALCULATION
                #     acc_cv += 2 * layer_edge_cts[nd.layer]
                #     layer_edge_cts[nd.layer] += 1
            if acc_cv > cutoff:
                cur_id += 1
                acc_cv = 0
                layer_edge_cts.clear()
            partition[nd.id] = cur_id
    return partition


def bfs_neighborhood(g: LayeredGraph, candidate_id: int, cutoff, nbhd_width=0):
    adj = g.get_adj_list()
    selected = [False] * g.n_nodes
    seen = [False] * g.n_nodes
    seen[candidate_id] = True
    if nbhd_width == 0:
        layer_edge_cts = [len(g.get_edges_by_layer()[i]) for i in range(g.n_layers - 1)]
    else:
        layer_edge_cts = find_vertical_nbhd_and_edge_counts(g, candidate_id, nbhd_width)
    bfsq = [candidate_id]
    acc_cv = 0
    while acc_cv <= cutoff and bfsq:
        cur_layer = bfsq.copy()
        bfsq.clear()
        for nd in cur_layer:
            if nbhd_width == 0 or g[nd].v_nbhd:
                for nd_adj in adj[nd]:
                    if selected[nd_adj]:
                        if g[nd_adj].layer < g[nd].layer:  # C-VAR CALCULATION
                            acc_cv += 2 * layer_edge_cts[g[nd_adj].layer]
                            layer_edge_cts[g[nd_adj].layer] -= 1
                        else:
                            acc_cv += 2 * layer_edge_cts[g[nd].layer]
                            layer_edge_cts[g[nd].layer] -= 1
                    elif not seen[nd_adj]:
                        bfsq.append(nd_adj)
                        seen[nd_adj] = True
                # print(acc_cv)
                if acc_cv <= cutoff:
                    selected[nd] = True
                else:
                    break
    return selected


def vertical_neighborhood(g: LayeredGraph, candidate_node: int, cutoff):
    """ greedy criterion: largest (incoming edges / (outgoing edges + 0.5)) """
    d_adj = g.get_double_adj_list()
    selected = [False] * g.n_nodes
    for n in g.layers[g[candidate_node].layer]:
        selected[n.id] = True
    layer_edge_cts = [len(g.get_edges_by_layer()[i]) for i in range(g.n_layers - 1)]
    layer_node_cts = defaultdict(int)
    layer_node_cts[g[candidate_node].layer] = len(g.layers[g[candidate_node].layer])
    left_layer = g[candidate_node].layer - 1
    right_layer = g[candidate_node].layer + 1
    pqueue = []
    acc_cv = 0
    exp_left = 1
    exp_right = 1
    open_left = True
    open_right = True
    if right_layer in g.layers:
        for nd in g.layers[right_layer]:
            heapq.heappush(pqueue, (-(len(d_adj[nd.id][0]) / (len(d_adj[nd.id][1]) + 0.5)), nd.id))
    if left_layer in g.layers:
        for nd in g.layers[left_layer]:
            heapq.heappush(pqueue, (-(len(d_adj[nd.id][1]) / (len(d_adj[nd.id][0]) + 0.5)), nd.id))
    while pqueue:
        next_nd = heapq.heappop(pqueue)[1]
        next_layer = g[next_nd].layer
        nd_lr = 0 if next_layer == right_layer else 1
        # print(d_adj[next_nd][nd_lr], nd_lr, next_layer, next_nd, right_layer, g[d_adj[next_nd][nd_lr][0]].layer - 1 + nd_lr)
        # print(layer_edge_cts[g[d_adj[next_nd][nd_lr][0]].layer - 1 + nd_lr])
        for i in range(2):  # C-VAR CALCULATION
            for nd_adj in d_adj[next_nd][i]:
                if selected[nd_adj]:
                    acc_cv += 2 * layer_edge_cts[next_layer - 1 + i]
                    layer_edge_cts[next_layer - 1 + i] -= 1
        if acc_cv > cutoff:
            break
        selected[next_nd] = True
        layer_node_cts[next_layer] += 1
        if layer_node_cts[next_layer] == len(g.layers[next_layer]):
            # convoluted impl but pauses expansion if one side gets ahead by 2 layers or more
            if exp_right - exp_left == 1 and nd_lr == 0:
                open_right = False
            elif exp_left - exp_right == 1 and nd_lr == 1:
                open_left = False
            right_layer += 1 if (nd_lr == 0) else 0
            exp_right += 1 if (nd_lr == 0) else 0
            left_layer -= 1 if (nd_lr == 1) else 0
            exp_left += 1 if (nd_lr == 1) else 0
            if right_layer in g.layers and nd_lr == 0 and open_right:
                for nd in g.layers[right_layer]:
                    heapq.heappush(pqueue, (-(len(d_adj[nd.id][0]) / (len(d_adj[nd.id][1]) + 0.5)), nd.id))
            if left_layer in g.layers and nd_lr == 1 and open_left:
                for nd in g.layers[left_layer]:
                    heapq.heappush(pqueue, (-(len(d_adj[nd.id][1]) / (len(d_adj[nd.id][0]) + 0.5)), nd.id))
        if exp_right - exp_left < 2 and not open_right:
            open_right = True
            if right_layer in g.layers:
                for nd in g.layers[right_layer]:
                    heapq.heappush(pqueue, (-(len(d_adj[nd.id][0]) / (len(d_adj[nd.id][1]) + 0.5)), nd.id))
        elif exp_left - exp_right < 2 and not open_left:
            open_left = True
            if left_layer in g.layers:
                for nd in g.layers[left_layer]:
                    heapq.heappush(pqueue, (-(len(d_adj[nd.id][1]) / (len(d_adj[nd.id][0]) + 0.5)), nd.id))
    return selected


def vertical_re_neighborhood(g: LayeredGraph, candidate_id: int, cutoff, nbhd_width=0):
    """ greedy criterion: largest (incoming edges / (outgoing edges + 0.5))
        this procedure randomly chooses an adjacent layer and adds from that layer until cutoff.
    """
    d_adj = g.get_double_adj_list()
    selected = [False] * g.n_nodes
    if nbhd_width == 0:
        layer_edge_cts = [len(g.get_edges_by_layer()[i]) for i in range(g.n_layers - 1)]
    else:
        layer_edge_cts = find_vertical_nbhd_and_edge_counts(g, candidate_id, nbhd_width)
    for n in g.layers[g[candidate_id].layer]:
        if nbhd_width == 0 or n.v_nbhd:
            selected[n.id] = True
    layer_node_cts = defaultdict(int)
    layer_node_max_cts = [len([v for v in lay if nbhd_width == 0 or v.v_nbhd]) for lay in g.layers.values()]
    layer_node_cts[g[candidate_id].layer] = len(g.layers[g[candidate_id].layer])
    next_layers = [g[candidate_id].layer - 1] if g[candidate_id].layer > 0 else [] + [g[candidate_id].layer + 1] if g[candidate_id].layer < g.n_layers - 1 else []
    next_layer = random.choice(next_layers)
    lr_minmax = sorted([g[candidate_id].layer, next_layer])
    pqueue = []
    acc_cv = 0
    for nd in g.layers[next_layer]:
        if nbhd_width == 0 or nd.v_nbhd:
            heapq.heappush(pqueue, (-(len(d_adj[nd.id][1]) / (len(d_adj[nd.id][0]) + 1)), nd.id))
    while pqueue:
        next_nd = heapq.heappop(pqueue)[1]
        for i in range(2):  # C-VAR CALCULATION
            for nd_adj in d_adj[next_nd][i]:
                if selected[nd_adj]:
                    acc_cv += 2 * layer_edge_cts[g[next_nd].layer - 1 + i]
                    layer_edge_cts[g[next_nd].layer - 1 + i] -= 1
        if acc_cv > cutoff:
            break
        selected[next_nd] = True
        layer_node_cts[next_layer] += 1
        if layer_node_cts[next_layer] == layer_node_max_cts[next_layer]:
            if lr_minmax[0] != 0 or lr_minmax[1] != g.n_layers - 1:
                if lr_minmax[0] == 0:
                    next_layer = lr_minmax[1] + 1
                elif lr_minmax[1] == g.n_layers - 1:
                    next_layer = lr_minmax[0] - 1
                else:
                    next_layer = random.choice([lr_minmax[0] - 1, lr_minmax[1] + 1])
            if next_layer > lr_minmax[1]:
                lr_minmax[1] = next_layer
                for nd in g.layers[next_layer]:
                    if nbhd_width == 0 or nd.v_nbhd:
                        heapq.heappush(pqueue, (-(len(d_adj[nd.id][0]) / (len(d_adj[nd.id][1]) + 1)), nd.id))
            if next_layer < lr_minmax[0]:
                lr_minmax[0] = next_layer
                for nd in g.layers[next_layer]:
                    if nbhd_width == 0 or nd.v_nbhd:
                        heapq.heappush(pqueue, (-(len(d_adj[nd.id][1]) / (len(d_adj[nd.id][0]) + 1)), nd.id))
    return selected


def degree_ratio_neighborhood(g: LayeredGraph, candidate_id: int, cutoff, nbhd_width=0):
    selected = [False] * g.n_nodes
    selected[candidate_id] = True
    if nbhd_width == 0:
        layer_edge_cts = [len(g.get_edges_by_layer()[i]) for i in range(g.n_layers - 1)]
    else:
        layer_edge_cts = find_vertical_nbhd_and_edge_counts(g, candidate_id, nbhd_width)
    print(layer_edge_cts)
    adj = g.get_adj_list()
    degree_ratios = {}
    included_counts = {}
    for n_adj in adj[candidate_id]:
        if nbhd_width == 0 or g[n_adj].v_nbhd:
            degree_ratios[n_adj] = 1 / len(adj[n_adj])
            included_counts[n_adj] = 1
    acc_cv = 0
    n_selected = 0
    while acc_cv <= cutoff and n_selected < g.n_nodes and len(degree_ratios) > 0:
        next_node = max(degree_ratios, key=lambda x: (degree_ratios[x], included_counts[x]))
        del degree_ratios[next_node]
        del included_counts[next_node]
        for nd_adj in adj[next_node]:
            if selected[nd_adj]:
                if g[nd_adj].layer < g[next_node].layer:  # C-VAR CALCULATION
                    acc_cv += 2 * layer_edge_cts[g[nd_adj].layer]
                    layer_edge_cts[g[nd_adj].layer] -= 1
                else:
                    acc_cv += 2 * layer_edge_cts[g[next_node].layer]
                    layer_edge_cts[g[next_node].layer] -= 1
            elif nbhd_width == 0 or g[nd_adj].v_nbhd:
                if nd_adj not in included_counts:
                    included_counts[nd_adj] = 0
                included_counts[nd_adj] += 1
                degree_ratios[nd_adj] = included_counts[nd_adj] / len(adj[nd_adj])
        if acc_cv <= cutoff:
            selected[next_node] = True
            n_selected += 1
        else:
            break
    return selected


def random_neighborhood(g: LayeredGraph, candidate_id: int, cutoff, nbhd_width=0):
    if nbhd_width == 0:
        layer_edge_cts = [len(g.get_edges_by_layer()[i]) for i in range(g.n_layers - 1)]
    else:
        layer_edge_cts = find_vertical_nbhd_and_edge_counts(g, candidate_id, nbhd_width)
    print(layer_edge_cts)
    selected = [False] * g.n_nodes
    selected[candidate_id] = True
    adj = g.get_adj_list()
    acc_cv = 0
    n_selected = 1
    while acc_cv <= cutoff and n_selected < g.n_nodes:
        next_node = random.choice([nd for nd in g.nodes if not selected[nd.id] and (nbhd_width == 0 or nd.v_nbhd)])
        for nd_adj in adj[next_node.id]:
            if selected[nd_adj]:
                if g[nd_adj].layer < next_node.layer:  # C-VAR CALCULATION
                    acc_cv += 2 * layer_edge_cts[g[nd_adj].layer]
                    layer_edge_cts[g[nd_adj].layer] -= 1
                else:
                    acc_cv += 2 * layer_edge_cts[next_node.layer]
                    layer_edge_cts[next_node.layer] -= 1
        if acc_cv <= cutoff:
            selected[next_node.id] = True
            n_selected += 1
        else:
            break
    return selected


def in_vertical_nbhd(yv, bounds, size):
    if bounds[0] <= yv <= bounds[1] or (yv < bounds[0] and bounds[1] - yv <= size) or (yv > bounds[1] and yv - bounds[0] <= size):
        return True
    return False


def find_vertical_nbhd_and_edge_counts(g: LayeredGraph, candidate, percent):
    # ASSUMPTION: all nodes nd have nd.y unique integer in [0, len(nd_layer)]
    edge_counts = [0] * (g.n_layers - 1)
    cnd = g[candidate]
    d_adj = g.get_double_adj_list()
    max_layer_size = max((len(lay) for lay in g.layers.values()))
    size = max_layer_size * percent
    cnd_rel_y = cnd.y / len(g.layers[cnd.layer])
    for lid, lay in g.layers.items():
        scaled_y = cnd_rel_y * len(lay)
        lb = max(0, scaled_y - (size / 2))
        ub = min(len(lay), lb + size)
        lb = min(ub - size, lb)
        for nd in lay:
            if lb <= nd.y <= ub:
                nd.v_nbhd = True
                # if lid > 0:
                #     edge_counts[lid - 1] += len(d_adj[0])
                # if lid < g.n_layers - 1:
                #     edge_counts[lid] += len(d_adj[1])
            else:
                nd.v_nbhd = False
            if lid > 0:
                for nd_adj in d_adj[nd.id][0]:
                    if nd.v_nbhd or g[nd_adj].v_nbhd:
                        edge_counts[lid - 1] += 1
    return edge_counts


def random_candidate(g: LayeredGraph, init=False):
    if init:
        for nd in g:
            nd.energy = 10
    else:
        max_e = max(x.energy for x in g)
        return random.choice([nd.id for nd in g if nd.energy == max_e])


def degree_candidate(g: LayeredGraph, init=False):
    if init:
        adj = g.get_adj_list()
        for nd in g:
            nd.energy = len(adj[nd.id])
    else:
        return max(g.node_ids.keys(), key=lambda x: g[x].energy)


def betweenness_candidate(g: LayeredGraph, init=False):
    if init:
        nxg = g.get_networkx_graph()
        b_c = nx.betweenness_centrality(nxg, normalized=False)
        for nd in g:
            nd.energy = b_c[nd.id]
    else:
        return max(g.node_ids.keys(), key=lambda x: g[x].energy)


def biconnected_candidate(g: LayeredGraph, init=False):
    if init:
        nxg = g.get_networkx_graph()
        try:
            unxg = g.unx_graph
        except AttributeError:
            g.unx_graph = nxg.to_undirected()
            unxg = g.unx_graph
        bicon = nx.biconnected_components(unxg)
        for subset in bicon:
            for nd in subset:
                g[nd].energy += 1
    else:
        return max(g.node_ids.keys(), key=lambda x: g[x].energy)


def avg_edge_length_candidate(g: LayeredGraph, init=False):
    if init:
        for nd in g:
            nd.tabu = False
    else:
        adj = g.get_adj_list()
        # max_seen, max_nd = 0, -1
        for nd in g:
            nd.energy = sum((abs(nd.y - g[a_nd].y) for a_nd in adj[nd.id])) / len(adj[nd.id])
            # if nd.energy > max_seen:  # calculate candidate after penalty but before update?
            #     max_seen = nd.energy
            #     max_nd = nd.id
        return max(g.node_ids.keys(), key=lambda x: g[x].energy if not g[x].tabu else 0)


def crossings_candidate(g: LayeredGraph, init=False):
    if init:
        for nd in g:
            nd.tabu = False
    else:
        e_b_l = g.get_edges_by_layer()
        for edge_list in e_b_l.values():
            for e1, e2 in itertools.combinations(edge_list, 2):
                if len({e1.n1, e1.n2, e2.n1, e2.n2}) == 4:
                    if (e1.n1.y > e2.n1.y and e1.n2.y < e2.n2.y) or (e1.n1.y < e2.n1.y and e1.n2.y > e2.n2.y):
                        g[e1.n1.id].energy += 1
                        g[e1.n2.id].energy += 1
                        g[e2.n1.id].energy += 1
                        g[e2.n2.id].energy += 1
        return max(g.node_ids.keys(), key=lambda x: g[x].energy if not g[x].tabu else 0)


def next_candidate(g: LayeredGraph):
    return max(g.node_ids.keys(), key=lambda x: g[x].energy)


def penalty_fn(g: LayeredGraph, neighborhood, candidate, movement, iteration, no_repeats=False):
    var_1, var_2 = 8, 4
    var_moved = 2  # factor for penalty on a neighborhood node that DID move (higher = less penalty)
    end_penalty_iter = 20  # subsequent iterations only penalize nodes that don't move
    # var_1 = (var_1 - 1) * (math.e ** (-1/end_penalty_iter * iteration)) + 1
    var_1 = (var_1 - 1) * (1 / (iteration // 10 + 1)) + 1
    # var_2 = (var_2 - 1) * (math.e ** (-1/end_penalty_iter * iteration)) + 1
    var_2 = (var_2 - 1) * (1 / (iteration // 10 + 1)) + 1
    adj = g.get_adj_list()
    g[candidate].tabu = True
    if no_repeats:
        g[candidate].energy = 0
    for i, nd in enumerate(neighborhood):
        if movement[i] == 0:
            if nd == candidate:
                g[nd].energy /= var_1
            else:
                g[nd].energy /= var_2
        else:
            if nd == candidate:
                g[nd].energy /= max(var_1 / var_moved, 1)
            else:
                g[nd].energy /= max(var_2 / var_moved, 1)
            g[nd].tabu = False
            for nd_adj in adj[nd]:
                g[nd_adj].tabu = False
