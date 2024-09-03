import time
from src.graph import LayeredGraph
from copy import deepcopy
import numpy as np


def tabu(g: LayeredGraph, nonimp=50, maxiter=-1, timelimit=-1, improvement_data=False):
    """
    :param g: required, LayeredGraph object to perform tabu search on
    :param nonimp: Number of non-improving iterations before program terminates, default 50
    :param maxiter: Number of second-level diversification swaps each iteration, if left as -1 will use 25*|V|
    :param timelimit: If set, will instead run until the specified time in seconds has elapsed
    :param improvement_data: If set, will record and return the best solution values and the time of improvement
    :return: None. positions of nodes of g will be updated according to the method of Laguna et al. in the 1997 article, 'Arc Crossing Minimization in Heirarchical Digraphs with Tabu Search'
    """
    assert g.is_proper(), "graph is not properly layered, try calling g.add_anchors() and then g.relayer()"
    g.y_val_setup()  # 1-indexes y-values and sorts layer
    cur_solution_pi = np.array([g[i].y - 1 for i in range(g.n_nodes)], np.int_)
    cur_solution_phi = {layer_id: np.array([nd.id for nd in g.layers[layer_id]], np.int_) for layer_id in g.layers}
    best_solution_pi = deepcopy(cur_solution_pi)
    d_adj = g.get_double_adj_list()
    attractiveness = np.array([sum(len(d_adj[v.id][0]) + len(d_adj[v.id][1]) for v in g.layers[l_i]) for l_i in range(g.n_layers)], np.int_)
    tabu_list = np.zeros(g.n_layers, np.int_)
    if maxiter == -1:
        maxiter = 25 * g.n_nodes
    lastimp = 0
    cur_iter = 0
    rng = np.random.default_rng()
    best_cr_val = g.num_edge_crossings()
    start_time = time.time()
    if improvement_data:
        time_values = [0]
        crossing_values = [best_cr_val]

    # begin tabu search loop
    while lastimp < nonimp or time.time() - start_time < timelimit:
        cur_iter += 1
        reduced_this_iter = False
        print("Cur iter:", cur_iter)

        # 1. perform initial left-right intensification sweep
        for layer_id in range(g.n_layers):
            crossings_improved = intensify_layer(g, cur_solution_pi, cur_solution_phi, layer_id, d_adj)
            if crossings_improved:
                tabu_list[layer_id] = cur_iter
                reduced_this_iter = True
            else:
                tabu_list[layer_id] = -cur_iter

        # 2. inner while loop, use first-level diversification and intensify accordingly
        next_layer = tabu_diversify_l1(tabu_list, attractiveness)
        inner_iter = 1
        while next_layer != -1:
            crossings_improved = intensify_layer(g, cur_solution_pi, cur_solution_phi, next_layer, d_adj)
            if crossings_improved:
                tabu_list[next_layer] = inner_iter
            else:
                tabu_list[next_layer] = -inner_iter
            inner_iter += 1
            next_layer = tabu_diversify_l1(tabu_list, attractiveness)

        # 3. check if solution is improved
        for nd, rk in enumerate(cur_solution_pi):
            g[nd].y = rk
        cur_crossings = g.num_edge_crossings()
        print(f"Improved best of {best_cr_val}?", f"yes, to {cur_crossings}" if cur_crossings < best_cr_val else "no")

        if cur_crossings < best_cr_val:
            best_solution_pi = deepcopy(cur_solution_pi)
            lastimp = 0
            best_cr_val = cur_crossings
            if improvement_data:
                time_values.append(time.time() - start_time)
                crossing_values.append(best_cr_val)
        else:
            lastimp += 1

        # 4. second-level diversification
        n_swaps_made = 0
        while n_swaps_made < maxiter:
            n_swaps_made += tabu_diversify_l2(g, cur_solution_pi, cur_solution_phi, rng)

    for nd, rk in enumerate(best_solution_pi):
        g[nd].y = rk
    print(f"\n... cut off solve after {nonimp} non-improving iterations" if timelimit == -1 else f"... cut off at {timelimit} seconds")

    if improvement_data:
        return crossing_values, time_values


def intensify_layer(g: LayeredGraph, cur_ranking, cur_phi, layer_id, d_adj):
    # 1. compute K(u,v) for all u,v in L_i
    crossing_matrix, crossing_map = get_crossing_matrix_layer_i(d_adj, cur_ranking, g.layers[layer_id])
    # 2. for each node, get positions that minimize the number of crossings, and swap to best. Repeat until no swaps remain
    crossings_reduced = False
    reduced_this_iter = True
    while reduced_this_iter:
        reduced_this_iter = False
        for nd in g.layers[layer_id]:
            max_indices, improves = efficient_movement_calculation(crossing_matrix, crossing_map, cur_ranking, cur_phi, nd.id, layer_id)
            # print(improves)
            if max(improves) == 0:  # don't perform swap
                continue
            elif len(max_indices) > 1:  # use barycenter to break ties
                nd_barycenter = sum(cur_ranking[v] for v in d_adj[nd.id][0] + d_adj[nd.id][1]) / (len(d_adj[nd.id][0]) + len(d_adj[nd.id][1]))
                best_idx = min(max_indices, key=lambda x: abs(x - nd_barycenter))
                if cur_ranking[nd.id] not in max_indices:
                    crossings_reduced, reduced_this_iter = True, True
                    # print(nd.id, max_indices, improves)
            else:
                best_idx = max_indices[0]
                crossings_reduced, reduced_this_iter = True, True
                # print(nd.id, max_indices, improves)
            if cur_ranking[nd.id] != best_idx:
                insert_at_position(nd.id, best_idx, cur_ranking, cur_phi[layer_id])
    # 3. Perform final pass, inserting each node at its barycenter if it does not increase crossings
    for nd in g.layers[layer_id]:
        _, e_list = efficient_movement_calculation(crossing_matrix, crossing_map, cur_ranking, cur_phi, nd.id, layer_id)
        nd_barycenter = sum(cur_ranking[v] for v in d_adj[nd.id][0] + d_adj[nd.id][1]) / (len(d_adj[nd.id][0]) + len(d_adj[nd.id][1]))
        bary_val = round(nd_barycenter + 0.0001)  # avoid rounding weirdness for numbers ending in .5
        if bary_val >= len(e_list):
            bary_val = len(e_list) - 1
        bary_cr = e_list[bary_val]
        if bary_cr >= 0:
            insert_at_position(nd.id, bary_val, cur_ranking, cur_phi[layer_id])
            if bary_cr > 0:
                crossings_reduced = True
    g.layers[layer_id].sort(key=lambda x: cur_ranking[x.id])
    return crossings_reduced


def get_crossing_matrix_layer_i(d_adj, rank, layer_i_nodes):
    """
    :param d_adj: double adjacency list (split by forward/backward edges)
    :param rank: ordering of all nodes (node id -> position in layer)
    :param layer_i_nodes: list of node objects in layer i
    :return: matrix of K(u,v) values for layer i
    """
    crossing_matrix = np.zeros((len(layer_i_nodes), len(layer_i_nodes)), np.int_)
    cr_map = {nd.id: i for i, nd in enumerate(layer_i_nodes)}
    for j, u in enumerate(layer_i_nodes):
        for v in layer_i_nodes[j + 1:]:
            for p in range(2):
                for u_adj in d_adj[u.id][p]:
                    for v_adj in d_adj[v.id][p]:
                        if u_adj != v_adj:
                            if rank[u_adj] < rank[v_adj]:
                                crossing_matrix[cr_map[v.id]][cr_map[u.id]] += 1
                                # crossing_matrix[u.id][rank[v.id]] -= 1
                            else:
                                crossing_matrix[cr_map[u.id]][cr_map[v.id]] += 1
                                # crossing_matrix[v.id][rank[u.id]] -= 1
    return crossing_matrix, cr_map


def efficient_movement_calculation(cr_matrix, cr_map, rank, phi, nd_id, layer_id):
    """
    :param cr_matrix: K(u,v) matrix for layer i
    :param cr_map:
    :param rank:
    :param phi:
    :param nd_id:
    :param layer_id:
    :return: list of indices for a node phi(j) such that E(phi(j),j') is maximal, using the efficient calculation of Laguna et al.
    """
    index_list = [rank[nd_id]]
    e_list = np.zeros(len(cr_matrix), np.int_)
    max_move_val = 0
    cur_val = 0
    for i in range(rank[nd_id] + 1, len(cr_matrix)):
        cur_val += cr_matrix[cr_map[nd_id]][cr_map[phi[layer_id][i]]] - cr_matrix[cr_map[phi[layer_id][i]]][cr_map[nd_id]]
        e_list[i] = cur_val
        if cur_val > max_move_val:
            index_list = [i]
            max_move_val = cur_val
        elif cur_val == max_move_val:
            index_list.append(i)
    cur_val = 0
    for i in range(rank[nd_id] - 1, -1, -1):
        cur_val += cr_matrix[cr_map[phi[layer_id][i]]][cr_map[nd_id]] - cr_matrix[cr_map[nd_id]][cr_map[phi[layer_id][i]]]
        e_list[i] = cur_val
        if cur_val > max_move_val:
            index_list = [i]
            max_move_val = cur_val
        elif cur_val == max_move_val:
            index_list.append(i)
    return index_list, e_list


def insert_at_position(node_id, target_pos, ranking, phi_layer) -> None:
    if ranking[node_id] < target_pos:
        for idx in range(ranking[node_id] + 1, target_pos + 1):
            ranking[phi_layer[idx]] -= 1
            phi_layer[idx - 1] = phi_layer[idx]
    elif ranking[node_id] > target_pos:
        for idx in range(ranking[node_id] - 1, target_pos - 1, -1):
            ranking[phi_layer[idx]] += 1
            phi_layer[idx + 1] = phi_layer[idx]
    ranking[node_id] = target_pos
    phi_layer[target_pos] = node_id


def tabu_diversify_l1(tabu_list, attractiveness):
    unnormalized_probs = np.array([attractiveness[i] if
                                   ((i < len(tabu_list) - 1 and tabu_list[i + 1] >= abs(tabu_list[i])) or
                                    (i > 0 and tabu_list[i - 1] >= abs(tabu_list[i])))
                                   else 0
                                   for i in range(len(tabu_list))
                                   ])
    if np.any(unnormalized_probs):
        return np.random.choice(len(tabu_list), p=unnormalized_probs / np.sum(unnormalized_probs))
    else:
        # no non-tabu layer exists
        return -1


def tabu_diversify_l2(g: LayeredGraph, cur_pi, cur_phi, rng):
    nd_to_swap = int(rng.integers(low=0, high=g.n_nodes))
    swap_dir = rng.choice([-1, 1])
    lid = g[nd_to_swap].layer
    if (cur_pi[nd_to_swap] == 0 and swap_dir == -1) or (cur_pi[nd_to_swap] == len(cur_phi[lid]) - 1 and swap_dir == 1):
        return 0
    else:
        cur_pos = cur_pi[nd_to_swap]
        swap_target = cur_phi[lid][cur_pos + swap_dir]
        cur_phi[lid][cur_pos], cur_phi[lid][cur_pos + swap_dir] = cur_phi[lid][cur_pos + swap_dir], cur_phi[lid][cur_pos]
        cur_pi[nd_to_swap], cur_pi[swap_target] = cur_pi[swap_target], cur_pi[nd_to_swap]
        return 1
