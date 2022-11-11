import random
import itertools
from src import graph


def normal_c_vars(g: graph.LayeredGraph, edges_by_layer):
    c_vars = []
    constants = []
    for i, edge_list in edges_by_layer.items():
        for pr in itertools.combinations(edge_list, 2):
            if pr[0][0] != pr[1][0] and pr[0][1] != pr[1][1] and pr[0][0] != pr[1][1] and pr[0][1] != pr[1][0]:
                c_vars.append(pr)
                constants.append(g.edge_names[pr[0]].weight * g.edge_names[pr[1]].weight)
    return c_vars, constants


def c_vars_with_crossing_var_sum_reduction(g: graph.LayeredGraph, nodes_by_layer, edges_by_layer, multiedge_threshold):
    adjacency = g.create_double_adj_list(forward_only=True)
    c_vars = []
    special_c_vars = []
    for i, edge_list in edges_by_layer.items():
        special_nodes = []
        for node in nodes_by_layer[i]:
            if len(adjacency[node]) >= multiedge_threshold:
                special_nodes.append(tuple((node, adj) for adj in adjacency[node]))
        for edge_gp in special_nodes:
            for edge in edge_list:
                if edge not in edge_gp:
                    special_c_vars.append((edge, tuple(edge_gp)))
        plain_edges = [edge for edge in edge_list if edge not in {item for subl in special_nodes for item in subl}]
        for pr in itertools.combinations(plain_edges, 2):
            if pr[0][0] != pr[1][0] and pr[0][1] != pr[1][1] and pr[0][0] != pr[1][1] and pr[0][1] != pr[1][0]:
                c_vars.append(pr)
    return c_vars, special_c_vars


def kargers_algo_cut_finder(g: graph.LayeredGraph, n_iter):
    min_cut_max_size = [[], []]
    corresponding_cross = []
    for i in range(n_iter):
        edges = [[e.n1.name, e.n2.name] for e in g.edges]
        random.shuffle(edges)
        combined = {n.name: {n.name} for n in g.nodes}
        num_nodes = g.n_nodes
        while num_nodes > 2:
            contract = edges.pop()
            if combined[contract[0]] is not combined[contract[1]]:
                new_node = combined[contract[0]].union(combined[contract[1]])
                for v in new_node:
                    combined[v] = new_node
                num_nodes -= 1
        cross_edges = []
        for edge in edges:
            if combined[edge[0]] is not combined[edge[1]]:
                cross_edges.append(edge)
        if len(cross_edges) <= 2:
            cut = [list(combined[cross_edges[0][0]]), list(combined[cross_edges[0][1]])]
            cut.sort(key=lambda x: len(x))
            if len(cut[0]) > len(min_cut_max_size[0]):
                min_cut_max_size = cut
                corresponding_cross = cross_edges
        # print(combined[cross_edges[0][0]], combined[cross_edges[0][1]])
        # print("cut size:", len(cross_edges))
    print("cut size: ", len(corresponding_cross), min_cut_max_size)
    return min_cut_max_size, corresponding_cross


def stack_edges_from_cutset(g: graph.LayeredGraph, cutset):
    by_layer = {}
    cutnodes = set(cutset[0])
    for edge in g.edges:
        if edge.n1 in cutnodes and edge.n2 in cutnodes:
            if edge.n1.layer not in by_layer:
                by_layer[edge.n1.layer] = []
            by_layer[edge.n1.layer].append(edge)
    for cutlist in by_layer.values():
        g.stack_subgraph(cutlist)
