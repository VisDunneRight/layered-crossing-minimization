import random
import itertools
from src import graph


def normal_c_vars(g: graph.LayeredGraph, edges_by_layer, mirror_vars):
    c_vars = []
    constants = []
    for i, edge_list in edges_by_layer.items():
        if mirror_vars:
            for pr in itertools.permutations(edge_list, 2):
                if pr[0][0] != pr[1][0] and pr[0][1] != pr[1][1] and pr[0][0] != pr[1][1] and pr[0][1] != pr[1][0]:
                    c_vars.append(pr)
                    constants.append(g.edge_ids[pr[0]].weight * g.edge_ids[pr[1]].weight)
        else:
            for pr in itertools.combinations(edge_list, 2):
                if pr[0][0] != pr[1][0] and pr[0][1] != pr[1][1] and pr[0][0] != pr[1][1] and pr[0][1] != pr[1][0]:
                    c_vars.append(pr)
                    constants.append(g.edge_ids[pr[0]].weight * g.edge_ids[pr[1]].weight)
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
    graph_cuts = []
    corresponding_cross = []
    for i in range(n_iter):
        edges = [[e.n1.id, e.n2.id] for e in g.edges]
        random.shuffle(edges)
        combined = {n.id: {n.id} for n in g.nodes}
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
        if len(cross_edges) <= 3:
            cut = [list(combined[cross_edges[0][0]]), list(combined[cross_edges[0][1]])]
            cut.sort(key=lambda x: len(x))
            if is_large_subgraph(g, cut[0]):
                if all((set(cut[0]).isdisjoint(other_cut) for other_cut in graph_cuts)):
                    graph_cuts.append(cut[0])
                    corresponding_cross.append(cross_edges)
                else:
                    ind = [i for i in range(len(graph_cuts)) if not set(cut[0]).isdisjoint(graph_cuts[i])]
                    if sum((len(graph_cuts[i]) for i in ind)) < len(cut[0]):
                        graph_cuts = [val for j, val in enumerate(graph_cuts) if j not in ind]
                        graph_cuts.append(cut[0])
                        corresponding_cross = [val for j, val in enumerate(corresponding_cross) if j not in ind]
                        corresponding_cross.append(cross_edges)
        # print(combined[cross_edges[0][0]], combined[cross_edges[0][1]])
        # print("cut size:", len(cross_edges))
    print("Graph cutsets:", graph_cuts, "sizes:", [len(edl) for edl in corresponding_cross])
    return graph_cuts, corresponding_cross


def is_large_subgraph(g: graph.LayeredGraph, subg_nodes):
    if len(subg_nodes) >= 5:
        if len([g.node_ids[n].layer for n in subg_nodes]) != len(set((g.node_ids[n].layer for n in subg_nodes))):
            return True
    return False


def expand_cut(subg_nodes, contacts, normal_adjacency, allotted_n):
    """
    Idea: run DFS on graph w/o connections to subgraph. Store the size s of the subtree under every node.
    """
    n_left = allotted_n - len(subg_nodes)
    for node in contacts:
        exit_cond = False
        prev_node = node
        for node_adj in normal_adjacency[node]:
            if node_adj not in subg_nodes:
                cur_node = node_adj
        while n_left > 0 and not exit_cond:
            if len(normal_adjacency[cur_node]) == 2:
                subg_nodes.append()


def branching_path(subg_nodes: set, adjacency, contact_node, branch_nodes, alotted_n):
    to_scoop = set(branch_nodes)
    cur_nodes = []
    # for b_node in branch_nodes:
    #     for adj in adjacency[b_node]:
    #         if adj not in subg_nodes and adj not in to_scoop:


