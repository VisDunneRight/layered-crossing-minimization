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


def get_groups(g: graph.LayeredGraph):
    if "groups" not in g.node_data:
        print("Need to add group assignments on nodes, use g.add_groups. Proceeding without groups.")
        return [], []
    groups = g.node_data["groups"]
    has_nested_groups = False if all(type(v) == int or type(v) == float for v in groups.values()) else True
    gp_vals = set()
    for v in groups.values():
        if type(v) == list:
            gp_vals.update(v)
        else:
            gp_vals.add(v)
    all_groups = [[nd for nd, val in groups.items() if val == gp_id] for gp_id in gp_vals]
    is_sl_group = [True if all(g[gp[i]].layer == g[gp[i+1]].layer for i in range(len(gp) - 1)) else False for gp in all_groups]
    sl_groups = [gp for i, gp in enumerate(all_groups) if is_sl_group[i]]
    ml_groups = [gp for i, gp in enumerate(all_groups) if not is_sl_group[i]]
    if not has_nested_groups:
        ml_idx = 0
        for i, gp in enumerate(all_groups):
            gp_id = groups[gp[0]]
            gp_layers = sorted(set(g[nd].layer for nd in gp))
            if not all(gp_layers[ix] + 1 == gp_layers[ix + 1] for ix in range(len(gp_layers) - 1)):
                print("Some groups are not continuous across layers. Proceeding without groups.")
                return [], []
            gp_by_layer = [[nd for nd in gp if g[nd].layer == lid] for lid in gp_layers]
            if not is_sl_group[i]:
                max_gp_layer_sz = max(len(lay) for lay in gp_by_layer)
                for layer_subgp in gp_by_layer:
                    lay_id = g[layer_subgp[0]].layer
                    for _ in range(max_gp_layer_sz - len(layer_subgp)):
                        new_nd = g.add_node(lay_id, data={"groups": gp_id, "filler": True})
                        ml_groups[ml_idx].append(new_nd.id)
                ml_idx += 1

        # filler nodes for dir. trans. ver.
        layers_with_ml_gps = set(g[nd].layer for gp in ml_groups for nd in gp)
        if layers_with_ml_gps:
            max_graph_layer_sz = max(len(g.layers[lay]) for lay in layers_with_ml_gps)
            for lid in layers_with_ml_gps:
                for _ in range(max_graph_layer_sz - len(g.layers[lid])):
                    g.add_node(lid, data={"filler": True})
    else:
        print("Multiple group membership not yet supported. Proceeding without groups")
    return sl_groups, ml_groups


def remove_filler_nodes(g: graph.LayeredGraph, x_vars):
    if "filler" not in g.node_data:
        return
    n_removed = 0
    n_b_l = g.get_ids_by_layer()
    for i in range(g.n_nodes - 1, -1, -1):
        nid = g.nodes[i].id
        if nid in g.node_data["filler"] and nid not in g.node_data["groups"]:
            nd_pop = g.nodes.pop(i)
            g.layers[nd_pop.layer].remove(nd_pop)
            g.node_ids.pop(nd_pop.id)
            n_removed += 1
            for nd in n_b_l[nd_pop.layer]:
                if nd != nid:
                    if (nd, nid) in x_vars:
                        del x_vars[nd, nid]
                    elif (nid, nd) in x_vars:
                        del x_vars[nid, nd]
    g.n_nodes -= n_removed
    g.invalidate_data()


def remove_filler_group_nodes(g: graph.LayeredGraph, x_vars):
    if "filler" not in g.node_data:
        return
    n_removed = 0
    n_b_l = g.get_ids_by_layer()
    for i in range(g.n_nodes - 1, -1, -1):
        nid = g.nodes[i].id
        if nid in g.node_data["filler"]:
            g.node_data["groups"].pop(nid)
            nd_pop = g.nodes.pop(i)
            g.layers[nd_pop.layer].remove(nd_pop)
            g.node_ids.pop(nd_pop.id)
            n_removed += 1
            for nd in n_b_l[nd_pop.layer]:
                if nd != nid:
                    if (nd, nid) in x_vars:
                        del x_vars[nd, nid]
                    elif (nid, nd) in x_vars:
                        del x_vars[nid, nd]
    g.node_data.pop("filler")
    g.n_nodes -= n_removed
    g.invalidate_data()
