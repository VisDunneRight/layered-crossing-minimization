from src import graph


def crossing_var_sum_reduction_finder(g: graph.LayeredGraph, nodes_by_layer, c_vars_layers):
    double_adj = g.create_double_adj_list()
    constraints_to_remove = []
    key_edges = []                  # key edges should instead be an item which can easily be used to form a constraint
    for i in range(1, len(nodes_by_layer)+1):
        for node in nodes_by_layer[i]:
            if len(double_adj[node][1]) >= 3:
                key_edges.append([])
                for adj in double_adj[node][1]:
                    key_edges[-1].append((node, adj))
                    for c_var in c_vars_layers[i]:
                        if c_var[0] == (node, adj) or c_var[1] == (node, adj):
                            constraints_to_remove.append(c_var)
    return key_edges, constraints_to_remove
