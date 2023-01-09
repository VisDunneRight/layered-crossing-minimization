from src import graph


def count_butterflies(g: graph.LayeredGraph):
    adjacency = g.create_double_adj_list()
    butterfly_count = 0
    for i in range(1, g.n_layers):
        wedges = []
        for node_1 in g.layers[i]:
            for adj in adjacency[node_1.name][1]:
                for node_2 in adjacency[adj][0]:
                    if node_1.name != node_2:
                        wedges.append((node_1.name, adj, node_2))
        for j, wedge_1 in enumerate(wedges):
            for wedge_2 in wedges[j+1:]:
                if wedge_1[0] == wedge_2[0] and wedge_1[2] == wedge_2[2] and wedge_1[1] != wedge_2[1]:
                    butterfly_count += 1
    return butterfly_count


def get_butterflies(g: graph.LayeredGraph):
    adjacency = g.create_double_adj_list()
    butterflies = []
    for i in range(1, g.n_layers):
        wedges = []
        for node_1 in g.layers[i]:
            for adj in adjacency[node_1.name][1]:
                for node_2 in adjacency[adj][0]:
                    if node_1.name != node_2:
                        wedges.append((node_1.name, adj, node_2))
        for j, wedge_1 in enumerate(wedges):
            for wedge_2 in wedges[j + 1:]:
                if wedge_1[0] == wedge_2[0] and wedge_1[2] == wedge_2[2] and wedge_1[1] != wedge_2[1]:
                    butterflies.append((wedge_1, wedge_2))
    return butterflies
