from src import graph


def count_butterflies(g: graph.LayeredGraph):
    adjacency = g.get_double_adj_list()
    butterfly_count = 0
    for i in range(1, g.n_layers):
        wedges = []
        for node_1 in g.layers[i]:
            for adj in adjacency[node_1.id][1]:
                for node_2 in adjacency[adj][0]:
                    if node_1.id < node_2:
                        wedges.append((node_1.id, adj, node_2))
        for j, wedge_1 in enumerate(wedges):
            for wedge_2 in wedges[j+1:]:
                if wedge_1[0] == wedge_2[0] and wedge_1[2] == wedge_2[2] and wedge_1[1] != wedge_2[1]:
                    butterfly_count += 1
    return butterfly_count


def get_butterflies(g: graph.LayeredGraph):
    adjacency = g.get_double_adj_list()
    butterflies = []
    for i in range(1, g.n_layers):
        wedges = []
        for node_1 in g.layers[i]:
            for adj in adjacency[node_1.id][1]:
                for node_2 in adjacency[adj][0]:
                    if node_1.id < node_2:
                        wedges.append((node_1.id, adj, node_2))
        for j, wedge_1 in enumerate(wedges):
            for wedge_2 in wedges[j + 1:]:
                if wedge_1[0] == wedge_2[0] and wedge_1[2] == wedge_2[2] and wedge_1[1] != wedge_2[1]:
                    butterflies.append((wedge_1[0], wedge_1[1], wedge_1[2], wedge_2[1]))
    return butterflies


def get_3claws(g: graph.LayeredGraph):
    bearclaws = []
    dbl_adj_list = g.get_double_adj_list()
    for i in range(1, g.n_layers):
        for node_1 in g.layers[i]:
            wedges = []
            for adj in dbl_adj_list[node_1.id][1]:
                for node_2 in dbl_adj_list[adj][0]:
                    if node_2 != node_1.id:
                        wedges.append((node_2, adj))
            for j, wedge_1 in enumerate(wedges):
                for k, wedge_2 in enumerate(wedges[j + 1:]):
                    for wedge_3 in wedges[j + k + 2:]:
                        if wedge_1[0] != wedge_2[0] and wedge_1[0] != wedge_3[0] and wedge_2[0] != wedge_3[0] and wedge_1[1] != wedge_2[1] and wedge_1[1] != wedge_3[1] and wedge_2[1] != wedge_3[1]:
                            bearclaws.append(((node_1.id, wedge_1[1]), (node_1.id, wedge_2[1]), (node_1.id, wedge_3[1]), wedge_1, wedge_2, wedge_3))
            wedges.clear()
            for adj in dbl_adj_list[node_1.id][0]:
                for node_2 in dbl_adj_list[adj][1]:
                    if node_2 != node_1.id:
                        wedges.append((adj, node_2))
            for j, wedge_1 in enumerate(wedges):
                for k, wedge_2 in enumerate(wedges[j + 1:]):
                    for wedge_3 in wedges[j + k + 2:]:
                        if wedge_1[0] != wedge_2[0] and wedge_1[0] != wedge_3[0] and wedge_2[0] != wedge_3[0] and wedge_1[1] != wedge_2[1] and wedge_1[1] != wedge_3[1] and wedge_2[1] != wedge_3[1]:
                            bearclaws.append(((wedge_1[0], node_1.id), (wedge_2[0], node_1.id), (wedge_3[0], node_1.id), wedge_1, wedge_2, wedge_3))
    return bearclaws


def get_domepaths(g: graph.LayeredGraph):
    domes = []
    dbl_adj_list = g.get_double_adj_list()
    for i in range(1, g.n_layers):
        for node_1 in g.layers[i]:
            wedges = []
            for adj in dbl_adj_list[node_1.id][1]:
                for node_2 in dbl_adj_list[adj][0]:
                    if node_2 != node_1.id:
                        wedges.append((node_2, adj))
            for j, wedge_1 in enumerate(wedges):
                for k, wedge_2 in enumerate(wedges[j + 1:]):
                    if wedge_1[0] != wedge_2[0] and wedge_1[1] != wedge_2[1]:
                        domes.append(((node_1.id, wedge_1[1]), (node_1.id, wedge_2[1]), wedge_1, wedge_2))
            wedges.clear()
            for adj in dbl_adj_list[node_1.id][0]:
                for node_2 in dbl_adj_list[adj][1]:
                    if node_2 != node_1.id:
                        wedges.append((adj, node_2))
            for j, wedge_1 in enumerate(wedges):
                for k, wedge_2 in enumerate(wedges[j + 1:]):
                    if wedge_1[0] != wedge_2[0] and wedge_1[1] != wedge_2[1]:
                        domes.append(((wedge_1[0], node_1.id), (wedge_2[0], node_1.id), wedge_1, wedge_2))
    return domes


# def get_2layer_cycles(g: graph.LayeredGraph):
#     cycles = []
#     domes = get_domepaths(g)
#     butterfiles = get_butterflies(g)
#     for

