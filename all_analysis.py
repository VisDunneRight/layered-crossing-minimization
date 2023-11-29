from src import optimization, graph, heuristics, vis
import json


if __name__ == '__main__':
    # with open("data storage/All/AK_Distribution_LastMileNetwork_All_2022.json", 'r') as fd:
    #     data = json.load(fd)
    with open("data storage/All/Distribution_reformatted.json", 'r') as fd:
        data = json.load(fd)
        # json.dump(data, fd, indent=4)
    g = graph.LayeredGraph()
    seen_subs = set()
    merge_edges = False  # Merge duplicate edges, totaling the supply delivered
    for n1 in data:
        x1 = g.add_node(0, name=data[n1]["DistributionCenterStatistics"]["NAME"] + "_" + str(data[n1]["DistributionCenterStatistics"]["OBJECTID"]))
        # x1 = g.add_node(0, name=str(data[n1]["DistributionCenterStatistics"]["OBJECTID"]))
        for adj in data[n1]["Destinations"]:
            if data[n1]["Destinations"][adj]["Substances"] not in seen_subs:
                seen_subs.add(data[n1]["Destinations"][adj]["Substances"])
            if data[n1]["Destinations"][adj]["Name"] not in g:
                x2 = g.add_node(1, name=data[n1]["Destinations"][adj]["Name"])
            else:
                x2 = g[data[n1]["Destinations"][adj]["Name"]]
            if merge_edges and (x1.id, x2.id) in g.edge_ids:
                g.edge_ids[(x1.id, x2.id)].weight += float(data[n1]["Destinations"][adj]["Delivered"])
            elif merge_edges:
                g.add_edge(x1.id, x2.id, weight=float(data[n1]["Destinations"][adj]["Delivered"]))
            else:
                g.add_edge(x1.id, x2.id, weight=float(data[n1]["Destinations"][adj]["Delivered"]), data={"substance": data[n1]["Destinations"][adj]["Substances"]})
    print(seen_subs)

    # Restrict to only one type
    gasoline_g = graph.LayeredGraph()
    for nd in g.nodes:
        gasoline_g.add_node(nd.layer, idx=nd.id, name=nd.name)
    for ed in g.edges:
        if g.edge_data["substance"][(ed.n1.id, ed.n2.id)] == "Gasoline":
            gasoline_g.add_edge(ed.n1.id, ed.n2.id, weight=ed.weight)
    g = gasoline_g

    # Collapse nodes with same adjacency
    subgs = [0] * g.n_nodes
    adj = g.get_adj_list()
    seen_gps = {}
    n_gps = 1
    for nd in g.layers[1]:
        nd_gps = tuple(sorted(adj[nd.id]))
        if nd_gps not in seen_gps:
            seen_gps[nd_gps] = n_gps
            subgs[nd.id] = n_gps
            n_gps += 1
        else:
            subgs[nd.id] = seen_gps[nd_gps]
    print(seen_gps)
    for nd in g.layers[0]:
        if not adj[nd.id]:
            subgs[nd.id] = n_gps
    for gp in range(n_gps):
        if sum((1 for s in subgs if s == gp)) == 1:
            subgs = [0 if s == gp else s for s in subgs]
    subgs_ord = sorted(list(set(subgs)))
    for i, x in enumerate(subgs.copy()):
        sum_missing = sum((1 for v in subgs_ord if v < x))
        subgs[i] = sum_missing
    g = g.stacked_graph_from_subgraph_nodes(subgs, only_subgraphs=True)

    heuristics.split(g, n_iter=30)
    heuristics.improved_sifting(g, n_iter=100)
    opt = optimization.LayeredOptimizer(g)
    opt.name = "DISTRIBUTION_gas"
    # opt.symmetry_breaking = True
    # opt.collapse_leaves = True
    # opt.vertical_transitivity = True
    # opt.direct_transitivity = False
    # opt.heuristic_start = True
    # opt.cycle_constraints = True
    # opt.butterfly_reduction = True
    # opt.draw_graph = True
    # opt.bendiness_reduction = True
    # opt.optimize_layout()
    opt.just_bendiness_reduction()
    vis.draw_graph(g, opt.name, node_x_distance=700)
