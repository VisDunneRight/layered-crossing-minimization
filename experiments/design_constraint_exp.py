import collections
import csv
import random
from src.benchmark import generate_benchmark
from src.optimization import LayeredOptimizer
from src.heuristics import improved_sifting


def add_groups(g):
    groups = {}
    gp_id, gp_nd_sum, gp_target = 0, 0, max(round(0.04 * g.n_nodes), 2)
    while gp_nd_sum < 0.2 * g.n_nodes:
        if random.random() < 0.5:  # sl-group
            success = False
            r_layer = random.choice([lid for lid in g.layers if len(g.layers[lid]) >= 2])
            for _ in range(random.randint(gp_target, max(len(g.layers[r_layer]) // 2, gp_target))):
                rns = [v.id for v in g.layers[r_layer] if v.id not in groups]
                if rns:
                    r_n = random.choice(rns)
                    groups[r_n] = gp_id
                    gp_nd_sum += 1
                    success = True
            if success:
                gp_id += 1
        else:  # BFS ml-group
            r_n = random.choice([v.id for v in g.nodes if v.id not in groups])
            a_l = g.get_adj_list()
            bfsq = [r_n]
            seen = [False] * g.n_nodes
            seen[r_n] = True
            gp_nds = []
            for _ in range(random.randint(gp_target, 3 * gp_target)):
                if bfsq:
                    next_nd = bfsq.pop(0)
                    if next_nd not in groups:
                        gp_nds.append(next_nd)
                    for v_adj in a_l[next_nd]:
                        if not seen[v_adj]:
                            seen[v_adj] = True
                            bfsq.append(v_adj)
            gp_layers = set((g[nd].layer for nd in gp_nds))
            # print(gp_nds, gp_layers)
            if len(gp_nds) > len(gp_layers) > max(gp_layers) - min(gp_layers):
                for gpnd in gp_nds:
                    groups[gpnd] = gp_id
                    gp_nd_sum += 1
                gp_id += 1
    # print(groups)
    g.add_groups(groups)
    # vis.draw_graph(g, "testing")


def add_sl_edges(g):
    sltoremove = round(0.1 * len(g.edges))
    nremoved = 0
    while nremoved < sltoremove:
        r_e = random.choice(g.edges)
        while r_e.n1.is_anchor_node or r_e.n2.is_anchor_node:
            r_e = random.choice(g.edges)
        ot_nodes = [nd.id for nd in g.layers[r_e.n1.layer] if nd.id != r_e.n1.id]
        if len(ot_nodes) > 0:
            nd_ot = random.choice(ot_nodes)
            if (r_e.n1.id, nd_ot) not in g.edge_ids:
                nremoved += 1
                g.edges.remove(r_e)
                del g.edge_ids[r_e.n1.id, r_e.n2.id]
                g.add_edge(r_e.n1.id, nd_ot)


def run_func(combo_idx, data_path):
    opt = LayeredOptimizer(data_path)
    tlimit = 300
    if combo_idx == 0:  # Node weight spacing for CR
        nweights = [random.randint(1, 5) for _ in range(opt.g.n_nodes)]
        opt.g.add_node_weights(nweights)
        opt.m_val *= 3
        res = opt.optimize_layout(cutoff_time=tlimit, crossing_minimization=True, apply_node_weight_spacing=True)
    elif combo_idx == 1:  # Edge weight spacing for CR à la Zarate et al.
        eweights = {ed: random.randint(1, 5) for ed in opt.g.edge_ids}
        opt.g.add_edge_weights(eweights)
        res = opt.optimize_layout(cutoff_time=tlimit, crossing_minimization=True, apply_edge_weight=True)
    elif combo_idx == 2:  # Edge weight mult on BR
        improved_sifting(opt.g)
        eweights = {ed: random.randint(1, 5) for ed in opt.g.edge_ids}
        opt.g.add_edge_weights(eweights)
        res = opt.optimize_layout(cutoff_time=tlimit, bendiness_reduction=True, apply_edge_weight=True, fix_x_vars=True)
    elif combo_idx == 3:  # Node Emphasis + CR
        enodes = random.sample(list(opt.g.node_ids.keys()), max(round(0.05 * opt.g.n_nodes), 1))
        opt.g.add_node_emphasis(enodes)
        res = opt.optimize_layout(cutoff_time=tlimit, crossing_minimization=True, node_emphasis=True)
    elif combo_idx == 4:  # Node Emphasis + BR
        improved_sifting(opt.g)
        enodes = random.sample(list(opt.g.node_ids.keys()), max(round(0.05 * opt.g.n_nodes), 1))
        opt.g.add_node_emphasis(enodes)
        res = opt.optimize_layout(cutoff_time=tlimit, bendiness_reduction=True, fix_x_vars=True, node_emphasis=True)
    elif combo_idx == 5:  # Same Layer Edges + CR
        add_sl_edges(opt.g)
        res = opt.optimize_layout(cutoff_time=tlimit, crossing_minimization=True)
    elif combo_idx == 6:  # Same layer edges + BR
        improved_sifting(opt.g)
        add_sl_edges(opt.g)
        res = opt.optimize_layout(cutoff_time=tlimit, bendiness_reduction=True, fix_x_vars=True)
    elif combo_idx == 7:  # Straight edges + CR
        res = opt.optimize_layout(cutoff_time=tlimit, constrain_straight_long_arcs=True, crossing_minimization=True, vertical_transitivity=True)
    elif combo_idx == 8:  # Streamlining + BR
        improved_sifting(opt.g)
        res = opt.optimize_layout(cutoff_time=tlimit, streamline=True, bendiness_reduction=True, fix_x_vars=True)
    elif combo_idx == 9:  # Groups + CR + BR
        add_groups(opt.g)
        res = opt.optimize_layout(cutoff_time=tlimit, grouping_constraints=True, crossing_minimization=True, bendiness_reduction=True, y_based_group_constraints=True)
    elif combo_idx == 10:  # Groups + CR then Groups + BR
        add_groups(opt.g)
        res1 = opt.optimize_layout(cutoff_time=tlimit, grouping_constraints=True, crossing_minimization=True)
        res2 = opt.optimize_layout(cutoff_time=tlimit, grouping_constraints=True, bendiness_reduction=True, y_based_group_constraints=True, fix_x_vars=True)
        retval = collections.namedtuple("retval", "runtime objval status")
        res = retval(res1.runtime + res2.runtime, res1.objval + res2.objval, 2 if res1.status == res2.status == 2 else 9)
    elif combo_idx == 11:  # Fixed nodes + CR
        fix_nds = random.sample(list(opt.g.node_ids.keys()), max(round(0.2 * opt.g.n_nodes), 2))
        fix_top = random.sample(fix_nds, random.randint(0, max(round(0.2 * opt.g.n_nodes), 2)))
        fix_bottom = [v for v in fix_nds if v not in fix_top]
        opt.g.add_node_fix(fix_top, "top")
        opt.g.add_node_fix(fix_bottom, "bottom")
        res = opt.optimize_layout(cutoff_time=tlimit, crossing_minimization=True, fix_nodes=True)
    else:
        raise Exception("No experiment with that idx")

    gcr = opt.g.num_edge_crossings()
    bnds = sum(abs(e.n1.y - e.n2.y) for e in opt.g.edges)
    return sum(1 for nd in opt.g.nodes if not nd.is_anchor_node), opt.g.n_nodes, res.objval, gcr, bnds, res.runtime, res.status


def calculate_cutoff(csv_file, num_nodes, files_per_bucket):
    with open(csv_file, 'r') as fd:
        rdr = csv.reader(fd)
        first_line = next(rdr)
        fl_idx, st_idx = 1, first_line.index("Status")
        nfls = 0
        n_cutoff = 0
        for ln in rdr:
            if int(ln[fl_idx].split('_')[1]) == num_nodes:
                nfls += 1
                if int(ln[st_idx]) != 2:
                    n_cutoff += 1
        if nfls != files_per_bucket:
            raise Exception(f"Wrong num files in bucket {num_nodes}, {files_per_bucket} != {nfls}")
    return n_cutoff / nfls


generate_benchmark(["combo_idx"], {"combo_idx": list(range(12))}, run_func, "../random graphs/networkx2", name="design_constraints", csv_header=["Nodes", "TotalNodes", "ObjVal", "Crossings", "EdgeLength", "Runtime", "Status"], class_dependencies=["src/optimization.LayeredOptimizer"], project_root="/Users/connorwilson/PycharmProjects/stratisfimal-python")
