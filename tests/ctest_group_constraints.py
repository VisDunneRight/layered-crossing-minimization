import os
import random
import unittest
from src import read_data, optimization, vis


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
            for _ in range(random.randint(gp_target, 2 * gp_target)):
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
    print(groups)
    g.add_groups(groups)
    # vis.draw_graph(g, "testing")


class TestOptimizationWithGroups(unittest.TestCase):
    def setUp(self) -> None:
        self.g1 = read_data.read("../Rome-Lib/graficon71nodi/grafo6545.71")
        self.g1.add_groups([[30, 37, 1, 25], [57, 43], [52, 9, 24]])
        self.g2 = read_data.read("../Rome-Lib/graficon71nodi/grafo6545.71")
        self.g2.add_groups([[31, 38, 18, 27], [2, 65, 12, 28]])
        self.g3 = read_data.read("../Rome-Lib/graficon71nodi/grafo6545.71")
        self.g3.add_groups([[30, 37, 1, 25], [52, 9, 24], [31, 38, 18, 27], [2, 65, 12, 28], [57, 43, 19], [41, 60]])
        self.g4 = read_data.read("../Rome-Lib/graficon71nodi/grafo6545.71")

    def test_sl_groups(self):
        opt = optimization.LayeredOptimizer(self.g1)
        opt.optimize_layout(grouping_constraints=True, crossing_minimization=True)
        vis.draw_graph(self.g1, "GRP_TEST")

    def test_ml_groups(self):
        opt = optimization.LayeredOptimizer(self.g2)
        opt.optimize_layout(grouping_constraints=True, crossing_minimization=True)
        vis.draw_graph(self.g2, "GRP_TEST_2")

    def test_all_groups_hard(self):
        opt = optimization.LayeredOptimizer(self.g3)
        opt.optimize_layout(grouping_constraints=True, crossing_minimization=True)
        vis.draw_graph(self.g3, "GRP_TEST_3")

    def test_all_groups_hard_with_y_vals(self):
        opt = optimization.LayeredOptimizer(self.g3)
        opt.optimize_layout(grouping_constraints=True, crossing_minimization=True, y_based_group_constraints=True)
        vis.draw_graph(self.g3, "GRP_TEST_3b")

    def test_all_groups_hard_with_bendiness(self):
        opt = optimization.LayeredOptimizer(self.g3)
        opt.optimize_layout(grouping_constraints=True, crossing_minimization=True)
        vis.draw_graph(self.g3, "GRP_TEST_3bendA")
        opt.optimize_layout(bendiness_reduction=True, grouping_constraints=True, y_based_group_constraints=True, fix_x_vars=True)
        vis.draw_graph(self.g3, "GRP_TEST_3bendB")

    def test_group_correctness(self):
        opt = optimization.LayeredOptimizer(self.g4)
        res1 = opt.optimize_layout(crossing_minimization=True)
        vis.draw_graph(self.g4, "GRP_TEST_4a")
        grouped = set()
        groups = []
        for _ in range(10):
            rn = random.choice(self.g4.nodes)
            lv, yv = rn.layer, rn.y
            if rn.id not in grouped:
                groups.append([rn.id])
                grouped.add(rn.id)
                for lid, py in [(lv, yv - 1), (lv, yv + 1), (lv - 1, yv - 1), (lv - 1, yv), (lv - 1, yv + 1), (lv + 1, yv - 1), (lv + 1, yv), (lv + 1, yv + 1)]:
                    if lid in self.g4.layers:
                        for nd in self.g4.layers[lid]:
                            if nd.y == py and nd.id not in grouped:
                                groups[-1].append(nd.id)
                                grouped.add(nd.id)
                                continue
        self.g4.add_groups(groups)
        opt2 = optimization.LayeredOptimizer(self.g4)
        res2 = opt2.optimize_layout(grouping_constraints=True, crossing_minimization=True)
        vis.draw_graph(self.g4, "GRP_TEST_4b")
        self.assertEqual(res1.objval, res2.objval)

    def test_stress_crbend(self):
        for i in range(5, 15):
            for j in range(50):
                print(f"graph_{i}_{j}")
                opt = optimization.LayeredOptimizer(f"../random graphs/networkx2/graph_{i}_{j}")
                add_groups(opt.g)
                opt.optimize_layout(grouping_constraints=True, crossing_minimization=True)
                vis.draw_graph(opt.g, "GRP_TEST_5bendA")
                opt.optimize_layout(bendiness_reduction=True, grouping_constraints=True, y_based_group_constraints=True, fix_x_vars=True)
                vis.draw_graph(opt.g, "GRP_TEST_5bendB")
        self.assertEqual(0, 0)

    def test_crbend_specific(self):
        opt = optimization.LayeredOptimizer("../random graphs/networkx2/graph_8_9")
        opt.g.add_groups({5: 0, 0: 0, 3: 0, 7: 0})
        opt.optimize_layout(grouping_constraints=True, crossing_minimization=True)
        vis.draw_graph(opt.g, "GRP_TEST_5bendA")
        print([(nd.id, nd.y) for nd in opt.g.nodes])
        opt.optimize_layout(bendiness_reduction=True, grouping_constraints=True, y_based_group_constraints=True,
                            fix_x_vars=True)
        vis.draw_graph(opt.g, "GRP_TEST_5bendB")
