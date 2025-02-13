import unittest
from src import read_data, optimization, vis


class TestOptimizationWithFairness(unittest.TestCase):
    def setUp(self) -> None:
        self.g1 = read_data.read("../Rome-Lib/graficon41nodi/grafo2019.41")
        self.fair_groups = {nid: 0 if nid <= self.g1.n_nodes // 2 else 1 for nid in self.g1.node_ids}
        self.g1.add_fairness_values(self.fair_groups)
        self.g2 = read_data.read("../Rome-Lib/graficon41nodi/grafo2531.41")
        self.fair_groups_2 = {nid: 0 if nid % 2 == 0 else 1 for nid in self.g2.node_ids}
        self.g2.add_fairness_values(self.fair_groups_2)

    def test_edge_length_fairness(self):
        opt = optimization.LayeredOptimizer(self.g1)
        opt.optimize_layout(fairness_constraints=True, fairness_metric="bends", crossing_minimization=False)
        vis.draw_graph(self.g1, "FAIR_BEND_TEST", groups=self.fair_groups)

        g_bends = [0, 0]
        len_g0, len_g1 = sum((1 for v in self.fair_groups.values() if v == 0)), sum((1 for v in self.fair_groups.values() if v == 1))
        for nd, v in self.fair_groups.items():
            for nd_ot in self.g1.get_adj_list()[nd]:
                if self.fair_groups[nd_ot] == v == 0:
                    g_bends[v] += abs(self.g1[nd].y - self.g1[nd_ot].y)
                elif self.fair_groups[nd_ot] == v == 1:
                    g_bends[v] += abs(self.g1[nd].y - self.g1[nd_ot].y)
        print(g_bends[0] / len_g0, g_bends[1] / len_g1)

    def test_edge_length_fairness_with_minimization(self):
        opt = optimization.LayeredOptimizer(self.g1)
        opt.optimize_layout(crossing_minimization=True)
        # opt.fairness_constraints, opt.fairness_metric, opt.crossing_minimization, opt.fix_x_vars, opt.bendiness_reduction = True, "bends", False, True, True
        opt.optimize_layout(fairness_constraints=True, fairness_metric="bends", fix_x_vars=True, bendiness_reduction=True)
        vis.draw_graph(self.g1, "FAIR_BEND_TEST_2", groups=self.fair_groups)

        g_bends = [0, 0]
        len_g0, len_g1 = sum((1 for v in self.fair_groups.values() if v == 0)), sum((1 for v in self.fair_groups.values() if v == 1))
        seen = set()
        for nd, v in self.fair_groups.items():
            for nd_ot in self.g1.get_adj_list()[nd]:
                if nd_ot not in seen:
                    if self.fair_groups[nd_ot] == v == 0:
                        g_bends[v] += abs(self.g1[nd].y - self.g1[nd_ot].y)
                    elif self.fair_groups[nd_ot] == v == 1:
                        g_bends[v] += abs(self.g1[nd].y - self.g1[nd_ot].y)
            seen.add(nd)
        print(g_bends)
        print(g_bends[0] / len_g0, g_bends[1] / len_g1)

    def test_crossing_fairness(self):
        opt = optimization.LayeredOptimizer(self.g2)
        opt.optimize_layout(fairness_constraints=True, fairness_metric="crossings", crossing_minimization=True, gamma_fair=5)
        opt.optimize_layout(bendiness_reduction=True, fix_x_vars=True)
        vis.draw_graph(self.g2, "FAIR_CR_TEST", groups=self.fair_groups_2)

    def test_crossing_fairness_no_min(self):
        opt = optimization.LayeredOptimizer(self.g2)
        opt.optimize_layout(fairness_constraints=True, fairness_metric="crossings")
        opt.optimize_layout(bendiness_reduction=True, fix_x_vars=True)
        vis.draw_graph(self.g2, "FAIR_CR_TEST_2", groups=self.fair_groups_2)

        g_cr = [0, 0]
        len_g0, len_g1 = sum((1 for v in self.fair_groups_2.values() if v == 0)), sum((1 for v in self.fair_groups_2.values() if v == 1))
        e_b_l = self.g2.get_edge_ids_by_layer()
        for l in range(self.g2.n_layers - 1):
            for i, e1 in enumerate(e_b_l[l]):
                for e2 in e_b_l[l][i+1:]:
                    if (self.g2[e1[0]].y > self.g2[e2[0]].y and self.g2[e1[1]].y < self.g2[e2[1]].y) or (self.g2[e1[0]].y < self.g2[e2[0]].y and self.g2[e1[1]].y > self.g2[e2[1]].y):
                        fsum = self.fair_groups_2[e1[0]] + self.fair_groups_2[e1[1]] + self.fair_groups_2[e2[0]] + self.fair_groups_2[e2[1]]
                        if fsum <= 1:
                            g_cr[0] += 1
                        elif fsum >= 3:
                            g_cr[1] += 1

        print(g_cr)
        print(g_cr[0], g_cr[1] * len_g0 / len_g1)

    def test_fairness_hybrid(self):
        opt = optimization.LayeredOptimizer(self.g2)
        res1 = opt.optimize_layout(crossing_minimization=True)
        opt.optimize_layout(fairness_constraints=True, fairness_metric="crossings", hybrid_constraints=[("crossings", res1.objval + 4)])
        res2 = opt.optimize_layout(bendiness_reduction=True, fix_x_vars=True)
        opt.optimize_layout(fairness_constraints=True, fix_x_vars=True, fairness_metric="bends", hybrid_constraints=[("bends", res2.objval + 5)])
        vis.draw_graph(self.g2, "FAIR_BOTH_TEST", groups=self.fair_groups_2)
