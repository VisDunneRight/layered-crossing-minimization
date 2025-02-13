import unittest
from src import read_data, optimization, vis


class TestCorrectnessAllDesigns(unittest.TestCase):
    def setUp(self) -> None:
        self.g1 = read_data.read("../Rome-Lib/graficon71nodi/grafo6545.71")
        self.g2 = read_data.read("../Rome-Lib/graficon71nodi/grafo6545.71", remove_sl=False)
        self.g2.add_edge(63, 64)
        self.g2.add_edge(8, 64)
        self.g2.add_edge(8, 34)
        self.g3 = read_data.read("../Rome-Lib/graficon71nodi/grafo6545.71")
        self.g3.add_groups([[30, 37, 1, 25], [52, 9, 24], [31, 38, 18, 27], [2, 65, 12, 28], [57, 43, 19], [41, 60]])
        self.g_bend = read_data.read("../Rome-Lib/graficon41nodi/grafo2019.41")
        self.g4 = read_data.read("../Rome-Lib/graficon18nodi/grafo198.18")
        self.g5 = read_data.read("../Rome-Lib/graficon13nodi/grafo201.13")
        self.g6 = read_data.read("../Rome-Lib/graficon41nodi/grafo2019.41")
        self.fair_groups = {nid: 0 if nid <= self.g1.n_nodes // 2 else 1 for nid in self.g6.node_ids}
        self.g6.add_fairness_values(self.fair_groups)

    def test_crmin(self):
        opt = optimization.LayeredOptimizer(self.g1)
        res = opt.optimize_layout(crossing_minimization=True)
        self.assertEqual(27, res.objval)
        self.assertEqual(res.objval, self.g1.num_edge_crossings())

    def test_crmin_plus_bend_combo(self):
        opt = optimization.LayeredOptimizer(self.g_bend)
        res = opt.optimize_layout(crossing_minimization=True, bendiness_reduction=True, sequential_bendiness=False, gamma_1=3)
        self.assertAlmostEqual(48, res.objval, delta=0.0001)
        self.assertAlmostEqual(res.objval, self.g_bend.calculate_stratisfimal_objective(3, 1), delta=0.0001)

    def test_crmin_plus_bend_sequential(self):
        opt_0 = optimization.LayeredOptimizer(self.g1)
        opt_0.optimize_layout(crossing_minimization=True)
        basic_objective = self.g1.calculate_stratisfimal_objective(3, 1)
        opt = optimization.LayeredOptimizer(self.g1)
        res1 = opt.optimize_layout(crossing_minimization=True)
        opt.optimize_layout(bendiness_reduction=True, fix_x_vars=True)
        self.assertEqual(27, res1.objval)
        self.assertGreater(basic_objective, self.g1.calculate_stratisfimal_objective(3, 1))

    def test_same_layer_edges(self):
        opt = optimization.LayeredOptimizer(self.g2)
        res = opt.optimize_layout(crossing_minimization=True)
        self.assertEqual(28, res.objval)

    def test_all_groups_hard(self):
        opt = optimization.LayeredOptimizer(self.g3)
        res = opt.optimize_layout(crossing_minimization=True, grouping_constraints=True)
        self.assertEqual(41, res.objval)
        self.assertEqual(res.objval, self.g3.num_edge_crossings())

    def test_all_groups_hard_with_bendiness(self):
        opt = optimization.LayeredOptimizer(self.g3)
        opt.optimize_layout(crossing_minimization=True, grouping_constraints=True)
        res = opt.optimize_layout(grouping_constraints=True, fix_x_vars=True, bendiness_reduction=True, sequential_grouping_constraints=True)
        vis.draw_graph(self.g3, "CORR_HARD")
        self.assertAlmostEqual(344, res.objval, delta=0.0001)
        self.assertEqual(41, self.g3.num_edge_crossings())

    def test_angle_opt(self):
        opt = optimization.LayeredOptimizer(self.g4)
        opt.optimize_layout(crossing_minimization=True)
        res = opt.optimize_layout(angular_resolution=True, fix_x_vars=True, m_val=12)
        self.assertEqual(0, res.objval)
        for e1 in self.g4.edges:
            for e2 in self.g4.edges:
                if e1.n1.layer == e2.n1.layer:
                    if (e1.n1.y > e2.n1.y and e1.n2.y < e2.n2.y) or (e1.n1.y < e2.n1.y and e1.n2.y > e2.n2.y):
                        self.assertEqual(-1, round((e1.n2.y - e1.n1.y)/2 * (e2.n2.y - e2.n1.y)/2, 6))

    def test_sym_nodes_and_edges(self):
        opt = optimization.LayeredOptimizer(self.g5)
        opt.optimize_layout(crossing_minimization=True)
        res = opt.optimize_layout(symmetry_maximization=True, symmetry_maximization_edges=True, fix_x_vars=True)
        sum_sym = 0
        for ndlist in self.g5.get_ids_by_layer().values():
            for i, nd1 in enumerate(ndlist):
                for nd2 in ndlist[i:]:
                    if round(self.g5[nd1].y + self.g5[nd2].y, 4) != opt.m_val:
                        sum_sym += 1
                    # else:
                    #     print(nd1, nd2, "are symmetric")
        for elist in self.g5.get_edges_by_layer().values():
            for i, e1 in enumerate(elist):
                for e2 in elist[i:]:
                    if e1 != e2 and (round(self.g5[e1.n1.id].y + self.g5[e2.n1.id].y, 4) != opt.m_val or round(self.g5[e1.n2.id].y + self.g5[e2.n2.id].y, 4) != opt.m_val):
                        sum_sym += 1
                    # else:
                    #     print(e1, e2, "are symmetric")
        self.assertEqual(res.objval, sum_sym)

    def test_edge_length_fairness_with_minimization(self):
        opt = optimization.LayeredOptimizer(self.g6)
        opt.optimize_layout(crossing_minimization=True)
        opt.optimize_layout(fairness_constraints=True, fairness_metric="bends", fix_x_vars=True, bendiness_reduction=True)
        g_bends = [0, 0]
        len_g0, len_g1 = sum((1 for v in self.fair_groups.values() if v == 0)), sum((1 for v in self.fair_groups.values() if v == 1))
        seen = set()
        for nd, v in self.fair_groups.items():
            for nd_ot in self.g6.get_adj_list()[nd]:
                if nd_ot not in seen:
                    if self.fair_groups[nd_ot] == v == 0:
                        g_bends[v] += abs(self.g6[nd].y - self.g6[nd_ot].y)
                    elif self.fair_groups[nd_ot] == v == 1:
                        g_bends[v] += abs(self.g6[nd].y - self.g6[nd_ot].y)
            seen.add(nd)
        self.assertEqual(round(g_bends[0] / len_g0, 6), round(g_bends[1] / len_g1, 6))
