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

    def test_crmin(self):
        opt = optimization.LayeredOptimizer(self.g1, symmetry_breaking=True, heuristic_start=True)
        res = opt.optimize_layout()
        self.assertEqual(27, res[1])
        self.assertEqual(res[1], self.g1.num_edge_crossings())

    def test_crmin_plus_bend_combo(self):
        opt = optimization.LayeredOptimizer(self.g_bend, symmetry_breaking=True, bendiness_reduction=True, sequential_bendiness=False, gamma_1=3, stratisfimal_y_vars=False)
        res = opt.optimize_layout()
        vis.draw_graph(self.g_bend, "Bend_is_wrong")
        self.assertEqual(48, res[1])
        self.assertEqual(res[1], self.g_bend.calculate_stratisfimal_objective(3, 1))

    def test_crmin_plus_bend_sequential(self):
        opt_0 = optimization.LayeredOptimizer(self.g1, symmetry_breaking=True)
        opt_0.optimize_layout()
        basic_objective = self.g1.calculate_stratisfimal_objective(3, 1)
        opt = optimization.LayeredOptimizer(self.g1, symmetry_breaking=True, bendiness_reduction=True, sequential_bendiness=True)
        res = opt.optimize_layout()
        self.assertEqual(27, res[1])
        self.assertGreater(basic_objective, self.g1.calculate_stratisfimal_objective(3, 1))

    def test_weird_crmin_sequential_no_bend(self):
        opt = optimization.LayeredOptimizer(self.g1, symmetry_breaking=True, bendiness_reduction=False, sequential_bendiness=True)
        res = opt.optimize_layout()
        self.assertEqual(27, res[1])
        self.assertEqual(res[1], self.g1.num_edge_crossings())

    def test_mirrored_vars_crmin(self):
        opt = optimization.LayeredOptimizer(self.g1, symmetry_breaking=True, mirror_vars=True)
        res = opt.optimize_layout()
        self.assertEqual(27, res[1])
        self.assertEqual(res[1], self.g1.num_edge_crossings())

    def test_same_layer_edges(self):
        opt = optimization.LayeredOptimizer(self.g2, symmetry_breaking=True)
        res = opt.optimize_layout()
        self.assertEqual(28, res[1])

    def test_all_groups_hard(self):
        opt = optimization.LayeredOptimizer(self.g3, symmetry_breaking=True, grouping_constraints=True)
        res = opt.optimize_layout()
        self.assertEqual(41, res[1])
        self.assertEqual(res[1], self.g3.num_edge_crossings())

    def test_all_groups_hard_with_bendiness(self):
        opt = optimization.LayeredOptimizer(self.g3, symmetry_breaking=True, grouping_constraints=True, sequential_bendiness=True, bendiness_reduction=True, sequential_grouping_constraints=True)
        res = opt.optimize_layout()
        self.assertEqual(41, res[1])
        self.assertEqual(res[1], self.g3.num_edge_crossings())
