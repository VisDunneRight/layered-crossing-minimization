import unittest
from layered_optimization import read_data, optimization, vis, metrics


class TestSymmetryMaximization(unittest.TestCase):
    def setUp(self) -> None:
        self.g1 = read_data.read("../datasets/Rome-Lib/graficon41nodi/grafo2019.41")
        self.g2 = read_data.read("../datasets/Rome-Lib/graficon13nodi/grafo201.13")
        self.g3 = read_data.read("../datasets/Rome-Lib/graficon13nodi/grafo2221.13")

    def test_cr_then_node_sym(self):
        opt = optimization.LayeredOptimizer(self.g1)
        opt.optimize_layout(crossing_minimization=True)
        res = opt.optimize_layout(horizontal_node_symmetry=True, fix_x_vars=True)
        self.assertEqual(res.metrics.horizontal_node_symmetry, metrics.horizontal_node_symmetry(self.g1))
        vis.draw_graph(self.g1, "SYM_TEST_1")

    def test_cr_then_edge_sym(self):
        opt = optimization.LayeredOptimizer(self.g2)
        opt.optimize_layout(crossing_minimization=True)
        res = opt.optimize_layout(horizontal_edge_symmetry=True, fix_x_vars=True)
        vis.draw_graph(self.g2, "SYM_TEST_2")
        self.assertEqual(res.metrics.horizontal_edge_symmetry, metrics.horizontal_edge_symmetry(self.g2))

    def test_cr_then_sym_plus_bend(self):
        opt = optimization.LayeredOptimizer(self.g2)
        opt.optimize_layout(crossing_minimization=True)
        opt.optimize_layout(horizontal_node_symmetry=True, horizontal_edge_symmetry=True, fix_x_vars=True, edge_length_minimization=True, streamline=True)
        vis.draw_graph(self.g2, "SYM_TEST_3")

    def test_only_edge_sym_bend(self):
        opt = optimization.LayeredOptimizer(self.g1)
        res = opt.optimize_layout(horizontal_edge_symmetry=True, streamline=True)
        self.assertEqual(res.metrics.horizontal_edge_symmetry, metrics.horizontal_edge_symmetry(self.g1))
        vis.draw_graph(self.g1, "SYM_TEST_4")

    def test_vertical_sym(self):
        opt = optimization.LayeredOptimizer(self.g1)
        res = opt.optimize_layout(vertical_node_symmetry=True, vertical_edge_symmetry=True, streamline=True, symmetry_vertical_axis=4)
        vis.draw_graph(self.g1, "SYM_TEST_5")
        self.assertEqual(res.metrics.vertical_node_symmetry, metrics.vertical_node_symmetry(self.g1))
        self.assertEqual(res.metrics.vertical_edge_symmetry, metrics.vertical_edge_symmetry(self.g1))

        opt2 = optimization.LayeredOptimizer(self.g2)
        res2 = opt2.optimize_layout(vertical_node_symmetry=True, vertical_edge_symmetry=True, streamline=True)
        self.assertEqual(res2.metrics.vertical_node_symmetry, metrics.vertical_node_symmetry(self.g2, v_axis=2))
        self.assertEqual(res2.metrics.vertical_edge_symmetry, metrics.vertical_edge_symmetry(self.g2))

        opt3 = optimization.LayeredOptimizer(self.g3)
        res3 = opt3.optimize_layout(vertical_node_symmetry=True, vertical_edge_symmetry=True, streamline=True)
        vis.draw_graph(self.g3, "SYM_TEST_5VB")
        self.assertEqual(res3.metrics.vertical_node_symmetry, metrics.vertical_node_symmetry(self.g3, v_axis=1.5))
        self.assertEqual(res3.metrics.vertical_edge_symmetry, metrics.vertical_edge_symmetry(self.g3))

    def test_rotation_sym(self):
        opt = optimization.LayeredOptimizer(self.g2)
        res = opt.optimize_layout(rotation_node_symmetry=True, streamline=True)
        self.assertEqual(res.metrics.rotation_node_symmetry, metrics.rotation_node_symmetry(self.g2))
        vis.draw_graph(self.g2, "SYM_TEST_6")

        opt2 = optimization.LayeredOptimizer(self.g3)
        res2 = opt2.optimize_layout(rotation_node_symmetry=True, streamline=True)
        self.assertEqual(res2.metrics.rotation_node_symmetry, metrics.rotation_node_symmetry(self.g3))

        opt3 = optimization.LayeredOptimizer(self.g3)
        res3 = opt3.optimize_layout(rotation_edge_symmetry=True, edge_length_minimization=False, streamline=True)
        vis.draw_graph(self.g3, "SYM_TEST_6VB")
        self.assertEqual(res3.metrics.actual_symmetries, metrics.rotation_edge_symmetry(self.g3, return_act_sym=True))
        self.assertEqual(res3.metrics.rotation_edge_symmetry, metrics.rotation_edge_symmetry(self.g3))

    def test_rotation_sym_both(self):
        opt = optimization.LayeredOptimizer(self.g2)
        res = opt.optimize_layout(rotation_node_symmetry=True, rotation_edge_symmetry=True, streamline=True)
        self.assertEqual(res.metrics.rotation_edge_symmetry, metrics.rotation_edge_symmetry(self.g2))
        self.assertEqual(res.metrics.rotation_node_symmetry, metrics.rotation_node_symmetry(self.g2))
        vis.draw_graph(self.g2, "SYM_TEST_7")

        opt2 = optimization.LayeredOptimizer(self.g3)
        res2 = opt2.optimize_layout(rotation_node_symmetry=True, rotation_edge_symmetry=True, streamline=True)
        self.assertEqual(res2.metrics.rotation_edge_symmetry, metrics.rotation_edge_symmetry(self.g3))
        self.assertEqual(res2.metrics.rotation_node_symmetry, metrics.rotation_node_symmetry(self.g3))
