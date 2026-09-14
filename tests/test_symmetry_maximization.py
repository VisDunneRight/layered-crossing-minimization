import unittest
from layered_optimization import read_data, optimization, vis


class TestSymmetryMaximization(unittest.TestCase):
    def setUp(self) -> None:
        self.g1 = read_data.read("../Rome-Lib/graficon41nodi/grafo2019.41")
        self.g2 = read_data.read("../Rome-Lib/graficon13nodi/grafo201.13")
        self.g3 = read_data.read("../Rome-Lib/graficon13nodi/grafo2221.13")

    def test_cr_then_node_sym(self):
        opt = optimization.LayeredOptimizer(self.g1)
        opt.optimize_layout(crossing_minimization=True)
        opt.optimize_layout(symmetry_maximization=True, fix_x_vars=True)
        vis.draw_graph(self.g1, "SYM_TEST_1")

    def test_cr_then_node_and_edge_sym(self):
        opt = optimization.LayeredOptimizer(self.g2)
        opt.optimize_layout(crossing_minimization=True)
        opt.optimize_layout(symmetry_maximization=True, symmetry_maximization_edges=True, fix_x_vars=True)
        vis.draw_graph(self.g2, "SYM_TEST_2")

    def test_cr_then_sym_plus_bend(self):
        opt = optimization.LayeredOptimizer(self.g2)
        opt.optimize_layout(crossing_minimization=True)
        opt.optimize_layout(symmetry_maximization=True, symmetry_maximization_edges=True, fix_x_vars=True, edge_length_minimization=True, streamline=True)
        vis.draw_graph(self.g2, "SYM_TEST_3")

    def test_only_edge_sym_bend(self):
        opt = optimization.LayeredOptimizer(self.g1)
        opt.optimize_layout(symmetry_maximization_edges=True, streamline=True)
        vis.draw_graph(self.g1, "SYM_TEST_4")

    def test_vertical_sym(self):
        opt = optimization.LayeredOptimizer(self.g1)
        opt.optimize_layout(symmetry_maximization=True, symmetry_maximization_edges=True, symmetry_maximization_type="vertical", streamline=True)
        vis.draw_graph(self.g1, "SYM_TEST_5")

        opt2 = optimization.LayeredOptimizer(self.g2)
        opt2.optimize_layout(symmetry_maximization=True, symmetry_maximization_edges=True, symmetry_maximization_type="vertical", streamline=True)

        opt3 = optimization.LayeredOptimizer(self.g3)
        opt3.optimize_layout(symmetry_maximization=True, symmetry_maximization_edges=True, symmetry_maximization_type="vertical", streamline=True)

    def test_rotation_sym(self):
        opt = optimization.LayeredOptimizer(self.g2)
        opt.optimize_layout(symmetry_maximization=True, symmetry_maximization_type="rotation", streamline=True)
        vis.draw_graph(self.g2, "SYM_TEST_6")

        opt2 = optimization.LayeredOptimizer(self.g3)
        opt2.optimize_layout(symmetry_maximization=True, symmetry_maximization_type="rotation", streamline=True)

    def test_rotation_sym_both(self):
        opt = optimization.LayeredOptimizer(self.g2)
        opt.optimize_layout(symmetry_maximization=True, symmetry_maximization_edges=True, symmetry_maximization_type="rotation", streamline=True)
        vis.draw_graph(self.g2, "SYM_TEST_7")

        opt2 = optimization.LayeredOptimizer(self.g3)
        opt2.optimize_layout(symmetry_maximization=True, symmetry_maximization_edges=True, symmetry_maximization_type="rotation", streamline=True)
