import unittest
from layered_optimization import read_data, optimization, vis


class TestOptimizationWithEdgeBundling(unittest.TestCase):
    def setUp(self) -> None:
        self.g1 = read_data.read("../datasets/Rome-Lib/graficon18nodi/grafo198.18")
        self.g2 = read_data.read("../datasets/Rome-Lib/graficon32nodi/grafo1053.32")

    def test1_edge_length_solo(self):
        opt = optimization.LayeredOptimizer(self.g1)
        opt.optimize_layout(edge_length_minimization=True)
        vis.draw_graph(self.g1, "EDGELEN_TEST_SOLO")

    def test2_edge_length_fixed(self):
        opt = optimization.LayeredOptimizer(self.g2)
        opt.optimize_layout(edge_length_minimization=True, fix_x_vars=True, streamline=True)
        vis.draw_graph(self.g2, "EDGELEN_TEST_FIXX")
