import unittest
from src import read_data, optimization, vis


class TestOptimizationWithEdgeBundling(unittest.TestCase):
    def setUp(self) -> None:
        self.g1 = read_data.read("../Rome-Lib/graficon41nodi/grafo2019.41")

    def test_bundling(self):
        opt = optimization.LayeredOptimizer(self.g1, symmetry_breaking=True, edge_bundling=True, gamma_1=4)
        opt.optimize_layout()
        vis.draw_graph(self.g1, "BUNDL_TEST")
