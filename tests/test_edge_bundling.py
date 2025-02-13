import unittest
from src import read_data, optimization, vis


class TestOptimizationWithEdgeBundling(unittest.TestCase):
    def setUp(self) -> None:
        self.g1 = read_data.read("../Rome-Lib/graficon41nodi/grafo2019.41")

    def test_bundling(self):
        opt = optimization.LayeredOptimizer(self.g1)
        opt.optimize_layout(crossing_minimization=True, edge_bundling=True, gamma_1=4)
        vis.draw_graph(self.g1, "BUNDL_TEST")

    def test_bundling_bland(self):
        opt = optimization.LayeredOptimizer(self.g1)
        opt.optimize_layout(edge_bundling=True)
        vis.draw_graph(self.g1, "BUNDL_TEST_0")

    def test_cr_bundle_hybrid(self):
        opt = optimization.LayeredOptimizer(self.g1)
        res1 = opt.optimize_layout(crossing_minimization=True)
        opt.optimize_layout(hybrid_constraints=[("crossings", res1.objval + 2)], edge_bundling=True, bendiness_reduction=True, streamline=True, anchor_proximity=0.2)
        vis.draw_graph(self.g1, "BUNDL_TEST_2")
