import unittest
from src import read_data, optimization, vis


class TestMinEdgesWithCrossings(unittest.TestCase):
    def setUp(self) -> None:
        self.g1 = read_data.read("../Rome-Lib/graficon71nodi/grafo6545.71")
        self.g2 = read_data.read("../Rome-Lib/graficon60nodi/grafo3199.60")

    def test_min_edges_only(self):
        opt = optimization.LayeredOptimizer(self.g2)
        opt.optimize_layout(min_edges_with_crossings=True, crossing_minimization=False)
        opt.optimize_layout(bendiness_reduction=True, fix_x_vars=True)
        vis.draw_graph(self.g2, "EWC_TEST")

    def test_min_edges_hybrid_cr(self):
        opt = optimization.LayeredOptimizer(self.g2)
        res = opt.optimize_layout(min_edges_with_crossings=True)
        opt.optimize_layout(crossing_minimization=True, hybrid_constraints=[("min_edges_with_crossings", res.objval)])
        opt.optimize_layout(bendiness_reduction=True, fix_x_vars=True)
        vis.draw_graph(self.g2, "EWC_TEST_3")
