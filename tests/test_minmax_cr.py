import unittest
from src import read_data, optimization, vis


class TestMinMaxCrossingOptimization(unittest.TestCase):
    def setUp(self) -> None:
        self.g1 = read_data.read("../Rome-Lib/graficon71nodi/grafo6545.71")
        self.g2 = read_data.read("../Rome-Lib/graficon60nodi/grafo3199.60")

    def test_minmax_only(self):
        opt = optimization.LayeredOptimizer(self.g2)
        opt.optimize_layout(min_max_crossings=True, crossing_minimization=False)
        vis.draw_graph(self.g2, "MINMAX_TEST")

    def test_minmax_plus_cr(self):
        opt = optimization.LayeredOptimizer(self.g1)
        opt.optimize_layout(min_max_crossings=True, crossing_minimization=True, gamma_min_max=3)
        opt.optimize_layout(bendiness_reduction=True, fix_x_vars=True)
        vis.draw_graph(self.g1, "MINMAX_TEST_2")

    def test_minmax_hybrid_cr(self):
        opt = optimization.LayeredOptimizer(self.g2)
        res = opt.optimize_layout(min_max_crossings=True)
        opt.optimize_layout(crossing_minimization=True, hybrid_constraints=[("min_max_crossings", res.objval)])
        opt.optimize_layout(bendiness_reduction=True, fix_x_vars=True)
        vis.draw_graph(self.g2, "MINMAX_TEST_3")
