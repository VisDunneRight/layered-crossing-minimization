import unittest
from layered_optimization import read_data, optimization, vis, metrics


class TestMinMaxCrossingOptimization(unittest.TestCase):
    def setUp(self) -> None:
        self.g1 = read_data.read("../datasets/Rome-Lib/graficon71nodi/grafo6545.71")
        self.g2 = read_data.read("../datasets/Rome-Lib/graficon60nodi/grafo3199.60")

    def test_minmax_only(self):
        opt = optimization.LayeredOptimizer(self.g2)
        res = opt.optimize_layout(min_max_crossings=True, crossing_minimization=False)
        print(res.metrics)
        self.assertEqual(res.metrics.minmax, metrics.min_max_crossings(self.g2))
        vis.draw_graph(self.g2, "MINMAX_TEST")

    def test_minmax_plus_cr(self):
        opt = optimization.LayeredOptimizer(self.g1)
        res = opt.optimize_layout(min_max_crossings=True, crossing_minimization=True, gamma_min_max=3)
        self.assertEqual(res.metrics.minmax, metrics.min_max_crossings(self.g1))
        opt.optimize_layout(edge_length_minimization=True, fix_x_vars=True)
        vis.draw_graph(self.g1, "MINMAX_TEST_2")

    def test_minmax_hybrid_cr(self):
        opt = optimization.LayeredOptimizer(self.g2)
        res = opt.optimize_layout(min_max_crossings=True)
        res2 = opt.optimize_layout(crossing_minimization=True, hybrid_constraints=[("min_max_crossings", res.objval)])
        res3 = opt.optimize_layout(edge_length_minimization=True, fix_x_vars=True)
        vis.draw_graph(self.g2, "MINMAX_TEST_3")
        self.assertEqual(res.metrics.minmax, res2.metrics.minmax, metrics.min_max_crossings(self.g2))
