import unittest
from src import read_data, optimization, vis


class TestSequentialOptimization(unittest.TestCase):
    def setUp(self) -> None:
        self.g1 = read_data.read("../Rome-Lib/graficon71nodi/grafo6545.71")

    def test_cr_then_bend(self):
        opt = optimization.LayeredOptimizer(self.g1)
        opt.optimize_layout(crossing_minimization=True)
        opt.optimize_layout(bendiness_reduction=True, fix_x_vars=True)
        vis.draw_graph(self.g1, "SEQ_TEST_1")

    def test_cr_then_bend_streamline(self):
        opt = optimization.LayeredOptimizer(self.g1)
        opt.optimize_layout(crossing_minimization=True)
        opt.optimize_layout(bendiness_reduction=True, fix_x_vars=True, streamline=True)
        vis.draw_graph(self.g1, "SEQ_TEST_2")

    def test_hybrid_crossing_then_bend(self):
        opt = optimization.LayeredOptimizer(self.g1)
        res = opt.optimize_layout(crossing_minimization=True)
        # opt.optimize_layout(bendiness_reduction=True, crossing_minimization=True, start_xy_vars=True)
        opt.optimize_layout(bendiness_reduction=True, hybrid_constraints=[("crossings", 4)], start_xy_vars=True)
        vis.draw_graph(self.g1, "SEQ_TEST_3")
