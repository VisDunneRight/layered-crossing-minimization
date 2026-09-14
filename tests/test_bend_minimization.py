import unittest
from layered_optimization import read_data, optimization, vis


class TestBendMinimization(unittest.TestCase):
    def setUp(self) -> None:
        self.g1 = read_data.read("../Rome-Lib/graficon71nodi/grafo6545.71")
        self.g2 = read_data.read("../Rome-Lib/graficon60nodi/grafo3199.60")

    def test_min_bends_only(self):
        opt = optimization.LayeredOptimizer(self.g2)
        opt.optimize_layout(crossing_minimization=True)
        opt.optimize_layout(bend_minimization=True, fix_x_vars=True, streamline=True)
        vis.draw_graph(self.g2, "BENDMIN_TEST", straighten_edges=False)

    def test_min_bends_with_edgelength(self):
        opt = optimization.LayeredOptimizer(self.g2)
        opt.optimize_layout(crossing_minimization=True)
        opt.optimize_layout(edge_length_minimization=True, bend_minimization=True, fix_x_vars=True)
        vis.draw_graph(self.g2, "BENDMIN_TEST_2")
