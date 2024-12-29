import unittest
from src import read_data, optimization, vis


class TestMinMaxCrossingOptimization(unittest.TestCase):
    def setUp(self) -> None:
        self.g1 = read_data.read("../Rome-Lib/graficon71nodi/grafo6545.71")

    def test_minmax_only(self):
        opt = optimization.LayeredOptimizer(self.g1, symmetry_breaking=True, min_max_crossings=True, only_min_max_crossings=True)
        opt.optimize_layout()
        vis.draw_graph(self.g1, "MINMAX_TEST")
