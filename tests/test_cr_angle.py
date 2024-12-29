import unittest
from src import read_data, optimization, vis


class TestAngularOptimization(unittest.TestCase):
    def setUp(self) -> None:
        self.g1 = read_data.read("../Rome-Lib/graficon71nodi/grafo6545.71")

    def test_angle_opt(self):
        opt = optimization.LayeredOptimizer(self.g1, symmetry_breaking=True, angular_resolution=True)
        opt.optimize_layout()
        vis.draw_graph(self.g1, "ANG_TEST")
