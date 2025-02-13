import unittest
from src import read_data, optimization, vis
from src.heuristics import improved_sifting


class TestAngularOptimization(unittest.TestCase):
    def setUp(self) -> None:
        self.g1 = read_data.read("../Rome-Lib/graficon18nodi/grafo198.18")
        self.g2 = read_data.read("../Rome-Lib/graficon32nodi/grafo1053.32")
        self.g3 = read_data.read("../Rome-Lib/graficon10nodi/grafo166.10")

    def test_angle_opt(self):
        opt = optimization.LayeredOptimizer(self.g1)
        opt.optimize_layout(crossing_minimization=True)
        opt.optimize_layout(angular_resolution=True, fix_x_vars=True, m_val=12)
        vis.draw_graph(self.g1, "ANG_TEST")

    def test_angle_opt_bland(self):
        opt = optimization.LayeredOptimizer(self.g3)
        opt.optimize_layout(angular_resolution=True)
        vis.draw_graph(self.g3, "ANG_TEST_0")

    def test_angle_opt_bigger(self):
        opt = optimization.LayeredOptimizer(self.g2)
        opt.optimize_layout(crossing_minimization=True)
        opt.optimize_layout(angular_resolution=True, fix_x_vars=True)
        vis.draw_graph(self.g2, "ANG_TEST_2")

    def test_angle_combined(self):
        opt = optimization.LayeredOptimizer(self.g1)
        opt.optimize_layout(crossing_minimization=True, angular_resolution=True, m_val=10)
        vis.draw_graph(self.g1, "ANG_TEST_3")

    def test_angle_combined_with_bend(self):
        opt = optimization.LayeredOptimizer(self.g1)
        opt.optimize_layout(crossing_minimization=True)
        opt.optimize_layout(angular_resolution=True, bendiness_reduction=True, gamma_2=0.1, fix_x_vars=True, m_val=12)
        vis.draw_graph(self.g1, "ANG_TEST_4")

    def test_fixed_x_sifting(self):
        opt = optimization.LayeredOptimizer(self.g2)
        improved_sifting(opt.g)
        vis.draw_graph(self.g2, "ANG_TEST_5a")
        opt.optimize_layout(fix_x_vars=True, angular_resolution=True)
        vis.draw_graph(self.g2, "ANG_TEST_5")
