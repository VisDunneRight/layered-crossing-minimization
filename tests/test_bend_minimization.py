import unittest
from layered_optimization import read_data, optimization, vis
from layered_optimization.heuristics import improved_sifting


class TestBendMinimization(unittest.TestCase):
    def setUp(self) -> None:
        self.g1 = read_data.read("../datasets/Rome-Lib/graficon71nodi/grafo6545.71")
        self.g2 = read_data.read("../datasets/Rome-Lib/graficon60nodi/grafo3199.60")
        self.g3 = read_data.read("../datasets/random_graphs/uniform_layered/n5/graph112.lgbin")

    def test1_min_bends_only(self):
        opt = optimization.LayeredOptimizer(self.g2)
        opt.optimize_layout(crossing_minimization=True)
        opt.optimize_layout(bend_minimization=True, fix_x_vars=True, streamline=True)
        vis.draw_graph(self.g2, "BENDMIN_TEST", straighten_edges=False)

    def test2_min_bends_with_edgelength(self):
        opt = optimization.LayeredOptimizer(self.g2)
        opt.optimize_layout(crossing_minimization=True)
        opt.optimize_layout(edge_length_minimization=True, bend_minimization=True, fix_x_vars=True, gamma_edgelength=0.5)
        vis.draw_graph(self.g2, "BENDMIN_TEST_2")

    def test3_fixed_x_solo(self):
        opt = optimization.LayeredOptimizer(self.g2)
        improved_sifting(opt.g)
        opt.optimize_layout(bend_minimization=True, fix_x_vars=True)
        vis.draw_graph(self.g2, "BENDMIN_TEST_3")

    def test4_problematic_graph(self):
        opt = optimization.LayeredOptimizer(self.g3)
        improved_sifting(opt.g)
        vis.draw_graph(self.g3, "BENDMIN_TEST_PROBa")
        opt.optimize_layout(bend_minimization=True, fix_x_vars=True, streamline=False)
        vis.draw_graph(self.g3, "BENDMIN_TEST_PROB")
