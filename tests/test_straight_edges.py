import unittest
from src import read_data, optimization, vis


class TestOptimizationWithStraightLongArcs(unittest.TestCase):
    def setUp(self) -> None:
        self.gsmall = read_data.read("../Rome-Lib/graficon41nodi/grafo2019.41")
        self.gsmall2 = read_data.read("../Rome-Lib/graficon41nodi/grafo2019.41")
        self.glarge = read_data.read("../Rome-Lib/graficon71nodi/grafo6545.71")

    def test_straightedge_basic(self):
        opt = optimization.LayeredOptimizer(self.gsmall)
        opt.optimize_layout(crossing_minimization=True, constrain_straight_long_arcs=True, m_val=30)
        vis.draw_graph(self.gsmall, "STRAIGHT_TEST")

        opt = optimization.LayeredOptimizer(self.gsmall2)
        opt.optimize_layout(crossing_minimization=True, constrain_straight_long_arcs=True, vertical_transitivity=True)
        vis.draw_graph(self.gsmall2, "STRAIGHT_TEST_2")

    def test_straightedge_large(self):
        opt = optimization.LayeredOptimizer(self.glarge)
        opt.optimize_layout(crossing_minimization=True, constrain_straight_long_arcs=True, vertical_transitivity=True)
        vis.draw_graph(self.glarge, "STRAIGHT_TEST_3")

    def test_straightedge_2bends(self):
        opt = optimization.LayeredOptimizer(self.gsmall2)
        opt.optimize_layout(crossing_minimization=True, constrain_straight_long_arcs=True, vertical_transitivity=True, long_arc_bend_limit=1)
        opt.optimize_layout(constrain_straight_long_arcs=True, long_arc_bend_limit=1, fix_x_vars=True, bendiness_reduction=True)
        vis.draw_graph(self.gsmall2, "STRAIGHT_TEST_4")
