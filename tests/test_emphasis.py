import unittest
from src import read_data, optimization, vis


class TestOptimizationWithNodeEmphasis(unittest.TestCase):
    def setUp(self) -> None:
        self.g1 = read_data.read("../Rome-Lib/graficon41nodi/grafo2019.41")
        self.g1.add_node_emphasis([2])

    def test_emphasis_basic(self):
        opt = optimization.LayeredOptimizer(self.g1)
        res1 = opt.optimize_layout(crossing_minimization=True, node_emphasis=True)
        vis.draw_graph(self.g1, "EMPH_TEST")

        opt2 = optimization.LayeredOptimizer(self.g1)
        res2 = opt2.optimize_layout(crossing_minimization=True)
        vis.draw_graph(self.g1, "EMPH_TEST_compr")

        self.assertGreater(res1.objval, res2.objval)
