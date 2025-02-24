import random
import unittest
from src import read_data, optimization, vis


class TestOptimizationWithNodeEmphasis(unittest.TestCase):
    def setUp(self) -> None:
        self.g1 = read_data.read("../Rome-Lib/graficon41nodi/grafo2019.41")
        self.g1.add_node_weights({9: 6})
        self.g1.add_edge_weights({(2, 17): 5, (2, 18): 5, (2, 0): 5, (25, 9): 4, (25, 10): 4})
        self.g2 = read_data.read("../Rome-Lib/graficon35nodi/grafo179.35")
        self.g2.add_node_weights([random.randint(1, 10) for _ in range(self.g2.n_nodes)])

    def test_node_weight_spacing_basic(self):
        opt = optimization.LayeredOptimizer(self.g1)
        opt.optimize_layout(crossing_minimization=True, apply_node_weight_spacing=True)
        vis.draw_graph(self.g1, "NDWT_TEST")

    def test_both_weights(self):
        opt = optimization.LayeredOptimizer(self.g1)
        opt.optimize_layout(bendiness_reduction=True, apply_node_weight_spacing=True, apply_edge_weight=True, hybrid_constraints=[("crossings", 4)])
        vis.draw_graph(self.g1, "EDWT_TEST")

    def test_lots_of_node_weights(self):
        opt = optimization.LayeredOptimizer(self.g2)
        opt.m_val *= 3
        opt.optimize_layout(crossing_minimization=True, apply_node_weight_spacing=True)
        vis.draw_graph(self.g2, "NDWT_TEST_2")
