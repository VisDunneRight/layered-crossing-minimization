import os
import unittest
from src import read_data, optimization, vis


class TestSameLayerEdgeOptimization(unittest.TestCase):
    def setUp(self) -> None:
        self.g1 = read_data.read("../Rome-Lib/graficon71nodi/grafo6545.71", remove_sl=False)
        self.g2 = read_data.read("../Rome-Lib/graficon71nodi/grafo6545.71", remove_sl=False)
        self.g2.add_edge(63, 64)
        self.g2.add_edge(8, 64)
        self.g2.add_edge(8, 34)

    def test_errorcheck_sl_edges(self):
        for grpath in os.listdir("../Rome-Lib/graficon71nodi"):
            the_g = read_data.read("../Rome-Lib/graficon71nodi/" + grpath, remove_sl=False)
            if not all((e1.n1.layer != e1.n2.layer for e1 in the_g.edges)):
                print(f"{grpath} Has SL Edge")
                opt = optimization.LayeredOptimizer(the_g, symmetry_breaking=True)
                opt.optimize_layout()
                break
            else:
                print("No SL Edge")

    def test_errorcheck_sl_edges_w_drawing(self):
        print("No SL Edge" if all((not e1.same_layer_edge for e1 in self.g1.edges)) else "Has SL Edge")
        opt = optimization.LayeredOptimizer(self.g1, symmetry_breaking=True)
        opt.optimize_layout()
        vis.draw_graph(self.g1, "SL_TEST")

    def test_errorcheck_double_sl_edges(self):
        print(list(filter(lambda x: x.same_layer_edge, self.g2.edges)))
        opt = optimization.LayeredOptimizer(self.g2, symmetry_breaking=True)
        opt.optimize_layout()
        vis.draw_graph(self.g2, "SL_TEST_2")

