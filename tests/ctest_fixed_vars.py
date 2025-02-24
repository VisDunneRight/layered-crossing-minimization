import random
import unittest
from src import read_data, optimization, vis


class TestOptimizationWithFixedNodes(unittest.TestCase):
    def setUp(self) -> None:
        self.g1 = read_data.read("../Rome-Lib/graficon71nodi/grafo6545.71")

    def test_fixnodes_topmiddlebottom(self):
        fix_nodes = random.sample(list(self.g1.node_ids.keys()), 20)
        opt = optimization.LayeredOptimizer(self.g1)
        self.g1.add_node_fix(fix_nodes, "top")
        opt.optimize_layout(crossing_minimization=True, fix_nodes=True)
        vis.draw_graph(self.g1, "FIX_TEST_1a", groups={nd: 1 if nd in fix_nodes else 0 for nd in self.g1.node_ids})

        self.g1.add_node_fix(fix_nodes, "middle")
        opt.optimize_layout(crossing_minimization=True, fix_nodes=True)
        vis.draw_graph(self.g1, "FIX_TEST_1b", groups={nd: 1 if nd in fix_nodes else 0 for nd in self.g1.node_ids})

        self.g1.add_node_fix(fix_nodes, "bottom")
        opt.optimize_layout(crossing_minimization=True, fix_nodes=True)
        vis.draw_graph(self.g1, "FIX_TEST_1c", groups={nd: 1 if nd in fix_nodes else 0 for nd in self.g1.node_ids})
