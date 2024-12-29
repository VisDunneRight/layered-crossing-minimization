import os
import random
import unittest
from src import read_data, optimization, vis


class TestOptimizationWithGroups(unittest.TestCase):
    def setUp(self) -> None:
        self.g1 = read_data.read("../Rome-Lib/graficon71nodi/grafo6545.71")
        self.g1.add_groups([[30, 37, 1, 25], [57, 43], [52, 9, 24]])
        self.g2 = read_data.read("../Rome-Lib/graficon71nodi/grafo6545.71")
        self.g2.add_groups([[31, 38, 18, 27], [2, 65, 12, 28]])
        self.g3 = read_data.read("../Rome-Lib/graficon71nodi/grafo6545.71")
        self.g3.add_groups([[30, 37, 1, 25], [52, 9, 24], [31, 38, 18, 27], [2, 65, 12, 28], [57, 43, 19], [41, 60]])
        self.g4 = read_data.read("../Rome-Lib/graficon71nodi/grafo6545.71")

    def test_sl_groups(self):
        opt = optimization.LayeredOptimizer(self.g1, symmetry_breaking=True, grouping_constraints=True)
        opt.optimize_layout()
        vis.draw_graph(self.g1, "GRP_TEST")

    def test_ml_groups(self):
        opt = optimization.LayeredOptimizer(self.g2, symmetry_breaking=True, grouping_constraints=True)
        opt.optimize_layout()
        vis.draw_graph(self.g2, "GRP_TEST_2")

    def test_all_groups_hard(self):
        opt = optimization.LayeredOptimizer(self.g3, symmetry_breaking=True, grouping_constraints=True)
        opt.optimize_layout()
        vis.draw_graph(self.g3, "GRP_TEST_3")

    def test_all_groups_hard_with_bendiness(self):
        opt = optimization.LayeredOptimizer(self.g3, symmetry_breaking=True, grouping_constraints=True, sequential_bendiness=True, bendiness_reduction=True, sequential_grouping_constraints=True)
        opt.optimize_layout()
        vis.draw_graph(self.g3, "GRP_TEST_3bend")

    def test_group_correctness(self):
        opt = optimization.LayeredOptimizer(self.g4, symmetry_breaking=True)
        res1 = opt.optimize_layout()
        vis.draw_graph(self.g4, "GRP_TEST_4a")
        grouped = set()
        groups = []
        for _ in range(10):
            rn = random.choice(self.g4.nodes)
            lv, yv = rn.layer, rn.y
            if rn.id not in grouped:
                groups.append([rn.id])
                grouped.add(rn.id)
                for lid, py in [(lv, yv - 1), (lv, yv + 1), (lv - 1, yv - 1), (lv - 1, yv), (lv - 1, yv + 1), (lv + 1, yv - 1), (lv + 1, yv), (lv + 1, yv + 1)]:
                    if lid in self.g4.layers:
                        for nd in self.g4.layers[lid]:
                            if nd.y == py and nd.id not in grouped:
                                groups[-1].append(nd.id)
                                grouped.add(nd.id)
                                continue
        self.g4.add_groups(groups)
        opt2 = optimization.LayeredOptimizer(self.g4, symmetry_breaking=True, grouping_constraints=True)
        res2 = opt2.optimize_layout()
        vis.draw_graph(self.g4, "GRP_TEST_4b")
        self.assertEqual(res1[1], res2[1])
