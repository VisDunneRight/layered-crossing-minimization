import unittest
from layered_optimization import read_data, optimization, vis, metrics
from layered_optimization.heuristics import improved_sifting
import os


class TestVariableFixing(unittest.TestCase):
    def setUp(self) -> None:
        self.g1 = read_data.read("../datasets/Rome-Lib/graficon18nodi/grafo198.18")
        self.g2 = read_data.read("../datasets/Rome-Lib/graficon32nodi/grafo1053.32")
        self.g3 = read_data.read("../datasets/Rome-Lib/graficon10nodi/grafo166.10")
        self.g4 = read_data.read("../datasets/random_graphs/uniform_layered/n5/graph144.lgbin")
        self.g5 = read_data.read("../datasets/random_graphs/uniform_layered/n40/graph118.lgbin")

    def test0_compare_fixing_x_with_c(self):
        x_times, c_times = [], []
        x_cr_obj, c_cr_obj, x_b_obj, c_b_obj = [], [], [], []
        n_graphs = 10
        for fl in os.listdir("../datasets/random_graphs/uniform_layered/n60")[:n_graphs]:
            gr = read_data.read(f"../datasets/random_graphs/uniform_layered/n60/{fl}")
            opt = optimization.LayeredOptimizer(gr)
            res0 = opt.optimize_layout(crossing_minimization=True)
            res1 = opt.optimize_layout(edge_length_minimization=True, fix_x_vars=True)
            x_times.append(res1.runtime)
            x_cr_obj.append(res0.metrics.crossings)
            x_b_obj.append(res1.metrics.length)
            res2 = opt.optimize_layout(edge_length_minimization=True, fix_c_vars=True)
            c_times.append(res2.runtime)
            c_cr_obj.append(res2.metrics.crossings)
            c_b_obj.append(res2.metrics.length)
        print(f"X-FIXING AVG TIME:{sum(x_times)/n_graphs} / OBJ:{(sum(x_cr_obj) + sum(x_b_obj))/n_graphs}")
        print(f"C-FIXING AVG TIME:{sum(c_times)/n_graphs} / OBJ:{(sum(c_cr_obj) + sum(c_b_obj))/n_graphs}")
        self.assertEqual(sum(c_cr_obj), sum(x_cr_obj))
        self.assertLess(sum(c_b_obj), sum(x_b_obj))

    def test_05_c_fixing_metrics(self):
        opt = optimization.LayeredOptimizer(self.g5)
        res1 = opt.optimize_layout(crossing_minimization=True)
        res2 = opt.optimize_layout(fix_c_vars=True, edge_length_minimization=True)
        self.assertEqual(res1.metrics.crossings, res2.metrics.crossings)

        res1 = opt.optimize_layout(planarization=True)
        res2 = opt.optimize_layout(fix_c_vars=True, edge_length_minimization=True, planarization=True)
        self.assertEqual(res1.metrics.planar, res2.metrics.planar)

        res1 = opt.optimize_layout(min_max_crossings=True)
        res2 = opt.optimize_layout(fix_c_vars=True, edge_length_minimization=True)
        self.assertEqual(res1.metrics.minmax, metrics.min_max_crossings(self.g5))

    def test1_fix_symmetry_vars(self):
        opt = optimization.LayeredOptimizer(self.g1)
        res1 = opt.optimize_layout(horizontal_node_symmetry=True)
        res2 = opt.optimize_layout(fix_horizontal_node_symmetry=True, vertical_node_symmetry=True)
        self.assertEqual(res1.metrics.horizontal_node_symmetry, res2.metrics.horizontal_node_symmetry)

        res1 = opt.optimize_layout(vertical_node_symmetry=True)
        res2 = opt.optimize_layout(fix_vertical_node_symmetry=True, edge_length_minimization=True)
        self.assertEqual(res1.metrics.vertical_node_symmetry, res2.metrics.vertical_node_symmetry)

        res1 = opt.optimize_layout(rotation_node_symmetry=True)
        res2 = opt.optimize_layout(fix_rotation_node_symmetry=True, edge_length_minimization=True)
        self.assertEqual(res1.metrics.rotation_node_symmetry, res2.metrics.rotation_node_symmetry)

        res1 = opt.optimize_layout(horizontal_edge_symmetry=True)
        res2 = opt.optimize_layout(fix_horizontal_edge_symmetry=True, vertical_edge_symmetry=True)
        self.assertEqual(res1.metrics.horizontal_edge_symmetry, res2.metrics.horizontal_edge_symmetry)

        res1 = opt.optimize_layout(vertical_edge_symmetry=True)
        res2 = opt.optimize_layout(fix_vertical_edge_symmetry=True, edge_length_minimization=True)
        self.assertEqual(res1.metrics.vertical_edge_symmetry, res2.metrics.vertical_edge_symmetry)

        res1 = opt.optimize_layout(rotation_edge_symmetry=True)
        res2 = opt.optimize_layout(fix_rotation_edge_symmetry=True, edge_length_minimization=True)
        self.assertEqual(res1.metrics.rotation_edge_symmetry, res2.metrics.rotation_edge_symmetry)

    def test2_fix_angle(self):
        opt = optimization.LayeredOptimizer(self.g3)
        res1 = opt.optimize_layout(crossing_angle=True)
        res2 = opt.optimize_layout(fix_crossing_angle=True, edge_length_minimization=True)
        vis.draw_graph(self.g3, "FIX_TEST_2_ANG")
        self.assertAlmostEqual(res1.metrics.angle, metrics.crossing_angle(self.g3))

    def test3_fix_bends(self):
        opt = optimization.LayeredOptimizer(self.g5)
        res1 = opt.optimize_layout(crossing_minimization=True)
        res2 = opt.optimize_layout(fix_c_vars=True, bend_minimization=True)
        res3 = opt.optimize_layout(fix_c_vars=True, fix_edge_bends=True, edge_length_minimization=True)
        vis.draw_graph(self.g5, "FIX_TEST_3_BEND")
        self.assertEqual(res1.metrics.crossings, res2.metrics.crossings, res3.metrics.crossings)
        self.assertEqual(res2.metrics.bend, metrics.bend_minimization(self.g5))

    def test4_fix_planar_edge_set(self):
        opt = optimization.LayeredOptimizer(self.g5)
        res00 = opt.optimize_layout(crossing_minimization=True)
        res1 = opt.optimize_layout(planarization=True)
        print("planar removal edges:", [k for k,v in opt.fix_planarization_removals.items() if v == 1])
        self.assertNotEqual(res00.metrics.crossings, res1.metrics.crossings)
        res2 = opt.optimize_layout(fix_planarization=True, crossing_minimization=True)
        res3 = opt.optimize_layout(fix_c_vars=True, fix_planarization=True, edge_length_minimization=True)
        vis.draw_graph(self.g5, "FIX_TEST_4_PLANAR")
        self.assertEqual(res1.metrics.planar, res2.metrics.planar, res3.metrics.planar)

    def test5_fix_min_max(self):
        opt = optimization.LayeredOptimizer(self.g5)
        res1 = opt.optimize_layout(min_max_crossings=True)
        print("counts per edge:", metrics.min_max_crossings(self.g5, return_counts_per_edge=True))
        vis.draw_graph(self.g5, "FIX_TEST_5a_MINMAX")
        res2 = opt.optimize_layout(fix_min_max_crossings=True, crossing_minimization=True)
        res3 = opt.optimize_layout(fix_c_vars=True, fix_min_max_crossings=True, edge_length_minimization=True)
        vis.draw_graph(self.g5, "FIX_TEST_5_MINMAX")
        self.assertEqual(res1.metrics.minmax, res2.metrics.minmax, res3.metrics.minmax)
