import unittest
from layered_optimization import read_data, optimization, vis, metrics
from layered_optimization.heuristics import improved_sifting
import os


class TestAngularOptimization(unittest.TestCase):
    def setUp(self) -> None:
        self.g1 = read_data.read("../datasets/Rome-Lib/graficon18nodi/grafo198.18")
        self.g2 = read_data.read("../datasets/Rome-Lib/graficon32nodi/grafo1053.32")
        self.g3 = read_data.read("../datasets/Rome-Lib/graficon10nodi/grafo166.10")
        self.g4 = read_data.read("../datasets/random_graphs/uniform_layered/n5/graph144.lgbin")

    def test1_angle_opt(self):
        opt = optimization.LayeredOptimizer(self.g1)
        opt.optimize_layout(crossing_minimization=True)
        res = opt.optimize_layout(crossing_angle=True, fix_x_vars=True, m_val=12)
        vis.draw_graph(self.g1, "ANG_TEST_CR_THEN_ANG")
        self.assertAlmostEqual(res.objval, metrics.crossing_angle(self.g1, use_approximation=True))

    def test2_angle_opt_bland(self):
        opt = optimization.LayeredOptimizer(self.g3)
        opt.optimize_layout(crossing_angle=True)
        vis.draw_graph(self.g3, "ANG_TEST_SOLO")

    def test3_angle_opt_bigger(self):
        opt = optimization.LayeredOptimizer(self.g2)
        opt.optimize_layout(crossing_minimization=True)
        opt.optimize_layout(crossing_angle=True, fix_x_vars=True)
        vis.draw_graph(self.g2, "ANG_TEST_CR_THEN_ANG_BIG")

    def test4_angle_combined(self):
        opt = optimization.LayeredOptimizer(self.g1)
        opt.optimize_layout(crossing_minimization=True, crossing_angle=True, m_val=10)
        vis.draw_graph(self.g1, "ANG_TEST_PLUS_CR")

    def test5_angle_combined_with_bend(self):
        opt = optimization.LayeredOptimizer(self.g1)
        opt.optimize_layout(crossing_minimization=True)
        opt.optimize_layout(crossing_angle=True, edge_length_minimization=True, gamma_edgelength=0.1, fix_x_vars=True, m_val=12)
        vis.draw_graph(self.g1, "ANG_TEST_W_BEND")

    def test6_fixed_x_sifting(self):
        opt = optimization.LayeredOptimizer(self.g2)
        improved_sifting(opt.g)
        vis.draw_graph(self.g2, "ANG_TEST_5a")
        opt.optimize_layout(fix_x_vars=True, crossing_angle=True)
        vis.draw_graph(self.g2, "ANG_TEST_5")

    def test7_angle_opt_hybrid(self):
        opt = optimization.LayeredOptimizer(self.g4)
        res0 = opt.optimize_layout(crossing_angle=True)
        print("RES0 OBJVAL", res0.objval)
        opt.optimize_layout(crossing_minimization=True, hybrid_constraints=[("crossing_angle", res0.objval + 0.001)])
        vis.draw_graph(self.g4, "ANG_TEST_HYBRD")

    def test8_true_angle_with_tangent(self):
        opt = optimization.LayeredOptimizer(self.g4)
        res0 = opt.optimize_layout(crossing_minimization=True)
        res = opt.optimize_layout(crossing_angle=True, ca_use_tangent=True, fix_x_vars=True)
        self.assertAlmostEqual(res.objval, metrics.crossing_angle(self.g4))
        vis.draw_graph(self.g4, "ANG_TEST_TAN")

    def test9_stress_test_compare(self):
        tru_times, apr_times = [], []
        tru_obj, apr_obj = [], []
        n1 = 10
        for i, fl in enumerate(os.listdir("../datasets/random_graphs/uniform_layered/n10")[:n1]):
            gr = read_data.read(f"../datasets/random_graphs/uniform_layered/n10/{fl}")
            opt = optimization.LayeredOptimizer(gr)
            res = opt.optimize_layout(crossing_angle=True, ca_use_tangent=True, fix_x_vars=False)
            vis.draw_graph(gr, f"ANG_TAN/gr{i}_TAN")
            tru_times.append(res.runtime)
            tru_obj.append(res.objval)
            res2 = opt.optimize_layout(crossing_angle=True, fix_x_vars=False)
            vis.draw_graph(gr, f"ANG_TAN/gr{i}_APX")
            apr_times.append(res2.runtime)
            apr_obj.append(res2.objval)
        n2 = 10
        for i, fl in enumerate(os.listdir("../datasets/random_graphs/uniform_layered/n30")[:n2]):
            gr = read_data.read(f"../datasets/random_graphs/uniform_layered/n30/{fl}")
            opt = optimization.LayeredOptimizer(gr)
            res0 = opt.optimize_layout(crossing_minimization=True)
            res = opt.optimize_layout(crossing_angle=True, ca_use_tangent=True, fix_x_vars=True)
            vis.draw_graph(gr, f"ANG_TAN/gr{i + n1}_TAN")
            tru_times.append(res.runtime)
            tru_obj.append(res.objval)
            res2 = opt.optimize_layout(crossing_angle=True, fix_x_vars=True)
            vis.draw_graph(gr, f"ANG_TAN/gr{i + n1}_APX")
            apr_times.append(res2.runtime)
            apr_obj.append(res2.objval)
        print("TRU TIME:", sum(tru_times)/(n1 + n2), tru_times)
        print("APRX TIME:", sum(apr_times)/(n1 + n2), apr_times)
        print("TRU OBJ:", sum(tru_obj)/(n1 + n2), tru_obj)
        print("APRX OBJ:", sum(apr_obj)/(n1 + n2), apr_obj)

    def test10_angle_with_minimum(self):
        opt = optimization.LayeredOptimizer(self.g2)
        res0 = opt.optimize_layout(crossing_minimization=True)
        res = opt.optimize_layout(crossing_angle=True, ca_use_minimum=True, edge_length_minimization=True, fix_x_vars=True)
        self.assertAlmostEqual(res.metrics.angle, metrics.crossing_angle(self.g2, use_min_approx=True), places=4)
        vis.draw_graph(self.g2, "ANG_TEST_MIN")
