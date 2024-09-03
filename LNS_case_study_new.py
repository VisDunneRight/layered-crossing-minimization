import csv
import time

from src.optimization import LayeredOptimizer
from src.neighborhood import degree_ratio_neighborhood, degree_candidate, random_candidate
from src.heuristics import barycenter, weighted_median
from src.tabu import tabu
from src.vis import draw_graph
from src.read_data import read


def write_to_csv(fpath, values):
    with open(fpath, 'a', newline='') as csvfd:
        wrt = csv.writer(csvfd)
        wrt.writerow(values)


if __name__ == '__main__':
    c1, c2 = 3, 1
    ctrlflow_files = ["head", "csplit", "ls"]

    # Crossings Only
    for control_flow_graph in ctrlflow_files:
        flname = f"control-flow-graphs/{control_flow_graph}/dbg.main.dot"

        # 1. LNS method
        gr1 = read(flname)
        barycenter(gr1)
        optim = LayeredOptimizer(gr1, cutoff_time=900)
        _, _, crossings, times = optim.local_opt_increment(1000, neighborhood_fn=degree_ratio_neighborhood, candidate_fn=random_candidate)
        optim.g.write_out(f"LNS_case_study/LNS_cr_{control_flow_graph}.lgbin")
        write_to_csv("LNS_case_study/results.csv", ["CR", control_flow_graph, "LNS", optim.g.num_edge_crossings(), crossings, times])

        # 2. Tabu
        gr2 = read(flname)
        crossings, times = tabu(gr2, timelimit=900, improvement_data=True)
        gr2.write_out(f"LNS_case_study/TABU_cr_{control_flow_graph}.lgbin")
        write_to_csv("LNS_case_study/results.csv", ["CR", control_flow_graph, "TABU", gr2.num_edge_crossings(), crossings, times])

        # 3. ILP
        gr3 = read(flname)
        optim = LayeredOptimizer(gr3, symmetry_breaking=True, cutoff_time=900, record_solution_data_over_time=True)
        crossings, times = optim.optimize_layout()
        optim.g.write_out(f"LNS_case_study/ILP_cr_{control_flow_graph}.lgbin")
        write_to_csv("LNS_case_study/results.csv", ["CR", control_flow_graph, "ILP", optim.g.num_edge_crossings(), crossings, times])

    # Crossings + Edge Bends
    for control_flow_graph in ctrlflow_files:
        flname = f"control-flow-graphs/{control_flow_graph}/dbg.main.dot"

        # 1. LNS method
        gr1 = read(flname)
        optim = LayeredOptimizer(gr1, vertical_transitivity=True, bendiness_reduction=True, sequential_bendiness=False, gamma_1=c1, gamma_2=c2, cutoff_time=900)
        _, _, crossings, times = optim.local_opt_increment(1000, neighborhood_fn=degree_ratio_neighborhood, candidate_fn=random_candidate)
        optim.g.write_out(f"LNS_case_study/LNS_cr+br_{control_flow_graph}.lgbin")
        write_to_csv("LNS_case_study/results.csv", ["CR+BR", control_flow_graph, "LNS", optim.g.calculate_stratisfimal_objective(c1, c2), crossings, times])

        # 2. Weighted Median
        gr2 = read(flname)
        gr2t1 = time.time()
        init_obj = gr2.calculate_stratisfimal_objective(c1, c2)
        weighted_median(gr2)
        gr2.write_out(f"LNS_case_study/WM_cr+br_{control_flow_graph}.lgbin")
        write_to_csv("LNS_case_study/results.csv", ["CR+BR", control_flow_graph, "WM", gr2.calculate_stratisfimal_objective(c1, c2), [init_obj, gr2.calculate_stratisfimal_objective(c1, c2)], [0, time.time() - gr2t1]])

        # 3. ILP
        gr3 = read(flname)
        optim2 = LayeredOptimizer(gr3, vertical_transitivity=True, bendiness_reduction=True, sequential_bendiness=False, symmetry_breaking=True, gamma_1=c1, gamma_2=c2, cutoff_time=900, record_solution_data_over_time=True)
        crossings, times = optim2.optimize_layout()
        optim2.g.write_out(f"LNS_case_study/ILP_cr+br_{control_flow_graph}.lgbin")
        write_to_csv("LNS_case_study/results.csv", ["CR+BR", control_flow_graph, "ILP", optim2.g.calculate_stratisfimal_objective(c1, c2), crossings, times])
