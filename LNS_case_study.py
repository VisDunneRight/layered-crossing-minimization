import csv
from src.optimization import LayeredOptimizer
from src.neighborhood import degree_ratio_neighborhood, degree_candidate
from src.heuristics import barycenter, weighted_median
from src.vis import draw_graph
from src.read_data import read


if __name__ == '__main__':
    # Note: this case study takes exactly 1 hour to run
    c1, c2 = 3, 1
    # ctrlflow_files = ["tail", "ptx", "sort", "sha256sum", "sleep", "ls"]
    # file_selection = -1
    control_flow_file = f"control-flow-graphs/ls/dbg.main.dot"
    gr = read(control_flow_file)
    # gr = read("casestudy_gstorenew_1.lgbin")
    # gr = read("casestudy_gstorenew_2.lgbin")
    # gr = read("casestudy_gstorenew_3.lgbin")

    barycenter(gr)
    # weighted_median(gr)
    print(gr.n_nodes, len(gr.edges))
    print(gr.calculate_stratisfimal_objective(c1, c2))

    optim = LayeredOptimizer(gr, vertical_transitivity=True, bendiness_reduction=True, sequential_bendiness=False, gamma_1=c1, gamma_2=c2, cutoff_time=600)
    optim.just_bendiness_reduction(streamline=False)

    # small neighborhoods
    res = optim.local_opt_increment(800, neighborhood_fn=degree_ratio_neighborhood, candidate_fn=degree_candidate)

    optim.g.write_out("casestudy_gstorenew_1.lgbin")
    with open("casestudy_data.csv", 'a') as fd:
        wrt = csv.writer(fd)
        wrt.writerow(res)
    optim.cutoff_time = 600

    # medium neighborhoods
    res = optim.local_opt_increment(1100, neighborhood_fn=degree_ratio_neighborhood, candidate_fn=degree_candidate)

    optim.g.write_out("casestudy_gstorenew_2.lgbin")
    with open("casestudy_data.csv", 'a') as fd:
        wrt = csv.writer(fd)
        wrt.writerow(res)
    optim.cutoff_time = 600

    # large neighborhoods
    res = optim.local_opt_increment(1400, neighborhood_fn=degree_ratio_neighborhood, candidate_fn=degree_candidate)

    optim.g.write_out("casestudy_gstorenew_3.lgbin")
    with open("casestudy_data.csv", 'a') as fd:
        wrt = csv.writer(fd)
        wrt.writerow(res)

    # global optimal
    gr2 = read(control_flow_file)
    optim2 = LayeredOptimizer(gr2, vertical_transitivity=True, bendiness_reduction=True, sequential_bendiness=False, gamma_1=c1, gamma_2=c2, cutoff_time=1800)
    res = optim2.optimize_layout()
    with open("casestudy_data.csv", 'a') as fd:
        wrt = csv.writer(fd)
        wrt.writerow(res)

    print(gr.calculate_stratisfimal_objective(c1, c2))

    draw_graph(gr, "case_study_ls_graph")
