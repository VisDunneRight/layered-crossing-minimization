from src.optimization import LayeredOptimizer
from src.neighborhood import degree_ratio_neighborhood, degree_candidate
from src.heuristics import barycenter, global_sifting, weighted_median
from src.vis import draw_graph
from src.read_data import read


if __name__ == '__main__':
    c1, c2 = 5, 1
    ctrlflow_files = ["tail", "ptx", "sort", "sha256sum", "sleep", "ls"]
    file_selection = -1
    control_flow_file = f"control-flow-graphs/{ctrlflow_files[file_selection]}/dbg.main.dot"
    # gr = read(control_flow_file)
    gr = read("casestudy_gstore.lgbin")

    # barycenter(gr)
    # global_sifting(gr)
    # weighted_median(gr)
    print(gr.calculate_stratisfimal_objective(c1, c2))

    # optim = LayeredOptimizer("casestudy_gstore.lgbin", vertical_transitivity=True, bendiness_reduction=True, sequential_bendiness=False, gamma_1=c1, gamma_2=c2, cutoff_time=600)
    # optim = LayeredOptimizer("casestudy_gstore_2.lgbin", vertical_transitivity=True, bendiness_reduction=True, sequential_bendiness=False, gamma_1=c1, gamma_2=c2, cutoff_time=600)
    # print(optim.g.calculate_stratisfimal_objective(c1, c2))
    # optim = LayeredOptimizer(gr, vertical_transitivity=True, bendiness_reduction=True, sequential_bendiness=False, gamma_1=c1, gamma_2=c2, cutoff_time=600)
    # barycenter(optim.g)
    # print(optim.g.n_layers)
    # print(optim.g.n_nodes, len(optim.g.edges))
    # optim.local_opt_increment(800, neighborhood_fn=degree_ratio_neighborhood)
    # optim.local_opt_increment(1300, neighborhood_fn=degree_ratio_neighborhood)
    # optim.g.write_out("casestudy_gstore_2.lgbin")

    draw_graph(gr, "big_test")
    # optim.local_opt_increment(1000, degree_ratio_neighborhood, degree_candidate)
