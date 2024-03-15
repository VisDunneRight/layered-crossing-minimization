from src.optimization import LayeredOptimizer
from src.neighborhood import degree_ratio_neighborhood, degree_candidate
from src.heuristics import barycenter
from src.vis import draw_graph


if __name__ == '__main__':
    c1, c2 = 5, 1
    ctrlflow_files = ["tail", "ptx", "sort", "sha256sum"]
    file_selection = 3
    control_flow_file = f"control-flow-graphs/{ctrlflow_files[file_selection]}/dbg.main.dot"

    optim = LayeredOptimizer(control_flow_file, stratisfimal_y_vars=True, bendiness_reduction=True, sequential_bendiness=False, gamma_1=c1, gamma_2=c2)
    # barycenter(optim.g)
    print(optim.g.n_layers)
    print(optim.g.n_nodes)
    draw_graph(optim.g, "big_test")
    # optim.local_opt_increment(1000, degree_ratio_neighborhood, degree_candidate)
