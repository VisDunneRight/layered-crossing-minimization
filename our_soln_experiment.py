from src.optimization import LayeredOptimizer
import experiments


def run_my_layout_algorithm(graph_file):
    results = []
    optimizer1 = LayeredOptimizer(graph_file)
    print(len(optimizer1.g.nodes), len(optimizer1.g.edges))
    optimizer1.direct_transitivity = True
    optimizer1.symmetry_breaking = True
    optimizer1.xvar_branch_priority = True
    optimizer1.mirror_vars = True
    optimizer1.return_experiment_data = True
    optimizer1.cutoff_time = 300
    results.append(optimizer1.optimize_layout())

    optimizer2 = LayeredOptimizer(graph_file)
    optimizer2.direct_transitivity = True
    optimizer2.symmetry_breaking = True
    optimizer2.return_experiment_data = True
    optimizer2.cutoff_time = 300
    results.append(optimizer2.optimize_layout())

    optimizer3 = LayeredOptimizer(graph_file)
    optimizer3.direct_transitivity = True
    optimizer3.symmetry_breaking = True
    optimizer3.cycle_constraints = True
    optimizer3.collapse_leaves = True
    optimizer3.return_experiment_data = True
    optimizer3.cutoff_time = 300
    results.append(optimizer3.optimize_layout())

    optimizer4 = LayeredOptimizer(graph_file)
    optimizer4.vertical_transitivity = True
    optimizer4.symmetry_breaking = True
    optimizer4.butterfly_reduction = True
    optimizer4.mirror_vars = True
    optimizer4.xvar_branch_priority = True
    optimizer4.return_experiment_data = True
    optimizer4.cutoff_time = 300
    results.append(optimizer4.optimize_layout())

    optimizer5 = LayeredOptimizer(graph_file)
    optimizer5.vertical_transitivity = True
    optimizer5.symmetry_breaking = True
    optimizer5.mirror_vars = True
    optimizer5.cycle_constraints = True
    optimizer5.collapse_leaves = True
    optimizer5.xvar_branch_priority = True
    optimizer5.return_experiment_data = True
    optimizer5.cutoff_time = 300
    results.append(optimizer5.optimize_layout())

    optimizer = LayeredOptimizer(graph_file)
    optimizer.direct_transitivity = True
    optimizer.symmetry_breaking = True
    optimizer.xvar_branch_priority = True
    optimizer.butterfly_reduction = True
    optimizer.cycle_constraints = True
    optimizer.heuristic_start = True
    optimizer.return_experiment_data = True
    optimizer.cutoff_time = 300
    results.append(optimizer.optimize_layout())

    return results


def run_best_layout(graph_file):
    optimizer = LayeredOptimizer(graph_file)
    optimizer.direct_transitivity = True
    optimizer.symmetry_breaking = True
    optimizer.xvar_branch_priority = True
    optimizer.butterfly_reduction = True
    optimizer.cycle_constraints = True
    optimizer.heuristic_start = True
    optimizer.return_experiment_data = True
    optimizer.cutoff_time = 300
    return optimizer.optimize_layout()


def run_stratisfimal_layout(graph_file):
    optimizer = LayeredOptimizer(graph_file)
    optimizer.vertical_transitivity = True
    # optimizer.direct_transitivity = True
    optimizer.return_experiment_data = True
    return optimizer.optimize_layout()


def run_optimal_sankey_layout(graph_file):
    optimizer = LayeredOptimizer(graph_file)
    optimizer.direct_transitivity = True
    optimizer.mirror_vars = True
    optimizer.butterfly_reduction = True
    optimizer.xvar_branch_priority = True
    optimizer.return_experiment_data = True
    return optimizer.optimize_layout()


def run_junger_polyhedral_layout(graph_file):
    optimizer = LayeredOptimizer(graph_file)
    optimizer.direct_transitivity = True
    optimizer.polyhedral_constraints = True
    optimizer.return_experiment_data = True
    return optimizer.optimize_layout()


def run_gange_planarity_approach_layout(graph_file):
    optimizer = LayeredOptimizer(graph_file)
    optimizer.direct_transitivity = True
    optimizer.symmetry_breaking = True
    optimizer.cycle_constraints = True
    optimizer.collapse_leaves = True
    optimizer.return_experiment_data = True
    return optimizer.optimize_layout()


def experiment(st_idx):
    with open("data storage/all_g_sorted.txt", 'r') as fd:
        graphs_files = []
        ct = 0
        for line in fd.readlines():
            if line[0] != "T":
                if ct >= st_idx:
                    graphs_files.append(line[:line.index(',')])
                ct += 1
    if st_idx == 0:
        experiments.insert_one(f"our_best_combination_10min.csv",
               ["Index", "File", "X-vars", "C-vars", "Total vars",
                "Total constraints", "Crossings", "Opttime", "Work", "Nodes visited", "Setup Time"])
    i = st_idx
    for file in graphs_files:
        print(f"Optimizing {file} ({i+1}/{ct})")
        res = run_my_layout_algorithm(file)
        formatted = [i, file] + [j for j in res]
        experiments.insert_one(f"our_best_combination_10min.csv", formatted)
        i += 1


if __name__ == '__main__':
    gname = "control-flow-graphs/chmod/dbg.main.dot"
    for i in range(5):
        xr = run_my_layout_algorithm(gname)
        print(min([x[5] for x in xr]))
        print(run_stratisfimal_layout(gname))
        print(run_junger_polyhedral_layout(gname))
        print(run_optimal_sankey_layout(gname))
        print(run_gange_planarity_approach_layout(gname))

