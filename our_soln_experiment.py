from src import optimization
import experiments


def run_my_layout_algorithm(graph_file):
    optimizer = optimization.LayeredOptimizer(graph_file)
    optimizer.vertical_transitivity = True
    optimizer.fix_one_var = True
    optimizer.mip_relax = True
    optimizer.xvar_branch_priority = True
    optimizer.aggro_presolve = True
    optimizer.mirror_vars = True
    optimizer.return_experiment_data = True
    optimizer.cutoff_time = 600
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
    experiment(7118)
