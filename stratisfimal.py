import os
import sys
import random
import cProfile
import pstats
from pstats import SortKey

import networkx as nx

from src import vis, layering, motifs, experiments
from src.graph import *
from src.read_data import *
from src.optimization import LayeredOptimizer


def run_all_rome_lib(num_nodes, num_graphs, num_drawings, bendiness_reduction, seq_bend, timelimit, save, savefile=None, shuffle=False, target=None, subgraph_reduction=False):
    i = 1
    outputs = ["Results:\n"]
    if target is not None and int(target) != 0:
        to_optimize = [f"grafo{target}.{num_nodes}"]
    elif shuffle:
        to_optimize = random.sample(os.listdir(f"Rome-Lib/graficon{num_nodes}nodi"), num_graphs)
    else:
        to_optimize = os.listdir(f"Rome-Lib/graficon{num_nodes}nodi")[:num_graphs]
    for file in to_optimize:
        g, tvert = layering.create_better_layered_graph(f"graficon{num_nodes}nodi/{file}", 4, 2)
        # g = layering.create_layered_graph(f"graficon{num_nodes}nodi/{file}")
        print(f"\n\n{file} ({i}/{num_graphs}):")
        # print("Number of butterflies:", motifs.count_butterflies(g))
        # print("Vertex promotion time:", tvert)
        run_optimizer(g, bendiness_reduction, seq_bend, timelimit, subgraph_reduction)
        # outputs.append(f"{num_nodes},{bendiness_reduction},{seq_bend},{timelimit}," + run_optimizer(g, bendiness_reduction, seq_bend, timelimit, subgraph_reduction) + '\n')
        i += 1
        if num_drawings > 0:
            num_drawings -= 1
            vis.draw_graph(g, f"Rome-Lib/{file}")
    if save:
        with open(savefile, 'a') as f:
            f.writelines(outputs)
    else:
        for output in outputs:
            print(output.replace('\n', ''))


def run_optimizer(g: LayeredGraph, bendiness_reduction, sequential, timelimit, subgraph):
    params = {"bendiness_raduction": bendiness_reduction, "sequential_bendiness": sequential, "do_subg_reduction": subgraph}
    if len(timelimit) > 0 and int(timelimit) > 0:
        params["cutoff_time"] = int(timelimit)
    optimizer = LayeredOptimizer(g, params)
    optimizer.optimize_layout()
    # return ','.join(str(e) for e in res[0] + res[1])


def run_test_pos_1_to_n():
    n_nodes = 57
    to_optimize = os.listdir(f"Rome-Lib/graficon{n_nodes}nodi")[:10]
    for to_opt in to_optimize:
        g = layering.create_better_layered_graph(f"graficon{n_nodes}nodi/{to_opt}", 4, 2)[0]
        save_vars, t_orig = constraints.optimize_layout(g, False, store_x_vars=True, cutoff_time=60)
        if t_orig < 55:
            times = []
            order = list(range(1, len(g.nodes)+1))
            random.seed(10)
            random.shuffle(order)
            current_vars_fixed = {}
            for i in range(len(g.nodes)):
                for k, v in save_vars.items():
                    if int(k[2:k.index(',')]) == order[i] or int(k[k.index(',') + 1:k.index(']')]) == order[i]:
                        current_vars_fixed[k] = v
                times.append(constraints.optimize_layout(g, False, assignment=current_vars_fixed))
            with open("1toNexperiment.txt", 'a') as f:
                f.write(str(n_nodes) + " " + str(to_opt[5:9]) + " " + str(t_orig) + " " + " ".join([str(i) for i in times]) + "\n")


def run_test_relative_1_to_n():
    n_nodes = 59
    to_optimize = os.listdir(f"Rome-Lib/graficon{n_nodes}nodi")[:10]
    for to_opt in to_optimize:
        g = layering.create_better_layered_graph(f"graficon{n_nodes}nodi/{to_opt}", 4, 2)[0]
        optimizer = LayeredOptimizer(g, {"return_x_vars": True})
        save_vars, t_orig = optimizer.optimize_layout()
        if t_orig < 55:
            times = []
            order = list(save_vars.keys())
            random.seed(10)
            random.shuffle(order)
            current_vars_fixed = {}
            for i in range(10):
                current_vars_fixed[order[i]] = save_vars[order[i]]
                # print(order[i], current_vars_fixed)
                times.append(constraints.optimize_layout(g, False, assignment=current_vars_fixed))
            for i in range(10, len(order), 10):
                for j in range(i, min(i+10, len(order))):
                    current_vars_fixed[order[j]] = save_vars[order[j]]
                times.append(constraints.optimize_layout(g, False, assignment=current_vars_fixed))
            with open("1toNexperiment.txt", 'a') as f:
                f.write(str(n_nodes) + " " + str(to_opt[5:9]) + " " + str(t_orig) + " " + " ".join([str(i) for i in times]) + "\n")


def run_test_start_assignments():
    n_nodes = 59
    to_optimize = os.listdir(f"Rome-Lib/graficon{n_nodes}nodi")[:10]
    for to_opt in to_optimize:
        g = layering.create_better_layered_graph(f"graficon{n_nodes}nodi/{to_opt}", 4, 2)[0]
        optimizer = LayeredOptimizer(g, {"return_x_vars": True, "name": to_opt})
        t_orig, save_vars = optimizer.optimize_layout()
        if t_orig < 50:
            times = []
            order = list(save_vars.keys())
            random.seed(22)
            random.shuffle(order)
            current_vars_started = {}
            i = 0
            while i < len(order) - 30:
                for j in range(i, min(i + 20, len(order))):
                    current_vars_started[order[j]] = save_vars[order[j]]
                optimizer = LayeredOptimizer(g, {})
                times.append(optimizer.optimize_with_starting_assignments(current_vars_started))
                i += 20
            for j in range(i, len(order)):
                current_vars_started[order[j]] = save_vars[order[j]]
                optimizer = LayeredOptimizer(g, {})
                times.append(optimizer.optimize_with_starting_assignments(current_vars_started))
            with open("1toN_varstart_experiment.txt", 'a') as f:
                f.write(str(n_nodes) + " " + str(to_opt[5:9]) + " " + str(t_orig) + " " + " ".join([str(i) for i in times]) + "\n")


def run_test_start_assignments_with_misleading():
    n_nodes = 59
    to_optimize = os.listdir(f"Rome-Lib/graficon{n_nodes}nodi")[:10]
    timeprint = []
    for to_opt in to_optimize:
        g = layering.create_better_layered_graph(f"graficon{n_nodes}nodi/{to_opt}", 4, 2)[0]
        optimizer = LayeredOptimizer(g, {"name": to_opt})
        t_orig = optimizer.optimize_layout()
        print(to_opt, t_orig)
    # for to_opt in to_optimize:
    #     g = layering.create_better_layered_graph(f"graficon{n_nodes}nodi/{to_opt}", 4, 2)[0]
    #     optimizer = LayeredOptimizer(g, {"return_x_vars": True, "name": to_opt})
    #     t_orig, save_vars = optimizer.optimize_layout()
    #     optimizer.return_x_vars = False
    #     optimizer.name = f"{to_opt} full start"
    #     t2 = optimizer.optimize_with_starting_assignments(save_vars)
    #     timeprint.append(f"{to_opt}: {t_orig}, {t2}")
    for prt in timeprint:
        print(prt)


def run_test_fix_x_vars():
    n_nodes = 59
    to_optimize = os.listdir(f"Rome-Lib/graficon{n_nodes}nodi")[:10]
    for to_opt in to_optimize:
        g = layering.create_better_layered_graph(f"graficon{n_nodes}nodi/{to_opt}", 4, 2)[0]
        optimizer = LayeredOptimizer(g, {"name": to_opt})
        t_orig, optval = optimizer.optimize_layout()
        times = []
        vals = []
        for i in range(25):
            pair = optimizer.optimize_with_starting_assignments(optimizer.generate_random_vars_to_fix(2))
            times.append(pair[0])
            vals.append(pair[1])
        with open("1toN_varstart_experiment.txt", 'a') as f:
            f.write(str(n_nodes) + " " + str(to_opt[5:9]) + "\n" + str(t_orig) + " " + " ".join(
                [str(i) for i in times]) + "\n" + str(optval) + " " + " ".join([str(i) for i in vals]) + "\n")


def run_my_algorithm():
    n_nodes = 77
    to_optimize = os.listdir(f"Rome-Lib/graficon{n_nodes}nodi")[:1]
    for to_opt in to_optimize:
        g = layering.create_better_layered_graph(f"graficon{n_nodes}nodi/{to_opt}", 4, 2)[0]
        optimizer = LayeredOptimizer(g, {"name": to_opt, "verbose": True, "do_subg_reduction": True})
        optimizer.optimize_layout()


def run_standard_version():
    n_nodes = 68
    to_optimize = os.listdir(f"Rome-Lib/graficon{n_nodes}nodi")[22:23]
    for to_opt in to_optimize:
        g = layering.create_better_layered_graph(f"graficon{n_nodes}nodi/{to_opt}", 4, 2)[0]
        optimizer = LayeredOptimizer(g, {"name": to_opt})
        optimizer.optimize_layout()


def butterfly_experiment():
    n_nodes = [68, 69, 69, 68, 68, 68]
    opt68 = os.listdir(f"Rome-Lib/graficon68nodi")
    opt69 = os.listdir(f"Rome-Lib/graficon69nodi")
    to_optimize = [opt68[5], opt69[5], opt69[8], opt68[16], opt68[41], opt68[56]]
    # to_optimize = os.listdir(f"Rome-Lib/graficon{n_nodes}nodi")[50: 60]
    times1 = []
    optvals1 = []
    times2 = []
    optvals2 = []
    for i, to_opt in enumerate(to_optimize):
        g = layering.create_better_layered_graph(f"graficon{n_nodes[i]}nodi/{to_opt}", 4, 2)[0]
        optimizer = LayeredOptimizer(g, {"name": to_opt, "butterfly_reduction": True, "verbose": False, "cutoff_time": 100})
        a, b = optimizer.optimize_layout()
        times1.append(str(a))
        optvals1.append(str(b))
        optimizer.butterfly_reduction = False
        a, b = optimizer.optimize_layout()
        times2.append(str(a))
        optvals2.append(str(b))
        print('-'*100)
    print("With butterfly reduction:")
    print('\t'.join(times1))
    print('\t'.join(optvals1))
    print("Without butterfly reduction:")
    print('\t'.join(times2))
    print('\t'.join(optvals2))


def transitivity_experiment():
    n_nodes = 67
    to_optimize = os.listdir(f"Rome-Lib/graficon{n_nodes}nodi")[:20]
    times1 = []
    optvals1 = []
    times2 = []
    optvals2 = []
    for i, to_opt in enumerate(to_optimize):
        g = layering.create_better_layered_graph(f"graficon{n_nodes}nodi/{to_opt}", 4, 2)[0]
        optimizer = LayeredOptimizer(g, {"name": to_opt, "transitivity_constraints": True, "verbose": False, "cutoff_time": 100})
        a, b = optimizer.optimize_layout()
        times1.append(str(a))
        optvals1.append(str(b))
        optimizer.transitivity_constraints = False
        a, b = optimizer.optimize_layout()
        times2.append(str(a))
        optvals2.append(str(b))
        print('-'*100)
    print("With transitivity constraints:")
    print('\t'.join(times1))
    print('\t'.join(optvals1))
    print("Without transitivity constraints:")
    print('\t'.join(times2))
    print('\t'.join(optvals2))


def fix_1_var_experiment():
    n_nodes = 67
    to_optimize = os.listdir(f"Rome-Lib/graficon{n_nodes}nodi")[:20]
    times1 = []
    optvals1 = []
    times2 = []
    optvals2 = []
    for to_opt in to_optimize:
        g = layering.create_better_layered_graph(f"graficon{n_nodes}nodi/{to_opt}", 4, 2)[0]
        optimizer = LayeredOptimizer(g, {"name": to_opt, "butterfly_reduction": False, "verbose": False, "cutoff_time": 100, "fix_one_var": True})
        a, b = optimizer.optimize_layout()
        times1.append(str(a))
        optvals1.append(str(b))
        optimizer.fix_one_var = False
        a, b = optimizer.optimize_layout()
        times2.append(str(a))
        optvals2.append(str(b))
        print('-' * 100)
    print("With fixing one reduction:")
    print('\t'.join(times1))
    print('\t'.join(optvals1))
    print("Without fixing one reduction:")
    print('\t'.join(times2))
    print('\t'.join(optvals2))


def write_file_name(filename, db):
    with open(f"data storage/{db}", 'a') as f:
        f.write(filename + '\n')


def randomly_select_files_for_exp(fname):
    dagmar_file_nums = random.sample(range(45), 8)
    dagmar_file_nums.extend(random.sample(range(20), 4))
    dagmar_files = [f"{1.6 if i < 8 else 2.6}/uniform_n{(k//10)*20+20}_e{(k//10)*32+32 if i < 8 else (k//10)*52+52}_i{k%10}.graphml" for i, k in enumerate(dagmar_file_nums)]
    rome_folds = random.choices(range(15, 90), k=25)
    rome_files = random.sample(range(50), 25)
    for file in dagmar_files:
        write_file_name(f"DAGmar/graphs/{file}", fname)
    for i, fold in enumerate(rome_folds):
        file = os.listdir(f"Rome-Lib/graficon{fold}nodi")[rome_files[i]]
        write_file_name(f"Rome-Lib/graficon{fold}nodi/{file}", fname)
    for file in random.sample(os.listdir("north"), 13):
        if file != "Graph.log":
            write_file_name(f"north/{file}", fname)
    with open(f"data storage/{fname}", 'r') as f:
        lines = list(f.readlines())
        random.shuffle(lines)
    with open(f"data storage/{fname}", 'w') as f:
        f.writelines(lines)


def make_altair_chart_for_ind_var():
    data = experiments.read_data_from_file("independent_var_study.csv", ',')
    data = [dat for dat in data if dat["opttime"] < 120 and dat["iterations"] > 0]
    print(len(data))
    for dat in data:
        dat['file'] = dat['file'][:dat['file'].index('/')]
        dat['xpc'] = dat['xvars'] + dat['cvars']
    vis.draw_altair_scatter(data, "xpc", "opttime", "file", "X-vars + c-vars", "Time (s)", "Decision variables vs time to optimize")
    vis.draw_altair_scatter(data, "xpc", "iterations", "file", "X-vars + c-vars", "Simplex iterations", "Decision variables vs iterations")


def starting_assignment_experiment():
    pass


def record_baseline_info(filename, start_idx):
    with open(f"data storage/{filename}", 'r') as f:
        i = start_idx - 1
        for line in f.readlines()[start_idx:]:
            i += 1
            g = read(line.removesuffix('\n'))
            opt = LayeredOptimizer(g, {})
            opt.return_full_data = True
            opt.fix_one_var = True
            opt.bendiness_reduction = False
            opt.aggro_presolve = True
            opt.xvar_branch_priority = True
            tup = opt.optimize_layout()
            with open(f"data storage/{filename}_info", 'a') as f2:
                f2.write(','.join(str(j) for j in [i, line.removesuffix('\n'), sum(1 for nd in g.nodes if not nd.is_anchor_node), len(g.nodes), len(g.edges), round(g.calculate_connectedness(), 3), tup[3], tup[4]]) + '\n')


if __name__ == '__main__':
    if len(sys.argv) < 2:
        # experiments.run_experiment(0, "data storage/experiment_set_50", exp_name="fix1var", param_to_set="fix_one_var", clear_files=True)
        # experiments.run_experiment(0, "data storage/experiment_set_50", exp_name="heuristic_start", param_to_set="heuristic_start", clear_files=True)
        # experiments.run_experiment(0, "data storage/experiment_set_50", exp_name="butterfly", param_to_set="butterfly_reduction", clear_files=True)
        # experiments.run_experiment(0, "data storage/experiment_set_50", exp_name="presolve", param_to_set="presolve", clear_files=True)
        # experiments.run_experiment(0, "data storage/experiment_set_50", exp_name="xvar_branch", param_to_set="priority", clear_files=True)
        # experiments.run_experiment(0, "data storage/experiment_set_50", exp_name="mip_relax", param_to_set="mip_relax", clear_files=True)
        # experiments.run_experiment(0, "data storage/experiment_set_50", exp_name="mirror_vars", param_to_set="mirror_vars", clear_files=True)

        experiments.make_altair_chart_all_three("baseline")

        # g = read("Rome-Lib/graficon17nodi/grafo224.17")
        # opt = LayeredOptimizer(g, {"priority": True, "strat_big_m": True})
        # opt.optimize_layout()

        # record_baseline_info("experiment_set_50", 50)

        # make_altair_chart_for_ind_var()

        # randomly_select_files_for_exp("experiment_set_50")

        # run_my_algorithm()

        # transitivity_experiment()

        # experiments.independent_var_experiment("independent_var_study_files", 80)

        # run_standard_version()
        # run_test_fix_x_vars()
        # cProfile.run("run_my_algorithm()", "my_algo_stats")
        # p = pstats.Stats("my_algo_stats")
        # cProfile.run("run_standard_version()", "standard_stats")
        # p = pstats.Stats("standard_stats")
        # p.strip_dirs().sort_stats(SortKey.CUMULATIVE).print_stats()
    elif len(sys.argv) < 10:
        run_all_rome_lib(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), bool(int(sys.argv[4])), bool(int(sys.argv[5])), sys.argv[6], bool(int(sys.argv[7])), savefile=sys.argv[8])
    elif len(sys.argv) < 11:
        run_all_rome_lib(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), bool(int(sys.argv[4])), bool(int(sys.argv[5])), sys.argv[6], bool(int(sys.argv[7])), savefile=sys.argv[8], shuffle=bool(int(sys.argv[9])))
    elif len(sys.argv) < 12:
        run_all_rome_lib(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), bool(int(sys.argv[4])), bool(int(sys.argv[5])), sys.argv[6], bool(int(sys.argv[7])), savefile=sys.argv[8], shuffle=bool(int(sys.argv[9])), target=sys.argv[10])
    else:
        run_all_rome_lib(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), bool(int(sys.argv[4])), bool(int(sys.argv[5])), sys.argv[6], bool(int(sys.argv[7])), savefile=sys.argv[8], shuffle=bool(int(sys.argv[9])), target=sys.argv[10], subgraph_reduction=bool(int(sys.argv[11])))
