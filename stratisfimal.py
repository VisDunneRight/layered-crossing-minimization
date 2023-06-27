import os
import csv
import random
from src import vis
import experiments
from src.graph import *
from src.read_data import *
from src.optimization import LayeredOptimizer


random.seed(22)


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


def run_stratisfimal_layout(graph_file):
    optimizer = LayeredOptimizer(graph_file)
    optimizer.vertical_transitivity = True
    optimizer.direct_transitivity = True
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
    optimizer.mirror_vars = True
    # optimizer.symmetry_constraints = False
    optimizer.return_experiment_data = True
    return optimizer.optimize_layout()


def run_my_layout_algorithm(graph_file):
    optimizer = LayeredOptimizer(graph_file)
    optimizer.vertical_transitivity = True
    optimizer.fix_one_var = True
    # optimizer.butterfly_reduction = True
    # optimizer.heuristic_start = True
    optimizer.mip_relax = True
    optimizer.xvar_branch_priority = True
    optimizer.aggro_presolve = True
    optimizer.mirror_vars = True
    optimizer.return_experiment_data = True
    return optimizer.optimize_layout()


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
            pair = optimizer.optimize_with_starting_assignments(optimizer.__generate_random_vars_to_fix(2))
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
    to_optimize = os.listdir(f"Rome-Lib/graficon{n_nodes}nodi")[0:1]
    for to_opt in to_optimize:
        optimizer = LayeredOptimizer(f"Rome-Lib/graficon{n_nodes}nodi/{to_opt}")
        optimizer.fix_one_var = True
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


def randomly_select_50_files(fname):
    with open('data storage/direct_transitivity/baseline_60.csv') as fd:
        reader = csv.reader(fd)
        rome_nums = random.sample(list(range(1, 9855)), 35)
        dagmar_nums = random.sample(list(range(9855, 9900)), 5)
        north_nums = random.sample(list(range(9900, 11177)), 10)
        gnames = [row[1]+'\n' for idx, row in enumerate(reader) if idx in rome_nums+dagmar_nums+north_nums]
    with open(f"data storage/{fname}", 'w') as f:
        f.writelines(gnames)
def make_altair_chart_for_ind_var():
    data = experiments.read_data_from_file("independent_var_study.csv", ',')
    data = [dat for dat in data if dat["opttime"] < 120 and dat["iterations"] > 0]
    print(len(data))
    for dat in data:
        dat['file'] = dat['file'][:dat['file'].index('/')]
        dat['xpc'] = dat['xvars'] + dat['cvars']
    vis.draw_altair_scatter(data, "xpc", "opttime", "file", "X-vars + c-vars", "Time (s)", "Decision variables vs time to optimize")
    vis.draw_altair_scatter(data, "xpc", "iterations", "file", "X-vars + c-vars", "Simplex iterations", "Decision variables vs iterations")


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


def case_study_graph_experiment(command):
    my_vals = []
    strat_vals = []
    junger_vals = []
    sankey_vals = []
    control_flow_file = f"control-flow-graphs/{command}/dbg.main.dot"
    for i in range(5):
        m1 = run_my_layout_algorithm(control_flow_file)
        m2 = run_junger_polyhedral_layout(control_flow_file)
        m3 = run_optimal_sankey_layout(control_flow_file)
        m4 = run_stratisfimal_layout(control_flow_file)
        my_vals.append(m1[5] + m1[8])
        junger_vals.append(m2[5] + m2[8])
        sankey_vals.append(m3[5] + m3[8])
        strat_vals.append(m4[5] + m4[8])
    print(sum(my_vals)/5, my_vals)
    print(sum(strat_vals)/5, strat_vals)
    print(sum(junger_vals)/5, junger_vals)
    print(sum(sankey_vals)/5, sankey_vals)


def my_fn(s):
    if s[1][0] == "R":
        add_val = -100000
    elif s[1][0] == "D":
        add_val = 0
    else:
        add_val = 100000
    return add_val + int(s[3])


""" find bucket of files size n in sorted exp file, works regardless of sorted """
def bucket_lines_in_data(file, bucket_size):
    lines_in_file = []
    seen_files = set()
    with open(file, 'r') as fd1:
        rdr = csv.reader(fd1)
        next(rdr)
        for ln in rdr:
            if bucket_size <= int(ln[3]) < bucket_size + 10 and ln[1] not in seen_files:
                lines_in_file.append(ln)
                seen_files.add(ln[1])
    lines_in_file.sort(key=my_fn)
    return lines_in_file


def get_all_files_in_bucket(bucket_size):
    with open("data storage/all_g_sorted.txt", 'r') as fd1:
        collect_lines = False
        filenames = []
        for line in fd1.readlines():
            if line[0] == "T":
                if int(line[line.index('[') + 1:line.index(',')]) == bucket_size:
                    collect_lines = True
                else:
                    collect_lines = False
            elif collect_lines:
                filenames.append(line[:line.index(',')])
    return filenames


def calc_if_bucket_donezo(datapts):
    timedout = sum((1 for pt in datapts if float(pt[10]) > 60))
    return True if timedout / len(datapts) >= 0.25 else False


def run_thing():
    """ find missing entries, run experiment, write to new file, cut off once >50% in bucket timeout """
    key1 = ["fix_one_var", "butterfly_reduction", "heuristic_start", "presolve", "priority", "mip_relax", "mirror_vars"]
    key2 = ["fix1var_60", "butterfly_60", "heuristic_60", "presolve_60", "xvar_branch_60", "mip_relax_60", "symmetry_60"]
    for j, inp1 in enumerate(["junger_basic", "vertical_transitivity", "redundancy"]):
        for i, inp2 in enumerate(key2):
            fname = f"{inp1}/{inp2}"
            experiments.insert_one(f"{fname}_new.csv", ["Index", "File", "Nodes", "Total Nodes", "Butterflies", "X-vars", "C-vars", "Total vars", "Total constraints", "Crossings", "Opttime", "Work", "Nodes visited", "Setup Time"])
            curindex = 0
            for bsize in range(10, 17141, 10):
                bfiles = bucket_lines_in_data("data storage/"+fname+".csv", bsize)
                bfnames = [bfl[1] for bfl in bfiles]
                all_bfiles = get_all_files_in_bucket(bsize)
                if len(all_bfiles) > 0:
                    for check_file in all_bfiles:
                        if check_file not in bfnames:
                            parameters = [key1[i], "direct_transitivity" if j % 2 == 0 else "baseline", "vertical_transitivity" if j > 0 else "baseline"]
                            experiments.run_one_graph(check_file, f"{fname[fname.index('/') + 1:]}_new", 60, parameters, curindex)
                        else:
                            bfiles[bfnames.index(check_file)][0] = curindex
                            experiments.insert_one(f"{fname}_new.csv", bfiles[bfnames.index(check_file)])
                        curindex += 1
                    if calc_if_bucket_donezo(bfiles):
                        print(f"{inp1} with switch {inp2} cutoff at bucket size {bsize}")
                        break


def run_select_multi_param():
    key = ["", "fix_one_var", "butterfly_reduction", "heuristic_start", "presolve", "priority", "mip_relax", "mirror_vars"]
    for form in [("vertical_transitivity", "12356"), ("direct_transitivity", "1357"), ("redundancy", "13467"), ("vertical_transitivity", "12567"), ("vertical_transitivity", "12456")]:
        fname = f"multi_param_results/{form[0]}_{form[1]}"
        experiments.insert_one(f"{fname}.csv",
                               ["Index", "File", "Nodes", "Total Nodes", "Butterflies", "X-vars", "C-vars",
                                "Total vars", "Total constraints", "Crossings", "Opttime", "Work", "Nodes visited",
                                "Setup Time"])
        curindex = 0
        for bsize in range(10, 17141, 10):
            all_bfiles = get_all_files_in_bucket(bsize)
            if len(all_bfiles) > 0:
                bucket_results = []
                for check_file in all_bfiles:
                    parameters = [key[int(i)] for i in form[1]] + ["direct_transitivity" if form[0] == "direct_transitivity" or form[0] == "redundancy" else "", "vertical_transitivity" if form[0] == "vertical_transitivity" or form[0] == "redundancy" else ""]
                    bucket_results.append(experiments.run_one_graph(check_file, f"{fname}", 60, parameters, curindex))
                    curindex += 1
                if calc_if_bucket_donezo(bucket_results):
                    print(f"Cutoff at bucket size {bsize}")
                    break


def calculate_success_rate_by_bucket(file):
    with open(file, 'r') as fd:
        rdr = csv.reader(fd)
        next(rdr)
        success_rates = {}
        bucket_size = 10
        running_count = 0
        running_success = 0
        for line in rdr:
            if int(line[3]) >= bucket_size + 10:
                success_rates[bucket_size] = running_success / running_count
                bucket_size += 10
                running_count = 0
                running_success = 0
            if int(line[10]) < 60:
                running_success += 1
            running_count += 1
        print(sorted([(hj, success_rates[hj]) for hj in success_rates], key=lambda x: x[0]))


def sample_10_percent():
    with open("data storage/5_percent_all_g_sorted.txt", 'w') as fd:
        estimate = 0
        for i in range(10, 1110, 10):
            with open("data storage/all_g_sorted.txt", 'r') as fd1:
                collect_lines = False
                files = []
                for line in fd1.readlines():
                    if line[0] == "T":
                        if int(line[line.index('[') + 1:line.index(',')]) == i:
                            collect_lines = True
                        else:
                            collect_lines = False
                    elif collect_lines:
                        files.append(line)
            if len(files) >= 5:
                new_files = random.sample(files, round(len(files)/20))
                for file in new_files:
                    fd.write(file)
                    with open("data storage/multi_param_results/junger_trans_1357.csv", 'r') as f:
                        for row in csv.reader(f):
                            if row[1] == file[:file.index(',')]:
                                estimate += float(row[10])
    print(estimate/60, "minutes")



def fix_27():
    for fil in os.listdir("data storage/vertical_transitivity/all_combos_5percent"):
        with open(f"data storage/vertical_transitivity/all_combos_5percent/{fil}", 'r') as fd1:
            lines = []
            for lin in csv.reader(fd1):
                lines.append(lin)
        with open(f"data storage/vertical_transitivity/all_combos_5percent/{fil}", 'w') as fd1:
            wrt = csv.writer(fd1)
            wrt.writerow(["Index", "File", "Nodes", "Total Nodes", "Butterflies", "X-vars", "C-vars", "Total vars", "Total constraints", "Crossings", "Opttime", "Work", "Nodes visited", "Setup Time"])
            if "2" in fil and "7" in fil:
                for lin in lines[:-1]:
                    if int(lin[4]) > 0:
                        lin[9] = str(int(lin[9]) - int(lin[4]))
                    wrt.writerow(lin)
            else:
                for lin in lines[:-1]:
                    wrt.writerow(lin)


def tfix():
    for fil in os.listdir("data storage/both_combined/all_combos_5percent"):
        if fil == "exp0.csv":
            with open(f"data storage/both_combined/all_combos_5percent/{fil}", 'r') as fd1:
                lines = []
                for lin in fd1.readlines():
                    lines.append(lin)
                lines = [ln for ln in lines if ln != "" and ln != "\n"]
            with open(f"data storage/both_combined/all_combos_5percent/{fil}", 'w') as fd1:
                for row in lines:
                    fd1.write(row)


def find_best():
    avgs = []
    meds = []
    for iv in ["direct_transitivity", "vertical_transitivity", "both_combined"]:
        for fl in os.listdir(f"data storage/{iv}/all_combos_5percent"):
            medlist = []
            with open(f"data storage/{iv}/all_combos_5percent/{fl}", 'r') as fd:
                rd = csv.reader(fd)
                next(rd)
                avsum = 0
                ct = 0
                success = 0
                cr = 0
                for ln in rd:
                    if len(ln) > 1:
                        avsum += float(ln[10]) if float(ln[10]) < 60 else 60
                        medlist.append(float(ln[10]) if float(ln[10]) < 60 else 60)
                        if float(ln[10]) < 60:
                            success += 1
                        ct += 1
                        cr += int(ln[9])
                medlist.sort()
                meds.append((f"{iv}/{fl}", medlist[len(medlist)//2] / ct, cr / ct, success))
                avgs.append((f"{iv}/{fl}", avsum / ct, cr / ct, success))
    avgs.sort(key=lambda x: (-x[3], x[1]))
    meds.sort(key=lambda x: (-x[3], x[1]))
    print(avgs)
    print(meds)
    # return avgs
    return meds


def find_best_same_files():
    avgs = []
    min_nsolved = 309
    # for iv in ["junger_basic", "vertical_transitivity"]:
    #     for fl in os.listdir(f"data storage/{iv}/all_combos_5percent"):
    #         with open(f"data storage/{iv}/all_combos_5percent/{fl}", 'r') as fd:
    #             fdl = len(fd.readlines())
    #             if fdl < min_nsolved:
    #                 # print(iv, fl, fdl)
    #                 min_nsolved = fdl
    print(min_nsolved)
    for iv in ["junger_basic", "vertical_transitivity"]:
        for fl in os.listdir(f"data storage/{iv}/all_combos_5percent"):
            with open(f"data storage/{iv}/all_combos_5percent/{fl}", 'r') as fd:
                rd = csv.reader(fd)
                next(rd)
                avsum = 0
                ct = 0
                success = 0
                cr = 0
                for ln in rd:
                    if len(ln) > 1 and ct < min_nsolved:
                        avsum += float(ln[10]) if float(ln[10]) < 60 else 60
                        if float(ln[10]) < 60:
                            success += 1
                        ct += 1
                        cr += int(ln[9])
                for i in range(ct, min_nsolved):
                    avsum += 60
                    ct += 1
                avgs.append((f"{iv}/{fl}", avsum / ct, cr / ct, success))
    avgs.sort(key=lambda x: (-x[3], x[1]))
    print(avgs)
    return avgs


def find_best_timeout_uncalculated_graphs():
    avgs = []
    max_nsolved = 322
    print(max_nsolved)
    for iv in ["junger_basic", "vertical_transitivity"]:
        for fl in os.listdir(f"data storage/{iv}/all_combos_5percent"):
            with open(f"data storage/{iv}/all_combos_5percent/{fl}", 'r') as fd:
                rd = csv.reader(fd)
                next(rd)
                avsum = 0
                ct = 0
                success = 0
                cr = 0
                for ln in rd:
                    if len(ln) > 1:
                        avsum += float(ln[10]) if float(ln[10]) < 60 else 60
                        if float(ln[10]) < 60:
                            success += 1
                        ct += 1
                        cr += int(ln[9])
                for i in range(ct, max_nsolved):
                    avsum += 60
                    ct += 1
                avgs.append((f"{iv}/{fl}", avsum / ct, cr / ct, success))
    avgs.sort(key=lambda x: x[1])
    print(avgs)
    return avgs


def rank_switches_all_combos():
    # avgs = find_best()
    avgs = find_best_same_files()
    # avgs = find_best_timeout_uncalculated_graphs()
    for switch in range(1, 8):
        swcrwins = 0
        swrunwins = 0
        swsuccwins = 0
        swct = 0
        for fl in avgs:
            if fl[0][:5] == "strat" and str(switch) in fl[0]:
                noswitch = fl[0].replace(str(switch), '0') if not fl[0][fl[0].index('.') - 2].isnumeric() else fl[0].replace(str(switch), '')
                nosentry = [f for f in avgs if f[0] == noswitch][0]
                if fl[1] < nosentry[1]:
                    swrunwins += 1
                if fl[2] < nosentry[2]:
                    swcrwins += 1
                # else:
                #     print(fl, nosentry)
                if fl[3] >= nosentry[3]:
                    swsuccwins += 1
                swct += 1
        print(switch, swsuccwins, swrunwins, swcrwins, swct, swsuccwins / swct, swrunwins / swct, swcrwins / swct)


def latex_print_ranking():
    all_results = find_best()
    all_r_vert = [entry for entry in all_results if entry[0][0] == "v"]
    all_r_direct = [entry for entry in all_results if entry[0][0] == "d"]
    all_r_both = [entry for entry in all_results if entry[0][0] == "b"]
    alphad = {1: "C", 2: "D", 3: "I", 4: "E", 5: "F", 6: "G", 7: "H"}
    alphae = {1: "1", 2: "2", 3: "7", 4: "3", 5: "4", 6: "5", 7: "6"}
    switch_conversion = {i: "\hyperref[sec:switchX]{\switchsboxonY}".replace("X", str(i)).replace("Y", alphad[i]) for i in range(1, 8)}
    switch_conversion.update({0: "baseline"})
    for i in range(len(all_r_vert)):
        to_conv1 = []
        to_conv2 = []
        to_conv3 = []
        for j in range(1,8):
            if alphae[j] in all_r_vert[i][0]:
                to_conv1.append(switch_conversion[j])
            else:
                to_conv1.append(' ')
            if alphae[j] in all_r_direct[i][0]:
                to_conv2.append(switch_conversion[j])
            else:
                to_conv2.append(' ')
            if alphae[j] in all_r_both[i][0]:
                to_conv3.append(switch_conversion[j])
            else:
                to_conv3.append(' ')
        # if "7" in all_r_vert[i][0]:
        #     to_conv1.insert(2, switch_conversion[7])
        # else:
        #     to_conv1.insert(2, ' ')
        # if "7" in all_r_direct[i][0]:
        #     to_conv2.insert(2, switch_conversion[7])
        # else:
        #     to_conv2.insert(2, ' ')
        # if "7" in all_r_both[i][0]:
        #     to_conv3.insert(2, switch_conversion[7])
        # else:
        #     to_conv3.insert(2, ' ')
        # if len(to_conv1) > 4:
        #     to_conv1.insert(4, '\\\\')
        # if len(to_conv2) > 4:
        #     to_conv2.insert(4, '\\\\')
        # if len(to_conv3) > 4:
        #     to_conv3.insert(4, '\\\\')
        # name = '{' + ''.join(to_conv1) + '}'
        # name2 = '{' + ''.join(to_conv2) + '}'
        # name3 = '{' + ''.join(to_conv3) + '}'
        name = ' & '.join(to_conv1)
        name2 = ' & '.join(to_conv2)
        name3 = ' & '.join(to_conv3)
        print(str(i+1), ' & ', name, ' & ', str(all_r_vert[i][3]), ' & & ', str(i+1), ' & ', name2, ' & ', str(all_r_direct[i][3]), ' & & ', str(i+1), ' & ', name3, ' & ', str(all_r_both[i][3]), ' \\\\')


def avg_comb_num_switch_ranking():
    all_results = find_best()
    all_r_vert = [entry for entry in all_results if entry[0][0] == "s"]
    avg_rank = [0]*8
    cts = [0]*8
    zeropos = 0
    for i in range(8):
        if i == 0:
            for j, val in enumerate(all_r_vert):
                if "0" in val[0]:
                    avg_rank[i] = j
                    cts[i] += 1
                    zeropos = j
        else:
            for j, val in enumerate(all_r_vert):
                if "1" not in val[0] and j != zeropos and len(val[0][val[0].index('p')+1:val[0].index('.')]) == i:
                    cts[i] += 1
                    avg_rank[i] += j
    print(avg_rank)
    print(cts)
    print([a/cts[i] for i, a in enumerate(avg_rank[:7])])


def group_success_calc(transitivity):
    with open(f"data storage/{transitivity}/all_combos_5percent/exp0.csv", 'r') as fd:
        rdr = csv.reader(fd)
        next(rdr)
        basevals = []
        for row in rdr:
            basevals.append(float(row[10]))
    all_cutoffs = {120: 0, 130: 0, 140: 0, 150: 0, 160: 0, 170: 0, 180: 0, 190: 0, 200: 0, 210: 0, 220: 0, 230: 0}
    for fil in os.listdir(f"data storage/{transitivity}/all_combos_5percent"):
        with open(f"data storage/{transitivity}/all_combos_5percent/{fil}", 'r') as fd:
            rdr = csv.reader(fd)
            next(rdr)
            ftimein = 0
            curbin = 10
            binct = 0
            for i, row in enumerate(rdr):
                if int(row[3])//20*20 > curbin:
                    # print(curbin, ftimein, binct, ftimein/binct)
                    if ftimein / binct < 0.75:
                        break
                    curbin = int(row[3])//20*20
                    ftimein = 0
                    binct = 0
                if i < len(basevals) and basevals[i] < 60 and float(row[10]) < 60:
                    ftimein += 1
                binct += 1
            all_cutoffs[curbin] += 1
    print(all_cutoffs)


def get_plotted_points_in_ind_switch_figures(formulation, switch):
    with open(f"data storage/{formulation}/baseline_60_new.csv", 'r') as fd:
        rdr = csv.reader(fd)
        next(rdr)
        basevals = []
        for row in rdr:
            basevals.append(float(row[10]))
    with open(f"data storage/{formulation}/{switch}_60_new.csv", 'r') as fd:
        rdr = csv.reader(fd)
        next(rdr)
        cur_tdiffs = []
        curbin = 10
        binct = 0
        for i, row in enumerate(rdr):
            if int(row[3])//10*10 > curbin:
                cur_tdiffs.sort()
                top_25 = cur_tdiffs[3*len(cur_tdiffs)//4:]
                bottom_25 = cur_tdiffs[:len(cur_tdiffs)//4]
                print(curbin, binct, sum(cur_tdiffs)/len(cur_tdiffs), sum(top_25)/len(top_25), sum(bottom_25)/len(bottom_25))
                curbin = int(row[3])//10*10
                binct = 0
            if 0 < basevals[i] < 60 and 0 < float(row[10]) < 60:
                cur_tdiffs.append(math.log(float(row[10])/basevals[i]))
                binct += 1


if __name__ == '__main__':
    # rank_switches_all_combos()

    # latex_print_ranking()
    # group_success_calc("junger_basic")
    # avg_comb_num_switch_ranking()

    # with open("data storage/5_percent_all_g_sorted.txt12", 'r') as fd:
    #     romemax, romemin = 0, 100
    #     dagmarmax, dagmarmin = 0, 100
    #     northmax, northmin = 0, 100
    #     northct, dagmarct, romect = 0, 0, 0
    #     for line in fd.readlines():
    #         if line[0] != "T":
    #             lns = line.split(',')
    #             tnodes = int(lns[1])
    #             col = lns[0][:lns[0].index('/')]
    #             if col == "Rome-Lib":
    #                 if tnodes > romemax:
    #                     romemax = tnodes
    #                 if tnodes < romemin:
    #                     romemin = tnodes
    #                 romect += 1
    #             if col == "DAGmar":
    #                 if tnodes > dagmarmax:
    #                     dagmarmax = tnodes
    #                 if tnodes < dagmarmin:
    #                     dagmarmin = tnodes
    #                 dagmarct += 1
    #             if col == "north":
    #                 if tnodes > northmax:
    #                     northmax = tnodes
    #                 if tnodes < northmin:
    #                     northmin = tnodes
    #                 northct += 1
    # print(romemax, romemin, dagmarmax, dagmarmin, northmax, northmin)
    # print(romect, dagmarct, northct)

    # case_study_graph_experiment("df")
    # run_select_multi_param()
    # sample_10_percent()

    # with open("data storage/5_percent_all_g_sorted.txt", 'r') as fd:
    #     gfiles = []
    #     for line in fd.readlines():
    #         gfiles.append(line[:line.index(',')])
    # experiments.all_combinations_experiment(gfiles, "all_combos_5percent")

    # fix_27()
    # tfix()

    # experiments.run_experiment((1,0), cutoff_time=60, exp_name="baseline", param_to_set="baseline", clear_files=False, max_timeout=15)
    # experiments.run_experiment((2,58), cutoff_time=60, exp_name="fix1var", param_to_set="fix_one_var", clear_files=False, max_timeout=5)
    # experiments.run_experiment((0,0), cutoff_time=60, exp_name="butterfly", param_to_set="butterfly_reduction", clear_files=True, max_timeout=3)
    # experiments.run_experiment((0,0), cutoff_time=60, exp_name="heuristic", param_to_set="heuristic_start", clear_files=True, max_timeout=3)
    # experiments.run_experiment((0,0), cutoff_time=60, exp_name="presolve", param_to_set="presolve", clear_files=True, max_timeout=3)
    # experiments.run_experiment((0,0), cutoff_time=60, exp_name="xvar_branch", param_to_set="priority", clear_files=True, max_timeout=3)
    # experiments.run_experiment((0,0), cutoff_time=60, exp_name="mip_relax", param_to_set="mip_relax", clear_files=True, max_timeout=3)
    # experiments.run_experiment((0,0), cutoff_time=60, exp_name="symmetry", param_to_set="mirror_vars", clear_files=True, max_timeout=3)

    # experiments.run_experiment(0, "data storage/experiment_set_50", exp_name="fix1var", param_to_set="fix_one_var", clear_files=True)
    # experiments.run_experiment(0, "data storage/experiment_set_50", exp_name="heuristic_start", param_to_set="heuristic_start", clear_files=True)
    # experiments.run_experiment(0, "data storage/experiment_set_50", exp_name="butterfly", param_to_set="butterfly_reduction", clear_files=True)
    # experiments.run_experiment(0, "data storage/experiment_set_50", exp_name="presolve", param_to_set="presolve", clear_files=True)
    # experiments.run_experiment(0, "data storage/experiment_set_50", exp_name="xvar_branch", param_to_set="priority", clear_files=True)
    # experiments.run_experiment(0, "data storage/experiment_set_50", exp_name="mip_relax", param_to_set="mip_relax", clear_files=True)
    # experiments.run_experiment(0, "data storage/experiment_set_50", exp_name="mirror_vars", param_to_set="mirror_vars", clear_files=True)

    # record_baseline_info("experiment_set_50", 50)

    # randomly_select_files_for_exp("experiment_set_50")

    # experiments.independent_var_experiment("independent_var_study_files", 80)

    run_standard_version()
    # run_test_fix_x_vars()
    # cProfile.run("run_my_algorithm()", "my_algo_stats")
    # p = pstats.Stats("my_algo_stats")
    # cProfile.run("run_standard_version()", "standard_stats")
    # p = pstats.Stats("standard_stats")
    # p.strip_dirs().sort_stats(SortKey.CUMULATIVE).print_stats()