import os
import itertools
import pickle
import csv
import random
import sys

from src import optimization, read_data, vis, motifs


def get_list_of_files(storage_file):
    filenames = []
    with open("data storage/" + storage_file, 'r') as f:
        for fname in f.readlines():
            filenames.append(fname.removesuffix('\n'))
    return filenames


def pickle_data(data, name):
    p_object = {"iter": len(data), "data": data}
    with open("data storage/" + name + "_pickled.bin", 'wb') as pickle_file:
        pickle.dump(p_object, pickle_file)


def get_pickled_data(name):
    with open("data storage/" + name + "_pickled.bin", 'rb') as fd:
        p_object = pickle.load(fd)
    return p_object["data"], p_object["iter"]


def insert_data(name, entries):
    with open("data storage/" + name, 'a', newline='') as f:
        wrt = csv.writer(f)
        for entry in entries:
            wrt.writerow(entry)


def insert_one(name, entry):
    with open("data storage/" + name, 'a', newline='') as f:
        wrt = csv.writer(f)
        wrt.writerow(entry)


def clear_file(name):
    f = open(f"data storage/{name}", 'w')
    f.close()


def basic_info(g):
    return [sum(1 for node in g if not node.is_anchor_node), len(g.nodes), motifs.count_butterflies(g)]


def read_data_from_file(name, split_char, header_select=None):
    data = []
    with open("data storage/" + name, 'r') as f:
        l1 = f.readline().removesuffix('\n')

        headers = [head for head in l1.split(split_char)]
        # data.append([int(val) if val.isnumeric() else float(val) if not val.isalpha() else val])
        for entry in f.readlines():
            if header_select is None:
                data.append({headers[i]: int(val) if val.isnumeric() else float(val) if val.replace('.', '').isnumeric() else val for i, val in enumerate(entry.removesuffix('\n').split(split_char))})
            else:
                data.append({headers[i]: int(val) if val.isnumeric() else float(val) if val.replace('.', '').isnumeric() else val for i, val in enumerate(entry.removesuffix('\n').split(split_char)) if i in header_select})
    return data


def tag_data(data_list, tagname, tag):
    for data in data_list:
        data[tagname] = tag


def baseline_experiment(start_idx, filename):
    with open(filename, 'r') as f:
        gfiles = [gname.removesuffix('\n') for gname in f.readlines()]
    for to_opt in gfiles[start_idx:]:
        g = read_data.read(to_opt)
        optimizer = optimization.LayeredOptimizer(g, {"name": to_opt, "cutoff_time": 600, "return_experiment_data": True, "stratisfimal_yvars": True})
        result = optimizer.optimize_layout()
        insert_data("strat_baseline.csv", [result])
        optimizer.junger_ec, optimizer.stratisfimal_y_vars = True, False
        result = optimizer.optimize_layout()
        insert_data("junger_baseline.csv", [result])
        optimizer.mirror_vars, optimizer.junger_ec = True, False
        result = optimizer.optimize_layout()
        insert_data("sankey_baseline.csv", [result])


def independent_var_experiment(file_name, start_ind):
    results = []
    for i, file in enumerate(get_list_of_files(file_name)[start_ind:]):
        print(file)
        result = [i+start_ind, file]
        g = read_data.read(file)
        result.extend((sum(1 for n in g.nodes if not n.is_anchor_node), len(g.nodes), len(g.edges)))
        opt = optimization.LayeredOptimizer(g, {"return_full_data": True, "cutoff_time": 120})
        result.extend(opt.optimize_layout())
        results.append(result)
        if i % 10 == 0:
            insert_data("independent_var_study.csv", results)
            results.clear()
    insert_data("independent_var_study.csv", results)


def fix_1_var_experiment(start_idx, filename):
    with open(filename, 'r') as f:
        gfiles = [gname.removesuffix('\n') for gname in f.readlines()]
    for to_opt in gfiles[start_idx:]:
        g = read_data.read(to_opt)
        optimizer = optimization.LayeredOptimizer(g, {"name": to_opt, "cutoff_time": 600, "symmetry_breaking": True, "return_experiment_data": True, "stratisfimal_yvars": True})
        result = optimizer.optimize_layout()
        insert_data("strat_fix1.csv", [result])
        optimizer.junger_ec, optimizer.stratisfimal_y_vars = True, False
        result = optimizer.optimize_layout()
        insert_data("junger_fix1.csv", [result])
        optimizer.mirror_vars, optimizer.junger_ec = True, False
        result = optimizer.optimize_layout()
        insert_data("sankey_fix1.csv", [result])


def get_all_graphs():
    all_g = []
    for i in range(10, 101):
        for file in os.listdir(f"Rome-Lib/graficon{i}nodi"):
            all_g.append(f"Rome-Lib/graficon{i}nodi/" + file)
    north_gs = sorted(list(os.listdir("north")), key=lambda fil: int(fil[2:(5 if fil[4] == '0' else 4)]))
    for i in range(len(north_gs)):
        north_gs[i] = "north/" + north_gs[i]
    north_gs.remove("north/g.57.26.graphml")   # skip this one graph that takes over an hour to insert the variables and constraints
    all_g.extend(north_gs)
    return all_g


def run_experiment(start_idx, cutoff_time, exp_name, param_to_set, clear_files, max_timeout):
    # with open("", 'r') as f:
    #     gfiles = [gname.removesuffix('\n') for gname in f.readlines()]
    gfiles = get_all_graphs()
    if clear_files:
        if start_idx != (0, 0):
            print("Something's wrong here...")
            return
        clear_file(f"junger_basic/{exp_name}_{cutoff_time}.csv")
        clear_file(f"vertical_transitivity/{exp_name}_{cutoff_time}.csv")
        clear_file(f"redundancy/{exp_name}_{cutoff_time}.csv")
        insert_one(f"junger_basic/{exp_name}_{cutoff_time}.csv", ["Index", "File", "Nodes", "Total Nodes", "Butterflies", "X-vars", "C-vars", "Total vars", "Total constraints", "Crossings", "Opttime", "Status", "Nodes visited", "Setup Time"])
        insert_one(f"vertical_transitivity/{exp_name}_{cutoff_time}.csv", ["Index", "File", "Nodes", "Total Nodes", "Butterflies", "X-vars", "C-vars", "Total vars", "Total constraints", "Crossings", "Opttime", "Status", "Nodes visited", "Setup Time"])
        insert_one(f"redundancy/{exp_name}_{cutoff_time}.csv", ["Index", "File", "Nodes", "Total Nodes", "Butterflies", "X-vars", "C-vars", "Total vars", "Total constraints", "Crossings", "Opttime", "Status", "Nodes visited", "Setup Time"])
    for j in range(start_idx[0], 3):
        junger_timedout = 0
        strat_timedout = 0
        redundant_timedout = 0
        for i, to_opt in enumerate(gfiles[j][start_idx[1]:]):
            if junger_timedout >= max_timeout and strat_timedout >= max_timeout and redundant_timedout >= max_timeout:
                break
            print(f"{i+start_idx[1] + 1} / {len(gfiles[j])}")
            g = read_data.read(to_opt)
            base_info = basic_info(g)
            optimizer = optimization.LayeredOptimizer(g, {"cutoff_time": cutoff_time, "return_experiment_data": True, "direct_transitivity": True, param_to_set: True})
            result = optimizer.optimize_layout()
            insert_one(f"junger_basic/{exp_name}_{cutoff_time}.csv", [i+start_idx[1], to_opt] + base_info + [j for j in result])
            if result[5] >= cutoff_time or junger_timedout >= max_timeout:
                junger_timedout += 1
            else:
                junger_timedout = 0
            optimizer.direct_transitivity, optimizer.vertical_transitivity = False, True
            result = optimizer.optimize_layout()
            insert_one(f"vertical_transitivity/{exp_name}_{cutoff_time}.csv", [i+start_idx[1], to_opt] + base_info + [j for j in result])
            if result[5] >= cutoff_time or strat_timedout >= max_timeout:
                strat_timedout += 1
            else:
                strat_timedout = 0
            optimizer.direct_transitivity = True
            result = optimizer.optimize_layout()
            insert_one(f"redundancy/{exp_name}_{cutoff_time}.csv", [i+start_idx[1], to_opt] + base_info + [j for j in result])
            if result[5] >= cutoff_time or redundant_timedout >= max_timeout:
                redundant_timedout += 1
            else:
                redundant_timedout = 0


def run_multi_param_experiment(start_idx, graphs_file, cutoff_time, exp_name, params_to_set, clear_files):
    with open(graphs_file, 'r') as f:
        gfiles = [gname.removesuffix('\n') for gname in f.readlines()]
    if clear_files:
        if start_idx != 0:
            print("Something's wrong here...")
            return
        clear_file(f"junger_basic/{exp_name}_{cutoff_time}.csv")
        clear_file(f"vertical_transitivity/{exp_name}_{cutoff_time}.csv")
        clear_file(f"redundancy/{exp_name}_{cutoff_time}.csv")
        insert_one(f"junger_basic/{exp_name}_{cutoff_time}.csv", ["Index", "File", "Nodes", "Total Nodes", "Butterflies", "X-vars", "C-vars", "Total vars", "Total constraints", "Crossings", "Opttime", "Status", "Nodes visited", "Setup Time"])
        insert_one(f"vertical_transitivity/{exp_name}_{cutoff_time}.csv", ["Index", "File", "Nodes", "Total Nodes", "Butterflies", "X-vars", "C-vars", "Total vars", "Total constraints", "Crossings", "Opttime", "Status", "Nodes visited", "Setup Time"])
        insert_one(f"redundancy/{exp_name}_{cutoff_time}.csv", ["Index", "File", "Nodes", "Total Nodes", "Butterflies", "X-vars", "C-vars", "Total vars", "Total constraints", "Crossings", "Opttime", "Status", "Nodes visited", "Setup Time"])
    for i, to_opt in enumerate(gfiles[start_idx:]):
        print(f"{i + start_idx + 1} / {len(gfiles)}")
        g = read_data.read(to_opt)
        base_info = basic_info(g)
        params = {param: True for param in params_to_set}
        params.update({"cutoff_time": cutoff_time, "return_experiment_data": True, "direct_transitivity": True})
        optimizer = optimization.LayeredOptimizer(g, params)
        result = optimizer.optimize_layout()
        insert_one(f"junger_basic/{exp_name}_{cutoff_time}.csv", [i + start_idx, to_opt] + base_info + [j for j in result])
        optimizer.direct_transitivity, optimizer.vertical_transitivity = False, True
        result = optimizer.optimize_layout()
        insert_one(f"vertical_transitivity/{exp_name}_{cutoff_time}.csv", [i + start_idx, to_opt] + base_info + [j for j in result])
        optimizer.direct_transitivity = True
        result = optimizer.optimize_layout()
        insert_one(f"redundancy/{exp_name}_{cutoff_time}.csv", [i + start_idx, to_opt] + base_info + [j for j in result])


def all_combinations_experiment(folder_to_write):
    key = ["symmetry_breaking", "butterfly_reduction", "heuristic_start", "presolve", "priority", "mip_relax", "mirror_vars"]
    with open("data storage/5_percent_all_g_sorted.txt", 'r') as fd:
        graphs_files = []
        for line in fd.readlines():
            graphs_files.append(line[:line.index(',')])
    for form in ["vertical_transitivity", "direct_transitivity", "both_combined"]:
        for combo in list(itertools.chain.from_iterable(itertools.combinations(key, r) for r in range(len(key)+1))):
            cur_tnodes = 10
            cur_success = 0
            cur_ct = 0
            for i, file in enumerate(graphs_files):
                print(f"{i} / {len(graphs_files)}")
                parameters = ["direct_transitivity" if form == "direct_transitivity" or form == "both_combined" else "",
                              "vertical_transitivity" if form == "vertical_transitivity" or form == "both_combined" else ""] + list(combo)
                folder = f"{form}/{folder_to_write}"
                if len(combo) == 0:
                    res = run_one_graph(file, f'{folder}/exp0', 60, parameters, i)
                else:
                    res = run_one_graph(file, f'{folder}/exp' + ''.join([str(ind+1) for ind, val in enumerate(key) if val in combo]), 60, parameters, i)
                if int(res[3]) >= cur_tnodes + 10:
                    if cur_success / cur_ct < 0.5:
                        print("-"*70, f"\n{form} experiment {combo} DONE\n", "-"*70)
                        break
                    cur_success = 0
                    cur_ct = 0
                    cur_tnodes += 10
                else:
                    if float(res[10]) < 60:
                        cur_success += 1
                    cur_ct += 1


def run_one_experiment(start_idx, graphs_file, exp_name, params_to_set, clear_files):
    with open(graphs_file, 'r') as f:
        gfiles = [gname.removesuffix('\n') for gname in f.readlines()]
    if clear_files:
        clear_file(f"{exp_name}.csv")
        insert_one(f"{exp_name}.csv", ["Index", "File", "X-vars", "C-vars", "Total vars", "Total constraints", "Crossings", "Opttime", "Nodes visited"])
    for i, to_opt in enumerate(gfiles[start_idx:]):
        print(f"{i+start_idx+1} / {len(gfiles)}")
        g = read_data.read(to_opt)
        params = {param: True for param in params_to_set}
        params.update({"cutoff_time": 300, "return_experiment_data": True})
        optimizer = optimization.LayeredOptimizer(g, params)
        result = optimizer.optimize_layout()
        insert_one(f"{exp_name}.csv", [i+start_idx+1, to_opt] + [j for j in result])


def run_one_graph(gfile, exp_name, cutoff_time, params_to_set, idx):
    g = read_data.read(gfile)
    base_info = basic_info(g)
    params = {param: True for param in params_to_set}
    params.update({"cutoff_time": cutoff_time, "return_experiment_data": True})
    optimizer = optimization.LayeredOptimizer(g, params)
    result = optimizer.optimize_layout()
    formatted = [idx, gfile] + base_info + [j for j in result]
    if int(formatted[11]) != 11:
        insert_one(f"{exp_name}.csv", formatted)
    return formatted


def sort_by_collection_and_tnodes(s):
    """ Sort primarily by Rome-Lib, then DAGmar, then North (AT&T), secondarily by total nodes """
    if s[1][0] == "R":
        add_val = -100000
    elif s[1][0] == "D":
        add_val = 0
    else:
        add_val = 100000
    return add_val + int(s[3])


# def bucket_lines_in_data(file, bucket_size):
#     """  """
#     lines_in_file = []
#     seen_files = set()
#     with open(file, 'r') as fd1:
#         rdr = csv.reader(fd1)
#         next(rdr)
#         for ln in rdr:
#             if bucket_size <= int(ln[3]) < bucket_size + 10 and ln[1] not in seen_files:
#                 lines_in_file.append(ln)
#                 seen_files.add(ln[1])
#     lines_in_file.sort(key=sort_by_collection_and_tnodes)
#     return lines_in_file


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


def get_all_files_by_bucket(n_g_per_bucket):
    with open(f"data storage/g{n_g_per_bucket}_sorted.txt", 'r') as fd1:
        filenames = []
        buckets = []
        for line in fd1.readlines():
            if line[0] == "T":
                filenames.append([])
                buckets.append(int(line[line.index('[')+1:line.index(',')]))
            else:
                filenames[-1].append(line[:line.index(',')])
    return filenames, buckets


def get_start_position(results_file, all_graphs, threshhold) -> tuple[int, int, int, int, bool]:
    x, y, cur_bucket, cur_success = 0, 0, 10, 0 
    with open(f"data storage/{results_file}", 'r') as fd:
        rdr = csv.reader(fd)
        for line in rdr:
            if line[1][0] == 'R' or line[1][0] == 'n':
                if int(line[3])//10*10 != cur_bucket:
                    cur_bucket = int(line[3])//10*10
                    x += 1
                    y = 0
                    cur_success = 0
                y += 1
                if int(line[11]) == 2:
                    cur_success += 1
    if y == len(all_graphs[x]) and cur_success / len(all_graphs[x]) < threshhold:
        return x + 1, 0, cur_success, x, True
    elif y == len(all_graphs[x]):
        return x + 1, 0, cur_success, x, False
    else:
        return x, y, cur_success, x, False


def individual_switch_cutoff(datapoints):
    timedout = sum((1 for pt in datapoints if float(pt[10]) > 60))
    return True if timedout / len(datapoints) >= 0.25 else False


def individual_switch_experiment(switch_num, transitivity_num, num_per_bucket):
    # This function is checkpoint-safe
    if "individual switch" not in os.listdir("data storage"):
        os.mkdir("data storage/individual switch")
    if "direct_transitivity" not in os.listdir("data storage/individual switch"):
        os.mkdir("data storage/individual switch/direct_transitivity")
    if "vertical_transitivity" not in os.listdir("data storage/individual switch"):
        os.mkdir("data storage/individual switch/vertical_transitivity")

    key1 = ["baseline", "symmetry_breaking", "butterfly_reduction", "mirror_vars", "cycle_constraints", "collapse_subgraphs", "heuristic_start", "priority", "mip_relax"]
    key2 = ["baseline_5m", "symmetry_breaking_5m", "butterfly_reduction_5m", "mirror_vars_5m", "cycle_constraints_5m", "collapse_subgraphs_5m", "heuristic_start_5m", "xvar_branch_priority_5m", "mip_relax_5m"]
    transitivity = "direct_transitivity" if transitivity_num == 0 else "vertical_transitivity"
    all_files, buckets = get_all_files_by_bucket(num_per_bucket)
    fname = f"individual switch/{transitivity}/{key2[switch_num]}"

    if switch_num != 0:
        if os.path.exists(f"data storage/{fname}.csv"):
            x, y, success_ct, furthest_bucket, is_complete = get_start_position(f"{fname}.csv", all_files, 0.75)
            if is_complete:
                print(f"{key1[switch_num]} with {transitivity} was cut off at bucket {furthest_bucket}")
                return
        else:
            insert_one(f"{fname}.csv", ["Index", "File", "Nodes", "Total Nodes", "Butterflies", "X-vars", "C-vars", "Total vars", "Total constraints", "Crossings", "Opttime", "Status", "Nodes visited", "Setup Time"])
            x, y, success_ct = 0, 0, 0
        success_condition = True
        parameters = [key1[switch_num], transitivity]
        while success_condition and x < len(all_files):
            for gidx in range(y, len(all_files[x])):
                print(f"\n{all_files[x][gidx]}  [bucket {buckets[x]}: {gidx}/{len(all_files[x])}]")
                res = run_one_graph(all_files[x][gidx], fname, 300, parameters, x * len(all_files[0]) + gidx)
                if int(res[11]) == 2:  # model status = 2 means solved to optimality w/o hitting cutoff
                    success_ct += 1
            if success_ct / len(all_files[x]) < 0.75:
                print(f"{transitivity} with switch {key1[switch_num]} cutoff at bucket size {buckets[x]}")
                success_condition = False
            else:
                x += 1
                y = 0
                success_ct = 0
    else:
        # baseline calculation
        furthest_reached = 0
        for key2v in key2[1:]:
            if os.path.exists(f"data storage/individual switch/{transitivity}/{key2v}.csv"):
                x, y, success_ct, furthest_bucket, is_complete = get_start_position(f"individual switch/{transitivity}/{key2v}.csv", all_files, 0.75)
                # if not is_complete:
                #     raise Exception("Wait until all other experiments are done to run the baseline")
                if furthest_bucket > furthest_reached:
                    furthest_reached = furthest_bucket
            # else:
            #     raise Exception("Wait until all other experiments are done to run the baseline")
        if os.path.exists(f"data storage/{fname}.csv"):
            x, y, dummy, dummy2, dummy3 = get_start_position(f"{fname}.csv", all_files, 0.75)
        else:
            insert_one(f"{fname}.csv", ["Index", "File", "Nodes", "Total Nodes", "Butterflies", "X-vars", "C-vars", "Total vars", "Total constraints", "Crossings", "Opttime", "Status", "Nodes visited", "Setup Time"])
            x, y = 0, 0
        parameters = [transitivity]
        while x <= furthest_reached:
            for gidx in range(y, len(all_files[x])):
                print(f"\n{all_files[x][gidx]}  [bucket {buckets[x]}: {gidx}/{len(all_files[x])}]")
                run_one_graph(all_files[x][gidx], fname, 300, parameters, x * len(all_files[0]) + gidx)
            x += 1
            y = 0


def all_combinations_experiment_checkpoint_safe(combination, num_per_bucket):
    if "all switches" not in os.listdir("data storage"):
        os.mkdir("data storage/all switches")
    if "direct_transitivity" not in os.listdir("data storage/all switches"):
        os.mkdir("data storage/all switches/direct_transitivity")
    if "vertical_transitivity" not in os.listdir("data storage/all switches"):
        os.mkdir("data storage/all switches/vertical_transitivity")

    key = ["symmetry_breaking", "butterfly_reduction", "mirror_vars", "cycle_constraints", "collapse_subgraphs", "heuristic_start", "priority", "mip_relax"]
    all_files, buckets = get_all_files_by_bucket(num_per_bucket)
    transitivity = "direct_transitivity" if combination[-1] else "vertical_transitivity"
    is_baseline = all((not j for j in combination[:-1]))
    fname = f"all switches/{transitivity}/exp_{''.join([str(j) for j, x in enumerate(combination) if x])}"

    if not is_baseline:
        if os.path.exists(f"data storage/{fname}.csv"):
            x, y, success_ct, furthest_bucket, is_complete = get_start_position(f"{fname}.csv", all_files, 0.5)
            if is_complete:
                print(f"Combination {combination} with {transitivity} was cut off at bucket {buckets[furthest_bucket]}")
                return
        else:
            insert_one(f"{fname}.csv", ["Index", "File", "Nodes", "Total Nodes", "Butterflies", "X-vars", "C-vars", "Total vars", "Total constraints", "Crossings", "Opttime", "Status", "Nodes visited", "Setup Time"])
            x, y, success_ct = 0, 0, 0
        success_condition = True
        parameters = [key[i] for i, val in enumerate(combination[:-1]) if val] + [transitivity]
        while success_condition and x < len(all_files):
            for gidx in range(y, len(all_files[x])):
                print(f"\n{all_files[x][gidx]}  [bucket {buckets[x]}: {gidx}/{len(all_files[x])}]")
                res = run_one_graph(all_files[x][gidx], fname, 60, parameters, x * len(all_files[0]) + gidx)
                if int(res[11]) == 2:  # model status = 2 means solved to optimality w/o hitting cutoff
                    success_ct += 1
            if success_ct / len(all_files[x]) < 0.5:
                print(f"{transitivity} with switches {', '.join([str(j) for j, x in enumerate(combination) if x])} cutoff at bucket size {buckets[x]}")
                success_condition = False
            else:
                x += 1
                y = 0
                success_ct = 0
    else:
        furthest_reached = 0
        for key2v in list(itertools.chain.from_iterable(itertools.combinations("12345678", r) for r in range(9))):
            key2 = ''.join(key2v)
            if os.path.exists(f"data storage/all switches/{transitivity}/exp_{key2}.csv"):
                x, y, success_ct, furthest_bucket, is_complete = get_start_position(f"all switches/{transitivity}/exp_{key2}.csv", all_files, 0.5)
                # if not is_complete:
                #     raise Exception("Wait until all other experiments are done to run the baseline")
                if furthest_bucket > furthest_reached:
                    furthest_reached = furthest_bucket
            # else:
            #     raise Exception("Wait until all other experiments are done to run the baseline")
        if os.path.exists(f"data storage/{fname}.csv"):
            x, y, dummy, dummy2, dummy3 = get_start_position(f"{fname}.csv", all_files, 0.5)
        else:
            insert_one(f"{fname}.csv", ["Index", "File", "Nodes", "Total Nodes", "Butterflies", "X-vars", "C-vars", "Total vars", "Total constraints", "Crossings", "Opttime", "Status", "Nodes visited", "Setup Time"])
            x, y = 0, 0
        parameters = [transitivity]
        while x <= furthest_reached:
            for gidx in range(y, len(all_files[x])):
                print(f"\n{all_files[x][gidx]}  [bucket {buckets[x]}: {gidx}/{len(all_files[x])}]")
                run_one_graph(all_files[x][gidx], fname, 60, parameters, x * len(all_files[0]) + gidx)
            x += 1
            y = 0


def sample_experiment_dataset(n_per_bucket):
    all_g = get_all_graphs()
    all_g_tnode = []
    for gr in all_g:
        my_g = read_data.read(gr)
        all_g_tnode.append((gr, len(my_g.nodes)))
    all_g_tnode.sort(key=lambda x: x[1])
    this_bucket = -10
    linestore = []
    bucket_sizes = []
    for gname in all_g_tnode:
        if gname[1] // 10 * 10 > this_bucket:
            bucket_sizes.append(gname[1] // 10 * 10)
            linestore.append([])
            this_bucket = gname[1] // 10 * 10
        linestore[-1].append(f"{gname[0]},{gname[1]}\n")
    for i, bucket in enumerate(linestore):
        if len(linestore[i]) > n_per_bucket:
            linestore[i] = random.sample(bucket, n_per_bucket)
    if not os.path.exists("data storage"):
        os.mkdir("data storage")
    with open(f"data storage/g{n_per_bucket}_sorted.txt", 'w') as fd1:
        for i, bucket in enumerate(linestore):
            fd1.write(f"Total nodes in [{bucket_sizes[i]},{bucket_sizes[i] + 10}):\n")
            for ln in bucket:
                fd1.write(ln)


# def sample_5_percent():
#     with open("data storage/5_percent_all_g_sorted.txt", 'w') as fd:
#         for i in range(10, 1110, 10):
#             with open("data storage/all_g_sorted.txt", 'r') as fd1:
#                 collect_lines = False
#                 files = []
#                 for line in fd1.readlines():
#                     if line[0] == "T":
#                         if int(line[line.index('[') + 1:line.index(',')]) == i:
#                             collect_lines = True
#                         else:
#                             collect_lines = False
#                     elif collect_lines:
#                         files.append(line)
#             if len(files) >= 5:
#                 new_files = random.sample(files, round(len(files)/20))
#                 for file in new_files:
#                     fd.write(file)


def count_complete(num_per_bucket):
    all_files, buckets = get_all_files_by_bucket(num_per_bucket)
    complete = 0
    incomplete = 0
    for file in os.listdir("data storage/all switches/vertical_transitivity"):
        x, y, success_ct, furthest_bucket, is_complete = get_start_position(f"/data storage/all switches/vertical_transitivity/{file}", all_files, 0.5)
        if is_complete:
            complete += 1
        else:
            incomplete += 1
    print(f"{complete} complete, {incomplete} incomplete")


if __name__ == '__main__':
    random.seed(22)
    num_g_per_bucket = 200
    if len(sys.argv) >= 2:
        exp_choice = int(sys.argv[1])
        switch_to_test = int(sys.argv[2])
    else:
        exp_choice = 2
        switch_to_test = 1

    """ Sample all data to generate experiment dataset, as described in our paper. Generates ./data storage/g[n per bucket]_sorted.txt """
    if exp_choice == 1 and not os.path.exists(f"data storage/g{num_g_per_bucket}_sorted.txt"):
        sample_experiment_dataset(num_g_per_bucket)
    elif exp_choice == 2 and not os.path.exists(f"data storage/g{num_g_per_bucket // 2}_sorted.txt"):
        sample_experiment_dataset(num_g_per_bucket // 2)

    """ Individual switch evaluation experiment """
    if exp_choice == 1:
        individual_switch_experiment(switch_to_test // 2, switch_to_test % 2, num_g_per_bucket)  # performs the individual switch experiment, writing all data to csv files in ./data storage/individual switch

    """ All combinations of switches/formulations experiment """
    if exp_choice == 2:
        combo_to_test = [bool(int(j)) for j in f'{switch_to_test:09b}']
        all_combinations_experiment_checkpoint_safe(combo_to_test, num_g_per_bucket // 2)  # performs all combinations experiment, writing all data to csv files in ./data storage/all switches
        # Files are named by switches used according to this key:
        # 1=symmetry breaking, 2=butterfly reduction, 3=mirrored vars, 4=heuristic start, 5=x-var priority, 6=mip relax, 7=cycle constraints, 8=collapse subgraphs
