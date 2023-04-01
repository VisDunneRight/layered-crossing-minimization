import os
import itertools
import pickle
import csv
import random
from src import optimization, read_data, vis, motifs


random.seed(22)


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
        optimizer = optimization.LayeredOptimizer(g, {"name": to_opt, "cutoff_time": 600, "fix_one_var": True, "return_experiment_data": True, "stratisfimal_yvars": True})
        result = optimizer.optimize_layout()
        insert_data("strat_fix1.csv", [result])
        optimizer.junger_ec, optimizer.stratisfimal_y_vars = True, False
        result = optimizer.optimize_layout()
        insert_data("junger_fix1.csv", [result])
        optimizer.mirror_vars, optimizer.junger_ec = True, False
        result = optimizer.optimize_layout()
        insert_data("sankey_fix1.csv", [result])


def get_all_graphs():
    all_g = [[], [], []]
    for i in range(10, 101):
        for file in os.listdir(f"Rome-Lib/graficon{i}nodi"):
            all_g[0].append(f"Rome-Lib/graficon{i}nodi/" + file)
    daglist = []
    for i in range(1, 8):
        for file in os.listdir(f"DAGmar/graphs/{i}.6"):
            daglist.append(f"DAGmar/graphs/{i}.6/" + file)
    g_t_nds = {}
    for file in daglist:
        g = read_data.read(file)
        g_t_nds[file] = len(g.nodes)
    all_g[1].extend(sorted(daglist, key=lambda x: g_t_nds[x]))
    north_gs = sorted(list(os.listdir("north")), key=lambda fil: int(fil[2:(5 if fil[4] == '0' else 4)]))
    for i in range(len(north_gs)):
        north_gs[i] = "north/" + north_gs[i]
    north_gs.remove("north/g.57.26.graphml")   # skip this one graph that takes over an hour to insert the variables and constraints
    all_g[2].extend(north_gs)
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
        insert_one(f"junger_basic/{exp_name}_{cutoff_time}.csv", ["Index", "File", "Nodes", "Total Nodes", "Butterflies", "X-vars", "C-vars", "Total vars", "Total constraints", "Crossings", "Opttime", "Work", "Nodes visited", "Setup Time"])
        insert_one(f"vertical_transitivity/{exp_name}_{cutoff_time}.csv", ["Index", "File", "Nodes", "Total Nodes", "Butterflies", "X-vars", "C-vars", "Total vars", "Total constraints", "Crossings", "Opttime", "Work", "Nodes visited", "Setup Time"])
        insert_one(f"redundancy/{exp_name}_{cutoff_time}.csv", ["Index", "File", "Nodes", "Total Nodes", "Butterflies", "X-vars", "C-vars", "Total vars", "Total constraints", "Crossings", "Opttime", "Work", "Nodes visited", "Setup Time"])
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
        insert_one(f"junger_basic/{exp_name}_{cutoff_time}.csv", ["Index", "File", "Nodes", "Total Nodes", "Butterflies", "X-vars", "C-vars", "Total vars", "Total constraints", "Crossings", "Opttime", "Work", "Nodes visited", "Setup Time"])
        insert_one(f"vertical_transitivity/{exp_name}_{cutoff_time}.csv", ["Index", "File", "Nodes", "Total Nodes", "Butterflies", "X-vars", "C-vars", "Total vars", "Total constraints", "Crossings", "Opttime", "Work", "Nodes visited", "Setup Time"])
        insert_one(f"redundancy/{exp_name}_{cutoff_time}.csv", ["Index", "File", "Nodes", "Total Nodes", "Butterflies", "X-vars", "C-vars", "Total vars", "Total constraints", "Crossings", "Opttime", "Work", "Nodes visited", "Setup Time"])
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
    key = ["fix_one_var", "butterfly_reduction", "heuristic_start", "presolve", "priority", "mip_relax", "mirror_vars"]
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
    # folder = "redundancy" if "direct_transitivity" in params_to_set and "vertical_transitivity" in params_to_set else "junger_basic" if "direct_transitivity" in params_to_set else "vertical_transitivity"
    formatted = [idx, gfile] + base_info + [j for j in result]
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


def individual_switch_cutoff(datapoints):
    timedout = sum((1 for pt in datapoints if float(pt[10]) > 60))
    return True if timedout / len(datapoints) >= 0.25 else False


def individual_switch_experiment():
    key1 = ["fix_one_var", "butterfly_reduction", "heuristic_start", "presolve", "priority", "mip_relax", "mirror_vars"]
    key2 = ["fix1var_60", "butterfly_60", "heuristic_60", "presolve_60", "xvar_branch_60", "mip_relax_60", "symmetry_60"]
    for j, inp1 in enumerate(["direct_transitivity", "vertical_transitivity", "both_combined"]):
        furthest_bucket_reached = 0
        for i, inp2 in enumerate(key2):
            fname = f"{inp1}/{inp2}"
            insert_one(f"{fname}_new.csv", ["Index", "File", "Nodes", "Total Nodes", "Butterflies", "X-vars", "C-vars", "Total vars", "Total constraints", "Crossings", "Opttime", "Work", "Nodes visited", "Setup Time"])
            curindex = 0
            for bsize in range(10, 17141, 10):
                all_bfiles = get_all_files_in_bucket(bsize)
                success_ct = 0
                if len(all_bfiles) > 0:
                    for check_file in all_bfiles:
                        parameters = [key1[i], "direct_transitivity" if j % 2 == 0 else "baseline",
                                      "vertical_transitivity" if j > 0 else "baseline"]
                        res = run_one_graph(check_file, f"{fname[fname.index('/') + 1:]}_new", 60, parameters, curindex)
                        curindex += 1
                        if res[10] < 60:
                            success_ct += 1
                if success_ct / len(all_bfiles) < 0.75:
                    print(f"{inp1} with switch {inp2} cutoff at bucket size {bsize}")
                    if bsize > furthest_bucket_reached:
                        furthest_bucket_reached = bsize
                    break
        # baseline calculation
        fname = f"{inp1}/baseline"
        insert_one(f"{fname}_new.csv",
                   ["Index", "File", "Nodes", "Total Nodes", "Butterflies", "X-vars", "C-vars", "Total vars",
                    "Total constraints", "Crossings", "Opttime", "Work", "Nodes visited", "Setup Time"])
        curindex = 0
        for bsize in range(10, furthest_bucket_reached + 10, 10):
            all_bfiles = get_all_files_in_bucket(bsize)
            if len(all_bfiles) > 0:
                for check_file in all_bfiles:
                    parameters = ["direct_transitivity" if j % 2 == 0 else "baseline", "vertical_transitivity" if j > 0 else "baseline"]
                    run_one_graph(check_file, f"{fname[fname.index('/') + 1:]}_new", 60, parameters, curindex)
                    curindex += 1


def sample_experiment_dataset():
    all_g = get_all_graphs()
    rome_g = []
    for gr in all_g[0]:
        my_g = read_data.read(gr)
        rome_g.append((gr, len(my_g.nodes)))
    # rome_g.sort(key=lambda x: x[1])
    print(len(rome_g))
    rome_g_sample = random.sample(rome_g, len(rome_g)//2)
    rome_g_sample.sort(key=lambda x: x[1])
    print(len(rome_g_sample))
    with open("data storage/rome_sorted.txt", 'w') as fd:
        fd.write("Total nodes in [10,20):\n")
        for i in range(len(rome_g_sample)):
            if i > 0 and rome_g_sample[i][1]//10 > rome_g_sample[i-1][1]//10:
                fd.write(f"Total nodes in [{rome_g_sample[i][1]//10*10},{rome_g_sample[i][1]//10*10+10}):\n")
            fd.write(f"{rome_g_sample[i][0]},{rome_g_sample[i][1]}\n")
    dagmar_g = []
    for gr in all_g[1]:
        my_g = read_data.read(gr)
        dagmar_g.append((gr, len(my_g.nodes)))
    dagmar_g.sort(key=lambda x: x[1])
    with open("data storage/dagmar_sorted.txt", 'w') as fd:
        for i in range(len(dagmar_g)):
            if i > 0 and dagmar_g[i][1] // 10 > dagmar_g[i - 1][1] // 10:
                fd.write(f"Total nodes in [{dagmar_g[i][1] // 10 * 10},{dagmar_g[i][1] // 10 * 10 + 10}):\n")
            fd.write(f"{dagmar_g[i][0]},{dagmar_g[i][1]}\n")
    north_g = []
    for gr in all_g[2]:
        my_g = read_data.read(gr)
        north_g.append((gr, len(my_g.nodes)))
    north_g.sort(key=lambda x: x[1])
    with open("data storage/north_sorted.txt", 'w') as fd:
        fd.write("Total nodes in [0,10):\n")
        for i in range(len(north_g)):
            if i > 0 and north_g[i][1] // 10 > north_g[i - 1][1] // 10:
                fd.write(f"Total nodes in [{north_g[i][1] // 10 * 10},{north_g[i][1] // 10 * 10 + 10}):\n")
            fd.write(f"{north_g[i][0]},{north_g[i][1]}\n")
    with open("data storage/all_g_sorted.txt", 'w') as fd1:
        tval_hash = {}  # tnode count -> index in linestore
        n_seen = 0
        linestore = []
        for fname in ["rome_sorted", "dagmar_sorted", "north_sorted"]:
            with open(f"data storage/{fname}.txt", 'r') as fd2:
                cur_tnodes = 0
                for line in fd2.readlines():
                    if line[0] == "T":
                        cur_tnodes = int(line[line.index('[') + 1:line.index(',')])
                        if cur_tnodes not in tval_hash:
                            tval_hash[int(line[line.index('[') + 1:line.index(',')])] = n_seen
                            n_seen += 1
                            linestore.append([])
                    else:
                        linestore[tval_hash[cur_tnodes]].append(line)
        for i in sorted(list(tval_hash.keys())):
            cur_idx = next((tval_hash[j] for j in tval_hash if j == i))  # index in linestore for next tnode val
            fd1.write(f"Total nodes in [{i},{i + 10}):\n")
            for ln in linestore[cur_idx]:
                fd1.write(ln)


def sample_5_percent():
    with open("data storage/5_percent_all_g_sorted.txt", 'w') as fd:
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


if __name__ == '__main__':
    """ NOTE: Running this file will take a very long time, at minimum 2-3 weeks. """

    """ Individual switch evaluation experiment """
    sample_experiment_dataset()  # sample all data to generate experiment dataset, as described in our paper. Generates data storage/all_g_sorted.txt
    individual_switch_experiment()  # performs the individual switch experiment, writing all data to csv files in /data storage/{formulation}

    """ All combinations of switches/formulations experiment """
    sample_5_percent()  # sample 5% of the experiment dataset (per 10-node bin), write to data storage/5_percent_all_g_sorted.txt
    all_combinations_experiment("all_combos_5percent")  # performs all combinations experiment, writing all data to csv files in /data storage/{formulation}/all_combos_5percent
    # Files are named by switches used according to this key: 0=baseline, 1=fix one var, 2=butterfly reduction, 3=heuristic start, 4=presolve, 5=xvar priority, 6=mip relax, 7=mirror vars

