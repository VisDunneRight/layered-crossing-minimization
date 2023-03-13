import os
import pickle
import csv
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
    for i in range(1, 11):
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
        clear_file(f"strat_big_m/{exp_name}_{cutoff_time}.csv")
        clear_file(f"redundancy/{exp_name}_{cutoff_time}.csv")
        insert_one(f"junger_basic/{exp_name}_{cutoff_time}.csv", ["Index", "File", "Nodes", "Total Nodes", "Butterflies", "X-vars", "C-vars", "Total vars", "Total constraints", "Crossings", "Opttime", "Work", "Nodes visited", "Setup Time"])
        insert_one(f"strat_big_m/{exp_name}_{cutoff_time}.csv", ["Index", "File", "Nodes", "Total Nodes", "Butterflies", "X-vars", "C-vars", "Total vars", "Total constraints", "Crossings", "Opttime", "Work", "Nodes visited", "Setup Time"])
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
            optimizer = optimization.LayeredOptimizer(g, {"cutoff_time": cutoff_time, "return_experiment_data": True, "junger_trans": True, param_to_set: True})
            result = optimizer.optimize_layout()
            insert_one(f"junger_basic/{exp_name}_{cutoff_time}.csv", [i+start_idx[1], to_opt] + base_info + [j for j in result])
            if result[5] >= cutoff_time or junger_timedout >= max_timeout:
                junger_timedout += 1
            else:
                junger_timedout = 0
            optimizer.junger_trans, optimizer.strat_big_m = False, True
            result = optimizer.optimize_layout()
            insert_one(f"strat_big_m/{exp_name}_{cutoff_time}.csv", [i+start_idx[1], to_opt] + base_info + [j for j in result])
            if result[5] >= cutoff_time or strat_timedout >= max_timeout:
                strat_timedout += 1
            else:
                strat_timedout = 0
            optimizer.junger_trans = True
            result = optimizer.optimize_layout()
            insert_one(f"redundancy/{exp_name}_{cutoff_time}.csv", [i+start_idx[1], to_opt] + base_info + [j for j in result])
            if result[5] >= cutoff_time or redundant_timedout >= max_timeout:
                redundant_timedout += 1
            else:
                redundant_timedout = 0


def run_multi_param_experiment(start_idx, graphs_file, cutoff_time, exp_name, params_to_set, clear_files, max_timeout):
    with open(graphs_file, 'r') as f:
        gfiles = [gname.removesuffix('\n') for gname in f.readlines()]
    if clear_files:
        if start_idx != (0, 0):
            print("Something's wrong here...")
            return
        clear_file(f"junger_basic/{exp_name}_{cutoff_time}.csv")
        clear_file(f"strat_big_m/{exp_name}_{cutoff_time}.csv")
        clear_file(f"redundancy/{exp_name}_{cutoff_time}.csv")
        insert_one(f"junger_basic/{exp_name}_{cutoff_time}.csv", ["Index", "File", "Nodes", "Total Nodes", "Butterflies", "X-vars", "C-vars", "Total vars", "Total constraints", "Crossings", "Opttime", "Work", "Nodes visited", "Setup Time"])
        insert_one(f"strat_big_m/{exp_name}_{cutoff_time}.csv", ["Index", "File", "Nodes", "Total Nodes", "Butterflies", "X-vars", "C-vars", "Total vars", "Total constraints", "Crossings", "Opttime", "Work", "Nodes visited", "Setup Time"])
        insert_one(f"redundancy/{exp_name}_{cutoff_time}.csv", ["Index", "File", "Nodes", "Total Nodes", "Butterflies", "X-vars", "C-vars", "Total vars", "Total constraints", "Crossings", "Opttime", "Work", "Nodes visited", "Setup Time"])
    for j in range(start_idx[0], 3):
        junger_timedout = 0
        strat_timedout = 0
        redundant_timedout = 0
        for i, to_opt in enumerate(gfiles[j][start_idx[1]:]):
            if junger_timedout >= max_timeout and strat_timedout >= max_timeout and redundant_timedout >= max_timeout:
                break
            print(f"{i + start_idx[1] + 1} / {len(gfiles[j])}")
            g = read_data.read(to_opt)
            base_info = basic_info(g)
            params = {param: True for param in params_to_set}
            params.update({"cutoff_time": cutoff_time, "return_experiment_data": True, "junger_trans": True})
            optimizer = optimization.LayeredOptimizer(g, params)
            result = optimizer.optimize_layout()
            insert_one(f"junger_basic/{exp_name}_{cutoff_time}.csv",
                       [i + start_idx[1], to_opt] + base_info + [j for j in result])
            if result[5] >= cutoff_time or junger_timedout >= max_timeout:
                junger_timedout += 1
            else:
                junger_timedout = 0
            optimizer.junger_trans, optimizer.strat_big_m = False, True
            result = optimizer.optimize_layout()
            insert_one(f"strat_big_m/{exp_name}_{cutoff_time}.csv",
                       [i + start_idx[1], to_opt] + base_info + [j for j in result])
            if result[5] >= cutoff_time or strat_timedout >= max_timeout:
                strat_timedout += 1
            else:
                strat_timedout = 0
            optimizer.junger_trans = True
            result = optimizer.optimize_layout()
            insert_one(f"redundancy/{exp_name}_{cutoff_time}.csv",
                       [i + start_idx[1], to_opt] + base_info + [j for j in result])
            if result[5] >= cutoff_time or redundant_timedout >= max_timeout:
                redundant_timedout += 1
            else:
                redundant_timedout = 0


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
