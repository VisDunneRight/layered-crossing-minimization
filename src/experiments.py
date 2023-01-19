import os
import pickle
import csv
from src import optimization, read_data


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
            # f.write(','.join(str(e) for e in entry) + '\n')


def read_data_from_file(name):
    data = []
    with open("data storage/" + name, 'r') as f:
        l1 = f.readline().removesuffix('\n')
        for val in l1.split(' '):
            data.append([int(val) if val.isnumeric() else float(val) if not val.isalpha() else val])
        for entry in f.readlines():
            for i, val in enumerate(entry.removesuffix('\n').split(' ')):
                data[i].append(int(val) if val.isnumeric() else float(val) if not val.isalpha() else val)
    return data


def baseline_experiment(start_idx, list_of_files):
    pass


def independent_var_experiment(file_name):
    results = []
    for i, file in enumerate(get_list_of_files(file_name)):
        result = [i, file]
        g = read_data.read(file)
        result.extend((sum(1 for n in g.nodes if not n.is_anchor_node), len(g.nodes), len(g.edges)))
        opt = optimization.LayeredOptimizer(g, {"return_full_data": True})
        result.extend(opt.optimize_layout())
        results.append(result)
        if i % 10 == 0:
            insert_data("independent_var_study.csv", results)
            results.clear()


def fix_1_var_experiment(start_idx, list_of_files):
    n_nodes = 67
    to_optimize = os.listdir(f"Rome-Lib/graficon{n_nodes}nodi")[:10]
    times1 = []
    optvals1 = []
    times2 = []
    optvals2 = []
    e_by_n = []
    for to_opt in list_of_files:
        g = read_data.read(to_opt)
        e_by_n.append(len(g.edges) / len(g.nodes))
        optimizer = optimization.LayeredOptimizer(g, {"name": to_opt, "butterfly_reduction": False, "verbose": False, "cutoff_time": 100, "fix_one_var": True})
        a, b = optimizer.optimize_layout()
        times1.append(a)
        optvals1.append(b)
        optimizer.fix_one_var = False
        a, b = optimizer.optimize_layout()
        times2.append(a)
        optvals2.append(b)
    values = []
    for i, v in enumerate(times1):
        values.append({'#edges/#nodes': e_by_n[i], 'Runtime': v, 'Fix 1 Var': 'Yes'})
        values.append({'#edges/#nodes': e_by_n[i], 'Runtime': times2[i], 'Fix 1 Var': 'No'})
    # data = alt.Data(values=values)
    # chart = alt.Chart(data).mark_circle(size=60).encode(
    #     x='#edges/#nodes:Q',
    #     y=alt.Y('Runtime:Q', scale=alt.Scale(type='log')),
    #     color=alt.Color('Fix 1 Var:N', scale=alt.Scale(scheme='dark2'))
    # )
    # save(chart, "testing1234.svg")
