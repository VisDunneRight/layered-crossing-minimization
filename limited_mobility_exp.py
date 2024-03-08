import os.path
import sys
import csv
from src.optimization import LayeredOptimizer
from src import heuristics
from src.neighborhood import *
from os import mkdir, listdir


def insert_one(filename, entry):
    with open(filename, 'a', newline='') as f:
        wrt = csv.writer(f)
        wrt.writerow(entry)


def run_one(filename, gfile, fn, index):
    g = LayeredOptimizer(gfile)
    output = fn(g)
    insert_one(filename, [index, gfile] + output)


def get_start_position(filename):
    if os.path.exists(filename):
        fset = set()
        with open(filename, 'r') as fd:
            rdr = csv.reader(fd)
            for ln in rdr:
                fset.add(ln[1])
        return fset
    else:
        return set()


def get_csv_lines(file, exp):
    with open(file, 'r') as fd:
        rdr = csv.reader(fd)
        return [row for idx, row in enumerate(rdr) if exp in row[1]]


def get_closest_cv(file_path, graph_size: str, target_avg: float):
    with open(file_path, 'r') as fd:
        rdr = csv.reader(fd)
        best_bnd = 0
        best_avg = sys.maxsize
        for ln in rdr:
            if ln[0] == graph_size:
                if abs(float(ln[2]) - target_avg) < abs(best_avg - target_avg):
                    best_bnd = int(ln[1])
                    best_avg = float(ln[2])
    return best_bnd


def run_experiment(neighborhood_fn, nbhd_size, initial_layout_fn, path_to_dataset, subdirs, num_graphs):
    if "results" not in listdir(path_to_dataset):
        mkdir(path_to_dataset + "/results")
    # if neighborhood_fn.__name__ not in listdir(path_to_dataset + "/results"):
    #     (path_to_dataset + "/results/" + neighborhood_fn.__name__.replace("_neighborhood", "") + "+" + candidate_fn.__name__.replace("_candidate", ""))
    nbhd_name = neighborhood_fn.__name__.replace("_neighborhood", "")
    fname = path_to_dataset + "/results/" + nbhd_name + "+" + str(nbhd_size) + ".csv"
    files_run = get_start_position(fname)
    if len(files_run) == 0 and (not os.path.isfile(fname) or os.path.getsize(fname) == 0):
        insert_one(fname, ["Index", "File", "OptTime", "CrFinal", "Cr1", "T1", "Cr2", "T2..."])
    cur_idx = 0
    for subdir in subdirs:
        n_cvs = get_closest_cv(f"{path_to_dataset}/bounds_results/{nbhd_name}_bounds.csv", subdir, nbhd_size)
        for fl_num in range(num_graphs):
            fl = f"graph{fl_num}.lgbin"
            file_path = f"{path_to_dataset}/{subdir}/{fl}"
            print(file_path)
            if file_path not in files_run:
                optim = LayeredOptimizer(file_path, cutoff_time=300)
                initial_layout_fn(optim.g)
                output = optim.local_opt_increment(n_cvs, neighborhood_fn=neighborhood_fn, candidate_fn=random_candidate)
                reordered = [v for i in range(len(output[2])) for v in (output[2][i], output[3][i])]
                insert_one(fname, [cur_idx, file_path, output[0], output[1]] + reordered)
            cur_idx += 1


if __name__ == '__main__':
    # sandbox()

    dataset_path = "random graphs/ratio_d3"
    subdirectories = ["r1.5k18n12", "r1.5k24n16", "r1.5k30n20", "r1.5k36n24", "r1.5k42n28"]
    num_graphs_in_subdir = 50
    nbhd_sizes = [10, 50, 100]
    nbhd_fns = [bfs_neighborhood, vertical_re_neighborhood, degree_ratio_neighborhood, random_neighborhood]
    graph_types = ["triangle", "big_layer"]
    if len(sys.argv) >= 2:
        nbhd_idx = int(sys.argv[1]) // (len(graph_types) * len(nbhd_sizes))
        type_idx = (int(sys.argv[1]) % (len(graph_types) * len(nbhd_sizes))) % len(graph_types)
        cv_idx = (int(sys.argv[1]) % (len(graph_types) * len(nbhd_sizes))) // len(graph_types)
    else:
        nbhd_idx, type_idx, cv_idx = 0, 1, 1
    run_experiment(nbhd_fns[nbhd_idx], nbhd_sizes[cv_idx], heuristics.barycenter, f"random graphs/{graph_types[type_idx]}", subdirectories, num_graphs_in_subdir)
