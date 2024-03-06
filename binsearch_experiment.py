import os.path
import random
import sys
import csv
from src.optimization import LayeredOptimizer
from src import heuristics, vis
from src.neighborhood import *


def insert_one(filename, entry):
    with open(filename, 'a', newline='') as f:
        wrt = csv.writer(f)
        wrt.writerow(entry)


def get_starting_bounds(file_path, graph_size: str, target_avg: float):
    with open(file_path, 'r') as fd:
        rdr = csv.reader(fd)
        cur_bnds = [0, sys.maxsize]
        cur_avgs = [sys.maxsize, -1]
        n_times_searched = 0
        for ln in rdr:
            if ln[0] == graph_size:
                if int(ln[1]) > cur_bnds[0] and float(ln[2]) >= target_avg:
                    cur_bnds[0] = int(ln[1])
                    cur_avgs[0] = float(ln[2])
                elif int(ln[1]) < cur_bnds[1] and float(ln[2]) <= target_avg:
                    cur_bnds[1] = int(ln[1])
                    cur_avgs[1] = float(ln[2])
                n_times_searched += 1
        if cur_bnds[1] == sys.maxsize:
            cur_bnds[1] = 2 * cur_bnds[0]
    return cur_bnds, cur_avgs, n_times_searched // 3


def get_start_position_binsearch(filename, cur_cv):
    if os.path.exists(filename):
        fset = {}
        with open(filename, 'r') as fd:
            rdr = csv.reader(fd)
            next(rdr)
            for ln in rdr:
                if ln[1] in fset:
                    fset.clear()
                if ln[2] == str(cur_cv):
                    fset[ln[1]] = (len(ln) - 7) / (2 * float(ln[3])) * 60
        return fset
    else:
        return {}


def run_experiment(neighborhood_fn, target_avg: int, graph_size: str, initial_layout_fn, path_to_dataset: str, n_graph_copies: int, depth=10, time_per_graph=60):
    """
    :param neighborhood_fn: neighborhood aggregation function, e.g. bfs_neighborhood from src/neighborhood.py
    :param target_avg: target num of opts/5mins
    :param graph_size: string corresponding to graph set, e.g. 'r1.5k12n8'
    :param initial_layout_fn: e.g. barycenter from src/heuristics.py
    :param path_to_dataset: path to dataset, full path will be path_to_dataset/graph_size
    :param n_graph_copies: num graphs in above directory, will optimize graph[0...n].lgbin
    :param depth: binary search depth to stop at
    :param time_per_graph: time spent optimizing each graph (seconds)
    :return: runs binary search to zero in on relative val for size calculation on this neighborhood, graph size
    """

    nbhd = neighborhood_fn.__name__.replace('_neighborhood', '')
    fname = f"{path_to_dataset}/bounds_results/{nbhd}+{graph_size}+{str(target_avg)}.csv"
    binfname = f"{path_to_dataset}/bounds_results/{nbhd}_bounds.csv"
    starting_bounds, starting_avgs, n_searches = get_starting_bounds(binfname, graph_size, target_avg)
    files_run = get_start_position_binsearch(fname, sum(starting_bounds)//2)
    if len(files_run) == 0 and (not os.path.isfile(fname) or os.path.getsize(fname) == 0):
        insert_one(fname, ["Index", "File", "CVcalc", "OptTime", "CrFinal", "Cr1", "T1", "Cr2", "T2..."])
    if len(files_run) == n_graph_copies:
        files_run.clear()
    cur_idx = 0
    while n_searches < depth:
        if n_searches == 0:
            cv = random.randint(1000, 3000)
            for i in range(n_graph_copies):
                fl = f"graph{i}.lgbin"
                print(fl)
                print(files_run)
                if f"{path_to_dataset}/{graph_size}/{fl}" not in files_run:
                    optim = LayeredOptimizer(f"{path_to_dataset}/{graph_size}/{fl}", cutoff_time=time_per_graph)
                    initial_layout_fn(optim.g)
                    output = optim.local_opt_increment(cv, neighborhood_fn=neighborhood_fn, candidate_fn=random_candidate)
                    reordered = [v for i in range(len(output[2])) for v in (output[2][i], output[3][i])]
                    files_run[f"{path_to_dataset}/{graph_size}/{fl}"] = (len(reordered) - 2) / (2 * output[0]) * 60
                    insert_one(fname, [cur_idx, f"{path_to_dataset}/{graph_size}/{fl}", cv, output[0], output[1]] + reordered)
                cur_idx += 1
            avg_opts = sorted(list(files_run.values()))[n_graph_copies // 2]
            insert_one(binfname, [graph_size, cv, avg_opts])
            starting_bounds, starting_avgs, n_searches = get_starting_bounds(binfname, graph_size, target_avg)
            n_searches -= 1
        else:
            print(starting_bounds, starting_avgs, n_searches)
            for i in range(n_graph_copies):
                fl = f"graph{i}.lgbin"
                print(fl)
                print(files_run)
                if f"{path_to_dataset}/{graph_size}/{fl}" not in files_run:
                    optim = LayeredOptimizer(f"{path_to_dataset}/{graph_size}/{fl}", cutoff_time=time_per_graph)
                    initial_layout_fn(optim.g)
                    output = optim.local_opt_increment(sum(starting_bounds) // 2, neighborhood_fn=neighborhood_fn, candidate_fn=random_candidate)
                    reordered = [v for i in range(len(output[2])) for v in (output[2][i], output[3][i])]
                    files_run[f"{path_to_dataset}/{graph_size}/{fl}"] = (len(reordered) - 2) / (2 * output[0]) * 60
                    insert_one(fname, [cur_idx, f"{path_to_dataset}/{graph_size}/{fl}", sum(starting_bounds) // 2, output[0], output[1]] + reordered)
                cur_idx += 1
            if len(files_run) == n_graph_copies:
                # avg_opts = sum(files_run.values()) / n_graph_copies
                avg_opts = sorted(list(files_run.values()))[n_graph_copies // 2]
                insert_one(binfname, [graph_size, sum(starting_bounds) // 2, avg_opts])
                if starting_avgs[1] == -1 or avg_opts < starting_avgs[1]:  # top priority is finding valid upper bound
                    if avg_opts < target_avg:  # upper bound found
                        starting_bounds[1] = sum(starting_bounds) // 2
                        starting_avgs[1] = avg_opts
                    else:  # need to set new upper bound target to enlarge search
                        starting_bounds[0] = sum(starting_bounds) // 2
                        starting_avgs[0] = avg_opts
                        starting_bounds[1] *= 2
                elif avg_opts >= starting_avgs[1]:  # don't care to separate when cur > avgs[0] since lower bounded by 0
                    if avg_opts > target_avg:  # cur_avg > target > st_avg[1] but bounds[0] < midpt < bounds[1]
                        starting_bounds[0] = sum(starting_bounds) // 2
                        starting_avgs[0] = avg_opts
                    else:
                        starting_bounds[1] = sum(starting_bounds) // 2
                        starting_avgs[1] = avg_opts
        n_searches += 1
        files_run.clear()
        cur_idx = 0


if __name__ == '__main__':
    n_graphs_in_bin = 50
    opts_targets = [10, 50, 100]
    nbhd_fns = [bfs_neighborhood, vertical_re_neighborhood, degree_ratio_neighborhood, random_neighborhood]
    graph_sizes = ["r1.5k18n12", "r1.5k24n16", "r1.5k30n20", "r1.5k36n24", "r1.5k42n28"]
    if len(sys.argv) >= 2:
        target_idx = int(sys.argv[1]) // (len(nbhd_fns) * len(graph_sizes))
        nbhd_idx = (int(sys.argv[1]) % (len(nbhd_fns) * len(graph_sizes))) % len(nbhd_fns)
        graph_idx = (int(sys.argv[1]) % (len(nbhd_fns) * len(graph_sizes))) // len(nbhd_fns)
    else:
        target_idx, nbhd_idx, graph_idx = 0, 0, 0
    run_experiment(nbhd_fns[nbhd_idx], opts_targets[target_idx], graph_sizes[graph_idx], heuristics.barycenter, "random graphs/ratio_d3", n_graphs_in_bin, time_per_graph=60)
