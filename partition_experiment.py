import os.path
import sys
import csv
from src.optimization import LayeredOptimizer
from src.helpers import *
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
        with open(filename, 'r') as fd:
            lines = fd.readlines()
            last = lines[-1].split(',')
            if last[0] == "Index":
                return 0
            else:
                return int(last[0]) + 1
    else:
        return -1


def small_test():
    opt = LayeredOptimizer("random graphs/n_by_n/n25/graph0.lgbin")
    opt.cutoff_time = 300
    opt.vertical_transitivity = True
    # opt = LayeredOptimizer("Rome-Lib/graficon96nodi/grafo3510.96")
    n_cv = opt.g.c_vars_count()
    n_p = 1
    heuristics.global_sifting(opt.g)
    print(n_cv)
    print(optimization_time_estimate(n_cv / n_p))
    while n_p * optimization_time_estimate(n_cv / n_p) > 60 and n_p < 100:
        n_p += 1
    print(n_p, optimization_time_estimate(n_cv / n_p), n_cv / n_p)
    if n_p > 1:
        ub_n_cv = n_cv / (n_p - 1)
        lb_n_cv = n_cv / n_p
        while abs(ub_n_cv - lb_n_cv) > 1:
            mid = (ub_n_cv + lb_n_cv) / 2
            if calc_time_taken_for_partition_size(n_p - 1, mid, n_cv) < 60:
                lb_n_cv = mid
            else:
                ub_n_cv = mid
        print(lb_n_cv)
        print(calc_time_taken_for_partition_size(n_p - 1, lb_n_cv, n_cv))
        print("individual partition", optimization_time_estimate(lb_n_cv))
        opt.optimize_layout(bucket_size=10000, local_opt=True, pct=1)
    # heuristics.global_sifting(opt.g)
    print(opt.g.num_edge_crossings())


def run_experiment(neighborhood_fn, candidate_fn, n_cvs, initial_layout_fn, path_to_dataset):
    if "results" not in listdir(path_to_dataset):
        mkdir(path_to_dataset + "/results")
    # if neighborhood_fn.__name__ not in listdir(path_to_dataset + "/results"):
    #     (path_to_dataset + "/results/" + neighborhood_fn.__name__.replace("_neighborhood", "") + "+" + candidate_fn.__name__.replace("_candidate", ""))
    fname = path_to_dataset + "/results/" + neighborhood_fn.__name__.replace("_neighborhood", "") + "+" + candidate_fn.__name__.replace("_candidate", "") + ".csv"
    fidx = get_start_position(fname)
    if fidx == -1:
        insert_one(fname, ["Index", "File", "CrFinal", "Cr1", "T1", "Cr2", "T2..."])
    cur_idx = 0
    for root, dirs, files in os.walk(path_to_dataset):
        dirs.sort()
        for fl in sorted(files):
            if os.path.splitext(fl)[1] == ".lgbin":
                if cur_idx >= fidx:
                    optim = LayeredOptimizer(root + "/" + fl, {"cutoff_time": 300, "vertical_transitivity": True})
                    initial_layout_fn(optim.g)
                    output = optim.optimize_layout(local_opt=True, bucket_size=n_cvs)
                    reordered = [(output[0][i], output[1][i]) for i in range(len(output[0]))]
                    insert_one(fname, [cur_idx, root + "/" + fl] + reordered)
                cur_idx += 1


if __name__ == '__main__':
    # small_test()
    cv_sizes = [2000, 4000, 6000, 8000]
    cand_fns = [degree_candidate, biconnected_candidate, betweenness_candidate, avg_adge_length_candidate, crossings_candidate]
    nbhd_fns = [bfs_neighborhood, vertical_neighborhood]
    if len(sys.argv) >= 2:
        nbhd_idx = int(sys.argv[1]) // (len(cand_fns) * len(cv_sizes))
        cand_idx = (int(sys.argv[1]) % (len(cand_fns) * len(cv_sizes))) % len(cand_fns)
        cv_idx = (int(sys.argv[1]) % (len(cand_fns) * len(cv_sizes))) // len(cand_fns)
    else:
        nbhd_idx, cand_idx, cv_idx = 0, 0, 0
    run_experiment(nbhd_fns[nbhd_idx], cand_fns[cand_idx], cv_sizes[cv_idx], heuristics.global_sifting, "random graphs/rectangles")
    # opt = LayeredOptimizer("random graphs/n_by_n/n25/graph0.lgbin")
    # crossings_candidate(opt.g, init=True)
    # print([nd.energy for nd in opt.g])
