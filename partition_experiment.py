import os.path
import sys
import csv
from src.optimization import LayeredOptimizer
from src.helpers import *
from src import heuristics, vis
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


def draw_line_charts():
    for fl in os.listdir("random graphs/rectangles/results"):
        with open(f"random graphs/rectangles/results/{fl}", 'r') as fd:
            dat = []
            rdr = csv.reader(fd)
            next(rdr)
            all_vals = []
            for i, ln in enumerate(rdr):
                st_cr = int(ln[4])
                for j in range(6, len(ln), 2):
                    all_vals.append((st_cr - int(ln[j]), float(ln[j+1]), i % 10))
                if i % 10 == 9:
                    cur_vals = [0] * 10
                    dat.append({"Graph": ln[1].split("/")[2], "Time": 0, "Crossings": 0})
                    all_vals.sort(key=lambda x: x[1])
                    for cr, t, gid in all_vals:
                        cur_vals[gid] = cr
                        dat.append({"Graph": ln[1].split("/")[2], "Time": t, "Crossings": sum(cur_vals)/10})
                    all_vals.clear()
            vis.draw_altair_simple_line_chart(dat, "Time", "Crossings", "Graph", "Time", "Crossings Improvement", os.path.splitext(fl)[0])


def small_test():
    # opt = LayeredOptimizer("random graphs/n_by_n/n25/graph0.lgbin")
    opt = LayeredOptimizer("random graphs/ratio_d3/r1.5k30n20/graph0.lgbin")
    opt.cutoff_time = 120
    # opt.create_video = True
    opt.name = "r1.5k30n20g0"
    # opt.vertical_transitivity = True
    # opt = LayeredOptimizer("Rome-Lib/graficon96nodi/grafo3510.96")
    # n_cv = opt.g.c_vars_count()

    heuristics.barycenter(opt.g)
    # opt.just_bendiness_reductiont()
    #     # vis.draw_graph(opt.g, "heurisic_bend", groups=[0] * opt.g.n_nodes, label_nodes=False)
    # vis.draw_graph(opt.g, "solution_heuristic")

    opt.local_opt_increment(2000, neighborhood_fn=vertical_neighborhood, candidate_fn=avg_edge_length_candidate)

    # opt.m_val *= 2
    # opt.just_bendiness_reduction()
    # vis.draw_graph(opt.g, "solution_bend", groups=[0] * opt.g.n_nodes, label_nodes=False)
    print(opt.g.num_edge_crossings())


def run_experiment(neighborhood_fn, candidate_fn, n_cvs, initial_layout_fn, path_to_dataset):
    if "results" not in listdir(path_to_dataset):
        mkdir(path_to_dataset + "/results")
    # if neighborhood_fn.__name__ not in listdir(path_to_dataset + "/results"):
    #     (path_to_dataset + "/results/" + neighborhood_fn.__name__.replace("_neighborhood", "") + "+" + candidate_fn.__name__.replace("_candidate", ""))
    fname = path_to_dataset + "/results/" + neighborhood_fn.__name__.replace("_neighborhood", "") + "+" + candidate_fn.__name__.replace("_candidate", "") + "+" + str(n_cvs) + ".csv"
    fidx = get_start_position(fname)
    if fidx == -1:
        insert_one(fname, ["Index", "File", "OptTime", "CrFinal", "Cr1", "T1", "Cr2", "T2..."])
    cur_idx = 0
    for root, dirs, files in os.walk(path_to_dataset):
        dirs.sort()
        for fl in sorted(files):
            if os.path.splitext(fl)[1] == ".lgbin":
                if cur_idx >= fidx:
                    optim = LayeredOptimizer(root + "/" + fl, cutoff_time=300, vertical_transitivity=True)
                    initial_layout_fn(optim.g)
                    output = optim.local_opt_increment(n_cvs, neighborhood_fn=neighborhood_fn, candidate_fn=candidate_fn)
                    reordered = [v for i in range(len(output[2])) for v in (output[2][i], output[3][i])]
                    insert_one(fname, [cur_idx, root + "/" + fl, output[0], output[1]] + reordered)
                cur_idx += 1


if __name__ == '__main__':
    # small_test()
    # draw_line_charts()

    cv_sizes = [1000, 2000, 3000]
    cand_fns = [degree_candidate, random_candidate, betweenness_candidate, avg_edge_length_candidate, crossings_candidate]  # biconnected candidate
    nbhd_fns = [bfs_neighborhood, vertical_neighborhood, degree_ratio_neighborhood]
    if len(sys.argv) >= 2:
        nbhd_idx = int(sys.argv[1]) // (len(cand_fns) * len(cv_sizes))
        cand_idx = (int(sys.argv[1]) % (len(cand_fns) * len(cv_sizes))) % len(cand_fns)
        cv_idx = (int(sys.argv[1]) % (len(cand_fns) * len(cv_sizes))) // len(cand_fns)
    else:
        nbhd_idx, cand_idx, cv_idx = 1, 1, 1
    run_experiment(nbhd_fns[nbhd_idx], cand_fns[cand_idx], cv_sizes[cv_idx], heuristics.barycenter, "random graphs/ratio_d3")

    # opt = LayeredOptimizer("random graphs/n_by_n/n25/graph0.lgbin")
    # crossings_candidate(opt.g, init=True)
    # print([nd.energy for nd in opt.g])
