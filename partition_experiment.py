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
        fset = set()
        with open(filename, 'r') as fd:
            rdr = csv.reader(fd)
            for ln in rdr:
                fset.add(ln[1])
            # lines = fd.readlines()
            # last = lines[-1].split(',')
            # if last[0] == "Index":
            #     return 0
            # else:
            #     return int(last[0]) + 1
        return fset
    else:
        return set()


def get_csv_lines(file, exp):
    with open(file, 'r') as fd:
        rdr = csv.reader(fd)
        return [row for idx, row in enumerate(rdr) if exp in row[1]]


def draw_line_charts(folder):
    n_repeats = 10
    max_time = 300
    nbhds = ["bfs", "degree_ratio", "vertical_re", "random"]
    cands = ["betweenness", "crossings", "degree", "avg_edge_length", "random"]
    sizes = [1000, 2000, 3000]
    ratios = [1.5, 1, 2]
    exps = ["r1.5k", "r1k", "r2k"]
    for size in sizes:
        for ridx, ratio in enumerate(ratios):
            all_data = {}
            for nbhd in nbhds:
                all_data[nbhd] = {}
                for cand in cands:
                    print(size, ratio, nbhd, cand)
                    all_data[nbhd][cand] = {}
                    fl = f"{folder}/{nbhd}+{cand}+{size}.csv"
                    t_avgs = []
                    for i, ln in enumerate(get_csv_lines(fl, exps[ridx])):
                        st_cr = int(ln[4])
                        all_vals = []
                        for j in range(4, len(ln), 2):
                            if int(ln[j]) == 0 and j == len(ln) - 2:
                                # print(ln)
                                print(f"zero problem, line {ln[0]}")
                                ln[j] = ln[j - 2]
                            all_vals.append((st_cr / int(ln[j]), float(ln[j + 1])))
                        ptr = 0
                        ot_ptr = 0
                        t_avgs.append([])
                        for t1 in range(max_time):
                            while all_vals[ot_ptr][1] <= t1:
                                ot_ptr += 1
                            if ptr == ot_ptr:
                                t_avgs[-1].append(t_avgs[-1][-1])
                            else:
                                sum_x = sum((all_vals[i][0] for i in range(ptr, ot_ptr)))
                                t_avgs[-1].append(sum_x / (ot_ptr - ptr))
                            ptr = ot_ptr
                        if i % n_repeats == n_repeats - 1:
                            gid = ln[1].split("/")[2]
                            # if gid == "r1.5k36n24" or gid == "r1k30n30" or gid == "r2k40n20":  # med-large graphs
                            # if gid == "r1.5k18n12" or gid == "r1k15n15" or gid == "r2k20n10":  # small graphs
                            if gid == "r1.5k12n8":
                                all_data[nbhd][cand][gid] = []
                                for t1 in range(max_time):
                                    avg_cr = sum((t_avgs[j][t1] for j in range(10))) / 10
                                    all_data[nbhd][cand][gid].append(avg_cr)
                                    # dat.append({"Graph": ln[1].split("/")[2], "Time": t1, "Crossings": avg_cr})

                            t_avgs.clear()
            dat = []
            for nbhd in nbhds:
                for cand in cands:
                    for gid, tvals in all_data[nbhd][cand].items():  # this is for specific sized graphs
                        for t1, v in enumerate(tvals):
                            dat.append({"Graph": gid, "Method": f"{nbhd}+{cand}", "Time": t1, "Crossings": v})

                    # for t1 in range(300):  # this is for full averaging
                    #     avg_inc = 0
                    #     for _, tvals in all_data[nbhd][cand].items():
                    #         avg_inc += tvals[t1]
                    #     dat.append({"Method": f"{nbhd}+{cand}", "Time": t1, "Crossings": avg_inc / 6})

            if dat:
                ymax = max((datx["Crossings"] for datx in dat))
                vis.draw_altair_simple_line_chart(dat, "Time", "Crossings", "Method", "Time", "Crossings Improvement", f"feb20/r1.5k12n8/{ratio}x{size}", ydom=[1, ymax + 0.25])

    # for fl in os.listdir(folder):
    #     with open(f"{folder}/{fl}", 'r') as fd:
    #         dat = []
    #         rdr = csv.reader(fd)
    #         next(rdr)
    #         t_avgs = []
    #         for i, ln in enumerate(rdr):
    #             st_cr = int(ln[4])
    #             all_vals = []
    #             for j in range(4, len(ln), 2):
    #                 all_vals.append((1 - int(ln[j])/st_cr, float(ln[j+1]), i % n_repeats))
    #             ptr = 0
    #             ot_ptr = 0
    #             t_avgs.append([])
    #             for t1 in range(max_time):
    #                 while all_vals[ot_ptr][1] <= t1:
    #                     ot_ptr += 1
    #                 if ptr == ot_ptr:
    #                     t_avgs[-1].append(t_avgs[-1][-1])
    #                 else:
    #                     sum_x = sum((all_vals[i][0] for i in range(ptr, ot_ptr)))
    #                     t_avgs[-1].append(sum_x/(ot_ptr - ptr))
    #                 ptr = ot_ptr
    #
    #             if i % n_repeats == n_repeats - 1:
    #                 # dat.append({"Graph": ln[1].split("/")[2], "Time": 0, "Crossings": 0})
    #                 print(t_avgs)
    #                 for t1 in range(max_time):
    #                     avg_cr = sum((t_avgs[j][t1] for j in range(10))) / 10
    #                     dat.append({"Graph": ln[1].split("/")[2], "Time": t1, "Crossings": avg_cr})
    #                 t_avgs.clear()
    #
    #         vis.draw_altair_simple_line_chart(dat, "Time", "Crossings", "Graph", "Time", "Crossings Improvement", os.path.splitext(fl)[0])


def print_exp_optcounts(results_folder_path):
    by_nbhd = {"bfs": [[0] * 3 for _ in range(6)], "degree_ratio": [[0] * 3 for _ in range(6)], "random": [[0] * 3 for _ in range(6)], "vertical_re": [[0] * 3 for _ in range(6)]}
    for root, dirs, files in os.walk(results_folder_path):
        for fl in sorted(files):
            if os.path.splitext(fl)[1] == ".csv":
                with open(root + "/" + fl, 'r') as fd:
                    rdr = csv.reader(fd)
                    next(rdr)
                    total_cts = 0
                    nlines = 0
                    grouped_by_size = {}
                    grpcounts = {}
                    for ln in rdr:
                        grpsize = ln[1].split('/')[2]
                        total_cts += (len(ln) - 6) / 2
                        nlines += 1
                        if grpsize not in grouped_by_size:
                            grouped_by_size[grpsize] = 0
                            grpcounts[grpsize] = 0
                        grouped_by_size[grpsize] += (len(ln) - 6) / (2 * float(ln[2])) * 60
                        grpcounts[grpsize] += 1
                    for k in grouped_by_size:
                        grouped_by_size[k] /= grpcounts[k]
                    print(fl, total_cts / nlines)
                    nsize = "1000" if "1000" in fl else ("2000" if "2000" in fl else "3000")
                    for k, v in grouped_by_size.items():
                        print(f"{k},{nsize},{v}")
                    # nsize = 0 if "1000" in fl else (1 if "2000" in fl else 2)
                    # gsize =
                    # if "bfs" in fl:
                    #     by_nbhd["bfs"][][nsize] += total_cts / nlines
    # print(by_nbhd)
    # for nh in by_nbhd:
    #     by_nbhd[nh] = [v / 5 for i in range(len(by_nbhd[nh])) for v in by_nbhd[nh][i]]
    # print(by_nbhd)


def sandbox():
    # opt = LayeredOptimizer("random graphs/ratio_d3/r1k10n10/graph0.lgbin")
    opt = LayeredOptimizer("random graphs/ratio_d3/r1.5k30n20/graph12.lgbin")
    # opt = LayeredOptimizer("random graphs/ratio_d3/r1.5k42n28/graph0.lgbin")
    opt.cutoff_time = 30
    # opt.create_video = True
    opt.name = "r1.5k12n8.deg"
    # opt.vertical_transitivity = True
    # opt = LayeredOptimizer("Rome-Lib/graficon96nodi/grafo3510.96")
    # n_cv = opt.g.c_vars_count()

    heuristics.barycenter(opt.g)
    # opt.just_bendiness_reductiont()
    #     # vis.draw_graph(opt.g, "heurisic_bend", groups=[0] * opt.g.n_nodes, label_nodes=False)
    # vis.draw_graph(opt.g, "solution_heuristic")

    for lay in opt.g.layers.values():  # central gravity, for making gifs
        min_y = min((nd.y for nd in lay))
        for nd in lay:
            nd.y -= min_y

    print([idx for idx, v in enumerate(bfs_neighborhood(opt.g, 140, 0)) if v])
    print(opt.g.n_nodes)

    opt.local_opt_increment(0, neighborhood_fn=bfs_neighborhood, candidate_fn=random_candidate, vertical_width=0)
    # opt.optimize_layout()

    # opt.m_val *= 2
    # opt.just_bendiness_reduction()
    # vis.draw_graph(opt.g, "solution_bend", groups=[0] * opt.g.n_nodes, label_nodes=False)
    # print(opt.g.num_edge_crossings())


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


def run_experiment(neighborhood_fn, candidate_fn, nbhd_size, initial_layout_fn, path_to_dataset, subdirs, num_graphs):
    if "results" not in listdir(path_to_dataset):
        mkdir(path_to_dataset + "/results")
    # if neighborhood_fn.__name__ not in listdir(path_to_dataset + "/results"):
    #     (path_to_dataset + "/results/" + neighborhood_fn.__name__.replace("_neighborhood", "") + "+" + candidate_fn.__name__.replace("_candidate", ""))
    nbhd_name = neighborhood_fn.__name__.replace("_neighborhood", "")
    fname = path_to_dataset + "/results/" + nbhd_name + "+" + candidate_fn.__name__.replace("_candidate", "") + "+" + str(nbhd_size) + ".csv"
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
                output = optim.local_opt_increment(n_cvs, neighborhood_fn=neighborhood_fn, candidate_fn=candidate_fn)
                reordered = [v for i in range(len(output[2])) for v in (output[2][i], output[3][i])]
                insert_one(fname, [cur_idx, file_path, output[0], output[1]] + reordered)
            cur_idx += 1


if __name__ == '__main__':
    # sandbox()
    # draw_line_charts("random graphs/ratio_d3/results")
    # print_exp_optcounts("./random graphs/ratio_d3/results")

    dataset_path = "random graphs/ratio_d3"
    subdirectories = ["r1.5k18n12", "r1.5k24n16", "r1.5k30n20", "r1.5k36n24", "r1.5k42n28"]
    num_graphs_in_subdir = 50
    nbhd_sizes = [10, 50, 100]
    cand_fns = [degree_candidate, random_candidate, betweenness_candidate, avg_edge_length_candidate, crossings_candidate]  # biconnected candidate
    nbhd_fns = [bfs_neighborhood, vertical_re_neighborhood, degree_ratio_neighborhood, random_neighborhood]
    if len(sys.argv) >= 2:
        nbhd_idx = int(sys.argv[1]) // (len(cand_fns) * len(nbhd_sizes))
        cand_idx = (int(sys.argv[1]) % (len(cand_fns) * len(nbhd_sizes))) % len(cand_fns)
        cv_idx = (int(sys.argv[1]) % (len(cand_fns) * len(nbhd_sizes))) // len(cand_fns)
    else:
        nbhd_idx, cand_idx, cv_idx = 0, 1, 1
    run_experiment(nbhd_fns[nbhd_idx], cand_fns[cand_idx], nbhd_sizes[cv_idx], heuristics.barycenter, dataset_path, subdirectories, num_graphs_in_subdir)
