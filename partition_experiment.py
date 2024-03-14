import os.path
import sys
import csv
from src.optimization import LayeredOptimizer
from src.helpers import *
from src import heuristics, vis
from src.neighborhood import *
from os import mkdir, listdir
import json


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
            l1 = next(rdr)
            if l1[1] != "File":
                fset.add(l1[1])
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
    n_repeats = 50
    max_time = 300
    nbhds = ["bfs", "degree_ratio", "vertical_re", "random"]
    cands = ["betweenness", "crossings", "degree", "avg_edge_length", "random"]
    sizes = [10, 50, 100]
    ratios = [1.5]  # , 1, 2]
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
                            if "r1.5k" in gid:
                                all_data[nbhd][cand][gid] = []
                                for t1 in range(max_time):
                                    avg_cr = sum((t_avgs[j][t1] for j in range(n_repeats))) / n_repeats
                                    all_data[nbhd][cand][gid].append(avg_cr)
                                    # dat.append({"Graph": ln[1].split("/")[2], "Time": t1, "Crossings": avg_cr})

                            t_avgs.clear()
            dat = [[] for _ in range(5)]
            for nbhd in nbhds:
                for cand in cands:

                    for szidx, (gid, tvals) in enumerate(all_data[nbhd][cand].items()):  # this is for specific sized graphs
                        for t1, v in enumerate(tvals):
                            dat[szidx].append({"Graph": gid, "Method": f"{nbhd}+{cand}", "Time": t1, "Crossings": v})

                    # for t1 in range(300):  # this is for full averaging
                    #     avg_inc = 0
                    #     for _, tvals in all_data[nbhd][cand].items():
                    #         avg_inc += tvals[t1]
                    #     dat.append({"Method": f"{nbhd}+{cand}", "Time": t1, "Crossings": avg_inc / 6})

            if dat:
                for datblock in dat:
                    # ymax = max((datx["Crossings"] for datx in dat))
                    vis.draw_altair_simple_line_chart(datblock, "Time", "Crossings", "Method", "Time", "Crossings Improvement", f"mar8/{datblock[0]['Graph']}/{ratio}x{size}", ydom=[1, 2.4])

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


def create_averaged_json_file(folder):
    n_repeats = 50
    n_graph_sizes = 5
    max_time = 300
    nbhds = ["bfs", "degree_ratio", "vertical_re", "random"]
    cands = ["betweenness", "crossings", "degree", "avg_edge_length", "random"]
    sizes = [10, 50, 100]
    ratios = [1.5]  # , 1, 2]
    exps = ["r1.5k", "r1k", "r2k"]
    json_object = []
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
                            if "r1.5k" in gid:
                                all_data[nbhd][cand][gid] = []
                                for t1 in range(max_time):
                                    avg_cr = sum((t_avgs[j][t1] for j in range(n_repeats))) / n_repeats
                                    all_data[nbhd][cand][gid].append(avg_cr)
                                    # dat.append({"Graph": ln[1].split("/")[2], "Time": t1, "Crossings": avg_cr})

                            t_avgs.clear()

            for nbhd in nbhds:
                for cand in cands:

                    for gid, tvals in all_data[nbhd][cand].items():  # this is for specific sized graphs
                        nlay = int(gid[gid.index("k")+1:gid.index("k") + 3])
                        nnod = int(gid[gid.index("k")+4:])
                        json_object.append({"Graph": gid, "NodesPerLayer": nnod, "Layers": nlay, "Neighborhood": nbhd, "Candidate": cand, "NbhdSize": size, "CrossingImprovement": tvals})

                    # for t1 in range(300):  # this is for full averaging
                    #     avg_inc = 0
                    #     for _, tvals in all_data[nbhd][cand].items():
                    #         avg_inc += tvals[t1]
                    #     dat.append({"Method": f"{nbhd}+{cand}", "Time": t1, "Crossings": avg_inc / 6})

    with open(f"{folder}/averaged_results.json", "w") as fd:
        json.dump(json_object, fd, indent=4)


def create_raw_results_json(folder):
    n_repeats = 50
    nbhds = ["bfs", "degree_ratio", "vertical_re", "random"]
    cands = ["betweenness", "crossings", "degree", "avg_edge_length", "random"]
    gsizes = ["r1.5k18n12", "r1.5k24n16", "r1.5k30n20", "r1.5k36n24", "r1.5k42n28"]
    nsizes = [10, 50, 100]
    json_object = {"r1.5k18n12": {}, "r1.5k24n16": {}, "r1.5k30n20": {}, "r1.5k36n24": {}, "r1.5k42n28": {}}
    for nsize in nsizes:
        for kv in json_object:
            json_object[kv][nsize] = {}
        for nbhd in nbhds:
            for cand in cands:
                ncid = f"{nbhd}+{cand}"
                for kv in json_object:
                    json_object[kv][nsize][ncid] = []
                with open(f"{folder}/movement/{ncid}+{nsize}.csv", 'r') as fdm:
                    rdr = csv.reader(fdm)
                    move_lines = list(rdr)[1:]
                with open(f"{folder}/{ncid}+{nsize}.csv", 'r') as fd1:
                    rdr = csv.reader(fd1)
                    next(rdr)
                    gsize_idx = 0
                    for ln in rdr:
                        idx = int(ln[0])
                        crossins = [int(ln[v]) for v in range(5, len(ln), 2)]
                        if crossins[-1] == 0:
                            crossins[-1] = crossins[-2]
                        json_object[gsizes[gsize_idx]][nsize][ncid].append({
                            "gid": idx % n_repeats,
                            "sizeCalc": int(ln[2]),
                            "nOpts": int(move_lines[idx][3]),
                            "crInit": int(ln[5]),
                            "crFinal": int(ln[4]),
                            "candTimesMoved": int(move_lines[idx][4]),
                            "totalMoves": int(move_lines[idx][5]),
                            "timesteps": [float(ln[v]) for v in range(6, len(ln), 2)],
                            "crossings": crossins,
                            "nodeMovement": json.loads(move_lines[idx][6])
                        })
                        if idx % 50 == 49:
                            gsize_idx += 1

    with open(f"{folder}/nbhd+cand_results.json", 'w') as fdR:
        json.dump(json_object, fdR)
    # js1 = {"r1.5k18n12": json_object["r1.5k18n12"], "r1.5k24n16": json_object["r1.5k24n16"], "r1.5k30n20": json_object["r1.5k30n20"]}
    # js2 = {"r1.5k36n24": json_object["r1.5k36n24"], "r1.5k42n28": json_object["r1.5k42n28"]}


def print_exp_optcounts(results_folder_path):
    nb_idxs = {"bfs": 0, "degree_ratio": 1, "random": 2, "vertical_re": 3}
    cd_idxs = {"degree": 0, "betweenness": 1, "crossings": 2, "avg_edge_length": 3, "random": 4}
    all_dat = {}
    for fl in sorted(os.listdir(results_folder_path)):
        if os.path.splitext(fl)[1] == ".csv" and "000" not in fl:
            with open(results_folder_path + "/" + fl, 'r') as fd:
                rdr = csv.reader(fd)
                next(rdr)
                total_cts = 0
                nlines = 0
                grouped_by_size = {}
                grpcounts = {}
                for ln in rdr:
                    grpsize = ln[1].split('/')[2]
                    total_cts += (len(ln) - 7) / 2
                    nlines += 1
                    if grpsize not in grouped_by_size:
                        grouped_by_size[grpsize] = 0  # []
                        grpcounts[grpsize] = 0
                    grouped_by_size[grpsize] += (len(ln) - 6) / (2 * float(ln[3])) * 60  # .append
                    grpcounts[grpsize] += 1
                for k in grouped_by_size:  # mean not median
                    grouped_by_size[k] /= grpcounts[k]
                # print(fl, total_cts / nlines)
                nsize = "10" if "10." in fl else ("50" if "50." in fl else "100")
                nbhd = fl.split("+")[0]
                cand = fl.split("+")[1]
                for k, v in grouped_by_size.items():
                    if f"{k} : {nsize}" not in all_dat:
                        all_dat[f"{k} : {nsize}"] = [["0.0000000000"] * 4 for _ in range(5)]
                    all_dat[f"{k} : {nsize}"][cd_idxs[cand]][nb_idxs[nbhd]] = str(v)[:12]
                    # print(f"{k},{nsize},{v}")  # {v[len(v)//2]}")
                # nsize = 0 if "1000" in fl else (1 if "2000" in fl else 2)
                # gsize =
                # if "bfs" in fl:
                #     by_nbhd["bfs"][][nsize] += total_cts / nlines
    for k, v in all_dat.items():
        print(k)
        print("="*82)
        print(" " * 21, "bfs\t\t\t  degree_ratio\t  random\t\t  vertical_re")
        for idx, row in enumerate(v):
            the_cand = [cc for cc in cd_idxs if cd_idxs[cc] == idx][0]
            print(the_cand, " " * (20 - len(the_cand)), '    '.join(row))
        print("="*82, "\n\n")


def sandbox():
    # opt = LayeredOptimizer("random graphs/ratio_d3/r1k10n10/graph0.lgbin")
    opt = LayeredOptimizer("random graphs/triangle/r1.5k30n20/graph12.lgbin")
    # opt = LayeredOptimizer("random graphs/ratio_d3/r1.5k42n28/graph0.lgbin")
    opt.cutoff_time = 30
    opt.create_video = True
    opt.name = "r1.5k30n20.tri+deg"
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

    # print([idx for idx, v in enumerate(bfs_neighborhood(opt.g, 140, 0)) if v])
    # print(opt.g.n_nodes)

    opt.local_opt_increment(1500, neighborhood_fn=degree_ratio_neighborhood, candidate_fn=random_candidate, vertical_width=0)
    # opt.optimize_layout()

    # opt.m_val *= 2
    # opt.just_bendiness_reduction()
    # vis.draw_graph(opt.g, "solution_bend", groups=[0] * opt.g.n_nodes, label_nodes=False)
    # print(opt.g.num_edge_crossings())


def get_closest_cv(file_path, graph_size: str, target_avg: float, print_best=False):
    with open(file_path, 'r') as fd:
        rdr = csv.reader(fd)
        best_bnd = 0
        best_avg = sys.maxsize
        for ln in rdr:
            if ln[0] == graph_size:
                if abs(float(ln[2]) - target_avg) < abs(best_avg - target_avg):
                    best_bnd = int(ln[1])
                    best_avg = float(ln[2])
    if print_best:
        print(f"{graph_size} target {target_avg}: bound {best_bnd}, avg={best_avg}")
    return best_bnd


def add_cvar_to_csv():
    datapath = "random graphs/ratio_d3_results"
    bounds_path = "random graphs/ratio_d3/bounds_results"
    sizes = ["r1.5k18n12", "r1.5k24n16", "r1.5k30n20", "r1.5k36n24", "r1.5k42n28"]
    cands = ["degree", "random", "betweenness", "avg_edge_length", "crossings"]
    nbhd_size = [10, 50, 100]
    func_names = ["bfs", "vertical_re", "degree_ratio", "random"]
    for nbfn in func_names:
        for cdfn in cands:
            for nbsz in nbhd_size:
                fpath = f"{datapath}/{nbfn}+{cdfn}+{nbsz}.csv"
                cur_gsize = sizes[0]
                cur_idx = 0
                cv = get_closest_cv(f"{bounds_path}/{nbfn}_bounds.csv", cur_gsize, nbsz)
                with open(fpath, 'r') as fd:
                    rdr = csv.reader(fd)
                    rows = list(rdr)
                rows[0].insert(2, "SizeCalc")
                for row in rows[1:]:
                    if cur_gsize not in row[1]:
                        cur_idx += 1
                        cur_gsize = sizes[cur_idx]
                        cv = get_closest_cv(f"{bounds_path}/{nbfn}_bounds.csv", cur_gsize, nbsz)
                    row.insert(2, cv)
                with open(fpath, 'w') as fd:
                    wrt = csv.writer(fd)
                    wrt.writerows(rows)


def merge_csvs(movement=False):
    path_to_merge_from = "random graphs/ratio_d3/results" + "/movement" if movement else ""
    path_to_merge_into = "random graphs/temporary rd3 storage copy" + "/movement" if movement else ""
    # sizes = ["r1.5k18n12", "r1.5k24n16", "r1.5k30n20", "r1.5k36n24", "r1.5k42n28"]
    cands = ["degree", "random", "betweenness", "avg_edge_length", "crossings"]
    nbhd_size = [10, 50, 100]
    func_names = ["bfs", "vertical_re", "degree_ratio", "random"]
    for nbfn in func_names:
        for cdfn in cands:
            for nbsz in nbhd_size:
                fpath = f"{path_to_merge_from}/{nbfn}+{cdfn}+{nbsz}.csv"
                ipath = f"{path_to_merge_into}/{nbfn}+{cdfn}+{nbsz}.csv"
                with open(fpath, 'r') as fd1:
                    rdr = csv.reader(fd1)
                    lines_from = list(rdr)
                if os.path.isfile(ipath):
                    with open(ipath, 'r') as fd2:
                        rdr = csv.reader(fd2)
                        lines_into = list(rdr)
                    if movement:
                        last_idx = int(lines_from[-1][0])
                        lines_into = [ln for ln in lines_into if int(ln[0]) > last_idx]
                    else:
                        lines_into = lines_into[len(lines_from):]
                else:
                    lines_into = []
                new_lines = lines_from + lines_into
                with open(ipath, 'w') as fd:
                    wrt = csv.writer(fd)
                    wrt.writerows(new_lines)


def print_binsearch_results(path_to_results, restriction=0):
    sizes = ["r1.5k18n12", "r1.5k24n16", "r1.5k30n20", "r1.5k36n24", "r1.5k42n28"]
    nbhd_size = [10, 50, 100]
    func_names = ["bfs", "vertical_re", "degree_ratio", "random"]
    for nbfn in func_names:
        print(nbfn, "\n=============")
        for sz in sizes:
            for nbsz in nbhd_size:
                get_closest_cv(f"{path_to_results}/{nbfn}_bounds{'' if restriction == 0 else str(restriction)}.csv", sz, nbsz, print_best=True)
        print("==============\n")


def run_experiment(neighborhood_fn, candidate_fn, nbhd_size, initial_layout_fn, path_to_dataset, subdirs, num_graphs, restriction, movement=True):
    # if "results" not in listdir(path_to_dataset):
    #     mkdir(path_to_dataset + "/results")
    #     mkdir(path_to_dataset + "/results/movement")
    # if neighborhood_fn.__name__ not in listdir(path_to_dataset + "/results"):
    #     (path_to_dataset + "/results/" + neighborhood_fn.__name__.replace("_neighborhood", "") + "+" + candidate_fn.__name__.replace("_candidate", ""))
    nbhd_name = neighborhood_fn.__name__.replace("_neighborhood", "")
    fname = path_to_dataset + "/results/" + nbhd_name + "+" + candidate_fn.__name__.replace("_candidate", "") + "+" + str(nbhd_size) + (f"+{restriction}" if restriction != 0 else "") + ".csv"
    movefname = f"{path_to_dataset}/results/movement/{nbhd_name}+{candidate_fn.__name__.replace('_candidate', '')}+{str(nbhd_size)}{f'+{restriction}' if restriction != 0 else ''}.csv"
    files_run = get_start_position(fname)
    if len(files_run) == 0 and (not os.path.isfile(fname) or os.path.getsize(fname) == 0):
        insert_one(fname, ["Index", "File", "SizeCalc", "OptTime", "CrFinal", "Cr1", "T1", "Cr2", "T2..."])
    if movement and len(files_run) == 0 and (not os.path.isfile(movefname) or os.path.getsize(movefname) == 0):
        insert_one(movefname, ["Index", "File", "SizeCalc", "nOpts", "CandTimesMoved", "TotalMoves", "TimesMoved"])
    cur_idx = 0
    for subdir in subdirs:
        n_cvs = get_closest_cv(f"{path_to_dataset}/bounds_results/{nbhd_name}_bounds{str(restriction) if restriction != 0 else ''}.csv", subdir, nbhd_size)
        for fl_num in range(num_graphs):
            fl = f"graph{fl_num}.lgbin"
            file_path = f"{path_to_dataset}/{subdir}/{fl}"
            print(file_path)
            if file_path not in files_run:
                optim = LayeredOptimizer(file_path, cutoff_time=300)
                initial_layout_fn(optim.g)
                output = optim.local_opt_increment(n_cvs, neighborhood_fn=neighborhood_fn, candidate_fn=candidate_fn, movement_data=movement, vertical_width=restriction)
                reordered = [v for i in range(len(output[2])) for v in (output[2][i], output[3][i])]
                insert_one(fname, [cur_idx, file_path, n_cvs, output[0], output[1]] + reordered)
                if movement:
                    insert_one(movefname, [cur_idx, file_path, n_cvs, len(output[2]), output[5], output[6], output[4]])
            cur_idx += 1


if __name__ == '__main__':
    # sandbox()
    # draw_line_charts("random graphs/ratio_d3/results")
    # create_raw_results_json("random graphs/temporary rd3 storage copy")
    # print_exp_optcounts("./random graphs/temporary rd3 storage copy")
    # print_binsearch_results("random graphs/big_layer/bounds_results", restriction=0.5)
    # add_cvar_to_csv()
    # merge_csvs(movement=True)

    dataset_path = "random graphs/ratio_d3"
    subdirectories = ["r1.5k18n12", "r1.5k24n16", "r1.5k30n20", "r1.5k36n24", "r1.5k42n28"]
    num_graphs_in_subdir = 50
    nbhd_sizes = [10, 50, 100]
    cand_fns = [degree_candidate, random_candidate, betweenness_candidate, avg_edge_length_candidate, crossings_candidate]  # biconnected candidate
    nbhd_fns = [bfs_neighborhood, vertical_re_neighborhood, degree_ratio_neighborhood, random_neighborhood]
    datasets = ["ratio_d3", "big_layer", "triangle"]
    restrictions = [0, 0.75, 0.5]
    if 2 <= len(sys.argv) <= 3:
        nbhd_idx = int(sys.argv[1]) // (len(cand_fns) * len(nbhd_sizes))
        cand_idx = (int(sys.argv[1]) % (len(cand_fns) * len(nbhd_sizes))) % len(cand_fns)
        cv_idx = (int(sys.argv[1]) % (len(cand_fns) * len(nbhd_sizes))) // len(cand_fns)
        dataset_idx = int(sys.argv[2]) if len(sys.argv) > 2 else 0
        restrict_idx = 0
    elif len(sys.argv) > 3:  # Limited mobility experiment
        cand_idx = 1  # use random candidate
        nbhd_idx = int(sys.argv[1]) // len(nbhd_sizes)
        cv_idx = int(sys.argv[1]) % len(nbhd_sizes)
        dataset_idx = int(sys.argv[2])
        restrict_idx = int(sys.argv[3])
    else:
        nbhd_idx, cand_idx, cv_idx, dataset_idx, restrict_idx = 0, 1, 0, 1, 2

    if dataset_idx == 1:
        subdirectories.insert(0, "r1.5k12n8")
        subdirectories.pop()
    run_experiment(nbhd_fns[nbhd_idx], cand_fns[cand_idx], nbhd_sizes[cv_idx], heuristics.barycenter, f"random graphs/{datasets[dataset_idx]}", subdirectories, num_graphs_in_subdir, restrictions[restrict_idx], movement=True)
