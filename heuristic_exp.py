import csv
from src.heuristics import *
from src import read_data
import os
import sys
import time


def insert_one(filename, entry):
    with open(filename, 'a', newline='') as f:
        wrt = csv.writer(f)
        wrt.writerow(entry)


def run_one(filename, gfile, fn, index):
    g = read_data.read(gfile)
    t = time.time()
    output = fn(g)
    insert_one(filename, [index, gfile, output, time.time() - t])


def get_start_pos(dataset):
    heuristic_names = ["bc", "med", "sift", "dot", "dwb", "ins", "sw", "sp", "bc_sift", "bc_sw", "bc_sp"]
    data_names = ["5_by_n", "10_by_10_density", "k_by_10", "rome"]
    data_lengths = [1000, 1700, 2400, 3983]
    for hname in heuristic_names[::-1]:
        if f"{data_names[dataset]}.csv" in os.listdir(f"data storage/heuristic/{hname}"):
            with open(f"data storage/heuristic/{hname}/{data_names[dataset]}.csv", 'r') as fd:
                lines = fd.readlines()
                last = lines[-1].split(',')
                if last[1] == "File":
                    continue
                elif last[0] == data_lengths[dataset] - 1:
                    return -1, -1
                else:
                    if dataset == 0:
                        if last[0][-2:] == "99":
                            return heuristic_names.index(hname), int(last[0]) // 100 * 10 + 20, 0
                        else:
                            return heuristic_names.index(hname), int(last[0]) // 100 * 10 + 10, int(last[0]) % 100 + 1
                    elif dataset == 1:
                        if last[0][-2:] == "99":
                            return heuristic_names.index(hname), int(last[0]) // 100 * 5 + 20, 0
                        else:
                            return heuristic_names.index(hname), int(last[0]) // 100 * 5 + 15, int(last[0]) % 100 + 1
                    elif dataset == 2:
                        if last[0][-2:] == "99":
                            return heuristic_names.index(hname), int(last[0]) // 100 + 3, 0
                        else:
                            return heuristic_names.index(hname), int(last[0]) // 100 + 2, int(last[0]) % 100 + 1
                    else:
                        return heuristic_names.index(hname), -1, int(last[0]) + 1
    if dataset == 0:
        return 0, 10, 0
    elif dataset == 1:
        return 0, 15, 0
    elif dataset == 2:
        return 0, 2, 0
    else:
        return 0, -1, 0


def run_heuristic_exp_checkpoint(dataset):
    heuristic_names = ["bc", "med", "sift", "dot", "dwb", "ins", "sw", "sp", "bc_sift", "bc_sw", "bc_sp"]
    heuristic_fns = [barycenter, median, global_sifting, weighted_median, degree_weighted_barycenter, greedy_insert, greedy_switching, split, improved_sifting, switching_with_preprocessing, barycenter_split]
    data_names = ["5_by_n", "10_by_10_density", "k_by_10", "rome"]
    if dataset == 3:
        rome_names = []
        with open("rome_crossings.csv", 'r') as fd:
            rdr = csv.reader(fd)
            next(rdr)
            for ln in rdr:
                rome_names.append(ln[0])
    if "heuristic" not in os.listdir("data storage"):
        os.mkdir("data storage/heuristic")
        os.mkdir("data storage/heuristic/bc")
        os.mkdir("data storage/heuristic/med")
        os.mkdir("data storage/heuristic/sift")
        os.mkdir("data storage/heuristic/dot")
        os.mkdir("data storage/heuristic/dwb")
        os.mkdir("data storage/heuristic/ins")
        os.mkdir("data storage/heuristic/sw")
        os.mkdir("data storage/heuristic/sp")
        os.mkdir("data storage/heuristic/bc_sift")
        os.mkdir("data storage/heuristic/bc_sw")
        os.mkdir("data storage/heuristic/bc_sp")

    hidx, fold, idx = get_start_pos(dataset)
    print(hidx, fold, idx)
    for i in range(hidx, len(heuristic_names)):
        hname = heuristic_names[i]
        print(hname)
        if (fold == 0 and idx == 0) or (dataset == 3 and idx == 0):
            insert_one(f"data storage/heuristic/{hname}/{data_names[dataset]}.csv", ["Index", "File", "Crossings", "Runtime"])
        if dataset == 0:
            if fold == 10 and idx == 0:
                insert_one(f"data storage/heuristic/{hname}/{data_names[dataset]}.csv", ["Index", "File", "Crossings", "Runtime"])
            for j in range(fold, 101, 10):
                for fileid in range(idx, 100):
                    print(hname, f"n={j} file", fileid)
                    run_one(f"data storage/heuristic/{hname}/{data_names[dataset]}.csv", f"random graphs/matuszewski/{data_names[dataset]}/n{j}/graph{fileid}.lgbin", heuristic_fns[i], 100 * ((j-1) // 10) + fileid)
                idx = 0
            fold = 10
        elif dataset == 1:
            if fold == 15 and idx == 0:
                insert_one(f"data storage/heuristic/{hname}/{data_names[dataset]}.csv", ["Index", "File", "Crossings", "Runtime"])
            for j in range(fold, 96, 5):
                for fileid in range(idx, 100):
                    print(hname, f"d={j} file", fileid)
                    run_one(f"data storage/heuristic/{hname}/{data_names[dataset]}.csv", f"random graphs/matuszewski/{data_names[dataset]}/d{j}/graph{fileid}.lgbin", heuristic_fns[i], 100 * ((j - 15) // 5) + fileid)
                idx = 0
            fold = 15
        elif dataset == 2:
            if fold == 2 and idx == 0:
                insert_one(f"data storage/heuristic/{hname}/{data_names[dataset]}.csv", ["Index", "File", "Crossings", "Runtime"])
            for j in range(fold, 26):
                for fileid in range(idx, 100):
                    print(hname, f"k={j} file", fileid)
                    run_one(f"data storage/heuristic/{hname}/{data_names[dataset]}.csv", f"random graphs/matuszewski/{data_names[dataset]}/k{j}/graph{fileid}.lgbin", heuristic_fns[i], 100 * (j - 2) + fileid)
                idx = 0
            fold = 2
        elif dataset == 3:
            if idx == 0:
                insert_one(f"data storage/heuristic/{hname}/{data_names[dataset]}.csv", ["Index", "File", "Crossings", "Runtime"])
            for fileid in range(3983):
                print(hname, "file", fileid)
                run_one(f"data storage/heuristic/{hname}/{data_names[dataset]}.csv", rome_names[fileid], heuristic_fns[i], fileid)
            idx = 0


if __name__ == '__main__':
    if len(sys.argv) >= 2:
        data_choice = int(sys.argv[1])
    else:
        data_choice = 0
    # run_heuristic_exp_checkpoint(data_choice)
    gr = read_data.read("random graphs/matuszewski/5_by_n/n10/graph17.lgbin")
    print(sorted([(v.n1.id, v.n2.id) for v in gr.edges], key=lambda x: x[0]))
    print(weighted_median(gr))
