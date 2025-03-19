import pandas as pd
import random
from aesthetics_exp import run_func
import os
import sys
import json
from src import read_data


def sample_random_rome_graphs(n_graphs):
    df = pd.read_csv("experiments/aesthetics/results/0_results.csv")
    gr1 = random.sample(df[df["Status"] == 2][["Data", "TotalNodes"]].values.tolist(), n_graphs // 2 + 2)
    df_sym = pd.read_csv("experiments/aesthetics/results/11_results.csv")
    gr2 = random.sample(df_sym[df_sym["Status"] == 2][["Data", "TotalNodes"]].values.tolist(), n_graphs // 2 + 2)
    grcombo = [','.join([v[0], str(v[1])]) + '\n' for v in gr1] + [','.join([v[0], str(v[1])]) + '\n' for v in gr2 if v not in gr1]
    with open("rome_selected.txt", "w") as fd:
        fd.writelines(grcombo)


def create_json_object_for_all_files():
    """ TODO: Add graphs that were rerun after cutoff and still completed """
    map_exp_names = ["Crossings", "EdgeLength", "CrossingAngle", "CrossingFairness", "LengthFairness", "NodeSymmetry", "EdgeSymmetry", "EdgeBundling", "MinMaxCrossings", "skip", "Crossings + EdgeLength", "Crossings + CrossingAngle", "Crossings + CrossingFairness", "Crossings + LengthFairness", "Crossings + NodeSymmetry", "Crossings + EdgeSymmetry", "Crossings + EdgeBundling", "Crossings + MinMaxCrossings", "skip", "EdgeLength (fixed x)", "CrossingAngle (fixed x)", "LengthFairness (fixed x)", "NodeSymmetry (fixed x)", "EdgeSymmetry (fixed x)", "EdgeLength + CrossingAngle (fixed x)", "EdgeLength + LengthFairness (fixed x)", "EdgeLength + NodeSymmetry (fixed x)", "Bend + EdgeSymmetry (fixed x)", "MaxPlanar", "Crossings + MaxPlanar"]
    full_data_obj = []
    the_rome_graphs = []
    with open("rome_selected.txt") as fd:
        for ln in fd.readlines():
            the_rome_graphs.append(ln.removesuffix('\n'))
    for expid in range(30):
        if map_exp_names[expid] != "skip":
            files_that_completed, tnodes, objvals = get_data_that_completed(expid)
            for fl in the_rome_graphs:
                if fl in files_that_completed:
                    idx = files_that_completed.index(fl)
                    grname = fl.split('/')[1].split('.')[0].replace("grafo", "graph")
                    full_data_obj.append({"src": f"images/exp{expid}/{grname}.png", "nodes": tnodes[idx], "exp": map_exp_names[expid], "obj": objvals[idx], "alt": f"{grname} / {tnodes[idx]} nodes / {map_exp_names[expid]} / Obj={objvals[idx]}", "stream": False})
                else:
                    gr = read_data.read(fl)
                    grname = fl.split('/')[1].split('.')[0].replace("grafo", "graph")
                    full_data_obj.append({"src": f"images/incomplete.png", "nodes": gr.n_nodes, "exp": map_exp_names[expid], "obj": -1, "alt": f"{grname} / {gr.n_nodes} nodes / {map_exp_names[expid]}", "stream": False})
    full_data_obj.sort(key=lambda x: (x["nodes"], map_exp_names.index(x["exp"])))
    with open("exp_json_data.json", "w") as fd:
        json.dump(full_data_obj, fd, indent=4)


def get_data_that_completed(exp_id, add_cutoff_graphs_up_to_2x=False):
    df = pd.read_csv(f"experiments/aesthetics/results/{exp_id}_results.csv")
    set_of_rome_graphs = []
    with open("rome_selected.txt") as fd:
        for ln in fd.readlines():
            set_of_rome_graphs.append(ln.removesuffix('\n').split(','))
    max_val = int(df.tail(1)["TotalNodes"].values[0])
    df = df[df["Data"].isin([v[0] for v in set_of_rome_graphs])]
    set_of_unrun_graphs = [[v[0], int(v[1])] for v in set_of_rome_graphs if v[0] not in df["Data"].values and int(v[1]) < 2 * max_val]
    set_of_unrun_graphs.sort(key=lambda xi: xi[1])
    df = df[df["Status"] == 2]
    expected_rtime = df["Runtime"].sum()
    print(f"Expected runtime, {len(df)} graphs: {expected_rtime / 60} minutes\nMax runtime if also running all {len(set_of_unrun_graphs)} graphs that were cutoff, 300sec timeout: {(expected_rtime + 300*len(set_of_unrun_graphs)) / 3600} hours")
    if add_cutoff_graphs_up_to_2x:
        return df["Data"].tolist() + [v[0] for v in set_of_unrun_graphs], df["TotalNodes"].tolist() + [v[1] for v in set_of_unrun_graphs], df["ObjVal"].tolist() + [-1] * len(set_of_unrun_graphs)
    else:
        return df["Data"].tolist(), df["TotalNodes"].tolist(), df["ObjVal"].tolist()


def run_all_data_that_completed(exp_id):
    if f"exp{exp_id}" not in os.listdir("Images/DrawingsForWebsite"):
        os.mkdir(f"Images/DrawingsForWebsite/exp{exp_id}")
    rfiles, _, _ = get_data_that_completed(exp_id, add_cutoff_graphs_up_to_2x=True)
    for i, romefl in enumerate(rfiles):
        grname = romefl.split('/')[1].split('.')[0].replace("grafo", "graph")
        if grname + ".png" not in os.listdir(f"Images/DrawingsForWebsite/exp{exp_id}"):
            run_func(exp_id, "Rome-Lib/" + romefl, drawing_filepath=f"DrawingsForWebsite/exp{exp_id}/{grname}", time_limit=200)


def run_all_data_that_completed_streamline(exp_id):
    if f"exp{exp_id}stream" not in os.listdir("Images/DrawingsForWebsite"):
        os.mkdir(f"Images/DrawingsForWebsite/exp{exp_id}stream")
    rfiles, _, _ = get_data_that_completed(exp_id, add_cutoff_graphs_up_to_2x=True)
    for i, romefl in enumerate(rfiles):
        grname = romefl.split('/')[1].split('.')[0].replace("grafo", "graph")
        if grname + ".png" not in os.listdir(f"Images/DrawingsForWebsite/exp{exp_id}stream"):
            run_func(exp_id, "Rome-Lib/" + romefl, drawing_filepath=f"DrawingsForWebsite/exp{exp_id}stream/{grname}", add_streamline=True, time_limit=200)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        for pid in list(range(9)) + list(range(10, 18)) + list(range(19, 30)):
            run_all_data_that_completed(pid)
    else:
        pid = int(sys.argv[1])
        if pid < 30:
            if pid != 9 and pid != 18:
                run_all_data_that_completed(pid)
        else:
            if pid - 30 != 9 and pid - 30 != 18:
                run_all_data_that_completed_streamline(pid - 30)
    # create_json_object_for_all_files()
    # sample_random_rome_graphs(150)
    # x = get_data_that_completed(24, add_cutoff_graphs_up_to_2x=True)
