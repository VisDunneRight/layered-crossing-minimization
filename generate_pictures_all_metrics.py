import csv
import pandas as pd
import random
from aesthetics_exp import run_func
import os


def sample_random_rome_graphs(n_graphs):
    df = pd.read_csv("experiments/aesthetics/results/0_results.csv")
    gr1 = random.sample(df["Data"].tolist(), n_graphs // 2 + 5)
    df_sym = pd.read_csv("experiments/aesthetics/results/11_results.csv")
    gr2 = random.sample(df_sym["Data"].tolist(), n_graphs // 2 + 5)
    grcombo = (v + '\n' for v in set(gr1).union(set(gr2)))
    with open("rome_selected.txt", "w") as fd:
        fd.writelines(grcombo)


def get_data_that_completed(exp_id):
    df = pd.read_csv(f"experiments/aesthetics/results/{exp_id}_results.csv")
    set_of_rome_graphs = set()
    with open("rome_selected.txt") as fd:
        for ln in fd.readlines():
            set_of_rome_graphs.add(ln.removesuffix('\n'))
    df = df[df["Data"].isin(set_of_rome_graphs)]
    df = df[df["Status"] == 2]
    expected_rtime = df["Runtime"].sum()
    print(f"Expected runtime: {expected_rtime / 60} minutes")
    return df["Data"].tolist()


def run_all_data_that_completed(exp_id):
    if f"exp{exp_id}" not in os.listdir("Images/DrawingsForWebsite"):
        os.mkdir(f"Images/DrawingsForWebsite/exp{exp_id}")
    for romefl in get_data_that_completed(exp_id):
        grname = romefl.split('/')[1]
        run_func(exp_id, "Rome-Lib/" + romefl, drawing_filepath=f"DrawingsForWebsite/exp{exp_id}/{grname}")


if __name__ == '__main__':
    run_all_data_that_completed(0)
