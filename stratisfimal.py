import os
import sys
from src import constraints, vis, read_rome_data, motifs
from src.graph import *


def run_all_rome_lib(num_nodes, num_graphs, num_drawings, bendiness_reduction, seq_bend, timelimit, save, savefile=None):
    i = 1
    outputs = []
    for file in os.listdir(f"Rome-Lib/graficon{num_nodes}nodi")[:num_graphs]:
        g = read_rome_data.create_better_layered_graph(f"graficon{num_nodes}nodi/{file}", 4, 2)
        # g = read_rome_data.create_layered_graph(f"graficon{num_nodes}nodi/{file}")
        print(f"\n\n{file} ({i}/{num_graphs}):")
        print("Number of butterflies: ", motifs.count_butterflies(g))
        outputs.append(f"{num_nodes},{bendiness_reduction},{seq_bend},{timelimit}," + run_optimizer(g, bendiness_reduction, seq_bend, timelimit) + '\n')
        i += 1
        if num_drawings > 0:
            num_drawings -= 1
            vis.draw(g, f"Rome-Lib/{file}")
    if save:
        with open(savefile, 'a') as f:
            f.writelines(outputs)
    else:
        for output in outputs:
            print(output.replace('\n', ''))


def run_optimizer(g: LayeredGraph, bendiness_reduction, sequential, timelimit):
    if len(timelimit) > 0 and int(timelimit) > 0:
        res = constraints.optimize_layout(g, 5, 1, bendiness_reduction, sequential_bendiness=sequential, cutoff_time=int(timelimit))
    else:
        res = constraints.optimize_layout(g, 5, 1, bendiness_reduction, sequential_bendiness=sequential)
    return ','.join(str(e) for e in res[0] + res[1])


if __name__ == '__main__':
    if len(sys.argv) < 2:
        n_nodes = int(input("Enter desired number of nodes in your graph (10 <= n <= 100): "))
        n_g = input("Enter desired number of graphs to optimize (1 <= n <= 100, default 1): ")
        n_graphs = 1 if len(n_g) == 0 else int(n_g)
        n_d = input("Enter desired number of svg drawings to create (<= previous input, default 0): ")
        n_drawings = 0 if len(n_d) == 0 else int(n_d)
        b_r = input("Perform bendiness reduction? Enter 0 or 1, default 1: ")
        b_reduction = 1 if len(b_r) == 0 else int(b_r)
        timeout = input("Set time limit in seconds (enter 0 for no time limit, default no time limit): ")
        run_all_rome_lib(n_nodes, n_graphs, n_drawings, b_reduction, True, timeout, False)
    else:
        run_all_rome_lib(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), bool(sys.argv[4]), bool(sys.argv[5]), sys.argv[6], bool(sys.argv[7]), savefile=sys.argv[8])
