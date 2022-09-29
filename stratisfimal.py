import os
import sys
from src import constraints, vis, read_rome_data, motifs
from src.graph import *


def make_graph(choice=1):
    g = LayeredGraph()
    if choice == 1:
        g.add_nodes([("1", 1), ("2", 1), ("3", 2), ("4", 2), ("5", 3), ("6", 4), ("7", 2), ("8", 4), ("9", 1)])
        g.add_edges([("1", "2"), ("1", "3"), ("2", "4"), ("3", "5"), ("4", "6"), ("7", "8"), ("9", "3")])
        print([n.name for n in g.nodes], [(e.n1.name, e.n2.name) for e in g.edges])
        g.add_anchors()
        print([n.name for n in g.nodes], [(e.n1.name, e.n2.name) for e in g.edges])
        print(g.get_names_by_layer())
    elif choice == 2:
        g.add_nodes([("1", 1), ("2", 2), ("3", 2), ("4", 2)])
        g.add_edges([("1", "2"), ("1", "3"), ("1", "4")])
        print([n.name for n in g.nodes], [(e.n1.name, e.n2.name) for e in g.edges])
        g.add_anchors()
        print([n.name for n in g.nodes], [(e.n1.name, e.n2.name) for e in g.edges])
        print(g.get_names_by_layer())
    return g


def run_all_rome_lib(num_nodes, num_graphs, num_drawings, bendiness_reduction):
    i = 1
    outputs = []
    for file in os.listdir(f"Rome-Lib/graficon{num_nodes}nodi")[:num_graphs]:
        g = read_rome_data.create_layered_graph(f"graficon{num_nodes}nodi/{file}")
        print(f"\n\n{file} ({i}/{num_graphs}):")
        print("Number of butterflies: ", motifs.count_butterflies(g))
        outputs.append(f"{num_nodes}\t{bendiness_reduction}\t" + run_optimizer(g, bendiness_reduction))
        i += 1
        if num_drawings > 0:
            num_drawings -= 1
            vis.draw(g, f"Rome-Lib/{file}")
    for output in outputs:
        print(output)


def run_optimizer(g: LayeredGraph, bendiness_reduction):
    res = constraints.optimize_layout(g, 5, 1, bendiness_reduction, cutoff_time=60)
    return '\t'.join(str(e) for e in res[0] + [res[1], res[2], res[3]])


if __name__ == '__main__':
    # if len(sys.argv) < 2:
    #     graph_to_draw = make_graph()
    # else:
    #     graph_to_draw = make_graph(choice=int(sys.argv[1]))
    # run_optimizer(graph_to_draw)
    # vis.draw(graph_to_draw, "eex")
    n_nodes = int(input("Enter desired number of nodes in your graph (10 <= n <= 100): "))
    n_graphs = int(input("Enter desired number of graphs to optimize (1 <= n <= 100): "))
    n_drawings = int(input("Enter desired number of svg drawings to create (<= previous input): "))
    b_reduction = int(input("Perform bendiness reduction? Enter 0 or 1: "))
    run_all_rome_lib(n_nodes, n_graphs, n_drawings, b_reduction)
