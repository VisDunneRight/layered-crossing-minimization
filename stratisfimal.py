import sys
import constraints
import graph
import vis
from graph import *


def make_graph(choice=1):
    g = LayeredGraph()
    if choice == 1:
        g.add_nodes([("1",1),("2",1),("3",2),("4",2),("5",3),("6",4),("7",2),("8",4),("9",1)])
        g.add_edges([("1","2"),("1","3"),("2","4"),("3","5"),("4","6"),("7","8"),("9","3")])
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


def run_optimizer(g: graph.LayeredGraph):
    constraints.optimize_layout(g)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        graph_to_draw = make_graph()
    else:
        graph_to_draw = make_graph(choice=int(sys.argv[1]))
    run_optimizer(graph_to_draw)
    vis.draw(graph_to_draw, "eex")
