import src.graph
import re
import random


def create_bfs_layered_graph(s_g):
    visited = {n: False for n in s_g}
    bfs_q = set()
    # first = 11
    first = random.randint(1, len(s_g))
    bfs_q.add(first)
    visited[first] = True
    g = src.graph.LayeredGraph()
    g.add_node(str(first), 1)
    layer = 1
    while bfs_q:
        to_explore = bfs_q.copy()
        bfs_q.clear()
        layer += 1
        for n in to_explore:
            for adj in s_g[n]:
                if not visited[adj]:
                    visited[adj] = True
                    g.add_node(str(adj), layer)
                    bfs_q.add(adj)
                    g.add_edge(str(n), str(adj))
                elif adj in bfs_q:
                    g.add_edge(str(n), str(adj))
                elif g.get_node(str(adj)).layer == g.get_node(str(n)).layer and g.get_edge(str(adj), str(n)) is None:
                    g.add_edge(str(n), str(adj))
    return g


def create_layered_graph(rome_file):
    with open(f"Rome-Lib/{rome_file}") as f:
        simple_g = {}
        n_e = True
        for line in f.readlines():
            if line[0] == '#':
                n_e = False
                continue
            elif n_e:
                simple_g[int(line.split(' ')[0])] = []
            else:
                e = re.split('[ \n]', line)
                simple_g[int(e[2])].append(int(e[3]))
                simple_g[int(e[3])].append(int(e[2]))
        return create_bfs_layered_graph(simple_g)
