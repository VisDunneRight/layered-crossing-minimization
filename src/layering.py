import copy
import time
from src import graph, vis
import re
import random
import networkx as nx


def create_bfs_layered_graph(s_g):
    visited = {n: False for n in s_g}
    bfs_q = set()
    # first = 1
    first = random.randint(1, len(s_g))
    bfs_q.add(first)
    visited[first] = True
    g = graph.LayeredGraph()
    g.add_node(1, name=first)
    layer = 1
    while bfs_q:
        to_explore = bfs_q.copy()
        bfs_q.clear()
        layer += 1
        for n in to_explore:
            for adj in s_g[n]:
                if not visited[adj]:
                    visited[adj] = True
                    g.add_node(layer, name=adj)
                    bfs_q.add(adj)
                    g.add_edge(n, adj)
                elif adj in bfs_q:
                    g.add_edge(n, adj)
                elif g.get_node(adj).layer == g.get_node(n).layer and g.get_edge(adj, n) is None:
                    g.add_edge(n, adj)
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


def create_better_layered_graph(rome_file, w, c):
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
        # print("s_g", simple_g)
        to_remove = cycle_removal(simple_g)
        for edge in to_remove:
            simple_g[edge[0]].remove(edge[1])
        g, tvert = min_width(simple_g, w, c)
        for edge in to_remove:
            g.add_edge(edge[0], edge[1])
        # print("removed edges", to_remove)
        for edge in g.edges:
            edge.update()
        g.add_anchors()
        g.relayer()
        g.y_val_setup()
        # vis.draw_graph(g, "example3")
        return g, tvert


def create_edge_list_layered_graph(filepath, w, c):
    with open(filepath) as f:
        simple_g = {}
        for line in f.readlines():
            if line == '\n':
                continue
            e = re.split('[ ,\n]', line)
            e = [entry for entry in e if entry.isnumeric()]
            assert len(e) == 2, "file format incorrect. Each line should be integers of the form: u, v"
            if int(e[0]) not in simple_g:
                simple_g[int(e[0])] = []
            if int(e[1]) not in simple_g:
                simple_g[int(e[1])] = []
            simple_g[int(e[0])].append(int(e[1]))
    # print("s_g", simple_g)
    to_remove = cycle_removal(simple_g)
    for edge in to_remove:
        simple_g[edge[0]].remove(edge[1])
    g, tvert = min_width(simple_g, w, c)
    for edge in to_remove:
        g.add_edge(edge[0], edge[1])
    # print("removed edges", to_remove)
    for edge in g.edges:
        edge.update()
    g.add_anchors()
    g.relayer()
    g.y_val_setup()
    # vis.draw_graph(g, "example3")
    return g, tvert


def create_edge_list_layered_graph_given_layering(filepath, layer_assign):
    with open(filepath) as f:
        g = graph.LayeredGraph()
        for line in f.readlines():
            if line == '\n':
                continue
            e = re.split('[ ,\n]', line)
            e = [entry for entry in e if entry.isnumeric()]
            assert len(e) == 2, "file format incorrect. Each line should be integers of the form: u, v"
            if int(e[0]) not in g.node_names:
                g.add_node(layer_assign[int(e[0])], name=int(e[0]))
            if int(e[1]) not in g.node_names:
                g.add_node(layer_assign[int(e[1])], name=int(e[1]))
            g.add_edge(int(e[0]), int(e[1]))
    g.add_anchors()
    g.y_val_setup()
    return g


def create_layered_graph_from_directed_nx_graph(nxg: nx.Graph, w, c):
    simple_g = {}
    names = {nname: i+1 for i, nname in enumerate(nxg)}
    for node in nxg:
        simple_g[names[node]] = []
        for adj in nxg[node]:
            simple_g[names[node]].append(names[adj])
    to_remove = cycle_removal(simple_g)
    for edge in to_remove:
        simple_g[edge[0]].remove(edge[1])
    g = min_width(simple_g, w, c)[0]
    for edge in to_remove:
        g.add_edge(edge[0], edge[1])
    for edge in g.edges:
        edge.update()
    g.add_anchors()
    g.relayer()
    g.y_val_setup()
    return g


def cycle_removal(s_g):
    s_l = []
    s_r = []
    s_gp = copy.deepcopy(s_g)
    in_deg = {n: 0 for n in s_g}
    for adj_list in s_g.values():
        for j in adj_list:
            in_deg[j] += 1
    while s_gp:
        s_r_add = []
        s_l_add = []
        for v in s_gp:
            if not len(s_gp[v]):
                s_r_add.append(v)
                del in_deg[v]
                s_r.append(v)
        for v in s_r_add:
            del s_gp[v]
            for each_list in s_gp.values():
                if v in each_list:
                    each_list.remove(v)
        for u in in_deg:
            if not in_deg[u]:
                for adj in s_gp[u]:
                    in_deg[adj] -= 1
                s_l_add.append(u)
                del s_gp[u]
                s_l.append(u)
        for u in s_l_add:
            del in_deg[u]
        if s_gp:
            w = max(s_gp, key=lambda x: len(s_gp[x])-in_deg[x])
            for each_list in s_gp.values():
                if w in each_list:
                    each_list.remove(w)
            for adj in s_gp[w]:
                in_deg[adj] -= 1
            del s_gp[w]
            del in_deg[w]
            s_l.append(w)
    s_r.reverse()
    seen = set()
    removed_edges = []
    for v in s_l + s_r:
        for adj in s_g[v]:
            if adj in seen:
                removed_edges.append((v, adj))
        seen.add(v)
    return removed_edges


def min_width(s_g, w, c):
    g = graph.LayeredGraph()
    v_minus_u = sorted(list(s_g.keys()), key=lambda x: -len(s_g[x]))
    z_list = {n: False for n in s_g.keys()}
    u_list = []
    in_deg = {n: 0 for n in s_g.keys()}
    for adj_list in s_g.values():
        for j in adj_list:
            in_deg[j] += 1
    cur_layer, width_cur, width_up = 1, 0, 0
    while len(v_minus_u) > 0:
        chosen = False
        for node in v_minus_u:
            if all(z_list[v] is True for v in s_g[node]):
                chosen = True
                selected = node
                g.add_node(cur_layer, name=node)
                v_minus_u.remove(node)
                u_list.append(node)
                width_cur -= len(s_g[node]) - 1
                width_up += in_deg[node]
                break

        if not chosen or (width_cur >= w and len(s_g[selected]) < 1) or width_up >= c*w:
            cur_layer += 1
            for u in u_list:
                z_list[u] = True
            u_list.clear()
            width_cur = width_up
            width_up = 0

    for k, adj_list in s_g.items():
        for j in adj_list:
            g.add_edge(k, j)
    t = time.time()
    vertex_promotion(g)
    return g, time.time() - t


def promote_vertex(layering, g_dl, v):
    dummy = 0
    for behind in g_dl[v][1]:
        if layering[v] + 1 == layering[behind]:
            dummy += promote_vertex(layering, g_dl, behind)
    layering[v] += 1
    dummy += len(g_dl[v][0]) - len(g_dl[v][1])
    return dummy


def vertex_promotion(g: graph.LayeredGraph):
    double_adj = g.create_double_adj_list()
    layering = [0] + [g.node_names[i].layer for i in range(1, g.n_nodes + 1)]
    layer_backup = layering.copy()
    for i in range(100):
        promotions = 0
        for v in g.nodes:
            if len(double_adj[v.name][1]) > 0:
                if promote_vertex(layering, double_adj, v.name) < 0:
                    promotions += 1
                    layer_backup = layering.copy()
                else:
                    layering = layer_backup.copy()
        if promotions == 0:
            break
    for i, l in enumerate(layering):
        if i == 0:
            continue
        g.node_names[i].layer = l
