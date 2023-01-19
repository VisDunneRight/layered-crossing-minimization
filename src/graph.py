import itertools
import math
from src.helpers import *


class LayeredNode:
    def __init__(self, node_name, node_layer, is_anchor=False, stacked=False):
        self.name = node_name
        self.layer = node_layer
        self.y = node_name
        self.is_anchor_node = is_anchor
        self.stacked = stacked


class LayeredEdge:
    def __init__(self, node1: LayeredNode, node2: LayeredNode, stacked=False, weight=1):
        self.n1 = node1
        self.n2 = node2
        self.length = abs(self.n1.layer - self.n2.layer)
        self.same_layer_edge = self.n1.layer == self.n2.layer
        self.stacked = stacked
        self.weight = weight

    def get_bendiness(self):
        return abs(self.n1.y - self.n2.y)

    def update(self):
        self.length = abs(self.n1.layer - self.n2.layer)
        self.same_layer_edge = self.n1.layer == self.n2.layer


def find_closest(val, taken_vals: set):
    if round(val, 1) not in taken_vals:
        return round(val, 1)
    i = 1
    if val < round(val, 1):
        while True:
            if round(val - i, 1) not in taken_vals and round(val - i, 1) > 0:
                return round(val - i, 1)
            if round(val + i) not in taken_vals:
                return round(val + i, 1)
            i += 1
    else:
        while True:
            if round(val + i, 1) not in taken_vals:
                return round(val + i, 1)
            if round(val - i, 1) not in taken_vals and round(val - i, 1) > 0:
                return round(val - i, 1)
            i += 1


class LayeredGraph:
    def __init__(self):
        self.n_layers = 0
        self.n_nodes = 0
        self.nodes = []
        self.layers = {}    # TODO convert to layernum->nodename instead of layernum->nodeobject
        self.edges = []
        self.stacked_edges = []
        self.stacked_nodes = []
        self.big_stack_nodes = []
        self.node_names = {}
        self.edge_names = {}
        self.adj_list = {}

    def __getitem__(self, item):
        return self.node_names[item]

    def __iter__(self):
        return iter(self.nodes)

    def __contains__(self, item):
        return item in self.node_names

    def get_node(self, node_name):
        if node_name in self.node_names:
            return self.node_names[node_name]
        return None

    def get_edge(self, node1_name, node2_name):
        if (node1_name, node2_name) in self.edge_names:
            return self.edge_names[node1_name, node2_name]
        return None

    def get_names(self):
        names = []
        for n in self.nodes:
            names.append(n.name)
        return names

    def get_names_by_layer(self):
        names = {}
        for n in self.nodes:
            if n.layer not in names:
                names[n.layer] = []
            names[n.layer].append(n.name)
        return names

    def get_edge_names_by_layer(self, only_diff_layer=False):
        edge_list = {}
        for edge in self.edges:
            if not only_diff_layer or edge.n1.layer != edge.n2.layer:
                if edge.n1.layer not in edge_list:
                    edge_list[edge.n1.layer] = []
                edge_list[edge.n1.layer].append((edge.n1.name, edge.n2.name))
        return edge_list

    def get_edges_by_layer(self, only_diff_layer=False):
        edge_list = {}
        for edge in self.edges:
            if not only_diff_layer or edge.n1.layer != edge.n2.layer:
                if edge.n1.layer not in edge_list:
                    edge_list[edge.n1.layer] = []
                edge_list[edge.n1.layer].append(edge)
        return edge_list

    def add_node(self, layer, name=None, is_anchor=False, stacked=False):
        if name is None:
            name = self.n_nodes + 1
        elif name in self.node_names:
            raise Exception(f"node {name} already exists in graph")
        x = LayeredNode(name, layer, is_anchor=is_anchor, stacked=stacked)
        if layer not in self.layers:
            self.layers[layer] = []
            self.n_layers += 1
        self.nodes.append(x)
        self.layers[layer].append(x)
        self.node_names[name] = x
        self.n_nodes += 1
        return x

    def add_nodes(self, names_and_layers):
        for nl in names_and_layers:
            self.add_node(nl[1], name=nl[0])

    def add_edge(self, n1, n2, stacked=False, weight=1):
        if n1 not in self.node_names or n2 not in self.node_names:
            print(f"failed to add edge ({n1}, {n2}): node DNE")
            return
        if self.node_names[n1].layer > self.node_names[n2].layer:
            e = LayeredEdge(self.node_names[n2], self.node_names[n1], stacked=stacked, weight=weight)
            self.edges.append(e)
            self.edge_names[n2, n1] = e
        else:
            e = LayeredEdge(self.node_names[n1], self.node_names[n2], stacked=stacked, weight=weight)
            self.edges.append(e)
            self.edge_names[n1, n2] = e
        return e

    def add_graph_by_edges(self, edge_list):
        for edge in edge_list:
            if edge[0] not in self.node_names:
                self.add_node(edge[1], name=edge[0])
            elif edge[2] not in self.node_names:
                self.add_node(edge[3], name=edge[2])
            self.add_edge(edge[0], edge[2])

    def add_edges(self, edge_list):
        for edge in edge_list:
            self.add_edge(edge[0], edge[1])

    def stack_subgraph(self, cut_nodes: set, crossing_edges):
        """
        Saves subgraph to new list in stacked_nodes/stacked_edges, deletes originals from nodes/edges/node_names
        Creates 1 new node, edge per layer with stacked=True and corresponding weight
        """
        level_seen = set()
        self.stacked_nodes.append([])
        self.stacked_edges.append([])
        self.big_stack_nodes.append([])
        for node in sorted(list(cut_nodes), key=lambda nod: self.node_names[nod].layer):
            if self.node_names[node].layer not in level_seen:
                level_seen.add(self.node_names[node].layer)
                x = self.add_node(self.node_names[node].layer, stacked=True)
                self.big_stack_nodes[-1].append(x)
            self.stacked_nodes[-1].append(self.node_names[node])
            self.nodes.remove(self.node_names[node])
            del self.node_names[node]
        cut_edges = []
        layer_counts = {}
        for edge in self.edges:
            if edge.n1.name in cut_nodes and edge.n2.name in cut_nodes:
                if edge.n1.layer not in layer_counts:
                    layer_counts[edge.n1.layer] = 0
                if not edge.same_layer_edge:
                    layer_counts[edge.n1.layer] += 1
                cut_edges.append(edge)
        for edge in cut_edges:
            self.stacked_edges[-1].append(edge)
            self.edges.remove(edge)
            del self.edge_names[edge.n1.name, edge.n2.name]
        for edge in crossing_edges:
            self.stacked_edges[-1].append(edge)
            self.edges.remove(edge)
            del self.edge_names[edge.n1.name, edge.n2.name]
            if edge.n1.name in cut_nodes:
                j = 0
                while self.big_stack_nodes[-1][j].layer != edge.n1.layer:
                    j += 1
                self.add_edge(self.big_stack_nodes[-1][j].name, edge.n2.name, stacked=True)
            else:
                j = 0
                while self.big_stack_nodes[-1][j].layer != edge.n2.layer:
                    j += 1
                self.add_edge(edge.n1.name, self.big_stack_nodes[-1][j].name, stacked=True)
        for i in range(len(self.big_stack_nodes[-1]) - 1):
            self.add_edge(self.big_stack_nodes[-1][i].name, self.big_stack_nodes[-1][i + 1].name, stacked=True, weight=layer_counts[self.big_stack_nodes[-1][i].layer])

    def stack_entire_graph(self, list_of_subgraphs):
        """
        Saves subgraphs to new lists in stacked_nodes/stacked_edges, deletes originals from nodes/edges/node_names
        Creates 1 new node, edge per layer for each subgraph with stacked=True and corresponding weight
        """
        min_max_layers = [(min((self.node_names[node].layer for node in subg_nodes)), max((self.node_names[node].layer for node in subg_nodes))) for subg_nodes in list_of_subgraphs]
        node_to_group = {}
        layer_counts = [{} for i in range(len(list_of_subgraphs))]
        for i, subg_list in enumerate(list_of_subgraphs):
            for node in subg_list:
                node_to_group[node] = i
        for i in range(len(list_of_subgraphs)):
            self.stacked_nodes.append([])
            self.stacked_edges.append([])
            self.big_stack_nodes.append([])
            for j in range(min_max_layers[i][0], min_max_layers[i][1]+1):
                x = self.add_node(j, stacked=True)
                self.big_stack_nodes[-1].append(x)
        for edge in list(self.edge_names).copy():
            if node_to_group[edge[0]] != node_to_group[edge[1]]:
                bsn1 = self.big_stack_nodes[node_to_group[edge[0]]][self.node_names[edge[0]].layer - min_max_layers[node_to_group[edge[0]]][0]]
                bsn2 = self.big_stack_nodes[node_to_group[edge[1]]][self.node_names[edge[1]].layer - min_max_layers[node_to_group[edge[1]]][0]]
                self.add_edge(bsn1.name, bsn2.name, stacked=True)
        for edge in self.edges:
            if not edge.stacked and node_to_group[edge.n1.name] == node_to_group[edge.n2.name]:
                if edge.n1.layer not in layer_counts[node_to_group[edge.n1.name]]:
                    layer_counts[node_to_group[edge.n1.name]][edge.n1.layer] = 0
                layer_counts[node_to_group[edge.n1.name]][edge.n1.layer] += 1
        print(layer_counts)
        for i, subg_nodestack in enumerate(self.big_stack_nodes):
            for j in range(len(subg_nodestack)-1):
                self.add_edge(subg_nodestack[j].name, subg_nodestack[j+1].name, stacked=True, weight=layer_counts[i][j+min_max_layers[i][0]])
        for node in self.nodes:
            if not node.stacked:
                self.stacked_nodes[node_to_group[node.name]].append(node)
                del self.node_names[node.name]
        self.nodes = [node for node in self.nodes if node.stacked]
        for edge in self.edges:
            if not edge.stacked:
                self.stacked_edges[node_to_group[edge.n1.name]].append(edge)
                del self.edge_names[edge.n1.name, edge.n2.name]
        self.edges = [edge for edge in self.edges if edge.stacked]

    def stacked_graph_from_subgraph_nodes(self, subgraph_assignments):
        """
        :param: subgraph_assignments: List mapping with index node ID minus 1 and elt subgraph assignment, integer in {0,...,#subgraphs-1}
        :return: new stacked graph object G', mapping of node ID -> stack node ID in G'
        """
        new_g = LayeredGraph()
        subgraphs = [[i + 1 for i, asgn in enumerate(subgraph_assignments[1:]) if asgn == j] for j in range(max(subgraph_assignments) + 1)]
        # node_to_group = {}
        node_to_stack_node = {}
        stack_node_to_nodelist = {}
        # crossing_edges = [[] for i in range(len(subgraphs))]
        crossing_edges = {}
        contact_nodes = [[] for i in range(len(subgraphs))]

        # for i, subg_list in enumerate(subgraphs):
        #     for node in subg_list:
        #         node_to_group[node] = i
        for subgraph in subgraphs:
            min_l = min((self.node_names[node].layer for node in subgraph))
            max_l = max((self.node_names[node].layer for node in subgraph))
            sn1 = new_g.add_node(min_l)
            stack_node_to_nodelist[sn1.name] = set()
            starting_node = new_g.n_nodes
            for level in range(min_l, max_l):
                sn2 = new_g.add_node(level + 1)
                stack_node_to_nodelist[sn2.name] = set()
                new_g.add_edge(sn1.name, sn2.name, weight=0)
                sn1 = sn2
            for node in subgraph:
                node_to_stack_node[node] = starting_node + (self.node_names[node].layer - min_l)
                stack_node_to_nodelist[starting_node + (self.node_names[node].layer - min_l)].add(node)
        for edge in self.edges:
            if subgraph_assignments[edge.n1.name] != subgraph_assignments[edge.n2.name]:
                if (node_to_stack_node[edge.n1.name], node_to_stack_node[edge.n2.name]) in new_g.edge_names:
                    new_g.edge_names[node_to_stack_node[edge.n1.name], node_to_stack_node[edge.n2.name]].weight += 1    # TODO fix handling of crossing edges between same stack nodes
                else:
                    new_g.add_edge(node_to_stack_node[edge.n1.name], node_to_stack_node[edge.n2.name])
                # crossing_edges[node_to_group[edge.n1.name]].append((edge.n1.name, edge.n2.name))
                # crossing_edges[node_to_group[edge.n2.name]].append((edge.n1.name, edge.n2.name))
                if edge.n1.name not in crossing_edges:
                    crossing_edges[edge.n1.name] = []
                crossing_edges[edge.n1.name].append(edge.n2.name)
                contact_nodes[subgraph_assignments[edge.n1.name]].append(edge.n1.name)
                contact_nodes[subgraph_assignments[edge.n2.name]].append(edge.n2.name)
            else:
                new_g.edge_names[node_to_stack_node[edge.n1.name], node_to_stack_node[edge.n2.name]].weight += 1

        return new_g, crossing_edges, contact_nodes, stack_node_to_nodelist, node_to_stack_node

    def unstack_graph_nodes(self, stack_index):
        self.nodes = [node for node in self.nodes if node not in self.big_stack_nodes[stack_index]]
        for s_node in self.big_stack_nodes[stack_index]:
            del self.node_names[s_node.name]
        for old_node in self.stacked_nodes[stack_index]:
            self.nodes.append(old_node)
            self.node_names[old_node.name] = old_node
        self.big_stack_nodes[stack_index].clear()

    def unstack_all_graph_edges(self):
        stacked_edge = [edge for edge in self.edges if edge.stacked]
        self.edges = [edge for edge in self.edges if not edge.stacked]
        for edge in stacked_edge:
            del self.edge_names[edge.n1.name, edge.n2.name]
        for stack_edge in self.stacked_edges:
            for edge in stack_edge:
                self.edges.append(edge)
                self.edge_names[edge.n1.name, edge.n2.name] = edge
        self.stacked_edges = []

    def adjacency_matrix(self):
        adj_list = self.create_normal_adj_list()
        new_matrix = [[0] * len(self.nodes) for i in range(len(self.nodes))]
        for v, l in adj_list.items():
            for u in l:
                new_matrix[v - 1][u - 1] = 1
                new_matrix[u - 1][v - 1] = 1
        return new_matrix

    def add_anchors(self):
        to_remove = []
        to_add = []
        for edge in self.edges:
            if edge.length > 1:
                to_remove.append(edge)
                last_node = self.add_node(edge.n1.layer + 1, is_anchor=True)
                to_add.append((edge.n1.name, last_node.name))
                for i in range(edge.length-2):
                    next_node = self.add_node(i + edge.n1.layer + 2, is_anchor=True)
                    to_add.append((last_node.name, next_node.name))
                    last_node = next_node
                to_add.append((last_node.name, edge.n2.name))
        for v in to_remove:
            self.edges.remove(v)
            del self.edge_names[v.n1.name, v.n2.name]
        for e in to_add:
            self.add_edge(e[0], e[1])

    def layer_counts(self):
        layer_counts = [len(lay) for lay in self.layers.values()]
        return f"Max height: {max(layer_counts)}", f"Layer counts: {layer_counts}"

    def create_double_adj_list(self, forward_only=False):
        for node in self.nodes:
            if forward_only:
                self.adj_list[node.name] = []
            else:
                self.adj_list[node.name] = [[], []]
        for edge in self.edges:
            if forward_only:
                self.adj_list[edge.n1.name].append(edge.n2.name)
            else:
                self.adj_list[edge.n1.name][1].append(edge.n2.name)
                self.adj_list[edge.n2.name][0].append(edge.n1.name)
        return self.adj_list

    def create_normal_adj_list(self):
        adjacency_list = {}
        for node in self.nodes:
            adjacency_list[node.name] = []
        for edge in self.edges:
            adjacency_list[edge.n1.name].append(edge.n2.name)
            adjacency_list[edge.n2.name].append(edge.n1.name)
        return adjacency_list

    def relayer(self):
        n_removals = min((n.layer for n in self.nodes)) - 1
        levels = sorted(list(self.layers.keys()))
        if n_removals > 0:
            for level in levels:
                self.layers[level - n_removals] = self.layers[level]
                for v in self.layers[level - n_removals]:
                    v.layer -= n_removals
                del self.layers[level]
        for node in self.nodes:
            if node.layer not in self.layers or node not in self.layers[node.layer]:
                for level in levels:
                    if level - n_removals in self.layers and node in self.layers[level - n_removals]:
                        self.layers[level - n_removals].remove(node)
                        if not self.layers[level - n_removals]:
                            del self.layers[level - n_removals]
                if node.layer not in self.layers:
                    self.layers[node.layer] = []
                self.layers[node.layer].append(node)
        for edge in self.edges:
            edge.update()
        self.n_layers = len(self.layers)
        # for i in self.layers.keys():
        #     print(i, [n.name for n in self.layers[i]], [n.layer for n in self.layers[i]])

    def y_val_setup(self):
        for level in self.layers:
            self.layers[level].sort(key=lambda x: x.y)
            for i, node in enumerate(self.layers[level]):
                node.y = i + 1

    def barycentric_reordering(self, n_iter):
        adjacency = self.create_normal_adj_list()
        min_y = min((n.y for n in self.nodes))
        for node_list in self.layers.values():
            min_l_y = min((n.y for n in node_list))
            if min_l_y > min_y:
                for n in node_list:
                    n.y -= min_l_y + min_y
        max_n_nodes = max((len(lay) for lay in self.layers.values()))
        for node_list in self.layers.values():
            for n in node_list:
                n.y += (max_n_nodes - len(node_list)) // 2
        # print(max(round(((max_n_nodes // 3) + 1) * (math.log10(n_iter)/2)), 1))
        for node in self.nodes:
            # node.y *= max(round(((max_n_nodes // 2) + 1) * (math.log(n_iter)/2)), 1)
            node.y *= 3
        # self.relayer()
        # self.y_val_setup()
        for level in self.layers:
            self.layers[level].sort(key=lambda x: -len(adjacency[x.name]))
        for i in range(n_iter):
            for j in range(1, len(self.layers) + 1):
                averages = [sum((self.node_names[m].y for m in adjacency[n.name]))/len(adjacency[n.name]) for n in self.layers[j]]
                taken = set()
                for k, node in enumerate(self.layers[j]):
                    insert_at = find_closest(averages[k], taken)
                    taken.add(insert_at)
                    node.y = insert_at
            for j in range(len(self.layers), 0, -1):
                averages = [sum((self.node_names[m].y for m in adjacency[n.name])) / len(adjacency[n.name]) for n in self.layers[j]]
                taken = set()
                for k, node in enumerate(self.layers[j]):
                    insert_at = find_closest(averages[k], taken)
                    taken.add(insert_at)
                    node.y = insert_at

    def num_edge_crossings(self):
        e_b_l = self.get_edges_by_layer()
        n_ec = 0
        for edge_list in e_b_l.values():
            for e1, e2 in itertools.combinations(edge_list, 2):
                if len({e1.n1, e1.n2, e2.n1, e2.n2}) == 4:
                    if e1.same_layer_edge and e2.same_layer_edge:
                        u1, w1, u2, w2 = e1.n1.y, e1.n2.y, e2.n1.y, e2.n2.y
                        if (u1 > u2 > w1 > w2) or (u1 > w2 > w1 > u2) or (w1 > u2 > u1 > w2) or (w1 > w2 > u1 > u2) or (u2 > u1 > w2 > w1) or (u2 > w1 > w2 > u1) or (w2 > u1 > u2 > w1) or (w2 > w1 > u2 > u1):
                            n_ec += 1
                    elif e1.same_layer_edge and ((e1.n1.y > e2.n1.y > e1.n1.y) or (e1.n2.y > e2.n1.y > e1.n1.y)):
                        n_ec += 1
                    elif e2.same_layer_edge and ((e2.n1.y > e1.n1.y > e2.n2.y) or (e2.n2.y > e1.n1.y > e2.n1.y)):
                        n_ec += 1
                    elif (e1.n1.y > e2.n1.y and e1.n2.y < e2.n2.y) or (e1.n1.y < e2.n1.y and e1.n2.y > e2.n2.y):
                        n_ec += 1
        return n_ec

    def num_edge_crossings_from_xvars_no_sl(self, x_vars):
        e_b_l = self.get_edges_by_layer()
        n_ec = 0
        for edge_list in e_b_l.values():
            for e1, e2 in itertools.combinations(edge_list, 2):
                if len({e1.n1, e1.n2, e2.n1, e2.n2}) == 4:
                    if (get_x_var(x_vars, e1.n1.name, e2.n1.name) and not get_x_var(x_vars, e1.n2.name, e2.n2.name)) or (not get_x_var(x_vars, e1.n1.name, e2.n1.name) and get_x_var(x_vars, e1.n2.name, e2.n2.name)):
                        n_ec += 1
        return n_ec

    def assign_y_vals_given_x_vars(self, x_vars):
        for nd in self.nodes:
            nd.y = 0
        for x_var, val in x_vars.items():
            self[x_var[val]].y += 1

    def wiggle_node(self, x_vars, edge_b_l, node, pos_or_neg):
        best_seen_n_cr = 0
        start_x_vars = x_vars.copy()
        x_vars_changed = x_vars
        relevant_x_vars = {}
        for n_other in self.layers[self[node].layer]:
            if n_other.name != node:
                relevant_x_vars[node, n_other.name] = get_x_var(x_vars, node, n_other.name)
