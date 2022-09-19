class Node:
    def __init__(self, node_name, node_layer):
        self.name = node_name
        self.layer = node_layer
        self.v = 0


class Edge:
    def __init__(self, node1: Node, node2: Node):
        self.n1 = node1
        self.n2 = node2
        self.length = abs(self.n1.layer - self.n2.layer)
        self.same_layer_edge = self.n1.layer == self.n2.layer

    def get_bendiness(self):
        return abs(self.n1.v - self.n2.v)


class LayeredGraph:
    def __init__(self):
        self.n_layers = 0
        self.nodes = []
        self.layers = {}
        self.edges = []
        self.node_names = {}
        self.edge_names = {}

    # def is_name_available(self, name):
    #     for n in self.nodes:
    #         if n.name == name:
    #             return False
    #     return True

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

    def add_node(self, name, layer):
        x = Node(name, layer)
        if layer not in self.layers:
            self.layers[layer] = []
            self.n_layers += 1
        self.nodes.append(x)
        self.layers[layer].append(x)
        self.node_names[name] = x
        return x

    def add_nodes(self, names_and_layers):
        for nl in names_and_layers:
            self.add_node(nl[0], nl[1])

    def add_edge(self, n1, n2):
        if n1 not in self.node_names or n2 not in self.node_names:
            print(f"failed to add edge ({n1}, {n2}): node DNE")
            return
        if self.node_names[n1].layer > self.node_names[n2].layer:
            e = Edge(self.node_names[n2], self.node_names[n1])
            self.edges.append(e)
            self.edge_names[n2, n1] = e
        else:
            e = Edge(self.node_names[n1], self.node_names[n2])
            self.edges.append(e)
            self.edge_names[n1, n2] = e
        return e

    def add_graph_by_edges(self, edge_list):
        for edge in edge_list:
            if edge[0] not in self.node_names:
                self.add_node(edge[0], edge[1])
            elif edge[2] not in self.node_names:
                self.add_node(edge[2], edge[3])
            self.add_edge(edge[0], edge[2])

    def add_edges(self, edge_list):
        for edge in edge_list:
            self.add_edge(edge[0], edge[1])

    def add_anchors(self):
        to_remove = []
        to_add = []
        for edge in self.edges:
            if edge.length > 1:
                to_remove.append(edge)
                last_node = self.add_node(f"{edge.n1.name}-{edge.n2.name}a1", edge.n1.layer + 1)
                to_add.append((edge.n1.name, last_node.name))
                for i in range(edge.length-2):
                    next_node = self.add_node(f"{edge.n1.name}-{edge.n2.name}a{i+2}", i + edge.n1.layer + 2)
                    to_add.append((last_node.name, next_node.name))
                    last_node = next_node
                to_add.append((last_node.name, edge.n2.name))
        for v in to_remove:
            self.edges.remove(v)
        for e in to_add:
            self.add_edge(e[0], e[1])
