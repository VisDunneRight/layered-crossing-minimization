import itertools
import pickle
import random

from src.helpers import *


class LayeredNode:
	def __init__(self, node_name, node_layer, is_anchor=False, stacked=False):
		self.name = node_name
		self.layer = node_layer
		self.y = node_name
		self.is_anchor_node = is_anchor
		self.stacked = stacked

	def __str__(self):
		return f"ID={self.name}/L={self.layer}"

	def __repr__(self):
		return f"ID={self.name}/L={self.layer}"


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

	def __str__(self):
		return f"({str(self.n1)}, {str(self.n2)})"

	def __repr__(self):
		return f"({str(self.n1)}, {str(self.n2)})"


def find_closest(val, taken_vals: set):
	if round(val) not in taken_vals:
		return round(val)
	i = 1
	if val < round(val):
		while True:
			if round(val - i) not in taken_vals and round(val - i) > 0:
				return round(val - i)
			if round(val + i) not in taken_vals:
				return round(val + i)
			i += 1
	else:
		while True:
			if round(val + i) not in taken_vals:
				return round(val + i)
			if round(val - i) not in taken_vals and round(val - i) > 0:
				return round(val - i)
			i += 1


class LayeredGraph:
	def __init__(self):
		self.n_layers = 0
		self.n_nodes = 0
		self.nodes = []
		self.layers = {}
		self.edges = []
		self.stacked_edges = []
		self.stacked_nodes = []
		self.big_stack_nodes = []
		self.node_names = {}
		self.edge_names = {}
		self.adj_list = {}
		self.double_adj_list = {}
		self.time = 0

	def __getitem__(self, item):
		return self.node_names[item]

	def __iter__(self):
		return iter(self.nodes)

	def __contains__(self, item):
		return item in self.node_names

	def get_node(self, node_id):
		if 0 <= node_id < len(self.nodes):
			return self.nodes[node_id]
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
			name = self.n_nodes
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

	def add_edge(self, n1_name, n2_name, stacked=False, weight=1):
		if n1_name not in self.node_names or n2_name not in self.node_names:
			print(f"failed to add edge ({n1_name}, {n2_name}): node DNE")
			return
		if self.node_names[n1_name].layer > self.node_names[n2_name].layer:
			e = LayeredEdge(self.node_names[n2_name], self.node_names[n1_name], stacked=stacked, weight=weight)
			self.edges.append(e)
			self.edge_names[n2_name, n1_name] = e
		else:
			e = LayeredEdge(self.node_names[n1_name], self.node_names[n2_name], stacked=stacked, weight=weight)
			self.edges.append(e)
			self.edge_names[n1_name, n2_name] = e
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

	def stack_subgraph(self, cut_nodes: set, crossing_edges):  # DEPRECATED
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

	def stack_entire_graph(self, list_of_subgraphs):  # DEPRECATED
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

	def stacked_graph_from_subgraph_nodes(self, subgraph_assignments, only_subgraphs=False):
		"""
		:param: subgraph_assignments: List mapping with index node ID and elt subgraph assignment, integer in {0,...,#subgraphs-1}
		:return: new stacked graph object G', mapping of node ID -> stack node ID in G'
		"""
		new_g = CollapsedGraph(self)
		new_g.subgraphs = [[i for i, asgn in enumerate(subgraph_assignments) if asgn == j] for j in range(1 if only_subgraphs else 0, max(subgraph_assignments) + 1)]
		# node_to_group = {}
		new_g.node_to_stack_node = {}
		new_g.stack_node_to_nodelist = {}
		# crossing_edges = [[] for i in range(len(subgraphs))]
		new_g.crossing_edges = {}
		new_g.contact_nodes = [[] for i in range(len(new_g.subgraphs))]

		# for i, subg_list in enumerate(subgraphs):
		#     for node in subg_list:
		#         node_to_group[node] = i
		if only_subgraphs:
			for i, nd_val in enumerate(subgraph_assignments):
				if nd_val == 0:
					new_g.add_node(self.nodes[i].layer, name=i)

		new_g.n_nodes = self.n_nodes
		for subgraph in new_g.subgraphs:
			starting_node = new_g.n_nodes
			min_l = min((self.node_names[node].layer for node in subgraph))
			max_l = max((self.node_names[node].layer for node in subgraph))
			sn1 = new_g.add_node(min_l, stacked=True)
			new_g.stack_node_to_nodelist[sn1.name] = set()
			for level in range(min_l, max_l):
				sn2 = new_g.add_node(level + 1, stacked=True)
				new_g.stack_node_to_nodelist[sn2.name] = set()
				new_g.add_edge(sn1.name, sn2.name, stacked=True, weight=0)
				sn1 = sn2
			for node in subgraph:
				new_g.node_to_stack_node[node] = starting_node + (self.node_names[node].layer - min_l)
				new_g.stack_node_to_nodelist[starting_node + (self.node_names[node].layer - min_l)].add(node)
		for edge in self.edges:
			if subgraph_assignments[edge.n1.name] != subgraph_assignments[edge.n2.name]:
				if only_subgraphs:
					if subgraph_assignments[edge.n1.name] != 0 and subgraph_assignments[edge.n2.name] != 0:
						new_g.add_edge(new_g.node_to_stack_node[edge.n1.name], new_g.node_to_stack_node[edge.n2.name])
					elif subgraph_assignments[edge.n1.name] != 0:
						new_g.add_edge(new_g.node_to_stack_node[edge.n1.name], edge.n2.name)
					elif subgraph_assignments[edge.n2.name] != 0:
						new_g.add_edge(edge.n1.name, new_g.node_to_stack_node[edge.n2.name])
				else:
					if (new_g.node_to_stack_node[edge.n1.name], new_g.node_to_stack_node[edge.n2.name]) in new_g.edge_names:
						new_g.edge_names[new_g.node_to_stack_node[edge.n1.name], new_g.node_to_stack_node[edge.n2.name]].weight += 1    # FIXME handling of cases with multiple crossing edges between the same two stack nodes
					else:
						new_g.add_edge(new_g.node_to_stack_node[edge.n1.name], new_g.node_to_stack_node[edge.n2.name])
				# crossing_edges[node_to_group[edge.n1.name]].append((edge.n1.name, edge.n2.name))
				# crossing_edges[node_to_group[edge.n2.name]].append((edge.n1.name, edge.n2.name))
				if edge.n1.name not in new_g.crossing_edges:
					new_g.crossing_edges[edge.n1.name] = []
				new_g.crossing_edges[edge.n1.name].append(edge.n2.name)
				if only_subgraphs:
					new_g.contact_nodes[subgraph_assignments[edge.n1.name] - 1].append(edge.n1.name)
					new_g.contact_nodes[subgraph_assignments[edge.n2.name] - 1].append(edge.n2.name)
				else:
					new_g.contact_nodes[subgraph_assignments[edge.n1.name]].append(edge.n1.name)
					new_g.contact_nodes[subgraph_assignments[edge.n2.name]].append(edge.n2.name)
			elif not only_subgraphs or subgraph_assignments[edge.n1.name] != 0:
				new_g.edge_names[new_g.node_to_stack_node[edge.n1.name], new_g.node_to_stack_node[edge.n2.name]].weight += 1
			elif only_subgraphs:
				new_g.add_edge(edge.n1.name, edge.n2.name)

		new_g.n_nodes = len(new_g.nodes)
		return new_g

	def unstack_graph_nodes(self, stack_index):  # DEPRECATED
		self.nodes = [node for node in self.nodes if node not in self.big_stack_nodes[stack_index]]
		for s_node in self.big_stack_nodes[stack_index]:
			del self.node_names[s_node.name]
		for old_node in self.stacked_nodes[stack_index]:
			self.nodes.append(old_node)
			self.node_names[old_node.name] = old_node
		self.big_stack_nodes[stack_index].clear()

	def unstack_all_graph_edges(self):  # DEPRECATED
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
				self.double_adj_list[node.name] = []
			else:
				self.double_adj_list[node.name] = [[], []]
		for edge in self.edges:
			if forward_only:
				self.double_adj_list[edge.n1.name].append(edge.n2.name)
			else:
				self.double_adj_list[edge.n1.name][1].append(edge.n2.name)
				self.double_adj_list[edge.n2.name][0].append(edge.n1.name)
		return self.double_adj_list

	def create_normal_adj_list(self):
		for node in self.nodes:
			self.adj_list[node.name] = []
		for edge in self.edges:
			self.adj_list[edge.n1.name].append(edge.n2.name)
			self.adj_list[edge.n2.name].append(edge.n1.name)
		return self.adj_list

	# Clean up graph by removing empty layers and making sure the first layer has label 1.
	def relayer(self, remove_sl=True):
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
		removals = []
		for edge in self.edges:
			edge.update()
			if edge.same_layer_edge:
				removals.append(edge)
		if remove_sl:
			for to_remove in removals:
				del self.edge_names[to_remove.n1.name, to_remove.n2.name]
				self.edges.remove(to_remove)
		self.n_layers = len(self.layers)
		self.nodes.sort(key=lambda x: x.name)
		# for i in self.layers.keys():
		#     print(i, [n.name for n in self.layers[i]], [n.layer for n in self.layers[i]])

	def y_val_setup(self):
		for level in self.layers:
			self.layers[level].sort(key=lambda x: x.y)
			for i, node in enumerate(self.layers[level]):
				node.y = i + 1

	def barycentric_reordering(self, n_iter):  # DEPRECATED: Use heuristics.py
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

	def calculate_connectedness(self):
		max_connectedness = 0
		for l_num in sorted(list(self.layers.keys()))[:-1]:
			max_connectedness += len(self.layers[l_num]) * len(self.layers[l_num + 1])
		return len(self.edges) / max_connectedness

	def collapse_ap_cases(self, leaves_only=False):
		# 1) Recursive DFS to label all articulation points
		idx = 0
		disc = [self.n_nodes] * self.n_nodes
		low = [self.n_nodes] * self.n_nodes
		visited = [False] * self.n_nodes
		parent = [-1] * self.n_nodes
		aps = [False] * self.n_nodes
		time = 0
		if self.adj_list == {}:
			self.create_normal_adj_list()

		def tarjan_ap(index, visit, discovery, lowest, apoints, parents):
			nonlocal time
			visit[index] = True
			discovery[index] = time
			lowest[index] = time
			time += 1
			for idx_adj in self.adj_list[index]:
				if not visit[idx_adj]:
					parents[idx_adj] = index
					tarjan_ap(idx_adj, visit, discovery, lowest, apoints, parents)
					lowest[index] = min(lowest[index], lowest[idx_adj])
					if parents[index] != -1 and lowest[idx_adj] >= discovery[index]:
						apoints[index] = True
			if parents[index] == -1 and len([v for v in parents if v == index]) > 1:
				apoints[index] = True

		tarjan_ap(idx, visited, disc, low, aps, parent)

		# 2) For each AP, calc if largest <|G|/2 subgraph is valid AP case, add to list + calculate ratio, ow check each combo in order of decreasing size
		#   If an AP is already inside valid subgraph, ignore.
		ap_idxs = [i for i, v in enumerate(aps) if v]
		subgraphs_selected = []
		subgraph_types = []
		type_names = ["articulation point with leaves", "extreme layer subgraph", "<=2 node/layer subgraph", "left L-tree", "right L-tree"]
		aps_covered = [False] * len(ap_idxs)
		aps_subg_assign = [-1] * len(ap_idxs)
		for i, ap in enumerate(ap_idxs):
			if not aps_covered[i]:
				subgs_adj = [[v for v in self.adj_list[ap] if low[v] == i] for i in {low[x] for x in self.adj_list[ap]}]
				if subgs_adj:
					for j in range(len(subgs_adj)):
						ap_bfsq = [v for v in subgs_adj[j]]
						bfs_visit = [False] * self.n_nodes
						bfs_visit[ap] = True
						for v in ap_bfsq:
							bfs_visit[v] = True
						while ap_bfsq:
							next_layer = ap_bfsq.copy()
							ap_bfsq.clear()
							for nd in next_layer:
								for nd_adj in self.adj_list[nd]:
									if not bfs_visit[nd_adj]:
										bfs_visit[nd_adj] = True
										ap_bfsq.append(nd_adj)
										subgs_adj[j].append(nd_adj)
					subgs_adj = [temp for temp in subgs_adj if len(temp) <= self.n_nodes / 2]
					combos = []
					if len(subgs_adj) < 10:
						for r_val in range(1, len(subgs_adj) + 1):
							for temp in itertools.combinations(subgs_adj, r_val):
								if sum(len(lis) for lis in temp) <= self.n_nodes / 2:
									combos.append([cidx for temp2 in temp for cidx in temp2])
					else:
						for j in range(100):
							kv = random.randint(1, len(subgs_adj))
							samp = random.sample(subgs_adj, kv)
							if sum(len(lis) for lis in samp) <= self.n_nodes / 2:
								combos.append([cidx for temp2 in samp for cidx in temp2])
						# combos += [[cidx for temp2 in temp for cidx in temp2] for temp in itertools.combinations(subgs_adj, r)]
					combos.sort(key=lambda x: -len(x))
					for cmb in combos:
						valid, subg_type = self.check_if_collapsible_subgraph(cmb, ap, leaves_only=leaves_only)
						if valid:
							overwrite_idx = -1
							for k, ap2 in enumerate(ap_idxs):
								if ap2 in cmb:
									if aps_covered[k]:  # already covered ap situation
										overwrite_idx = aps_subg_assign[k]
									aps_covered[k] = True
							aps_covered[i] = True
							if overwrite_idx == -1:
								subgraphs_selected.append(cmb + [ap])
								subgraph_types.append(subg_type)
								for k, ap2 in enumerate(ap_idxs):
									if ap2 in cmb + [ap]:
										aps_subg_assign[k] = len(subgraphs_selected) - 1
								print(f"Collapsible subgraph of type [{type_names[subg_type]}] found: {cmb + [ap]}")
							else:
								subgraphs_selected[overwrite_idx] = cmb + [ap]
								subgraph_types[overwrite_idx] = subg_type
								for k, ap2 in enumerate(ap_idxs):
									if ap2 in cmb + [ap]:
										aps_subg_assign[k] = overwrite_idx
								print(f"Updated collapsible subgraph of type [{type_names[subg_type]}] found: {cmb + [ap]}")
							break

		# 4) Mark each collapsible as a subgraph and call collapse_subgraphs
		subgraphs_marked = [0] * self.n_nodes
		subgraphs_selected = [subg for i, subg in enumerate(subgraphs_selected) if i in set(aps_subg_assign)]  # remove extraneous subgraphs, occurring when an updated subgraph encompasses multiple already found subgraphs
		for mkid in range(len(subgraphs_selected)):
			for nd in subgraphs_selected[mkid]:
				subgraphs_marked[nd] = mkid + 1
		new_g = self.stacked_graph_from_subgraph_nodes(subgraphs_marked, only_subgraphs=True)
		new_g.subgraph_types = subgraph_types
		return new_g

	def check_if_collapsible_subgraph(self, subgraph, joint_idx, leaves_only=False) -> tuple[bool, int]:
		subgraph_layers = [self.nodes[nd].layer for nd in subgraph] + [self.nodes[joint_idx].layer]
		if len(set(subgraph_layers)) != len(subgraph_layers):
			subgraph_layers.remove(self.nodes[joint_idx].layer)
			if joint_idx in subgraph:
				subgraph.remove(joint_idx)
			if all((self.adj_list[x] == [joint_idx] for x in subgraph)):
				# articulation point with leaves
				return True, 0
			if not leaves_only:
				joint_layer = self.nodes[joint_idx].layer
				subg_set = set(subgraph)
				outsiders = [nd for nd in self.adj_list[joint_idx] if nd not in subg_set]
				if (all((self.nodes[x].layer < joint_layer for x in outsiders)) and all((x >= joint_layer for x in subgraph_layers))) or (all((self.nodes[x].layer > joint_layer for x in outsiders)) and all((x <= joint_layer for x in subgraph_layers))):
					# extreme layer subgraph
					return True, 1
				if len(outsiders) <= 1:
					if len([x for x in subgraph if self.nodes[x].layer == joint_layer]) <= 1:
						# <2 node per layer subgraph
						return True, 2
					n_left_parents = [len([x for x in self.adj_list[y] if self.nodes[x].layer == self.nodes[y].layer - 1]) for y in subgraph]
					if all((x <= 1 for x in n_left_parents)):
						# left l-tree
						return True, 3
					n_right_parents = [len([x for x in self.adj_list[y] if self.nodes[x].layer == self.nodes[y].layer + 1]) for y in subgraph]
					if all((x <= 1 for x in n_right_parents)):
						# right l-tree
						return True, 4
		return False, 0

	def vertex_exchange_graph(self):
		veg = {}  # vertex exchange graph as adj list
		nd_set = {}  #
		nd_list = []
		ed_sign = {}
		nd_idx = 0
		for ln in range(1, self.n_layers + 1):
			for nd1, nd2 in itertools.combinations(self.layers[ln], 2):
				veg[nd_idx] = []
				nd_list.append((nd1.name, nd2.name))
				nd_set[nd1.name, nd2.name] = nd_idx
				nd_idx += 1
		ebl = self.get_edges_by_layer()
		for ln in range(1, self.n_layers):
			for ed1, ed2 in itertools.combinations(ebl[ln], 2):
				if len({ed1.n1.name, ed1.n2.name, ed2.n1.name, ed2.n2.name}) == 4:
					if (ed1.n1.y > ed2.n1.y and ed1.n2.y > ed2.n2.y) or (ed1.n1.y < ed2.n1.y and ed1.n2.y < ed2.n2.y):
						sign = 0
					else:
						sign = 1
					if (ed1.n1.name, ed2.n1.name) in nd_set:
						if (ed1.n2.name, ed2.n2.name) in nd_set:
							veg[nd_set[ed1.n1.name, ed2.n1.name]].append((nd_set[ed1.n2.name, ed2.n2.name], sign))
							veg[nd_set[ed1.n2.name, ed2.n2.name]].append((nd_set[ed1.n1.name, ed2.n1.name], sign))
							ed_sign[nd_set[ed1.n1.name, ed2.n1.name], nd_set[ed1.n2.name, ed2.n2.name]] = sign
						else:
							veg[nd_set[ed1.n1.name, ed2.n1.name]].append((nd_set[ed2.n2.name, ed1.n2.name], sign))
							veg[nd_set[ed2.n2.name, ed1.n2.name]].append((nd_set[ed1.n1.name, ed2.n1.name], sign))
							ed_sign[nd_set[ed1.n1.name, ed2.n1.name], nd_set[ed2.n2.name, ed1.n2.name]] = sign
					else:
						if (ed1.n2.name, ed2.n2.name) in nd_set:
							veg[nd_set[ed2.n1.name, ed1.n1.name]].append((nd_set[ed1.n2.name, ed2.n2.name], sign))
							veg[nd_set[ed1.n2.name, ed2.n2.name]].append((nd_set[ed2.n1.name, ed1.n1.name], sign))
							ed_sign[nd_set[ed2.n1.name, ed1.n1.name], nd_set[ed1.n2.name, ed2.n2.name]] = sign
						else:
							veg[nd_set[ed2.n1.name, ed1.n1.name]].append((nd_set[ed2.n2.name, ed1.n2.name], sign))
							veg[nd_set[ed2.n2.name, ed1.n2.name]].append((nd_set[ed2.n1.name, ed1.n1.name], sign))
							ed_sign[nd_set[ed2.n1.name, ed1.n1.name], nd_set[ed2.n2.name, ed1.n2.name]] = sign
		return veg, nd_list, ed_sign

	def wiggle_node(self, x_vars, edge_b_l, node, pos_or_neg):
		best_seen_n_cr = 0
		start_x_vars = x_vars.copy()
		x_vars_changed = x_vars
		relevant_x_vars = {}
		for n_other in self.layers[self[node].layer]:
			if n_other.name != node:
				relevant_x_vars[node, n_other.name] = get_x_var(x_vars, node, n_other.name)

	def write_out(self, path):
		with open(path, 'wb') as bfd:
			pickle.dump(self, bfd)


class CollapsedGraph(LayeredGraph):
	def __init__(self, g: LayeredGraph):
		super().__init__()
		self.old_g = g
		self.subgraphs = None
		self.subgraph_types = None
		self.node_to_stack_node = None
		self.stack_node_to_nodelist = None
		self.crossing_edges = None
		self.contact_nodes = None

	def create_layered_graphs_from_subgraphs(self):
		subgraph_lgs = []
		for subg in self.subgraphs:
			subg_obj = LayeredGraph()
			seen_nds = set()
			for nd in subg:
				subg_obj.add_node(self.old_g.nodes[nd].layer, name=nd)
				seen_nds.add(nd)
				for adj in self.old_g.adj_list[nd]:
					if adj in seen_nds:
						subg_obj.add_edge(nd, adj)
			subgraph_lgs.append(subg_obj)
		return subgraph_lgs

	def create_layered_graphs_from_subgraphs_dangling_nodes(self):
		# TODO (later): combine dangling nodes from same node into single node with higher edge weight, mark as node to fix
		subgraph_lgs = []
		for subg in self.subgraphs:
			subg_obj = LayeredGraph()
			seen_nds = set(subg)
			for nd in subg:
				subg_obj.add_node(self.old_g.nodes[nd].layer, name=nd)
				seen_nds.remove(nd)
				for adj in self.old_g.adj_list[nd]:
					if adj not in seen_nds:
						if adj not in subg_obj:
							subg_obj.add_node(self.old_g.nodes[adj].layer, name=adj)
						subg_obj.add_edge(subg_obj[nd], subg_obj[adj])
			subgraph_lgs.append(subg_obj)
		return subgraph_lgs
