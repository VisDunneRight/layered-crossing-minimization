import itertools
import pickle
import random
from typing import Tuple
from src.helpers import *
import networkx as nx


class LayeredNode:
	def __init__(self, node_id, node_layer, is_anchor=False, stacked=False, name=None, fix=0):
		self.id = node_id
		self.layer = node_layer
		self.y = node_id
		self.is_anchor_node = is_anchor
		self.stacked = stacked
		self.name = node_id if name is None else name
		self.fix = fix
		self.energy = 0
		self.tabu = False

	def __str__(self):
		return f"ID={self.id}/L={self.layer}"

	def __repr__(self):
		return f"ID={self.id}/L={self.layer}"


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
		self.node_ids = {}
		self.node_names = {}
		self.edge_ids = {}
		self.ids_by_layer = {}
		self.edge_ids_by_layer = {}
		self.edges_by_layer = {}
		self.name_to_id = {}
		self.adj_list = {}
		self.double_adj_list = {}
		self.node_data = {}
		self.edge_data = {}
		self.nx_graph = None

	def __getitem__(self, item):
		if type(item) == int:
			return self.node_ids[item]
		else:
			return self.node_names[item]

	def __iter__(self):
		return iter(self.nodes)

	def __contains__(self, item):
		if type(item) == int:
			return item in self.node_ids
		else:
			return item in self.node_names

	# def get_node(self, node_id):
	# 	if 0 <= node_id < len(self.nodes):
	# 		return self.nodes[node_id]
	# 	return None

	def get_edge(self, node1_name, node2_name):
		if (node1_name, node2_name) in self.edge_ids:
			return self.edge_ids[node1_name, node2_name]
		return None

	def get_ids(self):
		names = []
		for n in self.nodes:
			names.append(n.id)
		return names

	def get_ids_by_layer(self):
		if self.ids_by_layer == {}:
			for n in self.nodes:
				if n.layer not in self.ids_by_layer:
					self.ids_by_layer[n.layer] = []
				self.ids_by_layer[n.layer].append(n.id)
		return self.ids_by_layer

	def get_edge_ids_by_layer(self, only_diff_layer=False):
		if self.edge_ids_by_layer == {}:
			for edge in self.edges:
				if not only_diff_layer or edge.n1.layer != edge.n2.layer:
					if edge.n1.layer not in self.edge_ids_by_layer:
						self.edge_ids_by_layer[edge.n1.layer] = []
					self.edge_ids_by_layer[edge.n1.layer].append((edge.n1.id, edge.n2.id))
		return self.edge_ids_by_layer

	def get_edges_by_layer(self, only_diff_layer=False):
		if self.edges_by_layer == {}:
			for edge in self.edges:
				if not only_diff_layer or edge.n1.layer != edge.n2.layer:
					if edge.n1.layer not in self.edges_by_layer:
						self.edges_by_layer[edge.n1.layer] = []
					self.edges_by_layer[edge.n1.layer].append(edge)
		return self.edges_by_layer

	def get_name_to_id(self):
		if self.name_to_id == {}:
			for i, nd in enumerate(self.nodes):
				self.name_to_id[nd.name] = i
		return self.name_to_id

	def get_adj_list(self):
		if self.adj_list == {}:
			self.create_normal_adj_list()
		return self.adj_list

	def get_double_adj_list(self):
		if self.double_adj_list == {}:
			self.create_double_adj_list()
		return self.double_adj_list

	def invalidate_data(self):
		if self.adj_list != {}:
			self.adj_list = {}
		if self.double_adj_list != {}:
			self.double_adj_list = {}
		if self.ids_by_layer != {}:
			self.ids_by_layer = {}
		if self.edge_ids_by_layer != {}:
			self.edge_ids_by_layer = {}
		if self.edges_by_layer != {}:
			self.edges_by_layer = {}
		if self.name_to_id != {}:
			self.name_to_id = {}

	def add_node(self, layer, idx=None, name=None, is_anchor=False, stacked=False, data=None):
		if idx is None:
			idx = self.n_nodes
		elif idx in self.node_ids:
			raise Exception(f"node {idx} already exists in graph")
		elif name is not None and name in self.node_names:
			raise Exception(f"node with name {name} already exists in graph")
		x = LayeredNode(idx, layer, is_anchor=is_anchor, stacked=stacked, name=name)
		if layer not in self.layers:
			self.layers[layer] = []
			self.n_layers += 1
		if data is not None:
			for k, v in data.items():
				if k not in self.node_data:
					self.node_data[k] = {}
				self.node_data[k][idx] = v
		self.nodes.append(x)
		self.layers[layer].append(x)
		self.node_ids[idx] = x
		if name is not None:
			self.node_names[name] = x
		# else:
		# 	self.node_names[idx] = x
		self.n_nodes += 1
		self.invalidate_data()
		return x

	def add_nodes(self, ids_and_layers):
		for nl in ids_and_layers:
			self.add_node(nl[1], idx=nl[0])

	def add_edge(self, n1_id, n2_id, stacked=False, weight=1, data=None):
		if n1_id not in self.node_ids or n2_id not in self.node_ids:
			raise Exception(f"failed to add edge ({n1_id}, {n2_id}): node DNE")
		if (n1_id, n2_id) in self.edge_ids:
			raise Exception(f"attempted to create duplicate edge ({n1_id}, {n2_id})")
		if self.node_ids[n1_id].layer > self.node_ids[n2_id].layer:
			e = LayeredEdge(self.node_ids[n2_id], self.node_ids[n1_id], stacked=stacked, weight=weight)
			self.edges.append(e)
			self.edge_ids[n2_id, n1_id] = e
		else:
			e = LayeredEdge(self.node_ids[n1_id], self.node_ids[n2_id], stacked=stacked, weight=weight)
			self.edges.append(e)
			self.edge_ids[n1_id, n2_id] = e
		if data is not None:
			for k, v in data.items():
				if k not in self.edge_data:
					self.edge_data[k] = {}
				self.edge_data[k][n1_id, n2_id] = v
		self.invalidate_data()
		return e

	def add_graph_by_edges(self, edge_list):
		for edge in edge_list:
			if edge[0] not in self.node_ids:
				self.add_node(edge[1], idx=edge[0])
			elif edge[2] not in self.node_ids:
				self.add_node(edge[3], idx=edge[2])
			self.add_edge(edge[0], edge[2])

	def add_edges(self, edge_list):
		for edge in edge_list:
			self.add_edge(edge[0], edge[1])

	# def stack_subgraph(self, cut_nodes: set, crossing_edges):  # DEPRECATED
	# 	"""
	# 	Saves subgraph to new list in stacked_nodes/stacked_edges, deletes originals from nodes/edges/node_ids
	# 	Creates 1 new node, edge per layer with stacked=True and corresponding weight
	# 	"""
	# 	level_seen = set()
	# 	self.stacked_nodes.append([])
	# 	self.stacked_edges.append([])
	# 	self.big_stack_nodes.append([])
	# 	for node in sorted(list(cut_nodes), key=lambda nod: self.node_ids[nod].layer):
	# 		if self.node_ids[node].layer not in level_seen:
	# 			level_seen.add(self.node_ids[node].layer)
	# 			x = self.add_node(self.node_ids[node].layer, stacked=True)
	# 			self.big_stack_nodes[-1].append(x)
	# 		self.stacked_nodes[-1].append(self.node_ids[node])
	# 		self.nodes.remove(self.node_ids[node])
	# 		del self.node_ids[node]
	# 	cut_edges = []
	# 	layer_counts = {}
	# 	for edge in self.edges:
	# 		if edge.n1.name in cut_nodes and edge.n2.name in cut_nodes:
	# 			if edge.n1.layer not in layer_counts:
	# 				layer_counts[edge.n1.layer] = 0
	# 			if not edge.same_layer_edge:
	# 				layer_counts[edge.n1.layer] += 1
	# 			cut_edges.append(edge)
	# 	for edge in cut_edges:
	# 		self.stacked_edges[-1].append(edge)
	# 		self.edges.remove(edge)
	# 		del self.edge_ids[edge.n1.name, edge.n2.name]
	# 	for edge in crossing_edges:
	# 		self.stacked_edges[-1].append(edge)
	# 		self.edges.remove(edge)
	# 		del self.edge_ids[edge.n1.name, edge.n2.name]
	# 		if edge.n1.name in cut_nodes:
	# 			j = 0
	# 			while self.big_stack_nodes[-1][j].layer != edge.n1.layer:
	# 				j += 1
	# 			self.add_edge(self.big_stack_nodes[-1][j].name, edge.n2.name, stacked=True)
	# 		else:
	# 			j = 0
	# 			while self.big_stack_nodes[-1][j].layer != edge.n2.layer:
	# 				j += 1
	# 			self.add_edge(edge.n1.name, self.big_stack_nodes[-1][j].name, stacked=True)
	# 	for i in range(len(self.big_stack_nodes[-1]) - 1):
	# 		self.add_edge(self.big_stack_nodes[-1][i].name, self.big_stack_nodes[-1][i + 1].name, stacked=True, weight=layer_counts[self.big_stack_nodes[-1][i].layer])

	# def stack_entire_graph(self, list_of_subgraphs):  # DEPRECATED
	# 	"""
	# 	Saves subgraphs to new lists in stacked_nodes/stacked_edges, deletes originals from nodes/edges/node_ids
	# 	Creates 1 new node, edge per layer for each subgraph with stacked=True and corresponding weight
	# 	"""
	# 	min_max_layers = [(min((self.node_ids[node].layer for node in subg_nodes)), max((self.node_ids[node].layer for node in subg_nodes))) for subg_nodes in list_of_subgraphs]
	# 	node_to_group = {}
	# 	layer_counts = [{} for _ in range(len(list_of_subgraphs))]
	# 	for i, subg_list in enumerate(list_of_subgraphs):
	# 		for node in subg_list:
	# 			node_to_group[node] = i
	# 	for i in range(len(list_of_subgraphs)):
	# 		self.stacked_nodes.append([])
	# 		self.stacked_edges.append([])
	# 		self.big_stack_nodes.append([])
	# 		for j in range(min_max_layers[i][0], min_max_layers[i][1]+1):
	# 			x = self.add_node(j, stacked=True)
	# 			self.big_stack_nodes[-1].append(x)
	# 	for edge in list(self.edge_ids).copy():
	# 		if node_to_group[edge[0]] != node_to_group[edge[1]]:
	# 			bsn1 = self.big_stack_nodes[node_to_group[edge[0]]][self.node_ids[edge[0]].layer - min_max_layers[node_to_group[edge[0]]][0]]
	# 			bsn2 = self.big_stack_nodes[node_to_group[edge[1]]][self.node_ids[edge[1]].layer - min_max_layers[node_to_group[edge[1]]][0]]
	# 			self.add_edge(bsn1.name, bsn2.name, stacked=True)
	# 	for edge in self.edges:
	# 		if not edge.stacked and node_to_group[edge.n1.name] == node_to_group[edge.n2.name]:
	# 			if edge.n1.layer not in layer_counts[node_to_group[edge.n1.name]]:
	# 				layer_counts[node_to_group[edge.n1.name]][edge.n1.layer] = 0
	# 			layer_counts[node_to_group[edge.n1.name]][edge.n1.layer] += 1
	# 	print(layer_counts)
	# 	for i, subg_nodestack in enumerate(self.big_stack_nodes):
	# 		for j in range(len(subg_nodestack)-1):
	# 			self.add_edge(subg_nodestack[j].name, subg_nodestack[j+1].name, stacked=True, weight=layer_counts[i][j+min_max_layers[i][0]])
	# 	for node in self.nodes:
	# 		if not node.stacked:
	# 			self.stacked_nodes[node_to_group[node.name]].append(node)
	# 			del self.node_ids[node.name]
	# 	self.nodes = [node for node in self.nodes if node.stacked]
	# 	for edge in self.edges:
	# 		if not edge.stacked:
	# 			self.stacked_edges[node_to_group[edge.n1.name]].append(edge)
	# 			del self.edge_ids[edge.n1.name, edge.n2.name]
	# 	self.edges = [edge for edge in self.edges if edge.stacked]

	def stacked_graph_from_subgraph_nodes(self, subgraph_assignments, only_subgraphs=False, keep_indices_same=True):
		"""
		:param: subgraph_assignments: List with node IDs as indices and subgraph assignment as values, assignment is integer in {0,...,#subgraphs-1}
		:param: only_subgraphs: True if only the subgraph assignments >=1 are to be collapsed (assignment=0 for all uncollapsed nodes), else False to collapse all nodes
		:return: new stacked graph object G', mapping of node ID -> stack node ID in G'
		"""
		new_g = CollapsedGraph(self)
		new_g.subgraphs = [[i for i, asgn in enumerate(subgraph_assignments) if asgn == j] for j in range(1 if only_subgraphs else 0, max(subgraph_assignments) + 1)]
		new_g.node_to_stack_node = {}
		new_g.stack_node_to_nodelist = {}
		new_g.crossing_edges = {}
		new_g.contact_nodes = [[] for _ in range(len(new_g.subgraphs))]
		old_node_to_new_node = {}

		# for i, subg_list in enumerate(subgraphs):
		#     for node in subg_list:
		#         node_to_group[node] = i
		if only_subgraphs:
			for i, nd_val in enumerate(subgraph_assignments):
				if nd_val == 0 and keep_indices_same:
					old_node_to_new_node[i] = new_g.add_node(self.node_ids[i].layer, idx=self.node_ids[i].name)
				elif nd_val == 0:
					new_g.add_node(self.node_ids[i].layer, idx=i)

		if not keep_indices_same:
			new_g.n_nodes = self.n_nodes
		for subgraph in new_g.subgraphs:
			starting_node = new_g.n_nodes
			min_l = min((self.node_ids[node].layer for node in subgraph))
			max_l = max((self.node_ids[node].layer for node in subgraph))
			l_ids = {}
			for node in subgraph:
				if self.node_ids[node].layer not in l_ids:
					l_ids[self.node_ids[node].layer] = node
			sn1 = new_g.add_node(min_l, stacked=True, idx=l_ids[min_l])
			new_g.stack_node_to_nodelist[sn1.id] = set()
			for level in range(min_l, max_l):
				sn2 = new_g.add_node(level + 1, stacked=True, idx=l_ids[level + 1])
				new_g.stack_node_to_nodelist[sn2.id] = set()
				new_g.add_edge(sn1.id, sn2.id, stacked=True, weight=0)
				sn1 = sn2
			for node in subgraph:
				# new_g.node_to_stack_node[node] = starting_node + (self.node_ids[node].layer - min_l)
				new_g.node_to_stack_node[node] = l_ids[self.node_ids[node].layer]
				# new_g.stack_node_to_nodelist[starting_node + (self.node_ids[node].layer - min_l)].add(node)
				new_g.stack_node_to_nodelist[l_ids[self.node_ids[node].layer]].add(node)
		for edge in self.edges:
			if subgraph_assignments[edge.n1.id] != subgraph_assignments[edge.n2.id]:
				if only_subgraphs:
					if subgraph_assignments[edge.n1.id] != 0 and subgraph_assignments[edge.n2.id] != 0:
						new_g.add_edge(new_g.node_to_stack_node[edge.n1.id], new_g.node_to_stack_node[edge.n2.id])
					elif subgraph_assignments[edge.n1.id] != 0:
						new_g.add_edge(new_g.node_to_stack_node[edge.n1.id], old_node_to_new_node[edge.n2.id].id)
					elif subgraph_assignments[edge.n2.id] != 0:
						new_g.add_edge(old_node_to_new_node[edge.n1.id].id, new_g.node_to_stack_node[edge.n2.id])
				else:
					if (new_g.node_to_stack_node[edge.n1.id], new_g.node_to_stack_node[edge.n2.id]) in new_g.edge_ids:
						new_g.edge_ids[new_g.node_to_stack_node[edge.n1.id], new_g.node_to_stack_node[edge.n2.id]].weight += edge.weight    # FIXME handling of cases with multiple crossing edges between the same two stack nodes
					else:
						new_g.add_edge(new_g.node_to_stack_node[edge.n1.id], new_g.node_to_stack_node[edge.n2.id])
				# crossing_edges[node_to_group[edge.n1.name]].append((edge.n1.name, edge.n2.name))
				# crossing_edges[node_to_group[edge.n2.name]].append((edge.n1.name, edge.n2.name))
				if edge.n1.id not in new_g.crossing_edges:
					new_g.crossing_edges[edge.n1.id] = []
				new_g.crossing_edges[edge.n1.id].append(edge.n2.id)
				if only_subgraphs:
					new_g.contact_nodes[subgraph_assignments[edge.n1.id] - 1].append(edge.n1.id)
					new_g.contact_nodes[subgraph_assignments[edge.n2.id] - 1].append(edge.n2.id)
				else:
					new_g.contact_nodes[subgraph_assignments[edge.n1.id]].append(edge.n1.id)
					new_g.contact_nodes[subgraph_assignments[edge.n2.id]].append(edge.n2.id)
			elif not only_subgraphs or subgraph_assignments[edge.n1.id] != 0:
				new_g.edge_ids[new_g.node_to_stack_node[edge.n1.id], new_g.node_to_stack_node[edge.n2.id]].weight += edge.weight
			elif only_subgraphs:
				new_g.add_edge(old_node_to_new_node[edge.n1.id].id, old_node_to_new_node[edge.n2.id].id)

		new_g.n_nodes = len(new_g.nodes)
		return new_g

	# def unstack_graph_nodes(self, stack_index):  # DEPRECATED
	# 	self.nodes = [node for node in self.nodes if node not in self.big_stack_nodes[stack_index]]
	# 	for s_node in self.big_stack_nodes[stack_index]:
	# 		del self.node_ids[s_node.name]
	# 	for old_node in self.stacked_nodes[stack_index]:
	# 		self.nodes.append(old_node)
	# 		self.node_ids[old_node.name] = old_node
	# 	self.big_stack_nodes[stack_index].clear()

	# def unstack_all_graph_edges(self):  # DEPRECATED
	# 	stacked_edge = [edge for edge in self.edges if edge.stacked]
	# 	self.edges = [edge for edge in self.edges if not edge.stacked]
	# 	for edge in stacked_edge:
	# 		del self.edge_ids[edge.n1.name, edge.n2.name]
	# 	for stack_edge in self.stacked_edges:
	# 		for edge in stack_edge:
	# 			self.edges.append(edge)
	# 			self.edge_ids[edge.n1.name, edge.n2.name] = edge
	# 	self.stacked_edges = []

	def adjacency_matrix(self):
		adj_list = self.get_adj_list()
		new_matrix = [[0] * len(self.nodes) for _ in range(len(self.nodes))]
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
				to_add.append((edge.n1.id, last_node.id))
				for i in range(edge.length-2):
					next_node = self.add_node(i + edge.n1.layer + 2, is_anchor=True)
					to_add.append((last_node.id, next_node.id))
					last_node = next_node
				to_add.append((last_node.id, edge.n2.id))
		for v in to_remove:
			self.edges.remove(v)
			del self.edge_ids[v.n1.id, v.n2.id]
		for e in to_add:
			self.add_edge(e[0], e[1])

	def layer_counts(self):
		layer_counts = [len(lay) for lay in self.layers.values()]
		return f"Max height: {max(layer_counts)}", f"Layer counts: {layer_counts}"

	def create_double_adj_list(self, forward_only=False):
		for node in self.nodes:
			if forward_only:
				self.double_adj_list[node.id] = []
			else:
				self.double_adj_list[node.id] = [[], []]
		for edge in self.edges:
			if forward_only:
				self.double_adj_list[edge.n1.id].append(edge.n2.id)
			else:
				self.double_adj_list[edge.n1.id][1].append(edge.n2.id)
				self.double_adj_list[edge.n2.id][0].append(edge.n1.id)
		return self.double_adj_list

	def create_normal_adj_list(self):
		for node in self.nodes:
			self.adj_list[node.id] = []
		for edge in self.edges:
			self.adj_list[edge.n1.id].append(edge.n2.id)
			self.adj_list[edge.n2.id].append(edge.n1.id)
		return self.adj_list

	# Clean up graph by removing empty layers and making sure the first layer has label 1.
	def relayer(self, remove_sl=True):
		n_removals = min((n.layer for n in self.nodes))
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
				del self.edge_ids[to_remove.n1.id, to_remove.n2.id]
				self.edges.remove(to_remove)
				print(f"Same-layer edge {to_remove} removed")
		# print(self.nodes)
		# for i in range(1, len(self.layers) + 1):  # change ids to 0 index
		# 	self.layers[i - 1] = self.layers[i]
		# 	del self.layers[i]
		self.n_layers = len(self.layers)
		self.nodes.sort(key=lambda x: x.id)
		self.invalidate_data()

	def y_val_setup(self):
		for level in self.layers:
			self.layers[level].sort(key=lambda x: x.y)
			for i, node in enumerate(self.layers[level]):
				node.y = i + 1

	def barycentric_reordering(self, n_iter):  # DEPRECATED: Use heuristics.py
		adjacency = self.get_adj_list()
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
			self.layers[level].sort(key=lambda x: -len(adjacency[x.id]))
		for i in range(n_iter):
			for j in range(len(self.layers)):
				averages = [sum((self.node_ids[m].y for m in adjacency[n.id])) / len(adjacency[n.id]) for n in self.layers[j]]
				taken = set()
				for k, node in enumerate(self.layers[j]):
					insert_at = find_closest(averages[k], taken)
					taken.add(insert_at)
					node.y = insert_at
			for j in range(len(self.layers)-1, -1, -1):
				averages = [sum((self.node_ids[m].y for m in adjacency[n.id])) / len(adjacency[n.id]) for n in self.layers[j]]
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
					if (get_x_var(x_vars, e1.n1.id, e2.n1.id) and not get_x_var(x_vars, e1.n2.id, e2.n2.id)) or (not get_x_var(x_vars, e1.n1.id, e2.n1.id) and get_x_var(x_vars, e1.n2.id, e2.n2.id)):
						n_ec += 1
		return n_ec

	def assign_y_vals_given_x_vars(self, x_vars):
		# x_vars: dictionary of x-var to assignment, e.g. {(3, 4): 1} if node 3 is below node 4
		for nd in self.nodes:
			nd.y = 0
		for x_var, val in x_vars.items():
			if x_var[0] in self.node_ids and x_var[1] in self.node_ids and val != 2:
				self[x_var[val]].y += 1

	def check_position_validity(self):
		checks_out = True
		for lay in self.layers.values():
			xl = sorted([(nd.id, nd.y) for nd in lay], key=lambda x: x[1])
			for i in range(len(xl)-1):
				if xl[i][1] == xl[i+1][1]:
					print(f"Invalid position: node {xl[i][0]} and {xl[i+1][0]} both have position {xl[i][1]}")
					checks_out = False
		if checks_out:
			print("Nodes have valid positions.")

	def calculate_connectedness(self):
		max_connectedness = 0
		for l_num in sorted(list(self.layers.keys()))[:-1]:
			max_connectedness += len(self.layers[l_num]) * len(self.layers[l_num + 1])
		return len(self.edges) / max_connectedness

	def is_connected(self):
		self.get_adj_list()
		visited = [False] * self.n_nodes
		visited[0] = True
		bfsq = [0]
		while bfsq:
			next_layer = bfsq.copy()
			bfsq.clear()
			for nd in next_layer:
				for nd_adj in self.adj_list[nd]:
					if not visited[nd_adj]:
						bfsq.append(nd_adj)
						visited[nd_adj] = True
		return all(visited)

	def collapse_ap_cases(self, leaves_only=False):
		# 1) Recursive DFS to label all articulation points
		idx = 0
		disc = [self.n_nodes] * self.n_nodes
		low = [self.n_nodes] * self.n_nodes
		visited = [False] * self.n_nodes
		parent = [-1] * self.n_nodes
		aps = [False] * self.n_nodes
		time = 0
		self.get_adj_list()

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

	def collapse_leaves(self):
		adj_list = self.get_adj_list()
		leaf_subgs = {}
		for nd in self.nodes:
			if len(adj_list[nd.id]) == 1:
				if adj_list[nd.id][0] not in leaf_subgs:
					leaf_subgs[adj_list[nd.id][0]] = [[], []]
				leaf_subgs[adj_list[nd.id][0]][0 if nd.layer < self.node_ids[adj_list[nd.id][0]].layer else 1].append(nd.id)

		subgraphs_marked = [0] * self.n_nodes
		subg_identifier = 1
		for root, leaf_lists in leaf_subgs.items():
			if len(leaf_lists[0]) >= 2:
				subgraphs_marked[root] = subg_identifier
				for lnd in leaf_lists[0]:
					subgraphs_marked[lnd] = subg_identifier
				print(f"Leaf subgraph found: {[root] + [leaf_lists[0]]}")
			if len(leaf_lists[1]) >= 2:
				subgraphs_marked[root] = subg_identifier
				for lnd in leaf_lists[1]:
					subgraphs_marked[lnd] = subg_identifier
				print(f"Leaf subgraph found: {[root] + [leaf_lists[1]]}")
			if len(leaf_lists[0]) >= 2 or len(leaf_lists[1]) >= 2:
				subg_identifier += 1

		new_g = self.stacked_graph_from_subgraph_nodes(subgraphs_marked, only_subgraphs=True)
		new_g.subgraph_types = [0] * (subg_identifier - 1)
		return new_g

	def check_if_collapsible_subgraph(self, subgraph, joint_idx, leaves_only=False) -> Tuple[bool, int]:
		subgraph_layers = [self.node_ids[nd].layer for nd in subgraph] + [self.node_ids[joint_idx].layer]
		if len(set(subgraph_layers)) != len(subgraph_layers):
			subgraph_layers.remove(self.node_ids[joint_idx].layer)
			if joint_idx in subgraph:
				subgraph.remove(joint_idx)
			if all((self.adj_list[x] == [joint_idx] for x in subgraph)):
				# articulation point with leaves
				return True, 0
			if not leaves_only:
				joint_layer = self.node_ids[joint_idx].layer
				subg_set = set(subgraph)
				outsiders = [nd for nd in self.adj_list[joint_idx] if nd not in subg_set]
				if (all((self.node_ids[x].layer < joint_layer for x in outsiders)) and all((x >= joint_layer for x in subgraph_layers))) or (all((self.node_ids[x].layer > joint_layer for x in outsiders)) and all((x <= joint_layer for x in subgraph_layers))):
					# extreme layer subgraph
					return True, 1
				if len(outsiders) <= 1:
					if len([x for x in subgraph if self.node_ids[x].layer == joint_layer]) <= 1:
						# <2 node per layer subgraph
						return True, 2
					n_left_parents = [len([x for x in self.adj_list[y] if self.node_ids[x].layer == self.node_ids[y].layer - 1]) for y in subgraph]
					if all((x <= 1 for x in n_left_parents)):
						# left l-tree
						return True, 3
					n_right_parents = [len([x for x in self.adj_list[y] if self.node_ids[x].layer == self.node_ids[y].layer + 1]) for y in subgraph]
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
		for ln in range(self.n_layers):
			for nd1, nd2 in itertools.combinations(self.layers[ln], 2):
				veg[nd_idx] = []
				nd_list.append((nd1.id, nd2.id))
				nd_set[nd1.id, nd2.id] = nd_idx
				nd_idx += 1
		ebl = self.get_edges_by_layer()
		for ln in range(self.n_layers - 1):
			for ed1, ed2 in itertools.combinations(ebl[ln], 2):
				if len({ed1.n1.id, ed1.n2.id, ed2.n1.id, ed2.n2.id}) == 4:
					if (ed1.n1.y > ed2.n1.y and ed1.n2.y > ed2.n2.y) or (ed1.n1.y < ed2.n1.y and ed1.n2.y < ed2.n2.y):
						sign = 0
					else:
						sign = 1
					if (ed1.n1.id, ed2.n1.id) in nd_set:
						if (ed1.n2.id, ed2.n2.id) in nd_set:
							veg[nd_set[ed1.n1.id, ed2.n1.id]].append((nd_set[ed1.n2.id, ed2.n2.id], sign))
							veg[nd_set[ed1.n2.id, ed2.n2.id]].append((nd_set[ed1.n1.id, ed2.n1.id], sign))
							ed_sign[nd_set[ed1.n1.id, ed2.n1.id], nd_set[ed1.n2.id, ed2.n2.id]] = sign
						else:
							veg[nd_set[ed1.n1.id, ed2.n1.id]].append((nd_set[ed2.n2.id, ed1.n2.id], sign))
							veg[nd_set[ed2.n2.id, ed1.n2.id]].append((nd_set[ed1.n1.id, ed2.n1.id], sign))
							ed_sign[nd_set[ed1.n1.id, ed2.n1.id], nd_set[ed2.n2.id, ed1.n2.id]] = sign
					else:
						if (ed1.n2.id, ed2.n2.id) in nd_set:
							veg[nd_set[ed2.n1.id, ed1.n1.id]].append((nd_set[ed1.n2.id, ed2.n2.id], sign))
							veg[nd_set[ed1.n2.id, ed2.n2.id]].append((nd_set[ed2.n1.id, ed1.n1.id], sign))
							ed_sign[nd_set[ed2.n1.id, ed1.n1.id], nd_set[ed1.n2.id, ed2.n2.id]] = sign
						else:
							veg[nd_set[ed2.n1.id, ed1.n1.id]].append((nd_set[ed2.n2.id, ed1.n2.id], sign))
							veg[nd_set[ed2.n2.id, ed1.n2.id]].append((nd_set[ed2.n1.id, ed1.n1.id], sign))
							ed_sign[nd_set[ed2.n1.id, ed1.n1.id], nd_set[ed2.n2.id, ed1.n2.id]] = sign
		return veg, nd_list, ed_sign

	# def wiggle_node(self, x_vars, edge_b_l, node, pos_or_neg):
	# 	best_seen_n_cr = 0
	# 	start_x_vars = x_vars.copy()
	# 	x_vars_changed = x_vars
	# 	relevant_x_vars = {}
	# 	for n_other in self.layers[self[node].layer]:
	# 		if n_other.id != node:
	# 			relevant_x_vars[node, n_other.id] = get_x_var(x_vars, node, n_other.id)

	def get_networkx_graph(self, update=False):
		try:
			if self.nx_graph is None or update:
				self.nx_graph = nx.DiGraph()
				self.nx_graph.add_nodes_from(self.node_ids.keys())
				self.nx_graph.add_edges_from(self.edge_ids.keys())
		except AttributeError:
			self.nx_graph = nx.DiGraph()
			self.nx_graph.add_nodes_from(self.node_ids.keys())
			self.nx_graph.add_edges_from(self.edge_ids.keys())
		return self.nx_graph

	def c_vars_count(self):
		c_count = 0
		e_b_l = self.get_edges_by_layer()
		for elist in e_b_l.values():
			c_count += len(elist) ** 2
		return c_count

	def write_out(self, path):
		with open(path, 'wb') as bfd:
			pickle.dump(self, bfd)


class CollapsedGraph(LayeredGraph):
	def __init__(self, g: LayeredGraph, subgraphs=None):
		super().__init__()
		self.old_g = g
		self.subgraphs = None
		self.subgraph_types = None
		self.node_to_stack_node = None
		self.stack_node_to_nodelist = None
		self.crossing_edges = None
		self.contact_nodes = None
		# self.skeleton = None
		if subgraphs is not None:
			self.subgraphs = subgraphs

	def get_collapsed_node(self, layer, subg_id):
		for nd in self.layers[layer]:
			if self.subgraphs[self.stack_node_to_nodelist[nd.id][0]] == subg_id:
				return nd
		return -1

	def create_layered_graphs_from_subgraphs(self):
		subgraph_lgs = []
		old_adj = self.old_g.get_adj_list()
		subgs = [[idx for idx, v in enumerate(self.subgraphs) if v == i] for i in range(max(self.subgraphs) + 1)]
		for sid, subg in enumerate(subgs):
			subg_obj = LayeredGraph()
			seen_nds = set()
			for nd in subg:
				subg_obj.add_node(self.old_g.nodes[nd].layer, idx=nd)
				seen_nds.add(nd)
				for adj in old_adj[nd]:
					if adj in seen_nds:
						subg_obj.add_edge(nd, adj)
			subgraph_lgs.append(subg_obj)
			subg_obj.subg_id = sid
		return subgraph_lgs

	def create_layered_graphs_from_subgraphs_dangling_nodes(self):
		subgraph_lgs = []
		old_adj = self.old_g.get_adj_list()
		subgs = [[idx for idx, v in enumerate(self.subgraphs) if v == i] for i in range(max(self.subgraphs)+1)]
		already_connected = {}
		for sid, subg in enumerate(subgs):
			subg_obj = LayeredGraph()
			seen_nds = set(subg)
			for nd in subg:
				subg_obj.add_node(self.old_g[nd].layer, idx=nd)
				seen_nds.remove(nd)
				for adj in old_adj[nd]:
					if adj not in seen_nds:
						# the line below corrects for when a node has >1 adjacent contact nodes, but in diff layers
						nd_adj_val = nd + 0.5 if self.old_g[adj].layer > self.old_g[nd].layer else nd
						if adj not in subg_obj.node_ids and nd_adj_val not in already_connected:  # adj in diff subg
							subg_obj.add_node(self.old_g.nodes[adj].layer, idx=adj)
							already_connected[nd_adj_val] = adj
							subg_obj.add_edge(nd, adj)
						elif adj not in subg_obj.node_ids and nd_adj_val in already_connected:
							if (nd, already_connected[nd_adj_val]) in subg_obj.edge_ids:
								subg_obj.edge_ids[nd, already_connected[nd_adj_val]].weight += 1
							else:
								subg_obj.edge_ids[already_connected[nd_adj_val], nd].weight += 1
						else:
							subg_obj.add_edge(nd, adj)
			subgraph_lgs.append(subg_obj)
			subg_obj.subg_id = sid
		return subgraph_lgs

	def create_collapsed_graph_skeleton(self):
		self.node_to_stack_node = {}
		self.stack_node_to_nodelist = {}
		self.crossing_edges = []
		self.contact_nodes = []
		subg_layer_combos = {}

		for nd in self.old_g.nodes:
			if (nd.layer, self.subgraphs[nd.id]) not in subg_layer_combos:
				xnd = self.add_node(nd.layer)
				subg_layer_combos[nd.layer, self.subgraphs[nd.id]] = xnd.id
				self.stack_node_to_nodelist[xnd.id] = []
			else:
				xnd = self.node_ids[subg_layer_combos[nd.layer, self.subgraphs[nd.id]]]
			self.node_to_stack_node[nd.id] = xnd.id
			self.stack_node_to_nodelist[xnd.id].append(nd.id)

		for ed in self.old_g.edges:  # find crossing edges and nodes
			if self.subgraphs[ed.n1.id] != self.subgraphs[ed.n2.id]:
				self.crossing_edges.append((ed.n1.id, ed.n2.id))
				if ed.n1.id not in self.contact_nodes:
					self.contact_nodes.append(ed.n1.id)
				if ed.n2.id not in self.contact_nodes:
					self.contact_nodes.append(ed.n2.id)
			if (self.node_to_stack_node[ed.n1.id], self.node_to_stack_node[ed.n2.id]) not in self.edge_ids:
				self.add_edge(self.node_to_stack_node[ed.n1.id], self.node_to_stack_node[ed.n2.id], weight=0)
			self.edge_ids[self.node_to_stack_node[ed.n1.id], self.node_to_stack_node[ed.n2.id]].weight += 1
