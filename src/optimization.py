import itertools, re, time
import random
import gurobipy as gp
from gurobipy import GRB
from sklearn.cluster import SpectralClustering
from src import vis, reductions, motifs, type_conversions
from src.graph import *


class LayeredOptimizer:
	def __init__(self, layered_graph: LayeredGraph, parameters: dict):
		self.g = layered_graph
		self.x_var_assign = {x_v: 2 for n_l in self.g.get_names_by_layer().values() for x_v in itertools.combinations(n_l, 2)}
		self.bendiness_reduction = parameters["bendiness_reduction"] if "bendiness_reduction" in parameters else False
		self.gamma_1 = parameters["gamma_1"] if "gamma_1" in parameters else 1
		self.gamma_2 = parameters["gamma_2"] if "gamma_2" in parameters else 1
		self.m_val = parameters["self.m_val"] if "self.m_val" in parameters else 50
		self.sequential_bendiness = parameters["sequential_bendiness"] if "sequential_bendiness" in parameters else True
		self.cutoff_time = parameters["cutoff_time"] if "cutoff_time" in parameters else 0
		self.do_subg_reduction = parameters["do_subg_reduction"] if "do_subg_reduction" in parameters else False
		self.return_x_vars = parameters["return_x_vars"] if "return_x_vars" in parameters else False
		self.butterfly_reduction = parameters["butterfly_reduction"] if "butterfly_reduction" in parameters else False
		self.verbose = parameters["verbose"] if "verbose" in parameters else True
		self.subg_verbose = parameters["subg_verbose"] if "subg_verbose" in parameters else True
		self.name = parameters["name"] if "name" in parameters else "graph1"
		self.print_info = []

	def sequential_br(self, graph_arg=None, substitute_x_vars=None, subgraph_seq=False):
		g = self.g if graph_arg is None else graph_arg
		x_var_opt = self.x_var_assign if substitute_x_vars is None else substitute_x_vars
		y_vars = list(g.node_names)
		n_constr = 0
		m2 = gp.Model()
		y = m2.addVars(y_vars, vtype=GRB.CONTINUOUS, lb=0, ub=self.m_val, name="y")
		# y = m2.addVars(y_vars, vtype=GRB.INTEGER, lb=0, ub=mv, name="y")
		m2.update()
		for v in m2.getVars():
			v.start = g[int(v.varName[2:v.varName.index(']')])].y
		b_vars = list(g.edge_names.keys())
		b = m2.addVars(b_vars, vtype=GRB.CONTINUOUS, lb=0, ub=self.m_val, name="b")
		# b = m2.addVars(b_vars, vtype=GRB.INTEGER, lb=0, ub=mv, name="b")
		m2.setObjective(b.sum(), GRB.MINIMIZE)
		n_orig = sum((1 for node in g.nodes if not node.is_anchor_node))
		if subgraph_seq:
			n_orig = max((n.name for n in g.nodes))
		for var, val in x_var_opt.items():
			if val == 0:
				if int(var[0]) > n_orig and int(var[1]) > n_orig:
					m2.addConstr(y[var[0]] >= 0.15 + y[var[1]], f"vert{var}")
				elif int(var[0]) > n_orig or int(var[1]) > n_orig:
					m2.addConstr(y[var[0]] >= 0.3 + y[var[1]], f"vert{var}")
				else:
					m2.addConstr(y[var[0]] >= 1 + y[var[1]], f"vert{var}")
			else:
				if int(var[0]) > n_orig and int(var[1]) > n_orig:
					m2.addConstr(y[var[0]] + 0.15 <= y[var[1]], f"vert{var}")
				elif int(var[0]) > n_orig or int(var[1]) > n_orig:
					m2.addConstr(y[var[0]] + 0.3 <= y[var[1]], f"vert{var}")
				else:
					m2.addConstr(y[var[0]] + 1 <= y[var[1]], f"vert{var}")
		for b_var in b_vars:
			m2.addConstr(y[b_var[0]] - y[b_var[1]] <= b[b_var], f"1bend{b_var}")
			m2.addConstr(y[b_var[1]] - y[b_var[0]] <= b[b_var], f"2bend{b_var}")
			n_constr += 2
		m2.setParam("OutputFlag", 0)
		m2.optimize()
		for v in m2.getVars():
			if v.varName[:1] == "y":
				g[int(v.varName[2:v.varName.index(']')])].y = float(v.x)
		return n_constr

	def edge_crossings(self, model: gp.Model, c_vars, x, c, graph_arg=None, track_x_var_usage=False):
		g = self.g if graph_arg is None else graph_arg
		n_constr_0, n_constr_1, n_constr_2 = 0, 0, 0
		x_var_usage = {}
		for c_var in c_vars:
			x1_rev, x2_rev, x3_rev = 1, 1, 1
	
			# simple edge crossing
			if not g.edge_names[c_var[0]].same_layer_edge and not g.edge_names[c_var[1]].same_layer_edge:
				if (c_var[0][0], c_var[1][0]) in x:
					x1 = (c_var[0][0], c_var[1][0])
				else:
					x1 = (c_var[1][0], c_var[0][0])
					x1_rev = -1
				if (c_var[0][1], c_var[1][1]) in x:
					x2 = (c_var[0][1], c_var[1][1])
				else:
					x2 = (c_var[1][1], c_var[0][1])
					x2_rev = -1
				if track_x_var_usage:
					if x1 not in x_var_usage:
						x_var_usage[x1] = 0
					if x2 not in x_var_usage:
						x_var_usage[x2] = 0
					x_var_usage[x1] += 2
					x_var_usage[x2] += 2
				# print(f"(1-{x1_rev}*x[{x1}]) + {x2_rev}*x[{x2}] + c[c_var] + {(1-x1_rev)/2} >= 1")
				model.addConstr((1 - x1_rev * x[x1]) + x2_rev * x[x2] + c[c_var] - (1 - x1_rev) / 2 + (1 - x2_rev) / 2 >= 1, f"1se{c_var}")
				model.addConstr(x1_rev * x[x1] + (1 - x2_rev * x[x2]) + c[c_var] + (1 - x1_rev) / 2 - (1 - x2_rev) / 2 >= 1, f"2se{c_var}")
				n_constr_0 += 2
	
			# same-layer/2-layer edge crossing
			elif g.edge_names[c_var[0]].same_layer_edge and not g.edge_names[c_var[1]].same_layer_edge:
				if (c_var[0][0], c_var[1][0]) in x:
					x1 = (c_var[0][0], c_var[1][0])
				else:
					x1 = (c_var[1][0], c_var[0][0])
					x1_rev = -1
				if (c_var[0][1], c_var[1][0]) in x:
					x2 = (c_var[0][1], c_var[1][0])
				else:
					x2 = (c_var[1][0], c_var[0][1])
					x2_rev = -1
				if track_x_var_usage:
					if x1 not in x_var_usage:
						x_var_usage[x1] = 0
					if x2 not in x_var_usage:
						x_var_usage[x2] = 0
					x_var_usage[x1] += 2
					x_var_usage[x2] += 2
				model.addConstr(-x1_rev * x[x1] + x2_rev * x[x2] + c[c_var] - (1 - x1_rev) / 2 + (1 - x2_rev) / 2 >= 0, f"1hy{c_var}")
				model.addConstr(x1_rev * x[x1] - x2_rev * x[x2] + c[c_var] + (1 - x1_rev) / 2 - (1 - x2_rev) / 2 >= 0, f"2hy{c_var}")
				n_constr_1 += 2
	
			elif g.edge_names[c_var[1]].same_layer_edge and not g.edge_names[c_var[0]].same_layer_edge:
				if (c_var[1][0], c_var[0][0]) in x:
					x1 = (c_var[1][0], c_var[0][0])
				else:
					x1 = (c_var[0][0], c_var[1][0])
					x1_rev = -1
				if (c_var[1][1], c_var[0][0]) in x:
					x2 = (c_var[1][1], c_var[0][0])
				else:
					x2 = (c_var[0][0], c_var[1][1])
					x2_rev = -1
				if track_x_var_usage:
					if x1 not in x_var_usage:
						x_var_usage[x1] = 0
					if x2 not in x_var_usage:
						x_var_usage[x2] = 0
					x_var_usage[x1] += 2
					x_var_usage[x2] += 2
				model.addConstr(-x1_rev * x[x1] + x2_rev * x[x2] + c[c_var] - (1 - x1_rev) / 2 + (1 - x2_rev) / 2 >= 0, f"1hy{c_var}")
				model.addConstr(x1_rev * x[x1] - x2_rev * x[x2] + c[c_var] + (1 - x1_rev) / 2 - (1 - x2_rev) / 2 >= 0, f"2hy{c_var}")
				n_constr_1 += 2
	
			# same layer edge crossing
			else:
				x1_rev, x2_rev, x3_rev, x4_rev = 1, 1, 1, 1
				if (c_var[0][0], c_var[1][0]) in x:
					x1 = (c_var[0][0], c_var[1][0])
				else:
					x1 = (c_var[1][0], c_var[0][0])
					x1_rev = -1
				if (c_var[0][1], c_var[1][1]) in x:
					x2 = (c_var[0][1], c_var[1][1])
				else:
					x2 = (c_var[1][1], c_var[0][1])
					x2_rev = -1
				if (c_var[0][0], c_var[1][1]) in x:
					x3 = (c_var[0][0], c_var[1][1])
				else:
					x3 = (c_var[1][1], c_var[0][0])
					x3_rev = -1
				if (c_var[0][1], c_var[1][0]) in x:
					x4 = (c_var[0][1], c_var[1][0])
				else:
					x4 = (c_var[1][0], c_var[0][1])
					x4_rev = -1
				if track_x_var_usage:
					if x1 not in x_var_usage:
						x_var_usage[x1] = 0
					if x2 not in x_var_usage:
						x_var_usage[x2] = 0
					if x3 not in x_var_usage:
						x_var_usage[x3] = 0
					if x4 not in x_var_usage:
						x_var_usage[x4] = 0
					x_var_usage[x1] += 6
					x_var_usage[x2] += 6
					x_var_usage[x3] += 6
					x_var_usage[x4] += 6
				model.addConstr(
					c[c_var] + 1 - x1_rev * x[x1] + 1 - x4_rev * x[x4] + 1 - x2_rev * x[x2] - (1 - x1_rev) / 2 - (
							1 - x4_rev) / 2 - (1 - x2_rev) / 2 >= 1, f"1sl{c_var}")
				model.addConstr(c[c_var] + 1 - x3_rev * x[x3] + x2_rev * x[x2] + x4_rev * x[x4] - (1 - x3_rev) / 2 + (
						1 - x2_rev) / 2 + (1 - x4_rev) / 2 >= 1, f"2sl{c_var}")
				model.addConstr(c[c_var] + x1_rev * x[x1] + 1 - x3_rev * x[x3] + x2_rev * x[x2] + (1 - x1_rev) / 2 - (
						1 - x3_rev) / 2 + (1 - x2_rev) / 2 >= 1, f"3sl{c_var}")
				model.addConstr(c[c_var] + 1 - x4_rev * x[x4] + 1 - x2_rev * x[x2] + x3_rev * x[x3] - (1 - x4_rev) / 2 - (
						1 - x2_rev) / 2 + (1 - x3_rev) / 2 >= 1, f"4sl{c_var}")
				model.addConstr(c[c_var] + x4_rev * x[x4] + x1_rev * x[x1] + 1 - x3_rev * x[x3] + (1 - x4_rev) / 2 + (
						1 - x1_rev) / 2 - (1 - x3_rev) / 2 >= 1, f"5sl{c_var}")
				model.addConstr(c[c_var] + 1 - x2_rev * x[x2] + x3_rev * x[x3] + 1 - x1_rev * x[x1] - (1 - x2_rev) / 2 + (
						1 - x3_rev) / 2 - (1 - x1_rev) / 2 >= 1, f"6sl{c_var}")
				model.addConstr(c[c_var] + x2_rev * x[x2] + x4_rev * x[x4] + x1_rev * x[x1] + (1 - x2_rev) / 2 + (
						1 - x4_rev) / 2 + (1 - x1_rev) / 2 >= 1, f"7sl{c_var}")
				model.addConstr(c[c_var] + x3_rev * x[x3] + 1 - x1_rev * x[x1] + 1 - x4_rev * x[x4] + (1 - x3_rev) / 2 - (
						1 - x1_rev) / 2 - (1 - x4_rev) / 2 >= 1, f"8sl{c_var}")
				n_constr_2 += 8
		return n_constr_0, n_constr_1, n_constr_2, x_var_usage

	def compute_variable_assignments(self, x_vars, c_vars):
		x_assignments = {}
		c_assignments = {}
		for x_var in x_vars:
			if self.g[x_var[0]].y < self.g[x_var[1]].y:
				x_assignments[x_var] = 1
			else:
				x_assignments[x_var] = 0
		for c_var in c_vars:
			sl1 = self.g.edge_names[c_var[0]].same_layer_edge
			sl2 = self.g.edge_names[c_var[1]].same_layer_edge
			x1 = get_x_var(x_assignments, c_var[0][0], c_var[1][0])
			if sl1 and sl2:
				x2 = get_x_var(x_assignments, c_var[0][1], c_var[1][1])
				x3 = get_x_var(x_assignments, c_var[1][0], c_var[0][1])
				x4 = get_x_var(x_assignments, c_var[0][0], c_var[1][1])
				c_assignments[c_var] = x1 * x2 * x3 + x4 * (1 - x2) * (1 - x3) + (1 - x3) * (1 - x1) * x4 + x2 * (1 - x4) * x1 + (1 - x1) * x4 * (1 - x2) + x3 * x2 * (1 - x4) + (1 - x4) * x1 * x3 + (1 - x2) * (1 - x3) * (1 - x1)
			elif sl1:
				x2 = get_x_var(x_assignments, c_var[0][1], c_var[1][0])
				c_assignments[c_var] = x1 * (1 - x2) + (1 - x1) * x2
			elif sl2:
				x2 = get_x_var(x_assignments, c_var[0][0], c_var[1][1])
				c_assignments[c_var] = x1 * (1 - x2) + (1 - x1) * x2
			else:
				x2 = get_x_var(x_assignments, c_var[0][1], c_var[1][1])
				c_assignments[c_var] = x1 * (1 - x2) + (1 - x1) * x2
		return x_assignments, c_assignments

	def optimize_with_subgraph_reduction(self, n_partitions):
		self.print_info.append("")

		cluster = SpectralClustering(n_clusters=n_partitions, assign_labels="discretize", affinity="precomputed").fit(self.g.adjacency_matrix())
		subgraphs = [set(node.name for node in self.g.nodes if cluster.labels_[node.name - 1] == i) for i in range(n_partitions)]
		top_level_g, crosses, contacts, stack_to_nodeset, node_to_stack = self.g.stacked_graph_from_subgraph_nodes(cluster.labels_)
		top_level_optval, top_x_vars = self.optimize_layout_standard(graph_arg=top_level_g, return_x_vars=True, bendiness_reduction=False, is_subgraph=True, name="Collapsed graph", verbose=True)
		# TODO choose "clean" subgraph, merge the others
		t = time.time()
		vis.draw(top_level_g, "interim", gravity=True)
		vis.draw(self.g, "overall", groups=cluster.labels_)
		for x_var in top_x_vars:
			for low_n1 in stack_to_nodeset[x_var[0]]:
				for low_n2 in stack_to_nodeset[x_var[1]]:
					set_x_var(self.x_var_assign, low_n1, low_n2, top_x_vars[x_var])
		# x_var_assig = {}
		# for x_var in top_x_vars:
		# 	if next(iter(stack_to_nodeset[x_var[0]])) in subgraphs[0] or next(iter(stack_to_nodeset[x_var[1]])) in subgraphs[0]:
		# 		for low_n1 in stack_to_nodeset[x_var[0]]:
		# 			for low_n2 in stack_to_nodeset[x_var[1]]:
		# 				set_x_var(x_var_assig, low_n1, low_n2, top_x_vars[x_var])
		# self.optimize_with_starting_assignments(x_var_assig)
		# return

		sides = {}
		for contact_node in crosses:
			# if len(crosses[contact_node] > 1    <- TODO test if node has multiple contact edges
			# figure out contact_sides to pass to call of optimize. contact_sides=1 => node fixed at top
			for contact_other in crosses[contact_node]:
				if node_to_stack[contact_node] + 1 in stack_to_nodeset and cluster.labels_[next(iter(stack_to_nodeset[node_to_stack[contact_node] + 1])) - 1] == cluster.labels_[next(iter(stack_to_nodeset[node_to_stack[contact_node]])) - 1]:
					sides[contact_node] = 1 - get_x_var(top_x_vars, node_to_stack[contact_node] + 1, node_to_stack[contact_other])
				if node_to_stack[contact_other] - 1 in stack_to_nodeset and cluster.labels_[next(iter(stack_to_nodeset[node_to_stack[contact_other] - 1])) - 1] == cluster.labels_[next(iter(stack_to_nodeset[node_to_stack[contact_other]])) - 1]:
					sides[contact_other] = 1 - get_x_var(top_x_vars, node_to_stack[contact_other] - 1, node_to_stack[contact_node])

		self.g.adj_list = self.g.create_double_adj_list(forward_only=True)
		layered_subg_list = []
		x_var_opt_list = []
		unconstrained_opt_vals = []
		opt_vals = []
		for i, subg in enumerate(subgraphs):
			g_prime = LayeredGraph()
			layered_subg_list.append(g_prime)
			sides_prime = {}
			vars_to_fix = {}

			for subg_node in subg:
				g_prime.add_node(self.g[subg_node].layer, name=self.g[subg_node].name, is_anchor=self.g[subg_node].is_anchor_node)
			for subg_node in subg:
				if subg_node in sides:
					sides_prime[subg_node] = sides[subg_node]
				for adj_node in self.g.adj_list[subg_node]:
					if cluster.labels_[subg_node - 1] == cluster.labels_[adj_node - 1]:
						g_prime.add_edge(subg_node, adj_node)
			gp_layers = g_prime.get_names_by_layer()
			for contact_node, x_val in sides_prime.items():
				for node in gp_layers[g_prime[contact_node].layer]:
					if node != contact_node and len(self.g.adj_list[node]) >= 1:
						vars_to_fix[(contact_node, node)] = x_val

			unconstrained_opt_val = self.optimize_layout_standard(graph_arg=g_prime, bendiness_reduction=False)
			opt_val, x_vars_opt = self.optimize_layout_standard(graph_arg=g_prime, bendiness_reduction=False, is_subgraph=True, fix_x_vars=vars_to_fix, return_x_vars=True, name=f"subgraph {i+1}", verbose=self.subg_verbose)
			if unconstrained_opt_val != opt_val:
				self.print_info.append(f"\tSubgraph {i+1} has {opt_val} crossings but could have as few as {unconstrained_opt_val}")
			else:
				self.print_info.append(f"\tSubgraph {i+1} is optimal")
			x_var_opt_list.append(x_vars_opt)
			unconstrained_opt_vals.append(unconstrained_opt_val)
			opt_vals.append(opt_val)
			for x_var, val in x_vars_opt.items():
				set_x_var(self.x_var_assign, x_var[0], x_var[1], val)

			vis.draw(g_prime, f"interim_subg{i+1}")

		t = time.time() - t

		self.print_info.append(f"Total time to optimize and patch together subgraphs: {t}")
		self.print_info.append(f"Final edge crossing count: {self.g.num_edge_crossings()}")
		self.print_info.append(f"{top_level_optval} crossings in collapsed graph, {sum(opt_vals)} in subgraphs, {self.g.num_edge_crossings() - top_level_optval - sum(opt_vals)} from crossing edges")

		vis.draw(self.g, "endpoints_highlight", groups=cluster.labels_)

		# TODO return x_vars (or use self.x_var_assign?), n_crossings, subgraph LayeredGraph objects
		return 0

	def optimize_layout_standard(self, graph_arg=None, bendiness_reduction=True, assignment=None, return_x_vars=False, pool_solutions=False, name="graph1", fix_x_vars=None, start_x_vars=None, is_subgraph=False, verbose=False, use_top_level_params=False):
		g = self.g if graph_arg is None else graph_arg

		t1 = time.time()
		nodes_by_layer = g.get_names_by_layer()
		edges_by_layer = g.get_edge_names_by_layer()
		n_constraints_generated = [0] * 6  # simple edge, hybrid edge, same layer edge, vertical pos, bendiness, total
		pre_sym = '\t' if is_subgraph else ''
		if verbose or (use_top_level_params and self.verbose):
			self.print_info.extend(('-' * 70, f"{self.name if use_top_level_params else name}:", '-' * 70))
		m = gp.Model()
	
		""" Add all variables """
		x_vars = []
		x_vars_layers = {}
		# z_vars = []
		for i, name_list in nodes_by_layer.items():
			x_vars += list(itertools.combinations(name_list, 2))
			# z_vars += list(itertools.permutations(name_list, 2))
			x_vars_layers[i] = list(itertools.combinations(name_list, 2))
		# x = m.addVars(x_vars, vtype=GRB.CONTINUOUS, name="x")
		# k = m.addVars(x_vars, vtype=GRB.CONTINUOUS, name="k")
		x = m.addVars(x_vars, vtype=GRB.BINARY, name="x")
		# z = m.addVars(z_vars, vtype=GRB.CONTINUOUS, lb=0, ub=self.m_val, name="z")
		c_vars, c_consts = reductions.normal_c_vars(g, edges_by_layer)
		c = m.addVars(c_vars, vtype=GRB.CONTINUOUS, name="c")
		y_vars = [n.name for n in g]
		y = m.addVars(y_vars, vtype=GRB.CONTINUOUS, lb=0, ub=self.m_val, name="y")
		m.update()
	
		""" Fix variables """
		if assignment:  # TODO remove in place of fix_x_vars
			for k, v in assignment.items():
				m.getVarByName(k).ub, m.getVarByName(k).lb = v, v

		if fix_x_vars:
			for k, v in fix_x_vars.items():
				if k in x:
					m.getVarByName(f"x[{k[0]},{k[1]}]").lb, m.getVarByName(f"x[{k[0]},{k[1]}]").ub = v, v
				else:
					m.getVarByName(f"x[{k[1]},{k[0]}]").lb, m.getVarByName(f"x[{k[1]},{k[0]}]").ub = 1 - v, 1 - v

		""" Starting variable assignment """
		if start_x_vars:
			for v in m.getVars():
				v.Start = start_x_vars[v.varName]
			# for k, v in start_x_vars.items():
			# 	if k in x:
			# 		m.getVarByName(f"x[{k[0]},{k[1]}]").Start = v
			# 	else:
			# 		m.getVarByName(f"x[{k[1]},{k[0]}]").Start = 1 - v

		g_igraph = type_conversions.layered_graph_to_igraph(g)
		heuristic_layout = g_igraph.layout_sugiyama(layers=g_igraph.vs["layer"])
		print(heuristic_layout)

		# g.barycentric_reordering(10)  TODO igraph version
		# x_assign, c_assign = self.compute_variable_assignments(g, x_vars, c_vars)
		# for var in m.getVars():
		# 	if var.varName[:1] == "y":
		# 		var.start = g[int(var.varName[2:var.varName.index(']')])].y
		# 	elif var.varName[:1] == "x":
		# 		if x_assign[int(var.varName[2:var.varName.index(',')]), int(
		# 				var.varName[var.varName.index(',') + 1:var.varName.index(']')])]:
		# 			var.start = 1
		# 		else:
		# 			var.start = 0
		# 	elif var.varName[:1] == "c":
		# 		cnm = [int(x) for x in re.findall(r"[0-9]+", var.varName)]
		# 		if c_assign[(cnm[0], cnm[1]), (cnm[2], cnm[3])]:
		# 			var.start = 1
		# 		else:
		# 			var.start = 0
	
		""" Butterfly reduction """
		if self.butterfly_reduction:
			b_set_list = []
			for b_v in motifs.get_butterflies(g):
				b_set_list.append({b_v[0][0], b_v[0][1], b_v[1][0], b_v[1][1]})
			if b_set_list:
				self.print_info.append(f"{pre_sym}Butterflies found: {len(b_set_list)}")
				for var in m.getVars():
					if var.varName[:1] == 'c':
						if {int(x) for x in re.findall(r"[0-9]+", var.varName)} in b_set_list:
							var.lb, var.ub = 1, 1
	
		""" Set model objective function """
		if not self.sequential_bendiness:
			# z_vars = []
			# for i, name_list in nodes_by_layer.items():
			# 	z_vars += list(itertools.permutations(name_list, 2))
			# z = m.addVars(z_vars, vtype=GRB.INTEGER, lb=0, ub=self.m_val, name="z")
			if bendiness_reduction:
				b_vars = list(g.edge_names.keys())
				b = m.addVars(b_vars, vtype=GRB.INTEGER, lb=0, ub=self.m_val, name="b")
				m.setObjective(self.gamma_1 * c.sum() + self.gamma_2 * b.sum(), GRB.MINIMIZE)
			else:
				opt = gp.LinExpr()
				for i, c_var in enumerate(c_vars):
					opt += c_consts[i] * c[c_var]
				m.setObjective(opt, GRB.MINIMIZE)
		else:
			opt = gp.LinExpr()
			for i, c_var in enumerate(c_vars):
				opt += c_consts[i] * c[c_var]
			m.setObjective(opt, GRB.MINIMIZE)
	
		""" Long-version crossing reduction code """
		n_cs = self.edge_crossings(m, c_vars, x, c, graph_arg=g)
		for i, val in enumerate(n_cs[:-1]):
			n_constraints_generated[i] += val

		""" Fix key x-var """
		x_var_usage = n_cs[-1]
		most_used_x = max(x_var_usage, key=x_var_usage.get)
		m.getVarByName(f"x[{most_used_x[0]},{most_used_x[1]}]").lb = 0
		m.getVarByName(f"x[{most_used_x[0]},{most_used_x[1]}]").ub = 0

		""" Vertical position, implication version """
		for x_var in x_vars:
			m.addGenConstrIndicator(x[x_var], True, y[x_var[0]] + 1 <= y[x_var[1]])
			m.addGenConstrIndicator(x[x_var], False, y[x_var[0]] >= 1 + y[x_var[1]])
			n_constraints_generated[3] += 2
	
		""" Original vertical position constraints """
		# for x_var in x_vars:
		# 	m.addConstr(z[x_var] - self.m_val * x[x_var] <= 0, f"1.1z{x_var}")
		# 	m.addConstr(z[x_var] - y[x_var[0]] - (self.m_val * x[x_var]) >= -1 * self.m_val, f"1.2z{x_var}")
		# 	m.addConstr(y[x_var[1]] - z[x_var] - x[x_var] >= 0, f"1.3z{x_var}")
		# 	m.addConstr(z[x_var] <= y[x_var[0]], f"1.4z{x_var}")
		# 	# m.addConstr(z[x_var] >= 0, f"1.5z{x_var}")	# print(f"z{x_var}-{self.m_val}*x{x_var} <= 0"), print(f"y{x_var[1]} - z{x_var} - x{x_var} >= 0")
		#
		# 	m.addConstr(z[x_var[1], x_var[0]] - self.m_val * (1 - x[x_var]) <= 0, f"2.1z{x_var}")
		# 	m.addConstr(z[x_var[1], x_var[0]] - y[x_var[1]] - self.m_val * (1 - x[x_var]) >= -1 * self.m_val,
		# 				f"2.2z{x_var}")
		# 	m.addConstr(y[x_var[0]] - z[x_var[1], x_var[0]] - (1 - x[x_var]) >= 0, f"2.3z{x_var}")
		# 	m.addConstr(z[x_var[1], x_var[0]] <= y[x_var[1]], f"2.4z{x_var}")
		# 	# m.addConstr(z[x_var[1],x_var[0]] >= 0, f"2.5z{x_var}")	# print(f"z[{x_var[1],x_var[0]}]-{self.m_val}*(1-x{x_var}) <= 0"), print(f"y{x_var[0]} - z[{x_var[1],x_var[0]}] - (1-x{x_var}) >= 0")
		#
		# 	# m.addConstr(k[x_var] + x[x_var] - 1 >= 0)  # SOS-constraint version. Use with full LP relaxation
		# 	# m.addSOS(GRB.SOS_TYPE1, [k[x_var], x[x_var]])
		#
		# 	n_constraints_generated[3] += 10
	
		""" Non-sequential bendiness reduction, original Stratisfimal version"""
		if not self.sequential_bendiness and self.bendiness_reduction:
			for b_var in b_vars:
				m.addConstr(y[b_var[0]] - y[b_var[1]] <= b[b_var], f"1bend{b_var}")
				m.addConstr(y[b_var[1]] - y[b_var[0]] <= b[b_var], f"2bend{b_var}")
				n_constraints_generated[4] += 2
	
		""" If graph is a subgraph, fix vars to ensure contact nodes on the right side """
		# if is_subgraph:
		# 	adj_list = g.create_double_adj_list(forward_only=True)
		# 	print(contact_nodes, contact_sides)
		# 	for v in x_vars:
		# 		for i, contact_nd in enumerate(contact_nodes):
		# 			if contact_sides[i] >= 0:
		# 				if v[0] == contact_nd and contact_sides[i] == 1:
		# 					if len(adj_list[v[1]]) > 0:
		# 						m.getVarByName(f"x[{v[0]},{v[1]}]").lb, m.getVarByName(f"x[{v[0]},{v[1]}]").ub = 1, 1
		# 				elif v[0] == contact_nd:
		# 					if len(g.adj_list[v[1]]) > 0:
		# 						m.getVarByName(f"x[{v[0]},{v[1]}]").lb, m.getVarByName(f"x[{v[0]},{v[1]}]").ub = 0, 0
		# 				elif v[1] == contact_nd and contact_sides[i] == 1:
		# 					if len(adj_list[v[0]]) > 0:
		# 						m.getVarByName(f"x[{v[0]},{v[1]}]").lb, m.getVarByName(f"x[{v[0]},{v[1]}]").ub = 0, 0
		# 				elif v[1] == contact_nd:
		# 					if len(adj_list[v[0]]) > 0:
		# 						m.getVarByName(f"x[{v[0]},{v[1]}]").lb, m.getVarByName(f"x[{v[0]},{v[1]}]").ub = 1, 1
	
		""" Optimize model """
		if self.cutoff_time > 0:
			m.setParam("TimeLimit", self.cutoff_time)
		m.setParam("OutputFlag", 0)  # TODO measure time for individual branches/iterations
		n_constraints_generated[5] = sum(n_constraints_generated[:5])
		t1 = time.time() - t1
		if verbose or (use_top_level_params and self.verbose):
			self.print_info.append(f"{pre_sym}Time to input constraints: {t1}")
			self.print_info.append(f"{pre_sym}Constraint counts: {n_constraints_generated}")
			# self.print_info.append(f"{pre_sym}{g.layer_counts()[0]}\n{pre_sym}{g.layer_counts()[1]}")
		t2 = time.time()
		m.optimize()
		t2 = time.time() - t2
		if verbose or (use_top_level_params and self.verbose):
			self.print_info.append(f"{pre_sym}Objective: {m.objVal}")
			self.print_info.append(f"{pre_sym}Time to optimize: {t2}")
		x_vars_opt = {}
		for v in m.getVars():
			if v.varName[:1] == 'y':
				g[int(v.varName[2:v.varName.index(']')])].y = float(v.x)
			elif v.varName[:1] == "x" and use_top_level_params:
				self.x_var_assign[int(v.varName[2:v.varName.index(',')]), int(v.varName[v.varName.index(',') + 1:v.varName.index(']')])] = round(v.x)
			elif v.varName[:1] == "x":
				x_vars_opt[int(v.varName[2:v.varName.index(',')]), int(v.varName[v.varName.index(',') + 1:v.varName.index(']')])] = round(v.x)

		""" Draw pre-bendiness graph """
		# for v in m.getVars():
		# 	print('%s %g' % (v.varName, v.x))
		# 	if v.varName[:1] == "y":
		# 		g[int(v.varName[2:v.varName.index(']')])].y = round(v.x)
		# vis.draw(g, "interim")
	
		""" Sequantial bendiness reduction """
		if bendiness_reduction or (use_top_level_params and self.bendiness_reduction):
			t3 = time.time()
			# if not do_subg_reduction:
			# 	for v in m.getVars():
			# 		# print('%s %g' % (v.varName, v.x))
			# 		if v.varName[:1] == "x":
			# 			x_vars_opt[int(v.varName[2:v.varName.index(',')]), int(v.varName[v.varName.index(',') + 1:v.varName.index(']')])] = round(v.x)
			n_constraints_generated[4] += self.sequential_br(graph_arg=g, subgraph_seq=is_subgraph)
			t3 = time.time() - t3
			if not is_subgraph and verbose:
				self.print_info.append(f"{pre_sym}Time to perform bendiness reduction: {t3}")
		else:
			t3 = 0
			for v in m.getVars():
				if v.varName[:1] == "y":
					g[int(v.varName[2:v.varName.index(']')])].y = float(v.x)

		if verbose or (use_top_level_params and self.verbose):
			self.print_info.append(f"{pre_sym}Final edge crossing count: {g.num_edge_crossings()}")
			self.print_info.append(f"{pre_sym}Number of constraints: {n_constraints_generated}")
			self.print_info.append(f"{pre_sym}{round(t1, 3)}, {round(t2, 3)}, {round(t3, 3)}, {round(t1 + t2 + t3, 3)}")

		if return_x_vars or (use_top_level_params and self.return_x_vars):
			# vars_to_assign = {v.varName: v.x for v in m.getVars()}
			vars_to_assign = {}
			for v in m.getVars():
				if v.varName[:1] == "x":
					vars_to_assign[int(v.varName[2:v.varName.index(',')]), int(v.varName[v.varName.index(',') + 1:v.varName.index(']')])] = round(v.x)
			return round(t1 + t2 + t3 + t3, 3), vars_to_assign
	
		return round(t1 + t2 + t3 + t3, 3), int(m.objVal)

	""" Return min possible num crossings, given by my formula """
	def bound_on_optimality(self, list_of_subgraphs, crossing_num_tuple):
		return

	def optimize_full_subgraph_algorithm(self, t_allotted):
		t_start = time.time()
		n_partitions = 2
		best_so_far = self.x_var_assign.copy()
		best_n_cr = self.g.num_edge_crossings()
		lower_bound = 0
		while lower_bound != best_n_cr and time.time() - t_start < t_allotted:
			n_rep = max(6 - n_partitions, 1)
			cur_rep = 1
			while cur_rep <= n_rep and lower_bound != best_n_cr and time.time() - t_start < t_allotted:
				print(f"{n_partitions} partitions attempt {cur_rep}, {time.time()-t_start} elapased.")
				result_x_vars, n_cr = self.optimize_with_subgraph_reduction(n_partitions)
				if n_cr < best_n_cr:
					best_so_far = result_x_vars
					best_n_cr = n_cr
					print(f"New best, {n_cr} crossings")
					lower_bound = self.bound_on_optimality(result_x_vars)
				cur_rep += 1
			n_partitions += 1
		if lower_bound == best_n_cr:
			# we're optimal, baby
			print("wooooo")
		else:
			print("ran out of time ;(")

	def optimize_with_starting_assignments(self, assigments):
		out = self.optimize_layout_standard(use_top_level_params=True, fix_x_vars=assigments)
		if self.verbose:
			for string in self.print_info:
				print(string)
			self.print_info.clear()
		return out

	def generate_random_vars_to_fix(self, n_vars):
		r_vars = random.sample(list(self.x_var_assign.keys()), n_vars)
		assignments = {r_var: random.randint(0, 1) for r_var in r_vars}
		return assignments

	def normalized_loss(self, x_vars, correct_solution, g1=0, g2=1):
		correct_y_solution = self.find_y_assignment_given_x(correct_solution)
		my_y_solution = self.find_y_assignment_given_x(x_vars)
		loss = 0
		for nd, y_v in correct_y_solution.items():
			if y_v != my_y_solution[nd]:
				loss += g1 + g2 * (abs(y_v - my_y_solution[nd]) - 1)
		return loss

	def find_y_assignment_given_x(self, x_vars):
		l_to_vlist = {l: [v for v in x_vars if self.g[v[0]].layer == l] for l in range(1, self.g.n_layers + 1)}
		y_assign = {}
		for v_list in l_to_vlist.values():
			counts = {x_var[0]: 0 for x_var in v_list}
			counts.update({x_var[1]: 0 for x_var in v_list})
			for v in v_list:
				if x_vars[v] == 0:
					counts[v[0]] += 1
				else:
					counts[v[1]] += 1
			in_order = sorted(list(counts.values()))
			if all((in_order[i] < in_order[i+1] for i in range(len(in_order)-1))):
				for x_v, ct in counts.items():
					y_assign[x_v] = ct
			else:
				return "non-transitive x-values"
		return y_assign

	def optimize_layout(self):
		if self.do_subg_reduction:
			out = self.optimize_with_subgraph_reduction(2)
		else:
			out = self.optimize_layout_standard(use_top_level_params=True)
		if self.verbose:
			for string in self.print_info:
				print(string)
			self.print_info.clear()
		return out
