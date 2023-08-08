import numpy
from scipy.optimize import linprog
import numpy as np
import itertools, time
from src import motifs
from src.optimization import LayeredOptimizer


class HiGHSLayeredOptimizer(LayeredOptimizer):
	def __init__(self, layered_graph, parameters=None):
		LayeredOptimizer.__init__(self, layered_graph, parameters=parameters)
		if self.heuristic_start:
			print("Heuristic start only available with Gurobi implementation")
		del self.heuristic_start
		if self.xvar_branch_priority:
			print("Branching priority only available with Gurobi implementation")
		del self.xvar_branch_priority
		if self.vertical_transitivity:
			self.stratisfimal_y_vars = True

	def optimize_layout(self):
		out = self.__optimize_layout_standard_highs()
		if self.verbose:
			for string in self.print_info:
				print(string)
			self.print_info.clear()
		return out

	# TODO: implement. https://docs.scipy.org/doc/scipy/reference/optimize.linprog-highs.html#optimize-linprog-highs
	def __optimize_layout_standard_highs(self):
		g = self.g
		if not self.direct_transitivity and not self.vertical_transitivity:
			self.vertical_transitivity = True
		t1 = time.time()

		""" Collapse valid subgraphs """
		if self.collapse_subgraphs:
			g = g.collapse_ap_cases()
			self.x_var_assign = {x_v: 2 for n_l in g.get_names_by_layer().values() for x_v in itertools.combinations(n_l, 2)}

		nodes_by_layer = g.get_names_by_layer()
		edges_by_layer = g.get_edge_names_by_layer()

		""" Variables, constraints setup """
		x_vars = {}  # pointer to column in constraint matrix
		x_vars_eq = {}
		z_vars = {}
		y_vars = {}
		c_vars = {}
		c_vars_eq = {}
		c_var_constants = {}
		counter = 0
		counter_eq = 0
		x_var_counter_start = 0
		for name_list in nodes_by_layer.values():
			if self.mirror_vars:
				for x_pair in itertools.permutations(name_list, 2):
					x_vars[x_pair] = counter
					x_vars_eq[x_pair] = counter
					counter += 1
					counter_eq += 1
			else:
				for x_pair in itertools.combinations(name_list, 2):
					x_vars[x_pair] = counter
					counter += 1
		x_var_counter_end, z_var_counter_start = counter, counter
		for name_list in nodes_by_layer.values():
			if self.vertical_transitivity:
				for z_pair in itertools.permutations(name_list, 2):
					z_vars[z_pair] = counter
					counter += 1
		z_var_counter_end, c_var_counter_start = counter, counter
		for i, edge_list in edges_by_layer.items():
			if self.mirror_vars:
				for pr in itertools.permutations(edge_list, 2):
					if pr[0][0] != pr[1][0] and pr[0][1] != pr[1][1] and pr[0][0] != pr[1][1] and pr[0][1] != pr[1][0]:
						c_vars[pr] = counter
						c_var_constants[pr] = g.edge_names[pr[0]].weight * g.edge_names[pr[1]].weight
						c_vars_eq[pr] = counter
						counter += 1
			else:
				for pr in itertools.combinations(edge_list, 2):
					if pr[0][0] != pr[1][0] and pr[0][1] != pr[1][1] and pr[0][0] != pr[1][1] and pr[0][1] != pr[1][0]:
						c_vars[pr] = counter
						c_var_constants[pr] = g.edge_names[pr[0]].weight * g.edge_names[pr[1]].weight
						counter += 1
		c_var_counter_end, y_var_counter_start = counter, counter
		if self.vertical_transitivity or self.stratisfimal_y_vars:
			for nd in g.nodes:
				y_vars[nd.name] = counter
				counter += 1
		y_var_counter_end = counter
		model_constraints = []

		""" Butterfly reduction """
		butterfly_c_vars = set()
		butterfly_c_pairs = []
		if self.butterfly_reduction:
			b_set_list = []
			for b_v in motifs.get_butterflies(g):
				b_set_list.append(set(b_v))
			b_set_one_found = [False] * len(b_set_list)
			print("Butterfly set:", b_set_list)
			if b_set_list:
				for c_var in c_vars:
					c_set = {c_var[0][0], c_var[0][1], c_var[1][0], c_var[1][1]}
					if c_set in b_set_list:
						b_ind = b_set_list.index(c_set)
						if b_set_one_found[b_ind]:
							c_org = [c_v for c_v in butterfly_c_vars if {c_v[0][0], c_v[0][1], c_v[1][0], c_v[1][1]} == c_set][0]
							butterfly_c_pairs.append((c_org, c_var))
						else:
							b_set_one_found[b_ind] = True
						butterfly_c_vars.add(c_var)

		""" Transitivity constraints """
		# if self.direct_transitivity:
		#
		# elif self.vertical_transitivity:

		""" Long-version crossing reduction code """

		""" Symmetry constraints """

		""" Fix key x-var """
		# if self.symmetry_breaking:

		""" Non-sequential bendiness reduction, original Stratisfimal version"""
		# if not self.sequential_bendiness and self.bendiness_reduction:

		""" Create model objective """
		# model_objective = []
		# if not self.sequential_bendiness:
		# 	if self.bendiness_reduction:
		# 		b_vars = {}
		# 		for ed in g.edge_names:
		# 			b_vars[ed] = counter
		# 			counter += 1
		# 		b_vars = list(g.edge_names.keys())
		# 		b = m.addVars(b_vars, vtype=GRB.INTEGER, lb=0, ub=self.m_val, name="b")
		# 		m.setObjective(self.gamma_1 * c.sum() + self.gamma_2 * b.sum(), GRB.MINIMIZE)
		# 	else:
		# 		opt = gp.LinExpr()
		# 		for i, c_var in enumerate(c_vars):
		# 			opt += c_consts[i] * c[c_var]
		# 		m.setObjective(opt, GRB.MINIMIZE)
		# else:
		# 	opt = gp.LinExpr()
		# 	if self.mirror_vars and self.symmetry_constraints:
		# 		for i, c_var in enumerate(c_vars_orig):
		# 			opt += nc_consts[i] * c[c_var]
		# 	else:
		# 		for i, c_var in enumerate(c_vars):
		# 			opt += c_consts[i] * c[c_var]
		# 	m.setObjective(opt, GRB.MINIMIZE)

		""" Optimize model """

		""" Sequential bendiness reduction """

		""" Draw optimized graph """
		# if self.draw_graph:

		""" Add print info """
		# if self.verbose:
		# 	self.print_info.append(f"{pre_sym}Final edge crossing count: {g.num_edge_crossings()}")
		# 	self.print_info.append(f"{pre_sym}Number of constraints: {n_constraints_generated}")
		# 	self.print_info.append(f"{pre_sym}{round(t1, 3)}, {round(t2, 3)}, {round(t3, 3)}, {round(t1 + t2 + t3, 3)}")

		# print(f"Number of crossings: {round(m.objVal) if m.objVal != float('inf') else m.objVal}",
		# 	  f"\tOptimization time: {round(m.runtime, 3)}")

		""" Return data """
		n_cr = 0
		return n_cr

	def __edge_crossings_highs(self, sparse_matrix):
		pass  # TODO
