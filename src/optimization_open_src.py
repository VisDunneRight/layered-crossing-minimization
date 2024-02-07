import random
from scipy.optimize import linprog
import itertools
import time
from src import motifs, vis
from src.optimization import LayeredOptimizer
from src.graph import LayeredGraph
from src.helpers import *


class HiGHSLayeredOptimizer(LayeredOptimizer):
	def __init__(self, layered_graph, **kwargs):
		LayeredOptimizer.__init__(self, layered_graph, **kwargs)
		if self.heuristic_start:
			print("Heuristic start only available with Gurobi implementation")
		del self.heuristic_start
		if self.xvar_branch_priority:
			print("Branching priority only available with Gurobi implementation")
		del self.xvar_branch_priority
		if self.vertical_transitivity:
			self.stratisfimal_y_vars = True
		self.unbound_val = 12345
		self.constraints = []
		self.constraints_eq = []
		self.nvars = 0
		self.b_ub = []
		self.b_eq = []
		self.b_lb = []
		self.c_t = []
		self.bounds = []
		self.integrality = []

	def optimize_layout(self, fix_x_vars=None):
		out = self.__optimize_layout_standard_highs(fix_x_vars=fix_x_vars)
		return out

	# TODO: implement. https://docs.scipy.org/doc/scipy/reference/optimize.linprog-highs.html#optimize-linprog-highs
	def __optimize_layout_standard_highs(self, fix_x_vars=None):
		if not self.vertical_transitivity and not self.direct_transitivity:
			self.direct_transitivity = True
		if self.polyhedral_constraints:
			self.claw_constraints, self.dome_path_constraints = True, True
		g = self.g
		t1 = time.time()

		""" Collapse valid subgraphs """
		if self.collapse_leaves:
			g = g.collapse_leaves()
			self.x_var_assign = {x_v: 2 for n_l in g.get_ids_by_layer().values() for x_v in itertools.combinations(n_l, 2)}

		vis.draw_graph(g, "interim")

		""" Variables, constraints setup """
		x_vars = {}  # pointer to column in constraint matrix
		# if self.mip_relax:
		# 	x_vars_eq = x_vars
		z_vars = {}
		y_vars = {}
		c_vars = {}
		# if self.mip_relax:
		# 	c_vars_eq = c_vars
		c_var_constants = {}
		counter = 0
		xs = 0
		for name_list in g.get_ids_by_layer().values():
			combinatorics = itertools.permutations(name_list, 2) if self.mirror_vars else itertools.combinations(name_list, 2)
			for x_pair in combinatorics:
				x_vars[x_pair] = counter
				# x_vars_eq[x_pair] = counter_eq
				counter += 1
				# counter_eq += 1
		xf, zs = counter, counter
		if self.vertical_transitivity:
			for name_list in g.get_ids_by_layer().values():
				for z_pair in itertools.permutations(name_list, 2):
					z_vars[z_pair] = counter
					counter += 1
		zf, cs = counter, counter
		for i, edge_list in g.get_edge_ids_by_layer().items():
			if self.mirror_vars:
				for pr in itertools.permutations(edge_list, 2):
					if pr[0][0] != pr[1][0] and pr[0][1] != pr[1][1] and pr[0][0] != pr[1][1] and pr[0][1] != pr[1][0]:
						c_vars[pr] = counter
						if pr[0][0] < pr[1][0]:
							c_var_constants[pr] = g.edge_ids[pr[0]].weight * g.edge_ids[pr[1]].weight
						else:
							c_var_constants[pr] = 0
						# c_vars_eq[pr] = counter_eq
						counter += 1
						# counter_eq += 1
			else:
				for pr in itertools.combinations(edge_list, 2):
					if pr[0][0] != pr[1][0] and pr[0][1] != pr[1][1] and pr[0][0] != pr[1][1] and pr[0][1] != pr[1][0]:
						c_vars[pr] = counter
						c_var_constants[pr] = g.edge_ids[pr[0]].weight * g.edge_ids[pr[1]].weight
						# c_vars_eq[pr] = counter_eq
						counter += 1
						# counter_eq += 1
		cf, ys = counter, counter
		if self.vertical_transitivity:
			for nd in g.nodes:
				y_vars[nd.id] = counter
				counter += 1
		yf = counter
		self.nvars = counter
		self.bounds = [(0, None)] * self.nvars
		self.integrality = [0] * self.nvars
		for z_idx in z_vars.values():
			self.bounds[z_idx] = (0, self.m_val)
		for y_idx in y_vars.values():
			self.bounds[y_idx] = (0, self.m_val)

		""" Fix x-vars as necessary, if this is a collapsed subgraph """
		if fix_x_vars is not None:
			for x_v, val in fix_x_vars.items():
				if x_v in x_vars:
					self.bounds[x_vars[x_v]] = (val, val)
				else:
					self.bounds[x_vars[x_v[1], x_v[0]]] = (val, val)

		""" Butterfly reduction """
		butterfly_c_pairs = self.get_butterfly_cvars(g, c_vars)

		""" Transitivity constraints """
		self.__transitivity_matrix(x_vars, y_vars, z_vars, g)

		""" Edge crossing constraints """
		xvar_usage = self.__edge_crossings_matrix(c_vars, x_vars, butterfly_c_pairs, track_x_var_usage=self.symmetry_breaking)

		""" 3-claw constraints """
		self.__add_3claw_constraints_matrix(g, c_vars)

		""" Dome path constraints """
		self.__add_domepath_constraints_matrix(g, x_vars, c_vars)

		""" Symmetry constraints """
		self.__symmetry_constraints_matrix(c_vars, x_vars)

		""" Cycle constraints """
		self.__cycle_constraints_matrix(g, c_vars)

		""" Fix key x-var """
		self.__symmetry_breaking_matrix(x_vars, xvar_usage)

		""" Non-sequential bendiness reduction, original Stratisfimal version"""
		# if not self.sequential_bendiness and self.bendiness_reduction:

		""" Create model objective """
		self.c_t = [0] * self.nvars
		for c_v, const in c_var_constants.items():
			self.c_t[c_vars[c_v]] = const
		if self.mip_relax:
			for v in x_vars.values():
				self.integrality[v] = 1
		else:
			for v in x_vars.values():
				self.integrality[v] = 1
			for v in c_vars.values():
				self.integrality[v] = 1
			for v in y_vars.values():
				self.integrality[v] = 1
			for v in z_vars.values():
				self.integrality[v] = 1

		if len(self.c_t) == 0:  # Add dummy variables if there are none, to prevent errors with linprog
			self.__add_late_breaking_variable(integral=False)
			self.__add_late_breaking_variable(integral=False)
		if len(self.constraints) == 0:  # Add dummy constraint if there are none to prevent error with linprog
			self.__add_matrix_constraint([], [], 0)

		# print(x_vars)
		# print(z_vars)
		# print(c_vars)
		# print(y_vars)
		# print(self.constraints)
		# print(self.b_ub)
		# print(self.c_t)
		# print(self.integrality)

		""" Optimize model """
		t1 = time.time() - t1
		t2 = time.time()
		options = {"time_limit": self.cutoff_time} if self.cutoff_time > 0 else {}
		res = linprog(self.c_t, method='highs', A_ub=self.constraints, A_eq=self.constraints_eq if self.constraints_eq != [] else None, b_ub=self.b_ub, b_eq=self.b_eq, bounds=self.bounds, integrality=self.integrality, options=options)
		t2 = time.time() - t2
		n_cr = round(res.fun) if res.fun is not None else float('inf')
		if res.x is not None:
			for xv, idx in x_vars.items():
				self.x_var_assign[xv] = round(res.x[idx])
		else:
			for xv, idx in x_vars.items():
				self.x_var_assign[xv] = 0 if xv in self.x_var_assign else 1

		""" Optimize and merge collapsed subgraphs """
		g, t1, num_crossings = self.__optimize_subgraphs_matrix(g, x_vars, t1, n_cr)

		""" Sequential bendiness reduction """
		self.__sequential_bendiness_matrix()

		""" Draw optimized graph """
		if self.draw_graph:
			if not self.bendiness_reduction:
				if self.direct_transitivity:
					g.assign_y_vals_given_x_vars(self.x_var_assign)
				else:
					for y_v, y_idx in y_vars.items():
						g[y_v].y = res.x[y_idx]
			vis.draw_graph(g, self.name)

		""" Print output """
		print(res.message)
		print(f"Number of crossings: {n_cr}", f"\tOptimization time: {round(t2, 3)}")

		""" Return data """
		if self.return_experiment_data:
			return len(x_vars), len(c_vars), self.nvars, len(self.constraints), n_cr, t2, bool(res.success), int(res.nit), round(t1, 3)

		return round(t1 + t2, 3), n_cr

	def __add_matrix_constraint(self, col_idxs, coefficients, upper_bound, equality=False):
		constr = [0] * self.nvars
		for i in range(len(col_idxs)):
			constr[col_idxs[i]] = coefficients[i]
		if equality:
			# if self.mip_relax:
			# 	self.b_lb.append(upper_bound)
			# 	self.constraints.append(constr)
			# else:
			self.b_eq.append(upper_bound)
			self.constraints_eq.append(constr)
		else:
			# if self.mip_relax:
			# 	self.b_lb.append(-np.inf)
			self.b_ub.append(upper_bound)
			self.constraints.append(constr)

	def __add_late_breaking_variable(self, bound=None, integral=True):
		self.nvars += 1
		self.bounds.append(bound if bound is not None else (0, None))
		self.integrality.append(1 if integral else 0)
		self.c_t.append(0)
		for constr in self.constraints:
			constr.append(0)
		for constr in self.constraints_eq:
			constr.append(0)

	def __edge_crossings_matrix(self, c, x, butterflies, track_x_var_usage=False):
		x_var_usage = {}
		for c_var_pair in butterflies:
			self.__add_matrix_constraint([c[c_var_pair[0]], c[c_var_pair[1]]], [1, 1], 1, equality=True)
			# model.addConstr(c[c_var_pair[0]] + c[c_var_pair[1]] == 1)
		for c_var in c:
			if self.mirror_vars:
				self.__add_matrix_constraint([x[c_var[1][0], c_var[0][0]], x[c_var[0][1], c_var[1][1]], c[c_var]], [-1, -1, -1], -1)
				# model.addConstr(x[c_var[1][0], c_var[0][0]] + x[c_var[0][1], c_var[1][1]] + c[c_var] >= 1, f"1se{c_var}")
				self.__add_matrix_constraint([x[c_var[0][0], c_var[1][0]], x[c_var[1][1], c_var[0][1]], c[c_var]], [-1, -1, -1], -1)
				# model.addConstr(x[c_var[0][0], c_var[1][0]] + x[c_var[1][1], c_var[0][1]] + c[c_var] >= 1, f"2se{c_var}")
				if track_x_var_usage:
					if (c_var[1][0], c_var[0][0]) not in x_var_usage:
						x_var_usage[c_var[1][0], c_var[0][0]] = 0
					if (c_var[0][1], c_var[1][1]) not in x_var_usage:
						x_var_usage[c_var[0][1], c_var[1][1]] = 0
					if (c_var[0][0], c_var[1][0]) not in x_var_usage:
						x_var_usage[c_var[0][0], c_var[1][0]] = 0
					if (c_var[1][1], c_var[0][1]) not in x_var_usage:
						x_var_usage[c_var[1][1], c_var[0][1]] = 0
					x_var_usage[c_var[1][0], c_var[0][0]] += 1
					x_var_usage[c_var[0][1], c_var[1][1]] += 1
					x_var_usage[c_var[0][0], c_var[1][0]] += 1
					x_var_usage[c_var[1][1], c_var[0][1]] += 1
			else:
				x1_rev, x2_rev, x3_rev = 1, 1, 1
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
				self.__add_matrix_constraint([x[x1], x[x2], c[c_var]], [x1_rev, -x2_rev, -1], -(1 - x1_rev) / 2 + (1 - x2_rev) / 2)
				# model.addConstr((1 - x1_rev * x[x1]) + x2_rev * x[x2] + c[c_var] - (1 - x1_rev) / 2 + (1 - x2_rev) / 2 >= 1,f"1se{c_var}")
				self.__add_matrix_constraint([x[x1], x[x2], c[c_var]], [-x1_rev, x2_rev, -1], (1 - x1_rev) / 2 - (1 - x2_rev) / 2)
				# model.addConstr(x1_rev * x[x1] + (1 - x2_rev * x[x2]) + c[c_var] + (1 - x1_rev) / 2 - (1 - x2_rev) / 2 >= 1, f"2se{c_var}")
		return x_var_usage

	def __transitivity_matrix(self, x, y, z, g):
		if self.direct_transitivity:
			for x_vars_list in g.get_ids_by_layer().values():
				for x_1, x_2, x_3 in itertools.combinations(x_vars_list, 3):
					if self.mirror_vars:
						self.__add_matrix_constraint([x[x_1, x_2], x[x_2, x_3], x[x_1, x_3]], [-1, -1, 1], 0)
						# model.addConstr(x[x_1, x_2] + x[x_2, x_3] - x[x_1, x_3] >= 0)
						self.__add_matrix_constraint([x[x_1, x_3], x[x_1, x_2], x[x_2, x_3]], [-1, 1, 1], 1)
						# model.addConstr(x[x_1, x_3] - x[x_1, x_2] - x[x_2, x_3] >= -1)
					else:
						x1const, x11, x12 = get_x_var_consts(x, x_1, x_2)
						x2const, x21, x22 = get_x_var_consts(x, x_2, x_3)
						x3const, x31, x32 = get_x_var_consts(x, x_1, x_3)
						self.__add_matrix_constraint([x[x11, x12], x[x21, x22], x[x31, x32]], [-x1const, -x2const, x3const], (1 - x1const)//2 + (1 - x2const)//2 - (1 - x3const)//2)
						# model.addConstr(x1const * x[x11, x12] + x2const * x[x21, x22] - x3const * x[x31, x32] + (1 - x1const)//2 + (1 - x2const)//2 - (1 - x3const)//2 >= 0)
						self.__add_matrix_constraint([x[x11, x12], x[x21, x22], x[x31, x32]], [x1const, x2const, -x3const], -(1 - x1const)//2 - (1 - x2const)//2 + (1 - x3const)//2 + 1)
						# model.addConstr(-1 * x1const * x[x11, x12] - x2const * x[x21, x22] + x3const * x[x31, x32] - (1 - x1const)//2 - (1 - x2const)//2 + (1 - x3const)//2 >= -1)

		""" Original vertical position constraints as in Stratisfimal Layout """
		if self.vertical_transitivity:
			for x_var in x:
				self.__add_matrix_constraint([z[x_var], x[x_var]], [1, -self.m_val], 0)
				# model.addConstr(z[x_var] - self.m_val * x[x_var] <= 0, f"1.1z{x_var}")
				self.__add_matrix_constraint([z[x_var], y[x_var[0]], x[x_var]], [-1, 1, self.m_val], self.m_val)
				# model.addConstr(z[x_var] - y[x_var[0]] - (self.m_val * x[x_var]) >= -1 * self.m_val, f"1.2z{x_var}")
				self.__add_matrix_constraint([y[x_var[1]], z[x_var], x[x_var]], [-1, 1, 1], 0)
				# model.addConstr(y[x_var[1]] - z[x_var] - x[x_var] >= 0, f"1.3z{x_var}")
				self.__add_matrix_constraint([z[x_var], y[x_var[0]]], [1, -1], 0)
				# model.addConstr(z[x_var] <= y[x_var[0]], f"1.4z{x_var}")
				# m.addConstr(z[x_var] >= 0, f"1.5z{x_var}")	# already enforced by variable bounds

				self.__add_matrix_constraint([z[x_var[1], x_var[0]], x[x_var]], [1, self.m_val], self.m_val)
				# model.addConstr(z[x_var[1], x_var[0]] - self.m_val * (1 - x[x_var]) <= 0, f"2.1z{x_var}")
				self.__add_matrix_constraint([z[x_var[1], x_var[0]], y[x_var[1]], x[x_var]], [-1, 1, -self.m_val], 0)
				# model.addConstr(z[x_var[1], x_var[0]] - y[x_var[1]] - self.m_val * (1 - x[x_var]) >= -1 * self.m_val, f"2.2z{x_var}")
				self.__add_matrix_constraint([y[x_var[0]], z[x_var[1], x_var[0]], x[x_var]], [-1, 1, -1], -1)
				# model.addConstr(y[x_var[0]] - z[x_var[1], x_var[0]] - (1 - x[x_var]) >= 0, f"2.3z{x_var}")
				self.__add_matrix_constraint([z[x_var[1], x_var[0]], y[x_var[1]]], [1, -1], 0)
				# model.addConstr(z[x_var[1], x_var[0]] <= y[x_var[1]], f"2.4z{x_var}")
				# m.addConstr(z[x_var[1],x_var[0]] >= 0, f"2.5z{x_var}")	# already enforced by variable bounds

	def __symmetry_constraints_matrix(self, c, x):
		if self.mirror_vars and self.symmetry_constraints:
			x_v_seen = set()
			c_v_seen = set()
			for x_var in x:
				if (x_var[1], x_var[0]) not in x_v_seen:
					self.__add_matrix_constraint([x[x_var], x[x_var[1], x_var[0]]], [1, 1], 1, equality=True)
					# model.addConstr(x[x_var] + x[x_var[1], x_var[0]] == 1)
					x_v_seen.add(x_var)
			for c_var in c:
				if (c_var[1], c_var[0]) not in c_v_seen:
					self.__add_matrix_constraint([c[c_var], c[c_var[1], c_var[0]]], [1, -1], 0, equality=True)
					# model.addConstr(c[c_var] == c[c_var[1], c_var[0]])
					c_v_seen.add(c_var)

	def __symmetry_breaking_matrix(self, x, xvar_usage):
		if self.symmetry_breaking:
			if x:
				if xvar_usage != {}:
					most_used_x = max(xvar_usage, key=xvar_usage.get)
				else:
					most_used_x = random.choice(list(x.keys()))
				self.bounds[x[most_used_x]] = (0, 0)
				# model.getVarByName(f"x[{most_used_x[0]},{most_used_x[1]}]").lb = 0
				# model.getVarByName(f"x[{most_used_x[0]},{most_used_x[1]}]").ub = 0

	def __cycle_constraints_matrix(self, g: LayeredGraph, c):
		if self.cycle_constraints:
			# create vertex-exchange graph
			v_e_g, node_list, edge_signs = g.vertex_exchange_graph()

			# Find BFS tree, then calculate all fundamental cycles
			veg_bfsq = []
			veg_visit = [False] * len(node_list)
			n_visited = 0
			start_search_at = 0
			veg_parent = [-1] * len(node_list)
			fund_cycles = []
			while n_visited < len(node_list):
				while veg_visit[start_search_at]:
					start_search_at += 1
				veg_bfsq.append(start_search_at)
				veg_visit[start_search_at] = True
				n_visited += 1
				while veg_bfsq:
					next_layer = veg_bfsq.copy()
					veg_bfsq.clear()
					for nd in next_layer:
						for nd_adj in v_e_g[nd]:
							if not veg_visit[nd_adj[0]]:
								veg_visit[nd_adj[0]] = True
								veg_parent[nd_adj[0]] = nd
								veg_bfsq.append(nd_adj[0])
								n_visited += 1
							elif veg_parent[nd] != nd_adj[0]:
								branch1 = nd
								branch2 = nd_adj[0]
								fcycle_edges = [(nd, nd_adj[0], nd_adj[1])]
								while veg_parent[branch1] != branch2 and branch1 != veg_parent[branch2] and branch1 != branch2:
									if veg_parent[branch1] != -1:
										fcycle_edges.append((branch1, veg_parent[branch1], edge_signs[branch1, veg_parent[branch1]] if (branch1, veg_parent[branch1]) in edge_signs else edge_signs[veg_parent[branch1], branch1]))
										branch1 = veg_parent[branch1]
									if veg_parent[branch2] != -1:
										fcycle_edges.append((branch2, veg_parent[branch2], edge_signs[branch2, veg_parent[branch2]] if (branch2, veg_parent[branch2]) in edge_signs else edge_signs[veg_parent[branch2], branch2]))
										branch2 = veg_parent[branch2]
								if branch1 != branch2:
									fcycle_edges.append((branch1, branch2, edge_signs[branch1, branch2] if (branch1, branch2) in edge_signs else edge_signs[branch2, branch1]))
								fund_cycles.append(fcycle_edges)

			# Add cycle constraints from Healy and Kuusik, "The Vertex-Exchange Graph"
			cvars_set = set(c)
			for i, fcycle in enumerate(fund_cycles):
				if len(fcycle) <= 6:
					label = sum((fc[2] for fc in fcycle)) % 2
					c_idxs = []
					for edg in fcycle:
						u, v = node_list[edg[0]], node_list[edg[1]]
						if g.node_ids[u[0]].layer > g.node_ids[v[0]].layer:
							u, v = v, u
						if ((u[0], v[0]), (u[1], v[1])) in cvars_set or ((u[1], v[1]), (u[0], v[0])) in cvars_set:
							if (g.node_ids[u[0]].y < g.node_ids[u[1]].y) == (g.node_ids[v[0]].y < g.node_ids[v[1]].y):
								e1 = (u[0], v[0]) if edg[2] == 0 else (u[0], v[1])  # if sign doesn't match ypos comparison then we know this is a butterfly, and we got the wrong two edges in cvars
								e2 = (u[1], v[1]) if edg[2] == 0 else (u[1], v[0])
							else:
								e1 = (u[0], v[0]) if edg[2] == 1 else (u[0], v[1])
								e2 = (u[1], v[1]) if edg[2] == 1 else (u[1], v[0])
						elif ((u[0], v[1]), (u[1], v[0])) in cvars_set or ((u[1], v[0]), (u[0], v[1])) in cvars_set:
							if (g.node_ids[u[0]].y < g.node_ids[u[1]].y) == (g.node_ids[v[0]].y > g.node_ids[v[1]].y):  # XAND
								e1 = (u[0], v[1]) if edg[2] == 0 else (u[0], v[0])
								e2 = (u[1], v[0]) if edg[2] == 0 else (u[1], v[1])
							else:
								e1 = (u[0], v[1]) if edg[2] == 1 else (u[0], v[0])
								e2 = (u[1], v[0]) if edg[2] == 1 else (u[1], v[1])
						else:
							raise Exception("problem with vertex exchange graph")
						u, v = get_c_var(cvars_set, e1, e2)
						c_idxs.append(c[u, v])
					if label == 0:
						# model.addConstr(csum <= len(fcycle) - 1)
						self.__add_matrix_constraint(c_idxs, [1] * len(c_idxs), len(fcycle) - 1)
						# self.__add_late_breaking_variable(bound=(0, len(fcycle) / 2), integral=False)
						# self.__add_matrix_constraint(c_idxs + [self.nvars - 1], [1] * len(c_idxs) + [-2], 0, equality=True)
					else:
						# model.addConstr(csum >= 1)
						self.__add_matrix_constraint(c_idxs, [-1] * len(c_idxs), -1)
						# self.__add_late_breaking_variable(bound=(0, len(fcycle) / 2 - 1), integral=False)
						# self.__add_matrix_constraint(c_idxs + [self.nvars - 1], [1] * len(c_idxs) + [-2], 1, equality=True)
						# self.__add_matrix_constraint(c_idxs, [-1] * len(c_idxs), -1)

	def __optimize_subgraphs_matrix(self, graph, x, t1, num_crossings):
		if self.collapse_subgraphs or self.collapse_leaves:
			t4 = time.time()
			# new_g, stack_node_to_nodelist, node_to_stack_node, subgraphs_collapsed, subg_types
			subgraphs = graph.create_layered_graphs_from_subgraphs()
			for i, subgraph in enumerate(subgraphs):
				if graph.subgraph_types[i] != 0:
					var_to_fix = -1
					xvars_to_fix = {}
					if graph.subgraph_types[i] != 1:
						for nd in graph.contact_nodes[i]:
							if nd in subgraph:
								var_to_fix = nd
					if var_to_fix != -1:
						connect_nd = -1
						relative_xval = -1
						for nd, adj in graph.crossing_edges.items():
							if nd == var_to_fix:
								connect_nd = adj[0]
							elif var_to_fix in adj:
								connect_nd = nd
						if connect_nd not in graph:
							connect_nd = graph.node_to_stack_node[connect_nd]
						for nd_other in subgraph.nodes:
							if nd_other.layer == graph[connect_nd].layer:
								relative_xval = get_x_var(self.x_var_assign, graph.node_to_stack_node[nd_other.id], connect_nd)
						if relative_xval != -1:
							for other_nd in subgraph.layers[subgraph[var_to_fix].layer]:
								if other_nd.id != var_to_fix:
									xvars_to_fix[var_to_fix, other_nd.id] = 1 - relative_xval
					optim = HiGHSLayeredOptimizer(subgraph)
					# print(xvars_to_fix)
					if xvars_to_fix != {}:
						num_crossings += optim.optimize_layout(fix_x_vars=xvars_to_fix)[1]
					else:
						num_crossings += optim.optimize_layout()[1]
					for xv, assgn in optim.x_var_assign.items():
						self.x_var_assign[xv] = assgn
				else:
					for j, nd in enumerate(subgraph.nodes):
						for nd2 in subgraph.nodes[j+1:]:
							if nd.layer == nd2.layer:
								if (nd.id, nd2.id) in self.x_var_assign:
									self.x_var_assign[nd.id, nd2.id] = 0
								else:
									self.x_var_assign[nd2.id, nd.id] = 1
			to_remove = []
			for x_var in x:
				if graph[x_var[0]].stacked and graph[x_var[1]].stacked:
					for otv in graph.stack_node_to_nodelist[x_var[0]]:
						for otv2 in graph.stack_node_to_nodelist[x_var[1]]:
							set_x_var(self.x_var_assign, otv, otv2, get_x_var(self.x_var_assign, x_var[0], x_var[1]))
					to_remove.append(x_var)
				elif graph[x_var[0]].stacked:
					for otv in graph.stack_node_to_nodelist[x_var[0]]:
						set_x_var(self.x_var_assign, otv, x_var[1], get_x_var(self.x_var_assign, x_var[0], x_var[1]))
					to_remove.append(x_var)
				elif graph[x_var[1]].stacked:
					for otv in graph.stack_node_to_nodelist[x_var[1]]:
						set_x_var(self.x_var_assign, x_var[0], otv, get_x_var(self.x_var_assign, x_var[0], x_var[1]))
					to_remove.append(x_var)
			for xv in to_remove:
				if xv in self.x_var_assign:
					del self.x_var_assign[xv]
				else:
					del self.x_var_assign[xv[1], xv[0]]
			return graph.old_g, t1 + time.time() - t4, num_crossings
		else:
			return graph, t1, num_crossings

	def __add_3claw_constraints_matrix(self, g, c_vars):
		if self.claw_constraints:
			bearclaws = motifs.get_3claws(g)
			print(f"3-claws found: {len(bearclaws)}")
			for claw in bearclaws:
				claw_cvs = [
							# c_vars[get_c_var(c_vars, claw[3], claw[4])],
							c_vars[get_c_var(c_vars, claw[3], claw[1])],
							# c_vars[get_c_var(c_vars, claw[3], claw[5])],
							c_vars[get_c_var(c_vars, claw[3], claw[2])],
							c_vars[get_c_var(c_vars, claw[0], claw[4])],
							c_vars[get_c_var(c_vars, claw[0], claw[5])],
							# c_vars[get_c_var(c_vars, claw[4], claw[5])],
							c_vars[get_c_var(c_vars, claw[4], claw[2])],
							c_vars[get_c_var(c_vars, claw[1], claw[5])]]
				self.__add_matrix_constraint(claw_cvs, [-1] * 9, -1)

	def __add_domepath_constraints_matrix(self, g, x_vars, c_vars):
		if self.dome_path_constraints:
			domes = motifs.get_domepaths(g)
			print(f"Dome-paths found: {len(domes)}")
			for dome in domes:
				if dome[0][0] == dome[1][0]:
					klc, kl1, kl2 = get_x_var_consts(x_vars, dome[2][0], dome[0][0])
					kmc, km1, km2 = get_x_var_consts(x_vars, dome[2][0], dome[3][0])
					lmc, lm1, lm2 = get_x_var_consts(x_vars, dome[0][0], dome[3][0])
					cikjl1, cikjl2 = get_c_var(c_vars, dome[2], dome[1])
					ciljm1, ciljm2 = get_c_var(c_vars, dome[0], dome[3])
					# m.addConstr(klc * x[kl1, kl2] - 2 * kmc * x[km1, km2] + lmc * x[lm1, lm2] - c[cikjl1, cikjl2] - c[ciljm1, ciljm2] + (1 - klc)//2 - (1 - kmc) + (1 - lmc)//2 <= 0)
					self.__add_matrix_constraint([x_vars[kl1, kl2], x_vars[km1, km2], x_vars[lm1, lm2], c_vars[cikjl1, cikjl2], c_vars[ciljm1, ciljm2]], [klc, -2*kmc, lmc, -1, -1], 0 - (1 - klc)//2 + (1 - kmc) - (1 - lmc)//2)
					# m.addConstr(-klc * x[kl1, kl2] + 2 * kmc * x[km1, km2] - lmc * x[lm1, lm2] - c[cikjl1, cikjl2] - c[ciljm1, ciljm2] - (1 - klc)//2 + (1 - kmc) - (1 - lmc)//2 <= 0)
					self.__add_matrix_constraint([x_vars[kl1, kl2], x_vars[km1, km2], x_vars[lm1, lm2], c_vars[cikjl1, cikjl2], c_vars[ciljm1, ciljm2]], [-klc, 2*kmc, -lmc, -1, -1], (1 - klc)//2 - (1 - kmc) + (1 - lmc)//2)

	def __sequential_bendiness_matrix(self):
		if self.bendiness_reduction and self.sequential_bendiness:
			g = self.g
			x_var_opt = self.x_var_assign
			ct = 0
			y_vars_br = {}
			b_vars_br = {}
			for nd in g.nodes:
				y_vars_br[nd.id] = ct
				ct += 1
			for ed in g.edge_ids:
				b_vars_br[ed] = ct
				ct += 1
			a_ub = []
			b_ub = []
			c_t = [0] * ct
			for idx in b_vars_br.values():
				c_t[idx] = 1
			for var, val in x_var_opt.items():
				cnstr = [0] * ct
				cnstr[y_vars_br[var[0]]] = -1 if val == 0 else 1
				cnstr[y_vars_br[var[1]]] = 1 if val == 0 else -1
				a_ub.append(cnstr)
				if g[int(var[0])].is_anchor_node and g[int(var[1])].is_anchor_node:
					b_ub.append(-0.15)
					# m2.addConstr(y[var[0]] >= 0.15 + y[var[1]], f"vert{var}")
				elif g[int(var[0])].is_anchor_node or g[int(var[1])].is_anchor_node:
					b_ub.append(-0.3)
					# m2.addConstr(y[var[0]] >= 0.3 + y[var[1]], f"vert{var}")
				else:
					b_ub.append(-1)
					# m2.addConstr(y[var[0]] >= 1 + y[var[1]], f"vert{var}")
				# if g[int(var[0])].is_anchor_node and g[int(var[1])].is_anchor_node:
				# 	m2.addConstr(y[var[0]] + 0.15 <= y[var[1]], f"vert{var}")
				# elif g[int(var[0])].is_anchor_node or g[int(var[1])].is_anchor_node:
				# 	m2.addConstr(y[var[0]] + 0.3 <= y[var[1]], f"vert{var}")
				# else:
				# 	m2.addConstr(y[var[0]] + 1 <= y[var[1]], f"vert{var}")
			for b_v, b_idx in b_vars_br.items():
				cnstr1 = [0] * ct
				cnstr1[y_vars_br[b_v[0]]] = 1
				cnstr1[y_vars_br[b_v[1]]] = -1
				cnstr1[b_idx] = -1
				a_ub.append(cnstr1)
				b_ub.append(0)
				# m2.addConstr(y[b_var[0]] - y[b_var[1]] <= b[b_var], f"1bend{b_var}")
				cnstr2 = [0] * ct
				cnstr2[y_vars_br[b_v[0]]] = -1
				cnstr2[y_vars_br[b_v[1]]] = 1
				cnstr2[b_idx] = -1
				a_ub.append(cnstr2)
				b_ub.append(0)
				# m2.addConstr(y[b_var[1]] - y[b_var[0]] <= b[b_var], f"2bend{b_var}")

			seq_res = linprog(c_t, method="highs", A_ub=a_ub, b_ub=b_ub)
			for nd in g.nodes:
				nd.y = seq_res.x[y_vars_br[nd.id]]
