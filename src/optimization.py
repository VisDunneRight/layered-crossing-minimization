import os.path
import random
import time
import itertools
import gurobipy as gp
from gurobipy import GRB
from sklearn.cluster import SpectralClustering
from src import vis, reductions, motifs, type_conversions, read_data
from src.graph import LayeredGraph, CollapsedGraph
from src.helpers import *


class LayeredOptimizer:
	def __init__(self, layered_graph, parameters=None):
		if parameters is None:
			parameters = {}
		assert (type(layered_graph) == str and os.path.isfile(layered_graph)) or type(layered_graph) == LayeredGraph or type(layered_graph) == CollapsedGraph, "input needs to be a path to graph file or a LayeredGraph object"
		if type(layered_graph) == LayeredGraph or type(layered_graph) == CollapsedGraph:
			self.g = layered_graph
		else:
			self.g = read_data.read(layered_graph)
		self.x_var_assign = {x_v: 2 for n_l in self.g.get_ids_by_layer().values() for x_v in itertools.combinations(n_l, 2)}
		self.bendiness_reduction = parameters["bendiness_reduction"] if "bendiness_reduction" in parameters else False
		self.gamma_1 = parameters["gamma_1"] if "gamma_1" in parameters else 1
		self.gamma_2 = parameters["gamma_2"] if "gamma_2" in parameters else 1
		self.m_val = parameters["m_val"] if "m_val" in parameters else max(len(layer) for layer in self.g.layers.values())
		self.sequential_bendiness = parameters["sequential_bendiness"] if "sequential_bendiness" in parameters else True
		self.return_full_data = parameters["return_full_data"] if "return_full_data" in parameters else False
		self.cutoff_time = parameters["cutoff_time"] if "cutoff_time" in parameters else 0
		self.do_subg_reduction = parameters["do_subg_reduction"] if "do_subg_reduction" in parameters else False
		self.return_x_vars = parameters["return_x_vars"] if "return_x_vars" in parameters else False
		self.butterfly_reduction = parameters["butterfly_reduction"] if "butterfly_reduction" in parameters else False
		self.verbose = parameters["verbose"] if "verbose" in parameters else False
		self.subg_verbose = parameters["subg_verbose"] if "subg_verbose" in parameters else True
		self.draw_graph = parameters["draw_graph"] if "draw_graph" in parameters else False
		self.symmetry_breaking = parameters["symmetry_breaking"] if "symmetry_breaking" in parameters else False
		self.heuristic_start = parameters["heuristic_start"] if "heuristic_start" in parameters else False
		self.aggro_presolve = parameters["presolve"] if "presolve" in parameters else False
		self.mip_relax = parameters["mip_relax"] if "mip_relax" in parameters else False
		self.xvar_branch_priority = parameters["xvar_branch_priority"] if "xvar_branch_priority" in parameters else False
		self.direct_transitivity = parameters["direct_transitivity"] if "direct_transitivity" in parameters else False
		self.vertical_transitivity = parameters["vertical_transitivity"] if "vertical_transitivity" in parameters else False
		self.mirror_vars = parameters["mirror_vars"] if "mirror_vars" in parameters else False
		self.stratisfimal_y_vars = parameters["stratisfimal_yvars"] if "stratisfimal_yvars" in parameters else False
		self.symmetry_constraints = parameters["symmetry_constraints"] if "symmetry_constraints" in parameters else True
		self.cycle_constraints = parameters["cycle_constraints"] if "cycle_constraints" in parameters else False
		self.collapse_subgraphs = parameters["collapse_subgraphs"] if "collapse_subgraphs" in parameters else False
		self.collapse_leaves = parameters["collapse_leaves"] if "collapse_leaves" in parameters else False
		self.claw_constraints = parameters["claw_constraints"] if "claw_constraints" in parameters else False
		self.dome_path_constraints = parameters["dome_path_constraints"] if "dome_path_constraints" in parameters else False
		self.polyhedral_constraints = parameters["polyhedral_constraints"] if "polyhedral_constraints" in parameters else False
		self.return_experiment_data = parameters["return_experiment_data"] if "return_experiment_data" in parameters else False
		self.name = parameters["name"] if "name" in parameters else "graph1"
		self.print_info = []
		if self.polyhedral_constraints:
			self.claw_constraints, self.dome_path_constraints = True, True

	def __optimize_layout_standard(self, graph_arg=None, assignment=None, name="graph1", fix_x_vars=None, start_x_vars=None, is_subgraph=False, use_top_level_params=False):
		g = self.g if graph_arg is None else graph_arg
		if self.polyhedral_constraints:
			self.claw_constraints, self.dome_path_constraints = True, True
		if not self.direct_transitivity and not self.vertical_transitivity:
			self.direct_transitivity = True
		t1 = time.time()

		""" Collapse valid subgraphs """
		if self.collapse_leaves:
			g = g.collapse_leaves()
			self.x_var_assign = {x_v: 2 for n_l in g.get_ids_by_layer().values() for x_v in itertools.combinations(n_l, 2)}
			# vis.draw_graph(g, "after collapse", groups={nd.name: 1 if nd.stacked else 0 for nd in g.nodes})

		""" Create model and graph data """
		nodes_by_layer = g.get_ids_by_layer()
		edges_by_layer = g.get_edge_ids_by_layer()
		n_constraints_generated = [0] * 6  # simple edge, hybrid edge, same layer edge, vertical pos, bendiness, total
		pre_sym = '\t' if is_subgraph else ''
		if self.verbose:
			self.print_info.extend(('-' * 70, f"{self.name if use_top_level_params else name}:", '-' * 70))
		m = gp.Model()

		""" Add all variables """
		x_vars = []
		z_vars = []
		relax_type = GRB.INTEGER if not self.mip_relax else GRB.CONTINUOUS
		for i, name_list in nodes_by_layer.items():
			if self.mirror_vars:
				x_vars += list(itertools.permutations(name_list, 2))
			else:
				x_vars += list(itertools.combinations(name_list, 2))
			if self.stratisfimal_y_vars:
				z_vars += list(itertools.permutations(name_list, 2))
		x = m.addVars(x_vars, vtype=GRB.BINARY, name="x")
		if self.stratisfimal_y_vars:
			z = m.addVars(z_vars, vtype=relax_type, lb=0, ub=self.m_val, name="z")
		else:
			z = None
		c_vars, c_consts = reductions.normal_c_vars(g, edges_by_layer, self.mirror_vars)
		if self.mirror_vars:
			c_vars_orig, nc_consts = reductions.normal_c_vars(g, edges_by_layer, False)
		c = m.addVars(c_vars, vtype=relax_type, name="c")
		if self.vertical_transitivity or self.stratisfimal_y_vars:
			y_vars = [n.id for n in g]
			y = m.addVars(y_vars, vtype=relax_type, lb=0, ub=self.m_val, name="y")
		else:
			y = None
		m.update()

		""" Fix variables/set starting assignments, used in some old experiments """
		if assignment:  # TODO (later): can be removed, not used in full subg algo (remove in place of fix_x_vars)
			for k, v in assignment.items():
				m.getVarByName(k).ub, m.getVarByName(k).lb = v, v
		if fix_x_vars:
			for k, v in fix_x_vars.items():
				if k in x:
					m.getVarByName(f"x[{k[0]},{k[1]}]").lb, m.getVarByName(f"x[{k[0]},{k[1]}]").ub = v, v
				else:
					m.getVarByName(f"x[{k[1]},{k[0]}]").lb, m.getVarByName(f"x[{k[1]},{k[0]}]").ub = 1 - v, 1 - v
		if start_x_vars:
			for v in m.getVars():
				v.Start = start_x_vars[v.varName]

		""" Heuristic starting assignments """
		self.__heuristic_start(m, g)

		""" Set higher priority for branching on x-vars """
		if self.xvar_branch_priority:
			for v in m.getVars():
				if v.varName[:1] == "x":
					v.BranchPriority = 1
				else:
					v.BranchPriority = 0

		""" Butterfly reduction """
		butterfly_c_pairs = self.get_butterfly_cvars(g, c_vars)

		""" Set model objective function """
		if not self.sequential_bendiness:
			if self.bendiness_reduction:
				b_vars = list(g.edge_ids.keys())
				b = m.addVars(b_vars, vtype=GRB.INTEGER, lb=0, ub=self.m_val, name="b")
				m.setObjective(self.gamma_1 * c.sum() + self.gamma_2 * b.sum(), GRB.MINIMIZE)
			else:
				opt = gp.LinExpr()
				for i, c_var in enumerate(c_vars):
					opt += c_consts[i] * c[c_var]
				m.setObjective(opt, GRB.MINIMIZE)
		else:
			opt = gp.LinExpr()
			if self.mirror_vars and self.symmetry_constraints:
				for i, c_var in enumerate(c_vars_orig):
					opt += nc_consts[i] * c[c_var]
			else:
				for i, c_var in enumerate(c_vars):
					opt += c_consts[i] * c[c_var]
			m.setObjective(opt, GRB.MINIMIZE)

		""" Transitivity constraints """
		self.__transitivity(m, nodes_by_layer, x_vars, x, y, z)

		""" Edge crossing constraints """
		n_cs = self.__edge_crossings(m, c_vars, x, c, graph_arg=g, track_x_var_usage=self.symmetry_breaking, butterflies=butterfly_c_pairs)
		for i, val in enumerate(n_cs[:-1]):
			n_constraints_generated[i] += val

		""" 3-claw constraints """
		self.__add_3claw_constraints(m, g, c_vars, c)

		""" Dome path constraints """
		self.__add_dome_path_constraints(m, g, c_vars, c, x_vars, x)

		""" Symmetry constraints """
		self.__add_symmetry_constraints(m, x_vars, c_vars, x, c)

		""" Cycle constraints """
		self.__add_cycle_constraints(m, g, c_vars, c)

		""" Break symmetry by fixing key x-var """
		self.__symmetry_breaking(m, n_cs[-1], x_vars)

		""" Non-sequential bendiness reduction, original Stratisfimal version"""
		if not self.sequential_bendiness and self.bendiness_reduction:
			for b_var in b_vars:
				m.addConstr(y[b_var[0]] - y[b_var[1]] <= b[b_var], f"1bend{b_var}")
				m.addConstr(y[b_var[1]] - y[b_var[0]] <= b[b_var], f"2bend{b_var}")
				n_constraints_generated[4] += 2

		""" Optimize model """
		if self.cutoff_time > 0:
			m.setParam("TimeLimit", self.cutoff_time)
		m.setParam("OutputFlag", 0)
		if self.aggro_presolve:
			m.setParam("Presolve", 2)
		n_constraints_generated[5] = sum(n_constraints_generated[:5])
		t1 = time.time() - t1
		if self.verbose:
			self.print_info.append(f"{pre_sym}Time to input constraints: {t1}")
			self.print_info.append(f"{pre_sym}Constraint counts: {n_constraints_generated}")
		t2 = time.time()
		m.optimize()
		t2 = time.time() - t2
		if self.verbose:
			self.print_info.append(f"{pre_sym}Objective: {m.objVal}")
			self.print_info.append(f"{pre_sym}Time to optimize: {t2}")
		if (m.status != 2 and m.status != 9) or m.SolCount == 0:
			print("model returned status code:", m.status)
			print("4: model probably never found a feasible solution")
			print("11: solve interrupted")
			print("otherwise check https://www.gurobi.com/documentation/current/refman/optimization_status_codes.html")
			print("or, was cutoff without finding a solution. #solutions found:", m.SolCount)
			return 0, 0, 0, 0, 0, float('inf'), m.status, 0, 0, "INCORRECT STATUS"
		for v in m.getVars():
			if v.varName[:1] == "x":
				self.x_var_assign[int(v.varName[2:v.varName.index(',')]), int(v.varName[v.varName.index(',') + 1:v.varName.index(']')])] = round(v.x)
		num_crossings = round(m.objVal)

		""" Optimize and merge collapsed subgraphs """
		g, t1, num_crossings = self.__optimize_subgraphs(g, x_vars, t1, num_crossings)

		""" Draw pre-bendiness graph """
		# vis.draw_graph(g, "interim")

		""" Sequential bendiness reduction """
		if self.bendiness_reduction:
			t3 = time.time()
			n_constraints_generated[4] += self.__sequential_br(graph_arg=g, subgraph_seq=is_subgraph)
			t3 = time.time() - t3
			if not is_subgraph and self.verbose:
				self.print_info.append(f"{pre_sym}Time to perform bendiness reduction: {t3}")
		else:
			t3 = 0

		""" Draw the resulting graph """
		if self.draw_graph:
			if not self.bendiness_reduction:
				if self.direct_transitivity:
					g.assign_y_vals_given_x_vars(self.x_var_assign)
				else:
					for v in m.getVars():
						if v.varName[:1] == "y" and int(v.varName[2:v.varName.index(']')]) in g.node_ids:
							g[int(v.varName[2:v.varName.index(']')])].y = float(v.x)
			vis.draw_graph(g, self.name)

		if self.verbose:
			self.print_info.append(f"{pre_sym}Final edge crossing count: {g.num_edge_crossings()}")
			self.print_info.append(f"{pre_sym}Number of constraints: {n_constraints_generated}")
			self.print_info.append(f"{pre_sym}{round(t1, 3)}, {round(t2, 3)}, {round(t3, 3)}, {round(t1 + t2 + t3, 3)}")

		""" Return data """
		print(f"Number of crossings: {num_crossings}", f"\tOptimization time: {round(m.runtime, 3)}")
		# print(f"g calc: {g.num_edge_crossings()}")

		if self.return_x_vars:
			return num_crossings, self.x_var_assign

		if self.return_experiment_data:
			return len(x_vars), len(c_vars), m.numVars, m.numConstrs, num_crossings, m.runtime, m.status, int(m.nodeCount), round(t1, 3)

		return round(t1 + t2 + t3, 3), num_crossings

	def __sequential_br(self, graph_arg=None, substitute_x_vars=None, subgraph_seq=False):
		g = self.g if graph_arg is None else graph_arg
		x_var_opt = self.x_var_assign if substitute_x_vars is None else substitute_x_vars
		y_vars = list(g.node_ids)
		n_constr = 0
		m2 = gp.Model()
		y = m2.addVars(y_vars, vtype=GRB.CONTINUOUS, lb=0, ub=self.m_val, name="y")
		# y = m2.addVars(y_vars, vtype=GRB.INTEGER, lb=0, ub=mv, name="y")
		m2.update()
		for v in m2.getVars():
			v.start = g[int(v.varName[2:v.varName.index(']')])].y
		b_vars = list(g.edge_ids.keys())
		b = m2.addVars(b_vars, vtype=GRB.CONTINUOUS, lb=0, ub=self.m_val, name="b")
		# b = m2.addVars(b_vars, vtype=GRB.INTEGER, lb=0, ub=mv, name="b")
		m2.setObjective(b.sum(), GRB.MINIMIZE)
		n_orig = sum((1 for node in g.nodes if not node.is_anchor_node))
		if subgraph_seq:
			n_orig = max((n.id for n in g.nodes))
		for var, val in x_var_opt.items():
			if val == 0:
				if g[int(var[0])].is_anchor_node and int(var[1]) > n_orig:
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

	def __transitivity(self, model: gp.Model, names_by_layer, x_vars, x, y, z):
		if self.direct_transitivity:
			for x_vars_list in names_by_layer.values():
				for x_1, x_2, x_3 in itertools.combinations(x_vars_list, 3):
					if self.mirror_vars:
						model.addConstr(x[x_1, x_2] + x[x_2, x_3] - x[x_1, x_3] >= 0)
						model.addConstr(x[x_1, x_3] - x[x_1, x_2] - x[x_2, x_3] >= -1)
					else:
						x1const, x11, x12 = get_x_var_consts(x_vars, x_1, x_2)
						x2const, x21, x22 = get_x_var_consts(x_vars, x_2, x_3)
						x3const, x31, x32 = get_x_var_consts(x_vars, x_1, x_3)
						model.addConstr(x1const * x[x11, x12] + x2const * x[x21, x22] - x3const * x[x31, x32] + (1 - x1const)//2 + (1 - x2const)//2 - (1 - x3const)//2 >= 0)
						model.addConstr(-1 * x1const * x[x11, x12] - x2const * x[x21, x22] + x3const * x[x31, x32] - (1 - x1const)//2 - (1 - x2const)//2 + (1 - x3const)//2 >= -1)

		""" Vertical position, implication version """
		# Uses big-M method: https://support.gurobi.com/hc/en-us/articles/4414392016529-How-do-I-model-conditional-statements-in-Gurobi-
		if self.vertical_transitivity:
			for x_var in x_vars:
				model.addGenConstrIndicator(x[x_var], True, y[x_var[0]] + 1 <= y[x_var[1]])
				model.addGenConstrIndicator(x[x_var], False, y[x_var[0]] >= 1 + y[x_var[1]])

		""" Original vertical position constraints as in Stratisfimal Layout """
		if self.stratisfimal_y_vars:
			for x_var in x_vars:
				model.addConstr(z[x_var] - self.m_val * x[x_var] <= 0, f"1.1z{x_var}")
				model.addConstr(z[x_var] - y[x_var[0]] - (self.m_val * x[x_var]) >= -1 * self.m_val, f"1.2z{x_var}")
				model.addConstr(y[x_var[1]] - z[x_var] - x[x_var] >= 0, f"1.3z{x_var}")
				model.addConstr(z[x_var] <= y[x_var[0]], f"1.4z{x_var}")
				# m.addConstr(z[x_var] >= 0, f"1.5z{x_var}")	# print(f"z{x_var}-{self.m_val}*x{x_var} <= 0"), print(f"y{x_var[1]} - z{x_var} - x{x_var} >= 0")

				model.addConstr(z[x_var[1], x_var[0]] - self.m_val * (1 - x[x_var]) <= 0, f"2.1z{x_var}")
				model.addConstr(z[x_var[1], x_var[0]] - y[x_var[1]] - self.m_val * (1 - x[x_var]) >= -1 * self.m_val, f"2.2z{x_var}")
				model.addConstr(y[x_var[0]] - z[x_var[1], x_var[0]] - (1 - x[x_var]) >= 0, f"2.3z{x_var}")
				model.addConstr(z[x_var[1], x_var[0]] <= y[x_var[1]], f"2.4z{x_var}")
				# m.addConstr(z[x_var[1],x_var[0]] >= 0, f"2.5z{x_var}")	# print(f"z[{x_var[1],x_var[0]}]-{self.m_val}*(1-x{x_var}) <= 0"), print(f"y{x_var[0]} - z[{x_var[1],x_var[0]}] - (1-x{x_var}) >= 0")

				# m.addConstr(k[x_var] + x[x_var] - 1 >= 0)  # SOS-constraint version. Use with full LP relaxation
				# m.addSOS(GRB.SOS_TYPE1, [k[x_var], x[x_var]])

	def __edge_crossings(self, model: gp.Model, c_vars, x, c, graph_arg=None, track_x_var_usage=False, butterflies=None):
		g = self.g if graph_arg is None else graph_arg
		n_constr_0, n_constr_1, n_constr_2 = 0, 0, 0
		x_var_usage = {}
		if butterflies is not None:
			for c_var_pair in butterflies:
				model.addConstr(c[c_var_pair[0]] + c[c_var_pair[1]] == 1)
		for c_var in c_vars:
			if self.mirror_vars and not g.edge_ids[c_var[0]].same_layer_edge and not g.edge_ids[c_var[1]].same_layer_edge:
				model.addConstr(x[c_var[1][0], c_var[0][0]] + x[c_var[0][1], c_var[1][1]] + c[c_var] >= 1, f"1se{c_var}")
				model.addConstr(x[c_var[0][0], c_var[1][0]] + x[c_var[1][1], c_var[0][1]] + c[c_var] >= 1, f"2se{c_var}")
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

				# simple edge crossing
				if not g.edge_ids[c_var[0]].same_layer_edge and not g.edge_ids[c_var[1]].same_layer_edge:
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
				elif g.edge_ids[c_var[0]].same_layer_edge and not g.edge_ids[c_var[1]].same_layer_edge:
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

				elif g.edge_ids[c_var[1]].same_layer_edge and not g.edge_ids[c_var[0]].same_layer_edge:
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

	# def __edge_crossings_junger(self, model: gp.Model, c_vars, x_vars, x, c, graph_arg=None, track_x_var_usage=False, butterflies=None):
	# 	g = self.g if graph_arg is None else graph_arg
	# 	n_constraints = 0
	# 	x_var_usage = {}
	# 	for c_var in c_vars:
	# 		if (butterflies is None or c_var not in butterflies) and c_var[0][0] != c_var[1][0]:
	# 			if not self.mirror_vars:
	# 				if c_var[0][1] < c_var[1][1]:
	# 					x1const, x11, x12 = get_x_var_consts(x_vars, c_var[0][1], c_var[1][1])
	# 					x2const, x21, x22 = get_x_var_consts(x_vars, c_var[0][0], c_var[1][0])
	# 					model.addConstr(c[c_var] + x1const * x[x11, x12] - x2const * x[x21, x22] + (1 - x1const)//2 - (1 - x2const)//2 >= 0)
	# 					model.addConstr(c[c_var] - x1const * x[x11, x12] + x2const * x[x21, x22] - (1 - x1const)//2 + (1 - x2const)//2 >= 0)
	# 					if track_x_var_usage:
	# 						if (x11, x12) not in x_var_usage:
	# 							x_var_usage[x11, x12] = 0
	# 						if (x21, x22) not in x_var_usage:
	# 							x_var_usage[x21, x22] = 0
	# 						x_var_usage[x11, x12] += 1
	# 						x_var_usage[x21, x22] += 1
	# 				else:
	# 					x1const, x11, x12 = get_x_var_consts(x_vars, c_var[1][1], c_var[0][1])
	# 					x2const, x21, x22 = get_x_var_consts(x_vars, c_var[0][0], c_var[1][0])
	# 					model.addConstr(c[c_var] + x1const * x[x11, x12] + x2const * x[x21, x22] - 1 + (1 - x1const)//2 + (1 - x2const)//2 >= 0)
	# 					model.addConstr(c[c_var] - x1const * x[x11, x12] - x2const * x[x21, x22] + 1 - (1 - x1const)//2 - (1 - x2const)//2 >= 0)
	# 					if track_x_var_usage:
	# 						if (x11, x12) not in x_var_usage:
	# 							x_var_usage[x11, x12] = 0
	# 						if (x21, x22) not in x_var_usage:
	# 							x_var_usage[x21, x22] = 0
	# 						x_var_usage[x11, x12] += 1
	# 						x_var_usage[x21, x22] += 1
	# 			else:
	# 				if c_var[0][1] < c_var[1][1]:
	# 					model.addConstr(c[c_var] + x[c_var[0][1], c_var[1][1]] - x[c_var[0][0], c_var[1][0]] >= 0)
	# 					model.addConstr(c[c_var] - x[c_var[0][1], c_var[1][1]] + x[c_var[0][0], c_var[1][0]] >= 0)
	# 					if track_x_var_usage:
	# 						if (c_var[0][1], c_var[1][1]) not in x_var_usage:
	# 							x_var_usage[c_var[0][1], c_var[1][1]] = 0
	# 						if (c_var[0][0], c_var[1][0]) not in x_var_usage:
	# 							x_var_usage[c_var[0][0], c_var[1][0]] = 0
	# 						x_var_usage[c_var[0][1], c_var[1][1]] += 2
	# 						x_var_usage[c_var[0][0], c_var[1][0]] += 2
	# 				else:
	# 					model.addConstr(c[c_var] + x[c_var[1][1], c_var[0][1]] + x[c_var[0][0], c_var[1][0]] - 1 >= 0)
	# 					model.addConstr(c[c_var] - x[c_var[1][1], c_var[0][1]] - x[c_var[0][0], c_var[1][0]] + 1 >= 0)
	# 					if track_x_var_usage:
	# 						if (c_var[0][0], c_var[1][0]) not in x_var_usage:
	# 							x_var_usage[c_var[0][0], c_var[1][0]] = 0
	# 						if (c_var[1][1], c_var[0][1]) not in x_var_usage:
	# 							x_var_usage[c_var[1][1], c_var[0][1]] = 0
	# 						x_var_usage[c_var[0][0], c_var[1][0]] += 2
	# 						x_var_usage[c_var[1][1], c_var[0][1]] += 2
	# 			n_constraints += 2
	# 	return n_constraints, track_x_var_usage

	def __symmetry_breaking(self, model: gp.Model, x_var_usage, x_vars):
		if self.symmetry_breaking:
			if x_vars:
				if x_var_usage != {}:
					most_used_x = max(x_var_usage, key=x_var_usage.get)
				else:
					most_used_x = random.choice(x_vars)
				model.getVarByName(f"x[{most_used_x[0]},{most_used_x[1]}]").lb = 0
				model.getVarByName(f"x[{most_used_x[0]},{most_used_x[1]}]").ub = 0

	def get_butterfly_cvars(self, graph: LayeredGraph, c_vars):
		butterfly_c_vars = set()
		butterfly_c_pairs = []
		if self.butterfly_reduction:
			b_set_list = []
			for b_v in motifs.get_butterflies(graph):
				b_set_list.append(set(b_v))
			b_set_one_found = [0] * len(b_set_list)
			print("Butterfly set:", b_set_list)
			if b_set_list:
				for c_var in c_vars:
					c_set = {c_var[0][0], c_var[0][1], c_var[1][0], c_var[1][1]}
					if c_set in b_set_list:
						b_ind = b_set_list.index(c_set)
						if self.mirror_vars and b_set_one_found[b_ind] == 3:
							c_org = [c_v for c_v in butterfly_c_vars if {c_v[0][0], c_v[0][1], c_v[1][0], c_v[1][1]} == c_set]
							cv_p = [c_v for c_v in c_org if c_v[0][0] == c_var[0][0]][0]
							cv_r = [c_v for c_v in c_org if c_v[0][1] == c_var[0][1]][0]
							cv_q = [c_v for c_v in c_org if c_v[0] == c_var[1]][0]
							butterfly_c_pairs.append((c_var, cv_p))
							butterfly_c_pairs.append((c_var, cv_r))
							butterfly_c_pairs.append((cv_q, cv_p))
							butterfly_c_pairs.append((cv_q, cv_r))
						elif not self.mirror_vars and b_set_one_found[b_ind]:
							c_org = [c_v for c_v in butterfly_c_vars if {c_v[0][0], c_v[0][1], c_v[1][0], c_v[1][1]} == c_set][0]
							butterfly_c_pairs.append((c_org, c_var))
						else:
							b_set_one_found[b_ind] += 1
						butterfly_c_vars.add(c_var)
		return butterfly_c_pairs

	def __heuristic_start(self, model: gp.Model, graph: LayeredGraph):
		if self.heuristic_start:
			g_igraph = type_conversions.layered_graph_to_igraph(graph)
			heuristic_layout = g_igraph.layout_sugiyama(layers=g_igraph.vs["layer"])
			for i, coord in enumerate(heuristic_layout.coords[:graph.n_nodes]):
				graph.nodes[i].y = coord[0]
			for lay in graph.layers:
				graph.layers[lay].sort(key=lambda nd1: nd1.y)
				l_offset = 0
				for j in range(1, len(graph.layers[lay])):
					if graph.layers[lay][j].y <= graph.layers[lay][j - 1].y:
						l_offset += 1
					graph.layers[lay][j].y += l_offset
			for v in model.getVars():
				if v.varName[:1] == "y":
					v.Start = graph[int(v.varName[2:v.varName.index(']')])].y
				elif v.varName[:1] == "x":
					v.Start = 1 if graph[int(v.varName[2:v.varName.index(',')])].y < graph[int(v.varName[v.varName.index(',') + 1:v.varName.index(']')])].y else 0

	def __add_symmetry_constraints(self, model: gp.Model, x_vars, c_vars, x, c):
		if self.mirror_vars and self.symmetry_constraints:
			x_v_seen = set()
			c_v_seen = set()
			for x_var in x_vars:
				if (x_var[1], x_var[0]) not in x_v_seen:
					model.addConstr(x[x_var] + x[x_var[1], x_var[0]] == 1)
					x_v_seen.add(x_var)
			for c_var in c_vars:
				if (c_var[1], c_var[0]) not in c_v_seen:
					model.addConstr(c[c_var] == c[c_var[1], c_var[0]])
					c_v_seen.add(c_var)

	def __add_cycle_constraints(self, model: gp.Model, g: LayeredGraph, cvars, c):
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
			cvars_set = set(cvars)
			for i, fcycle in enumerate(fund_cycles):
				if len(fcycle) <= 6:
					label = sum((fc[2] for fc in fcycle)) % 2
					csum = gp.LinExpr()
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
						csum += c[u, v]
					if label == 0:
						# kc = model.addVar(0, len(fcycle) / 2, vtype=GRB.CONTINUOUS)
						# model.addConstr(2 * kc == csum)
						model.addConstr(csum <= len(fcycle) - 1)
					else:
						# kc = model.addVar(0, len(fcycle) / 2 - 1, vtype=GRB.CONTINUOUS)
						# model.addConstr(2 * kc + 1 == csum)
						model.addConstr(csum >= 1)

	def __optimize_subgraphs(self, graph, x_vars, t1, num_crossings):
		if self.collapse_leaves or self.collapse_subgraphs:
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
					optim = LayeredOptimizer(subgraph)
					# print(xvars_to_fix)
					if xvars_to_fix != {}:
						num_crossings += optim.optimize_layout(fix_xvars=xvars_to_fix)[1]
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
			for x_var in x_vars:
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

	def __add_3claw_constraints(self, m: gp.Model, g: LayeredGraph, cvars, c):
		if self.claw_constraints:
			cvset = set(cvars)
			bearclaws = motifs.get_3claws(g)
			print(f"3-claws found: {len(bearclaws)}")
			for claw in bearclaws:
				claw_cvs = gp.LinExpr()
				# claw_cvs += c[get_c_var(cvset, claw[3], claw[4])]
				claw_cvs += c[get_c_var(cvset, claw[3], claw[1])]
				# claw_cvs += c[get_c_var(cvset, claw[3], claw[5])]
				claw_cvs += c[get_c_var(cvset, claw[3], claw[2])]
				claw_cvs += c[get_c_var(cvset, claw[0], claw[4])]
				claw_cvs += c[get_c_var(cvset, claw[0], claw[5])]
				# claw_cvs += c[get_c_var(cvset, claw[4], claw[5])]
				claw_cvs += c[get_c_var(cvset, claw[4], claw[2])]
				claw_cvs += c[get_c_var(cvset, claw[1], claw[5])]
				m.addConstr(claw_cvs >= 1)

	def __add_dome_path_constraints(self, m: gp.Model, g: LayeredGraph, cvars, c, xvars, x):
		if self.dome_path_constraints:
			cvset = set(cvars)
			domes = motifs.get_domepaths(g)
			print(f"Dome-paths found: {len(domes)}")
			for dome in domes:
				if dome[0][0] == dome[1][0]:
					klc, kl1, kl2 = get_x_var_consts(xvars, dome[2][0], dome[0][0])
					kmc, km1, km2 = get_x_var_consts(xvars, dome[2][0], dome[3][0])
					lmc, lm1, lm2 = get_x_var_consts(xvars, dome[0][0], dome[3][0])
					cikjl1, cikjl2 = get_c_var(cvset, dome[2], dome[1])
					ciljm1, ciljm2 = get_c_var(cvset, dome[0], dome[3])
					m.addConstr(klc * x[kl1, kl2] - 2 * kmc * x[km1, km2] + lmc * x[lm1, lm2] - c[cikjl1, cikjl2] - c[ciljm1, ciljm2] + (1 - klc)//2 - (1 - kmc) + (1 - lmc)//2 <= 0)
					m.addConstr(-klc * x[kl1, kl2] + 2 * kmc * x[km1, km2] - lmc * x[lm1, lm2] - c[cikjl1, cikjl2] - c[ciljm1, ciljm2] - (1 - klc)//2 + (1 - kmc) - (1 - lmc)//2 <= 0)

	def __optimize_with_subgraph_reduction(self, n_partitions, cluster, top_level_g, crosses, contacts, stack_to_nodeset, node_to_stack):
		# TODO (later): remove all the parameters to standard_opt, replace non-top-level calls with creation of new optimizer object
		self.print_info.append("")

		subgraphs = [set(node.id for node in self.g.nodes if cluster[node.id] == i) for i in range(n_partitions)]
		top_level_subgraphs = {v: cluster[next(iter(stack_to_nodeset[v]))] for v in top_level_g.node_ids.keys()}
		# 2-subgraph optimization version
		# top_x_vars = {}
		# for node_list in top_level_g.layers.values():
		# 	for n1, n2 in itertools.combinations(node_list, 2):
		# 		set_x_var(top_x_vars, n1.name, n2.name, top_level_subgraphs[n1.name])
		# top_level_optval = top_level_g.num_edge_crossings_from_xvars_no_sl(top_x_vars)
		top_level_optval, top_x_vars = self.__optimize_layout_standard(graph_arg=top_level_g, return_x_vars=True, bendiness_reduction=False, is_subgraph=True, name="Collapsed graph", verbose=True)

		subg_x_var_colors = {}
		for top_x_var, val in top_x_vars.items():
			set_x_var(subg_x_var_colors, top_level_subgraphs[top_x_var[0]], top_level_subgraphs[top_x_var[1]], val)
			if len(subg_x_var_colors) == n_partitions * (n_partitions - 1) // 2:
				break

		# select clean graph (or random one, temporarily) - merge all other subgraphs
		fix_subg = random.choice(list(range(n_partitions)))
		if n_partitions > 2:
			for i in range(n_partitions):
				if not all(get_x_var(subg_x_var_colors, i, j) == 1 for j in itertools.chain(range(i), range(i+1, n_partitions))) and not all(get_x_var(subg_x_var_colors, i, j) == 0 for j in itertools.chain(range(i), range(i+1, n_partitions))):
					fix_subg = i
					break
		other_subg = set()
		for i in range(len(subgraphs)):
			if i != fix_subg:
				other_subg.update(subgraphs[i])
		subgraphs_merged = [other_subg, subgraphs[fix_subg]]

		t = time.time()
		vis.draw_graph(top_level_g, "interim", gravity=True, groups=top_level_subgraphs)
		vis.draw_graph(self.g, "overall", groups=cluster)
		f_subg_l = {node.layer: node.id for node in top_level_g.nodes if top_level_subgraphs[node.id] == fix_subg}
		fix_x_vars_for_merge = {}
		for x_var in top_x_vars:
			for low_n1 in stack_to_nodeset[x_var[0]]:
				for low_n2 in stack_to_nodeset[x_var[1]]:
					if cluster[low_n1] == fix_subg or cluster[low_n2] == fix_subg:
						set_x_var(self.x_var_assign, low_n1, low_n2, top_x_vars[x_var])
					elif top_level_g[x_var[0]].layer in f_subg_l and get_x_var(top_x_vars, f_subg_l[top_level_g[x_var[0]].layer], x_var[0]) != get_x_var(top_x_vars, f_subg_l[top_level_g[x_var[0]].layer], x_var[1]):
						set_x_var(fix_x_vars_for_merge, low_n1, low_n2, top_x_vars[x_var])

		sides = {}
		for contact_node in crosses:
			# if len(crosses[contact_node] > 1
			# figure out contact_sides to pass to call of optimize. contact_sides=1 => node fixed at top
			for contact_other in crosses[contact_node]:
				if cluster[contact_node] == fix_subg or cluster[contact_other] == fix_subg:
					if node_to_stack[contact_node] + 1 in top_level_subgraphs and top_level_subgraphs[node_to_stack[contact_node]] == top_level_subgraphs[node_to_stack[contact_node] + 1]:
						sides[contact_node] = 1 - get_x_var(top_x_vars, node_to_stack[contact_node] + 1, node_to_stack[contact_other])
					if node_to_stack[contact_other] - 1 in top_level_subgraphs and top_level_subgraphs[node_to_stack[contact_other]] == top_level_subgraphs[node_to_stack[contact_other] - 1]:
						sides[contact_other] = 1 - get_x_var(top_x_vars, node_to_stack[contact_other] - 1, node_to_stack[contact_node])

		# self.g.create_double_adj_list(forward_only=True)
		layered_subg_list = []
		x_vars_opt = None
		unconstrained_opt_vals = []
		opt_vals = []
		for i, subg in enumerate(subgraphs_merged):
			g_prime = LayeredGraph()
			layered_subg_list.append(g_prime)
			sides_prime = {}
			vars_to_fix = {}
			extra_node_closest_subg = {}
			if i == 0:
				for x_var in fix_x_vars_for_merge:
					vars_to_fix[x_var] = fix_x_vars_for_merge[x_var]

			for subg_node in subg:
				g_prime.add_node(self.g[subg_node].layer, idx=self.g[subg_node].id, is_anchor=self.g[subg_node].is_anchor_node)
			for subg_node in subg:
				if subg_node in sides:
					for cnode, clist in crosses.items():
						if cnode == subg_node:
							for cadj in clist:
								if cadj not in g_prime:
									g_prime.add_node(g_prime[subg_node].layer + 1, idx=cadj, stacked=True)
									extra_node_closest_subg[cadj] = cluster[subg_node]
								g_prime.add_edge(subg_node, cadj)
								sides_prime[cadj] = sides[subg_node]
						elif subg_node in clist:
							if cnode not in g_prime:
								g_prime.add_node(g_prime[subg_node].layer - 1, idx=cnode, stacked=True)
								extra_node_closest_subg[cnode] = cluster[subg_node]
							g_prime.add_edge(cnode, subg_node)
							sides_prime[cnode] = sides[subg_node]
				for adj_node in self.g.get_double_adj_list()[subg_node]:
					if adj_node in subg:
						g_prime.add_edge(subg_node, adj_node)
			gp_layers = g_prime.get_ids_by_layer()

			for nd, clr in extra_node_closest_subg.items():  # line the extra crossing nodes up with vars_to_fix
				for other in gp_layers[g_prime[nd].layer]:
					if nd != other and not g_prime[other].stacked and clr != cluster[other]:
						set_x_var(vars_to_fix, nd, other, get_x_var(subg_x_var_colors, clr, cluster[other]))
					elif nd != other and g_prime[other].stacked and clr != extra_node_closest_subg[other]:
						set_x_var(vars_to_fix, nd, other, get_x_var(subg_x_var_colors, clr, cluster[other]))

			for contact_node, x_val in sides_prime.items():
				for node in gp_layers[g_prime[contact_node].layer]:
					# if node != contact_node and len(self.g.double_adj_list[node]) >= 1 and (contact_node, node) not in vars_to_fix and (node, contact_node) not in vars_to_fix:
					if i == 0:
						if node != contact_node and (contact_node, node) not in vars_to_fix and (node, contact_node) not in vars_to_fix and (not g_prime[node].stacked or sides_prime[node] != sides_prime[contact_node]):
							vars_to_fix[contact_node, node] = x_val
					else:
						if node != contact_node and (g_prime[node].stacked and sides_prime[node] == sides_prime[contact_node]):
							vars_to_fix[contact_node, node] = get_x_var(x_vars_opt, contact_node, node)
						elif node != contact_node:
							vars_to_fix[contact_node, node] = x_val

			unconstrained_opt_val = self.__optimize_layout_standard(graph_arg=g_prime)[1]
			opt_val, x_vars_opt = self.__optimize_layout_standard(graph_arg=g_prime, is_subgraph=True, fix_x_vars=vars_to_fix, return_x_vars=True, name=f"subgraph {i + 1}", verbose=self.subg_verbose)
			if unconstrained_opt_val != opt_val:
				self.print_info.append(f"\tSubgraph {i+1} has {opt_val} crossings but could have as few as {unconstrained_opt_val}")
			else:
				self.print_info.append(f"\tSubgraph {i+1} is optimal")
			unconstrained_opt_vals.append(unconstrained_opt_val)
			opt_vals.append(opt_val)
			for x_var, val in x_vars_opt.items():
				set_x_var(self.x_var_assign, x_var[0], x_var[1], val)

			self.__sequential_br(graph_arg=g_prime, substitute_x_vars=x_vars_opt, subgraph_seq=True)
			vis.draw_graph(g_prime, f"interim_subg{i + 1}", groups=cluster)

		t = time.time() - t

		n_ec = self.g.num_edge_crossings_from_xvars_no_sl(self.x_var_assign)
		self.print_info.append(f"Total time to optimize and patch together subgraphs: {t}")
		self.print_info.append(f"Final edge crossing count: {n_ec}")
		self.print_info.append(f"{top_level_optval} crossings in collapsed graph, {sum(opt_vals)} in subgraphs, {n_ec - top_level_optval - sum(opt_vals)} from crossing edges")  # FIXME (later)

		self.__sequential_br()
		vis.draw_graph(self.g, "endpoints_highlight", groups=cluster)

		# LOOKAT return x_vars (or use self.x_var_assign?), n_crossings, subgraph LayeredGraph objects
		return n_ec, opt_vals, unconstrained_opt_vals, top_level_optval

	""" Return min possible num crossings, given by my formula """
	def __bound_on_optimality(self, subg_assign, top_lv_g, cr_num, cr_num_tuple, opt_cr_num_tuple, top_ov):
		# LOOKAT search reduced graph for butterflies which are inexcusable, to improve bound
		return sum(opt_cr_num_tuple)

	def __optimize_full_subgraph_algorithm(self):
		t_allotted = self.cutoff_time if self.cutoff_time > 0 else 60
		t_start = time.time()
		n_partitions = 2
		best_so_far = self.x_var_assign.copy()
		best_n_cr = self.g.num_edge_crossings()
		lower_bound = 0
		# while lower_bound != best_n_cr and time.time() - t_start < t_allotted:
		for i in range(3):
			# n_rep = max(6 - n_partitions, 1)
			n_rep = 1
			cur_rep = 1
			while cur_rep <= n_rep and lower_bound != best_n_cr and time.time() - t_start < t_allotted:
				print(f"{n_partitions} partitions attempt {cur_rep}, {time.time()-t_start} elapased.")
				cluster = list(SpectralClustering(n_clusters=n_partitions, assign_labels="discretize", affinity="precomputed").fit(self.g.adjacency_matrix()).labels_)
				cluster.insert(0, 0)
				top_level_g, crosses, contacts, stack_to_nodeset, node_to_stack = self.g.stacked_graph_from_subgraph_nodes(cluster)
				n_cr, o_v, uo_v, top_ov = self.__optimize_with_subgraph_reduction(n_partitions, cluster, top_level_g, crosses, contacts, stack_to_nodeset, node_to_stack)
				if n_cr < best_n_cr:
					best_so_far = self.x_var_assign
					best_n_cr = n_cr
					print(f"New best, {n_cr} crossings")
				new_lower_bound = self.__bound_on_optimality(cluster, top_level_g, n_cr, o_v, uo_v, top_ov)
				if new_lower_bound > lower_bound:
					lower_bound = new_lower_bound
					print(f"New lower bound on optimal #crossings: {lower_bound}")
				cur_rep += 1
			n_partitions += 1
		self.__sequential_br(substitute_x_vars=best_so_far)
		if lower_bound == best_n_cr:
			# we're optimal, baby
			print("wooooo")
		else:
			print("ran out of time ;(")
		return 0

	def optimize_with_starting_assignments(self, assigments):
		out = self.__optimize_layout_standard(use_top_level_params=True, fix_x_vars=assigments)
		if self.verbose:
			for string in self.print_info:
				print(string)
			self.print_info.clear()
		return out

	def __generate_random_vars_to_fix(self, n_vars):
		r_vars = random.sample(list(self.x_var_assign.keys()), n_vars)
		assignments = {r_var: random.randint(0, 1) for r_var in r_vars}
		return assignments

	def __normalized_loss(self, x_vars, correct_solution, g1=0, g2=1):
		correct_y_solution = self.__find_y_assignment_given_x(correct_solution)
		my_y_solution = self.__find_y_assignment_given_x(x_vars)
		loss = 0
		for nd, y_v in correct_y_solution.items():
			if y_v != my_y_solution[nd]:
				loss += g1 + g2 * (abs(y_v - my_y_solution[nd]) - 1)
		return loss

	def __find_y_assignment_given_x(self, x_vars):
		l_to_vlist = {lay: [v for v in x_vars if self.g[v[0]].layer == lay] for lay in range(1, self.g.n_layers + 1)}
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

	def __assign_x_given_y(self):
		for k in self.x_var_assign:
			if self.g[k[0]].y < self.g[k[1]].y:
				self.x_var_assign[k] = 0
			else:
				self.x_var_assign[k] = 1

	def optimize_layout(self, fix_xvars=None):
		if self.do_subg_reduction:
			out = self.__optimize_full_subgraph_algorithm()
		elif fix_xvars is not None:
			out = self.__optimize_layout_standard(use_top_level_params=True, fix_x_vars=fix_xvars)
		else:
			out = self.__optimize_layout_standard(use_top_level_params=True)
		if self.verbose:
			for string in self.print_info:
				print(string)
			self.print_info.clear()
		return out

	def just_bendiness_reduction(self):
		self.__assign_x_given_y()
		self.__sequential_br()
		vis.draw_graph(self.g, self.name)
