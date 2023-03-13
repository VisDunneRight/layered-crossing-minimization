import itertools, re, time
import os.path
import random
import gurobipy as gp
from gurobipy import GRB
from sklearn.cluster import SpectralClustering
from src import vis, reductions, motifs, type_conversions, read_data
from src.graph import *
from src.helpers import *


class LayeredOptimizer:
	def __init__(self, layered_graph, parameters=None):
		if parameters is None:
			parameters = {}
		assert os.path.isfile(layered_graph) or type(layered_graph) == LayeredGraph, "input needs to be a path to graph file or a LayeredGraph object"
		if type(layered_graph) == LayeredGraph:
			self.g = layered_graph
		else:
			self.g = read_data.read(layered_graph)
		self.x_var_assign = {x_v: 2 for n_l in self.g.get_names_by_layer().values() for x_v in itertools.combinations(n_l, 2)}
		self.bendiness_reduction = parameters["bendiness_reduction"] if "bendiness_reduction" in parameters else False
		self.gamma_1 = parameters["gamma_1"] if "gamma_1" in parameters else 1
		self.gamma_2 = parameters["gamma_2"] if "gamma_2" in parameters else 1
		self.m_val = parameters["m_val"] if "m_val" in parameters else max(len(layer) for layer in self.g.layers.values())
		self.sequential_bendiness = parameters["sequential_bendiness"] if "sequential_bendiness" in parameters else True
		# self.transitivity_constraints = parameters["transitivity_constraints"] if "transitivity_constraints" in parameters else False
		self.return_full_data = parameters["return_full_data"] if "return_full_data" in parameters else False
		self.cutoff_time = parameters["cutoff_time"] if "cutoff_time" in parameters else 0
		self.do_subg_reduction = parameters["do_subg_reduction"] if "do_subg_reduction" in parameters else False
		self.return_x_vars = parameters["return_x_vars"] if "return_x_vars" in parameters else False
		self.butterfly_reduction = parameters["butterfly_reduction"] if "butterfly_reduction" in parameters else False
		self.verbose = parameters["verbose"] if "verbose" in parameters else False
		self.subg_verbose = parameters["subg_verbose"] if "subg_verbose" in parameters else True
		self.draw_graph = parameters["draw_graph"] if "draw_graph" in parameters else False
		self.fix_one_var = parameters["fix_one_var"] if "fix_one_var" in parameters else False
		self.heuristic_start = parameters["heuristic_start"] if "heuristic_start" in parameters else False
		self.aggro_presolve = parameters["presolve"] if "presolve" in parameters else False
		self.mip_relax = parameters["mip_relax"] if "mip_relax" in parameters else False
		self.xvar_branch_priority = parameters["priority"] if "priority" in parameters else False
		# self.junger_ec = parameters["junger_ec"] if "junger_ec" in parameters else False
		self.junger_trans = parameters["junger_trans"] if "junger_trans" in parameters else False
		self.strat_big_m = parameters["strat_big_m"] if "strat_big_m" in parameters else False
		self.mirror_vars = parameters["mirror_vars"] if "mirror_vars" in parameters else False
		self.stratisfimal_y_vars = parameters["stratisfimal_yvars"] if "stratisfimal_yvars" in parameters else False
		self.symmetry_constraints = parameters["symmetry_constraints"] if "symmetry_constraints" in parameters else True
		# self.indicator_y_constraints = True
		self.return_experiment_data = parameters["return_experiment_data"] if "return_experiment_data" in parameters else False
		self.name = parameters["name"] if "name" in parameters else "graph1"
		self.print_info = []
		if not self.junger_trans and not self.strat_big_m and not self.stratisfimal_y_vars:
			self.strat_big_m = True

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

	def transitivity(self, model: gp.Model, names_by_layer, x_vars, x):
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

	def edge_crossings(self, model: gp.Model, c_vars, x, c, graph_arg=None, track_x_var_usage=False, butterflies=None):
		g = self.g if graph_arg is None else graph_arg
		n_constr_0, n_constr_1, n_constr_2 = 0, 0, 0
		x_var_usage = {}
		if butterflies == set():
			butterflies = None
		for c_var in c_vars:
			if butterflies is None or c_var not in butterflies:
				if self.mirror_vars and not g.edge_names[c_var[0]].same_layer_edge and not g.edge_names[c_var[1]].same_layer_edge:
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

	def edge_crossings_equal_butterfly(self, model: gp.Model, c_vars, x, c, graph_arg=None, track_x_var_usage=False, butterflies=None):
		g = self.g if graph_arg is None else graph_arg
		n_constr_0, n_constr_1, n_constr_2 = 0, 0, 0
		x_var_usage = {}
		if butterflies == set():
			butterflies = None
		for c_var_pair in butterflies:
			model.addConstr(c[c_var_pair[0]] + c[c_var_pair[1]] == 1)
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

	def edge_crossings_junger(self, model: gp.Model, c_vars, x_vars, x, c, graph_arg=None, track_x_var_usage=False, butterflies=None):
		g = self.g if graph_arg is None else graph_arg
		n_constraints = 0
		x_var_usage = {}
		for c_var in c_vars:
			if (butterflies is None or c_var not in butterflies) and c_var[0][0] != c_var[1][0]:
				if not self.mirror_vars:
					if c_var[0][1] < c_var[1][1]:
						x1const, x11, x12 = get_x_var_consts(x_vars, c_var[0][1], c_var[1][1])
						x2const, x21, x22 = get_x_var_consts(x_vars, c_var[0][0], c_var[1][0])
						model.addConstr(c[c_var] + x1const * x[x11, x12] - x2const * x[x21, x22] + (1 - x1const)//2  - (1 - x2const)//2 >= 0)
						model.addConstr(c[c_var] - x1const * x[x11, x12] + x2const * x[x21, x22] - (1 - x1const)//2  + (1 - x2const)//2 >= 0)
						if track_x_var_usage:
							if (x11, x12) not in x_var_usage:
								x_var_usage[x11, x12] = 0
							if (x21, x22) not in x_var_usage:
								x_var_usage[x21, x22] = 0
							x_var_usage[x11, x12] += 1
							x_var_usage[x21, x22] += 1
					else:
						x1const, x11, x12 = get_x_var_consts(x_vars, c_var[1][1], c_var[0][1])
						x2const, x21, x22 = get_x_var_consts(x_vars, c_var[0][0], c_var[1][0])
						model.addConstr(c[c_var] + x1const * x[x11, x12] + x2const * x[x21, x22] - 1  + (1 - x1const)//2  + (1 - x2const)//2 >= 0)
						model.addConstr(c[c_var] - x1const * x[x11, x12] - x2const * x[x21, x22] + 1  - (1 - x1const)//2  - (1 - x2const)//2 >= 0)
						if track_x_var_usage:
							if (x11, x12) not in x_var_usage:
								x_var_usage[x11, x12] = 0
							if (x21, x22) not in x_var_usage:
								x_var_usage[x21, x22] = 0
							x_var_usage[x11, x12] += 1
							x_var_usage[x21, x22] += 1
				else:
					if c_var[0][1] < c_var[1][1]:
						model.addConstr(c[c_var] + x[c_var[0][1], c_var[1][1]] - x[c_var[0][0], c_var[1][0]] >= 0)
						model.addConstr(c[c_var] - x[c_var[0][1], c_var[1][1]] + x[c_var[0][0], c_var[1][0]] >= 0)
						if track_x_var_usage:
							if (c_var[0][1], c_var[1][1]) not in x_var_usage:
								x_var_usage[c_var[0][1], c_var[1][1]] = 0
							if (c_var[0][0], c_var[1][0]) not in x_var_usage:
								x_var_usage[c_var[0][0], c_var[1][0]] = 0
							x_var_usage[c_var[0][1], c_var[1][1]] += 2
							x_var_usage[c_var[0][0], c_var[1][0]] += 2
					else:
						model.addConstr(c[c_var] + x[c_var[1][1], c_var[0][1]] + x[c_var[0][0], c_var[1][0]] - 1 >= 0)
						model.addConstr(c[c_var] - x[c_var[1][1], c_var[0][1]] - x[c_var[0][0], c_var[1][0]] + 1 >= 0)
						if track_x_var_usage:
							if (c_var[0][0], c_var[1][0]) not in x_var_usage:
								x_var_usage[c_var[0][0], c_var[1][0]] = 0
							if (c_var[1][1], c_var[0][1]) not in x_var_usage:
								x_var_usage[c_var[1][1], c_var[0][1]] = 0
							x_var_usage[c_var[0][0], c_var[1][0]] += 2
							x_var_usage[c_var[1][1], c_var[0][1]] += 2
				n_constraints += 2
		return n_constraints, track_x_var_usage

	def add_symmetry_constraints(self, model: gp.Model, x_vars, c_vars, x, c):
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

	def optimize_with_subgraph_reduction(self, n_partitions, cluster, top_level_g, crosses, contacts, stack_to_nodeset, node_to_stack):
		self.print_info.append("")

		subgraphs = [set(node.name for node in self.g.nodes if cluster[node.name] == i) for i in range(n_partitions)]
		top_level_subgraphs = {v: cluster[next(iter(stack_to_nodeset[v]))] for v in top_level_g.node_names.keys()}
		# 2-subgraph optimization version
		# top_x_vars = {}
		# for node_list in top_level_g.layers.values():
		# 	for n1, n2 in itertools.combinations(node_list, 2):
		# 		set_x_var(top_x_vars, n1.name, n2.name, top_level_subgraphs[n1.name])
		# top_level_optval = top_level_g.num_edge_crossings_from_xvars_no_sl(top_x_vars)
		top_level_optval, top_x_vars = self.optimize_layout_standard(graph_arg=top_level_g, return_x_vars=True, bendiness_reduction=False, is_subgraph=True, name="Collapsed graph", verbose=True)

		subg_x_var_colors = {}
		for top_x_var, val in top_x_vars.items():
			set_x_var(subg_x_var_colors, top_level_subgraphs[top_x_var[0]], top_level_subgraphs[top_x_var[1]], val)
			if len(subg_x_var_colors) == n_partitions * (n_partitions - 1) // 2:
				break

		# TODO select clean graph (or random one, temporarily)
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
		f_subg_l = {node.layer: node.name for node in top_level_g.nodes if top_level_subgraphs[node.name] == fix_subg}
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

		self.g.adj_list = self.g.create_double_adj_list(forward_only=True)
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
				g_prime.add_node(self.g[subg_node].layer, name=self.g[subg_node].name, is_anchor=self.g[subg_node].is_anchor_node)
			for subg_node in subg:
				if subg_node in sides:
					for cnode, clist in crosses.items():
						if cnode == subg_node:
							for cadj in clist:
								if cadj not in g_prime:
									g_prime.add_node(g_prime[subg_node].layer + 1, name=cadj, stacked=True)
									extra_node_closest_subg[cadj] = cluster[subg_node]
								g_prime.add_edge(subg_node, cadj)
								sides_prime[cadj] = sides[subg_node]
						elif subg_node in clist:
							if cnode not in g_prime:
								g_prime.add_node(g_prime[subg_node].layer - 1, name=cnode, stacked=True)
								extra_node_closest_subg[cnode] = cluster[subg_node]
							g_prime.add_edge(cnode, subg_node)
							sides_prime[cnode] = sides[subg_node]
				for adj_node in self.g.adj_list[subg_node]:
					if adj_node in subg:
						g_prime.add_edge(subg_node, adj_node)
			gp_layers = g_prime.get_names_by_layer()

			for nd, clr in extra_node_closest_subg.items():  # line the extra crossing nodes up with vars_to_fix
				for other in gp_layers[g_prime[nd].layer]:
					if nd != other and not g_prime[other].stacked and clr != cluster[other]:
						set_x_var(vars_to_fix, nd, other, get_x_var(subg_x_var_colors, clr, cluster[other]))
					elif nd != other and g_prime[other].stacked and clr != extra_node_closest_subg[other]:
						set_x_var(vars_to_fix, nd, other, get_x_var(subg_x_var_colors, clr, cluster[other]))

			for contact_node, x_val in sides_prime.items():
				for node in gp_layers[g_prime[contact_node].layer]:
					# if node != contact_node and len(self.g.adj_list[node]) >= 1 and (contact_node, node) not in vars_to_fix and (node, contact_node) not in vars_to_fix:
					if i == 0:
						if node != contact_node and (contact_node, node) not in vars_to_fix and (node, contact_node) not in vars_to_fix and (not g_prime[node].stacked or sides_prime[node] != sides_prime[contact_node]):
							vars_to_fix[contact_node, node] = x_val
					else:
						if node != contact_node and (g_prime[node].stacked and sides_prime[node] == sides_prime[contact_node]):
							vars_to_fix[contact_node, node] = get_x_var(x_vars_opt, contact_node, node)
						elif node != contact_node:
							vars_to_fix[contact_node, node] = x_val

			unconstrained_opt_val = self.optimize_layout_standard(graph_arg=g_prime, bendiness_reduction=False)[1]
			opt_val, x_vars_opt = self.optimize_layout_standard(graph_arg=g_prime, bendiness_reduction=False, is_subgraph=True, fix_x_vars=vars_to_fix, return_x_vars=True, name=f"subgraph {i+1}", verbose=self.subg_verbose)
			if unconstrained_opt_val != opt_val:
				self.print_info.append(f"\tSubgraph {i+1} has {opt_val} crossings but could have as few as {unconstrained_opt_val}")
			else:
				self.print_info.append(f"\tSubgraph {i+1} is optimal")
			unconstrained_opt_vals.append(unconstrained_opt_val)
			opt_vals.append(opt_val)
			for x_var, val in x_vars_opt.items():
				set_x_var(self.x_var_assign, x_var[0], x_var[1], val)

			self.sequential_br(graph_arg=g_prime, substitute_x_vars=x_vars_opt, subgraph_seq=True)
			vis.draw_graph(g_prime, f"interim_subg{i + 1}", groups=cluster)

		t = time.time() - t

		n_ec = self.g.num_edge_crossings_from_xvars_no_sl(self.x_var_assign)
		self.print_info.append(f"Total time to optimize and patch together subgraphs: {t}")
		self.print_info.append(f"Final edge crossing count: {n_ec}")
		self.print_info.append(f"{top_level_optval} crossings in collapsed graph, {sum(opt_vals)} in subgraphs, {n_ec - top_level_optval - sum(opt_vals)} from crossing edges")  # incorrect

		self.sequential_br()
		vis.draw_graph(self.g, "endpoints_highlight", groups=cluster)

		# TODO return x_vars (or use self.x_var_assign?), n_crossings, subgraph LayeredGraph objects
		return n_ec, opt_vals, unconstrained_opt_vals, top_level_optval

	def optimize_layout_standard(self, graph_arg=None, bendiness_reduction=False, assignment=None, return_x_vars=False, heuristic_start=False, transitivity=False, presolve=0, name="graph1", fix_x_vars=None, start_x_vars=None, fix_1_xvar=False, branch_on_x_vars=False, is_subgraph=False, verbose=False, use_top_level_params=False):
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
		# x_vars_layers = {}
		z_vars = []
		relax_type = GRB.INTEGER if not self.mip_relax else GRB.CONTINUOUS
		for i, name_list in nodes_by_layer.items():
			if self.mirror_vars:
				x_vars += list(itertools.permutations(name_list, 2))
			else:
				x_vars += list(itertools.combinations(name_list, 2))
			if self.stratisfimal_y_vars:
				z_vars += list(itertools.permutations(name_list, 2))
			# x_vars_layers[i] = list(itertools.combinations(name_list, 2))
		# x = m.addVars(x_vars, vtype=GRB.CONTINUOUS, name="x")
		# k = m.addVars(x_vars, vtype=GRB.CONTINUOUS, name="k")
		x = m.addVars(x_vars, vtype=GRB.BINARY, name="x")
		if self.stratisfimal_y_vars:
			z = m.addVars(z_vars, vtype=relax_type, lb=0, ub=self.m_val, name="z")
		c_vars, c_consts = reductions.normal_c_vars(g, edges_by_layer, self.mirror_vars)
		if self.mirror_vars:
			c_vars_orig, c_consts = reductions.normal_c_vars(g, edges_by_layer, False)
		c = m.addVars(c_vars, vtype=relax_type, name="c")
		if self.strat_big_m or self.stratisfimal_y_vars:
			y_vars = [n.name for n in g]
			y = m.addVars(y_vars, vtype=relax_type, lb=0, ub=self.m_val, name="y")
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

		""" Set variable starting values """
		if start_x_vars:
			for v in m.getVars():
				v.Start = start_x_vars[v.varName]
			# for k, v in start_x_vars.items():
			# 	if k in x:
			# 		m.getVarByName(f"x[{k[0]},{k[1]}]").Start = v
			# 	else:
			# 		m.getVarByName(f"x[{k[1]},{k[0]}]").Start = 1 - v

		if heuristic_start or (use_top_level_params and self.heuristic_start):
			g_igraph = type_conversions.layered_graph_to_igraph(g)
			heuristic_layout = g_igraph.layout_sugiyama(layers=g_igraph.vs["layer"])
			for i, coord in enumerate(heuristic_layout.coords[:g.n_nodes]):
				g[i + 1].y = coord[0]
			for v in m.getVars():
				if v.varName[:1] == "y":
					v.Start = g[int(v.varName[2:v.varName.index(']')])].y
				elif v.varName[:1] == "x":
					v.Start = 1 if g[int(v.varName[2:v.varName.index(',')])].y < g[int(v.varName[v.varName.index(',') + 1:v.varName.index(']')])].y else 0

		# g.barycentric_reordering(10)
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

		if branch_on_x_vars or (use_top_level_params and self.xvar_branch_priority):
			for v in m.getVars():
				if v.varName[:1] == "x":
					v.BranchPriority = 1
				else:
					v.BranchPriority = 0
	
		""" Butterfly reduction """
		butterfly_c_vars = set()
		butterfly_c_pairs = []
		if self.butterfly_reduction:
			b_set_list = []
			for b_v in motifs.get_butterflies(g):
				b_set_list.append(set(b_v))
			b_set_one_found = [False]*len(b_set_list)
			print("Butterfly set:", b_set_list)
			if b_set_list:
				self.print_info.append(f"{pre_sym}Num butterflies found: {len(b_set_list)}")
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
				opt += len(butterfly_c_vars) // 2
				m.setObjective(opt, GRB.MINIMIZE)
		else:
			opt = gp.LinExpr()
			if self.mirror_vars:
				for i, c_var in enumerate(c_vars_orig):
					opt += c_consts[i] * c[c_var]
			else:
				for i, c_var in enumerate(c_vars):
					opt += c_consts[i] * c[c_var]
			opt += len(butterfly_c_vars) // 2
			m.setObjective(opt, GRB.MINIMIZE)

		""" Transitivity constraints """
		# if transitivity or (use_top_level_params and self.transitivity_constraints):
		if self.junger_trans:
			self.transitivity(m, nodes_by_layer, x_vars, x)

		""" Long-version crossing reduction code """
		n_cs = self.edge_crossings(m, c_vars, x, c, graph_arg=g, track_x_var_usage=fix_1_xvar or (use_top_level_params and self.fix_one_var), butterflies=butterfly_c_vars)
		# n_cs = self.edge_crossings_2(m, c_vars, x, c, graph_arg=g, track_x_var_usage=fix_1_xvar or (use_top_level_params and self.fix_one_var), butterflies=butterfly_c_pairs)
		for i, val in enumerate(n_cs[:-1]):
			n_constraints_generated[i] += val

		""" Symmetry constraints """
		if self.mirror_vars and self.symmetry_constraints:
			self.add_symmetry_constraints(m, x_vars, c_vars, x, c)

		""" Fix key x-var """
		if fix_1_xvar or (use_top_level_params and self.fix_one_var):
			x_var_usage = n_cs[-1]
			if x_vars:
				if x_var_usage != {}:
					most_used_x = max(x_var_usage, key=x_var_usage.get)
				else:
					most_used_x = random.choice(x_vars)
				m.getVarByName(f"x[{most_used_x[0]},{most_used_x[1]}]").lb = 0
				m.getVarByName(f"x[{most_used_x[0]},{most_used_x[1]}]").ub = 0

		""" Vertical position, implication version """
		# Uses big-M method: https://support.gurobi.com/hc/en-us/articles/4414392016529-How-do-I-model-conditional-statements-in-Gurobi-
		if self.strat_big_m:
			for x_var in x_vars:
				m.addGenConstrIndicator(x[x_var], True, y[x_var[0]] + 1 <= y[x_var[1]])
				m.addGenConstrIndicator(x[x_var], False, y[x_var[0]] >= 1 + y[x_var[1]])
				n_constraints_generated[3] += 2
	
		""" Original vertical position constraints """
		if self.stratisfimal_y_vars:
			for x_var in x_vars:
				m.addConstr(z[x_var] - self.m_val * x[x_var] <= 0, f"1.1z{x_var}")
				m.addConstr(z[x_var] - y[x_var[0]] - (self.m_val * x[x_var]) >= -1 * self.m_val, f"1.2z{x_var}")
				m.addConstr(y[x_var[1]] - z[x_var] - x[x_var] >= 0, f"1.3z{x_var}")
				m.addConstr(z[x_var] <= y[x_var[0]], f"1.4z{x_var}")
				# m.addConstr(z[x_var] >= 0, f"1.5z{x_var}")	# print(f"z{x_var}-{self.m_val}*x{x_var} <= 0"), print(f"y{x_var[1]} - z{x_var} - x{x_var} >= 0")

				m.addConstr(z[x_var[1], x_var[0]] - self.m_val * (1 - x[x_var]) <= 0, f"2.1z{x_var}")
				m.addConstr(z[x_var[1], x_var[0]] - y[x_var[1]] - self.m_val * (1 - x[x_var]) >= -1 * self.m_val, f"2.2z{x_var}")
				m.addConstr(y[x_var[0]] - z[x_var[1], x_var[0]] - (1 - x[x_var]) >= 0, f"2.3z{x_var}")
				m.addConstr(z[x_var[1], x_var[0]] <= y[x_var[1]], f"2.4z{x_var}")
				# m.addConstr(z[x_var[1],x_var[0]] >= 0, f"2.5z{x_var}")	# print(f"z[{x_var[1],x_var[0]}]-{self.m_val}*(1-x{x_var}) <= 0"), print(f"y{x_var[0]} - z[{x_var[1],x_var[0]}] - (1-x{x_var}) >= 0")

				# m.addConstr(k[x_var] + x[x_var] - 1 >= 0)  # SOS-constraint version. Use with full LP relaxation
				# m.addSOS(GRB.SOS_TYPE1, [k[x_var], x[x_var]])

				n_constraints_generated[3] += 10

		""" Non-sequential bendiness reduction, original Stratisfimal version"""
		# if not self.sequential_bendiness and self.bendiness_reduction and not (transitivity or (use_top_level_params and self.transitivity_constraints)):
		if not self.sequential_bendiness and self.bendiness_reduction:
			for b_var in b_vars:
				m.addConstr(y[b_var[0]] - y[b_var[1]] <= b[b_var], f"1bend{b_var}")
				m.addConstr(y[b_var[1]] - y[b_var[0]] <= b[b_var], f"2bend{b_var}")
				n_constraints_generated[4] += 2

		""" Optimize model """
		if self.cutoff_time > 0:
			m.setParam("TimeLimit", self.cutoff_time)
		m.setParam("OutputFlag", 0)
		# m.setParam("LogFile", "standard_log")
		# m.setParam("LogToConsole", 0)
		if self.aggro_presolve:
			m.setParam("Presolve", 2)
		n_constraints_generated[5] = sum(n_constraints_generated[:5])
		t1 = time.time() - t1
		if verbose or (use_top_level_params and self.verbose):
			self.print_info.append(f"{pre_sym}Time to input constraints: {t1}")
			self.print_info.append(f"{pre_sym}Constraint counts: {n_constraints_generated}")
		t2 = time.time()
		m.optimize()
		t2 = time.time() - t2
		if verbose or (use_top_level_params and self.verbose):
			self.print_info.append(f"{pre_sym}Objective: {m.objVal}")
			self.print_info.append(f"{pre_sym}Time to optimize: {t2}")
		x_vars_opt = {}
		for v in m.getVars():
			if v.varName[:1] == "x" and use_top_level_params:
				self.x_var_assign[int(v.varName[2:v.varName.index(',')]), int(v.varName[v.varName.index(',') + 1:v.varName.index(']')])] = round(v.x)
			elif v.varName[:1] == "x":
				x_vars_opt[int(v.varName[2:v.varName.index(',')]), int(v.varName[v.varName.index(',') + 1:v.varName.index(']')])] = round(v.x)

		""" Draw pre-bendiness graph """
		# vis.draw_graph(g, "interim")
	
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
			# if transitivity or (use_top_level_params and self.transitivity_constraints):
			# 	if use_top_level_params:
			# 		g.assign_y_vals_given_x_vars(self.x_var_assign)
			# 	else:
			# 		g.assign_y_vals_given_x_vars(x_vars_opt)
			# for v in m.getVars():
			# 	if v.varName[:1] == "y":
			# 		g[int(v.varName[2:v.varName.index(']')])].y = float(v.x)

		if self.draw_graph:
			if not self.bendiness_reduction:
				if self.junger_trans:
					if use_top_level_params:
						g.assign_y_vals_given_x_vars(self.x_var_assign)
					else:
						g.assign_y_vals_given_x_vars(x_vars_opt)
				else:
					for v in m.getVars():
						if v.varName[:1] == "y":
							g[int(v.varName[2:v.varName.index(']')])].y = float(v.x)
			vis.draw_graph(g, "interim")

		if verbose or (use_top_level_params and self.verbose):
			self.print_info.append(f"{pre_sym}Final edge crossing count: {g.num_edge_crossings()}")
			self.print_info.append(f"{pre_sym}Number of constraints: {n_constraints_generated}")
			self.print_info.append(f"{pre_sym}{round(t1, 3)}, {round(t2, 3)}, {round(t3, 3)}, {round(t1 + t2 + t3, 3)}")

		print(m.objVal, round(m.runtime, 3))

		if use_top_level_params and self.return_x_vars:
			return int(m.objVal), self.x_var_assign
		elif return_x_vars:
			return int(m.objVal), x_vars_opt

		if self.return_experiment_data:
			if m.objVal == float('inf'):
				objv = "inf"
			else:
				objv = round(m.objVal)
			return len(x_vars), len(c_vars), m.numVars, m.numConstrs, objv, round(m.runtime, 4), round(m.work, 4), int(m.nodeCount), round(t1, 3)

		if self.return_full_data:
			return len(x_vars), len(c_vars), n_cs[0], motifs.count_butterflies(g), round(m.objVal), round(t1 + t2 + t3 + t3, 3), round(m.runtime, 3), round(m.work, 3), int(m.iterCount), round(t1, 3)

		return round(t1 + t2 + t3, 3), int(m.objVal)

	""" Return min possible num crossings, given by my formula """
	def bound_on_optimality(self, subg_assign, top_lv_g, cr_num, cr_num_tuple, opt_cr_num_tuple, top_ov):
		# TODO search reduced graph for butterflies which are inexcusable, to improve bound
		return sum(opt_cr_num_tuple)

	def optimize_full_subgraph_algorithm(self, t_allotted):
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
				n_cr, o_v, uo_v, top_ov = self.optimize_with_subgraph_reduction(n_partitions, cluster, top_level_g, crosses, contacts, stack_to_nodeset, node_to_stack)
				if n_cr < best_n_cr:
					best_so_far = self.x_var_assign
					best_n_cr = n_cr
					print(f"New best, {n_cr} crossings")
				new_lower_bound = self.bound_on_optimality(cluster, top_level_g, n_cr, o_v, uo_v, top_ov)
				if new_lower_bound > lower_bound:
					lower_bound = new_lower_bound
					print(f"New lower bound on optimal #crossings: {lower_bound}")
				cur_rep += 1
			n_partitions += 1
		self.sequential_br(substitute_x_vars=best_so_far)
		if lower_bound == best_n_cr:
			# we're optimal, baby
			print("wooooo")
		else:
			print("ran out of time ;(")
		return 0

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
			out = self.optimize_full_subgraph_algorithm(60)
		else:
			out = self.optimize_layout_standard(use_top_level_params=True)
		if self.verbose:
			for string in self.print_info:
				print(string)
			self.print_info.clear()
		return out
