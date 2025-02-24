import collections
import math
import os.path
import random
import time
import itertools
import gurobipy as gp
# import multiprocessing as mp
from gurobipy import GRB
# from sklearn.cluster import SpectralClustering
from src import vis, reductions, motifs, type_conversions, read_data
from src.graph import LayeredGraph, CollapsedGraph
from src.helpers import *
from src.neighborhood import *
from src.heuristics import improved_sifting


class LayeredOptimizer:
	def __init__(self, layered_graph, **kwargs):
		assert (type(layered_graph) == str and os.path.isfile(layered_graph)) or type(layered_graph) == LayeredGraph or type(layered_graph) == CollapsedGraph, "input needs to be a valid path or LayeredGraph object"
		if type(layered_graph) == LayeredGraph or type(layered_graph) == CollapsedGraph:
			self.g = layered_graph
		else:
			self.g = read_data.read(layered_graph)
		self.x_var_assign = {x_v: 2 for n_l in self.g.get_ids_by_layer().values() for x_v in itertools.combinations(n_l, 2)}
		self.crossing_minimization = kwargs.get("crossing_minimization", False)
		self.bendiness_reduction = kwargs.get("bendiness_reduction", False)
		self.gamma_1 = kwargs.get("gamma_1", 1)
		self.gamma_2 = kwargs.get("gamma_2", 1)
		self.m_val = kwargs.get("m_val", round(1.5 * max(len(lr) for lr in self.g.layers.values())))
		self.sequential_bendiness = kwargs.get("sequential_bendiness", False)
		self.local_opt = kwargs.get("local_opt", False)
		self.local_opt_heuristic = kwargs.get("local_opt_heuristic", "incremental")
		self.n_partitions = kwargs.get("n_partitions", -1)
		self.return_full_data = kwargs.get("return_full_data", False)
		self.cutoff_time = kwargs.get("cutoff_time", 0)
		self.do_subg_reduction = kwargs.get("do_subg_reduction", False)
		self.return_x_vars = kwargs.get("return_x_vars", False)
		self.butterfly_reduction = kwargs.get("butterfly_reduction", False)
		self.draw_graph = kwargs.get("draw_graph", False)
		self.symmetry_breaking = kwargs.get("symmetry_breaking", True)
		self.heuristic_start = kwargs.get("heuristic_start", False)
		self.aggro_presolve = kwargs.get("presolve", False)
		self.mip_relax = kwargs.get("mip_relax", False)
		self.xvar_branch_priority = kwargs.get("xvar_branch_priority", False)
		self.direct_transitivity = kwargs.get("direct_transitivity", False)
		self.vertical_transitivity = kwargs.get("vertical_transitivity", False)
		self.mirror_vars = kwargs.get("mirror_vars", False)
		self.stratisfimal_y_vars = kwargs.get("stratisfimal_y_vars", False)
		self.symmetry_constraints = kwargs.get("symmetry_constraints", True)
		self.cycle_constraints = kwargs.get("cycle_constraints", False)
		self.collapse_subgraphs = kwargs.get("collapse_subgraphs", False)
		self.collapse_leaves = kwargs.get("collapse_leaves", False)
		self.claw_constraints = kwargs.get("claw_constraints", False)
		self.dome_path_constraints = kwargs.get("dome_path_constraints", False)
		self.polyhedral_constraints = kwargs.get("polyhedral_constraints", False)
		self.grouping_constraints = kwargs.get("grouping_constraints", False)
		self.y_based_group_constraints = kwargs.get("y_based_group_constraints", False)
		self.node_emphasis = kwargs.get("node_emphasis", False)
		self.emphasis_cr_weight = kwargs.get("emphasis_cr_weight", 3)
		self.emphasis_br_weight = kwargs.get("emphasis_br_weight", 1)
		self.apply_node_weight_spacing = kwargs.get("apply_node_weight_spacing", False)
		self.apply_edge_weight = kwargs.get("apply_edge_weight", False)
		self.angular_resolution = kwargs.get("angular_resolution", False)
		self.symmetry_maximization = kwargs.get("symmetry_maximization", False)
		self.symmetry_maximization_edges = kwargs.get("symmetry_maximization_edges", False)
		self.min_max_crossings = kwargs.get("min_max_crossings", False)
		self.gamma_min_max = kwargs.get("gamma_min_max", 1)
		self.streamline = kwargs.get("streamline", False)
		self.anchor_proximity = kwargs.get("anchor_proximity", 0.3)
		self.fix_x_vars = kwargs.get("fix_x_vars", False)
		self.start_xy_vars = kwargs.get("start_xy_vars", False)
		self.fix_nodes = kwargs.get("fix_nodes", False)
		# self.only_min_max_crossings = kwargs.get("only_min_max_crossings", False)
		self.edge_bundling = kwargs.get("edge_bundling", False)
		self.gamma_bundle = kwargs.get("gamma_bundle", 1)
		self.fairness_constraints = kwargs.get("fairness_constraints", False)
		self.fairness_metric = kwargs.get("fairness_metric", "crossings")
		self.gamma_fair = kwargs.get("gamma_fair", 1)
		self.return_experiment_data = kwargs.get("return_experiment_data", False)
		self.create_video = kwargs.get("create_video", False)
		self.constrain_straight_long_arcs = kwargs.get("constrain_straight_long_arcs", False)
		self.long_arc_bend_limit = kwargs.get("long_arc_bend_limit", 0)
		self.record_solution_data_over_time = kwargs.get("record_solution_data_over_time", False)
		self.crossing_lower_constraints = kwargs.get("crossing_lower_constraints", False)
		self.name = kwargs.get("name", "graph1")
		self.nthreads = kwargs.get("nthreads", 0)
		self.hybrid_constraints = kwargs.get("hybrid_constraints", [])
		self.print_info = []
		if self.polyhedral_constraints:
			self.claw_constraints, self.dome_path_constraints = True, True

	def __optimize_layout_standard(self, graph_arg=None, fix_x_vars=None, start_x_vars=None):
		with gp.Env() as env, gp.Model(env=env) as m:
			g = self.g if graph_arg is None else graph_arg

			""" Collapse valid subgraphs """
			if self.collapse_leaves:
				g = g.collapse_leaves()
				self.x_var_assign = {x_v: 2 for n_l in g.get_ids_by_layer().values() for x_v in itertools.combinations(n_l, 2)}
			# vis.draw_graph(g, "after collapse", groups={nd.name: 1 if nd.stacked else 0 for nd in g.nodes}

			""" Get groups and add filler nodes if necessary """
			groups = self.__setup_filler_nodes(g)

			""" Create model """
			t1 = time.time()
			x_vars, c_vars = self.__create_optimization_model(m, g, fix_x_vars=fix_x_vars, start_x_vars=start_x_vars, groups=groups)
			# x_vars, c_vars = self.__crossing_reduction_model(m, g, fix_x_vars=fix_x_vars, start_x_vars=start_x_vars, groups=groups)
			t1 = time.time() - t1

			""" Optimize crossing minimization model """
			if self.record_solution_data_over_time:
				obj_val, t2, ncr_list, ntimes_list = self.__optimize_crossing_reduction_model(m, g, env, x_vars=x_vars)
			else:
				obj_val, t2 = self.__optimize_crossing_reduction_model(m, g, env, x_vars=x_vars)

			# vis.draw_graph(g, "interim")

			""" Sequential bendiness reduction """
			t3 = 0
			if not self.__we_need_y_vars():
				g.assign_y_vals_given_x_vars(self.x_var_assign)
			# if self.bendiness_reduction and self.sequential_bendiness:
			# 	t3 = time.time()
			# 	self.__sequential_br(graph_arg=g, env=env, groups=groups)
			# 	t3 = time.time() - t3
			# else:
			# 	t3 = 0
			# 	if self.direct_transitivity and not self.grouping_constraints:
			# 		g.assign_y_vals_given_x_vars(self.x_var_assign)
			# 	else:
			# 		for v in m.getVars():
			# 			if v.varName[:1] == "y" and int(v.varName[2:v.varName.index(']')]) in g.node_ids:
			# 				g[int(v.varName[2:v.varName.index(']')])].y = float(v.x)

			""" Adjust objects if leaves were collapsed/filler nodes added """
			self.__teardown_filler_nodes(g)
			if self.collapse_leaves:
				g = g.old_g

			""" Remove any additional filler group nodes """
			if self.grouping_constraints:
				reductions.remove_filler_group_nodes(g, self.x_var_assign)

			""" Draw the resulting graph """
			if self.draw_graph:
				if self.name == "collapsed_graph":
					vis.draw_graph(g, self.name, groups=[g.subgraphs[g.stack_node_to_nodelist[nod.id][0]] for nod in g.nodes])
				else:
					vis.draw_graph(g, self.name)

			""" Print results and return data """
			print(f"Optimization objective value: {obj_val}\t\tOptimization runtime: {round(m.runtime, 3)}s")
			print("g calculation \t", end='')
			g.calculate_stratisfimal_objective(1, 1)

			if self.return_x_vars:
				return obj_val, self.x_var_assign

			if self.return_experiment_data:
				return len(x_vars), len(c_vars), m.numVars, m.numConstrs, obj_val, m.runtime, m.status, int(m.nodeCount), round(t1, 3)

			if self.record_solution_data_over_time:
				return ncr_list, ntimes_list

			retval = collections.namedtuple("retval", "runtime objval status")
			return retval(t1 + t2 + t3, obj_val, m.status)

	def __crossing_reduction_model(self, m: gp.Model, g: LayeredGraph, fix_x_vars=None, start_x_vars=None, groups=None):
		if self.polyhedral_constraints:
			self.claw_constraints, self.dome_path_constraints = True, True
		if not self.direct_transitivity and not self.vertical_transitivity:
			self.direct_transitivity = True
		# if self.constrain_straight_long_arcs and not self.vertical_transitivity:
		# 	print("Using vertical transitivity instead to constrain long arcs")
		# 	self.vertical_transitivity, self.direct_transitivity = True, False
		if self.bendiness_reduction and self.direct_transitivity and not self.sequential_bendiness:
			self.vertical_transitivity, self.direct_transitivity = True, False

		nodes_by_layer = g.get_ids_by_layer()
		edges_by_layer = g.get_edge_ids_by_layer()

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
		c_vars_orig, nc_consts = None, None
		if self.mirror_vars:
			c_vars_orig, nc_consts = reductions.normal_c_vars(g, edges_by_layer, False)
		c = m.addVars(c_vars, vtype=relax_type, name="c")
		# if self.grouping_constraints:  # THIS IS FOR Y-VALUE BASED GROUP CONSTRAINTS
		# 	sl_groups, ml_groups = groups[0], groups[1]
		# 	grp_lb, grp_ub = [], []
		# 	if ml_groups:
		# 		grp_vars = list(range(len(ml_groups)))
		# 		grp_lb = m.addVars(grp_vars, vtype=GRB.CONTINUOUS, lb=0, ub=self.m_val, name="lb")
		# 		grp_ub = m.addVars(grp_vars, vtype=GRB.CONTINUOUS, lb=0, ub=self.m_val, name="ub")
		# 	if len(ml_groups) > 0 and not self.vertical_transitivity:
		# 		print("There are multilayer groups—swapping to vertical transitivity.")
		# 		self.vertical_transitivity, self.direct_transitivity = True, False
		y = None
		if self.vertical_transitivity or (self.bendiness_reduction and not self.sequential_bendiness):
			y = m.addVars([n.id for n in g], vtype=relax_type, lb=0, ub=self.m_val, name="y")
		b, b_vars = None, None
		if not self.sequential_bendiness and self.bendiness_reduction:
			b_vars = list(g.edge_ids.keys())
			b = m.addVars(b_vars, vtype=GRB.INTEGER, lb=0, ub=self.m_val, name="b")
		alpha, alpha_vars, a_aux_diff, bundle = None, None, None, None
		if self.edge_bundling:
			alpha_vars = [(a1, a2) for lid in g.layers for ix, a1 in enumerate(nodes_by_layer[lid]) if g[a1].is_anchor_node for a2 in nodes_by_layer[lid][ix+1:] if g[a2].is_anchor_node]
			alpha = m.addVars(alpha_vars, vtype=GRB.INTEGER, lb=0, ub=self.m_val, name="alpha")
			a_aux_diff = m.addVars(alpha_vars, vtype=GRB.INTEGER, lb=-self.m_val, ub=self.m_val, name="a_aux_diff")
			bundle = m.addVar(vtype=GRB.INTEGER, lb=0, name="bundle")
		e, e_vars = None, None
		if self.node_emphasis:
			require_graph_props(g, require_node_data=["emphasis"])
			e_vars = [(nd, nd_adj) if g[nd].layer < g[nd_adj].layer else (nd_adj, nd) for nd in g.node_data["emphasis"] for nd_adj in g.get_adj_list()[nd]]
			e = m.addVars(e_vars, vtype=GRB.INTEGER, lb=0, ub=self.m_val, name="e")
		cp, cp_vars, big_c = None, None, None
		if self.min_max_crossings:
			cp_vars = [(e.n1.id, e.n2.id) for e in g.edges]
			cp = m.addVars(cp_vars, vtype=GRB.INTEGER, lb=0, name="cp")
			big_c = m.addVar(lb=0, vtype=GRB.INTEGER, name="C")

		m.update()  # required after adding variables in order to use them in constraints

		""" Fix variables/set starting assignments """
		if fix_x_vars:
			for k, v in fix_x_vars.items():
				a1 = m.getVarByName(f"x[{k[0]},{k[1]}]")
				if a1:
					a1.lb, a1.ub = v, v
				else:
					a2 = m.getVarByName(f"x[{k[1]},{k[0]}]")
					a2.lb, a2.ub = 1 - v, 1 - v
		if start_x_vars:
			for v in m.getVars():
				v.Start = start_x_vars[v.varName]
		if not all((nd.fix == 0 for nd in g)):
			self.symmetry_breaking = False
			for nd in g:
				if nd.fix != 0:  # LOOKAT: make multiple fixed nodes in same layer not fix each other
					for nd2 in g.layers[nd.layer]:
						if (nd.id, nd2.id) in x:
							if nd.fix > nd2.fix:
								m.getVarByName(f"x[{nd.id},{nd2.id}]").lb, m.getVarByName(f"x[{nd.id},{nd2.id}]").ub = 0, 0
							elif nd.fix < nd2.fix:
								m.getVarByName(f"x[{nd.id},{nd2.id}]").lb, m.getVarByName(f"x[{nd.id},{nd2.id}]").ub = 1, 1
						elif (nd2.id, nd.id) in x:
							if nd.fix > nd2.fix:
								m.getVarByName(f"x[{nd2.id},{nd.id}]").lb, m.getVarByName(f"x[{nd2.id},{nd.id}]").ub = 1, 1
							elif nd.fix < nd2.fix:
								m.getVarByName(f"x[{nd2.id},{nd.id}]").lb, m.getVarByName(f"x[{nd2.id},{nd.id}]").ub = 0, 0

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

		""" Node emphasis constraints """
		e_consts = self.__emphasis_constraints(m, g, self.x_var_assign, x, e_vars, e, y, b, nodes_by_layer)

		""" Set model objective function """
		opt_func = self.__optimization_function(c, c_vars, c_vars_orig, b, c_consts, nc_consts, bundle, e, e_vars, e_consts, big_c)
		m.setObjective(opt_func, GRB.MINIMIZE)

		""" Transitivity constraints """
		self.__transitivity(m, nodes_by_layer, self.x_var_assign, x, y, z)

		""" Edge crossing constraints """
		xv_use = self.__edge_crossings(m, c_vars, x, c, graph_arg=g, track_x_var_usage=self.symmetry_breaking, butterflies=butterfly_c_pairs)

		""" 3-claw constraints """
		self.__add_3claw_constraints(m, g, c_vars, c)

		""" Dome path constraints """
		self.__add_dome_path_constraints(m, g, c_vars, c, x_vars, x)

		""" Symmetry constraints """
		self.__add_symmetry_constraints(m, x_vars, c_vars, x, c)

		""" Cycle constraints """
		self.__add_cycle_constraints(m, g, c_vars, c)

		""" Break symmetry by fixing key x-var """
		self.__symmetry_breaking(m, xv_use, x_vars)

		""" Long-edge constraints """
		self.__long_edge_constraints(m, g, self.x_var_assign, x, y, nodes_by_layer)

		""" Group constraints """
		self.__add_group_constraints(m, g, x_vars, x, groups, nodes_by_layer)

		""" Edge bundling constraints """
		self.__edge_bundling(m, g, alpha, a_aux_diff, alpha_vars, bundle, x, x_vars, nodes_by_layer)

		""" Min-max edge crossing constraints """
		self.__min_max_crossing_constraints(m, g, c_vars, c, cp_vars, cp, big_c, edges_by_layer)

		""" Non-sequential bendiness reduction, original Stratisfimal version"""
		if not self.sequential_bendiness and self.bendiness_reduction:
			for b_var in b_vars:
				m.addConstr(y[b_var[0]] - y[b_var[1]] <= b[b_var], f"1bend{b_var}")
				m.addConstr(y[b_var[1]] - y[b_var[0]] <= b[b_var], f"2bend{b_var}")

		return x_vars, c_vars

	def __optimize_crossing_reduction_model(self, m: gp.Model, g: LayeredGraph, env, x_vars=None):
		if self.cutoff_time > 0:
			m.setParam("TimeLimit", self.cutoff_time)
		m.setParam("OutputFlag", 0)
		if self.aggro_presolve:
			m.setParam("Presolve", 2)
		if self.nthreads > 0:
			m.setParam("Threads", self.nthreads)
			print("Threads =", m.Params.Threads)
		t2 = time.time()
		if self.record_solution_data_over_time:
			time_values, crossing_values = [], []
			def callbackfn(model, where):
				if where == GRB.Callback.MIPSOL:
					time_values.append(model.cbGet(GRB.Callback.RUNTIME))
					crossing_values.append(model.cbGet(GRB.Callback.MIPSOL_OBJBST))
			m.optimize(callbackfn)
		else:
			m.optimize()
		t2 = time.time() - t2
		if (m.status != 2 and m.status != 9) or m.SolCount == 0:
			print("model returned status code:", m.status)
			print("3: model unsolvable")
			print(f"4: model never found a feasible solution (#solutions found:{m.SolCount})")
			print("11: solve interrupted")
			print("otherwise check https://www.gurobi.com/documentation/current/refman/optimization_status_codes.html")
			return 0, 0, 0, 0, 0, float('inf'), m.status, 0, 0, "INCORRECT STATUS"
		# gs1, gs2 = 0, 0
		for v in m.getVars():
			if v.varName[:2] == "x[":
				xv1 = int(v.varName[2:v.varName.index(',')])
				xv2 = int(v.varName[v.varName.index(',') + 1:v.varName.index(']')])
				set_x_var(self.x_var_assign, xv1, xv2, round(v.x))
			elif v.varName[:2] == "y[":
				g[int(v.varName[2:v.varName.index(']')])].y = v.x
			# elif v.varName[:3] == "ang":
			# 	print(v.varName, v.x)
			# elif v.varName[:2] == "c[":
			# 	print(v.varName, v.x)
			# elif v.varName[:3] == "sym":
			# 	print(v.varName, v.x)
			# elif v.varName[:2] == "b[":
			# 	print(v.varName, v.x)
			# 	xv1 = int(v.varName[2:v.varName.index(',')])
			# 	xv2 = int(v.varName[v.varName.index(',') + 1:v.varName.index(']')])
			# 	if g.node_data["fairness"][xv1] == g.node_data["fairness"][xv2] == 0:
			# 		gs1 += v.x
			# 	elif g.node_data["fairness"][xv1] == g.node_data["fairness"][xv2] == 1:
			# 		gs2 += v.x
			# elif v.varName[:4] == "fair":
			# 	print(v.varName, v.x)
			# elif v.varName[:1] == "c" and round(v.x) != 0:
			# 	print(v.varName, v.x)
		# print(gs1, gs2)
		model_objval = m.objVal

		""" Optimize and merge collapsed subgraphs """
		g, t4, model_objval = self.__optimize_subgraphs(g, x_vars, model_objval)

		if self.record_solution_data_over_time:
			return model_objval, t2, crossing_values, time_values
		else:
			return model_objval, t2

	def __optimization_function(self, g, c, c_vars, c_vars_orig, b, b_vars, c_consts, nc_consts, alpha, e, e_vars, e_consts, big_c, fair_var, sym_n, sym_e, ang, ang_vars):
		opt = gp.LinExpr()
		c_mult, b_mult = self.gamma_1, self.gamma_2
		c_to_iter = c_vars_orig if self.mirror_vars and self.symmetry_constraints else c_vars
		use_consts = nc_consts if self.mirror_vars and self.symmetry_constraints else c_consts
		if self.crossing_minimization:
			for i, c_var in enumerate(c_to_iter):
				if self.node_emphasis:
					if c_var[0] in e_consts:
						use_consts[i] += self.emphasis_cr_weight
					if c_var[1] in e_consts:
						use_consts[i] += self.emphasis_cr_weight
				opt += c_mult * use_consts[i] * c[c_var]
		if not self.sequential_bendiness and self.bendiness_reduction:
			if self.apply_edge_weight:
				opt += sum(g.get_edge(b_v[0], b_v[1]).weight * b_mult * b[b_v] for b_v in b_vars)
			else:
				opt += b_mult * b.sum()
		if self.edge_bundling:
			opt += self.gamma_bundle * alpha.sum()
		if self.node_emphasis and e_vars is not None:
			for e_var in e_vars:
				opt += e_consts[e_var] * e[e_var]
		if self.min_max_crossings:
			opt += self.gamma_min_max * big_c
		if self.fairness_constraints:
			opt += self.gamma_fair * fair_var
		if self.symmetry_maximization:
			opt += sym_n.sum()
		if self.symmetry_maximization_edges:
			opt += sym_e.sum()
		if self.angular_resolution:
			opt += sum((ang[cv] * c[cv] for cv in ang_vars))
		return opt

	def __sequential_br(self, graph_arg=None, substitute_x_vars=None, env=None, streamline=True, groups=None):
		g = self.g if graph_arg is None else graph_arg
		x_var_opt = self.x_var_assign if substitute_x_vars is None else substitute_x_vars
		y_vars = list(g.node_ids)
		m2 = gp.Model(env=env) if env is not None else gp.Model()
		y = m2.addVars(y_vars, vtype=GRB.CONTINUOUS, lb=0, ub=self.m_val, name="y")
		if self.grouping_constraints and self.y_based_group_constraints:
			grp_vars = list(range(len(groups[0]) + len(groups[1])))
			y_t = m2.addVars(grp_vars, vtype=GRB.CONTINUOUS, lb=0, ub=self.m_val, name="t")
			y_b = m2.addVars(grp_vars, vtype=GRB.CONTINUOUS, lb=0, ub=self.m_val, name="b")
		m2.update()
		for v in m2.getVars():
			if v.varName[0] == "y":
				v.start = g[int(v.varName[2:v.varName.index(']')])].y
		b_vars = list(g.edge_ids.keys())
		b = m2.addVars(b_vars, vtype=GRB.CONTINUOUS, lb=0, ub=self.m_val, name="b")
		m2.setObjective(b.sum(), GRB.MINIMIZE)
		for var, val in x_var_opt.items():
			if val == 0:
				if streamline and (g[int(var[0])].is_anchor_node and g[int(var[1])].is_anchor_node):
					m2.addConstr(y[var[0]] >= 0.15 + y[var[1]], f"vert{var}")
				elif streamline and (g[int(var[0])].is_anchor_node or g[int(var[1])].is_anchor_node):
					m2.addConstr(y[var[0]] >= 0.3 + y[var[1]], f"vert{var}")
				else:
					m2.addConstr(y[var[0]] >= 1 + y[var[1]], f"vert{var}")
			else:
				if streamline and (g[int(var[0])].is_anchor_node and g[int(var[1])].is_anchor_node):
					m2.addConstr(y[var[0]] + 0.15 <= y[var[1]], f"vert{var}")
				elif streamline and (g[int(var[0])].is_anchor_node or g[int(var[1])].is_anchor_node):
					m2.addConstr(y[var[0]] + 0.3 <= y[var[1]], f"vert{var}")
				else:
					m2.addConstr(y[var[0]] + 1 <= y[var[1]], f"vert{var}")
		for b_var in b_vars:
			m2.addConstr(y[b_var[0]] - y[b_var[1]] <= b[b_var], f"1bend{b_var}")
			m2.addConstr(y[b_var[1]] - y[b_var[0]] <= b[b_var], f"2bend{b_var}")

		if self.grouping_constraints and self.y_based_group_constraints:
			for i, grp in enumerate(groups[0] + groups[1]):
				grp_layers = {}
				for nd in grp:
					m2.addConstr(y[nd] >= y_b[i])  # group nodes placed within boundaries
					m2.addConstr(y[nd] <= y_t[i])
					if g[nd].layer not in grp_layers:
						grp_layers[g[nd].layer] = []
					grp_layers[g[nd].layer].append(nd)
				grpmax = max((len(lylist) for lylist in grp_layers.values()))
				m2.addConstr(y_t[i] - y_b[i] + 1 == grpmax)  # enforce close groups
				for lid in grp_layers:  # find the closest node above & below group for each layer
					gp_lay_yv = g[grp_layers[lid][0]].y
					closest_below, closest_below_nd = -1, -1
					closest_above, closest_above_nd = self.m_val, -1
					for nd_ot in g.layers[lid]:
						if nd_ot.id not in grp_layers[lid]:
							if gp_lay_yv > nd_ot.y > closest_below:
								closest_below, closest_below_nd = nd_ot.y, nd_ot.id
							elif gp_lay_yv < nd_ot.y < closest_above:
								closest_above, closest_above_nd = nd_ot.y, nd_ot.id
					if closest_below_nd != -1:
						m2.addConstr(y_b[i] >= y[closest_below_nd] + 0.5)  # highest node below grp must be below lower bound, add buffer of 0.5
					if closest_above_nd != -1:
						m2.addConstr(y_t[i] <= y[closest_above_nd] - 0.5)  # lowest node above grp must be above upper bound, add buffer of 0.5

		m2.setParam("OutputFlag", 0)
		m2.optimize()
		for v in m2.getVars():
			if v.varName[:1] == "y":
				g[int(v.varName[2:v.varName.index(']')])].y = float(v.x)

	def __bendiness_constraints(self, m: gp.Model, y, b_vars, b, b_aux):
		if self.__we_need_b_vars():
			for b_var in b_vars:
				if b_aux is None:
					m.addConstr(y[b_var[0]] - y[b_var[1]] <= b[b_var], f"1bend{b_var}")
					m.addConstr(y[b_var[1]] - y[b_var[0]] <= b[b_var], f"2bend{b_var}")
				else:
					m.addConstr(b_aux[b_var] == y[b_var[0]] - y[b_var[1]], f"bend{b_var}")
					m.addGenConstrAbs(b[b_var], b_aux[b_var])

	def __transitivity(self, model: gp.Model, g: LayeredGraph, names_by_layer, x_vars, x, y, z):
		if self.direct_transitivity and not self.fix_x_vars:
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
			gap, gap_scale, max_gap = 1, 1, 3
			if self.apply_node_weight_spacing:
				max_wt = max((nd.weight for nd in g))
				if max_wt > 1:
					gap_scale = (max_gap - 1) / (max_wt - 1)
			if self.fix_x_vars:
				for var, val in self.x_var_assign.items():
					if self.apply_node_weight_spacing:
						gap = gap_scale * (g[int(var[0])].weight + g[int(var[1])].weight) / 2 - gap_scale + 1
					if val == 0:
						if self.streamline and (g[int(var[0])].is_anchor_node and g[int(var[1])].is_anchor_node):
							model.addConstr(y[var[0]] >= (self.anchor_proximity / 2) * gap + y[var[1]], f"vert{var}")
						elif self.streamline and (g[int(var[0])].is_anchor_node or g[int(var[1])].is_anchor_node):
							model.addConstr(y[var[0]] >= max(self.anchor_proximity, 0.3) * gap + y[var[1]], f"vert{var}")
						else:
							model.addConstr(y[var[0]] >= gap + y[var[1]], f"vert{var}")
					else:
						if self.streamline and (g[int(var[0])].is_anchor_node and g[int(var[1])].is_anchor_node):
							model.addConstr(y[var[0]] + (self.anchor_proximity / 2) * gap <= y[var[1]], f"vert{var}")
						elif self.streamline and (g[int(var[0])].is_anchor_node or g[int(var[1])].is_anchor_node):
							model.addConstr(y[var[0]] + max(self.anchor_proximity, 0.3) * gap <= y[var[1]], f"vert{var}")
						else:
							model.addConstr(y[var[0]] + gap <= y[var[1]], f"vert{var}")
			else:
				for x_var in x_vars:
					if self.apply_node_weight_spacing:
						gap = gap_scale * (g[int(x_var[0])].weight + g[int(x_var[1])].weight) / 2 - gap_scale + 1
					if self.streamline and (g[int(x_var[0])].is_anchor_node and g[int(x_var[1])].is_anchor_node):
						model.addGenConstrIndicator(x[x_var], True, y[x_var[0]] + (self.anchor_proximity / 2) * gap <= y[x_var[1]])
						model.addGenConstrIndicator(x[x_var], False, y[x_var[0]] >= (self.anchor_proximity / 2) * gap + y[x_var[1]])
					elif self.streamline and (g[int(x_var[0])].is_anchor_node or g[int(x_var[1])].is_anchor_node):
						model.addGenConstrIndicator(x[x_var], True, y[x_var[0]] + max(self.anchor_proximity, 0.3) * gap <= y[x_var[1]])
						model.addGenConstrIndicator(x[x_var], False, y[x_var[0]] >= max(self.anchor_proximity, 0.3) * gap + y[x_var[1]])
					else:
						model.addGenConstrIndicator(x[x_var], True, y[x_var[0]] + gap <= y[x_var[1]])
						model.addGenConstrIndicator(x[x_var], False, y[x_var[0]] >= gap + y[x_var[1]])

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

	def __edge_crossings(self, model: gp.Model, g: LayeredGraph, c_vars, x, c, track_x_var_usage=False, butterflies=None):
		x_var_usage = {}
		if c:
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
						if self.crossing_lower_constraints or (self.fairness_constraints and self.fairness_metric == "crossings") or any(ele[0] == "crossing_fairness" for ele in self.hybrid_constraints):
							model.addConstr(- x1_rev * x[x1] - x2_rev * x[x2] + c[c_var] - (1 - x1_rev) / 2 - (1 - x2_rev) / 2 <= 0, f"3se{c_var}")
							model.addConstr(-(1 - x1_rev * x[x1]) - (1 - x2_rev * x[x2]) + c[c_var] + (1 - x1_rev) / 2 + (1 - x2_rev) / 2 <= 0, f"4se{c_var}")

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
						model.addConstr(c[c_var] + (1 - x1_rev * x[x1]) + x4_rev * x[x4] + (1 - x2_rev * x[x2]) - (1 - x1_rev) / 2 + (1 - x4_rev) / 2 - (1 - x2_rev) / 2 >= 1, f"1sl{c_var}")
						model.addConstr(c[c_var] + (1 - x3_rev * x[x3]) + x2_rev * x[x2] + (1 - x4_rev * x[x4]) - (1 - x3_rev) / 2 + (1 - x2_rev) / 2 - (1 - x4_rev) / 2 >= 1, f"2sl{c_var}")
						model.addConstr(c[c_var] + x1_rev * x[x1] + (1 - x3_rev * x[x3]) + x2_rev * x[x2] + (1 - x1_rev) / 2 - (1 - x3_rev) / 2 + (1 - x2_rev) / 2 >= 1, f"3sl{c_var}")
						model.addConstr(c[c_var] + x4_rev * x[x4] + (1 - x2_rev * x[x2]) + x3_rev * x[x3] + (1 - x4_rev) / 2 - (1 - x2_rev) / 2 + (1 - x3_rev) / 2 >= 1, f"4sl{c_var}")
						model.addConstr(c[c_var] + (1 - x4_rev * x[x4]) + x1_rev * x[x1] + (1 - x3_rev * x[x3]) - (1 - x4_rev) / 2 + (1 - x1_rev) / 2 - (1 - x3_rev) / 2 >= 1, f"5sl{c_var}")
						model.addConstr(c[c_var] + (1 - x2_rev * x[x2]) + x3_rev * x[x3] + (1 - x1_rev * x[x1]) - (1 - x2_rev) / 2 + (1 - x3_rev) / 2 - (1 - x1_rev) / 2 >= 1, f"6sl{c_var}")
						model.addConstr(c[c_var] + x2_rev * x[x2] + (1 - x4_rev * x[x4]) + x1_rev * x[x1] + (1 - x2_rev) / 2 - (1 - x4_rev) / 2 + (1 - x1_rev) / 2 >= 1, f"7sl{c_var}")
						model.addConstr(c[c_var] + x3_rev * x[x3] + (1 - x1_rev * x[x1]) + x4_rev * x[x4] + (1 - x3_rev) / 2 - (1 - x1_rev) / 2 + (1 - x4_rev) / 2 >= 1, f"8sl{c_var}")
		return x_var_usage

	def __symmetry_breaking(self, model: gp.Model, x_var_usage, x_vars):
		if self.symmetry_breaking:
			if x_vars:
				if x_var_usage != {}:
					most_used_x = max(x_var_usage, key=x_var_usage.get)
				else:
					most_used_x = random.choice(x_vars)
				val_to_set = get_x_var(self.x_var_assign, most_used_x[0], most_used_x[1])
				if val_to_set != 0 and val_to_set != 1:
					val_to_set = 0
				model.getVarByName(f"x[{most_used_x[0]},{most_used_x[1]}]").lb = val_to_set
				model.getVarByName(f"x[{most_used_x[0]},{most_used_x[1]}]").ub = val_to_set

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
			improved_sifting(graph)
			# g_igraph = type_conversions.layered_graph_to_igraph(graph)
			# heuristic_layout = g_igraph.layout_sugiyama(layers=g_igraph.vs["layer"])
			# for i, coord in enumerate(heuristic_layout.coords[:graph.n_nodes]):
			# 	graph.nodes[i].y = coord[0]
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

	def __optimize_subgraphs(self, graph, x_vars, num_crossings):
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
			# to_remove = []
			for x_var in x_vars:
				if graph[x_var[0]].stacked and graph[x_var[1]].stacked:
					for otv in graph.stack_node_to_nodelist[x_var[0]]:
						for otv2 in graph.stack_node_to_nodelist[x_var[1]]:
							set_x_var(self.x_var_assign, otv, otv2, get_x_var(self.x_var_assign, x_var[0], x_var[1]))
					# to_remove.append(x_var)
				elif graph[x_var[0]].stacked:
					for otv in graph.stack_node_to_nodelist[x_var[0]]:
						if x_var[1] != otv:
							set_x_var(self.x_var_assign, otv, x_var[1], get_x_var(self.x_var_assign, x_var[0], x_var[1]))
					# to_remove.append(x_var)
				elif graph[x_var[1]].stacked:
					for otv in graph.stack_node_to_nodelist[x_var[1]]:
						if x_var[0] != otv:
							set_x_var(self.x_var_assign, x_var[0], otv, get_x_var(self.x_var_assign, x_var[0], x_var[1]))
					# to_remove.append(x_var)
			# for xv in to_remove:
			# 	if xv in self.x_var_assign:
			# 		del self.x_var_assign[xv]
			# 	else:
			# 		del self.x_var_assign[xv[1], xv[0]]
			# print(to_remove)
			# print(self.x_var_assign)

			return graph.old_g, time.time() - t4, num_crossings
		else:
			return graph, 0, num_crossings

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

	def __long_edge_constraints(self, m: gp.Model, g: LayeredGraph, long_e, x_vars, x, y, n_b_l):
		if self.constrain_straight_long_arcs:
			long_edges = g.get_long_edges()
			l_e_counts = {}
			for l_e in long_edges:
				if l_e[0] not in l_e_counts:
					l_e_counts[l_e[0]] = 0
				if l_e[-1] not in l_e_counts:
					l_e_counts[l_e[-1]] = 0
				l_e_counts[l_e[0]] += 1
				l_e_counts[l_e[-1]] += 1
			if y is not None:
				for l_e_list in long_edges:
					left, right = l_e_list[0], l_e_list[1]
					idx = 0
					l_e_sum = gp.LinExpr()
					while idx < len(l_e_list) - 2:
						m.addConstr(self.m_val * long_e[left, right] + y[left] - y[right] >= 0)
						m.addConstr(self.m_val * long_e[left, right] - y[left] + y[right] >= 0)
						if 0 < idx or l_e_counts[l_e_list[0]] <= 1 or self.long_arc_bend_limit >= 1:
							l_e_sum += long_e[left, right]
						idx += 1
						left = right
						right = l_e_list[idx + 1]
					m.addConstr(self.m_val * long_e[left, right] + y[left] - y[right] >= 0)
					m.addConstr(self.m_val * long_e[left, right] - y[left] + y[right] >= 0)
					if l_e_counts[l_e_list[-1]] <= 1 or self.long_arc_bend_limit >= 2:
						l_e_sum += long_e[left, right]
					m.addConstr(l_e_sum <= self.long_arc_bend_limit)
				# else:
				# 	for l_e_list in long_edges:
				# 		if l_e_counts[l_e_list[0]] <= 1:
				# 			m.addConstr(y[l_e_list[0]] == y[l_e_list[1]])
				# 		left, right = l_e_list[1], l_e_list[2]
				# 		idx = 1
				# 		while right != l_e_list[-1]:
				# 			idx += 1
				# 			m.addConstr(y[left] == y[right])
				# 			left = right
				# 			right = l_e_list[idx]
				# 		if l_e_counts[right] <= 1:
				# 			m.addConstr(y[left] == y[right])
			else:
				require_graph_props(g, require_node_data=["filler"])
				for l_e_list in long_edges:
					left = l_e_list[0]
					left_xsum = calc_x_var_sum(left, n_b_l[g[left].layer], x_vars, x)
					right = l_e_list[1]
					right_xsum = calc_x_var_sum(right, n_b_l[g[right].layer], x_vars, x)
					if l_e_counts[left] <= 1:
						m.addConstr(left_xsum == right_xsum)
					left, left_xsum = right, right_xsum
					right = l_e_list[2]
					idx = 2
					right_xsum = calc_x_var_sum(right, n_b_l[g[right].layer], x_vars, x)
					while right != l_e_list[-1]:
						idx += 1
						m.addConstr(left_xsum == right_xsum)
						left, left_xsum = right, right_xsum
						right = l_e_list[idx]
						right_xsum = calc_x_var_sum(right, n_b_l[g[right].layer], x_vars, x)
					if l_e_counts[right] <= 1:
						m.addConstr(left_xsum == right_xsum)

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

	def __edge_bundling(self, m: gp.Model, g: LayeredGraph, a, a_vars, x, x_vars, y, n_b_l):
		if self.edge_bundling or any(ele[0] == "edge_bundles" for ele in self.hybrid_constraints):
			if self.__we_need_y_vars():
				for a1, a2 in a_vars:
					m.addConstr(self.m_val * a[a1, a2] + y[a1] - y[a2] + self.anchor_proximity / 2 >= 0)
					m.addConstr(self.m_val * a[a1, a2] - y[a1] + y[a2] + self.anchor_proximity / 2 >= 0)
			else:
				for a1, a2 in a_vars:
					lid = g[a1].layer
					a1_xsum, a2_xsum = gp.LinExpr(), gp.LinExpr()
					for nd in n_b_l[lid]:
						if not g[nd].is_anchor_node:
							x_r, u1, u2 = get_x_var_consts(x_vars, a1, nd)
							x_r2, v1, v2 = get_x_var_consts(x_vars, a2, nd)
							a1_xsum += x_r * x[u1, u2] + (1 - x_r) // 2
							a2_xsum += x_r2 * x[v1, v2] + (1 - x_r2) // 2
					m.addConstr(self.m_val * a[a1, a2] >= a2_xsum - a1_xsum)
					m.addConstr(self.m_val * a[a1, a2] >= a1_xsum - a2_xsum)
				# 	m.addConstr(a_aux_diff[a1, a2] == a2_xsum - a1_xsum)
				# 	m.addGenConstrAbs(a[a1, a2], a_aux_diff[a1, a2])
				# m.addGenConstrNorm(bundle, a, 0)

	def __fairness_constraints(self, m: gp.Model, g: LayeredGraph, c_vars, c, b_vars, b, fair_var):
		if self.fairness_constraints or any(ele[0] == "crossing_fairness" or ele[0] == "bend_fairness" for ele in self.hybrid_constraints):
			if "fairness" not in g.node_data:
				raise Exception("Need to add fairness assignments on nodes, use g.add_fairness_values().")
			total_zeros = sum(1 for v in g.node_data["fairness"].values() if v == 0)
			total_ones = sum(1 for v in g.node_data["fairness"].values() if v == 1)
			f0_csum, f1_csum = gp.LinExpr(), gp.LinExpr()
			if self.fairness_metric == "crossings" or any(ele[0] == "crossing_fairness" for ele in self.hybrid_constraints):
				for c_var in c_vars:
					cv_total = sum(g.node_data["fairness"][nd] for nd in (c_var[0][0], c_var[0][1], c_var[1][0], c_var[1][1]))
					if cv_total <= 1:
						f0_csum += c[c_var]
					elif cv_total >= 3:
						f1_csum += c[c_var]
			elif self.fairness_metric == "bends" or any(ele[0] == "bend_fairness" for ele in self.hybrid_constraints):
				for b_var in b_vars:
					if g.node_data["fairness"][b_var[0]] == 0 and g.node_data["fairness"][b_var[1]] == 0:
						f0_csum += b[b_var]
					elif g.node_data["fairness"][b_var[0]] == 1 and g.node_data["fairness"][b_var[1]] == 1:
						f1_csum += b[b_var]
			else:
				raise Exception(f"'{self.fairness_metric}' not supported with fairness optimization.\nAllowed metrics: 'crossings', 'bends'\nProceeding without fairness.")
			f1_csum = f1_csum * total_zeros / total_ones
			m.addConstr(fair_var >= f0_csum - f1_csum)
			m.addConstr(fair_var >= f1_csum - f0_csum)

	def __emphasis_constraints(self, m: gp.Model, g: LayeredGraph, x_vars, x, e_vars, e, y, b, n_b_l):
		if self.node_emphasis:
			if "emphasis" not in g.node_data:
				raise Exception("Need to add indicate which nodes to emphasize, use g.add_node_emphasis().")
			emph_consts = {}
			if y is None:
				require_graph_props(g, require_node_data=["filler"])
				for nd in g.node_data["emphasis"]:
					nd_xsum = gp.LinExpr()
					for nd_ot in n_b_l[g[nd].layer]:
						if nd_ot != nd:
							x_r, u1, u2 = get_x_var_consts(x_vars, nd, nd_ot)
							nd_xsum += x_r * x[u1, u2] + (1 - x_r) // 2
					for nd_adj in g.get_adj_list()[nd]:
						adj_xsum = gp.LinExpr()
						for adj_ot in n_b_l[g[nd_adj].layer]:
							x_r, u1, u2 = get_x_var_consts(x_vars, nd_adj, adj_ot)
							adj_xsum += x_r * x[u1, u2] + (1 - x_r) // 2
						e1 = nd if g[nd].layer < g[nd_adj].layer else nd_adj
						e2 = nd_adj if g[nd].layer < g[nd_adj].layer else nd
						if (e1, e2) not in emph_consts:
							emph_consts[e1, e2] = 0
						emph_consts[e1, e2] += 1
						m.addConstr(e[e1, e2] >= nd_xsum - adj_xsum)
						m.addConstr(e[e1, e2] >= adj_xsum - nd_xsum)
			else:
				for nd in g.node_data["emphasis"]:
					for nd_adj in g.get_adj_list()[nd]:
						e1 = nd if g[nd].layer < g[nd_adj].layer else nd_adj
						e2 = nd_adj if g[nd].layer < g[nd_adj].layer else nd
						if (e1, e2) not in emph_consts:
							emph_consts[e1, e2] = 0
						emph_consts[e1, e2] += 1
						if b is None:
							m.addConstr(e[e1, e2] >= y[e1] - y[e2])
							m.addConstr(e[e1, e2] >= y[e2] - y[e1])
						else:
							m.addConstr(e[e1, e2] >= b[e1, e2])
							m.addConstr(e[e1, e2] >= b[e1, e2])
			return emph_consts
		return None

	def __fix_specific_nodes(self, m: gp.Model, g: LayeredGraph, x, n_b_l):
		if self.fix_nodes:
			if "fix_nodes" not in g.node_data:
				raise Exception("Need to add indicate which nodes to fix, use g.add_node_fix().")
			for nid, loc in g.node_data["fix_nodes"].items():
				if loc == "top" or loc == "bottom":
					fx_val = 1 if loc == "top" else 0
					for nd_ot in n_b_l[g[nid].layer]:
						if nd_ot not in g.node_data["fix_nodes"] or g.node_data["fix_nodes"][nd_ot] != loc:
							self.__fix_x_var(m, (nid, nd_ot), fx_val)
				# else:  # loc == 'middle'
				# 	sum_ot = sum(1 for v in n_b_l[g[nid].layer] if v not in g.node_data["fix_nodes"] or g.node_data["fix_nodes"][v] != "middle")
				# 	xsum = LinExpr()
				# 	for nd_ot in n_b_l[g[nid].layer]:
				# 		if nd_ot != nid and (nd_ot not in g.node_data["fix_nodes"] or g.node_data["fix_nodes"][nd_ot] != "middle"):
				# 			x_r, u1, u2 = get_x_var_consts(self.x_var_assign, nid, nd_ot)
				# 			xsum += x_r * x[u1, u2] + (1 - x_r) // 2
				# 	if sum_ot % 2 == 0:
				# 		m.addConstr(xsum == sum_ot // 2)
				# 	else:
				# 		m.addConstr(xsum >= sum_ot // 2)
				# 		m.addConstr(xsum <= sum_ot // 2 + 1)

	def __min_max_crossing_constraints(self, m: gp.Model, g: LayeredGraph, c_vars, c, cp_vars, cp, big_c_var, e_b_l):
		if self.min_max_crossings or any(ele[0] == "min_max_crossings" for ele in self.hybrid_constraints):
			cvset = set(c_vars)
			for cp_var in cp_vars:
				cp_esum = gp.LinExpr()
				nd1l = g[cp_var[0]].layer if g[cp_var[0]].layer < g[cp_var[1]].layer else g[cp_var[1]].layer
				for ed_ot in e_b_l[nd1l]:
					if len({ed_ot[0], ed_ot[1], cp_var[0], cp_var[1]}) == 4:
						e1, e2 = get_c_var(cvset, cp_var, ed_ot)
						cp_esum += c[e1, e2]
				m.addConstr(cp[cp_var] == cp_esum)
				m.addConstr(cp[cp_var] <= big_c_var)

	def __setup_filler_nodes(self, g: LayeredGraph):
		sl_groups, ml_groups = None, None
		add_regular_filler_nodes = not self.__we_need_y_vars()
		if self.grouping_constraints:
			sl_groups, ml_groups, recovered_nodes = reductions.get_groups(g, add_regular_filler_nodes)
			if recovered_nodes:
				if "filler_xvs" in g.node_data and add_regular_filler_nodes:
					for x_v, val in g.node_data["filler_xvs"]:
						self.x_var_assign[x_v] = val
				if "filler_group_xvs" in g.node_data:
					for x_v, val in g.node_data["filler_group_xvs"]:
						self.x_var_assign[x_v] = val
			else:
				self.x_var_assign = {x_v: 2 for n_l in self.g.get_ids_by_layer().values() for x_v in itertools.combinations(n_l, 2)}
		elif (self.constrain_straight_long_arcs or self.node_emphasis) and not self.__we_need_y_vars():
			recovered_nodes = reductions.add_filler_nodes(g)
			if recovered_nodes:
				if "filler_xvs" in g.node_data and add_regular_filler_nodes:
					for x_v, val in g.node_data["filler_xvs"]:
						self.x_var_assign[x_v] = val
				if "filler_group_xvs" in g.node_data:
					for x_v, val in g.node_data["filler_group_xvs"]:
						self.x_var_assign[x_v] = val
			else:
				self.x_var_assign = {x_v: 2 for n_l in self.g.get_ids_by_layer().values() for x_v in itertools.combinations(n_l, 2)}
		return sl_groups, ml_groups

	def __teardown_filler_nodes(self, g: LayeredGraph):
		if self.grouping_constraints:
			reductions.remove_filler_nodes(g, self.x_var_assign)
		elif (self.constrain_straight_long_arcs or self.node_emphasis) and not self.__we_need_y_vars():
			reductions.remove_filler_nodes(g, self.x_var_assign)

	def __add_group_constraints(self, m: gp.Model, g: LayeredGraph, xvars, x, y, y_top, y_bottom, groups, n_b_l):
		if self.grouping_constraints:
			sl_groups, ml_groups = groups[0], groups[1]
			if self.y_based_group_constraints:
				for i, grp in enumerate(groups[0] + groups[1]):
					grp_layers = {}
					for nd in grp:
						m.addConstr(y[nd] >= y_bottom[i])  # group nodes placed within boundaries
						m.addConstr(y[nd] <= y_top[i])
						if g[nd].layer not in grp_layers:
							grp_layers[g[nd].layer] = []
						grp_layers[g[nd].layer].append(nd)
					# grpmax = max((len(lylist) for lylist in grp_layers.values()))
					# m.addConstr(y_top[i] - y_bottom[i] + 1 == grpmax)  # enforce close groups
					if self.fix_x_vars:
						grpmax = max((len(lylist) for lylist in grp_layers.values()))
						m.addConstr(y_top[i] - y_bottom[i] + 1 == grpmax)
						for lid in grp_layers:  # find the closest node above & below group for each layer
							gp_lay_yv = g[grp_layers[lid][0]].y
							closest_below, closest_below_nd = -1, -1
							closest_above, closest_above_nd = self.m_val, -1
							for nd_ot in g.layers[lid]:
								if nd_ot.id not in grp_layers[lid]:
									if gp_lay_yv > nd_ot.y > closest_below:
										closest_below, closest_below_nd = nd_ot.y, nd_ot.id
									elif gp_lay_yv < nd_ot.y < closest_above:
										closest_above, closest_above_nd = nd_ot.y, nd_ot.id
							if closest_below_nd != -1:
								m.addConstr(y_bottom[i] >= y[closest_below_nd] + 0.5)  # highest node below grp must be below lower bound, add buffer of 0.5
							if closest_above_nd != -1:
								m.addConstr(y_top[i] <= y[closest_above_nd] - 0.5)  # lowest node above grp must be above upper bound, add buffer of 0.5
					else:
						for lid in grp_layers:
							gp_lay_eg = grp_layers[lid][0]
							for nd in g.layers[lid]:
								if nd.id not in grp_layers[lid]:
									x12_r, u1, u2 = get_x_var_consts(xvars, gp_lay_eg, nd.id)
									m.addConstr(y[nd.id] - self.m_val * x12_r * x[u1, u2] - self.m_val * (1 - x12_r) // 2 + 1 <= y_bottom[i])
									m.addConstr(-y[nd.id] + self.m_val * x12_r * x[u1, u2] + self.m_val * (1 - x12_r) // 2 + 1 <= self.m_val - y_top[i])
			else:
				for group in sl_groups:
					lid = g[group[0]].layer
					for nd1 in group:
						for nd2 in group:
							if nd1 != nd2:
								for nd3 in (v for v in n_b_l[lid] if v not in group):
									x13_r, u131, u132 = get_x_var_consts(xvars, nd1, nd3)
									x23_r, u231, u232 = get_x_var_consts(xvars, nd2, nd3)
									m.addConstr(x13_r * x[u131, u132] - x23_r * x[u231, u232] + (1 - x13_r)//2 - (1 - x23_r)//2 == 0)  # Eq. 8
				ml_gp_to_layers = [set(g[nd].layer for nd in grp) for grp in ml_groups]
				ml_grps_layerids = [[] for _ in range(g.n_layers)]
				for i, group in enumerate(ml_groups):
					for lid in ml_gp_to_layers[i]:
						ml_grps_layerids[lid].append(i)

				for i, group in enumerate(ml_groups):
					for lid in sorted(list(ml_gp_to_layers[i]))[:-1]:
						ar_gnd = 0
						for i_nd in group:
							if g[i_nd].layer == lid:
								ar_gnd = i_nd
								break
						ar_gnd_2 = 0
						for i_nd in group:
							if g[i_nd].layer == lid + 1:
								ar_gnd_2 = i_nd
								break
						ot_xsum = gp.LinExpr()
						for nd in n_b_l[lid]:
							if nd not in group:
								x_r, u1, u2 = get_x_var_consts(xvars, ar_gnd, nd)
								ot_xsum += x_r * x[u1, u2] + (1 - x_r)//2
						ot_l2_xsum = gp.LinExpr()
						for nd in n_b_l[lid + 1]:
							if nd not in group:
								x_r, u1, u2 = get_x_var_consts(xvars, ar_gnd_2, nd)
								ot_l2_xsum += x_r * x[u1, u2] + (1 - x_r)//2
						m.addConstr(ot_xsum == ot_l2_xsum)

				for i, group in enumerate(ml_groups):
					# if not self.vertical_transitivity:
					# 	raise Exception("Need vertical transitivity for multi-layer group constraints")
					for nd1 in group:
						lid1 = g[nd1].layer
						for nd2 in group:
							lid2 = g[nd2].layer
							if lid1 == lid2 and nd1 != nd2:
								for nd3 in (v for v in n_b_l[lid1] if v not in group):
									x13_r, u131, u132 = get_x_var_consts(xvars, nd1, nd3)
									x23_r, u231, u232 = get_x_var_consts(xvars, nd2, nd3)
									m.addConstr(x13_r * x[u131, u132] - x23_r * x[u231, u232] + (1 - x13_r) // 2 - (1 - x23_r) // 2 == 0)  # Eq. 8
				# 		m.addConstr(y[nd1] >= lbx[i])  # Eqs. 10
				# 		m.addConstr(y[nd1] <= ubx[i])
				# 	for lid in ml_gp_to_layers[i]:
				# 		ar_gnd = 0
				# 		for i_nd in group:
				# 			if g[i_nd].layer == lid:
				# 				ar_gnd = i_nd
				# 				break
				# 		for nd_ot in n_b_l[lid]:
				# 			if nd_ot not in group:
				# 				x12_r, u1, u2 = get_x_var_consts(xvars, ar_gnd, nd_ot)
				# 				m.addConstr(y[ar_gnd] - x12_r * self.m_val * x[u1, u2] - self.m_val * (1 - x12_r) // 2 <= lbx[i])  # Eq. 11
				# 				m.addConstr(- y[ar_gnd] + x12_r * self.m_val * x[u1, u2] + self.m_val * (1 - x12_r) // 2 <= self.m_val - ubx[i])  # Eq. 12
				# 		for grp_ot in ml_grps_layerids[lid]:
				# 			if grp_ot != i and (grp_ot, i) not in ml_grp_pairs_seen:
				# 				ar_gnd_ot = 0
				# 				ml_grp_pairs_seen.append((i, grp_ot))
				# 				for i_nd in ml_groups[grp_ot]:
				# 					if g[i_nd].layer == lid:
				# 						ar_gnd_ot = i_nd
				# 						break
				# 				x21_r, u2, u1 = get_x_var_consts(xvars, ar_gnd_ot, ar_gnd)
				# 				print(i, ar_gnd, grp_ot, ar_gnd_ot, x21_r, u2, u1)
				# 				m.addConstr(ubx[grp_ot] - x21_r * self.m_val * x[u2, u1] - lbx[i] - self.m_val * (1 - x21_r) // 2 <= 0)  # Eqs. 13
				# 				m.addConstr(-ubx[grp_ot] + x21_r * self.m_val * x[u2, u1] + lbx[i] + self.m_val * (1 - x21_r) // 2 <= self.m_val)

	def __add_symmetry_maximization_constraints(self, m: gp.Model, y, ysum, ys_vars, sym_n, sym_e, sym_e_vars):
		if self.symmetry_maximization or any(ele[0] == "node_symmetry" or ele[0] == "node+edge_symmetry" for ele in self.hybrid_constraints):
			for ysv in ys_vars:
				m.addConstr(ysum[ysv] == y[ysv[0]] + y[ysv[1]])
				m.addConstr(self.m_val * sym_n[ysv] + ysum[ysv] - self.m_val >= 0)
				m.addConstr(self.m_val * sym_n[ysv] + self.m_val - ysum[ysv] >= 0)
			if self.symmetry_maximization_edges or any(ele[0] == "node+edge_symmetry" for ele in self.hybrid_constraints):
				ys_set = set(ys_vars)
				for e_v in sym_e_vars:
					sv1 = (e_v[0][0], e_v[1][0])
					if sv1 not in ys_set:
						sv1 = (e_v[1][0], e_v[0][0])
					sv2 = (e_v[0][1], e_v[1][1])
					if sv2 not in ys_set:
						sv2 = (e_v[1][1], e_v[0][1])
					m.addConstr(2 * sym_e[e_v] - sym_n[sv1] - sym_n[sv2] >= 0)

	def __add_angular_resolution_constraints(self, m: gp.Model, m_vars, m_v, y, combo_vars, combo, final):
		if self.angular_resolution or any(ele[0] == "angular_resolution" for ele in self.hybrid_constraints):
			for m_var in m_vars:
				m.addConstr(m_v[m_var] == (y[m_var[0]] - y[m_var[1]]) / 2)
			for c_var in combo_vars:
				if combo is None:
					m.addConstr(m_v[c_var[0]] * m_v[c_var[1]] + 1 <= final[c_var])
					m.addConstr(- m_v[c_var[0]] * m_v[c_var[1]] - 1 <= final[c_var])
				else:
					m.addConstr(combo[c_var] == m_v[c_var[0]] * m_v[c_var[1]] + 1)
					m.addGenConstrAbs(final[c_var], combo[c_var])

	def __add_hybrid_constraints(self, m: gp.Model, c, c_vars, c_consts, b, alpha, big_c, fair_var, n_sym, e_sym, ang, ang_vars):
		if self.hybrid_constraints:
			for metric, bound in self.hybrid_constraints:
				if metric == "crossings":
					csum = gp.LinExpr()
					for i, c_var in enumerate(c_vars):
						csum += self.gamma_1 * c_consts[i] * c[c_var]
					m.addConstr(csum <= int(bound))
				elif metric == "bends":
					m.addConstr(self.gamma_2 * b.sum() <= int(bound))
				elif metric == "edge_bundles":
					m.addConstr(self.gamma_bundle * alpha.sum() <= int(bound))
				elif metric == "min_max_crossings":
					m.addConstr(big_c <= int(bound))
				elif metric == "crossing_fairness" or metric == "bend_fairness":
					m.addConstr(fair_var <= int(bound))
				elif metric == "node_symmetry":
					m.addConstr(n_sym.sum() <= int(bound))
				elif metric == "edge_symmetry":
					m.addConstr(e_sym.sum() <= int(bound))
				elif metric == "angular_resolution":
					m.addConstr(sum((ang[cv] * c[cv] for cv in ang_vars)) <= int(bound))
				else:
					raise Exception(f"No metric with name {metric}.\nAllowed metrics: [crossings, bends, edge_bundles, min_max_crossings, crossing_fairness, bend_fairness, node_symmetry, node+edge_symmetry, angular_resolution]")

	def __we_need_c_vars(self):
		if self.crossing_minimization or (self.fairness_constraints and self.fairness_metric == "crossings") or self.min_max_crossings or self.angular_resolution or any(ele[0] == "crossings" or ele[0] == "min_max_crossings" or ele[0] == "angular_resolution" or ele[0] == "crossing_fairness" for ele in self.hybrid_constraints):
			return True
		return False

	def __we_need_y_vars(self):
		if self.vertical_transitivity or (self.fairness_constraints and self.fairness_metric == "bends") or self.bendiness_reduction or self.streamline or self.stratisfimal_y_vars or self.symmetry_maximization or self.angular_resolution or self.apply_node_weight_spacing or self.node_emphasis or (self.grouping_constraints and self.y_based_group_constraints) or any(ele[0] == "bends" or ele[0] == "edge_bundles" or ele[0] == "angular_resolution" or ele[0] == "bend_fairness" or ele[0] == "node_symmetry" or ele[0] == "node+edge_symmetry" for ele in self.hybrid_constraints):
			return True
		return False

	def __we_need_b_vars(self):
		if self.bendiness_reduction or (self.fairness_constraints and self.fairness_metric == "bends") or any(ele[0] == "bends" or ele[0] == "bend_fairness" for ele in self.hybrid_constraints):
			return True
		return False

	def __create_optimization_model(self, m: gp.Model, g: LayeredGraph, fix_x_vars=None, start_x_vars=None, groups=None):
		if not self.direct_transitivity and not self.vertical_transitivity:
			self.direct_transitivity = True
		if self.__we_need_y_vars():
			self.vertical_transitivity, self.direct_transitivity = True, False

		nodes_by_layer = g.get_ids_by_layer()
		edges_by_layer = g.get_edge_ids_by_layer()

		""" Add all variables """
		x_vars = []
		z_vars = []
		for i, name_list in nodes_by_layer.items():
			x_vars += list(itertools.combinations(name_list, 2))
			if self.stratisfimal_y_vars:
				z_vars += list(itertools.permutations(name_list, 2))
		x = m.addVars(x_vars, vtype=GRB.BINARY, name="x")
		z = None
		if self.stratisfimal_y_vars:
			z = m.addVars(z_vars, vtype=GRB.CONTINUOUS, lb=0, ub=self.m_val, name="z")
		c_vars, c_consts, c = None, None, None
		if self.__we_need_c_vars():
			c_vars, c_consts = reductions.normal_c_vars(g, edges_by_layer, self.mirror_vars, use_e_weights=self.apply_edge_weight)
			c = m.addVars(c_vars, vtype=GRB.BINARY, name="c")
		# if self.grouping_constraints:  # THIS IS FOR Y-VALUE BASED GROUP CONSTRAINTS
		# 	sl_groups, ml_groups = groups[0], groups[1]
		# 	grp_lb, grp_ub = [], []
		# 	if ml_groups:
		# 		grp_vars = list(range(len(ml_groups)))
		# 		grp_lb = m.addVars(grp_vars, vtype=GRB.CONTINUOUS, lb=0, ub=self.m_val, name="lb")
		# 		grp_ub = m.addVars(grp_vars, vtype=GRB.CONTINUOUS, lb=0, ub=self.m_val, name="ub")
		# 	if len(ml_groups) > 0 and not self.vertical_transitivity:
		# 		print("There are multilayer groups—swapping to vertical transitivity.")
		# 		self.vertical_transitivity, self.direct_transitivity = True, False
		y_t, y_b = None, None
		if self.grouping_constraints and self.y_based_group_constraints:
			grp_vars = list(range(len(groups[0]) + len(groups[1])))
			y_t = m.addVars(grp_vars, vtype=GRB.CONTINUOUS, lb=0, ub=self.m_val, name="y_t")
			y_b = m.addVars(grp_vars, vtype=GRB.CONTINUOUS, lb=0, ub=self.m_val, name="y_b")
		y = None
		if self.__we_need_y_vars():
			y = m.addVars([n.id for n in g], vtype=GRB.CONTINUOUS, lb=0, ub=self.m_val, name="y")
		b, b_vars = None, None
		if self.__we_need_b_vars():
			b_vars = list(g.edge_ids.keys())
			b = m.addVars(b_vars, vtype=GRB.CONTINUOUS, lb=0, ub=self.m_val, name="b")
		alpha, alpha_vars = None, None
		if self.edge_bundling or any(ele[0] == "edge_bundles" for ele in self.hybrid_constraints):
			alpha_vars = [(a1, a2) for lid in g.layers for ix, a1 in enumerate(nodes_by_layer[lid]) if g[a1].is_anchor_node for a2 in nodes_by_layer[lid][ix + 1:] if g[a2].is_anchor_node]
			alpha = m.addVars(alpha_vars, vtype=GRB.BINARY, name="alpha")
			# 	alpha = m.addVars(alpha_vars, vtype=GRB.INTEGER, lb=0, ub=self.m_val, name="alpha")
			# 	a_aux_diff = m.addVars(alpha_vars, vtype=GRB.INTEGER, lb=-self.m_val, ub=self.m_val, name="a_aux_diff")
			# 	bundle = m.addVar(vtype=GRB.INTEGER, lb=0, name="bundle")
		e, e_vars = None, None
		if self.node_emphasis:
			require_graph_props(g, require_node_data=["emphasis"])
			e_vars = list(set([(nd, nd_adj) if g[nd].layer < g[nd_adj].layer else (nd_adj, nd) for nd in g.node_data["emphasis"] for nd_adj in g.get_adj_list()[nd]]))
			e = m.addVars(e_vars, vtype=GRB.CONTINUOUS, lb=0, ub=self.m_val, name="e")
		ste = None
		if self.constrain_straight_long_arcs:
			ste_vars = [(le[i], le[i+1]) for le in g.get_long_edges() for i in range(len(le) - 1)]
			ste = m.addVars(ste_vars, vtype=GRB.BINARY, name="st_edge")
		cp, cp_vars, big_c = None, None, None
		if self.min_max_crossings or any(ele[0] == "min_max_crossings" for ele in self.hybrid_constraints):
			cp_vars = [(e.n1.id, e.n2.id) for e in g.edges]
			cp = m.addVars(cp_vars, vtype=GRB.INTEGER, lb=0, name="cp")
			big_c = m.addVar(lb=0, vtype=GRB.INTEGER, name="C")
		fair_var, b_aux = None, None
		if self.fairness_constraints or any(ele[0] == "crossing_fairness" or ele[0] == "bend_fairness" for ele in self.hybrid_constraints):
			fair_var = m.addVar(lb=0, vtype=GRB.CONTINUOUS, name="fair")
			if self.fairness_metric == "bends" or any(ele[0] == "bend_fairness" for ele in self.hybrid_constraints):
				b_aux = m.addVars(b_vars, vtype=GRB.CONTINUOUS, lb=-self.m_val, ub=self.m_val, name="b_aux")
		ysum_vars, ysum, sym_n, sym_e_vars, sym_e = None, None, None, None, None
		if self.symmetry_maximization or any(ele[0] == "node_symmetry" or ele[0] == "node+edge_symmetry" for ele in self.hybrid_constraints):
			ysum_vars = [v for nlist in nodes_by_layer.values() for v in itertools.combinations(nlist, 2)] + [(v, v) for v in g.node_ids]
			ysum = m.addVars(ysum_vars, vtype=GRB.CONTINUOUS, lb=0, ub=2*self.m_val, name="ysum")
			sym_n = m.addVars(ysum_vars, vtype=GRB.BINARY, name="sym_n")
			if self.symmetry_maximization_edges or any(ele[0] == "node+edge_symmetry" for ele in self.hybrid_constraints):
				sym_e_vars = [v for elist in edges_by_layer.values() for v in itertools.combinations(elist, 2)]
				sym_e = m.addVars(sym_e_vars, vtype=GRB.BINARY, name="sym_e")
		ang_m_vars, ang_m, ang_combo_vars, ang_combo, ang_final = None, None, None, None, None
		if self.angular_resolution or any(ele[0] == "angular_resolution" for ele in self.hybrid_constraints):
			if self.fix_x_vars:
				ang_combo_vars = g.edge_crossing_edges()
				ang_m_vars = list(set([v[0] for v in ang_combo_vars] + [v[1] for v in ang_combo_vars]))
				# mv_e_b_l = [[mv for mv in ang_m_vars if g[mv[0]].layer == lid] for lid in range(g.n_layers)]
			else:
				ang_m_vars = list(g.edge_ids.keys())
				ang_combo_vars = c_vars
			ang_m = m.addVars(ang_m_vars, vtype=GRB.CONTINUOUS, lb=-self.m_val, ub=self.m_val)
			ang_final = m.addVars(ang_combo_vars, vtype=GRB.CONTINUOUS, lb=0, name="ang")
			if any(ele[0] == "angular_resolution" for ele in self.hybrid_constraints):
				ang_combo = m.addVars(ang_combo_vars, vtype=GRB.CONTINUOUS)
			m.setParam("NonConvex", 2)

		m.update()  # required after adding variables in order to use them in constraints

		""" Fix variables/set starting assignments """
		if self.start_xy_vars:
			for v in m.getVars():
				if v.varName[:2] == "x[":
					xv1 = int(v.varName[2:v.varName.index(',')])
					xv2 = int(v.varName[v.varName.index(',') + 1:v.varName.index(']')])
					v.Start = get_x_var(self.x_var_assign, xv1, xv2)
				elif v.varName[:2] == "y[":
					v.Start = g[int(v.varName[2:v.varName.index(']')])].y
		if self.fix_x_vars:
			if any(v == 2 for v in self.x_var_assign.values()):  # use current rel positions if no prior optimization
				self.__assign_x_given_y()
			for k, v in self.x_var_assign.items():
				if v == 0 or v == 1:
					a1 = m.getVarByName(f"x[{k[0]},{k[1]}]")
					if a1:
						a1.lb, a1.ub = v, v
					else:
						a2 = m.getVarByName(f"x[{k[1]},{k[0]}]")
						a2.lb, a2.ub = 1 - v, 1 - v
		if fix_x_vars:
			for k, v in fix_x_vars.items():
				a1 = m.getVarByName(f"x[{k[0]},{k[1]}]")
				if a1:
					a1.lb, a1.ub = v, v
				else:
					a2 = m.getVarByName(f"x[{k[1]},{k[0]}]")
					a2.lb, a2.ub = 1 - v, 1 - v
		if start_x_vars:
			for v in m.getVars():
				v.Start = start_x_vars[v.varName]

		""" Set higher priority for branching on x-vars """
		if self.xvar_branch_priority:
			for v in m.getVars():
				if v.varName[:1] == "x":
					v.BranchPriority = 1
				else:
					v.BranchPriority = 0

		""" Butterfly reduction """
		butterfly_c_pairs = self.get_butterfly_cvars(g, c_vars)

		""" Node emphasis constraints """
		e_consts = self.__emphasis_constraints(m, g, self.x_var_assign, x, e_vars, e, y, b, nodes_by_layer)

		""" Set model objective function """
		opt_func = self.__optimization_function(g, c, c_vars, None, b, b_vars, c_consts, None, alpha, e, e_vars, e_consts, big_c, fair_var, sym_n, sym_e, ang_final, ang_combo_vars)
		m.setObjective(opt_func, GRB.MINIMIZE)

		""" Transitivity constraints """
		self.__transitivity(m, g, nodes_by_layer, self.x_var_assign, x, y, z)

		""" Edge crossing constraints """
		xv_use = self.__edge_crossings(m, g, c_vars, x, c, track_x_var_usage=self.symmetry_breaking, butterflies=butterfly_c_pairs)

		""" Bendiness reduction constraints """
		self.__bendiness_constraints(m, y, b_vars, b, b_aux)

		""" Break symmetry by fixing key x-var """
		self.__symmetry_breaking(m, xv_use, x_vars)

		""" Long-edge constraints """
		self.__long_edge_constraints(m, g, ste, self.x_var_assign, x, y, nodes_by_layer)

		""" Group constraints """
		self.__add_group_constraints(m, g, x_vars, x, y, y_t, y_b, groups, nodes_by_layer)

		""" Fix nodes at top/bottom/middle """
		self.__fix_specific_nodes(m, g, x, nodes_by_layer)

		""" Edge bundling constraints """
		self.__edge_bundling(m, g, alpha, alpha_vars, x, x_vars, y, nodes_by_layer)

		""" Min-max edge crossing constraints """
		self.__min_max_crossing_constraints(m, g, c_vars, c, cp_vars, cp, big_c, edges_by_layer)

		""" Fairness constraints """
		self.__fairness_constraints(m, g, c_vars, c, b_vars, b, fair_var)

		""" Symmetry (aesthetic metric) constraints """
		self.__add_symmetry_maximization_constraints(m, y, ysum, ysum_vars, sym_n, sym_e, sym_e_vars)

		""" Angular resolution constraints, make crossings close to 90 degrees """
		self.__add_angular_resolution_constraints(m, ang_m_vars, ang_m, y, ang_combo_vars, ang_combo, ang_final)

		""" Hybrid model bounding constraints """
		self.__add_hybrid_constraints(m, c, c_vars, c_consts, b, alpha, big_c, fair_var, sym_n, sym_e, ang_final, ang_combo_vars)

		return x_vars, c_vars

	def __optimize_with_subgraph_reduction(self, n_partitions, cluster, top_level_g, crosses, contacts, stack_to_nodeset, node_to_stack):  # DEPRECATED
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
			opt_val, x_vars_opt = self.__optimize_layout_standard(graph_arg=g_prime, is_subgraph=True, fix_x_vars=vars_to_fix)
			if unconstrained_opt_val != opt_val:
				self.print_info.append(f"\tSubgraph {i+1} has {opt_val} crossings but could have as few as {unconstrained_opt_val}")
			else:
				self.print_info.append(f"\tSubgraph {i+1} is optimal")
			unconstrained_opt_vals.append(unconstrained_opt_val)
			opt_vals.append(opt_val)
			for x_var, val in x_vars_opt.items():
				set_x_var(self.x_var_assign, x_var[0], x_var[1], val)

			self.__sequential_br(graph_arg=g_prime, substitute_x_vars=x_vars_opt)
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

	# def __optimize_full_subgraph_algorithm(self):  # DEPRECATED
	# 	t_allotted = self.cutoff_time if self.cutoff_time > 0 else 60
	# 	t_start = time.time()
	# 	n_partitions = 2
	# 	best_so_far = self.x_var_assign.copy()
	# 	best_n_cr = self.g.num_edge_crossings()
	# 	lower_bound = 0
	# 	# while lower_bound != best_n_cr and time.time() - t_start < t_allotted:
	# 	for i in range(3):
	# 		# n_rep = max(6 - n_partitions, 1)
	# 		n_rep = 1
	# 		cur_rep = 1
	# 		while cur_rep <= n_rep and lower_bound != best_n_cr and time.time() - t_start < t_allotted:
	# 			print(f"{n_partitions} partitions attempt {cur_rep}, {time.time()-t_start} elapased.")
	# 			cluster = list(SpectralClustering(n_clusters=n_partitions, assign_labels="discretize", affinity="precomputed").fit(self.g.adjacency_matrix()).labels_)
	# 			cluster.insert(0, 0)
	# 			top_level_g, crosses, contacts, stack_to_nodeset, node_to_stack = self.g.stacked_graph_from_subgraph_nodes(cluster)
	# 			n_cr, o_v, uo_v, top_ov = self.__optimize_with_subgraph_reduction(n_partitions, cluster, top_level_g, crosses, contacts, stack_to_nodeset, node_to_stack)
	# 			if n_cr < best_n_cr:
	# 				best_so_far = self.x_var_assign
	# 				best_n_cr = n_cr
	# 				print(f"New best, {n_cr} crossings")
	# 			new_lower_bound = self.__bound_on_optimality(cluster, top_level_g, n_cr, o_v, uo_v, top_ov)
	# 			if new_lower_bound > lower_bound:
	# 				lower_bound = new_lower_bound
	# 				print(f"New lower bound on optimal #crossings: {lower_bound}")
	# 			cur_rep += 1
	# 		n_partitions += 1
	# 	self.__sequential_br(substitute_x_vars=best_so_far)
	# 	if lower_bound == best_n_cr:
	# 		# we're optimal, baby
	# 		print("wooooo")
	# 	else:
	# 		print("ran out of time ;(")
	# 	return 0

	# def __optimize_locally_optimal(self):
	# 	# step 1: calculate num partitions by increasing until size falls within predicted bound (self.cutoff_time)
	# 	# step 2: split up partitions into layeredgraph objects (add hanger nodes)
	# 	do_bendiness_reduction, self.bendiness_reduction = self.bendiness_reduction, False
	# 	do_draw_graph, self.draw_graph = self.draw_graph, False
	# 	cutoff_partitions, n_partitions = 100, 1
	# 	cg = CollapsedGraph(self.g)
	# 	if self.n_partitions == -1:
	# 		while n_partitions <= cutoff_partitions:
	# 			cluster = list(SpectralClustering(n_clusters=n_partitions, assign_labels="discretize", affinity="precomputed").fit(self.g.adjacency_matrix()).labels_)
	# 			cg.subgraphs = cluster
	# 			lgs = cg.create_layered_graphs_from_subgraphs_dangling_nodes()
	# 			if max((lg.optimization_time_estimate() for lg in lgs)) < self.cutoff_time:
	# 				print(f"Optimizing with {n_partitions} partition{'s' if n_partitions > 1 else ''}")
	# 				print("Runtime estimates:", [lg.optimization_time_estimate() for lg in lgs])
	# 				break
	# 	else:
	# 		cluster = list(SpectralClustering(n_clusters=self.n_partitions, assign_labels="discretize", affinity="precomputed").fit(self.g.adjacency_matrix()).labels_)
	# 		cg.subgraphs = cluster
	# 		lgs = cg.create_layered_graphs_from_subgraphs_dangling_nodes()
	# 	cg.create_collapsed_graph_skeleton()
	# 	# print([[idx for idx, v in enumerate(cg.subgraphs) if v == i] for i in range(max(cg.subgraphs)+1)])
	#
	# 	if cg.optimization_time_estimate() < self.cutoff_time / 2:
	# 		# note: need separate opt to preserve x_var_assign
	# 		skele_opt = LayeredOptimizer(cg)
	# 		skele_opt.set_reasonable_params()
	# 		skele_opt.draw_graph = True
	# 		skele_opt.name = "collapsed_graph"
	# 		skele_opt.optimize_layout()
	# 		y_vs = skele_opt.__find_y_assignment_given_x(skele_opt.x_var_assign)
	# 		for nd, y_v in y_vs.items():
	# 			cg[nd].y = y_v
	# 	else:
	# 		# call bary-sift
	# 		improved_sifting(cg)
	#
	# 	for cr_edge in cg.crossing_edges:  # fix nodes
	# 		#  cn11----cn12    <- subgraph 1
	# 		#     `----.       <- crossing edge (cn11, cn22)
	# 		#  cn21----cn22    <- subgraph 2
	# 		cn11y = cg[cg.node_to_stack_node[cr_edge[0]]].y
	# 		cn12 = cg.get_collapsed_node(self.g[cr_edge[1]].layer, cg.subgraphs[cr_edge[0]])
	# 		cn22y = cg[cg.node_to_stack_node[cr_edge[1]]].y
	# 		cn21 = cg.get_collapsed_node(self.g[cr_edge[0]].layer, cg.subgraphs[cr_edge[1]])
	# 		# print(cr_edge, cn12, cn21)
	# 		if cn12 != -1 and cr_edge[1] in lgs[cg.subgraphs[cr_edge[0]]]:
	# 			# fix node corresponding to cn22 in subg1
	# 			# issue: cr_edge[0] may not be a node in lg2
	# 			# solution: since iteration is over all cr_edges, can just skip if not there
	# 			lgs[cg.subgraphs[cr_edge[0]]][cr_edge[1]].fix = cn22y - cn12.y
	# 			# print(lgs[cg.subgraphs[cr_edge[0]]][cr_edge[1]], lgs[cg.subgraphs[cr_edge[0]]][cr_edge[1]].fix)
	# 		elif cn12 == -1 and cr_edge[1] in lgs[cg.subgraphs[cr_edge[0]]]:
	# 			lgs[cg.subgraphs[cr_edge[0]]][cr_edge[1]].fix = cn22y
	# 			# print(lgs[cg.subgraphs[cr_edge[0]]][cr_edge[1]], lgs[cg.subgraphs[cr_edge[0]]][cr_edge[1]].fix)
	# 		if cn21 != -1 and cr_edge[0] in lgs[cg.subgraphs[cr_edge[1]]]:
	# 			# fix node corresponding to cn11 in subg2
	# 			lgs[cg.subgraphs[cr_edge[1]]][cr_edge[0]].fix = cn11y - cn21.y
	# 			# print(lgs[cg.subgraphs[cr_edge[1]]][cr_edge[0]], lgs[cg.subgraphs[cr_edge[1]]][cr_edge[0]].fix)
	# 		elif cn21 == -1 and cr_edge[0] in lgs[cg.subgraphs[cr_edge[1]]]:
	# 			lgs[cg.subgraphs[cr_edge[1]]][cr_edge[0]].fix = cn11y
	# 			# print(lgs[cg.subgraphs[cr_edge[1]]][cr_edge[0]], lgs[cg.subgraphs[cr_edge[1]]][cr_edge[0]].fix)
	#
	# 	# step 3: multiprocessing pool, optimize normally all layeredgraphs
	# 	manager = mp.Manager()
	# 	lgs_out = manager.list()
	# 	jobs = []
	# 	for lg in lgs:
	# 		p = mp.Process(target=self.optimize_target, args=(lg, lgs_out))
	# 		jobs.append(p)
	# 		p.start()
	# 	for proc in jobs:
	# 		proc.join()
	#
	# 	# step 4: extract node positions from optimized subgraphs
	# 	print("before", self.g.num_edge_crossings())
	# 	max_height = max((v.y for lg in lgs_out for v in lg.nodes))
	# 	for lg in lgs_out:
	# 		for nd in lg.nodes:
	# 			if cg.subgraphs[nd.id] == lg.subg_id:
	# 				self.g[nd.id].y = cg[cg.node_to_stack_node[nd.id]].y * max_height + nd.y
	# 	# self.g.check_position_validity()
	# 	self.__assign_x_given_y()
	# 	print("after", self.g.num_edge_crossings())
	#
	# 	for i in range(len(lgs_out)):
	# 		vis.draw_graph(lgs_out[i], f"example{i + 1}", groups=[0 if nd.fix == 0 else 1 for nd in lgs_out[i].nodes])
	#
	# 	# step 5: apply neighborhood_sift
	#
	# 	if do_bendiness_reduction:
	# 		self.__sequential_br()
	# 	if do_draw_graph:
	# 		vis.draw_graph(self.g, self.name, groups=cg.subgraphs)
	#
	# 	return self.g.num_edge_crossings()

	def __fix_x_var(self, m: gp.Model, k, v, y_var=False, c_var=False):
		if y_var:
			a1 = m.getVarByName(f"y[{k}]")
			if a1:
				a1.lb, a1.ub = v, v
		elif c_var:
			a1 = m.getVarByName(f"c[({k[0][0]},{k[0][1]}),({k[1][0]},{k[1][1]})]")
			if a1:
				a1.lb, a1.ub = v, v
		else:
			a1 = m.getVarByName(f"x[{k[0]},{k[1]}]")
			if a1:
				a1.lb, a1.ub = v, v
			else:
				a2 = m.getVarByName(f"x[{k[1]},{k[0]}]")
				a2.lb, a2.ub = 1 - v, 1 - v

	def __unfix_x_var(self, m: gp.Model, k, v, y_var=False, c_var=False):
		if y_var:
			a1 = m.getVarByName(f"y[{k}]")
			if a1:
				a1.lb, a1.ub = 0, self.m_val
				a1.Start = v
		elif c_var:
			a1 = m.getVarByName(f"c[({k[0][0]},{k[0][1]}),({k[1][0]},{k[1][1]})]")
			if a1:
				a1.lb, a1.ub = 0, self.m_val
				a1.Start = v
		else:
			a1 = m.getVarByName(f"x[{k[0]},{k[1]}]")
			if a1:
				a1.lb, a1.ub = 0, 1
				a1.Start = v
			else:
				a2 = m.getVarByName(f"x[{k[1]},{k[0]}]")
				a2.lb, a2.ub = 0, 1
				a2.Start = v

	def __fix_everything(self, m: gp.Model):
		for v in m.getVars():
			if v.varName[:1] == "x":
				xv1 = int(v.varName[v.varName.index('[')+1: v.varName.index(',')])
				xv2 = int(v.varName[v.varName.index(',')+1: v.varName.index(']')])
				self.__fix_x_var(m, (xv1, xv2), self.x_var_assign[xv1, xv2])
			elif v.varName[:1] == "c":
				cvsp = v.varName.translate({ord(i): None for i in 'c[]()'}).split(',')
				e1 = self.g.edge_ids[(int(cvsp[0]), int(cvsp[1]))]
				e2 = self.g.edge_ids[(int(cvsp[2]), int(cvsp[3]))]
				if (e1.n1.y > e2.n1.y and e1.n2.y < e2.n2.y) or (e1.n1.y < e2.n1.y and e1.n2.y > e2.n2.y):
					self.__fix_x_var(m, ((e1.n1.id, e1.n2.id), (e2.n1.id, e2.n2.id)), 1, c_var=True)
				else:
					self.__fix_x_var(m, ((e1.n1.id, e1.n2.id), (e2.n1.id, e2.n2.id)), 0, c_var=True)

	def __incremetal_opt(self, graph: LayeredGraph, subgraph: list, m: gp.Model, env, cv_set, nbhd_width=0):
		# TODO (future): make it re-fix all vars used after optimizing, and only unfix subgraph node vars
		# Fix everything except subgraph and optimize.
		# subgraph: idx=node id, val=True if in subgraph else False
		# If all nodes already positioned, trust self.x_var_assign, else add the subgraph. Need crossing edges though
		if len(graph.nodes) != len(self.g.nodes):  # add new nodes if necessary
			for nid, v in enumerate(subgraph):
				if v and nid not in graph:
					graph.add_node(self.g[nid].layer, idx=nid)
			for edge in self.g.edges:
				if (subgraph[edge.n1.id] == 1 or subgraph[edge.n2.id] == 1) and (edge.n1.id, edge.n2.id) not in graph.edge_ids and edge.n1.id in graph.node_ids and edge.n2.id in graph.node_ids:
					graph.add_edge(edge.n1.id, edge.n2.id)
		# print("before", len(self.x_var_assign))
		# for xv in list(self.x_var_assign.keys()):

		subgraph_nodes = [nid for nid, v in enumerate(subgraph) if v]
		adj_list = graph.get_adj_list()
		ebl = graph.get_edges_by_layer()
		done = set()

		# for xv, val in self.x_var_assign.items():
		# 	if nbhd_width == 0:
		# 		if subgraph[xv[0]] or subgraph[xv[1]]:  # delete x-var assignements
		# 			self.__unfix_x_var(m, xv, val)
		# 			# del self.x_var_assign[xv]
		# 		else:
		# 			self.__fix_x_var(m, xv, val)
		# 	else:
		# 		if (subgraph[xv[0]] or subgraph[xv[1]]) and graph[xv[0]].v_nbhd and graph[xv[1]].v_nbhd:
		# 			self.__unfix_x_var(m, xv, val)
		# 		else:
		# 			self.__fix_x_var(m, xv, val)

		# if self.vertical_transitivity:  # Fix y-vars (probably not necessary since fixed by xvars)
		# 	subg_layers = {nd.layer for nd in graph if subgraph[nd.id]}
		# 	for nid, nd in graph.node_ids.items():
		# 		if nd.layer in subg_layers:
		# 			self.__unfix_x_var(m, nid, nd.y, y_var=True)
		# 		else:
		# 			self.__fix_x_var(m, nid, nd.y, y_var=True)

		for nd in subgraph_nodes:  # fix all vars associated with
			for ot_nd in graph.layers[graph[nd].layer]:  # fix x-vars
				if nd != ot_nd.id and (nbhd_width == 0 or ot_nd.v_nbhd):
					if (nd, ot_nd.id) in self.x_var_assign:
						self.__unfix_x_var(m, (nd, ot_nd.id), self.x_var_assign[(nd, ot_nd.id)])
					else:
						self.__unfix_x_var(m, (ot_nd.id, nd), self.x_var_assign[(ot_nd.id, nd)])

			for nd_adj in adj_list[nd]:  # fix c-vars
				if nd_adj not in done:
					edge_l = graph[nd].layer if graph[nd].layer < graph[nd_adj].layer else graph[nd_adj].layer
					the_edge = (nd, nd_adj) if graph[nd].layer < graph[nd_adj].layer else (nd_adj, nd)
					for ed_adj in ebl[edge_l]:
						if (the_edge, ed_adj) in cv_set:
							if ed_adj[0] != the_edge[0] and ed_adj[1] != the_edge[1]:
								e1 = graph.edge_ids[the_edge]
								e2 = graph.edge_ids[ed_adj]
								if (e1.n1.y > e2.n1.y and e1.n2.y < e2.n2.y) or (e1.n1.y < e2.n1.y and e1.n2.y > e2.n2.y):
									self.__unfix_x_var(m, (the_edge, ed_adj), 1, c_var=True)
								else:
									self.__unfix_x_var(m, (the_edge, ed_adj), 0, c_var=True)
						elif (ed_adj, the_edge) in cv_set:
							if ed_adj[0] != the_edge[0] and ed_adj[1] != the_edge[1]:
								e1 = graph.edge_ids[the_edge]
								e2 = graph.edge_ids[ed_adj]
								if (e1.n1.y > e2.n1.y and e1.n2.y < e2.n2.y) or (e1.n1.y < e2.n1.y and e1.n2.y > e2.n2.y):
									self.__unfix_x_var(m, (ed_adj, the_edge), 1, c_var=True)
								else:
									self.__unfix_x_var(m, (ed_adj, the_edge), 0, c_var=True)
			done.add(nd)

		# print([nd for nd in graph if nd.layer in subg_layers])
		# print([xv for xv in self.x_var_assign if subgraph[xv[0]] or subgraph[xv[1]]])
		# print("after", len(self.x_var_assign))
		# return self.__optimize_layout_standard(graph_arg=graph, fix_x_vars=self.x_var_assign)
		ret_v = self.__optimize_crossing_reduction_model(m, graph, env)

		done.clear()
		for nd in subgraph_nodes:  # fix all vars associated with
			for ot_nd in graph.layers[graph[nd].layer]:  # fix x-vars
				if nd != ot_nd.id and (nbhd_width == 0 or ot_nd.v_nbhd):
					if (nd, ot_nd.id) in self.x_var_assign:
						self.__fix_x_var(m, (nd, ot_nd.id), self.x_var_assign[(nd, ot_nd.id)])
					else:
						self.__fix_x_var(m, (ot_nd.id, nd), self.x_var_assign[(ot_nd.id, nd)])

			for nd_adj in adj_list[nd]:  # fix c-vars
				if nd_adj not in done:
					edge_l = graph[nd].layer if graph[nd].layer < graph[nd_adj].layer else graph[nd_adj].layer
					the_edge = (nd, nd_adj) if graph[nd].layer < graph[nd_adj].layer else (nd_adj, nd)
					for ed_adj in ebl[edge_l]:
						if (the_edge, ed_adj) in cv_set:
							if ed_adj[0] != the_edge[0] and ed_adj[1] != the_edge[1]:
								e1 = graph.edge_ids[the_edge]
								e2 = graph.edge_ids[ed_adj]
								if (e1.n1.y > e2.n1.y and e1.n2.y < e2.n2.y) or (e1.n1.y < e2.n1.y and e1.n2.y > e2.n2.y):
									self.__fix_x_var(m, (the_edge, ed_adj), 1, c_var=True)
								else:
									self.__fix_x_var(m, (the_edge, ed_adj), 0, c_var=True)
						elif (ed_adj, the_edge) in cv_set:
							if ed_adj[0] != the_edge[0] and ed_adj[1] != the_edge[1]:
								e1 = graph.edge_ids[the_edge]
								e2 = graph.edge_ids[ed_adj]
								if (e1.n1.y > e2.n1.y and e1.n2.y < e2.n2.y) or (e1.n1.y < e2.n1.y and e1.n2.y > e2.n2.y):
									self.__fix_x_var(m, (ed_adj, the_edge), 1, c_var=True)
								else:
									self.__fix_x_var(m, (ed_adj, the_edge), 0, c_var=True)
			done.add(nd)
		return ret_v

	def local_opt_increment(self, bucket_size, neighborhood_fn=bfs_neighborhood, candidate_fn=degree_candidate, vertical_width=0, movement_data=False):
		# opt_g = LayeredGraph()
		g = self.g
		if self.sequential_bendiness:
			do_bendiness_reduction, self.bendiness_reduction = self.bendiness_reduction, False
		do_draw_graph, self.draw_graph = self.draw_graph, False
		collapse_leaves, self.collapse_leaves, self.collapse_subgraphs = self.collapse_leaves, False, False
		if self.create_video and not os.path.isdir(f"Images/{self.name}"):
			os.mkdir(f"Images/{self.name}")
		if collapse_leaves:
			g = g.collapse_leaves()
			self.x_var_assign = {x_v: 2 for n_l in g.get_ids_by_layer().values() for x_v in itertools.combinations(n_l, 2)}
		cr_counts = [g.num_edge_crossings()]
		times = [0]
		times_moved = [0] * g.n_nodes
		candidate_moved = 0
		iterations_with_movement = 0
		opt_time = 0
		with gp.Env() as env, gp.Model(env=env) as m:
			x_vs, c_vs = self.__crossing_reduction_model(m, g)
			cv_set = set(c_vs)
			# TODO: Fix everything to start
			self.__assign_x_given_y()
			self.__fix_everything(m)
			if self.cutoff_time > 0:
				candidate_fn(g, init=True)
				iter_ct = 0
				frame_count = 0
				t_since_last = time.time()
				st_time = time.time()
				iter_without_improvement = 0
				while self.cutoff_time > 0:
					if self.create_video:
						vis.draw_graph(g, f"{self.name}/frame_{frame_count}", as_png=True, label_nodes=False, groups=[0] * g.n_nodes, gravity=True, copies=2)
						frame_count += 2
					candidate = candidate_fn(g)
					next_partition = neighborhood_fn(g, candidate, bucket_size, nbhd_width=vertical_width)
					neighborhood = [nid for nid, v in enumerate(next_partition) if v]
					y_save = [g[nd].y for nd in neighborhood]
					print(candidate, neighborhood)
					if self.create_video:
						vis.draw_graph(g, f"{self.name}/frame_{frame_count}", as_png=True, emphasize_nodes=[True if ind == candidate else False for ind in range(g.n_nodes)], groups=next_partition, gravity=True, label_nodes=False, copies=3)
						frame_count += 3
					out = self.__incremetal_opt(g, next_partition, m, env, cv_set, nbhd_width=vertical_width)
					self.cutoff_time -= time.time() - t_since_last
					t_since_last = time.time()
					opt_time += out[1]
					cr_counts.append(out[0])
					movement = [g[neighborhood[i]].y - y_save[i] for i in range(len(neighborhood))]
					for i, v in enumerate(movement):
						if v != 0:
							times_moved[neighborhood[i]] += 1
						if neighborhood[i] == candidate and v != 0:
							candidate_moved += 1
					if not all((v == 0 for v in movement)):
						iterations_with_movement += 1
					if cr_counts[-1] == cr_counts[-2]:
						iter_without_improvement += 1
					else:
						iter_without_improvement = 0
					if self.create_video and not all((v == 0 for v in movement)):
						gps = next_partition.copy()
						for i, v in enumerate(neighborhood):
							if movement[i] != 0:
								gps[v] = 2
						edges_moved = {(v, vp) for i, v in enumerate(neighborhood) if movement[i] != 0 for vp in g.get_adj_list()[v]}
						vis.draw_graph(g, f"{self.name}/frame_{frame_count}", as_png=True, emphasize_nodes=[True if ind == candidate else False for ind in range(g.n_nodes)], groups=gps, emphasize_edges=edges_moved, gravity=True, label_nodes=False, copies=3)
						frame_count += 3
					penalty_fn(g, neighborhood, candidate, movement, iter_ct, no_repeats=True)
					iter_ct += 1
					if iter_without_improvement == 5:
						iter_without_improvement = 0
						# TODO increase all blinds, refresh blinds. Full algorithm only
					times.append(time.time() - st_time)
					print("Iteration:", iter_ct, "\tTime left:", self.cutoff_time, "\tCrossings:", out[0])
			if collapse_leaves:
				self.collapse_leaves = True
				self.__optimize_subgraphs(g, x_vs, 0)
		if self.create_video:
			import imageio.v3 as iio
			from numpy import stack
			from pygifsicle import optimize
			frames = stack([iio.imread(f"Images/{self.name}/frame_{i}.png") for i in range(frame_count)], axis=0)
			iio.imwrite(f'Images/{self.name}.gif', frames, format="gif", fps=5)
			optimize(f'Images/{self.name}.gif')
			# with imageio.get_writer(f'Images/{self.name}.gif', mode='I') as writer:
			# 	for i in range(iter_ct * 3):
			# 		print(i)
			# 		image = imageio.v3.imread(f"Images/{self.name}/frame_{i}.png")
			# 		writer.append_data(image)

		if do_draw_graph:
			print(times_moved)
			vis.draw_graph(g, "solution_neighborhood", color_scale=times_moved)

		if movement_data:
			return opt_time, cr_counts[-1], cr_counts, times, times_moved, candidate_moved, iterations_with_movement
		else:
			return opt_time, cr_counts[-1], cr_counts, times

	def optimize_target(self, graph: LayeredGraph, return_list):
		self.__optimize_layout_standard(graph_arg=graph)
		return_list.append(graph)

	def optimize_with_starting_assignments(self, assigments):
		out = self.__optimize_layout_standard(fix_x_vars=assigments)
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
				for nd, ct in counts.items():
					y_assign[nd] = ct
			else:
				return "non-transitive x-values"
		return y_assign

	def __assign_x_given_y(self):
		for k in self.x_var_assign:
			if self.g[k[0]].y > self.g[k[1]].y:
				self.x_var_assign[k] = 0
			else:
				self.x_var_assign[k] = 1

	def __assign_y_given_x(self):
		# print(self.x_var_assign)
		# print([nd.y for nd in self.g.nodes])
		for nd in self.g:
			nd.y = 0
		for x_var, val in self.x_var_assign.items():
			if x_var[0] in self.g.node_ids and x_var[1] in self.g.node_ids and val != 2:
				self.g[x_var[val]].y += 1
		# print([nd.y for nd in self.g.nodes])

	def optimize_layout(self, fix_xvars=None, local_opt=False, force_optimal=False, bucket_size=1000, pct=1, **kwargs):
		self.crossing_minimization = kwargs.get("crossing_minimization", False)
		self.bendiness_reduction = kwargs.get("bendiness_reduction", False)
		self.gamma_1 = kwargs.get("gamma_1", 1)
		self.gamma_2 = kwargs.get("gamma_2", 1)
		# self.m_val = kwargs.get("m_val", round(1.5 * max(len(lr) for lr in self.g.layers.values())))
		self.sequential_bendiness = kwargs.get("sequential_bendiness", False)
		self.local_opt = kwargs.get("local_opt", False)
		self.local_opt_heuristic = kwargs.get("local_opt_heuristic", "incremental")
		self.n_partitions = kwargs.get("n_partitions", -1)
		self.return_full_data = kwargs.get("return_full_data", False)
		self.cutoff_time = kwargs.get("cutoff_time", 0)
		self.do_subg_reduction = kwargs.get("do_subg_reduction", False)
		self.return_x_vars = kwargs.get("return_x_vars", False)
		self.butterfly_reduction = kwargs.get("butterfly_reduction", False)
		self.draw_graph = kwargs.get("draw_graph", False)
		self.symmetry_breaking = kwargs.get("symmetry_breaking", True)
		self.heuristic_start = kwargs.get("heuristic_start", False)
		self.aggro_presolve = kwargs.get("presolve", False)
		self.mip_relax = kwargs.get("mip_relax", False)
		self.xvar_branch_priority = kwargs.get("xvar_branch_priority", False)
		self.direct_transitivity = kwargs.get("direct_transitivity", False)
		self.vertical_transitivity = kwargs.get("vertical_transitivity", False)
		self.mirror_vars = kwargs.get("mirror_vars", False)
		self.stratisfimal_y_vars = kwargs.get("stratisfimal_y_vars", False)
		self.symmetry_constraints = kwargs.get("symmetry_constraints", True)
		self.cycle_constraints = kwargs.get("cycle_constraints", False)
		self.collapse_subgraphs = kwargs.get("collapse_subgraphs", False)
		self.collapse_leaves = kwargs.get("collapse_leaves", False)
		self.claw_constraints = kwargs.get("claw_constraints", False)
		self.dome_path_constraints = kwargs.get("dome_path_constraints", False)
		self.polyhedral_constraints = kwargs.get("polyhedral_constraints", False)
		self.grouping_constraints = kwargs.get("grouping_constraints", False)
		self.y_based_group_constraints = kwargs.get("y_based_group_constraints", False)
		self.node_emphasis = kwargs.get("node_emphasis", False)
		self.emphasis_cr_weight = kwargs.get("emphasis_cr_weight", 3)
		self.emphasis_br_weight = kwargs.get("emphasis_br_weight", 1)
		self.apply_node_weight_spacing = kwargs.get("apply_node_weight_spacing", False)
		self.apply_edge_weight = kwargs.get("apply_edge_weight", False)
		self.angular_resolution = kwargs.get("angular_resolution", False)
		self.symmetry_maximization = kwargs.get("symmetry_maximization", False)
		self.symmetry_maximization_edges = kwargs.get("symmetry_maximization_edges", False)
		self.min_max_crossings = kwargs.get("min_max_crossings", False)
		self.gamma_min_max = kwargs.get("gamma_min_max", 1)
		self.streamline = kwargs.get("streamline", False)
		self.anchor_proximity = kwargs.get("anchor_proximity", 0.3)
		self.fix_x_vars = kwargs.get("fix_x_vars", False)
		self.start_xy_vars = kwargs.get("start_xy_vars", False)
		self.fix_nodes = kwargs.get("fix_nodes", False)
		self.edge_bundling = kwargs.get("edge_bundling", False)
		self.gamma_bundle = kwargs.get("gamma_bundle", 1)
		self.fairness_constraints = kwargs.get("fairness_constraints", False)
		self.fairness_metric = kwargs.get("fairness_metric", "crossings")
		self.gamma_fair = kwargs.get("gamma_fair", 1)
		self.return_experiment_data = kwargs.get("return_experiment_data", False)
		self.create_video = kwargs.get("create_video", False)
		self.constrain_straight_long_arcs = kwargs.get("constrain_straight_long_arcs", False)
		self.long_arc_bend_limit = kwargs.get("long_arc_bend_limit", 0)
		self.record_solution_data_over_time = kwargs.get("record_solution_data_over_time", False)
		self.name = kwargs.get("name", "graph1")
		self.nthreads = kwargs.get("nthreads", 0)
		self.hybrid_constraints = kwargs.get("hybrid_constraints", [])
		if self.polyhedral_constraints:
			self.claw_constraints, self.dome_path_constraints = True, True

		if force_optimal:
			self.local_opt, local_opt = False, False

		if self.local_opt or local_opt:
			if self.local_opt_heuristic == "partition":
				out = self.__optimize_locally_optimal()
			elif self.local_opt_heuristic == "incremental":
				# out = self.__optimize_incremental_local()
				out = self.local_opt_increment(bucket_size)
			else:
				raise Exception("no heuristic of that name")
		elif fix_xvars is not None:
			out = self.__optimize_layout_standard(fix_x_vars=fix_xvars)
		else:
			out = self.__optimize_layout_standard()
		return out

	def just_bendiness_reduction(self, streamline=True):
		self.__assign_x_given_y()
		self.__sequential_br(streamline=streamline)
		vis.draw_graph(self.g, self.name)

	def set_reasonable_params(self):
		self.symmetry_breaking = True
		self.xvar_branch_priority = True
		self.direct_transitivity = True
