import itertools
import re
import gurobipy as gp
from gurobipy import GRB
from src import graph, vis, reductions
import time


def sequential_br(g, y_vars, x_var_opt, mv, subgraph_seq=False):
	n_constr = 0
	m2 = gp.Model()
	y = m2.addVars(y_vars, vtype=GRB.CONTINUOUS, lb=0, ub=mv, name="y")
	# y = m2.addVars(y_vars, vtype=GRB.INTEGER, lb=0, ub=mv, name="y")
	m2.update()
	for v in m2.getVars():
		v.start = g.node_names[int(v.varName[2:v.varName.index(']')])].y
	b_vars = list(g.edge_names.keys())
	b = m2.addVars(b_vars, vtype=GRB.CONTINUOUS, lb=0, ub=mv, name="b")
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
			# print(v.varName, v.x)
			g.node_names[int(v.varName[2:v.varName.index(']')])].y = float(v.x)
	# g.node_names[int(v.varName[2:v.varName.index(']')])].y = round(v.x)
	return n_constr


def edge_crossings(model: gp.Model, c_vars, g, x, c):
	n_constr_0, n_constr_1, n_constr_2 = 0, 0, 0
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
	return n_constr_0, n_constr_1, n_constr_2


def compute_variable_assignments(g, x_vars, c_vars):
	x_assignments = {}
	c_assignments = {}
	for x_var in x_vars:
		if g.node_names[x_var[0]].y < g.node_names[x_var[1]].y:
			x_assignments[x_var] = 1
		else:
			x_assignments[x_var] = 0
	for c_var in c_vars:
		sl1 = g.edge_names[c_var[0]].same_layer_edge
		sl2 = g.edge_names[c_var[1]].same_layer_edge
		if (c_var[0][0], c_var[1][0]) in x_assignments:
			x1 = x_assignments[c_var[0][0], c_var[1][0]]
		else:
			x1 = 1 - x_assignments[c_var[1][0], c_var[0][0]]
		if sl1 and sl2:
			if (c_var[0][1], c_var[1][1]) in x_assignments:
				x2 = x_assignments[c_var[0][1], c_var[1][1]]
			else:
				x2 = 1 - x_assignments[c_var[1][1], c_var[0][1]]
			if (c_var[1][0], c_var[0][1]) in x_assignments:
				x3 = x_assignments[c_var[1][0], c_var[0][1]]
			else:
				x3 = 1 - x_assignments[c_var[0][1], c_var[1][0]]
			if (c_var[0][0], c_var[1][1]) in x_assignments:
				x4 = x_assignments[c_var[0][0], c_var[1][1]]
			else:
				x4 = 1 - x_assignments[c_var[1][1], c_var[0][0]]
			c_assignments[c_var] = x1 * x2 * x3 + x4 * (1 - x2) * (1 - x3) + (1 - x3) * (1 - x1) * x4 + x2 * (1 - x4) * x1 + (1 - x1) * x4 * (1 - x2) + x3 * x2 * (1 - x4) + (1 - x4) * x1 * x3 + (1 - x2) * (1 - x3) * (1 - x1)
		elif sl1:
			if (c_var[0][1], c_var[1][0]) in x_assignments:
				x2 = x_assignments[c_var[0][1], c_var[1][0]]
			else:
				x2 = 1 - x_assignments[c_var[1][0], c_var[0][1]]
			c_assignments[c_var] = x1 * (1 - x2) + (1 - x1) * x2
		elif sl2:
			if (c_var[0][0], c_var[1][1]) in x_assignments:
				x2 = x_assignments[c_var[0][0], c_var[1][1]]
			else:
				x2 = 1 - x_assignments[c_var[1][1], c_var[0][0]]
			c_assignments[c_var] = x1 * (1 - x2) + (1 - x1) * x2
		else:
			if (c_var[0][1], c_var[1][1]) in x_assignments:
				x2 = x_assignments[c_var[0][1], c_var[1][1]]
			else:
				x2 = 1 - x_assignments[c_var[1][1], c_var[0][1]]
			c_assignments[c_var] = x1 * (1 - x2) + (1 - x1) * x2
	return x_assignments, c_assignments


def optimize_layout(g: graph.LayeredGraph, bendiness_reduction, gamma_1=1, gamma_2=1, sequential_bendiness=True,
					cutoff_time=0, do_subg_reduction=False, is_subgraph=False, contact_nodes=None, contact_sides=None):
	t1 = time.time()
	if do_subg_reduction:
		cut, cross = reductions.kargers_algo_cut_finder(g, 10)
		cutnode_layers = [g.node_names[n].layer for n in cut[0]]
		if len(cutnode_layers) == len(set(cutnode_layers)):
			do_subg_reduction = False
			print("Failed to find large enough cut set")
		else:
			g.stack_subgraph(set(cut[0]), [g.edge_names[ed[0], ed[1]] for ed in cross])

	nodes_by_layer = g.get_names_by_layer()
	edges_by_layer = g.get_edge_names_by_layer()
	n_constraints_generated = [0] * 6  # simple edge, hybrid edge, same layer edge, vertical pos, bendiness, total
	m_val = 50

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
	# z = m.addVars(z_vars, vtype=GRB.CONTINUOUS, lb=0, ub=m_val, name="z")
	c_vars, c_consts = reductions.normal_c_vars(g, edges_by_layer)
	c = m.addVars(c_vars, vtype=GRB.CONTINUOUS, name="c")
	y_vars = [n.name for n in g.nodes]
	y = m.addVars(y_vars, vtype=GRB.CONTINUOUS, lb=0, ub=m_val, name="y")
	m.update()

	""" Starting variable assignment """
	# g.barycentric_reordering(10)
	# x_assign, c_assign = compute_variable_assignments(g, x_vars, c_vars)
	# for var in m.getVars():
	# 	if var.varName[:1] == "y":
	# 		var.start = g.node_names[int(var.varName[2:var.varName.index(']')])].y
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

	""" Non-sequential optimization function """
	if not sequential_bendiness:  # TODO convert to LinExp optimization func
		# z_vars = []
		# for i, name_list in nodes_by_layer.items():
		# 	z_vars += list(itertools.permutations(name_list, 2))
		# z = m.addVars(z_vars, vtype=GRB.INTEGER, lb=0, ub=m_val, name="z")
		if bendiness_reduction:
			b_vars = list(g.edge_names.keys())
			b = m.addVars(b_vars, vtype=GRB.INTEGER, lb=0, ub=m_val, name="b")
			m.setObjective(gamma_1 * c.sum() + gamma_2 * b.sum(), GRB.MINIMIZE)
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
	n_cs = edge_crossings(m, c_vars, g, x, c)
	for i, val in enumerate(n_cs):
		n_constraints_generated[i] += val

	""" Vertical position, implication version """
	for x_var in x_vars:
		m.addGenConstrIndicator(x[x_var], True, y[x_var[0]] + 1 <= y[x_var[1]])
		m.addGenConstrIndicator(x[x_var], False, y[x_var[0]] >= 1 + y[x_var[1]])
		n_constraints_generated[3] += 2

	""" Original vertical position constraints """
	# for x_var in x_vars:
	# 	# print(x_var)
	# 	m.addConstr(z[x_var] - m_val * x[x_var] <= 0, f"1.1z{x_var}")
	# 	m.addConstr(z[x_var] - y[x_var[0]] - (m_val * x[x_var]) >= -1 * m_val, f"1.2z{x_var}")
	# 	m.addConstr(y[x_var[1]] - z[x_var] - x[x_var] >= 0, f"1.3z{x_var}")
	# 	m.addConstr(z[x_var] <= y[x_var[0]], f"1.4z{x_var}")
	# 	# m.addConstr(z[x_var] >= 0, f"1.5z{x_var}")
	# 	# print(f"z{x_var}-{m_val}*x{x_var} <= 0")
	# 	# print(f"y{x_var[1]} - z{x_var} - x{x_var} >= 0")
	#
	# 	m.addConstr(z[x_var[1], x_var[0]] - m_val * (1 - x[x_var]) <= 0, f"2.1z{x_var}")
	# 	m.addConstr(z[x_var[1], x_var[0]] - y[x_var[1]] - m_val * (1 - x[x_var]) >= -1 * m_val,
	# 				f"2.2z{x_var}")
	# 	m.addConstr(y[x_var[0]] - z[x_var[1], x_var[0]] - (1 - x[x_var]) >= 0, f"2.3z{x_var}")
	# 	m.addConstr(z[x_var[1], x_var[0]] <= y[x_var[1]], f"2.4z{x_var}")
	# 	# m.addConstr(z[x_var[1],x_var[0]] >= 0, f"2.5z{x_var}")
	# 	# print(f"z[{x_var[1],x_var[0]}]-{m_val}*(1-x{x_var}) <= 0")
	# 	# print(f"y{x_var[0]} - z[{x_var[1],x_var[0]}] - (1-x{x_var}) >= 0")
	#
	# 	# m.addConstr(k[x_var] + x[x_var] - 1 >= 0)  # SOS-constraint version. Use with full LP relaxation
	# 	# m.addSOS(GRB.SOS_TYPE1, [k[x_var], x[x_var]])
	#
	# 	n_constraints_generated[3] += 10

	""" Non-sequential bendiness reduction, original Stratisfimal verseion"""
	if not sequential_bendiness and bendiness_reduction:
		for b_var in b_vars:
			m.addConstr(y[b_var[0]] - y[b_var[1]] <= b[b_var], f"1bend{b_var}")
			m.addConstr(y[b_var[1]] - y[b_var[0]] <= b[b_var], f"2bend{b_var}")
			n_constraints_generated[4] += 2

	""" Optimize model """
	t1 = time.time() - t1
	print("Time to input constraints: ", t1)
	n_constraints_generated[5] = sum(n_constraints_generated[:5])
	print(n_constraints_generated)
	g.print_layer_counts()
	t2 = time.time()
	if cutoff_time > 0:
		m.setParam("TimeLimit", cutoff_time)
	m.setParam("OutputFlag", 0)
	m.optimize()
	print('Obj: %g' % m.objVal)
	t2 = time.time() - t2
	print("Time to optimize: ", t2)
	for v in m.getVars():
		if v.varName[:1] == 'y':
			g.node_names[int(v.varName[2:v.varName.index(']')])].y = float(v.x)

	""" Optimize subgraphs and re-insert """
	x_vars_opt = {}
	if do_subg_reduction:
		t3 = time.time()
		g_prime = graph.LayeredGraph()
		for node in g.stacked_nodes[0]:
			g_prime.add_node(node.name, node.layer, is_anchor=node.is_anchor_node)
		for edge in g.stacked_edges[0]:
			if edge.n1.name in g_prime.node_names and edge.n2.name in g_prime.node_names:
				g_prime.add_edge(edge.n1.name, edge.n2.name)
		optimize_layout(g_prime, True, is_subgraph=True)

		vis.draw(g, "interim")
		vis.draw(g_prime, "interim_subg")

		for v in m.getVars():
			if v.varName[:1] == "x":
				x_vars_opt[int(v.varName[2:v.varName.index(',')]), int(v.varName[v.varName.index(',') + 1:v.varName.index(']')])] = round(v.x)

		gp_layers = g_prime.get_names_by_layer()
		for node_list in gp_layers.values():
			for nn1, nn2 in itertools.combinations(node_list, 2):
				if g_prime.node_names[nn1].y < g_prime.node_names[nn2].y:
					x_vars_opt[nn1, nn2] = 1
				else:
					x_vars_opt[nn1, nn2] = 0
		for level, node_list in g.get_names_by_layer().items():
			if level in gp_layers:
				stack_nodes = [node for node in node_list if g.node_names[node].stacked]
				for stack_node in stack_nodes:
					for gp_node in gp_layers[level]:
						for node in node_list:
							if node != stack_node:
								if (stack_node, node) in x_vars_opt:
									x_vars_opt[gp_node, node] = x_vars_opt[stack_node, node]
								else:
									x_vars_opt[node, gp_node] = x_vars_opt[node, stack_node]
					for node in node_list:
						if node != stack_node:
							if (stack_node, node) in x_vars_opt:
								del x_vars_opt[stack_node, node]
							else:
								del x_vars_opt[node, stack_node]
		g.unstack_graph()
		y_vars = [n.name for n in g.nodes]
		t3 = time.time() - t3
	else:
		t3 = 0

	""" Draw pre-bendiness graph """
	# for v in m.getVars():
	# 	print('%s %g' % (v.varName, v.x))
	# 	if v.varName[:1] == "y":
	# 		g.node_names[int(v.varName[2:v.varName.index(']')])].y = round(v.x)
	# vis.draw(g, "interim")

	""" Sequantial bendiness reduction """
	if bendiness_reduction:
		t4 = time.time()
		if not do_subg_reduction:
			for v in m.getVars():
				# print('%s %g' % (v.varName, v.x))
				if v.varName[:1] == "x":
					x_vars_opt[int(v.varName[2:v.varName.index(',')]), int(v.varName[v.varName.index(',') + 1:v.varName.index(']')])] = round(v.x)
		n_constraints_generated[4] += sequential_br(g, y_vars, x_vars_opt, m_val, subgraph_seq=not do_subg_reduction)
		t4 = time.time() - t4
		print("Time to perform bendiness reduction: ", t4)
	else:
		t4 = 0
		for v in m.getVars():
			if v.varName[:1] == "y":
				g.node_names[int(v.varName[2:v.varName.index(']')])].y = float(v.x)
	print("Final edge crossing count:", g.num_edge_crossings())

	return n_constraints_generated, [round(t1, 3), round(t2, 3), round(t3, 3), round(t4, 3), round(t1 + t2 + t3 + t4, 2)]
