import itertools
import gurobipy as gp
from gurobipy import GRB
from src import graph, vis, reductions
import time


def sequential_br(g, y_vars, x_var_opt, mv):
	n_constr = 0
	m2 = gp.Model()
	y = m2.addVars(y_vars, vtype=GRB.INTEGER, lb=0, ub=mv, name="y")
	b_vars = list(g.edge_names.keys())
	b = m2.addVars(b_vars, vtype=GRB.INTEGER, lb=0, ub=mv, name="b")
	m2.setObjective(b.sum(), GRB.MINIMIZE)
	for var, val in x_var_opt.items():
		if val == 0:
			m2.addConstr(y[var[0]] >= 1 + y[var[1]], f"vert{var}")
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
			g.node_names[int(v.varName[2:v.varName.index(']')])].y = int(v.x)
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
			model.addConstr((1 - x1_rev * x[x1]) + x2_rev * x[x2] + c[c_var] - (1 - x1_rev) / 2 + (1 - x2_rev) / 2 >= 1,
							f"1se{c_var}")
			model.addConstr(x1_rev * x[x1] + (1 - x2_rev * x[x2]) + c[c_var] + (1 - x1_rev) / 2 - (1 - x2_rev) / 2 >= 1,
							f"2se{c_var}")
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
			model.addConstr(-x1_rev * x[x1] + x2_rev * x[x2] + c[c_var] - (1 - x1_rev) / 2 + (1 - x2_rev) / 2 >= 0,
							f"1hy{c_var}")
			model.addConstr(x1_rev * x[x1] - x2_rev * x[x2] + c[c_var] + (1 - x1_rev) / 2 - (1 - x2_rev) / 2 >= 0,
							f"2hy{c_var}")
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
			model.addConstr(-x1_rev * x[x1] + x2_rev * x[x2] + c[c_var] - (1 - x1_rev) / 2 + (1 - x2_rev) / 2 >= 0,
							f"1hy{c_var}")
			model.addConstr(x1_rev * x[x1] - x2_rev * x[x2] + c[c_var] + (1 - x1_rev) / 2 - (1 - x2_rev) / 2 >= 0,
							f"2hy{c_var}")
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
	# for keylist in keys:
	return n_constr_0, n_constr_1, n_constr_2


def transitivity(model: gp.Model, x_var_layers, x, n_b_l):
	n_constr = 0
	for x_list in x_var_layers.values():
		if len(x_list) >= 3:
			for comb in itertools.combinations(x_list, 3):
				x2_rev, x3_rev = 1, 1
				if comb[0][1] == comb[1][0] or comb[0][1] == comb[1][1]:
					x2 = comb[1]
					x3 = comb[2]
					if comb[0][1] == comb[1][1]:
						x2_rev = -1
					if comb[0][0] == comb[2][1]:
						x3_rev = -1
				else:
					x2 = comb[2]
					x3 = comb[1]
					if comb[0][1] == comb[2][1]:
						x2_rev = -1
					if comb[0][0] == comb[1][1]:
						x3_rev = -1
				model.addConstr(x[comb[0]] + x2_rev * x[x2] - x3_rev * x[x3] + (1 - x2_rev) / 2 - (1 - x3_rev) / 2 >= 0,
								f"1t{comb[0]},{comb[1]},{comb[2]}")
				model.addConstr(
					-x[comb[0]] - x2_rev * x[x2] + x3_rev * x[x3] - (1 - x2_rev) / 2 + (1 - x3_rev) / 2 >= -1,
					f"2t{comb[0]},{comb[1]},{comb[2]}")
				n_constr += 2

	# def get_xval(i_x, j_x):
	# 	if i_x == j_x:
	# 		return 0
	# 	elif (i_x, j_x) in x:
	# 		return x[i_x, j_x]
	# 	return 1 - x[j_x, i_x]
	#
	# xsum_vars = []
	# xsum_checkvars = []
	# for name_list in n_b_l.values():
	# 	if len(name_list) > 2:
	# 		xsum_vars += name_list
	# 		xsum_checkvars += list(itertools.product(name_list, range(len(name_list))))
	# xsum = model.addVars(xsum_vars, vtype=GRB.INTEGER, name="xsum")
	# xsum_check = model.addVars(xsum_checkvars, vtype=GRB.BINARY, name="xsum_check")
	#
	# # new transitivity
	# for node_list in n_b_l.values():
	# 	if len(node_list) >= 3:
	# 		for v in node_list:
	# 			model.addConstr(xsum[v] == sum((get_xval(v, ot) for ot in node_list)), f"2t{v}")
	# for i, x_list in n_b_l.items():
	# 	if len(x_list) >= 3:
	# 		for x_var in x_var_layers[i]:
	# 			model.addConstr(xsum[x_var[0]] != xsum[x_var[1]], f"23t{x_var}")
	return n_constr


def optimize_layout(g: graph.LayeredGraph, gamma_1, gamma_2, bendiness_reduction, sequential_bendiness=True,
					cutoff_time=0):
	t1 = time.time()
	nodes_by_layer = g.get_names_by_layer()
	edges_by_layer = g.get_edge_names_by_layer()
	n_constraints_generated = [0] * 6  # simple edge, hybrid edge, same layer edge, vertical pos, bendiness, total
	m_val = 50

	m = gp.Model()

	""" Add all variables """
	x_vars = []
	x_vars_layers = {}
	for i, name_list in nodes_by_layer.items():
		x_vars += list(itertools.combinations(name_list, 2))
		x_vars_layers[i] = list(itertools.combinations(name_list, 2))
	x = m.addVars(x_vars, vtype=GRB.BINARY, name="x")
	c_vars = []
	c_vars_layers = {}
	for i, edge_list in edges_by_layer.items():
		c_vars_layers[i] = []
		for comb in itertools.combinations(edge_list, 2):
			if comb[0][0] != comb[1][0] and comb[0][1] != comb[1][1] and comb[0][0] != comb[1][1] and comb[0][1] != \
					comb[1][0]:
				c_vars.append(comb)
				c_vars_layers[i].append(comb)
	# c2_vars, to_remove = reductions.crossing_var_sum_reduction_finder(g, nodes_by_layer, c_vars_layers)
	# c_vars_update = [c for c in c_vars if c not in to_remove]
	c = m.addVars(c_vars, vtype=GRB.BINARY, name="c")
	y_vars = [n.name for n in g.nodes]
	y = m.addVars(y_vars, vtype=GRB.INTEGER, lb=0, ub=m_val, name="y")

	""" Non-sequential optimization function """
	if not sequential_bendiness:
		z_vars = []
		for i, name_list in nodes_by_layer.items():
			z_vars += list(itertools.permutations(name_list, 2))
		z = m.addVars(z_vars, vtype=GRB.INTEGER, lb=0, ub=m_val, name="z")
		if bendiness_reduction:
			b_vars = list(g.edge_names.keys())
			b = m.addVars(b_vars, vtype=GRB.INTEGER, lb=0, ub=m_val, name="b")
			m.setObjective(gamma_1 * c.sum() + gamma_2 * b.sum(), GRB.MINIMIZE)
		else:
			m.setObjective(c.sum(), GRB.MINIMIZE)
	else:
		m.setObjective(c.sum(), GRB.MINIMIZE)

	""" Transitivity attempts """
	# n_constraints_generated[0] += transitivity(m, x_vars_layers, x, nodes_by_layer)

	""" Long-version crossing reduction code """
	n_cs = edge_crossings(m, c_vars, g, x, c)
	for i, val in enumerate(n_cs):
		n_constraints_generated[i] += val

	""" Vertical position, implication version """
	for x_var in x_vars:
		m.addGenConstrIndicator(x[x_var], True, y[x_var[0]] + 1 <= y[x_var[1]])
		m.addGenConstrIndicator(x[x_var], False, y[x_var[0]] >= 1 + y[x_var[1]])
		n_constraints_generated[3] += 2

	""" Non-sequential bendiness reduction """
	if not sequential_bendiness and bendiness_reduction:
		for x_var in x_vars:
			# print(x_var)
			m.addConstr(z[x_var] - m_val * x[x_var] <= 0, f"1.1z{x_var}")
			m.addConstr(z[x_var] - y[x_var[0]] - (m_val * x[x_var]) >= -1 * m_val, f"1.2z{x_var}")
			m.addConstr(y[x_var[1]] - z[x_var] - x[x_var] >= 0, f"1.3z{x_var}")
			m.addConstr(z[x_var] <= y[x_var[0]], f"1.4z{x_var}")
			# m.addConstr(z[x_var] >= 0, f"1.5z{x_var}")
			# print(f"z{x_var}-{m_val}*x{x_var} <= 0")
			# print(f"y{x_var[1]} - z{x_var} - x{x_var} >= 0")

			m.addConstr(z[x_var[1], x_var[0]] - m_val * (1 - x[x_var]) <= 0, f"2.1z{x_var}")
			m.addConstr(z[x_var[1], x_var[0]] - y[x_var[1]] - m_val * (1 - x[x_var]) >= -1 * m_val,
						f"2.2z{x_var}")
			m.addConstr(y[x_var[0]] - z[x_var[1], x_var[0]] - (1 - x[x_var]) >= 0, f"2.3z{x_var}")
			m.addConstr(z[x_var[1], x_var[0]] <= y[x_var[1]], f"2.4z{x_var}")
			# m.addConstr(z[x_var[1],x_var[0]] >= 0, f"2.5z{x_var}")
			# print(f"z[{x_var[1],x_var[0]}]-{m_val}*(1-x{x_var}) <= 0")
			# print(f"y{x_var[0]} - z[{x_var[1],x_var[0]}] - (1-x{x_var}) >= 0")

			n_constraints_generated[3] += 10

		# bendiness
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
	# m.tune()
	# for i in range(m.tuneResultCount):
	#     m.getTuneResult(i)
	#     m.write('tune' + str(i) + '.prm')
	if cutoff_time > 0:
		m.setParam("TimeLimit", cutoff_time)
	m.setParam("OutputFlag", 0)
	m.optimize()
	print('Obj: %g' % m.objVal)
	t2 = time.time() - t2
	print("Time to optimize: ", t2)

	""" Draw pre-bendiness graph """
	# for v in m.getVars():
	#     # print('%s %g' % (v.varName, v.x))
	#     if v.varName[:1] == "y":
	#         g.node_names[v.varName[2:v.varName.index(']')]].y = int(v.x)
	# vis.draw(g, "interim")

	""" Sequantial bendiness reduction """
	if bendiness_reduction:
		t3 = time.time()
		x_vars_opt = {}
		for v in m.getVars():
			# print('%s %g' % (v.varName, v.x))
			if v.varName[:1] == "x":
				x_vars_opt[int(v.varName[2:v.varName.index(',')]), int(v.varName[v.varName.index(',') + 1:v.varName.index(']')])] = int(v.x)
		n_constraints_generated[4] += sequential_br(g, y_vars, x_vars_opt, m_val)
		t3 = time.time() - t3
		print("Time to perform bendiness reduction: ", t3)
	else:
		t3 = 0
		for v in m.getVars():
			# print('%s %g' % (v.varName, v.x))
			if v.varName[:1] == "y":
				g.node_names[int(v.varName[2:v.varName.index(']')])].y = int(v.x)
	# print('Obj: %g' % m.objVal)

	return n_constraints_generated, [round(t1, 3), round(t2, 3), round(t3, 3), round(t1+t2+t3, 2)]

# except gp.GurobiError as e:
# 	print('Error code ' + str(e.errno) + ': ' + str(e))
# except AttributeError:
#     print('Encountered an attribute error')
