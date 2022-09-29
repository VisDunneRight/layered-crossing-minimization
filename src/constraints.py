import itertools
import gurobipy as gp
from gurobipy import GRB
from src import graph, vis
import time


def optimize_layout(g: graph.LayeredGraph, gamma_1, gamma_2, bendiness_reduction, cutoff_time=0):
    t1 = time.time()
    nodes_by_layer = g.get_names_by_layer()
    edges_by_layer = g.get_edge_names_by_layer()
    n_constraints_generated = [0]*5  # simple edge, hybrid edge, same layer edge, vertical pos, bendiness
    m_val = 50

    try:
        m = gp.Model()

        # add all variables
        x_vars = []
        z_vars = []
        x_vars_layers = {}
        for i, name_list in nodes_by_layer.items():
            x_vars += list(itertools.combinations(name_list, 2))
            z_vars += list(itertools.permutations(name_list, 2))
            x_vars_layers[i] = list(itertools.combinations(name_list, 2))
        x = m.addVars(x_vars, vtype=GRB.BINARY, name="x")
        # z = m.addVars(z_vars, vtype=GRB.INTEGER, lb=0, ub=m_val, name="z")

        c_vars = []
        for edge_list in edges_by_layer.values():
            for comb in itertools.combinations(edge_list, 2):
                if comb[0][0] != comb[1][0] and comb[0][1] != comb[1][1] and comb[0][0] != comb[1][1] and comb[0][1] != comb[1][0]:
                    c_vars.append(comb)
        c = m.addVars(c_vars, vtype=GRB.BINARY, name="c")

        y_vars = [n.name for n in g.nodes]
        y = m.addVars(y_vars, vtype=GRB.INTEGER, lb=0, ub=m_val, name="y")

        # b_vars = list(g.edge_names.keys())
        # b = m.addVars(b_vars, vtype=GRB.INTEGER, lb=0, ub=m_val, name="b")

        # optimization function
        # if bendiness_reduction:
        #     m.setObjective(gamma_1*c.sum() + gamma_2*b.sum(), GRB.MINIMIZE)
        # else:
        #     m.setObjective(c.sum(), GRB.MINIMIZE)
        m.setObjective(c.sum(), GRB.MINIMIZE)

        # transitivity
        # for x_list in x_vars_layers.values():
        #     if len(x_list) >= 3:
        #         for comb in itertools.combinations(x_list, 3):
        #             x2_rev, x3_rev = 1, 1
        #             if comb[0][1] == comb[1][0] or comb[0][1] == comb[1][1]:
        #                 x2 = comb[1]
        #                 x3 = comb[2]
        #                 if comb[0][1] == comb[1][1]:
        #                     x2_rev = -1
        #                 if comb[0][0] == comb[2][1]:
        #                     x3_rev = -1
        #             else:
        #                 x2 = comb[2]
        #                 x3 = comb[1]
        #                 if comb[0][1] == comb[2][1]:
        #                     x2_rev = -1
        #                 if comb[0][0] == comb[1][1]:
        #                     x3_rev = -1
        #             m.addConstr(x[comb[0]] + x2_rev*x[x2] - x3_rev*x[x3] + (1-x2_rev)/2 - (1-x3_rev)/2 >= 0, f"1t{comb[0]},{comb[1]},{comb[2]}")
        #             m.addConstr(-x[comb[0]]-x2_rev*x[x2]+x3_rev*x[x3] - (1-x2_rev)/2 + (1-x3_rev)/2 >= -1, f"2t{comb[0]},{comb[1]},{comb[2]}")
        #             n_constraints_generated[0] += 2

        # def get_xval(i_x, j_x):
        #     if i_x == j_x:
        #         return 0
        #     elif (i_x, j_x) in x:
        #         return x[i_x, j_x]
        #     return 1 - x[j_x, i_x]
        #
        # xsum_vars = []
        # xsum_checkvars = []
        # for name_list in nodes_by_layer.values():
        #     if len(name_list) > 2:
        #         xsum_vars += name_list
        #         xsum_checkvars += list(itertools.product(name_list, range(len(name_list))))
        # xsum = m.addVars(xsum_vars, vtype=GRB.INTEGER, name="xsum")
        # xsum_check = m.addVars(xsum_checkvars, vtype=GRB.BINARY, name="xsum_check")
        #
        # # new transitivity
        # for node_list in nodes_by_layer.values():
        #     if len(node_list) >= 3:
        #         for v in node_list:
        #             m.addConstr(xsum[v] == sum((get_xval(v, ot) for ot in node_list)), f"2t{v}")
        # for i, x_list in nodes_by_layer.items():
        #     if len(x_list) >= 3:
        #         for x_var in x_vars_layers[i]:
        #             m.addConstr(xsum[x_var[0]] != xsum[x_var[1]], f"23t{x_var}")

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
                m.addConstr((1-x1_rev*x[x1]) + x2_rev*x[x2] + c[c_var] - (1-x1_rev)/2 + (1-x2_rev)/2 >= 1, f"1se{c_var}")
                m.addConstr(x1_rev*x[x1] + (1-x2_rev*x[x2]) + c[c_var] + (1-x1_rev)/2 - (1-x2_rev)/2 >= 1, f"2se{c_var}")
                n_constraints_generated[0] += 2

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
                m.addConstr(-x1_rev*x[x1] + x2_rev*x[x2] + c[c_var] - (1-x1_rev)/2 + (1-x2_rev)/2 >= 0, f"1hy{c_var}")
                m.addConstr(x1_rev*x[x1] - x2_rev*x[x2] + c[c_var] + (1-x1_rev)/2 - (1-x2_rev)/2 >= 0, f"2hy{c_var}")
                n_constraints_generated[1] += 2

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
                m.addConstr(-x1_rev*x[x1] + x2_rev*x[x2] + c[c_var] - (1-x1_rev)/2 + (1-x2_rev)/2 >= 0, f"1hy{c_var}")
                m.addConstr(x1_rev*x[x1] - x2_rev*x[x2] + c[c_var] + (1-x1_rev)/2 - (1-x2_rev)/2 >= 0, f"2hy{c_var}")
                n_constraints_generated[1] += 2

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

                m.addConstr(c[c_var] + 1-x1_rev*x[x1] + 1-x4_rev*x[x4] + 1-x2_rev*x[x2] - (1-x1_rev)/2 - (1-x4_rev)/2 - (1-x2_rev)/2 >= 1,f"1sl{c_var}")
                m.addConstr(c[c_var] + 1-x3_rev*x[x3] + x2_rev*x[x2] + x4_rev*x[x4] - (1-x3_rev)/2 + (1-x2_rev)/2 + (1-x4_rev)/2 >= 1,f"2sl{c_var}")
                m.addConstr(c[c_var] + x1_rev*x[x1] + 1-x3_rev*x[x3] + x2_rev*x[x2] + (1-x1_rev)/2 - (1-x3_rev)/2 + (1-x2_rev)/2 >= 1,f"3sl{c_var}")
                m.addConstr(c[c_var] + 1-x4_rev*x[x4] + 1-x2_rev*x[x2] + x3_rev*x[x3] - (1-x4_rev)/2 - (1-x2_rev)/2 + (1-x3_rev)/2 >= 1,f"4sl{c_var}")
                m.addConstr(c[c_var] + x4_rev*x[x4] + x1_rev*x[x1] + 1-x3_rev*x[x3] + (1-x4_rev)/2 + (1-x1_rev)/2 - (1-x3_rev)/2 >= 1,f"5sl{c_var}")
                m.addConstr(c[c_var] + 1-x2_rev*x[x2] + x3_rev*x[x3] + 1-x1_rev*x[x1] - (1-x2_rev)/2 + (1-x3_rev)/2 - (1-x1_rev)/2 >= 1,f"6sl{c_var}")
                m.addConstr(c[c_var] + x2_rev*x[x2] + x4_rev*x[x4] + x1_rev*x[x1] + (1-x2_rev)/2 + (1-x4_rev)/2 + (1-x1_rev)/2 >= 1,f"7sl{c_var}")
                m.addConstr(c[c_var] + x3_rev*x[x3] + 1-x1_rev*x[x1] + 1-x4_rev*x[x4] + (1-x3_rev)/2 - (1-x1_rev)/2 - (1-x4_rev)/2 >= 1,f"8sl{c_var}")

                n_constraints_generated[2] += 8

        # vertical position
        for x_var in x_vars:
            m.addGenConstrIndicator(x[x_var], True, y[x_var[0]] + 1 <= y[x_var[1]])
            m.addGenConstrIndicator(x[x_var], False, y[x_var[0]] >= 1 + y[x_var[1]])
            n_constraints_generated[3] += 2

        # if bendiness_reduction:
            # vertical position
            # for x_var in x_vars:
            #     # print(x_var)
            #     m.addConstr(z[x_var] - m_val * x[x_var] <= 0, f"1.1z{x_var}")
            #     m.addConstr(z[x_var] - y[x_var[0]] - (m_val * x[x_var]) >= -1 * m_val, f"1.2z{x_var}")
            #     m.addConstr(y[x_var[1]] - z[x_var] - x[x_var] >= 0, f"1.3z{x_var}")
            #     m.addConstr(z[x_var] <= y[x_var[0]], f"1.4z{x_var}")
            #     # m.addConstr(z[x_var] >= 0, f"1.5z{x_var}")
            #     # print(f"z{x_var}-{m_val}*x{x_var} <= 0")
            #     # print(f"y{x_var[1]} - z{x_var} - x{x_var} >= 0")
            #
            #     m.addConstr(z[x_var[1], x_var[0]] - m_val * (1 - x[x_var]) <= 0, f"2.1z{x_var}")
            #     m.addConstr(z[x_var[1], x_var[0]] - y[x_var[1]] - m_val * (1 - x[x_var]) >= -1 * m_val,
            #                 f"2.2z{x_var}")
            #     m.addConstr(y[x_var[0]] - z[x_var[1], x_var[0]] - (1 - x[x_var]) >= 0, f"2.3z{x_var}")
            #     m.addConstr(z[x_var[1], x_var[0]] <= y[x_var[1]], f"2.4z{x_var}")
            #     # m.addConstr(z[x_var[1],x_var[0]] >= 0, f"2.5z{x_var}")
            #     # print(f"z[{x_var[1],x_var[0]}]-{m_val}*(1-x{x_var}) <= 0")
            #     # print(f"y{x_var[0]} - z[{x_var[1],x_var[0]}] - (1-x{x_var}) >= 0")
            #
            #     n_constraints_generated[3] += 10

            # bendiness
            # for b_var in b_vars:
            #     m.addConstr(y[b_var[0]] - y[b_var[1]] <= b[b_var], f"1bend{b_var}")
            #     m.addConstr(y[b_var[1]] - y[b_var[0]] <= b[b_var], f"2bend{b_var}")
            #     n_constraints_generated[4] += 2

        t1 = time.time() - t1
        print("Time to input constraints: ", t1)
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
        t2 = time.time() - t2
        print("Time to optimize: ", t2)

        # for v in m.getVars():
        #     # print('%s %g' % (v.varName, v.x))
        #     if v.varName[:1] == "y":
        #         g.node_names[v.varName[2:v.varName.index(']')]].y = int(v.x)
        # vis.draw(g, "interim")

        if bendiness_reduction:
            t3 = time.time()
            x_vars_opt = {}
            for v in m.getVars():
                # print('%s %g' % (v.varName, v.x))
                if v.varName[:1] == "x":
                    x_vars_opt[v.varName[2:v.varName.index(',')], v.varName[v.varName.index(',') + 1:v.varName.index(']')]] = int(v.x)
            m2 = gp.Model()
            y = m2.addVars(y_vars, vtype=GRB.INTEGER, lb=0, ub=m_val, name="y")
            b_vars = list(g.edge_names.keys())
            b = m2.addVars(b_vars, vtype=GRB.INTEGER, lb=0, ub=m_val, name="b")
            m2.setObjective(b.sum(), GRB.MINIMIZE)
            for var, val in x_vars_opt.items():
                if val == 0:
                    m2.addConstr(y[var[0]] >= 1 + y[var[1]], f"vert{var}")
                else:
                    m2.addConstr(y[var[0]] + 1 <= y[var[1]], f"vert{var}")
            for b_var in b_vars:
                m2.addConstr(y[b_var[0]] - y[b_var[1]] <= b[b_var], f"1bend{b_var}")
                m2.addConstr(y[b_var[1]] - y[b_var[0]] <= b[b_var], f"2bend{b_var}")
                n_constraints_generated[4] += 2
            m2.setParam("OutputFlag", 0)
            m2.optimize()
            for v in m2.getVars():
                if v.varName[:1] == "y":
                    g.node_names[v.varName[2:v.varName.index(']')]].y = int(v.x)
            t3 = time.time() - t3
            print("Time to perform bendiness reduction: ", t3)
        else:
            t3 = 0
            for v in m.getVars():
                # print('%s %g' % (v.varName, v.x))
                if v.varName[:1] == "y":
                    g.node_names[v.varName[2:v.varName.index(']')]].y = int(v.x)
            # print('Obj: %g' % m.objVal)

        return n_constraints_generated, t1, t2, t3

    except gp.GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))

    # except AttributeError:
    #     print('Encountered an attribute error')
