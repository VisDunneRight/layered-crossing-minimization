import math
import warnings
import inspect
from gurobipy import LinExpr


def get_x_var(x_vars_dict, u1, u2):
	if (u1, u2) in x_vars_dict:
		return x_vars_dict[u1, u2]
	else:
		return 1 - x_vars_dict[u2, u1]


def set_x_var(x_vars_dict, u1, u2, val):
	if (u1, u2) in x_vars_dict:
		x_vars_dict[u1, u2] = val
	else:
		x_vars_dict[u2, u1] = 1 - val


def get_x_var_consts(x_vars_dict, u1, u2):
	if (u1, u2) in x_vars_dict:
		return 1, u1, u2
	else:
		return -1, u2, u1


def get_c_var(c_vars_set, e1, e2):
	if (e1, e2) in c_vars_set:
		return e1, e2
	elif (e2, e1) in c_vars_set:
		return e2, e1
	else:
		raise Exception(f"Edges {e1} and {e2} are not a crossing variable in this model")


def optimization_time_estimate(c_vars_count):
	return math.e ** (math.e ** (0.4615 * math.log(c_vars_count) - 1.542) - 6.9078) if c_vars_count > 0 else 0


def inv_opt_time_estimate(runtime):
	return max(math.e ** ((math.log(runtime) + 1.542) / 0.4615), 0)


def calc_time_taken_for_partition_size(n_p, n_cv, cv_total):
	assert n_p * n_cv <= cv_total
	return n_p * optimization_time_estimate(n_cv / n_p) + optimization_time_estimate(cv_total - n_p * n_cv)


def require(require_true: dict, require_false: dict = None, warn_true: dict = None, warn_false: dict = None, warning_message: str = None):
	if warn_true:
		for k, v in warn_true.items():
			if not v:
				warnings.warn(f"{inspect.stack()[1].function} should be using {k}")
				if warning_message:
					print(warning_message)
	if warn_false:
		for k, v in warn_false.items():
			if v:
				warnings.warn(f"{inspect.stack()[1].function} should be using {k}=False")
				if warning_message:
					print(warning_message)
	for k, v in require_true.items():
		if not v:
			raise Exception(f"{inspect.stack()[1].function} requires {k}")
	if require_false:
		for k, v in require_false.items():
			if v:
				raise Exception(f"{inspect.stack()[1].function} requires {k}=False")


def require_graph_props(g, require_node_data: list = None, warn_node_data: list = None):
	if require_node_data:
		for req in require_node_data:
			if req not in g.node_data or g.node_data[req] is None:
				raise Exception(f"{inspect.stack()[1].function} requires graph property [{req}]")
	if warn_node_data:
		for req in warn_node_data:
			if req not in g.node_data or g.node_data[req] is None:
				warnings.warn(f"{inspect.stack()[1].function} recommends graph property [{req}]")


def calc_x_var_sum(nd, layer_nodes, x_v_set, x):
	xsum = LinExpr()
	for nd_ot in layer_nodes:
		if nd_ot != nd:
			x_r, u1, u2 = get_x_var_consts(x_v_set, nd, nd_ot)
			xsum += x_r * x[u1, u2] + (1 - x_r) // 2
	return xsum
