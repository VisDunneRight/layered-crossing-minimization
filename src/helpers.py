import math


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
