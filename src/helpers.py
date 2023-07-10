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


def get_c_var(c_vars_dict, u, v):
	if (u, v) in c_vars_dict:
		return u, v
	elif (v, u) in c_vars_dict:
		return v, u
	else:
		raise Exception("Edges are not a crossing variable in this model")
