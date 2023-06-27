from scipy.optimize import linprog
import numpy as np
from src import motifs
from src.optimization import LayeredOptimizer


class HiGHSLayeredOptimizer(LayeredOptimizer):
	def __init__(self, layered_graph, parameters=None):
		LayeredOptimizer.__init__(self, layered_graph, parameters=parameters)
		if self.heuristic_start:
			print("Heuristic start only available in Gurobi implementation")
		del self.heuristic_start

	def optimize_layout(self):
		out = self.__optimize_layout_std()
		if self.verbose:
			for string in self.print_info:
				print(string)
			self.print_info.clear()
		return out

	# TODO: implement. https://docs.scipy.org/doc/scipy/reference/optimize.linprog-highs.html#optimize-linprog-highs
	def __optimize_layout_std(self):
		g = self.g
		if not self.direct_transitivity and not self.vertical_transitivity:
			self.vertical_transitivity = True
		nodes_by_layer = g.get_names_by_layer()
		edges_by_layer = g.get_edge_names_by_layer()
		pre_sym = ''
		if self.verbose:
			self.print_info.extend(('-' * 70, f"{self.name}:", '-' * 70))

		""" Variables, constraints setup """

		""" Butterfly reduction """

		""" Set model objective function """

		""" Transitivity constraints """
		# if self.direct_transitivity:

		""" Long-version crossing reduction code """

		""" Symmetry constraints """

		""" Fix key x-var """
		# if self.fix_one_var:

		""" Original vertical position constraints """
		# if self.vertical_transitivity:

		""" Non-sequential bendiness reduction, original Stratisfimal version"""
		# if not self.sequential_bendiness and self.bendiness_reduction:

		""" Optimize model """

		""" Sequential bendiness reduction """

		""" Draw optimized graph """
		# if self.draw_graph:

		""" Add print info """
		# if self.verbose:
		# 	self.print_info.append(f"{pre_sym}Final edge crossing count: {g.num_edge_crossings()}")
		# 	self.print_info.append(f"{pre_sym}Number of constraints: {n_constraints_generated}")
		# 	self.print_info.append(f"{pre_sym}{round(t1, 3)}, {round(t2, 3)}, {round(t3, 3)}, {round(t1 + t2 + t3, 3)}")

		# print(f"Number of crossings: {round(m.objVal) if m.objVal != float('inf') else m.objVal}",
		# 	  f"\tOptimization time: {round(m.runtime, 3)}")

		""" Return data """
		n_cr = 0
		return n_cr
