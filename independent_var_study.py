import csv
import math
import numpy as np
from src import read_data, vis


def calc_transform_1(dpt, a=1, b=1):
	gr = read_data.read(dpt[0])
	card_x = 0
	for lay in gr.layers.values():
		card_x += len(lay)**2
	card_c = 0
	e_b_l = gr.get_edges_by_layer()
	for elist in e_b_l.values():
		card_c += len(elist)**2
	return math.log(card_x)**a + math.log(card_c)**b


def calc_transform_2(dpt):
	gr = read_data.read(dpt[0])
	card_x = 0
	for lay in gr.layers.values():
		card_x += len(lay)**2
	card_c = 0
	e_b_l = gr.get_edges_by_layer()
	for elist in e_b_l.values():
		card_c += len(elist)**2
	adj = gr.create_normal_adj_list()
	sqss_deg = 0
	for v in adj.values():
		sqss_deg += len(v)**2
	sqss_deg = math.sqrt(sqss_deg)
	return sqss_deg * (math.log(card_x) + math.log(card_c))


def calc_transform_3(dpt):
	gr = read_data.read(dpt[0])
	card_c = 0
	e_b_l = gr.get_edges_by_layer()
	for elist in e_b_l.values():
		card_c += len(elist)**2
	return math.log(card_c)
	# return card_c


if __name__ == '__main__':
	with open(f"data storage/vertical_transitivity/fix1var_60.csv", 'r') as fd:
		rdr = csv.reader(fd)
		next(rdr)
		basevals = []
		for row in rdr:
			if 0.001 < float(row[10]) < 60:
				basevals.append((row[1], float(row[10]), row[1][:row[1].index('/')]))
	# datapoints = []
	# for bval in basevals:
	# 	datapoints.append({'x': calc_transform_1(bval, 1, 1), 'y': bval[1], 'c': bval[2]})
	# print(f"\nlog base = e")
	# print("#data points:", len(datapoints))
	# # vis.draw_altair_scatter(datapoints, "x", "y", "none", "ln(|X||C|)", "ln(t)", "transform1", True)
	# vis.draw_altair_scatter_with_regression_line(datapoints, "x", "y", "c", f"|X||C|", f"t", f"transform_times")
	# y_bar = sum((dp['y'] for dp in datapoints))/len(datapoints)
	# x_bar = sum((dp['x'] for dp in datapoints))/len(datapoints)
	# b_hat = sum(((dp['x'] - x_bar)*(dp['y'] - y_bar) for dp in datapoints)) / sum(((dp['x'] - x_bar)**2 for dp in datapoints))
	# a_hat = y_bar - (b_hat * x_bar)
	# print(f"y = {b_hat}x + {a_hat}")
	# xy_bar = sum((dp['x']*dp['y'] for dp in datapoints))/len(datapoints)
	# x2_bar = sum((dp['x']**2 for dp in datapoints))/len(datapoints)
	# y2_bar = sum((dp['y']**2 for dp in datapoints))/len(datapoints)
	# r_xy = (xy_bar - (x_bar * y_bar)) / math.sqrt((x2_bar - x_bar**2) * (y2_bar - y_bar**2))
	# print(f"Correlation coefficient: {r_xy}")

	""" exponential version """
	# for a in range(1,5):
	# 	for b in range(1,5):
	# 		datax = []
	# 		datay = []
	# 		colors = []
	# 		for bval in basevals:
	# 			datax.append(calc_transform_1(bval, a=a, b=b))
	# 			datay.append(math.log(bval[1]) + 8)
	# 			colors.append(bval[2])
	# 		datax = np.array(datax)
	# 		datay = np.array(datay)
	#
	# 		results = {}
	# 		coeffs = np.polyfit(datax, np.log(datay), 1, w=np.sqrt(datay))
	# 		p = np.poly1d(coeffs)
	# 		datay = datay - 8
	#
	# 		yhat = np.e**p(datax) - 8
	# 		ybar = np.sum(datay) / len(datay)
	# 		ssreg = np.sum((yhat - ybar) ** 2)
	# 		sstot = np.sum((datay - ybar) ** 2)
	# 		results['r_squared'] = ssreg / sstot
	# 		results['r'] = math.sqrt(results['r_squared'])
	#
	# 		max_v, min_v = max(datax), min(datax)
	# 		line_xv = np.linspace(min_v, max_v, 500)
	# 		line_dicts = [{'x': line_xv[i], 'y': np.e**(p(line_xv[i])) - 8} for i in range(50)]
	# 		scatter_dicts = [{'x': datax[i], 'y': datay[i], 'c': colors[i]} for i in range(len(datax))]
	# 		vis.draw_altair_scatter_with_custom_line(scatter_dicts, line_dicts, 'x', 'y', 'c', f'ln(|X|)^{a} + ln(|C|)^{b}', 'ln(t)', f'transform_{a}_{b}')
	#
	# 		print(p)
	# 		print(results)

	""" linear version compare """
	# for a in range(1, 5):
	# 	for b in range(1, 5):
	# 		datax = []
	# 		datay = []
	# 		colors = []
	# 		for bval in basevals:
	# 			datax.append(calc_transform_1(bval, a=a, b=b))
	# 			datay.append(math.log(bval[1]))
	# 			colors.append(bval[2])
	# 		datax = np.array(datax)
	# 		datay = np.array(datay)
	#
	# 		results = {}
	# 		coeffs = np.polyfit(datax, datay, 1)
	# 		p = np.poly1d(coeffs)
	#
	# 		yhat = p(datax)
	# 		ybar = np.sum(datay) / len(datay)
	# 		ssreg = np.sum((yhat - ybar) ** 2)
	# 		sstot = np.sum((datay - ybar) ** 2)
	# 		results['r_squared'] = ssreg / sstot
	# 		results['r'] = math.sqrt(results['r_squared'])
	#
	# 		max_v, min_v = max(datax), min(datax)
	# 		line_xv = np.linspace(min_v, max_v, 500)
	# 		line_dicts = [{'x': line_xv[i], 'y': p(line_xv[i])} for i in range(500)]
	# 		scatter_dicts = [{'x': datax[i], 'y': datay[i], 'c': colors[i]} for i in range(len(datax))]
	# 		vis.draw_altair_scatter_with_custom_line(scatter_dicts, line_dicts, 'x', 'y', 'c', f'ln(|X|)^{a} + ln(|C|)^{b}', 'ln(t)', f'transform_{a}_{b}')
	#
	# 		print(p)
	# 		print(f"a={a} b={b} results:", results)

	""" linear version only |C| """
	datax = []
	datay = []
	colors = []
	for bval in basevals:
		datax.append(calc_transform_3(bval))
		datay.append(math.log(1000*bval[1]))
		colors.append(bval[2])
	datax = np.array(datax)
	datay = np.array(datay)

	results = {}
	coeffs = np.polyfit(datax, np.log(datay), 1)
	p = np.poly1d(coeffs)

	yhat = np.e**p(datax)
	ybar = np.sum(datay) / len(datay)
	ssreg = np.sum((yhat - ybar) ** 2)
	sstot = np.sum((datay - ybar) ** 2)
	results['r_squared'] = ssreg / sstot
	results['r'] = math.sqrt(results['r_squared'])

	max_v, min_v = max(datax), min(datax)
	line_xv = np.linspace(min_v, max_v, 500)
	line_dicts = [{'x': line_xv[i], 'y': np.e**(p(line_xv[i]))} for i in range(500)]
	scatter_dicts = [{'x': datax[i], 'y': datay[i], 'c': colors[i]} for i in range(len(datax))]
	vis.draw_altair_scatter_with_custom_line(scatter_dicts, line_dicts, 'x', 'y', 'c', 'ln(|C|)', 'ln(t)', f'Cscore_fix1vert')

	print(p)
	print(f"results:", results)
