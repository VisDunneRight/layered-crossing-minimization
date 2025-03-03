import csv
import random
from src.benchmark import generate_benchmark
from src.optimization import LayeredOptimizer
from src.heuristics import improved_sifting
import pickle


def run_func(combo_idx, data_path):
    opt = LayeredOptimizer(data_path)
    tlimit = 600
    combo_map = ["CR", "Bend", "Angle", "CRFair", "BendFair", "SymN", "SymNE", "Bundle", "MinMax", "MinEwCr", "CR+Bend", "CR+Angle", "CR+CRFair", "CR+BendFair", "CR+SymN", "CR+SymNE", "CR+Bundle", "CR+MinMax", "CR+MinEwCr", "Bend (fixed x)", "Angle (fixed x)", "BendFair (fixed x)", "SymN (fixed x)", "SymNE (fixed x)", "Angle+Bend (fixed x)", "BendFair+Bend (fixed x)", "SymN+Bend (fixed x)", "SymNE+Bend (fixed x)"]
    combo_choice = combo_map[combo_idx]
    if combo_choice == "CR":
        res = opt.optimize_layout(cutoff_time=tlimit, crossing_minimization=True)
    elif combo_choice == "Bend":
        res = opt.optimize_layout(cutoff_time=tlimit, bendiness_reduction=True)
    elif combo_choice == "Angle":
        res = opt.optimize_layout(cutoff_time=tlimit, angular_resolution=True)
    elif combo_choice == "CRFair":
        fair_vals = random.choices([0, 1], k=opt.g.n_nodes)
        while 0 not in fair_vals or 1 not in fair_vals:
            fair_vals = random.choices([0, 1], k=opt.g.n_nodes)
        fair_nds = [[v for v in opt.g.node_ids if fair_vals[v] == 0], [v for v in opt.g.node_ids if fair_vals[v] == 1]]
        opt.g.add_fairness_values(fair_nds)
        res = opt.optimize_layout(cutoff_time=tlimit, fairness_constraints=True, fairness_metric="crossings")
    elif combo_choice == "BendFair":
        fair_vals = random.choices([0, 1], k=opt.g.n_nodes)
        while 0 not in fair_vals or 1 not in fair_vals:
            fair_vals = random.choices([0, 1], k=opt.g.n_nodes)
        fair_nds = [[v for v in opt.g.node_ids if fair_vals[v] == 0], [v for v in opt.g.node_ids if fair_vals[v] == 1]]
        opt.g.add_fairness_values(fair_nds)
        res = opt.optimize_layout(cutoff_time=tlimit, fairness_constraints=True, fairness_metric="bends")
    elif combo_choice == "SymN":
        res = opt.optimize_layout(cutoff_time=tlimit, symmetry_maximization=True, symmetry_maximization_edges=False)
    elif combo_choice == "SymNE":
        res = opt.optimize_layout(cutoff_time=tlimit, symmetry_maximization=True, symmetry_maximization_edges=True)
    elif combo_choice == "Bundle":
        res = opt.optimize_layout(cutoff_time=tlimit, edge_bundling=True)
    elif combo_choice == "MinMax":
        res = opt.optimize_layout(cutoff_time=tlimit, min_max_crossings=True)
    elif combo_choice == "MinEwCr":
        res = opt.optimize_layout(cutoff_time=tlimit, min_edges_with_crossings=True)
    elif combo_choice == "CR+Bend":
        res = opt.optimize_layout(cutoff_time=tlimit, crossing_minimization=True, bendiness_reduction=True)
    elif combo_choice == "CR+Angle":
        res = opt.optimize_layout(cutoff_time=tlimit, crossing_minimization=True, angular_resolution=True)
    elif combo_choice == "CR+CRFair":
        fair_vals = random.choices([0, 1], k=opt.g.n_nodes)
        while 0 not in fair_vals or 1 not in fair_vals:
            fair_vals = random.choices([0, 1], k=opt.g.n_nodes)
        fair_nds = [[v for v in opt.g.node_ids if fair_vals[v] == 0], [v for v in opt.g.node_ids if fair_vals[v] == 1]]
        opt.g.add_fairness_values(fair_nds)
        res = opt.optimize_layout(cutoff_time=tlimit, crossing_minimization=True, fairness_constraints=True,
                                  fairness_metric="crossings", gamma_fair=5)
    elif combo_choice == "CR+BendFair":
        fair_vals = random.choices([0, 1], k=opt.g.n_nodes)
        while 0 not in fair_vals or 1 not in fair_vals:
            fair_vals = random.choices([0, 1], k=opt.g.n_nodes)
        fair_nds = [[v for v in opt.g.node_ids if fair_vals[v] == 0], [v for v in opt.g.node_ids if fair_vals[v] == 1]]
        opt.g.add_fairness_values(fair_nds)
        res = opt.optimize_layout(cutoff_time=tlimit, crossing_minimization=True, fairness_constraints=True,
                                  fairness_metric="bends", gamma_fair=5)
    elif combo_choice == "CR+SymN":
        res = opt.optimize_layout(cutoff_time=tlimit, crossing_minimization=True, symmetry_maximization=True,
                                  symmetry_maximization_edges=False)
    elif combo_choice == "CR+SymNE":
        res = opt.optimize_layout(cutoff_time=tlimit, crossing_minimization=True, symmetry_maximization=True,
                                  symmetry_maximization_edges=True)
    elif combo_choice == "CR+Bundle":
        res = opt.optimize_layout(cutoff_time=tlimit, crossing_minimization=True, edge_bundling=True)
    elif combo_choice == "CR+MinMax":
        res = opt.optimize_layout(cutoff_time=tlimit, crossing_minimization=True, min_max_crossings=True,
                                  gamma_min_max=5)
    elif combo_choice == "CR+MinEwCr":
        res = opt.optimize_layout(cutoff_time=tlimit, crossing_minimization=True, min_edges_with_crossings=True)
    elif combo_choice == "Bend (fixed x)":
        improved_sifting(opt.g)
        res = opt.optimize_layout(cutoff_time=tlimit, fix_x_vars=True, bendiness_reduction=True)
    elif combo_choice == "Angle (fixed x)":
        improved_sifting(opt.g)
        res = opt.optimize_layout(cutoff_time=tlimit, fix_x_vars=True, angular_resolution=True)
    elif combo_choice == "BendFair (fixed x)":
        improved_sifting(opt.g)
        fair_vals = random.choices([0, 1], k=opt.g.n_nodes)
        while 0 not in fair_vals or 1 not in fair_vals:
            fair_vals = random.choices([0, 1], k=opt.g.n_nodes)
        fair_nds = [[v for v in opt.g.node_ids if fair_vals[v] == 0], [v for v in opt.g.node_ids if fair_vals[v] == 1]]
        opt.g.add_fairness_values(fair_nds)
        res = opt.optimize_layout(cutoff_time=tlimit, fix_x_vars=True, fairness_constraints=True,
                                  fairness_metric="bends")
    elif combo_choice == "SymN (fixed x)":
        improved_sifting(opt.g)
        res = opt.optimize_layout(cutoff_time=tlimit, fix_x_vars=True, symmetry_maximization=True,
                                  symmetry_maximization_edges=False)
    elif combo_choice == "SymNE (fixed x)":
        improved_sifting(opt.g)
        res = opt.optimize_layout(cutoff_time=tlimit, fix_x_vars=True, symmetry_maximization=True,
                                  symmetry_maximization_edges=True)
    elif combo_choice == "Angle+Bend (fixed x)":
        improved_sifting(opt.g)
        res = opt.optimize_layout(cutoff_time=tlimit, fix_x_vars=True, bendiness_reduction=True,
                                  angular_resolution=True)
    elif combo_choice == "BendFair+Bend (fixed x)":
        improved_sifting(opt.g)
        fair_vals = random.choices([0, 1], k=opt.g.n_nodes)
        while 0 not in fair_vals or 1 not in fair_vals:
            fair_vals = random.choices([0, 1], k=opt.g.n_nodes)
        fair_nds = [[v for v in opt.g.node_ids if fair_vals[v] == 0], [v for v in opt.g.node_ids if fair_vals[v] == 1]]
        opt.g.add_fairness_values(fair_nds)
        res = opt.optimize_layout(cutoff_time=tlimit, fix_x_vars=True, bendiness_reduction=True,
                                  fairness_constraints=True, fairness_metric="bends")
    elif combo_choice == "SymN+Bend (fixed x)":
        improved_sifting(opt.g)
        res = opt.optimize_layout(cutoff_time=tlimit, fix_x_vars=True, bendiness_reduction=True,
                                  symmetry_maximization=True, symmetry_maximization_edges=False)
    elif combo_choice == "SymNE+Bend (fixed x)":
        improved_sifting(opt.g)
        res = opt.optimize_layout(cutoff_time=tlimit, fix_x_vars=True, bendiness_reduction=True,
                                  symmetry_maximization=True, symmetry_maximization_edges=True)

    gcr = opt.g.num_edge_crossings()
    bnds = sum(abs(e.n1.y - e.n2.y) for e in opt.g.edges)
    return sum(1 for nd in opt.g.nodes if not nd.is_anchor_node), sum(1 for ed in opt.g.edges if not ed.n1.is_anchor_node), opt.g.n_nodes, len(opt.g.edges), res.objval, gcr, bnds, res.runtime, res.status


def calculate_cutoff(csv_file, num_nodes, files_per_bucket):
    with open(csv_file, 'r') as fd:
        rdr = csv.reader(fd)
        first_line = next(rdr)
        fl_idx, st_idx = 1, first_line.index("Status")
        nfls = 0
        n_cutoff = 0
        for ln in rdr:
            if int(ln[fl_idx].split('_')[1]) == num_nodes:
                nfls += 1
                if int(ln[st_idx]) != 2:
                    n_cutoff += 1
        if nfls != files_per_bucket:
            raise Exception(f"Wrong num files in bucket {num_nodes}, {files_per_bucket} != {nfls}")
    return n_cutoff / nfls


def calculate_cutoff_rome(csv_file, num_nodes_b, files_in_bucket):
    with open(csv_file, 'r') as fd:
        rdr = csv.reader(fd)
        first_line = next(rdr)
        tnodes_idx, st_idx = first_line.index("TotalNodes"), first_line.index("Status")
        nfls = 0
        n_cutoff = 0
        for ln in rdr:
            if int(ln[tnodes_idx]) // 5 * 5 == num_nodes_b:
                nfls += 1
                if int(ln[st_idx]) != 2:
                    n_cutoff += 1
        if nfls != files_in_bucket:
            raise Exception(f"Wrong num files in bucket {num_nodes_b}, {files_in_bucket} != {nfls}")
    return n_cutoff / nfls


def get_prev_result_if_exists(file_path, graph_name):
    with open(file_path, 'r') as fd:
        for line in fd:
            if line.split(',')[1] == graph_name:
                if line.split(',')[11] == "2":
                    return line
                else:
                    return []


fls_per_bucket = {10: 229, 15: 343, 20: 404, 25: 355, 30: 301, 35: 258, 40: 273, 45: 303, 50: 297, 55: 360, 60: 371, 65: 351, 70: 345, 75: 359, 80: 317, 85: 301, 90: 317, 95: 253, 100: 252, 105: 292, 110: 200, 115: 277, 120: 216, 125: 201, 130: 204, 135: 193, 140: 204, 145: 149, 150: 151, 155: 161, 160: 162, 165: 143, 170: 182, 175: 162, 180: 120, 185: 158, 190: 142, 195: 143, 200: 141, 205: 123, 210: 144, 215: 126, 220: 121, 225: 121, 230: 118, 235: 108, 240: 105, 245: 97, 250: 81}
fls_to_exclude = ['grafo3836.53', 'grafo3821.55', 'grafo4872.60', 'grafo7230.61', 'grafo9815.62', 'grafo8076.62', 'grafo4428.64', 'grafo9321.64', 'grafo9374.65', 'grafo9433.65', 'grafo4609.65', 'grafo4820.66', 'grafo1300.66', 'grafo9154.66', 'grafo7905.67', 'grafo9169.67', 'grafo5230.69', 'grafo7869.69', 'grafo1233.70', 'grafo9799.70', 'grafo8973.70', 'grafo9117.71', 'grafo8672.71', 'grafo4422.71', 'grafo5144.71', 'grafo8467.71', 'grafo8523.71', 'grafo4925.71', 'grafo9279.72', 'grafo7461.72', 'grafo4830.72', 'grafo4819.72', 'grafo9063.73', 'grafo4771.73', 'grafo8142.73', 'grafo8678.73', 'grafo7970.74', 'grafo8545.74', 'grafo9153.74', 'grafo5902.74', 'grafo7864.75', 'grafo9502.75', 'grafo7993.75', 'grafo9290.75', 'grafo4674.75', 'grafo8961.75', 'grafo4725.75', 'grafo7973.75', 'grafo4665.76', 'grafo4917.76', 'grafo8530.76', 'grafo5866.76', 'grafo8194.76', 'grafo8085.76', 'grafo9631.76', 'grafo4493.76', 'grafo4672.76', 'grafo8109.77', 'grafo7829.77', 'grafo9076.77', 'grafo9371.77', 'grafo7849.77', 'grafo4968.77', 'grafo5606.77', 'grafo9786.77', 'grafo4778.77', 'grafo4613.77', 'grafo4574.77', 'grafo5439.78', 'grafo8321.78', 'grafo8263.78', 'grafo4264.78', 'grafo8131.78', 'grafo5996.78', 'grafo8312.78', 'grafo9789.78', 'grafo4148.78', 'grafo5383.78', 'grafo8772.79', 'grafo9482.79', 'grafo9781.79', 'grafo4395.79', 'grafo9770.79', 'grafo9721.79', 'grafo5681.79', 'grafo8141.79', 'grafo9557.79', 'grafo9453.80', 'grafo7908.80', 'grafo9047.80', 'grafo4437.80', 'grafo7839.80', 'grafo7879.80', 'grafo4735.80', 'grafo8805.80', 'grafo4696.80', 'grafo4507.80', 'grafo5051.80', 'grafo8057.80', 'grafo9747.80', 'grafo8247.81', 'grafo8290.81', 'grafo8086.81', 'grafo9228.81', 'grafo9846.81', 'grafo8429.81', 'grafo8266.81', 'grafo9674.81', 'grafo4745.81', 'grafo9486.81', 'grafo6903.82', 'grafo8073.82', 'grafo8275.82', 'grafo9787.82', 'grafo7450.82', 'grafo9436.82', 'grafo8428.82', 'grafo9412.82', 'grafo7451.82', 'grafo4137.82', 'grafo8316.82', 'grafo8357.82', 'grafo9212.82', 'grafo9307.82', 'grafo9237.82', 'grafo8574.83', 'grafo9230.83', 'grafo8118.83', 'grafo9685.83', 'grafo7942.83', 'grafo4202.83', 'grafo7940.83', 'grafo4608.83', 'grafo8272.83', 'grafo5028.83', 'grafo6550.83', 'grafo9306.83', 'grafo7918.84', 'grafo6670.84', 'grafo9355.84', 'grafo8440.84', 'grafo8370.84', 'grafo9450.84', 'grafo7917.84', 'grafo9726.84', 'grafo8284.84', 'grafo8528.84', 'grafo8600.84', 'grafo8537.84', 'grafo9252.84', 'grafo9171.84', 'grafo9331.85', 'grafo8486.85', 'grafo8551.85', 'grafo8666.85', 'grafo9550.85', 'grafo9824.85', 'grafo9464.85', 'grafo7981.85', 'grafo8024.85', 'grafo7944.85', 'grafo8281.85', 'grafo7824.85', 'grafo9004.85', 'grafo8976.85', 'grafo6766.85', 'grafo8567.85', 'grafo8485.85', 'grafo9661.86', 'grafo8448.86', 'grafo8620.86', 'grafo8725.86', 'grafo8432.86', 'grafo9316.86', 'grafo8765.86', 'grafo8687.86', 'grafo8495.86', 'grafo7100.86', 'grafo7101.86', 'grafo8191.86', 'grafo8768.86', 'grafo8260.86', 'grafo7671.86', 'grafo8310.86', 'grafo8578.87', 'grafo9542.87', 'grafo8764.87', 'grafo8601.87', 'grafo8993.87', 'grafo8804.87', 'grafo8602.87', 'grafo8726.87', 'grafo9027.87', 'grafo8009.87', 'grafo8729.87', 'grafo8355.87', 'grafo3871.88', 'grafo8014.88', 'grafo8083.88', 'grafo8021.88', 'grafo8877.88', 'grafo8813.88', 'grafo8810.88', 'grafo8733.88', 'grafo8830.88', 'grafo8897.88', 'grafo8330.88', 'grafo8177.88', 'grafo7378.88', 'grafo8352.89', 'grafo7411.89', 'grafo8909.89', 'grafo8818.89', 'grafo3139.89', 'grafo9432.89', 'grafo8394.89', 'grafo8477.89', 'grafo8664.89', 'grafo9241.89', 'grafo8643.89', 'grafo8747.89', 'grafo8500.89', 'grafo7936.89', 'grafo8351.89', 'grafo9676.89', 'grafo11648.90', 'grafo8778.90', 'grafo8474.90', 'grafo10076.90', 'grafo8712.90', 'grafo8850.90', 'grafo8737.90', 'grafo11647.90', 'grafo8851.90', 'grafo8744.90', 'grafo8923.90', 'grafo8711.90', 'grafo8848.90', 'grafo9823.90', 'grafo11390.90', 'grafo9936.90', 'grafo8938.90', 'grafo7780.90', 'grafo3756.90', 'grafo8226.90', 'grafo9854.91', 'grafo11300.91', 'grafo9986.91', 'grafo11350.91', 'grafo8002.91', 'grafo11488.91', 'grafo9987.91', 'grafo9008.91', 'grafo9993.91', 'grafo9983.91', 'grafo11315.91', 'grafo11534.91', 'grafo11487.91', 'grafo8179.91', 'grafo8801.91', 'grafo11470.91', 'grafo10639.91', 'grafo8662.91', 'grafo11520.91', 'grafo11471.91', 'grafo8954.91', 'grafo11505.91', 'grafo8647.91', 'grafo11616.91', 'grafo11472.91', 'grafo10419.91', 'grafo8786.91', 'grafo8650.91', 'grafo10312.91', 'grafo3184.91', 'grafo11655.91', 'grafo11537.91', 'grafo8802.91', 'grafo11453.91', 'grafo11600.91', 'grafo11507.91', 'grafo10317.91', 'grafo10502.91', 'grafo9994.91', 'grafo11448.91', 'grafo8808.91', 'grafo10624.91', 'grafo9937.91', 'grafo9985.91', 'grafo9894.91', 'grafo8922.92', 'grafo11635.92', 'grafo11379.92', 'grafo8491.92', 'grafo8651.92', 'grafo8218.92', 'grafo11439.92', 'grafo10912.92', 'grafo8377.92', 'grafo9832.92', 'grafo10572.92', 'grafo9576.92', 'grafo11381.92', 'grafo10740.92', 'grafo10614.92', 'grafo11267.92', 'grafo11366.92', 'grafo9634.92', 'grafo9299.92', 'grafo10470.92', 'grafo10893.92', 'grafo11190.92', 'grafo11365.92', 'grafo8215.92', 'grafo6307.92', 'grafo11220.92', 'grafo11393.92', 'grafo8659.92', 'grafo10713.92', 'grafo11435.92', 'grafo10999.92', 'grafo10245.92', 'grafo11389.92', 'grafo10326.93', 'grafo7354.93', 'grafo11580.93', 'grafo10291.93', 'grafo8823.93', 'grafo8981.93', 'grafo10408.93', 'grafo11480.93', 'grafo8715.93', 'grafo11591.93', 'grafo8675.93', 'grafo11269.93', 'grafo10967.93', 'grafo10813.93', 'grafo10953.93', 'grafo3464.93', 'grafo11186.93', 'grafo11569.93', 'grafo10413.93', 'grafo8878.93', 'grafo11114.93', 'grafo10527.93', 'grafo11303.93', 'grafo11290.93', 'grafo10517.93', 'grafo10507.93', 'grafo10903.93', 'grafo11222.93', 'grafo11297.93', 'grafo10860.93', 'grafo10541.93', 'grafo7985.93', 'grafo10681.93', 'grafo10883.93', 'grafo11274.93', 'grafo11007.93', 'grafo11270.93', 'grafo11215.93', 'grafo11184.93', 'grafo11095.93', 'grafo11195.93', 'grafo10834.93', 'grafo8719.93', 'grafo6242.93', 'grafo11133.93', 'grafo11080.93', 'grafo11251.93', 'grafo10410.93', 'grafo11275.93', 'grafo7529.93', 'grafo10270.93', 'grafo10588.93', 'grafo8951.93', 'grafo8800.93', 'grafo10589.93', 'grafo10162.93', 'grafo10107.93', 'grafo10608.93', 'grafo11069.93', 'grafo10360.93', 'grafo10669.93', 'grafo3143.94', 'grafo10266.94', 'grafo8953.94', 'grafo8512.94', 'grafo8833.94', 'grafo10518.94', 'grafo10151.94', 'grafo11219.94', 'grafo10548.94', 'grafo8863.94', 'grafo10391.94', 'grafo10362.94', 'grafo11249.94', 'grafo11189.94', 'grafo10154.94', 'grafo10198.94', 'grafo11025.94', 'grafo10349.94', 'grafo10547.94', 'grafo8809.94', 'grafo8949.94', 'grafo10994.94', 'grafo10099.94', 'grafo11070.94', 'grafo10973.94', 'grafo8969.94', 'grafo10594.94', 'grafo11071.94', 'grafo11045.94', 'grafo10917.94', 'grafo10832.94', 'grafo10933.94', 'grafo10630.94', 'grafo10797.94', 'grafo10309.94', 'grafo10249.94', 'grafo10992.94', 'grafo10870.94', 'grafo10904.94', 'grafo10955.94', 'grafo11205.94', 'grafo10850.94', 'grafo10493.94', 'grafo10637.94', 'grafo10482.94', 'grafo10901.94', 'grafo7990.94', 'grafo10785.94', 'grafo11185.94', 'grafo10606.94', 'grafo8214.94', 'grafo8036.94', 'grafo10254.94', 'grafo10382.94', 'grafo8821.94', 'grafo8571.94', 'grafo8900.94', 'grafo10301.94', 'grafo8791.94', 'grafo10113.94', 'grafo10127.94', 'grafo6817.94', 'grafo8476.95', 'grafo11557.95', 'grafo10306.95', 'grafo10390.95', 'grafo10949.95', 'grafo9928.95', 'grafo10819.95', 'grafo10114.95', 'grafo8490.95', 'grafo10849.95', 'grafo10121.95', 'grafo10243.95', 'grafo10150.95', 'grafo8890.95', 'grafo10688.95', 'grafo8655.95', 'grafo10543.95', 'grafo10807.95', 'grafo11104.95', 'grafo10379.95', 'grafo3712.95', 'grafo10427.95', 'grafo10782.95', 'grafo10943.95', 'grafo11100.95', 'grafo11141.95', 'grafo11217.95', 'grafo10635.95', 'grafo10446.95', 'grafo10697.95', 'grafo10610.95', 'grafo10675.95', 'grafo10927.95', 'grafo10867.95', 'grafo11101.95', 'grafo9630.95', 'grafo10480.95', 'grafo11105.95', 'grafo10228.95', 'grafo10537.95', 'grafo8789.95', 'grafo10661.95', 'grafo10880.95', 'grafo10431.95', 'grafo8334.95', 'grafo10388.95', 'grafo10149.95', 'grafo10451.95', 'grafo10109.95', 'grafo11181.95', 'grafo10288.95', 'grafo11117.95', 'grafo10993.95', 'grafo10851.95', 'grafo10925.95', 'grafo10684.95', 'grafo10139.95', 'grafo10636.95', 'grafo10560.95', 'grafo8514.95', 'grafo10489.95', 'grafo10371.95', 'grafo10898.95', 'grafo11178.95', 'grafo10224.95', 'grafo10241.95', 'grafo11112.96', 'grafo10811.96', 'grafo10627.96', 'grafo10722.96', 'grafo10961.96', 'grafo10757.96', 'grafo10914.96', 'grafo11057.96', 'grafo10501.96', 'grafo10944.96', 'grafo11002.96', 'grafo10531.96', 'grafo10694.96', 'grafo10987.96', 'grafo7792.96', 'grafo10492.96', 'grafo11012.96', 'grafo10028.96', 'grafo8811.96', 'grafo10122.96', 'grafo10499.96', 'grafo10779.96', 'grafo8780.96', 'grafo10330.96', 'grafo11039.96', 'grafo8893.96', 'grafo8997.96', 'grafo8667.96', 'grafo10123.96', 'grafo3879.96', 'grafo10246.96', 'grafo11491.96', 'grafo10409.96', 'grafo11369.96', 'grafo10367.96', 'grafo10104.96', 'grafo10252.96', 'grafo9347.96', 'grafo3315.96', 'grafo10303.96', 'grafo10160.96', 'grafo10578.96', 'grafo11098.96', 'grafo10823.96', 'grafo10720.96', 'grafo9984.96', 'grafo11021.96', 'grafo11097.96', 'grafo10687.96', 'grafo9973.96', 'grafo11087.96', 'grafo10952.96', 'grafo10634.96', 'grafo10966.96', 'grafo10416.96', 'grafo10686.96', 'grafo10764.96', 'grafo10744.96', 'grafo11469.96', 'grafo10760.96', 'grafo11313.96', 'grafo10815.97', 'grafo10716.97', 'grafo10118.97', 'grafo11046.97', 'grafo11017.97', 'grafo10179.97', 'grafo11076.97', 'grafo10405.97', 'grafo10420.97', 'grafo10810.97', 'grafo10974.97', 'grafo10404.97', 'grafo10905.97', 'grafo10486.97', 'grafo11006.97', 'grafo3995.97', 'grafo10921.97', 'grafo10915.97', 'grafo10361.97', 'grafo8790.97', 'grafo6585.97', 'grafo10194.97', 'grafo10205.97', 'grafo10365.97', 'grafo10351.97', 'grafo10240.97', 'grafo10173.97', 'grafo10072.97', 'grafo10275.97', 'grafo10146.97', 'grafo10387.97', 'grafo10156.97', 'grafo10176.97', 'grafo10220.97', 'grafo10321.97', 'grafo10084.97', 'grafo10130.97', 'grafo8856.97', 'grafo10322.97', 'grafo10216.97', 'grafo10273.97', 'grafo10828.97', 'grafo11546.97', 'grafo10206.97', 'grafo10346.97', 'grafo10226.97', 'grafo10092.97', 'grafo10212.97', 'grafo10859.97', 'grafo10919.97', 'grafo8867.97', 'grafo10327.97', 'grafo11001.97', 'grafo10491.97', 'grafo8699.97', 'grafo10683.97', 'grafo10495.97', 'grafo10947.97', 'grafo10822.97', 'grafo11092.97', 'grafo10806.97', 'grafo10705.97', 'grafo10877.97', 'grafo11096.97', 'grafo11010.97', 'grafo10523.97', 'grafo10913.97', 'grafo10894.97', 'grafo10406.97', 'grafo10825.98', 'grafo10695.98', 'grafo10128.98', 'grafo11143.98', 'grafo11102.98', 'grafo7793.98', 'grafo10685.98', 'grafo10924.98', 'grafo11056.98', 'grafo11072.98', 'grafo10945.98', 'grafo10805.98', 'grafo10623.98', 'grafo10831.98', 'grafo10672.98', 'grafo10960.98', 'grafo10441.98', 'grafo10814.98', 'grafo11016.98', 'grafo10434.98', 'grafo10647.98', 'grafo10496.98', 'grafo10667.98', 'grafo11067.98', 'grafo11152.98', 'grafo3979.98', 'grafo11009.98', 'grafo6562.98', 'grafo8864.98', 'grafo10648.98', 'grafo10297.98', 'grafo10112.98', 'grafo10126.98', 'grafo10988.98', 'grafo10260.98', 'grafo6278.98', 'grafo8921.98', 'grafo8673.98', 'grafo10006.98', 'grafo10488.98', 'grafo8875.98', 'grafo10709.98', 'grafo10559.98', 'grafo10689.98', 'grafo11198.98', 'grafo10207.98', 'grafo10380.98', 'grafo10093.98', 'grafo10908.98', 'grafo8943.98', 'grafo10918.98', 'grafo10144.98', 'grafo10479.98', 'grafo10529.98', 'grafo10948.98', 'grafo10808.98', 'grafo8566.98', 'grafo10519.98', 'grafo10978.98', 'grafo3224.98', 'grafo10693.98', 'grafo10631.98', 'grafo8035.98', 'grafo10833.98', 'grafo10792.98', 'grafo8478.98', 'grafo10650.98', 'grafo10936.98', 'grafo10876.98', 'grafo11050.98', 'grafo10318.98', 'grafo10866.98', 'grafo10700.98', 'grafo10580.98', 'grafo11187.98', 'grafo11140.98', 'grafo10581.98', 'grafo11055.98', 'grafo3482.98', 'grafo10596.99', 'grafo10780.99', 'grafo10854.99', 'grafo11106.99', 'grafo3532.99', 'grafo11250.99', 'grafo11037.99', 'grafo10520.99', 'grafo10864.99', 'grafo10736.99', 'grafo11254.99', 'grafo11142.99', 'grafo11271.99', 'grafo10737.99', 'grafo11204.99', 'grafo10663.99', 'grafo10293.99', 'grafo11169.99', 'grafo10274.99', 'grafo10386.99', 'grafo11008.99', 'grafo6625.99', 'grafo10354.99', 'grafo10383.99', 'grafo3382.99', 'grafo10213.99', 'grafo10353.99', 'grafo9357.99', 'grafo10161.99', 'grafo8198.99', 'grafo10468.99', 'grafo10217.99', 'grafo11188.99', 'grafo10175.99', 'grafo11248.99', 'grafo11099.99', 'grafo9054.99', 'grafo10100.99', 'grafo10549.99', 'grafo10307.99', 'grafo11358.99', 'grafo8378.99', 'grafo8857.99', 'grafo8812.99', 'grafo7221.99', 'grafo10366.99', 'grafo3898.99', 'grafo8348.99', 'grafo11223.99', 'grafo11110.99', 'grafo10219.99', 'grafo11035.99', 'grafo10532.99', 'grafo10803.99', 'grafo10674.99', 'grafo10926.99', 'grafo10209.99', 'grafo11316.99', 'grafo11182.99', 'grafo8223.99', 'grafo11145.99', 'grafo10536.99', 'grafo10605.99', 'grafo11155.99', 'grafo8000.99', 'grafo10907.99', 'grafo10484.99', 'grafo10269.99', 'grafo3680.99', 'grafo11000.99', 'grafo3861.99', 'grafo10696.99', 'grafo10641.99', 'grafo10432.99', 'grafo10735.99', 'grafo11150.99', 'grafo11230.100', 'grafo10503.100', 'grafo10265.100', 'grafo10305.100', 'grafo10316.100', 'grafo11020.100', 'grafo10262.100', 'grafo10538.100', 'grafo8793.100', 'grafo9738.100', 'grafo10671.100', 'grafo10467.100', 'grafo11157.100', 'grafo10248.100', 'grafo10116.100', 'grafo10603.100', 'grafo10429.100', 'grafo8720.100', 'grafo11657.100', 'grafo11093.100', 'grafo10204.100', 'grafo11324.100', 'grafo8905.100', 'grafo8865.100', 'grafo10629.100', 'grafo11335.100', 'grafo10638.100', 'grafo10772.100', 'grafo10820.100', 'grafo11138.100', 'grafo10970.100', 'grafo11065.100', 'grafo8674.100', 'grafo10184.100', 'grafo11073.100', 'grafo8892.100', 'grafo10998.100', 'grafo10556.100', 'grafo10426.100', 'grafo10750.100', 'grafo8510.100', 'grafo10962.100', 'grafo10786.100', 'grafo10584.100', 'grafo8882.100', 'grafo8855.100', 'grafo10633.100', 'grafo10223.100', 'grafo10237.100', 'grafo8513.100', 'grafo11458.100', 'grafo10418.100', 'grafo10632.100', 'grafo7797.100', 'grafo10183.100', 'grafo10331.100', 'grafo10292.100', 'grafo8549.100', 'grafo10250.100', 'grafo10865.100', 'grafo10124.100', 'grafo9053.100', 'grafo8826.100', 'grafo11038.100', 'grafo10508.100', 'grafo10534.100', 'grafo10911.100', 'grafo11011.100', 'grafo10469.100', 'grafo10294.100', 'grafo10928.100', 'grafo11014.100', 'grafo10444.100', 'grafo10652.100', 'grafo10646.100', 'grafo10644.100', 'grafo11162.100', 'grafo11003.100']
with open("../localfiles/fl_to_tnodes", 'rb') as fd:
    sortkey = pickle.load(fd)

generate_benchmark(["combo_idx"], {"combo_idx": list(range(28))}, run_func, "../Rome-Lib", name="aesthetics", exclude_files=fls_to_exclude, csv_header=["Nodes", "TotalNodes", "ObjVal", "Crossings", "EdgeLength", "Runtime", "Status"], class_dependencies=["src/optimization.LayeredOptimizer"], project_root="/Users/connorwilson/PycharmProjects/stratisfimal-python", file_sort_key=sortkey)
