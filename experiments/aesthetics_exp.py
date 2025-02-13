import csv

from src.benchmark import generate_benchmark
from src.optimization import LayeredOptimizer
from src.heuristics import improved_sifting


def run_func(combo_idx, data_path):
    opt = LayeredOptimizer(data_path)
    tlimit = 300
    if combo_idx == 0:  # CR
        res = opt.optimize_layout(cutoff_time=tlimit, crossing_minimization=True)
    elif combo_idx == 1:  # Bend
        res = opt.optimize_layout(cutoff_time=tlimit, bendiness_reduction=True)
    elif combo_idx == 2:  # Angle
        res = opt.optimize_layout(cutoff_time=tlimit, angular_resolution=True)
    elif combo_idx == 3:  # CR Fair
        res = opt.optimize_layout(cutoff_time=tlimit, fairness_constraints=True, fairness_metric="crossings")
    elif combo_idx == 4:  # Bend Fair
        res = opt.optimize_layout(cutoff_time=tlimit, fairness_constraints=True, fairness_metric="bends")
    elif combo_idx == 5:  # SymN
        res = opt.optimize_layout(cutoff_time=tlimit, symmetry_maximization=True, symmetry_maximization_edges=False)
    elif combo_idx == 6:  # SymN+E
        res = opt.optimize_layout(cutoff_time=tlimit, symmetry_maximization=True, symmetry_maximization_edges=True)
    elif combo_idx == 7:  # Bundle
        res = opt.optimize_layout(cutoff_time=tlimit, edge_bundling=True)
    elif combo_idx == 8:  # MinMax
        res = opt.optimize_layout(cutoff_time=tlimit, min_max_crossings=True)
    elif combo_idx == 9:  # CR+Bend
        res = opt.optimize_layout(cutoff_time=tlimit, crossing_minimization=True, bendiness_reduction=True)
    elif combo_idx == 10:  # CR+Angle
        res = opt.optimize_layout(cutoff_time=tlimit, crossing_minimization=True, angular_resolution=True)
    elif combo_idx == 11:  # CR+CRFair
        res = opt.optimize_layout(cutoff_time=tlimit, crossing_minimization=True, fairness_constraints=True, fairness_metric="crossings", gamma_fair=5)
    elif combo_idx == 12:  # CR+FairBend
        res = opt.optimize_layout(cutoff_time=tlimit, crossing_minimization=True, fairness_constraints=True, fairness_metric="bends", gamma_fair=5)
    elif combo_idx == 13:  # CR+SymN
        res = opt.optimize_layout(cutoff_time=tlimit, crossing_minimization=True, symmetry_maximization=True, symmetry_maximization_edges=False)
    elif combo_idx == 14:  # CR+SymN+E
        res = opt.optimize_layout(cutoff_time=tlimit, crossing_minimization=True, symmetry_maximization=True, symmetry_maximization_edges=True)
    elif combo_idx == 15:  # CR+Bundle
        res = opt.optimize_layout(cutoff_time=tlimit, crossing_minimization=True, edge_bundling=True)
    elif combo_idx == 16:  # CR+MinMax
        res = opt.optimize_layout(cutoff_time=tlimit, crossing_minimization=True, min_max_crossings=True, gamma_min_max=5)
    elif combo_idx == 17:  # Angle (fixed x)
        improved_sifting(opt.g)
        res = opt.optimize_layout(cutoff_time=tlimit, fix_x_vars=True, angular_resolution=True)
    elif combo_idx == 18:  # FairBend (fixed x)
        improved_sifting(opt.g)
        res = opt.optimize_layout(cutoff_time=tlimit, fix_x_vars=True, fairness_constraints=True, fairness_metric="bends")
    elif combo_idx == 19:  # SymN (fixed x)
        improved_sifting(opt.g)
        res = opt.optimize_layout(cutoff_time=tlimit, fix_x_vars=True, symmetry_maximization=True, symmetry_maximization_edges=False)
    elif combo_idx == 20:  # SymN+E (fixed x)
        improved_sifting(opt.g)
        res = opt.optimize_layout(cutoff_time=tlimit, fix_x_vars=True, symmetry_maximization=True, symmetry_maximization_edges=True)
    elif combo_idx == 21:  # Angle (fixed x) + Bend
        improved_sifting(opt.g)
        res = opt.optimize_layout(cutoff_time=tlimit, fix_x_vars=True, bendiness_reduction=True, angular_resolution=True)
    elif combo_idx == 22:  # FairBend (fixed x) + Bend
        improved_sifting(opt.g)
        res = opt.optimize_layout(cutoff_time=tlimit, fix_x_vars=True, bendiness_reduction=True, fairness_constraints=True, fairness_metric="bends")
    elif combo_idx == 23:  # SymN (fixed x) + Bend
        improved_sifting(opt.g)
        res = opt.optimize_layout(cutoff_time=tlimit, fix_x_vars=True, bendiness_reduction=True, symmetry_maximization=True, symmetry_maximization_edges=False)
    elif combo_idx == 24:  # SymN+E (fixed x) + Bend
        improved_sifting(opt.g)
        res = opt.optimize_layout(cutoff_time=tlimit, fix_x_vars=True, bendiness_reduction=True, symmetry_maximization=True, symmetry_maximization_edges=True)

    gcr = opt.g.num_edge_crossings()
    bnds = sum(abs(e.n1.y - e.n2.y) for e in opt.g.edges)
    return combo_idx, sum(1 for nd in opt.g.nodes if not nd.is_anchor_node), opt.g.n_nodes, res.objval, gcr, bnds, res.runtime, res.status


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


generate_benchmark(["combo_idx"], {"combo_idx": list(range(25))}, run_func, "../random graphs/networkx", name="aesthetics", csv_header=["ComboIdx", "Nodes", "TotalNodes", "ObjVal", "Crossings", "EdgeLength", "Runtime", "Status"], class_dependencies=["src/optimization.LayeredOptimizer"], project_root="/Users/connorwilson/PycharmProjects/stratisfimal-python")
