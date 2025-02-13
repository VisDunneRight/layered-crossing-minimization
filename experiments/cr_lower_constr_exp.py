from src.benchmark import generate_benchmark
from src.optimization import LayeredOptimizer


def run_func(data_path):
    opt1 = LayeredOptimizer(data_path, crossing_lower_constraints=False)
    res1 = opt1.optimize_layout(cutoff_time=300, crossing_minimization=True)
    opt2 = LayeredOptimizer(data_path, crossing_lower_constraints=True)
    res2 = opt2.optimize_layout(cutoff_time=300, crossing_minimization=True)
    return res1.runtime, res2.runtime


exclude_dirs = [f"graficon{i}nodi" for i in range(10, 75)] + [f"graficon{i}nodi" for i in range(77, 101)]
generate_benchmark([], {}, run_func, "../Rome-Lib", name="cr_lower_constraint", exclude_dirs=exclude_dirs, csv_header=["Tstandard", "Twithconstraint"], class_dependencies=["src/optimization.LayeredOptimizer"], project_root="/Users/connorwilson/PycharmProjects/stratisfimal-python")
