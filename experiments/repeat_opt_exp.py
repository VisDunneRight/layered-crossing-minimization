from src.benchmark import generate_benchmark
from src.optimization import LayeredOptimizer


def run_func(data_path):
    opt1 = LayeredOptimizer(data_path)
    res1 = opt1.optimize_layout(cutoff_time=60, crossing_minimization=True)
    res2 = opt1.optimize_layout(cutoff_time=60, crossing_minimization=True)
    opt2 = LayeredOptimizer(data_path)
    res3 = opt2.optimize_layout(cutoff_time=60, crossing_minimization=True)
    return res1.runtime, res2.runtime, res3.runtime


exclude_dirs = [f"graficon{i}nodi" for i in range(10, 101) if i not in (20, 21, 30, 31, 40, 41, 51, 61)]
generate_benchmark([], {}, run_func, "../Rome-Lib", name="repeat_opt", exclude_dirs=exclude_dirs, csv_header=["T1", "T2", "T3"], class_dependencies=["src/optimization.LayeredOptimizer"], project_root="/Users/connorwilson/PycharmProjects/stratisfimal-python")
