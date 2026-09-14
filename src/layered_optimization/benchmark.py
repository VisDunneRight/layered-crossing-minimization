import os
import sys
import csv
import inspect
import math
from natsort import natsorted


def generate_benchmark(conditions, condition_map, run_function, path_to_data, name="my_benchmark", in_conds=None, nested_data=-1, exclude_files=None, exclude_dirs=None, allowed_filetypes=None, csv_header=None, dependencies=None, local_dependencies=None, class_dependencies=None, project_root=None, file_sort_key=None):
    """
    :param conditions: List of conditions as strings
    :param condition_map: Dictionary mapping all conditions to a list of values for that condition. Values may be any hashable type
    :param run_function: Run function for the benchmark experiemnt. Must take as input:
            A value for each condition. Parameter name must match name of condition
            Path to a data file. Parameter name must be 'data_path'
        return: List of data to be recorded for that run of the experiment. This will be written to a csv
    :param path_to_data: String. Can be either a flat directory of data files or split into directories each containing data files
    :param in_conds: Optional, list of interior conditions. These conditions will not be parallelized. For each combination of outer conditions, one csv file will be generated containing data for all combinations of in-conditions
    :param nested_data: Optional, specify 0 for flat directory of files or 1 for directories of data, will infer otherwise
    :param exclude_files: Optional, [string | list] of files to exclude from dataset. Dotfiles are already excluded by default
    :param exclude_dirs: Optional, [string | list] of directories to exclude from dataset. Dot-directories are already excluded by default
    :param allowed_filetypes: Optional [string | list of strings] of allowed file extensions for the dataset. Only files of the specified extension will be collected
    :param csv_header: Optional, list of strings which will be added to the header row of all results files, for the results returned by run_function. An index column and columns for each condition will be automatically added
    :param dependencies: Optional, list of strings for all libraries that your run function requires
    :param project_root: String, required only if using local or class dependencies. System absolute path to project root
    :param local_dependencies: Optional, list of strings for local files your run function requires. Specify as path from project root
    :param class_dependencies: Optional, list of strings for local classes your run function requires. Specify as path from project root with the class at the end, e.g. 'src/myfile.MyClass'
    :param file_sort_key: Optional, dictionary for sorting the dataset files. Otherwise will use natural sort on file names
    :param name: String, name for your benchmark. Will be applied to created files and directories
    :return:
    """

    # self.in_conds = in_conds
    # self.out_conds = out_conds
    # self.condition_map = condition_map
    # self.path_to_data = path_to_data
    # self.nested_data = nested_data
    # self.exclude_files = exclude_files
    # self.exclude_dirs = exclude_dirs

    in_conds = [] if in_conds is None else in_conds
    out_conds = [cond for cond in conditions if cond not in in_conds]

    # collect data files and directory info
    n_flat_files = 0
    n_dirs = 0
    dir_paths = []
    file_paths = []

    if exclude_files is None:
        exclude_files = []
    elif type(exclude_files) == str:
        exclude_files = [exclude_files]
    elif type(exclude_files) != list:
        raise TypeError("exclude_files is not a list or string")

    if exclude_dirs is None:
        exclude_dirs = []
    elif type(exclude_dirs) == str:
        exclude_dirs = [exclude_dirs]
    elif type(exclude_dirs) != list:
        raise TypeError("exclude_dirs is not a list or string")

    if allowed_filetypes is None:
        pass
    elif type(allowed_filetypes) == str:
        allowed_filetypes = [allowed_filetypes]
    elif type(allowed_filetypes) != list:
        raise TypeError("allowed_filetypes is not a list or string")

    with os.scandir(path_to_data) as ita:
        for entry in ita:
            if entry.is_dir() and entry.name[0] != '.' and entry.name not in exclude_dirs:
                n_dirs += 1
                dir_paths.append(entry.name)
            elif entry.is_file() and entry.name[0] != '.' and entry.name not in exclude_files and (allowed_filetypes is None or entry.name.split('.')[-1] in allowed_filetypes):
                n_flat_files += 1
                file_paths.append(entry.name)
    dir_paths.sort()

    # determine if nested data directory
    if nested_data == -1:
        data_is_nested = True if 2 ** n_dirs >= n_flat_files else False
    else:
        data_is_nested = bool(nested_data)

    # recollect files if data is nested
    if data_is_nested:
        file_paths.clear()
        for direc in dir_paths:
            files_in_dir = []
            for flname in os.listdir(path_to_data + '/' + direc):
                if flname[0] != '.' and flname not in exclude_files and (allowed_filetypes is None or flname.split('.')[-1] in allowed_filetypes):
                    files_in_dir.append(direc + '/' + flname)
            files_in_dir.sort()
            file_paths += files_in_dir

    if file_sort_key is None:
        file_paths = natsorted(file_paths)
    elif data_is_nested:
        file_paths.sort(key=lambda x: file_sort_key[x[x.rindex('/') + 1:]])
    else:
        file_paths.sort(key=lambda x: file_sort_key[x])
    print("Dataset size:", len(file_paths), "files")

    if not os.path.isdir(f"./{name}"):
        os.mkdir(f"./{name}")
    if not os.path.isdir(f"./{name}/results"):
        os.mkdir(f"./{name}/results")

    # build python executable
    tab = ' ' * 4
    build_code = ["import sys", "import csv", "import os"]

    # import all dependencies
    if dependencies is not None:
        if type(dependencies) == list:
            build_code += [f"import {dep}" for dep in dependencies]
        elif type(dependencies) == str:
            build_code.append(f"import {dependencies}")
        else:
            raise TypeError("Dependencies is not a list or string")
    if local_dependencies is not None:
        if project_root is None:
            raise Exception("Need to specify project root")
        if type(local_dependencies) == list:
            build_code += [f"from {'.'.join(ldep.split('/')[:-1])} import {ldep.split('/')[-1]}" for ldep in local_dependencies]
        elif type(local_dependencies) == str:
            build_code.append(f"from {'.'.join(local_dependencies.split('/')[:-1])} import {local_dependencies.split('/')[-1]}")
        else:
            raise TypeError("Local dependencies is not a list or string")
    if class_dependencies is not None:
        if project_root is None:
            raise Exception("Need to specify project root")
        if type(class_dependencies) == list:
            build_code += [f"from {'.'.join(cdep.replace('.py', '').split('.')[0].split('/'))} import {cdep.replace('.py', '').split('.')[1]}" for cdep in class_dependencies]
        elif type(class_dependencies) == str:
            build_code.append(f"from {'.'.join(class_dependencies.replace('.py', '').split('.')[0].split('/'))} import {class_dependencies.replace('.py', '').split('.')[1]}")
        else:
            raise TypeError("Class dependencies is not a list or string")
    if project_root is not None:
        assert type(project_root) == str
        build_code.append(f"\nsys.path.append('{project_root}')")

    # build helper functions
    build_code += [
        "\n", inspect.getsource(get_start_position), "",
        inspect.getsource(insert_one), "",
        inspect.getsource(run_function), ""
    ]

    # start building run method, initialize variables
    build_code.append("if __name__ == '__main__':")
    build_code += [
        f"{tab}files = {file_paths}",
        f"{tab}prefix = '../{path_to_data}/'",
        f"{tab}conditions = {conditions}",
        f"{tab}csv_header = ['Index', 'Data'] + conditions + {csv_header if csv_header is not None else []}",
        # f"{tab}out_conditions = {out_conds}",
        f"{tab}condition_map = {condition_map}\n",
    ]

    # build local non-poarallelized version first
    build_code.append(f"{tab}if len(sys.argv) != 2:")
    out_c_var_names = [out_c + '_iter' for out_c in out_conds]
    all_c_var_names = [cond + '_iter' for cond in conditions]

    for i, out_c in enumerate(out_conds):
        tabs = tab * (i + 2)
        build_code.append(f"{tabs}for {out_c + '_iter'} in condition_map['{out_c}']:")
    build_code.append(f"{tab * (len(out_conds) + 2)}out_vars = [{', '.join(out_c_var_names)}]")
    build_code.append(f"{tab * (len(out_conds) + 2)}csvpath = './results/' + '_'.join(str(x) for x in out_vars) + '_results.csv'")
    build_code.append(f"{tab * (len(out_conds) + 2)}idx = get_start_position(csvpath, csv_header)")
    build_code.append(f"{tab * (len(out_conds) + 2)}ct = 1")
    for i, in_c in enumerate(in_conds):
        tabs = tab * (i + 2 + len(out_conds))
        build_code.append(f"{tabs}for {in_c + '_iter'} in condition_map['{in_c}']:")
    build_code.append(f"{tab * (len(conditions) + 2)}cur_conds = [{', '.join([str(cond) + '_iter' for cond in conditions])}]")
    build_code.append(f"{tab * (len(conditions) + 2)}print('Combination', ' + '.join([str(cond) for cond in cur_conds]))")
    build_code.append(f"{tab * (len(conditions) + 2)}for file in files:")
    build_code.append(f"{tab * (len(conditions) + 3)}if idx == ct:")
    i_tab = tab * (len(conditions) + 4)
    build_code.append(f"{i_tab}print(len(files) if idx % len(files) == 0 else idx % len(files), '/', len(files))")
    build_code.append(f"{i_tab}out = {run_function.__name__}({'_iter, '.join([str(cond) + '=' + cond for cond in conditions])}{'_iter, ' if len(conditions) >= 1 else ''}data_path=prefix+file)")
    build_code.append(f"{i_tab}insert_one(csvpath, [idx, file] + cur_conds + list(out))")
    build_code.append(f"{i_tab}idx += 1")
    build_code.append(f"{tab * (len(conditions) + 3)}ct += 1")

    # build parallelized code
    build_code.append(f"{tab}else:")

    out_val_code = []
    len_out_conds = [len(condition_map[otc]) for otc in out_conds]
    cumulative_products = [math.prod(len_out_conds[i:]) for i in range(1, len(len_out_conds))] + [1]
    total_parallel_combos = math.prod(len_out_conds)
    for i, out_c in enumerate(out_conds):
        out_val_code.append(f"{tab * 2}{out_c}_val = condition_map['{out_c}'][(pid{''.join([' % ' + str(prd) for prd in cumulative_products[:i]])}) // {cumulative_products[i]}]")
    out_c_var_names = [out_c + '_val' for out_c in out_conds]
    all_c_var_names = [cond + '_val' for cond in conditions]

    build_code += [
        f"{tab * 2}pid = int(sys.argv[1])",
        '\n'.join(out_val_code),
        f"{tab * 2}out_vars = [{', '.join(out_c_var_names)}]",
        f"{tab * 2}csvpath = './results/' + '_'.join(str(x) for x in out_vars) + '_results.csv'",
        f"{tab * 2}idx = get_start_position(csvpath, csv_header)",
        f"{tab * 2}ct = 1",
        f"{tab * 2}start_pos = get_start_position(csvpath, csv_header)"
    ]

    for i, in_c in enumerate(in_conds):
        tabs = tab * (i + 2)
        build_code.append(f"{tabs}for {in_c + '_val'} in condition_map['{in_c}']:")
    i_tab = tab * (len(in_conds) + 4)

    build_code += [
        f"{tab * (len(in_conds) + 2)}cur_conds = [{', '.join([str(cond) + '_val' for cond in conditions])}]",
        f"{tab * (len(in_conds) + 2)}print('Combination', ' + '.join([str(cond) for cond in cur_conds]))",
        f"{tab * (len(in_conds) + 2)}for file in files:",
        f"{tab * (len(in_conds) + 3)}if idx == ct:",
        f"{i_tab}print(len(files) if idx % len(files) == 0 else idx % len(files), '/', len(files))",
        f"{i_tab}out = {run_function.__name__}({'_val, '.join([str(cond) + '=' + cond for cond in conditions])}{'_val, ' if len(conditions) >= 1 else ''}data_path=prefix+file)",
        f"{i_tab}insert_one(csvpath, [idx, file] + [{', '.join(all_c_var_names)}] + list(out))",
        f"{i_tab}idx += 1",
        f"{tab * (len(in_conds) + 3)}ct += 1"
    ]

    with open(f"./{name}/benchmark_exec.py", 'w') as fw:
        fw.writelines(line + '\n' for line in build_code)

    # build bash file for running slurm job in parallel
    build_script = [
        "#!/bin/bash",
        f"#SBATCH -J {name}               # Job name",
        "#SBATCH --partition=short        # Use short partition (24hrs max)",
        "#SBATCH -N 1                   # Number of nodes",
        "#SBATCH -n 1                  # Number of tasks",
        "#SBATCH --time 12:00:00           # Request 12 hours of compute time",
        "#SBATCH --mem=8G           # Request 8GB memory",
        "#SBATCH --constraint=broadwell           # Run job on the broadwell hardware cluster",
        f"#SBATCH -o outputs/{name}_%A_%a.txt       # Standard output file",
        f"#SBATCH -e errors/{name}_%A_%a.txt        # Standard error file",
        "",
        "srun python benchmark_exec.py $SLURM_ARRAY_TASK_ID"
    ]

    with open(f"./{name}/benchmark_slurm.sh", 'w') as fw:
        fw.writelines(line + '\n' for line in build_script)

    # Build local parallelized script
    build_parallel = [
        "# !/bin/bash\n",
        f"for value in {'{'}0..{total_parallel_combos}{'}'}",
        "do",
        "\tpython3 benchmark_exec.py $value &",
        "done"
    ]

    with open(f"./{name}/benchmark_parallel.sh", 'w') as fw:
        fw.writelines(line + '\n' for line in build_parallel)

    # Build README file
    build_readme = [
        "# Benchmarking Instructions",
        "## Running the Benchmark in Series",
        "To run the benchmark on your local machine:",
        f"1. Navigate to the {name} directory in your terminal",
        "2. Run `python3 benchmark_exec.py`\n",
        "Or, open `benchmark_exec.py` in your python IDE of choice and use the run feature\n",
        "## Running the Benchmark in Parallel",
        "To run the benchmark in parallel on your local machine:",
        f"1. Navigate to the {name} directory in your terminal",
        "2. Run `bash benchmark_parallel.sh`",
        "## Running the Benchmark Remotely using SLURM",
        "This is recommended for large benchmark experiments with many combinations of conditions. To run the benchmark in parallel on a remote cluster that uses SLURM:",
        f"1. Access the cluster and navigate to the {name} directory",
        "2. If necessary, activate a python virtual environment using as conda or pip with all required packages installed",
        "3. If necessary, modify `benchmark_slurm.sh` to load any needed modules prior to the `srun` call, e.g. `module load gurobi`"
        f"4. Run the following: `sbatch --array=0-{total_parallel_combos} benchmark_slurm.sh`\n",
        "*If you need to run the benchmark in series using SLURM, use the steps above with the following modifications:",
        "- delete $SLURM_ARRAY_TASK_ID from benchmark_slurm.sh",
        "- run `sbatch benchmark_slurm.sh`"
    ]

    with open(f"./{name}/README.md", 'w') as fw:
        fw.writelines(line + '\n' for line in build_readme)


def get_start_position(csv_filepath, header):
    if os.path.exists(csv_filepath):
        with open(csv_filepath, 'r') as fd:
            for line in fd:
                pass
            last_line = line
            if last_line.split(',')[0] != "Index":
                return int(last_line.split(',')[0]) + 1
            else:
                return 1
    else:
        insert_one(csv_filepath, header)
        return 1


def insert_one(filename, entry):
    with open(filename, 'a', newline='') as f:
        wrt = csv.writer(f)
        wrt.writerow(entry)
