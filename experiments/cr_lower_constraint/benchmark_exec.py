import sys
import csv
import os
from src.optimization import LayeredOptimizer

sys.path.append('/Users/connorwilson/PycharmProjects/stratisfimal-python')


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


def run_func(data_path):
    opt1 = LayeredOptimizer(data_path, crossing_lower_constraints=False)
    res1 = opt1.optimize_layout(cutoff_time=300, crossing_minimization=True)
    opt2 = LayeredOptimizer(data_path, crossing_lower_constraints=True)
    res2 = opt2.optimize_layout(cutoff_time=300, crossing_minimization=True)
    return res1.runtime, res2.runtime


if __name__ == '__main__':
    files = ['graficon75nodi/grafo1451.75', 'graficon75nodi/grafo4156.75', 'graficon75nodi/grafo4187.75', 'graficon75nodi/grafo4203.75', 'graficon75nodi/grafo4233.75', 'graficon75nodi/grafo4253.75', 'graficon75nodi/grafo4274.75', 'graficon75nodi/grafo4276.75', 'graficon75nodi/grafo4292.75', 'graficon75nodi/grafo4355.75', 'graficon75nodi/grafo4383.75', 'graficon75nodi/grafo4449.75', 'graficon75nodi/grafo4519.75', 'graficon75nodi/grafo4537.75', 'graficon75nodi/grafo4640.75', 'graficon75nodi/grafo4669.75', 'graficon75nodi/grafo4674.75', 'graficon75nodi/grafo4683.75', 'graficon75nodi/grafo4724.75', 'graficon75nodi/grafo4725.75', 'graficon75nodi/grafo4797.75', 'graficon75nodi/grafo4802.75', 'graficon75nodi/grafo4870.75', 'graficon75nodi/grafo4903.75', 'graficon75nodi/grafo4989.75', 'graficon75nodi/grafo4992.75', 'graficon75nodi/grafo5037.75', 'graficon75nodi/grafo5060.75', 'graficon75nodi/grafo5169.75', 'graficon75nodi/grafo5263.75', 'graficon75nodi/grafo5280.75', 'graficon75nodi/grafo5363.75', 'graficon75nodi/grafo5434.75', 'graficon75nodi/grafo5540.75', 'graficon75nodi/grafo5588.75', 'graficon75nodi/grafo5594.75', 'graficon75nodi/grafo5859.75', 'graficon75nodi/grafo5933.75', 'graficon75nodi/grafo5938.75', 'graficon75nodi/grafo6005.75', 'graficon75nodi/grafo6016.75', 'graficon75nodi/grafo6486.75', 'graficon75nodi/grafo6681.75', 'graficon75nodi/grafo6700.75', 'graficon75nodi/grafo6772.75', 'graficon75nodi/grafo6954.75', 'graficon75nodi/grafo6995.75', 'graficon75nodi/grafo7803.75', 'graficon75nodi/grafo7864.75', 'graficon75nodi/grafo7926.75', 'graficon75nodi/grafo7937.75', 'graficon75nodi/grafo7964.75', 'graficon75nodi/grafo7973.75', 'graficon75nodi/grafo7993.75', 'graficon75nodi/grafo8029.75', 'graficon75nodi/grafo8056.75', 'graficon75nodi/grafo8084.75', 'graficon75nodi/grafo8099.75', 'graficon75nodi/grafo8120.75', 'graficon75nodi/grafo8162.75', 'graficon75nodi/grafo8166.75', 'graficon75nodi/grafo8174.75', 'graficon75nodi/grafo8197.75', 'graficon75nodi/grafo8235.75', 'graficon75nodi/grafo8343.75', 'graficon75nodi/grafo8381.75', 'graficon75nodi/grafo8406.75', 'graficon75nodi/grafo8426.75', 'graficon75nodi/grafo8433.75', 'graficon75nodi/grafo8458.75', 'graficon75nodi/grafo8459.75', 'graficon75nodi/grafo8466.75', 'graficon75nodi/grafo8493.75', 'graficon75nodi/grafo8507.75', 'graficon75nodi/grafo8517.75', 'graficon75nodi/grafo8524.75', 'graficon75nodi/grafo8546.75', 'graficon75nodi/grafo8548.75', 'graficon75nodi/grafo8622.75', 'graficon75nodi/grafo8633.75', 'graficon75nodi/grafo8634.75', 'graficon75nodi/grafo8638.75', 'graficon75nodi/grafo8639.75', 'graficon75nodi/grafo8663.75', 'graficon75nodi/grafo8906.75', 'graficon75nodi/grafo8939.75', 'graficon75nodi/grafo8961.75', 'graficon75nodi/grafo9024.75', 'graficon75nodi/grafo9050.75', 'graficon75nodi/grafo9105.75', 'graficon75nodi/grafo9147.75', 'graficon75nodi/grafo9161.75', 'graficon75nodi/grafo9188.75', 'graficon75nodi/grafo9226.75', 'graficon75nodi/grafo9281.75', 'graficon75nodi/grafo9290.75', 'graficon75nodi/grafo9353.75', 'graficon75nodi/grafo9414.75', 'graficon75nodi/grafo9502.75', 'graficon75nodi/grafo9505.75', 'graficon75nodi/grafo9545.75', 'graficon75nodi/grafo9645.75', 'graficon75nodi/grafo9735.75', 'graficon75nodi/grafo9736.75', 'graficon75nodi/grafo9768.75', 'graficon75nodi/grafo9782.75', 'graficon75nodi/grafo9802.75', 'graficon75nodi/grafo9845.75', 'graficon76nodi/grafo752.76', 'graficon76nodi/grafo2762.76', 'graficon76nodi/grafo4120.76', 'graficon76nodi/grafo4168.76', 'graficon76nodi/grafo4174.76', 'graficon76nodi/grafo4226.76', 'graficon76nodi/grafo4272.76', 'graficon76nodi/grafo4287.76', 'graficon76nodi/grafo4354.76', 'graficon76nodi/grafo4379.76', 'graficon76nodi/grafo4477.76', 'graficon76nodi/grafo4484.76', 'graficon76nodi/grafo4492.76', 'graficon76nodi/grafo4493.76', 'graficon76nodi/grafo4614.76', 'graficon76nodi/grafo4617.76', 'graficon76nodi/grafo4650.76', 'graficon76nodi/grafo4658.76', 'graficon76nodi/grafo4665.76', 'graficon76nodi/grafo4672.76', 'graficon76nodi/grafo4680.76', 'graficon76nodi/grafo4692.76', 'graficon76nodi/grafo4694.76', 'graficon76nodi/grafo4708.76', 'graficon76nodi/grafo4743.76', 'graficon76nodi/grafo4917.76', 'graficon76nodi/grafo4962.76', 'graficon76nodi/grafo5090.76', 'graficon76nodi/grafo5119.76', 'graficon76nodi/grafo5285.76', 'graficon76nodi/grafo5293.76', 'graficon76nodi/grafo5361.76', 'graficon76nodi/grafo5423.76', 'graficon76nodi/grafo5611.76', 'graficon76nodi/grafo5725.76', 'graficon76nodi/grafo5841.76', 'graficon76nodi/grafo5866.76', 'graficon76nodi/grafo5874.76', 'graficon76nodi/grafo5898.76', 'graficon76nodi/grafo6000.76', 'graficon76nodi/grafo6002.76', 'graficon76nodi/grafo6041.76', 'graficon76nodi/grafo6065.76', 'graficon76nodi/grafo6078.76', 'graficon76nodi/grafo6083.76', 'graficon76nodi/grafo6602.76', 'graficon76nodi/grafo6958.76', 'graficon76nodi/grafo7027.76', 'graficon76nodi/grafo7058.76', 'graficon76nodi/grafo7319.76', 'graficon76nodi/grafo7361.76', 'graficon76nodi/grafo7807.76', 'graficon76nodi/grafo7821.76', 'graficon76nodi/grafo7837.76', 'graficon76nodi/grafo7846.76', 'graficon76nodi/grafo7853.76', 'graficon76nodi/grafo7903.76', 'graficon76nodi/grafo7929.76', 'graficon76nodi/grafo7934.76', 'graficon76nodi/grafo7943.76', 'graficon76nodi/grafo7967.76', 'graficon76nodi/grafo7982.76', 'graficon76nodi/grafo8013.76', 'graficon76nodi/grafo8044.76', 'graficon76nodi/grafo8078.76', 'graficon76nodi/grafo8080.76', 'graficon76nodi/grafo8085.76', 'graficon76nodi/grafo8095.76', 'graficon76nodi/grafo8097.76', 'graficon76nodi/grafo8124.76', 'graficon76nodi/grafo8151.76', 'graficon76nodi/grafo8163.76', 'graficon76nodi/grafo8164.76', 'graficon76nodi/grafo8168.76', 'graficon76nodi/grafo8194.76', 'graficon76nodi/grafo8196.76', 'graficon76nodi/grafo8200.76', 'graficon76nodi/grafo8233.76', 'graficon76nodi/grafo8267.76', 'graficon76nodi/grafo8282.76', 'graficon76nodi/grafo8301.76', 'graficon76nodi/grafo8306.76', 'graficon76nodi/grafo8324.76', 'graficon76nodi/grafo8326.76', 'graficon76nodi/grafo8366.76', 'graficon76nodi/grafo8367.76', 'graficon76nodi/grafo8374.76', 'graficon76nodi/grafo8412.76', 'graficon76nodi/grafo8436.76', 'graficon76nodi/grafo8442.76', 'graficon76nodi/grafo8446.76', 'graficon76nodi/grafo8456.76', 'graficon76nodi/grafo8460.76', 'graficon76nodi/grafo8462.76', 'graficon76nodi/grafo8483.76', 'graficon76nodi/grafo8522.76', 'graficon76nodi/grafo8526.76', 'graficon76nodi/grafo8530.76', 'graficon76nodi/grafo8616.76', 'graficon76nodi/grafo9044.76', 'graficon76nodi/grafo9208.76', 'graficon76nodi/grafo9253.76', 'graficon76nodi/grafo9263.76', 'graficon76nodi/grafo9277.76', 'graficon76nodi/grafo9438.76', 'graficon76nodi/grafo9631.76', 'graficon76nodi/grafo9636.76', 'graficon76nodi/grafo9718.76', 'graficon76nodi/grafo9767.76', 'graficon76nodi/grafo9774.76', 'graficon76nodi/grafo9798.76', 'graficon76nodi/grafo9848.76']
    prefix = '../../Rome-Lib/'
    conditions = []
    csv_header = ['Index', 'Data'] + conditions + ['Tstandard', 'Twithconstraint']
    condition_map = {}

    if len(sys.argv) != 2:
        out_vars = []
        csvpath = './results/' + '_'.join(str(x) for x in out_vars) + '_results.csv'
        idx = get_start_position(csvpath, csv_header)
        ct = 1
        cur_conds = []
        print('Combination', ' + '.join([str(cond) for cond in cur_conds]))
        for file in files:
            if idx == ct:
                print(len(files) if idx % len(files) == 0 else idx % len(files), '/', len(files))
                out = run_func(data_path=prefix+file)
                insert_one(csvpath, [idx, file] + cur_conds + list(out))
                idx += 1
            ct += 1
    else:
        pid = int(sys.argv[1])

        out_vars = []
        csvpath = './results/' + '_'.join(str(x) for x in out_vars) + '_results.csv'
        idx = get_start_position(csvpath, csv_header)
        ct = 1
        start_pos = get_start_position(csvpath, csv_header)
        cur_conds = []
        print('Combination', ' + '.join([str(cond) for cond in cur_conds]))
        for file in files:
            if idx == ct:
                print(len(files) if idx % len(files) == 0 else idx % len(files), '/', len(files))
                out = run_func(data_path=prefix+file)
                insert_one(csvpath, [idx, file] + [] + list(out))
                idx += 1
            ct += 1
