import multiprocessing as mp
import subprocess
import sys
import os

dataset_path = '/scratch/users/syt3/clustering/CM_paper_LFR_data/'
software_path = '/scratch/users/syt3/clustering/res_limit_exps.py'


def run_consensus(comb):
    net, res = comb
    tr = 0.2 # threshold
    data_path = dataset_path + net + '_' + res + '_lfr/'
    network_path = data_path + 'network.dat'
    gt_path = data_path + 'community.dat'
    output_path = data_path + 'sc_' + str(tr) + '.dat'
    log_path = data_path + 'sc_' + str(tr) + '.log'
    if res == 'modularity':
        cmd = 'python3 ' + software_path + ' -n ' + network_path + ' -gt ' + gt_path + ' --relabel -t ' + str(tr) + \
              ' --algorithm leiden-mod -p 10 -o ' + output_path + ' > ' + log_path
    else:
        cmd = 'python3 ' + software_path + ' -n ' + network_path + ' -gt ' + gt_path + ' --relabel -t ' + str(tr) +\
              ' --algorithm leiden-cpm -p 10 --resolution ' + res + ' -o ' + output_path + ' > ' + log_path
    full_cmd = '/usr/bin/time -v -o ' + output_path + '.stat' + ' -f "QR*\t%e\t%M" ' + cmd
    p = subprocess.Popen(full_cmd, stdout=subprocess.PIPE, shell=True)
    _, _ = p.communicate()


if __name__ == '__main__':
    combinations = []
    networks = ['cit_patents', 'cen', 'oc', 'cit_hepph', 'wiki_talk', 'wiki_topcats']
    res_values = ['modularity']#['0.0001', '0.001', '0.01', '0.1', '0.5']
    for n in networks:
        for r in res_values:
            combinations.append((n, r))
    with mp.Pool(mp.cpu_count()) as p:
        p.map(run_consensus, combinations)

