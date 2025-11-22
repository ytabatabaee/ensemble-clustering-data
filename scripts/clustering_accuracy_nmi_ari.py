from sklearn.metrics.cluster import normalized_mutual_info_score, adjusted_mutual_info_score, adjusted_rand_score
import numpy as np
import argparse


def membership_to_partition(membership):
    part_dict = {}
    for value in membership:
        if value in part_dict:
            part_dict[value] += 1
        else:
            part_dict[value] = 1
    return list(part_dict.values())


def get_membership_list_shared_nodes(gt_path, f1_path):
    gt_membership = dict()
    membership1 = dict()
    with open(gt_path) as fgt:
        for line in fgt:
            i, m = line.strip().split()
            gt_membership[int(i)] = m
    with open(f1_path) as f1:
        for line in f1:
            i, m = line.strip().split()
            membership1[int(i)] = m
    print('#nodes in ground-truth partition:', len(gt_membership.keys()))
    print('#nodes in estimated partition:', len(membership1.keys()))
    keys = list(set(membership1.keys()) & set(gt_membership.keys()))
    keys.sort()
    print('common nodes between partitions:', len(keys))
    mem1 = {i: membership1[i] for i in keys}
    memgt = {i: gt_membership[i] for i in keys}
    return list(memgt.values()), list(mem1.values())


def get_membership_list_add_singletons(gt_path, f1_path):
    gt_membership = dict()
    membership1 = dict()
    with open(gt_path) as fgt:
        for line in fgt:
            i, m = line.strip().split()
            gt_membership[int(i)] = m
    with open(f1_path) as f1:
        for line in f1:
            i, m = line.strip().split()
            membership1[int(i)] = m
    print('#nodes in ground-truth partition:', len(gt_membership.keys()))
    print('#nodes in estimated partition:', len(membership1.keys()))
    keys = list(gt_membership.keys())
    keys.sort()
    mem_gt = {i: gt_membership[i] for i in keys}
    mem1 = dict()
    for i in keys:
        if i in membership1:
            mem1[i] = membership1[i]
        else:
            mem1[i] = i
    # print('#singletons added to estimated clustering:', len(set(mem1.values())) - len(set(membership1.values())))
    return list(mem_gt.values()), list(mem1.values())


def measure_accuracy(mem_true, mem_est):
    nmi = normalized_mutual_info_score(mem_true, mem_est)
    ari = adjusted_rand_score(mem_true, mem_est)
    ami = adjusted_mutual_info_score(mem_true, mem_est)

    return nmi, ari, ami


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="LFR accuracy for pre-CM vs post-CM partitions")
    parser.add_argument("-gt", "--groundtruth", type=str, required=True,
                        help="File containing ground-truth community membership")
    parser.add_argument("-p", "--partition", type=str, required=True,
                        help="File containing estimated community membership")
    args = parser.parse_args()
    gt_membership, est_membership = get_membership_list_add_singletons(args.groundtruth, args.partition)

    nmi, ari, ami = measure_accuracy(gt_membership, est_membership)

    print('\nAccuracy:')
    print("Normalized mutual information (NMI): ", nmi)
    print("Adjusted mutual information (AMI): ", ami)
    print("Adjusted rand index (ARI): ", ari)

