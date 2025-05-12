import argparse
import networkx as nx
import powerlaw
import numpy as np
import json
from collections import defaultdict
from networkx.algorithms.community import modularity
import matplotlib.pyplot as plt


def membership_to_partition(membership):
    part_dict = {}
    for index, value in membership.items():
        if value in part_dict:
            part_dict[value].append(index)
        else:
            part_dict[value] = [index]
    return part_dict.values()


def get_membership_list_from_file(net, file_name):
    membership = dict()
    with open(file_name) as f:
        for line in f:
            i, m = line.strip().split()
            if int(i) in net.nodes:
                membership[int(i)] = m
    return membership


def write_membership_list_to_file(file_name, membership):
    with open(file_name, 'w') as f:
        f.write('\n'.join(str(i)+' '+str(membership[i]+1) for i in range(len(membership))))


def compute_mixing_param(net, membership):
    n = net.number_of_nodes()
    in_degree = defaultdict(int)
    out_degree = defaultdict(int)
    for n1, n2 in net.edges:
        if membership[n1] == membership[n2]: # nodes are co-clustered
            in_degree[n1] += 1
            in_degree[n2] += 1
        else:
            out_degree[n1] += 1
            out_degree[n2] += 1
    mus = [out_degree[i]/(out_degree[i]+in_degree[i]) if (out_degree[i]+in_degree[i]) != 0 else 0 for i in net.nodes]

    mixing_param = np.mean(mus)
    return mus



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Estimating properties of a network/clustering pair.')
    parser.add_argument('-n', metavar='net', type=str, required=True,
                        help='network edge-list path')
    parser.add_argument('-c', metavar='clustering', type=str, required=True,
                        help='clustering membership path')
    parser.add_argument('-o', metavar='output', type=str, required=True,
                        help='mixing parameter distribution')
    args = parser.parse_args()

    net = nx.read_edgelist(args.n, nodetype=int)
    membership = get_membership_list_from_file(net, args.c)

    mus=compute_mixing_param(net, membership)

    with open(args.o, "w") as f:
        for mu in mus:
            f.write(str(mu)+'\n')
