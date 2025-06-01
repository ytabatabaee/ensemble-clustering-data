import networkx as nx
from collections import defaultdict
import numpy as np
import argparse
from networkx.algorithms.community import modularity
from sklearn.metrics.cluster import adjusted_mutual_info_score
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.cluster import normalized_mutual_info_score
import csv


def partition_statistics(G, partition, show_cluster_size_dist=False):
    cluster_num = len(partition)
    cluster_sizes = [len(c) for c in partition]
    min_size, max_size, mean_size, median_size = np.min(cluster_sizes), np.max(cluster_sizes), np.mean(
        cluster_sizes), np.median(cluster_sizes)
    singletons = [c for c in partition if len(c) == 1]
    singletons_num = len(singletons)
    non_singleton_num = len(partition) - len(singletons)
    modularity_score = modularity(G, partition)
    coverage = (G.number_of_nodes() - len(singletons)) / G.number_of_nodes()

    print('#clusters in partition:', cluster_num)
    if show_cluster_size_dist:
        print(sorted(cluster_sizes, reverse=True))
    print('min, max, mean, median cluster sizes:', min_size, max_size, mean_size, median_size)
    print('number of singletons:', singletons_num)
    print('number of non-singleton clusters:', non_singleton_num)
    print('modularity:', modularity_score)
    print('coverage:', coverage)
    return cluster_sizes, coverage

    # return cluster_num, min_size, max_size, mean_size, median_size, singletons_num, non_singleton_num, modularity_score, coverage


def group_to_partition(partition):
    part_dict = {}
    for index, value in partition.items():
        if value in part_dict:
            part_dict[value].append(index)
        else:
            part_dict[value] = [index]
    return part_dict.values()


def membership_list_to_dict(membership_list):
    membership_dict = {}
    for i in range(len(membership_list)):
        membership_dict[i] = membership_list[i]
    return membership_dict


def get_membership_list_from_dict(membership_dict, n):
    memberships = [0]*n
    for i in range(len(membership_dict)):
        for x in membership_dict[i]:
            memberships[x] = i
    return memberships


def communities_to_dict(communities):
    result = {}
    community_index = 0
    for c in communities:
        community_mapping = ({node: community_index for index, node in enumerate(c)})
        result = {**result, **community_mapping}
        community_index += 1
    return result


def measure_accuracy(mem_true, mem_est, community_strength = None, threshold=1.0):
    weak_communities = []
    for i in range(len(community_strength)): # assuming memberships are consecutive integers starting from 0
        if community_strength[i] < threshold:
            weak_communities.append(i)

    mem_true_new, mem_est_new = dict(), dict()
    #print(mem_true)
    for node in mem_true.keys():
        if mem_true[node] in weak_communities: # only pick the nodes that are in clusters that belong to the weak communities
            mem_true_new[node] = mem_true[node]
            mem_est_new[node] = mem_est[node]

    mem_true_new = list(mem_true_new.values())
    mem_est_new = list(mem_est_new.values())
    n = len(mem_true_new)

    nmi = normalized_mutual_info_score(mem_true_new, mem_est_new)
    ari = adjusted_rand_score(mem_true_new, mem_est_new)
    ami = adjusted_mutual_info_score(mem_true_new, mem_est_new)

    tn, tp, fn, fp = 0, 0, 0, 0
    for i in range(n):
        for j in range(i + 1, n):
            if mem_true_new[i] == mem_true_new[j]:
                if mem_est_new[i] == mem_est_new[j]:
                    tp += 1
                else:
                    fn += 1
            else:
                if mem_est_new[i] == mem_est_new[j]:
                    fp += 1
                else:
                    tn += 1
    precision = tp / (tp + fp)
    recall = tp / (tp + fn)
    f1_score = 2 * precision * recall / (precision + recall)
    fnr = fn / (fn + tp)
    fpr = fp / (fp + tn)

    #return nmi, ami, ari, 0, 0, 0, 0, 0
    return nmi, ami, ari, precision, recall, f1_score, fnr, fpr


def read_network_partition(net_path, gt_path):
    graph = nx.read_edgelist(net_path, nodetype=int)
    membership = [''] * graph.number_of_nodes()
    with open(gt_path) as fgt:
        for line in fgt:
            i, m = line.strip().split()
            membership[int(i)] = int(m)
    #n = graph.number_of_nodes()
    '''for i in range(len(membership)):
        if membership[i] == '':
            membership[i] = n + i'''
    return graph, membership


def get_membership_list_from_file(net, file_name):
    membership = dict()
    with open(file_name) as f:
        for line in f:
            i, m = line.strip().split()
            if int(i) in net.nodes:
                membership[int(i)] = int(m)
    return membership


def write_membership_to_file(file_path, membership):
    with open(file_path, 'w') as out_file:
        writer = csv.writer(out_file, delimiter=' ')
        for i in range(len(membership)):
            writer.writerow([i] + [membership[i]])


def calculate_stats(graph, partition, n, k, res, method, partition_name, ground_truth_membership, acc_df, df):
    cluster_sizes, node_coverage = partition_statistics(graph, partition)
    for i in range(len(cluster_sizes)):
        df.loc[len(df.index)] = [n, k, res, method, partition_name, i, cluster_sizes[i]]
    membership = get_membership_list_from_dict(partition, n * k)
    nmi, ami, ari, precision, recall, f1_score, fnr, fpr = measure_accuracy(ground_truth_membership, membership)
    print("Normalized mutual information (NMI): ", nmi)
    print("Adjusted rand index (ARI): ", ari)
    print("Adjusted mutual information (AMI): ", ami)
    print("False positive rate (FPR), False negative rate (FNR):", fpr, fnr)
    print("Precision, Recall, F1-score:", precision, recall, f1_score)
    acc_df.loc[len(acc_df.index)] = [n, k, res, method, partition_name, 'NMI', nmi]
    acc_df.loc[len(acc_df.index)] = [n, k, res, method, partition_name, 'AMI', ami]
    acc_df.loc[len(acc_df.index)] = [n, k, res, method, partition_name, 'ARI', ari]
    acc_df.loc[len(acc_df.index)] = [n, k, res, method, partition_name, 'FPR', fpr]
    acc_df.loc[len(acc_df.index)] = [n, k, res, method, partition_name, 'FNR', fnr]
    acc_df.loc[len(acc_df.index)] = [n, k, res, method, partition_name, 'Precision', precision]
    acc_df.loc[len(acc_df.index)] = [n, k, res, method, partition_name, 'Recall', recall]
    acc_df.loc[len(acc_df.index)] = [n, k, res, method, partition_name, 'F1-Score', f1_score]
    return df, acc_df, membership


# calculates community strength (ratio of external to total degree) according to equation (9.2) in Barabasi (2016)
# communities where this ratio is below the threshold t=0.5 are considered weak
def calculate_community_strength(net, membership):
    n = net.number_of_nodes()
    in_degree = defaultdict(int)
    ex_degree = defaultdict(int)

    for n1, n2 in net.edges:
        if membership[n1] == membership[n2]:  # nodes are co-clustered
            in_degree[membership[n1]] += 2
        else:
            ex_degree[membership[n1]] += 1
            ex_degree[membership[n2]] += 1

    community_strength = [ex_degree[c] / (ex_degree[c] + in_degree[c]) for c in range(len(in_degree))]
    return community_strength


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Estimating properties of a network/clustering pair.')
    parser.add_argument('-n', metavar='net', type=str, required=True,
                        help='network edge-list path')
    parser.add_argument('-gt', metavar='groundtruth', type=str, required=True,
                        help='ground-truth membership path')
    parser.add_argument('-c', metavar='clustering', type=str, required=True,
                        help='clustering membership path')
    parser.add_argument('-t', metavar='threshold', type=float, required=False,
                        help='weak community threshold', default=0.5)
    args = parser.parse_args()
    net = nx.read_edgelist(args.n, nodetype=int)
    groundtruth = get_membership_list_from_file(net, args.gt)
    membership = get_membership_list_from_file(net, args.c)
    community_strength = calculate_community_strength(net, groundtruth)

    nmi, ami, ari, precision, recall, f1_score, fnr, fpr = measure_accuracy(groundtruth, membership, community_strength=community_strength, threshold=args.t)
    print("Normalized mutual information (NMI): ", nmi)
    print("Adjusted rand index (ARI): ", ari)
    print("Adjusted mutual information (AMI): ", ami)
    print("False positive rate (FPR), False negative rate (FNR):", fpr, fnr)
    print("Precision, Recall, F1-score:", precision, recall, f1_score)

