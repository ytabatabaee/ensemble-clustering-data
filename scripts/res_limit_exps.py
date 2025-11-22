import leidenalg
import networkx as nx
import igraph as ig
import community.community_louvain as cl
import numpy as np
import pandas as pd
import argparse
import os
import sys
import seaborn as sns
from networkx.algorithms.community import modularity
from sklearn.metrics.cluster import adjusted_mutual_info_score
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.cluster import normalized_mutual_info_score
import matplotlib.pyplot as plt
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


def initialize(graph, value):
    for u, v in graph.edges():
        graph[u][v]['weight'] = value
    return graph

def community_ecg(self, weights=None, ens_size=16, min_weight=0.05):
    W = [0]*self.ecount()
    ## Ensemble of level-1 Louvain
    for i in range(ens_size):
        p = np.random.permutation(self.vcount()).tolist()
        g = self.permute_vertices(p)
        l = g.community_multilevel(weights=weights, return_levels=True)[0].membership
        b = [l[p[x.tuple[0]]]==l[p[x.tuple[1]]] for x in self.es]
        W = [W[i]+b[i] for i in range(len(W))]
    W = [min_weight + (1-min_weight)*W[i]/ens_size for i in range(len(W))]
    part = self.community_multilevel(weights=W)
    part._modularity_params['weights'] = weights
    part.recalculate_modularity()
    ## Force min_weight outside 2-core
    core = self.shell_index()
    ecore = [min(core[x.tuple[0]],core[x.tuple[1]]) for x in self.es]
    part.W = [W[i] if ecore[i]>1 else min_weight for i in range(len(ecore))]
    part.CSI = 1-2*np.sum([min(1-i,i) for i in part.W])/len(part.W)
    return part
ig.Graph.community_ecg = community_ecg


def get_communities(graph, algorithm, seed, res_val=0.001):
    if algorithm == 'louvain':
        return cl.best_partition(graph, random_state=seed, weight='weight')
    elif algorithm == 'leiden-cpm':
        return communities_to_dict(leidenalg.find_partition(ig.Graph.from_networkx(graph),
                                                  leidenalg.CPMVertexPartition,
                                                  resolution_parameter=res_val,
                                                  weights='weight',
                                                  n_iterations=2,
                                                  seed=seed).as_cover())
    elif algorithm == 'leiden-mod':
        return communities_to_dict(leidenalg.find_partition(ig.Graph.from_networkx(graph),
                                        leidenalg.ModularityVertexPartition,
                                        weights='weight',
                                        seed=seed).as_cover())

def check_convergence(G, n_p, delta):
    count = 0
    for wt in nx.get_edge_attributes(G, 'weight').values():
        if wt != 0 and wt != n_p:
            count += 1
    if count > delta * G.number_of_edges():
        return False
    return True

def thresholding(graph, n_p, thresh):
    remove_edges = []
    for u, v in graph.edges():
        if graph[u][v]['weight'] < thresh * n_p:
            remove_edges.append((u, v))
    graph.remove_edges_from(remove_edges)
    return graph


def simple_consensus(G, algorithm='leiden-mod', n_p=10, thresh=0.8, delta=0.02, max_iter=10, res_value=0.01):
    graph = G.copy()
    graph = initialize(graph, 1.0)
    iter_count = 0

    while True:
        iter_count += 1
        if iter_count > max_iter:
            iter_count -= 1
            break
        nextgraph = graph.copy()
        nextgraph = initialize(nextgraph, 0.0)
        partitions = [get_communities(graph, algorithm, i, res_val=res_value) for i in range(n_p)]

        # print('edges', len(graph.edges()))
        for i in range(n_p):
            # print('np: ', i)
            c = partitions[i]
            for node, nbr in graph.edges():
                if graph[node][nbr]['weight'] not in (0, n_p):
                    if c[node] == c[nbr]:
                        nextgraph[node][nbr]['weight'] += 1
                else:
                    nextgraph[node][nbr]['weight'] = graph[node][nbr]['weight']
                # print(node, nbr, nextgraph[node][nbr]['weight'])

        nextgraph = thresholding(nextgraph, n_p, thresh)
        if check_convergence(nextgraph, n_p, delta=delta):
            break
        graph = nextgraph.copy()

    print('number of iterations:', iter_count)
    return group_to_partition(get_communities(graph, algorithm, 0, res_val=res_value))


def strict_consensus(G, algorithm='leiden-cpm', n_p=20, res_val=0.01):
    graph = G.copy()
    graph = initialize(graph, 1.0)
    iter_count = 0

    partitions = [get_communities(graph, algorithm, i, res_val) for i in range(n_p)]

    for i in range(n_p):
        c = partitions[i]
        for node, nbr in graph.edges():
            if c[node] != c[nbr]:
                graph[node][nbr]['weight'] = 0

    graph = thresholding(graph, 1, 1)

    return group_to_partition(get_communities(graph, algorithm, 0, res_val))


def fast_ensemble(G, algorithm='leiden-mod', n_p=10, tr=0.8, res_value=0.01):
    graph = G.copy()
    graph = initialize(graph, 1)
    partitions = [get_communities(graph, algorithm, i, res_val=res_value) for i in range(n_p)]

    remove_edges = []
    for node, nbr in graph.edges():
        for i in range(n_p):
            if partitions[i][node] != partitions[i][nbr]:
                graph[node][nbr]['weight'] -= 1/n_p
            if graph[node][nbr]['weight'] < tr:
                remove_edges.append((node, nbr))
                break
    graph.remove_edges_from(remove_edges)

    return group_to_partition(get_communities(graph, algorithm, 0, res_val=res_value))



def normal_clustering(graph, algorithm, res_val=0.01):
    if algorithm == 'leiden-cpm':
        return leidenalg.find_partition(ig.Graph.from_networkx(graph),
                                              leidenalg.CPMVertexPartition,
                                              resolution_parameter=res_val,
                                              n_iterations=2).as_cover()
    elif algorithm == 'leiden-mod':
        return leidenalg.find_partition(ig.Graph.from_networkx(graph),
                                         leidenalg.ModularityVertexPartition,
                                         seed=1234).as_cover()
    elif algorithm == 'louvain':
        #return cl.best_partition(graph)
        return list(group_to_partition(cl.best_partition(graph)))


def gen_tree_of_cliques(k, n):
  '''
  k: size of clique
  n: number of cliques
  '''
  cliques = [nx.complete_graph(k) for _ in range(n)]
  tree = nx.random_tree(n)
  tree_of_cliques = nx.disjoint_union_all(cliques)
  for s, d in tree.edges():
    tree_of_cliques.add_edge(s*k+4, d*k)
  return tree_of_cliques


def measure_accuracy(mem_true, mem_est):
    n = len(mem_true)
    tn, tp, fn, fp = 0, 0, 0, 0
    for i in range(n):
        for j in range(i + 1, n):
            if mem_true[i] == mem_true[j]:
                if mem_est[i] == mem_est[j]:
                    tp += 1
                else:
                    fn += 1
            else:
                if mem_est[i] == mem_est[j]:
                    fp += 1
                else:
                    tn += 1
    precision = tp / (tp + fp)
    recall = tp / (tp + fn)
    f1_score = 2 * precision * recall / (precision + recall)
    fnr = fn / (fn + tp)
    fpr = fp / (fp + tn)

    nmi = normalized_mutual_info_score(mem_true, mem_est)
    ari = adjusted_rand_score(mem_true, mem_est)
    ami = adjusted_mutual_info_score(mem_true, mem_est)

    #return nmi, ami, ari, 0, 0, 0, 0, 0
    return nmi, ami, ari, precision, recall, f1_score, fnr, fpr


def read_network_partition(net_path, gt_path):
    graph = nx.read_edgelist(net_path, nodetype=int)
    #print(graph.edges)

    # for E-R experiments
    '''
    for i in range(2000):
        if not graph.has_node(i):
            graph.add_node(i)
    print(graph.number_of_nodes())
    '''

    membership = [''] * graph.number_of_nodes()
    with open(gt_path) as fgt:
        for line in fgt:
            i, m = line.strip().split()
            membership[int(i)] = int(m)
    n = graph.number_of_nodes()
    for i in range(len(membership)):
        if membership[i] == '':
            membership[i] = n + i
    #print(membership)
    return graph, membership


def write_membership_to_file(file_path, memberhsip):
    with open(file_path, 'w') as out_file:
        writer = csv.writer(out_file, delimiter=' ')
        for i in range(len(memberhsip)):
            writer.writerow([i] + [memberhsip[i]])


def lfr_exp_mu():
    acc_df = pd.DataFrame(columns=["n", "mu", "partition", "acc_measure", "acc_value"])
    method = 'leiden-mod'
    n = 10000
    for mu in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
         ### Generating the graph
        graph = nx.generators.community.LFR_benchmark_graph(n=10000, tau1=3, tau2=1.5, mu=mu, average_degree=10,
                                                            min_community=10, seed=1932)
        if not os.path.exists('lfr_training_new/mu_'+str(mu)):
            os.mkdir('lfr_training_new/mu_'+str(mu))
        nx.write_edgelist(graph, 'lfr_training_new/mu_'+str(mu) + '/network.dat', data=False)

        communities = {frozenset(graph.nodes[v]["community"]) for v in graph}
        print(len(communities))
        cluster_sizes = [len(c) for c in communities]
        print(min(cluster_sizes), max(cluster_sizes))
        communities = communities_to_dict(communities)
        keys = list(communities.keys())
        keys.sort()
        membership = {i: communities[i] for i in keys}
        ground_truth_membership = list(membership.values())

        write_membership_to_file('lfr_training_new/mu_'+str(mu) + '/community.dat', ground_truth_membership)

        #graph, ground_truth_membership = read_network_ground_truth('lfr_training/mu_'+str(mu) + '/network.dat', 'lfr_training/mu_'+str(mu) + '/community.dat')

        print('Leiden-MOD')
        partition = normal_clustering(graph, method)
        acc_df, membership = calculate_stats_mu(list(partition), n, mu, 'Leiden-MOD', ground_truth_membership, acc_df)
        write_membership_to_file('lfr_training_new/mu_' + str(mu) + '/leiden_mod.dat', membership)

        print('ECG')
        g = ig.Graph.from_networkx(graph)
        ec = g.community_ecg(ens_size=10)
        acc_df, membership = calculate_stats_mu(list(ec), n, mu, 'ECG', ground_truth_membership, acc_df)
        write_membership_to_file('lfr_training_new/mu_' + str(mu) + '/ecg.dat', membership)


        '''partition = fast_ensemble(graph, method, n_p=10, tr=0.2)
        acc_df, membership = calculate_stats_mu(list(partition), n, mu, 'FastEnsemble(tr=0.2,Leiden-MOD)', ground_truth_membership, acc_df)
        write_membership_to_file('lfr_training_new/mu_' + str(mu) + '/fec_leiden_mod_0.2.dat', membership)

        partition = fast_ensemble(graph, method, n_p=10, tr=0.5)
        acc_df, membership = calculate_stats_mu(list(partition), n, mu, 'FastEnsemble(tr=0.5,Leiden-MOD)',
                                                ground_truth_membership, acc_df)
        write_membership_to_file('lfr_training_new/mu_' + str(mu) + '/fec_leiden_mod_0.5.dat', membership)'''

        print('FastEnsemble(tr=0.8,Leiden-MOD)')
        partition = fast_ensemble(graph, method, n_p=10, tr=0.8)
        acc_df, membership = calculate_stats_mu(list(partition), n, mu, 'FastEnsemble(tr=0.8,Leiden-MOD)',
                                                ground_truth_membership, acc_df)

        print('FastEnsemble(tr=0.9,Leiden-MOD)')
        partition = fast_ensemble(graph, method, n_p=10, tr=0.9)
        acc_df, membership = calculate_stats_mu(list(partition), n, mu, 'FastEnsemble(tr=0.9,Leiden-MOD)',
                                                 ground_truth_membership, acc_df)
        write_membership_to_file('lfr_training_new/mu_' + str(mu) + '/fec_leiden_mod_0.9.dat', membership)

    acc_df.to_csv('training_exp_tandon.csv')


def lfr_exp_deg():
    acc_df = pd.DataFrame(columns=["deg", "mu", "partition", "acc_measure", "acc_value"])
    method = 'leiden-mod'
    n = 10000
    data_dir = 'lfr_training_deg'
    for d in [5, 10, 20]:
        if not os.path.exists(data_dir+'/d_'+str(d)):
            os.mkdir(data_dir+'/d_'+str(d))
        for mu in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
            ### Generating the graph
            dir=data_dir+'/d_'+str(d)+'/mu_'+str(mu)
            graph = nx.generators.community.LFR_benchmark_graph(n=10000, tau1=3, tau2=1.5, mu=mu, average_degree=d,
                                                                min_community=10, seed=1932)
            if not os.path.exists(dir):
                os.mkdir(dir)
            nx.write_edgelist(graph, dir + '/network.dat', data=False)

            communities = {frozenset(graph.nodes[v]["community"]) for v in graph}
            print(len(communities))
            cluster_sizes = [len(c) for c in communities]
            print(min(cluster_sizes), max(cluster_sizes))
            communities = communities_to_dict(communities)
            keys = list(communities.keys())
            keys.sort()
            membership = {i: communities[i] for i in keys}
            ground_truth_membership = list(membership.values())

            write_membership_to_file(dir + '/community.dat', ground_truth_membership)

            #graph, ground_truth_membership = read_network_ground_truth('lfr_training/mu_'+str(mu) + '/network.dat', 'lfr_training/mu_'+str(mu) + '/community.dat')

            print('Leiden-MOD')
            partition = normal_clustering(graph, method)
            acc_df, membership = calculate_stats_deg(list(partition), n, d, mu, 'Leiden-MOD', ground_truth_membership, acc_df)
            write_membership_to_file(dir + '/leiden_mod.dat', membership)

            print('ECG')
            g = ig.Graph.from_networkx(graph)
            ec = g.community_ecg(ens_size=10)
            acc_df, membership = calculate_stats_deg(list(ec), n, d, mu, 'ECG', ground_truth_membership, acc_df)
            write_membership_to_file(dir + '/ecg.dat', membership)

            for tr in [0.2, 0.5, 0.8, 0.9]:
                partition = fast_ensemble(graph, method, n_p=10, tr=tr)
                acc_df, membership = calculate_stats_deg(list(partition), n, d, mu, 'FastEnsemble(tr='+str(tr)+',Leiden-MOD)', ground_truth_membership, acc_df)
                write_membership_to_file(dir + '/fec_leiden_mod_'+str(tr)+'.dat', membership)

    acc_df.to_csv('training_exp_degree.csv')


def lfr_exp_size():
    acc_df = pd.DataFrame(columns=["n", "mu", "partition", "acc_measure", "acc_value"])
    method = 'leiden-mod'
    data_dir = 'lfr_training_size'
    d = 10
    for n in [1000, 10000, 100000]:
        if not os.path.exists(data_dir+'/n_'+str(n)):
            os.mkdir(data_dir+'/n_'+str(n))
        for mu in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
            ### Generating the graph
            dir=data_dir+'/n_'+str(n)+'/mu_'+str(mu)
            graph = nx.generators.community.LFR_benchmark_graph(n=n, tau1=3, tau2=1.5, mu=mu, average_degree=d,
                                                                min_community=10, seed=1932)
            if not os.path.exists(dir):
                os.mkdir(dir)
            nx.write_edgelist(graph, dir + '/network.dat', data=False)

            communities = {frozenset(graph.nodes[v]["community"]) for v in graph}
            print(len(communities))
            cluster_sizes = [len(c) for c in communities]
            print(min(cluster_sizes), max(cluster_sizes))
            communities = communities_to_dict(communities)
            keys = list(communities.keys())
            keys.sort()
            membership = {i: communities[i] for i in keys}
            ground_truth_membership = list(membership.values())

            write_membership_to_file(dir + '/community.dat', ground_truth_membership)

            #graph, ground_truth_membership = read_network_ground_truth('lfr_training/mu_'+str(mu) + '/network.dat', 'lfr_training/mu_'+str(mu) + '/community.dat')

            print('Leiden-MOD')
            partition = normal_clustering(graph, method)
            acc_df, membership = calculate_stats_mu(list(partition), n, mu, 'Leiden-MOD', ground_truth_membership, acc_df)
            write_membership_to_file(dir + '/leiden_mod.dat', membership)

            print('ECG')
            g = ig.Graph.from_networkx(graph)
            ec = g.community_ecg(ens_size=10)
            acc_df, membership = calculate_stats_mu(list(ec), n, mu, 'ECG', ground_truth_membership, acc_df)
            write_membership_to_file(dir + '/ecg.dat', membership)

            for tr in [0.2, 0.5, 0.8, 0.9]:
                partition = fast_ensemble(graph, method, n_p=10, tr=tr)
                acc_df, membership = calculate_stats_mu(list(partition), n, mu, 'FastEnsemble(tr='+str(tr)+',Leiden-MOD)', ground_truth_membership, acc_df)
                write_membership_to_file(dir + '/fec_leiden_mod_'+str(tr)+'.dat', membership)

    acc_df.to_csv('training_exp_size.csv')


def lfr_exp_tau():
    mu = 0.5
    graph = nx.generators.community.LFR_benchmark_graph(n=10000, tau1=3, tau2=1.5, mu=mu, average_degree=10,
                                                        min_community=10, seed=1932)
    communities = {frozenset(graph.nodes[v]["community"]) for v in graph}
    print(len(communities))
    cluster_sizes = [len(c) for c in communities]
    print(min(cluster_sizes), max(cluster_sizes))
    communities = communities_to_dict(communities)
    keys = list(communities.keys())
    keys.sort()
    membership = {i: communities[i] for i in keys}
    ground_truth_membership = list(membership.values())

    acc_df = pd.DataFrame(columns=["mu", "tr", "partition", "acc_measure", "acc_value"])
    method = 'leiden-mod'

    for tr in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]:
        print('tr:', tr)
        partition = normal_clustering(graph, method)
        membership = get_membership_list_from_dict(partition, 10000)
        nmi, ami, ari, precision, recall, f1_score, fnr, fpr = measure_accuracy(ground_truth_membership,
                                                                                membership)
        print("Normalized mutual information (NMI): ", nmi)
        print("Adjusted rand index (ARI): ", ari)
        print("Adjusted mutual information (AMI): ", ami)
        print("False positive rate (FPR), False negative rate (FNR):", fpr, fnr)
        print("Precision, Recall, F1-score:", precision, recall, f1_score)
        acc_df.loc[len(acc_df.index)] = [mu, tr, 'Leiden-MOD', 'NMI', nmi]
        acc_df.loc[len(acc_df.index)] = [mu, tr, 'Leiden-MOD', 'AMI', ami]
        acc_df.loc[len(acc_df.index)] = [mu, tr, 'Leiden-MOD', 'ARI', ari]
        acc_df.loc[len(acc_df.index)] = [mu, tr, 'Leiden-MOD', 'FPR', fpr]
        acc_df.loc[len(acc_df.index)] = [mu, tr, 'Leiden-MOD', 'FNR', fnr]
        acc_df.loc[len(acc_df.index)] = [mu, tr, 'Leiden-MOD', 'Precision', precision]
        acc_df.loc[len(acc_df.index)] = [mu, tr, 'Leiden-MOD', 'Recall', recall]
        acc_df.loc[len(acc_df.index)] = [mu, tr, 'Leiden-MOD', 'F1-Score', f1_score]

        partition = fast_ensemble(graph, method, n_p=10, tr=tr)
        partition = communities_to_dict(partition)
        keys = list(partition.keys())
        keys.sort()
        membership = {i: partition[i] for i in keys}
        membership = list(membership.values())
        nmi, ami, ari, precision, recall, f1_score, fnr, fpr = measure_accuracy(ground_truth_membership,
                                                                                membership)
        print("Normalized mutual information (NMI): ", nmi)
        print("Adjusted rand index (ARI): ", ari)
        print("Adjusted mutual information (AMI): ", ami)
        print("False positive rate (FPR), False negative rate (FNR):", fpr, fnr)
        print("Precision, Recall, F1-score:", precision, recall, f1_score)
        acc_df.loc[len(acc_df.index)] = [mu, tr, 'TC(np=10)+Leiden-MOD', 'NMI', nmi]
        acc_df.loc[len(acc_df.index)] = [mu, tr, 'TC(np=10)+Leiden-MOD', 'AMI', ami]
        acc_df.loc[len(acc_df.index)] = [mu, tr, 'TC(np=10)+Leiden-MOD', 'ARI', ari]
        acc_df.loc[len(acc_df.index)] = [mu, tr, 'TC(np=10)+Leiden-MOD', 'FPR', fpr]
        acc_df.loc[len(acc_df.index)] = [mu, tr, 'TC(np=10)+Leiden-MOD', 'FNR', fnr]
        acc_df.loc[len(acc_df.index)] = [mu, tr, 'TC(np=10)+Leiden-MOD', 'Precision', precision]
        acc_df.loc[len(acc_df.index)] = [mu, tr, 'TC(np=10)+Leiden-MOD', 'Recall', recall]
        acc_df.loc[len(acc_df.index)] = [mu, tr, 'TC(np=10)+Leiden-MOD', 'F1-Score', f1_score]

    acc_df.to_csv('res_limit_exps_leiden_mod_lfr_tr_acc.csv')


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


def calculate_stats_mu(partition, n, mu, partition_name, ground_truth_membership, acc_df):
    membership = get_membership_list_from_dict(partition, n)
    nmi, ami, ari, precision, recall, f1_score, fnr, fpr = measure_accuracy(ground_truth_membership, membership)
    print("Normalized mutual information (NMI): ", nmi)
    print("Adjusted rand index (ARI): ", ari)
    acc_df.loc[len(acc_df.index)] = [n, mu, partition_name, 'NMI', nmi]
    acc_df.loc[len(acc_df.index)] = [n, mu, partition_name, 'ARI', ari]
    return acc_df, membership

def calculate_stats_deg(partition, n, d, mu, partition_name, ground_truth_membership, acc_df):
    membership = get_membership_list_from_dict(partition, n)
    nmi, ami, ari, precision, recall, f1_score, fnr, fpr = measure_accuracy(ground_truth_membership, membership)
    print("Normalized mutual information (NMI): ", nmi)
    print("Adjusted rand index (ARI): ", ari)
    acc_df.loc[len(acc_df.index)] = [d, mu, partition_name, 'NMI', nmi]
    acc_df.loc[len(acc_df.index)] = [d, mu, partition_name, 'ARI', ari]
    return acc_df, membership


def ring_exp_cpm():
    df = pd.DataFrame(columns=['k', "n", "res", "method", "partition", "cluster_id", "cluster_size"])
    acc_df = pd.DataFrame(columns=['k', "n", "res", "method", "partition", "acc_measure", "acc_value"])
    n = 1000
    graph = nx.ring_of_cliques(num_cliques=n, clique_size=10)
    if not os.path.exists('ring_cpm_res'):
        os.mkdir('ring_cpm_res')
    nx.write_edgelist(graph, 'ring_cpm_res/network.dat', data=False)
    ground_truth_membership = sum([[i] * 10 for i in range(n)], [])
    write_membership_to_file('ring_cpm_res/community.dat', ground_truth_membership)
    # graph = gen_tree_of_cliques(k, n)
    for res in [0.00001, 0.0001, 0.001, 0.01, 0.1, 0.5]:
    #for n in [90, 100, 500, 1000, 5000]:#, 10000]:
        print('n:', n)
        for k in [10]:
            for method in ['leiden-cpm']:
                print('original')
                partition = normal_clustering(graph, method, res_val=res)
                df, acc_df, membership = calculate_stats(graph, partition, n, k, res, method,
                                                             'Leiden-CPM', ground_truth_membership, acc_df, df)
                write_membership_to_file('ring_cpm_res/leiden_cpm_r'+str(res)+'.dat', membership)

                print('TC n = 10')
                partition = fast_ensemble(graph, method, n_p=10, res_value=res)
                df, acc_df, membership = calculate_stats(graph, list(partition), n, k, res, method,
                                             'FastEnsemble(Leiden-CPM)', ground_truth_membership, acc_df, df)
                write_membership_to_file('ring_cpm_res/fec_leiden_cpm_r' + str(res) + '.dat', membership)

                print('TC n = 10')
                partition = fast_ensemble(graph, method, n_p=10, res_value=res, tr=0.9)
                df, acc_df, membership = calculate_stats(graph, list(partition), n, k, res, method,
                                             'FastEnsemble(Leiden-CPM,tr=0.9)', ground_truth_membership, acc_df, df)
                write_membership_to_file('ring_cpm_res/fec_leiden_cpm_r' + str(res) + '.dat', membership)

                print('SC n = 10')
                partition = strict_consensus(graph, method, n_p=10, res_val=res)
                df, acc_df, membership = calculate_stats(graph, list(partition), n, k, res, method,
                                             'Strict(np=10,Leiden-CPM)', ground_truth_membership, acc_df, df)
                write_membership_to_file('ring_cpm_res/strict_np10_leiden_cpm_r' + str(res) + '.dat', membership)

                print('SC n = 50')
                partition = strict_consensus(graph, method, n_p=50, res_val=res)
                df, acc_df, membership = calculate_stats(graph, list(partition), n, k, res, method,
                                             'Strict(np=50,Leiden-CPM)', ground_truth_membership, acc_df, df)
                write_membership_to_file('ring_cpm_res/strict_np50_leiden_cpm_r' + str(res) + '.dat', membership)

    df.to_csv('res_limit_exps_leiden_cpm_vary_res.csv')
    acc_df.to_csv('res_limit_exps_leiden_cpm_vary_res_acc.csv')


def ring_exp_cpm_k():
    df = pd.DataFrame(columns=['k', "n", "res", "method", "partition", "cluster_id", "cluster_size"])
    acc_df = pd.DataFrame(columns=['k', "n", "res", "method", "partition", "acc_measure", "acc_value"])
    res = 0.001
    for n in [1000]:
        print('n:', n)
        for k in [3, 4, 5, 6, 7, 8, 9, 10]:
            graph = nx.ring_of_cliques(num_cliques=n, clique_size=k)
            ground_truth_membership = sum([[i]*k for i in range(n)], [])

            if not os.path.exists('ring_cpm/k_' + str(k)):
                os.mkdir('ring_cpm/k_' + str(k))
            nx.write_edgelist(graph, 'ring_cpm/k_' + str(k) + '/network.dat', data=False)

            write_membership_to_file('ring_cpm/k_' + str(k) + '/community.dat', ground_truth_membership)

            #graph = gen_tree_of_cliques(k, n)
            for method in ['leiden-cpm']:
                print('original')
                partition = normal_clustering(graph, method, res_val=res)
                df, acc_df, membership = calculate_stats(graph, partition, n, k, res, method,
                                                             'Leiden-CPM', ground_truth_membership, acc_df, df)
                write_membership_to_file('ring_cpm/k_' + str(k) +'/leiden_cpm_r' + str(res) + '.dat', membership)

                print('fast ensemble')
                partition = fast_ensemble(graph, method, n_p=10, res_value=res)
                df, acc_df, membership = calculate_stats(graph, list(partition), n, k, res, method,
                                            'FastEnsemble(Leiden-CPM)', ground_truth_membership, acc_df, df)
                write_membership_to_file('ring_cpm/k_' + str(k) +'/fec_leiden_cpm_r' + str(res) + '.dat', membership)

                print('fast ensemble(tr=0.9)')
                partition = fast_ensemble(graph, method, n_p=10, res_value=res, tr=0.9)
                df, acc_df, membership = calculate_stats(graph, list(partition), n, k, res, method,
                                            'FastEnsemble(Leiden-CPM,tr=0.9)', ground_truth_membership, acc_df, df)
                write_membership_to_file('ring_cpm/k_' + str(k) +'/fec_leiden_cpm_r' + str(res) + '.dat', membership)

                print('SC n = 10')
                partition = strict_consensus(graph, method, n_p=10, res_val=res)
                df, acc_df, membership = calculate_stats(graph, list(partition), n, k, res, method,
                                             'Strict(np=10,Leiden-CPM)', ground_truth_membership, acc_df, df)
                write_membership_to_file('ring_cpm/k_' + str(k) +'strict_np10_leiden_cpm_r' + str(res) + '.dat', membership)

                print('SC n = 50')
                partition = strict_consensus(graph, method, n_p=50, res_val=res)
                df, acc_df, membership = calculate_stats(graph, list(partition), n, k, res, method,
                                             'Strict(np=50,Leiden-CPM)', ground_truth_membership, acc_df, df)
                write_membership_to_file('ring_cpm/k_' + str(k) +'/strict_np50_leiden_cpm_r' + str(res) + '.dat', membership)


    df.to_csv('res_limit_exps_leiden_cpm_k.csv')
    acc_df.to_csv('res_limit_exps_leiden_cpm_acc_k.csv')


def ring_exp_cpm_n():
    df = pd.DataFrame(columns=['k', "n", "res", "method", "partition", "cluster_id", "cluster_size"])
    acc_df = pd.DataFrame(columns=['k', "n", "res", "method", "partition", "acc_measure", "acc_value"])
    res = 0.001
    dir = 'ring_cpm_DC/'
    for n in [90, 100, 500, 1000, 5000, 10000]:
        print('n:', n)
        for k in [10]:
            graph = nx.ring_of_cliques(num_cliques=n, clique_size=k)
            ground_truth_membership = sum([[i]*k for i in range(n)], [])

            if not os.path.exists(dir+'n_' + str(n)):
                os.mkdir(dir+'n_' + str(n))
            nx.write_edgelist(graph, dir+'n_' + str(n) + '/network.dat', data=False)

            write_membership_to_file(dir+'n_' + str(n) + '/community.dat', ground_truth_membership)

            #graph = gen_tree_of_cliques(k, n)
            for method in ['leiden-cpm']:
                print('original')
                partition = normal_clustering(graph, method, res_val=res)
                df, acc_df, membership = calculate_stats(graph, partition, n, k, res, method,
                                                             'Leiden-CPM', ground_truth_membership, acc_df, df)
                write_membership_to_file(dir+'n_' + str(n) +'/leiden_cpm_r' + str(res) + '.dat', membership)

                print('fast ensemble')
                partition = fast_ensemble(graph, method, n_p=10, res_value=res)
                df, acc_df, membership = calculate_stats(graph, list(partition), n, k, res, method,
                                            'FastEnsemble(Leiden-CPM)', ground_truth_membership, acc_df, df)
                write_membership_to_file(dir+'n_' + str(n) +'/fec_leiden_cpm_r' + str(res) + '.dat', membership)

                '''
                print('fast ensemble(tr=0.9)')
                partition = fast_ensemble(graph, method, n_p=10, res_value=res, tr=0.9)
                df, acc_df, membership = calculate_stats(graph, list(partition), n, k, res, method,
                                            'FastEnsemble(Leiden-CPM,tr=0.9)', ground_truth_membership, acc_df, df)
                write_membership_to_file(dir+'n_' + str(n) +'/fec_leiden_cpm_r' + str(res) + '.dat', membership)
                '''

                print('SC n = 10')
                partition = strict_consensus(graph, method, n_p=10, res_val=res)
                df, acc_df, membership = calculate_stats(graph, list(partition), n, k, res, method,
                                             'Strict(np=10,Leiden-CPM)', ground_truth_membership, acc_df, df)
                write_membership_to_file(dir+'n_' + str(n) +'strict_np10_leiden_cpm_r' + str(res) + '.dat', membership)

                print('SC n = 50')
                partition = strict_consensus(graph, method, n_p=50, res_val=res)
                df, acc_df, membership = calculate_stats(graph, list(partition), n, k, res, method,
                                             'Strict(np=50,Leiden-CPM)', ground_truth_membership, acc_df, df)
                write_membership_to_file(dir+'n_' + str(n) +'/strict_np50_leiden_cpm_r' + str(res) + '.dat', membership)


    df.to_csv('res_limit_exps_leiden_cpm_DC.csv')
    acc_df.to_csv('res_limit_exps_leiden_cpm_acc_DC.csv')


def ring_exp():
    df = pd.DataFrame(columns=['k', "n", "res", "method", "partition", "cluster_id", "cluster_size"])
    acc_df = pd.DataFrame(columns=['k', "n", "res", "method", "partition", "acc_measure", "acc_value"])
    for n in [90, 100, 500, 1000, 5000, 10000]:
        print('n:', n)
        for k in [10]:
            #graph = nx.ring_of_cliques(num_cliques=n, clique_size=k)
            # ground_truth_membership = sum([[i]*k for i in range(n)], [])
            # #graph = gen_tree_of_cliques(k, n)
            #
            # if not os.path.exists('ring_mod/n_' + str(n)):
            #     os.mkdir('ring_mod/n_' + str(n))
            # nx.write_edgelist(graph, './network.dat', data=False)#, data=False)
            #
            # write_membership_to_file('ring_mod/n_' + str(n) + '/community.dat', ground_truth_membership)

            # reading networkX graphs don't match when is written to file
            graph, ground_truth_membership = read_network_partition('../ring_mod/n_' + str(n) + '/network.dat',
                                                                       '../ring_mod/n_' + str(n) + '/community.dat')
            #graph = nx.read_edgelist('./network.dat', nodetype=int)
            graph = nx.ring_of_cliques(num_cliques=n, clique_size=k)
            #print(nx.is_isomorphic(graph, graph2))
            #print(graph.edges)

            for method in ['leiden-mod']:
                # print('original')
                # partition = normal_clustering(graph, method)
                # df, acc_df, membership = calculate_stats(graph, partition, n, k, 'mod', method,
                #                              'Leiden-MOD', ground_truth_membership, acc_df, df)
                # write_membership_to_file('ring_mod/n_' + str(n) + '/leiden_mod.dat', membership)

                '''print('FEC n = 10')
                partition = simple_consensus(graph, method, n_p=10)
                df, acc_df = calculate_stats(graph, list(partition), n, k, 'mod', method,
                                             'FastEnsemble(np=10)', ground_truth_membership, acc_df, df)'''

                print('FEC n = 5')
                partition = fast_ensemble(graph, 'leiden-mod', n_p=5, tr=0.8)
                df, acc_df, membership = calculate_stats(graph, list(partition), n, k, 'mod', method,
                                             'FastEnsemble(np=5,Leiden-MOD)', ground_truth_membership, acc_df, df)
                write_membership_to_file('../ring_mod/n_' + str(n) + '/fec_5_leiden_mod.dat', membership)

                print('FEC n = 50')
                partition = fast_ensemble(graph, 'leiden-mod', n_p=50, tr=0.8)
                df, acc_df, membership = calculate_stats(graph, list(partition), n, k, 'mod', method,
                                                         'FastEnsemble(np=50,Leiden-MOD)', ground_truth_membership, acc_df,
                                                         df)
                write_membership_to_file('../ring_mod/n_' + str(n) + '/fec_50_leiden_mod.dat', membership)

                print('FEC n = 100')
                partition = fast_ensemble(graph, 'leiden-mod', n_p=100, tr=0.8)
                df, acc_df, membership = calculate_stats(graph, list(partition), n, k, 'mod', method,
                                                         'FastEnsemble(np=100,Leiden-MOD)', ground_truth_membership, acc_df,
                                                         df)
                write_membership_to_file('../ring_mod/n_' + str(n) + '/fec_100_leiden_mod.dat', membership)

                '''print('FEC(louvain) n = 10')
                partition = fast_ensemble(graph, 'louvain', n_p=10)
                df, acc_df = calculate_stats(graph, list(partition), n, k, 'mod', method,
                                             'FastEnsemble(Louvain)', ground_truth_membership, acc_df, df)'''

                # print('ECG n = 10')
                # g = ig.Graph.from_networkx(graph)
                # ec = g.community_ecg(ens_size=10)
                # df, acc_df, membership = calculate_stats(graph, list(ec), n, k, 'mod', method,
                #                              'ECG', ground_truth_membership, acc_df, df)
                # write_membership_to_file('ring_mod/n_' + str(n) + '/ecg.dat', membership)

                print('SC n = 5')
                partition = strict_consensus(graph, method, n_p=5)
                df, acc_df, membership = calculate_stats(graph, list(partition), n, k, 'mod', method,
                                                         'Strict(np=5,Leiden-MOD)', ground_truth_membership, acc_df,
                                                         df)
                write_membership_to_file('../ring_mod/n_' + str(n) + '/strict_np5_leiden_mod.dat', membership)

                print('SC n = 100')
                partition = strict_consensus(graph, method, n_p=100)
                df, acc_df, membership = calculate_stats(graph, list(partition), n, k, 'mod', method,
                                                         'Strict(np=100,Leiden-MOD)', ground_truth_membership, acc_df,
                                                         df)
                write_membership_to_file('../ring_mod/n_' + str(n) + '/strict_np100_leiden_mod.dat', membership)

                '''print('SC n = 10')
                partition = strict_consensus(graph, method, n_p=10)
                df, acc_df, membership = calculate_stats(graph, list(partition), n, k, 'mod', method,
                                             'Strict(np=10,Leiden-MOD)', ground_truth_membership, acc_df, df)
                write_membership_to_file('ring_mod/n_' + str(n) + '/strict_np10_leiden_mod.dat', membership)

                print('SC n = 50')
                partition = strict_consensus(graph, method, n_p=50)
                df, acc_df, membership = calculate_stats(graph, list(partition), n, k, 'mod', method,
                                             'Strict(np=50,Leiden-MOD)', ground_truth_membership, acc_df, df)
                write_membership_to_file('ring_mod/n_' + str(n) + '/strict_np50_leiden_mod.dat', membership)'''

    df.to_csv('res_limit_exps_leiden_mod_ring_vary_np.csv')
    acc_df.to_csv('res_limit_exps_leiden_mod_ring_acc_vary_np.csv')


def tree_mod_exp():
    df = pd.DataFrame(columns=['k', "n", "res", "method", "partition", "cluster_id", "cluster_size"])
    acc_df = pd.DataFrame(columns=['k', "n", "res", "method", "partition", "acc_measure", "acc_value"])
    for n in [90, 100, 500, 1000, 5000]:
        print('n:', n)
        for k in [10]:
            graph = gen_tree_of_cliques(k, n)
            ground_truth_membership = sum([[i]*k for i in range(n)], [])

            if not os.path.exists('tree_mod/n_' + str(n)):
                os.mkdir('tree_mod/n_' + str(n))
            nx.write_edgelist(graph, 'tree_mod/n_' + str(n) + '/network.dat', data=False)

            write_membership_to_file('tree_mod/n_' + str(n) + '/community.dat', ground_truth_membership)

            for method in ['leiden-mod']:
                print('leiden-mod')
                partition = normal_clustering(graph, method)
                df, acc_df, membership = calculate_stats(graph, partition, n, k, 'mod', method,
                                             'Leiden-mod', ground_truth_membership, acc_df, df)
                write_membership_to_file('tree_mod/n_' + str(n) + '/leiden_mod.dat', membership)

                print('louvain')
                partition = normal_clustering(graph, 'louvain')
                df, acc_df, membership = calculate_stats(graph, partition, n, k, 'mod', 'louvain',
                                             'Louvain', ground_truth_membership, acc_df, df)
                write_membership_to_file('tree_mod/n_' + str(n) + '/louvain.dat', membership)

                '''print('FEC n = 10')
                partition = simple_consensus(graph, method, n_p=10)
                df, acc_df = calculate_stats(graph, list(partition), n, k, 'mod', method,
                                             'FastEnsemble(np=10)', ground_truth_membership, acc_df, df)'''

                print('TC n = 10')
                partition = fast_ensemble(graph, 'leiden-mod', n_p=10, tr=0.8)
                df, acc_df, membership = calculate_stats(graph, list(partition), n, k, 'mod', method,
                                             'FastEnsemble(Leiden-mod)', ground_truth_membership, acc_df, df)
                write_membership_to_file('tree_mod/n_' + str(n) + '/fec_leiden_mod.dat', membership)

                print('FEC(louvain) n = 10')
                partition = fast_ensemble(graph, 'louvain', n_p=10)
                df, acc_df, membership = calculate_stats(graph, list(partition), n, k, 'mod', 'louvain',
                                             'FastEnsemble(Louvain)', ground_truth_membership, acc_df, df)
                write_membership_to_file('tree_mod/n_' + str(n) + '/fec_louvain.dat', membership)


                print('ECG n = 10')
                g = ig.Graph.from_networkx(graph)
                ec = g.community_ecg(ens_size=10)
                df, acc_df, membership = calculate_stats(graph, list(ec), n, k, 'mod', method,
                                             'ECG', ground_truth_membership, acc_df, df)
                write_membership_to_file('tree_mod/n_' + str(n) + '/ecg.dat', membership)

                print('SC n = 10')
                partition = strict_consensus(graph, method, n_p=10)
                df, acc_df, membership = calculate_stats(graph, list(partition), n, k, 'mod', method,
                                             'Strict(np=10,Leiden-mod)', ground_truth_membership, acc_df, df)
                write_membership_to_file('tree_mod/n_' + str(n) + '/strict_np10_leiden_mod.dat', membership)

                print('SC n = 50')
                partition = strict_consensus(graph, method, n_p=50)
                df, acc_df, membership = calculate_stats(graph, list(partition), n, k, 'mod', method,
                                             'Strict(np=50,Leiden-mod)', ground_truth_membership, acc_df, df)
                write_membership_to_file('tree_mod/n_' + str(n) + '/strict_np50_leiden_mod.dat', membership)

    df.to_csv('res_limit_exps_leiden_mod_tree_final.csv')
    acc_df.to_csv('res_limit_exps_leiden_mod_tree_acc_final.csv')


'''if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Threshold Consensus")
    parser.add_argument("-n", "--edgelist", type=str,  required=True,
                        help="Network edge-list file")
    parser.add_argument("-alg", "--algorithm", type=str, required=False,
                        help="Clustering algorithm (leiden-mod or leiden-cpm)", default='leiden-cpm')
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Output community file")
    parser.add_argument("-gt", "--groundtruth", type=str, required=False,
                        help="Ground-truth")
    parser.add_argument("-t", "--threshold", type=float, required=False,
                        help="Threshold value", default=0.8)
    parser.add_argument("-r", "--resolution", type=float, required=False,
                        help="Resolution value for leiden-cpm", default=0.01)
    parser.add_argument("-p", "--partitions", type=int, required=False,
                        help="Number of partitions in consensus clustering", default=10)
    parser.add_argument("-d", "--delta", type=float, required=False,
                        help="Convergence parameter", default=0.02)
    parser.add_argument("-m", "--maxiter", type=int, required=False,
                        help="Maximum number of iterations in simple consensus", default=3)
    parser.add_argument("-rl", "--relabel", required=False, action='store_true',
                        help="Relabel network nodes from 0 to #nodes-1.", default=False)

    args = parser.parse_args()
    net = nx.read_edgelist(args.edgelist, nodetype=int)

    # relabeling nodes from 0 to n-1
    if args.relabel:
        mapping = dict(zip(sorted(net), range(0, net.number_of_nodes())))
        net = nx.relabel_nodes(net, mapping)
        reverse_mapping = {y: x for x, y in mapping.items()}


    #ecg_out = get_communities(net, 'ecg', 1234)
    #print(ecg_out)

    n_p = args.partitions
    #leiden_alg_list = ['leiden-mod'] * n_p
    #lou_lei_alg_list = ['leiden-mod'] * int(n_p / 2) + ['leiden-cpm'] * (n_p - int(n_p / 2))
    sc = simple_consensus(net, algorithm='leiden-mod', n_p=args.partitions, thresh=args.threshold, max_iter=3)
    if args.algorithm == 'leiden-cpm':
        sc = simple_consensus(net, algorithm='leiden-cpm', n_p=args.partitions, thresh=args.threshold, max_iter=3, res_value=args.resolution)
    elif args.algorithm == 'leiden-mod':
        sc = simple_consensus(net, algorithm='leiden-mod', n_p=args.partitions, thresh=args.threshold, max_iter=3)

    #sc = get_communities(initialize(net, 1.0), 'leiden-cpm', 1234, res_val=0.0001)
    #sc = threshold_consensus(net, algorithm='leiden-cpm', n_p=args.partitions, tr=args.threshold, res_value=args.resolution)
    #partition = communities_to_dict(sc)
    keys = list(sc.keys())
    keys.sort()
    #print(keys)
    membership_dict = {i: sc[i] for i in keys}
    #print(membership)
    membership = list(membership_dict.values())
    #print(membership)

    with open(args.output, 'w') as out_file:
        writer = csv.writer(out_file, delimiter=' ')
        for i in range(len(membership)):
            if args.relabel:
                writer.writerow([reverse_mapping[i]] + [membership[i]])
            else:
                writer.writerow([i] + [membership[i]])

    gt_membership = [''] * net.number_of_nodes()
    with open(args.groundtruth) as fgt:
        for line in fgt:
            i, m = line.strip().split()
            gt_membership[int(i)-1] = int(m)-1

    #print(gt_membership)

    #print(len(membership), len(gt_membership))

    nmi, ami, ari, precision, recall, f1_score, fnr, fpr = measure_accuracy(gt_membership, membership)
    print("Normalized mutual information (NMI): ", nmi)
    print("Adjusted rand index (ARI): ", ari)
    print("Adjusted mutual information (AMI): ", ami)
    print("False positive rate (FPR), False negative rate (FNR):", fpr, fnr)
    print("Precision, Recall, F1-score:", precision, recall, f1_score)'''


def read_fast_consensus_results():
    # this doesn't check for node labeling and might be incorrect for calculating accuracy
    #dir_path = '/projects/tallis/shared/min_consensus_clustering/'
    #dir_path = '/projects/illinois/eng/cs/warnow/shared/min_consensus_clustering/'
    dir_path = '../'

    ## ring-mod exp
    df = pd.DataFrame(columns=['k', "n", "res", "method", "partition", "cluster_id", "cluster_size"])
    for n in [90, 100, 500, 1000, 5000]:
        graph, partition = read_network_partition(dir_path+'tree_mod/n_' + str(n) + '/network.dat', dir_path+'tree_mod/n_' + str(n) + '/fastconsensus_t0.2_d0.02_np10_louvain.tsv')
        partition = group_to_partition(membership_list_to_dict(partition))
        cluster_sizes, node_coverage = partition_statistics(graph, partition)
        for i in range(len(cluster_sizes)):
            df.loc[len(df.index)] = [n, 10, 'mod', 'louvain', 'FastConsensus(Louvain)', i,
                                     cluster_sizes[i]]

    df.to_csv('fast_consensus_mod_tree.csv')

    ## erdos-renyi exp
    '''df = pd.DataFrame(columns=["n", "m", "p", "method", "partition", "cluster_id", "cluster_size"])
    for p in [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1]:
        graph, partition = read_network_partition(dir_path+'erdos_renyi/p_' + str(p) + '/network.dat',
                                                  dir_path+'erdos_renyi/p_' + str(p) + '/fec_0.9_leiden_mod.dat') # 'original_fastconsensus_louvain_np10_t0.2_d0.02.dat'
        partition = group_to_partition(membership_list_to_dict(partition))
        cluster_sizes, node_coverage = partition_statistics(graph, partition)
        for i in range(len(cluster_sizes)):
             df.loc[len(df.index)] = [1000, graph.number_of_edges(), p, 'leiden-mod', 'FastEnsemble(Leiden-MOD,tr=0.9)', i,
                                         cluster_sizes[i]]
            #df.loc[len(df.index)] = [1000, graph.number_of_edges(), p, 'louvain', 'FastConsensus(Louvain)', i,
            #                             cluster_sizes[i]]

    df.to_csv('fast_ensemble_0.9_erdos_renyi.csv')'''

    '''
    ## erdos-renyi+lfr exp
    df = pd.DataFrame(columns=["n", "m", "p", "method", "partition", "cluster_id", "cluster_size"])
    for p in [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1]:
        graph, partition = read_network_partition(dir_path + 'erdos_renyi_lfr/p_' + str(p) + '/network.dat',
                                                  dir_path + 'erdos_renyi_lfr/p_' + str(
                                                      p) + '/fec_0.9_leiden_mod.dat')
        partition = group_to_partition(membership_list_to_dict(partition))
        cluster_sizes, node_coverage = partition_statistics(graph, partition)
        for i in range(len(cluster_sizes)):
            df.loc[len(df.index)] = [1000, graph.number_of_edges(), p, 'leiden-mod', 'FastEnsemble(Leiden-mod,tr=0.9)', i,
                                     cluster_sizes[i]]

    df.to_csv('fast_ensemble_0.9_erdos_renyi_lfr.csv')
    '''
    '''
    ## erdos-renyi+lfr exp
    df = pd.DataFrame(columns=["n", "m", "p", "method", "partition", "cluster_id", "cluster_size"])
    for p in [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1]:
        graph, partition = read_network_partition(dir_path + 'erdos_renyi_ring/p_' + str(p) + '/network.dat',
                                                  dir_path + 'erdos_renyi_ring/p_' + str(
                                                      p) + '/original_fastconsensus_louvain_np10_t0.2_d0.02.dat')
        partition = group_to_partition(membership_list_to_dict(partition))
        cluster_sizes, node_coverage = partition_statistics(graph, partition)
        for i in range(len(cluster_sizes)):
            df.loc[len(df.index)] = [2000, graph.number_of_edges(), p, 'louvain', 'FastConsensus(Louvain)', i,
                                     cluster_sizes[i]]

    df.to_csv('fast_ensemble_erdos_renyi_ring.csv')
    '''
    '''
    ## lfr training experiment
    n=10000
    acc_df = pd.DataFrame(columns=["n", "mu", "partition", "acc_measure", "acc_value"])
    for mu in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
        graph, partition = read_network_partition(dir_path + 'lfr_training/mu_' + str(mu) + '/network.dat',
                                                  dir_path + 'lfr_training/mu_' + str(
                                                      mu) + '/louvain.dat')
        graph, groundtruth = read_network_partition(dir_path + 'lfr_training/mu_' + str(mu) + '/network.dat',
                                                  dir_path + 'lfr_training/mu_' + str(
                                                      mu) + '/community.dat')
        #partition = group_to_partition(membership_list_to_dict(partition))
        #groundtruth = group_to_partition(membership_list_to_dict(groundtruth))
        print(list(partition))
        acc_df, membership = calculate_stats_mu(list(partition), n, mu, 'Louvain', list(groundtruth), acc_df)
        acc_df.to_csv('lfr_training_louvain.csv')
    '''


if __name__ == "__main__":
    #read_fast_consensus_results()
    #tree_mod_exp()
    #ring_exp_cpm()
    ring_exp()
    #ring_exp_cpm_n()
    #ring_exp_cpm_k()
    #ring_exp_cpm()
    #lfr_exp_mu()
    #lfr_exp_deg()
    #lfr_exp_size()
