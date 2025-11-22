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


def fast_ensemble(G, algorithm='leiden-mod', n_p=10, tr=0.8, res_value=0.01, random_seed=1):
    graph = G.copy()
    graph = initialize(graph, 1)
    partitions = [get_communities(graph, algorithm, i*random_seed, res_val=res_value) for i in range(n_p)]

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



def normal_clustering(graph, algorithm, res_val=0.01, seed=1234):
    if algorithm == 'leiden-cpm':
        return leidenalg.find_partition(ig.Graph.from_networkx(graph),
                                              leidenalg.CPMVertexPartition,
                                              resolution_parameter=res_val,
                                              n_iterations=2,
                                              seed=seed).as_cover()
    elif algorithm == 'leiden-mod':
        return leidenalg.find_partition(ig.Graph.from_networkx(graph),
                                         leidenalg.ModularityVertexPartition,
                                         seed=seed).as_cover()
    elif algorithm == 'louvain':
        #return cl.best_partition(graph)
        return list(group_to_partition(cl.best_partition(graph)))



def measure_accuracy(mem_true, mem_est):
    '''n = len(mem_true)
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
    fpr = fp / (fp + tn)'''

    nmi = normalized_mutual_info_score(mem_true, mem_est)
    ari = adjusted_rand_score(mem_true, mem_est)
    ami = adjusted_mutual_info_score(mem_true, mem_est)

    return nmi, ami, ari, 0, 0, 0, 0, 0
    #return nmi, ami, ari, precision, recall, f1_score, fnr, fpr


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


def calculate_stats_mu(partition, n, mu, rep, partition_name, ground_truth_membership, acc_df):
    membership = partition# get_membership_list_from_dict(partition, n)
    #print(membership)
    nmi, ami, ari, precision, recall, f1_score, fnr, fpr = measure_accuracy(ground_truth_membership, membership)
    print("Normalized mutual information (NMI): ", nmi)
    print("Adjusted rand index (ARI): ", ari)
    print("Adjusted mutual information (AMI): ", ami)
    acc_df.loc[len(acc_df.index)] = [n, mu, rep, partition_name, 'NMI', nmi]
    acc_df.loc[len(acc_df.index)] = [n, mu, rep, partition_name, 'ARI', ari]
    acc_df.loc[len(acc_df.index)] = [n, mu, rep, partition_name, 'AMI', ami]
    return acc_df, membership


def calculate_stats_mu_pairwise(membership1, n, mu, rep1, rep2, partition_name, membership2, acc_df):
    nmi, ami, ari, precision, recall, f1_score, fnr, fpr = measure_accuracy(membership1, membership2)
    print("Normalized mutual information (NMI): ", nmi)
    print("Adjusted rand index (ARI): ", ari)
    print("Adjusted mutual information (AMI): ", ami)
    acc_df.loc[len(acc_df.index)] = [n, mu, rep1, rep2, partition_name, 'NMI', nmi]
    acc_df.loc[len(acc_df.index)] = [n, mu, rep1, rep2, partition_name, 'ARI', ari]
    acc_df.loc[len(acc_df.index)] = [n, mu, rep1, rep2, partition_name, 'AMI', ami]
    return acc_df, membership1


def calculate_stats_deg(partition, n, d, mu, partition_name, ground_truth_membership, acc_df):
    membership = get_membership_list_from_dict(partition, n)
    nmi, ami, ari, precision, recall, f1_score, fnr, fpr = measure_accuracy(ground_truth_membership, membership)
    print("Normalized mutual information (NMI): ", nmi)
    print("Adjusted rand index (ARI): ", ari)
    acc_df.loc[len(acc_df.index)] = [d, mu, partition_name, 'NMI', nmi]
    acc_df.loc[len(acc_df.index)] = [d, mu, partition_name, 'ARI', ari]
    return acc_df, membership


def stability_exp():
    acc_df = pd.DataFrame(columns=['n', "mu", "rep", "method", "acc_measure", "acc_value"])
    acc_df_pairwise = pd.DataFrame(columns=['n', "mu", "rep1", "rep2", "method", "acc_measure", "acc_value"])
    dir_path = '../lfr_training/default/'
    n = 10000
    for mu in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
        print(mu)
        _, groundtruth = read_network_partition(dir_path + 'mu_' + str(mu) + '/network.dat',
                                                                   dir_path + 'mu_' + str(mu) + '/community.dat')
        graph = nx.generators.community.LFR_benchmark_graph(n=10000, tau1=3, tau2=1.5, mu=mu, average_degree=10,
                                                            min_community=10, seed=1932)
        #groundtruth = list(group_to_partition(membership_list_to_dict(groundtruth)))
        #print(groundtruth)

        # for rep in range(1, 3):
        #     print('ECG')
        #     g = ig.Graph.from_networkx(graph)
        #     ec = g.community_ecg(ens_size=10)
        #     acc_df, membership = calculate_stats_mu(list(ec), n, mu, rep, 'ECG', groundtruth, acc_df)
        #     write_membership_to_file(dir_path + 'mu_' + str(mu) + '/ecg_r'+str(rep).zfill(2)+'.dat', membership)

            # print('Leiden-MOD')
            # partition = normal_clustering(graph, 'leiden-mod', seed=rep)
            # acc_df, membership = calculate_stats_mu(list(partition), n, mu, rep, 'Leiden-MOD', groundtruth, acc_df)
            # write_membership_to_file(dir_path + 'mu_' + str(mu) + '/leiden_mod_r'+str(rep).zfill(2)+'.dat', membership)
            #
            # print('FastEnsemble')
            # partition = fast_ensemble(graph, 'leiden-mod', n_p=10, tr=0.8)
            # acc_df, membership = calculate_stats_mu(list(partition), n, mu, rep, 'FastEnsemble', groundtruth, acc_df)
            # write_membership_to_file(dir_path + 'mu_' + str(mu) + '/fec_rs_leiden_mod_r'+str(rep).zfill(2)+'.dat', membership)

        leiden_partitions = []
        fast_ensemble_partitions = []
        ecg_partitions=[]
        for rep in range(1, 3):
            _, partition = read_network_partition(dir_path + 'mu_' + str(mu) + '/network.dat',
                                                        dir_path + 'mu_' + str(mu) + '/leiden_mod_r'+str(rep).zfill(2)+'.dat')
            leiden_partitions.append(list(partition))
            acc_df, _ = calculate_stats_mu(list(partition), n, mu, rep, 'Leiden-mod', groundtruth, acc_df)

            _, partition = read_network_partition(dir_path + 'mu_' + str(mu) + '/network.dat',
                                                  dir_path + 'mu_' + str(mu) + '/fec_leiden_mod_r' + str(rep).zfill(2) + '.dat')
            fast_ensemble_partitions.append(list(partition))
            acc_df, _ = calculate_stats_mu(list(partition), n, mu, rep, 'FastEnsemble', groundtruth, acc_df)

            _, partition = read_network_partition(dir_path + 'mu_' + str(mu) + '/network.dat',
                                                  dir_path + 'mu_' + str(mu) + '/ecg_r'+str(rep).zfill(2)+'.dat')
            ecg_partitions.append(list(partition))
            acc_df, _ = calculate_stats_mu(list(partition), n, mu, rep, 'ECG', groundtruth, acc_df)

        for rep1 in range(0, 2):
            for rep2 in range(rep1+1, 2):
                print(rep1, rep2)
                print('ECG')
                acc_df_pairwise, _ = calculate_stats_mu_pairwise(ecg_partitions[rep1], n, mu, rep1, rep2,
                                                        'ECG', ecg_partitions[rep2],
                                                        acc_df_pairwise)
                print('Leiden-mod')
                acc_df_pairwise, _ = calculate_stats_mu_pairwise(leiden_partitions[rep1], n, mu, rep1, rep2,
                                                                 'Leiden-mod', leiden_partitions[rep2],
                                                                acc_df_pairwise)
                print('FastEnsemble')
                acc_df_pairwise, _ = calculate_stats_mu_pairwise(fast_ensemble_partitions[rep1], n, mu, rep1, rep2,
                                                                 'FastEnsemble', fast_ensemble_partitions[rep2],
                                                                 acc_df_pairwise)

    acc_df.to_csv('stability_exp.csv')
    acc_df_pairwise.to_csv('stability_exp_pairwise.csv')




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
    stability_exp()
