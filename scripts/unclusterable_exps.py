import leidenalg
import networkx as nx
import igraph as ig
import community.community_louvain as cl
import numpy as np
import pandas as pd
import os
import csv
import seaborn as sns
from networkx.algorithms.community import modularity
from sklearn.metrics.cluster import adjusted_mutual_info_score
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.cluster import normalized_mutual_info_score
import matplotlib.pyplot as plt


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

    #return cluster_num, min_size, max_size, mean_size, median_size, singletons_num, non_singleton_num, modularity_score, coverage


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

def get_membership_list_from_dict(membership_dict):
    memberships = []
    for i in range(len(membership_dict)):
        memberships.append(membership_dict[i])
    return memberships


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

def get_communities(graph, algorithm, seed, res_val=0.01):
    if algorithm == 'louvain':
        return cl.best_partition(graph, random_state=seed, weight='weight')
    elif algorithm == 'leiden-cpm':
        return communities_to_dict(leidenalg.find_partition(ig.Graph.from_networkx(graph),
                                                  leidenalg.CPMVertexPartition,
                                                  resolution_parameter=res_val,
                                                  n_iterations=2).as_cover())
    elif algorithm == 'leiden-mod':
        return communities_to_dict(leidenalg.find_partition(ig.Graph.from_networkx(graph),
                                        leidenalg.ModularityVertexPartition,
                                        seed=seed).as_cover())
    elif algorithm == 'leiden-mod-weights':
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


def simple_consensus(G, algorithm='leiden-mod-weights', n_p=10, thresh=0.5, delta=0.02, max_iter=10):
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
        partitions = [get_communities(graph, algorithm, i) for i in range(n_p)]

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
    return group_to_partition(get_communities(graph, algorithm, 0))

def false_negative():
    return

def false_positive():
    return

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

def fast_ensemble(G, algorithm='leiden-cpm', n_p=10, tr=0.8, res_value=0.01):
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
        return group_to_partition(cl.best_partition(graph))

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
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1_score = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
    fnr = fn / (fn + tp) if (fn + tp) > 0 else 0
    fpr = fp / (fp + tn) if (fp + tn) > 0 else 0'''

    nmi = normalized_mutual_info_score(mem_true, mem_est)
    ari = adjusted_rand_score(mem_true, mem_est)
    #ami = adjusted_mutual_info_score(mem_true, mem_est)

    #print(tp, tn, fp, fn)

    #return nmi, ami, ari, precision, recall, f1_score, fnr, fpr
    return nmi, 0, ari, 0, 0, 0, 0, 0


def er_ring():
    n = 2000
    k = 10
    df = pd.DataFrame(columns=["n", "m", "p", "method", "partition", "cluster_id", "cluster_size"])
    df_coverage = pd.DataFrame(columns=["n", "m", "p", "method", "partition", "coverage"])
    acc_df = pd.DataFrame(columns=["n", "m", "p", "method", "partition", "acc_measure", "acc_value"])
    ring = nx.ring_of_cliques(num_cliques=100, clique_size=10)

    ground_truth_membership = sum([[i] * k for i in range(100)], [])
    ground_truth_membership = [(1000 + i) for i in range(1, 1001)] + ground_truth_membership
    #print(ground_truth_membership)

    for p in [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1]:
        graph = nx.erdos_renyi_graph(n=1000, p=p)
        graph = nx.disjoint_union_all([graph, ring])
        graph.add_edge(0, graph.number_of_nodes() - 1)
        m = graph.number_of_edges()
        if not os.path.exists('erdos_renyi_ring/p_'+str(p)):
            os.mkdir('erdos_renyi_ring/p_'+str(p))
        nx.write_edgelist(graph, 'erdos_renyi_ring/p_'+str(p) + '/network.dat', data=False)
        write_membership_to_file('erdos_renyi_ring/p_' + str(p) + '/community.dat', ground_truth_membership)
        print(p, m)

        for method in ['leiden-mod']:
            partition = normal_clustering(graph, method)
            df, df_coverage, acc_df, membership = calculate_stats_er(graph, partition, df_coverage, n, m, p, method,
                                                         'Leiden-MOD', ground_truth_membership, acc_df, df)
            write_membership_to_file('erdos_renyi_ring/p_' + str(p) + '/leiden_mod.dat', membership)

            partition = fast_ensemble(graph, method, n_p=10)
            df, df_coverage, acc_df, membership = calculate_stats_er(graph, list(partition), df_coverage, n, m, p, method,
                                                         'FastEnsemble(Leiden-MOD)', ground_truth_membership, acc_df,
                                                         df)
            write_membership_to_file('erdos_renyi_ring/p_' + str(p) + '/fec_leiden_mod.dat', membership)

            g = ig.Graph.from_networkx(graph)
            ec = g.community_ecg(ens_size=10)
            df, df_coverage, acc_df, membership = calculate_stats_er(graph, list(ec), df_coverage, n, m, p, method,
                                                         'ECG', ground_truth_membership, acc_df, df)
            write_membership_to_file('erdos_renyi_ring/p_' + str(p) + '/ecg.dat', membership)

            partition = strict_consensus(graph, method, n_p=10)
            df, df_coverage, acc_df, membership = calculate_stats_er(graph, list(partition), df_coverage, n, m, p, method,
                                                         'Strict(np=10)+Leiden-MOD', ground_truth_membership, acc_df,
                                                         df)
            write_membership_to_file('erdos_renyi_ring/p_' + str(p) + '/strict_np10_leiden_mod.dat', membership)


            partition = strict_consensus(graph, method, n_p=50)
            df, df_coverage, acc_df, membership = calculate_stats_er(graph, list(partition), df_coverage, n, m, p, method,
                                                         'Strict(np=50)+Leiden-MOD', ground_truth_membership, acc_df,
                                                         df)
            write_membership_to_file('erdos_renyi_ring/p_' + str(p) + '/strict_np50_leiden_mod.dat', membership)

    df.to_csv('erdos_renyi_ring_exps_leiden_mod.csv')
    df_coverage.to_csv('erdos_renyi_ring_exps_leiden_mod_coverage.csv')
    acc_df.to_csv('erdos_renyi_ring_exps_leiden_acc.csv')


def calculate_stats_er(graph, partition, df_coverage, n, m, p, method, partition_name, ground_truth_membership, acc_df, df):
    #print(partition)
    cluster_sizes, node_coverage = partition_statistics(graph, partition)
    df_coverage.loc[len(df_coverage.index)] = [n, m, p, method, partition_name, node_coverage]
    for i in range(len(cluster_sizes)):
        df.loc[len(df.index)] = [n, m, p, method, partition_name, i, cluster_sizes[i]]
    membership = get_membership_list_from_dict(partition, n)
    nmi, ami, ari, precision, recall, f1_score, fnr, fpr = measure_accuracy(ground_truth_membership, membership)
    print("Normalized mutual information (NMI): ", nmi)
    print("Adjusted rand index (ARI): ", ari)
    #print("Adjusted mutual information (AMI): ", ami)
    #print("False positive rate (FPR), False negative rate (FNR):", fpr, fnr)
    #print("Precision, Recall, F1-score:", precision, recall, f1_score)
    acc_df.loc[len(acc_df.index)] = [n, m, p, method, partition_name, 'NMI', nmi]
    #acc_df.loc[len(acc_df.index)] = [n, m, p, method, partition_name, 'AMI', ami]
    acc_df.loc[len(acc_df.index)] = [n, m, p, method, partition_name, 'ARI', ari]
    #acc_df.loc[len(acc_df.index)] = [n, m, p, method, partition_name, 'FPR', fpr]
    #acc_df.loc[len(acc_df.index)] = [n, m, p, method, partition_name, 'FNR', fnr]
    #acc_df.loc[len(acc_df.index)] = [n, m, p, method, partition_name, 'Precision', precision]
    #acc_df.loc[len(acc_df.index)] = [n, m, p, method, partition_name, 'Recall', recall]
    #acc_df.loc[len(acc_df.index)] = [n, m, p, method, partition_name, 'F1-Score', f1_score]
    return df, df_coverage, acc_df, membership


def write_membership_to_file(file_path, memberhsip):
    with open(file_path, 'w') as out_file:
        writer = csv.writer(out_file, delimiter=' ')
        for i in range(len(memberhsip)):
            writer.writerow([i] + [memberhsip[i]])


def erdos_renyi():
    df = pd.DataFrame(columns=["n", "m", "p", "method", "partition", "cluster_id", "cluster_size"])
    df_coverage = pd.DataFrame(columns=["n", "m", "p", "method", "partition", "coverage"])
    acc_df = pd.DataFrame(columns=["n", "m", "p", "method", "partition", "acc_measure", "acc_value"])
    n = 1000
    ground_truth_membership = [i for i in range(n)]
    for p in [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1]:
        graph = nx.erdos_renyi_graph(n=1000, p=p)

        if not os.path.exists('erdos_renyi/p_'+str(p)):
            os.mkdir('erdos_renyi/p_'+str(p))
        nx.write_edgelist(graph, 'erdos_renyi/p_'+str(p) + '/network.dat', data=False)

        write_membership_to_file('erdos_renyi/p_' + str(p) + '/community.dat', ground_truth_membership)

        m = graph.number_of_edges()
        print(p, m)
        for method in ['leiden-mod']:
            partition = normal_clustering(graph, method)
            df, df_coverage, acc_df, membership = calculate_stats_er(graph, partition, df_coverage, n, m, p, method,
                                                         'Leiden-MOD', ground_truth_membership, acc_df, df)
            write_membership_to_file('erdos_renyi/p_' + str(p) + '/leiden_mod.dat', membership)

            '''print('TC n = 10')
            partition = fast_ensemble(graph, method, n_p=10)
            df, acc_df = calculate_stats(graph, list(partition), n, k, 'mod', method,
                                         'FastEnsemble(Leiden-MOD)', ground_truth_membership, acc_df, df)'''

            partition = fast_ensemble(graph, method, n_p=10)
            df, df_coverage, acc_df, membership = calculate_stats_er(graph, list(partition), df_coverage, n, m, p, method,
                                                         'FastEnsemble(Leiden-MOD)', ground_truth_membership, acc_df,df)
            write_membership_to_file('erdos_renyi/p_' + str(p) + '/fec_leiden_mod.dat', membership)

            g = ig.Graph.from_networkx(graph)
            ec = g.community_ecg(ens_size=10)
            df, df_coverage, acc_df, membership = calculate_stats_er(graph, list(ec), df_coverage, n, m, p, method,
                                                         'ECG', ground_truth_membership, acc_df,df)
            write_membership_to_file('erdos_renyi/p_' + str(p) + '/ecg.dat', membership)

            partition = strict_consensus(graph, method, n_p=10)
            df, df_coverage, acc_df, membership = calculate_stats_er(graph, list(partition), df_coverage, n, m, p, method,
                                                         'Strict(np=10,Leiden-MOD)', ground_truth_membership, acc_df, df)
            write_membership_to_file('erdos_renyi/p_' + str(p) + '/strict_np10_leiden_mod.dat', membership)

            partition = strict_consensus(graph, method, n_p=50)
            df, df_coverage, acc_df, membership = calculate_stats_er(graph, list(partition), df_coverage, n, m, p, method,
                                                         'Strict(np=50,Leiden-MOD)', ground_truth_membership, acc_df, df)
            write_membership_to_file('erdos_renyi/p_' + str(p) + '/strict_np50_leiden_mod.dat', membership)

    df.to_csv('erdos_renyi_exps_leiden_mod.csv')
    df_coverage.to_csv('erdos_renyi_exps_leiden_mod_coverage.csv')
    acc_df.to_csv('erdos_renyi_exps_leiden_acc.csv')


def er_lfr():
    n = 2000
    df = pd.DataFrame(columns=["n", "m", "p", "method", "partition", "cluster_id", "cluster_size"])
    df_coverage = pd.DataFrame(columns=["n", "m", "p", "method", "partition", "coverage"])
    acc_df = pd.DataFrame(columns=["n", "m", "p", "method", "partition", "acc_measure", "acc_value"])
    lfr = nx.generators.community.LFR_benchmark_graph(n=1000, tau1=3, tau2=1.5, mu=0.1, average_degree=9.198,
                                                      min_community=45, seed=10)
    communities = {frozenset(lfr.nodes[v]["community"]) for v in lfr}
    communities = communities_to_dict(communities)
    keys = list(communities.keys())
    keys.sort()
    membership = {i: communities[i] for i in keys}
    ground_truth_membership = list(membership.values())
    ground_truth_membership = [(len(keys)+i) for i in range(1, 1001)] + ground_truth_membership
    for p in [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1]:
        graph = nx.erdos_renyi_graph(n=1000, p=p)
        graph = nx.disjoint_union_all([graph, lfr])
        graph.add_edge(0, graph.number_of_nodes() - 1)
        if not os.path.exists('erdos_renyi_lfr/p_'+str(p)):
            os.mkdir('erdos_renyi_lfr/p_'+str(p))
        nx.write_edgelist(graph, 'erdos_renyi_lfr/p_'+str(p) + '/network.dat', data=False)
        write_membership_to_file('erdos_renyi_lfr/p_' + str(p) + '/community.dat', ground_truth_membership)
        m = graph.number_of_edges()
        print(p, m)
        for method in ['leiden-mod']:
            partition = normal_clustering(graph, method)

            df, df_coverage, acc_df, membership = calculate_stats_er(graph, partition, df_coverage, n, m, p, method,
                                                         'Leiden-MOD', ground_truth_membership, acc_df, df)
            write_membership_to_file('erdos_renyi_lfr/p_' + str(p) + '/leiden_mod.dat', membership)

            partition = fast_ensemble(graph, method, n_p=10)
            df, df_coverage, acc_df, membership = calculate_stats_er(graph, list(partition), df_coverage, n, m, p, method,
                                                         'FastEnsemble(Leiden-MOD)', ground_truth_membership, acc_df, df)
            write_membership_to_file('erdos_renyi_lfr/p_' + str(p) + '/fec_leiden_mod.dat', membership)

            g = ig.Graph.from_networkx(graph)
            ec = g.community_ecg(ens_size=10)
            df, df_coverage, acc_df, membership = calculate_stats_er(graph, list(ec), df_coverage, n, m, p, method,
                                                         'ECG', ground_truth_membership, acc_df, df)
            write_membership_to_file('erdos_renyi_lfr/p_' + str(p) + '/ecg.dat', membership)

            partition = strict_consensus(graph, method, n_p=10)
            df, df_coverage, acc_df, membership = calculate_stats_er(graph, list(partition), df_coverage, n, m, p, method,
                                                         'Strict(np=10,Leiden-MOD)', ground_truth_membership, acc_df, df)
            write_membership_to_file('erdos_renyi_lfr/p_' + str(p) + '/strict_np10_leiden_mod.dat', membership)

            partition = strict_consensus(graph, method, n_p=50)
            df, df_coverage, acc_df, membership = calculate_stats_er(graph, list(partition), df_coverage, n, m, p, method,
                                                         'Strict(np=50,Leiden-MOD)', ground_truth_membership, acc_df, df)
            write_membership_to_file('erdos_renyi_lfr/p_' + str(p) + '/strict_np50_leiden_mod.dat', membership)

    df.to_csv('erdos_renyi_lfr_exps_leiden_mod.csv')
    df_coverage.to_csv('erdos_renyi_lfr_exps_leiden_mod_coverage.csv')
    acc_df.to_csv('erdos_renyi_lfr_exps_leiden_acc.csv')


if __name__ == "__main__":
    erdos_renyi()
    #er_lfr()
    #er_ring()








