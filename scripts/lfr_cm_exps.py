import leidenalg
import networkx as nx
import igraph as ig
#import community.community_louvain as cl
import numpy as np
import pandas as pd
import argparse
import seaborn as sns
from networkx.algorithms.community import modularity
from sklearn.metrics.cluster import adjusted_mutual_info_score
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.cluster import normalized_mutual_info_score
import matplotlib.pyplot as plt
import csv


def partition_statistics(G, partition, show_cluster_size_dist=True):
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
    """
    Runs an ensemble of single-level randomized Louvain;
    each member of the ensemble gets a "vote" to determine if the edges
    are intra-community or not;
    the votes are aggregated into an "ECG edge weights" in range [0,1];
    a final (full depth) Louvain (using the louvain package) is run
    using those edge weights;

    Parameters
    ----------
    self : :class:`igraph.Graph`
      Graph to define the partition on.

    weights: list of double, optional
      the edge weights

    ens_size: int
      the size of the ensemble of single-level Louvain.

    min_weight: double in range [0,1]
      the ECG edge weight for edges with zero votes from the ensemble.

    Returns
    -------
    partition
      The optimised partition, of class `igraph.clustering.VertexClustering`.

    partition.W
      The ECG edge weights

    partition.CSI
      The community strength index

    Notes
    -----
    The ECG edge weight function is defined as:

      min_weight + ( 1 - min_weight ) x (#votes_in_ensemble) / ens_size

    The weights are linear in terms of the #votes, in the range [min_weight,1].


    Examples
    --------
    >>> g = igraph.Graph.Famous('Zachary')
    >>> part = g.community_ecg(ens_size=25, min_weight = .1)
    """
    W = [0] * self.ecount()
    ## Ensemble of level-1 Louvain
    for i in range(ens_size):
        p = np.random.permutation(self.vcount()).tolist()
        g = self.permute_vertices(p)
        l1 = g.community_multilevel(weights=weights, return_levels=True)[0].membership
        b = [l1[p[x.tuple[0]]] == l1[p[x.tuple[1]]] for x in self.es]
        W = [W[i] + b[i] for i in range(len(W))]
    W = [min_weight + (1 - min_weight) * W[i] / ens_size for i in range(len(W))]
    ## Force min_weight outside 2-core
    core = self.shell_index()
    ecore = [min(core[x.tuple[0]], core[x.tuple[1]]) for x in self.es]
    w = [W[i] if ecore[i] > 1 else min_weight for i in range(len(ecore))]
    part = self.community_multilevel(weights=w)
    part._modularity_params['weights'] = weights
    part.recalculate_modularity()
    part.W = w
    part.CSI = 1 - 2 * np.sum([min(1 - i, i) for i in w]) / len(w)
    return part


ig.Graph.community_ecg = community_ecg


def get_communities(graph, algorithm, seed, res_val=0.01):
    #if algorithm == 'louvain':
    #    return cl.best_partition(graph, random_state=seed, weight='weight')
    if algorithm == 'ecg':
        relabelled_graph = ig.Graph.from_networkx(graph)
        networkx_node_id_dict = {}
        igraph_node_id_dict = relabelled_graph.community_ecg(ens_size=10)
        print(igraph_node_id_dict)
        for igraph_index, vertex in enumerate(relabelled_graph.vs):
            vertex_attributes = vertex.attributes()
            original_id = int(vertex_attributes["_nx_name"])
            relabelled_id = int(igraph_index)
            networkx_node_id_dict[original_id] = igraph_node_id_dict[relabelled_id]
        return networkx_node_id_dict
    if algorithm == 'leiden-cpm':
        relabelled_graph = ig.Graph.from_networkx(graph)
        networkx_node_id_dict = {}
        igraph_node_id_dict = leidenalg.find_partition(relabelled_graph, leidenalg.CPMVertexPartition,
                                                       resolution_parameter=res_val, n_iterations=1, seed=seed).membership
        print(igraph_node_id_dict)
        for igraph_index, vertex in enumerate(relabelled_graph.vs):
            vertex_attributes = vertex.attributes()
            original_id = int(vertex_attributes["_nx_name"])
            relabelled_id = int(igraph_index)
            networkx_node_id_dict[original_id] = igraph_node_id_dict[relabelled_id]
        return networkx_node_id_dict
        #return communities_to_dict(leidenalg.find_partition(ig.Graph.from_networkx(graph),
        #                                          leidenalg.CPMVertexPartition,
        #                                          resolution_parameter=res_val,
        #                                          n_iterations=2).as_cover())
    elif algorithm == 'leiden-mod':
        relabelled_graph = ig.Graph.from_networkx(graph)
        networkx_node_id_dict = {}
        igraph_node_id_dict = leidenalg.find_partition(relabelled_graph, leidenalg.ModularityVertexPartition,
                                                       weights='weight', n_iterations=-1, seed=seed).membership
        for igraph_index, vertex in enumerate(relabelled_graph.vs):
            vertex_attributes = vertex.attributes()
            original_id = int(vertex_attributes["_nx_name"])
            relabelled_id = int(igraph_index)
            networkx_node_id_dict[original_id] = igraph_node_id_dict[relabelled_id]
        return networkx_node_id_dict
        #return communities_to_dict(leidenalg.find_partition(ig.Graph.from_networkx(graph),
        #                                leidenalg.ModularityVertexPartition,
        #                                weights='weight',
        #                                seed=seed).as_cover())


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


def simple_consensus(G, algorithm='leiden-mod', n_p=10, thresh=0.9, delta=0.02, max_iter=10, res_value=0.01):
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
    final_comm = get_communities(graph, algorithm, 0, res_val=res_value)
    #print(final_comm)
    #final_comm = group_to_partition(final_comm)
    #print(final_comm)
    return final_comm


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

    return get_communities(graph, algorithm, 0, res_val) # group_to_partition(


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
    precision = tp / (tp + fp)
    recall = tp / (tp + fn)
    f1_score = 2 * precision * recall / (precision + recall)
    fnr = fn / (fn + tp)
    fpr = fp / (fp + tn)'''

    nmi = normalized_mutual_info_score(mem_true, mem_est)
    ari = adjusted_rand_score(mem_true, mem_est)
    #ami = adjusted_mutual_info_score(mem_true, mem_est)

    #print(tp, tn, fp, fn)

    #return nmi, ami, ari, precision, recall, f1_score, fnr, fpr
    return nmi, 0, ari, 0, 0, 0, 0, 0


def threshold_consensus(G, algorithm='leiden-cpm', n_p=20, tr=0.2, res_value=0.01):
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

    #return group_to_partition(get_communities(graph, algorithm, 0, res_val=res_value))
    return get_communities(graph, algorithm, 0, res_val=res_value)


def lfr_exp_mu():
    acc_df = pd.DataFrame(columns=["mu", "partition", "acc_measure", "acc_value"])
    method = 'leiden-mod'

    for mu in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
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

        print('mu:', mu)
        partition = normal_clustering(graph, method)
        membership = get_membership_list_from_dict(partition, 10000)
        nmi, ami, ari, precision, recall, f1_score, fnr, fpr = measure_accuracy(ground_truth_membership,
                                                                                membership)
        print("Normalized mutual information (NMI): ", nmi)
        print("Adjusted rand index (ARI): ", ari)
        print("Adjusted mutual information (AMI): ", ami)
        print("False positive rate (FPR), False negative rate (FNR):", fpr, fnr)
        print("Precision, Recall, F1-score:", precision, recall, f1_score)
        acc_df.loc[len(acc_df.index)] = [mu, 'Leiden-MOD', 'NMI', nmi]
        acc_df.loc[len(acc_df.index)] = [mu, 'Leiden-MOD', 'AMI', ami]
        acc_df.loc[len(acc_df.index)] = [mu, 'Leiden-MOD', 'ARI', ari]
        acc_df.loc[len(acc_df.index)] = [mu, 'Leiden-MOD', 'FPR', fpr]
        acc_df.loc[len(acc_df.index)] = [mu, 'Leiden-MOD', 'FNR', fnr]
        acc_df.loc[len(acc_df.index)] = [mu, 'Leiden-MOD', 'Precision', precision]
        acc_df.loc[len(acc_df.index)] = [mu, 'Leiden-MOD', 'Recall', recall]
        acc_df.loc[len(acc_df.index)] = [mu, 'Leiden-MOD', 'F1-Score', f1_score]

        partition = threshold_consensus(graph, method, n_p=10, tr=0.2)
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
        acc_df.loc[len(acc_df.index)] = [mu, 'TC(np=10,tr=0.2)+Leiden-MOD', 'NMI', nmi]
        acc_df.loc[len(acc_df.index)] = [mu, 'TC(np=10,tr=0.2)+Leiden-MOD', 'AMI', ami]
        acc_df.loc[len(acc_df.index)] = [mu, 'TC(np=10,tr=0.2)+Leiden-MOD', 'ARI', ari]
        acc_df.loc[len(acc_df.index)] = [mu, 'TC(np=10,tr=0.2)+Leiden-MOD', 'FPR', fpr]
        acc_df.loc[len(acc_df.index)] = [mu, 'TC(np=10,tr=0.2)+Leiden-MOD', 'FNR', fnr]
        acc_df.loc[len(acc_df.index)] = [mu, 'TC(np=10,tr=0.2)+Leiden-MOD', 'Precision', precision]
        acc_df.loc[len(acc_df.index)] = [mu, 'TC(np=10,tr=0.2)+Leiden-MOD', 'Recall', recall]
        acc_df.loc[len(acc_df.index)] = [mu, 'TC(np=10,tr=0.2)+Leiden-MOD', 'F1-Score', f1_score]

        partition = threshold_consensus(graph, method, n_p=10, tr=0.5)
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
        acc_df.loc[len(acc_df.index)] = [mu, 'TC(np=10,tr=0.5)+Leiden-MOD', 'NMI', nmi]
        acc_df.loc[len(acc_df.index)] = [mu, 'TC(np=10,tr=0.5)+Leiden-MOD', 'AMI', ami]
        acc_df.loc[len(acc_df.index)] = [mu, 'TC(np=10,tr=0.5)+Leiden-MOD', 'ARI', ari]
        acc_df.loc[len(acc_df.index)] = [mu, 'TC(np=10,tr=0.5)+Leiden-MOD', 'FPR', fpr]
        acc_df.loc[len(acc_df.index)] = [mu, 'TC(np=10,tr=0.5)+Leiden-MOD', 'FNR', fnr]
        acc_df.loc[len(acc_df.index)] = [mu, 'TC(np=10,tr=0.5)+Leiden-MOD', 'Precision', precision]
        acc_df.loc[len(acc_df.index)] = [mu, 'TC(np=10,tr=0.5)+Leiden-MOD', 'Recall', recall]
        acc_df.loc[len(acc_df.index)] = [mu, 'TC(np=10,tr=0.5)+Leiden-MOD', 'F1-Score', f1_score]

        partition = threshold_consensus(graph, method, n_p=10, tr=0.8)
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
        acc_df.loc[len(acc_df.index)] = [mu, 'TC(np=10,tr=0.8)+Leiden-MOD', 'NMI', nmi]
        acc_df.loc[len(acc_df.index)] = [mu, 'TC(np=10,tr=0.8)+Leiden-MOD', 'AMI', ami]
        acc_df.loc[len(acc_df.index)] = [mu, 'TC(np=10,tr=0.8)+Leiden-MOD', 'ARI', ari]
        acc_df.loc[len(acc_df.index)] = [mu, 'TC(np=10,tr=0.8)+Leiden-MOD', 'FPR', fpr]
        acc_df.loc[len(acc_df.index)] = [mu, 'TC(np=10,tr=0.8)+Leiden-MOD', 'FNR', fnr]
        acc_df.loc[len(acc_df.index)] = [mu, 'TC(np=10,tr=0.8)+Leiden-MOD', 'Precision', precision]
        acc_df.loc[len(acc_df.index)] = [mu, 'TC(np=10,tr=0.8)+Leiden-MOD', 'Recall', recall]
        acc_df.loc[len(acc_df.index)] = [mu, 'TC(np=10,tr=0.8)+Leiden-MOD', 'F1-Score', f1_score]

    acc_df.to_csv('res_limit_exps_leiden_mod_lfr_mu_acc.csv')

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

        partition = threshold_consensus(graph, method, n_p=10, tr=tr)
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
    return df, acc_df


def ring_exp_cpm():
    df = pd.DataFrame(columns=['k', "n", "res", "method", "partition", "cluster_id", "cluster_size"])
    acc_df = pd.DataFrame(columns=['k', "n", "res", "method", "partition", "acc_measure", "acc_value"])
    n = 1000
    for res in [0.00001, 0.0001, 0.001, 0.01, 0.1, 0.5]:
    #for n in [90, 100, 500, 1000, 5000]:#, 10000]:
        print('n:', n)
        for k in [10]:
            graph = nx.ring_of_cliques(num_cliques=n, clique_size=k)
            ground_truth_membership = sum([[i]*k for i in range(n)], [])
            #graph = gen_tree_of_cliques(k, n)
            for method in ['leiden-cpm']:
                print('original')
                partition = normal_clustering(graph, method, res_val=res)
                df, acc_df = calculate_stats(graph, partition, n, k, res, method,
                                                             'Leiden-CPM', ground_truth_membership, acc_df, df)

                print('SC n = 10')
                partition = strict_consensus(graph, method, n_p=10, res_val=res)
                df, acc_df = calculate_stats(graph, list(partition), n, k, res, method,
                                             'Strict(np=10)+Leiden-CPM', ground_truth_membership, acc_df, df)

                print('SC n = 50')
                partition = strict_consensus(graph, method, n_p=50, res_val=res)
                df, acc_df = calculate_stats(graph, list(partition), n, k, res, method,
                                             'Strict(np=50)+Leiden-CPM', ground_truth_membership, acc_df, df)

    df.to_csv('res_limit_exps_leiden_cpm_vary_res.csv')
    acc_df.to_csv('res_limit_exps_leiden_cpm_vary_res_acc.csv')


def ring_exp():
    # method is Leiden, Louvain, Leiden-mod
    # partition is SC, original, ground-truth
    #df = pd.DataFrame(columns=['k', "n", "method", "partition", "cluster_id", "cluster_size"])
    # acc_df = pd.DataFrame(columns=['k', "n", "method", "partition", "nmi", "ami", "ari", "fpr", "fnr", "precision", "recall", "f1_score"])
    acc_df = pd.DataFrame(columns=['k', "n", "method", "partition", "acc_measure", "acc_value"])
    #for n in [90, 100, 500, 1000, 5000, 10000]:
    #n = 1000
    #for res in [0.00001, 0.0001, 0.001, 0.01, 0.1, 0.5]:
    for n in [90, 100, 500, 1000, 5000]:#, 10000]:
        print('n:', n)
        for k in [10]:
            graph = nx.ring_of_cliques(num_cliques=n, clique_size=k)
            ground_truth_membership = sum([[i]*k for i in range(n)], [])
            #graph = gen_tree_of_cliques(k, n)
            for method in ['leiden-mod']: # , 'leiden-mod', 'louvain'
                print('original')
                partition = normal_clustering(graph, method, res_val=0.0001)
                membership = get_membership_list_from_dict(partition, n * k)
                nmi, ami, ari, precision, recall, f1_score, fnr, fpr = measure_accuracy(ground_truth_membership, membership)
                print("Normalized mutual information (NMI): ", nmi)
                print("Adjusted rand index (ARI): ", ari)
                print("Adjusted mutual information (AMI): ", ami)
                print("False positive rate (FPR), False negative rate (FNR):", fpr, fnr)
                print("Precision, Recall, F1-score:", precision, recall, f1_score)
                acc_df.loc[len(acc_df.index)] = [k, n, method, 'Leiden-MOD', 'NMI', nmi]
                acc_df.loc[len(acc_df.index)] = [k, n, method, 'Leiden-MOD', 'AMI', ami]
                acc_df.loc[len(acc_df.index)] = [k, n, method, 'Leiden-MOD', 'ARI', ari]
                acc_df.loc[len(acc_df.index)] = [k, n, method, 'Leiden-MOD', 'FPR', fpr]
                acc_df.loc[len(acc_df.index)] = [k, n, method, 'Leiden-MOD', 'FNR', fnr]
                acc_df.loc[len(acc_df.index)] = [k, n, method, 'Leiden-MOD', 'Precision', precision]
                acc_df.loc[len(acc_df.index)] = [k, n, method, 'Leiden-MOD', 'Recall', recall]
                acc_df.loc[len(acc_df.index)] = [k, n, method, 'Leiden-MOD', 'F1-Score', f1_score]

                #cluster_sizes = partition_statistics(graph, partition)
                #for i in range(len(cluster_sizes)):
                #    df.loc[len(df.index)] = [k, n, method, 'Leiden-MOD', i, cluster_sizes[i]]

                print('SC n = 10')
                partition = strict_consensus(graph, method, n_p=10, res_val=0.0001)
                membership = list(partition.values())
                # print(membership)
                #partition = group_to_partition(partition)
                nmi, ami, ari, precision, recall, f1_score, fnr, fpr = measure_accuracy(ground_truth_membership, membership)
                print("Normalized mutual information (NMI): ", nmi)
                print("Adjusted rand index (ARI): ", ari)
                print("Adjusted mutual information (AMI): ", ami)
                print("False positive rate (FPR), False negative rate (FNR):", fpr, fnr)
                print("Precision, Recall, F1-score:", precision, recall, f1_score)
                acc_df.loc[len(acc_df.index)] = [k, n, method, 'SC(np=10)+Leiden-MOD', 'NMI', nmi]
                acc_df.loc[len(acc_df.index)] = [k, n, method, 'SC(np=10)+Leiden-MOD', 'AMI', ami]
                acc_df.loc[len(acc_df.index)] = [k, n, method, 'SC(np=10)+Leiden-MOD', 'ARI', ari]
                acc_df.loc[len(acc_df.index)] = [k, n, method, 'SC(np=10)+Leiden-MOD', 'FPR', fpr]
                acc_df.loc[len(acc_df.index)] = [k, n, method, 'SC(np=10)+Leiden-MOD', 'FNR', fnr]
                acc_df.loc[len(acc_df.index)] = [k, n, method, 'SC(np=10)+Leiden-MOD', 'Precision', precision]
                acc_df.loc[len(acc_df.index)] = [k, n, method, 'SC(np=10)+Leiden-MOD', 'Recall', recall]
                acc_df.loc[len(acc_df.index)] = [k, n, method, 'SC(np=10)+Leiden-MOD', 'F1-Score', f1_score]
                #acc_df.loc[len(acc_df.index)] = [k, n, method, 'SC(np=10)+Leiden-MOD', nmi, ami, ari, fpr, fnr, precision, recall, f1_score]
                #cluster_sizes = partition_statistics(graph, partition)
                #for i in range(len(cluster_sizes)):
                #    df.loc[len(df.index)] = [k, n, method, 'SC(np=10)+Leiden-MOD', i, cluster_sizes[i]]

                print('SC n = 50')
                partition = strict_consensus(graph, method, n_p=50, res_val=0.0001)
                membership = list(partition.values())
                #partition = group_to_partition(partition)
                nmi, ami, ari, precision, recall, f1_score, fnr, fpr = measure_accuracy(ground_truth_membership, membership)
                print("Normalized mutual information (NMI): ", nmi)
                print("Adjusted rand index (ARI): ", ari)
                print("Adjusted mutual information (AMI): ", ami)
                print("False positive rate (FPR), False negative rate (FNR):", fpr, fnr)
                print("Precision, Recall, F1-score:", precision, recall, f1_score)
                #acc_df.loc[len(acc_df.index)] = [k, n, method, 'SC(np=50)+Leiden-MOD', nmi, ami, ari, fpr, fnr, precision, recall, f1_score]
                acc_df.loc[len(acc_df.index)] = [k, n, method, 'SC(np=50)+Leiden-MOD', 'NMI', nmi]
                acc_df.loc[len(acc_df.index)] = [k, n, method, 'SC(np=50)+Leiden-MOD', 'AMI', ami]
                acc_df.loc[len(acc_df.index)] = [k, n, method, 'SC(np=50)+Leiden-MOD', 'ARI', ari]
                acc_df.loc[len(acc_df.index)] = [k, n, method, 'SC(np=50)+Leiden-MOD', 'FPR', fpr]
                acc_df.loc[len(acc_df.index)] = [k, n, method, 'SC(np=50)+Leiden-MOD', 'FNR', fnr]
                acc_df.loc[len(acc_df.index)] = [k, n, method, 'SC(np=50)+Leiden-MOD', 'Precision', precision]
                acc_df.loc[len(acc_df.index)] = [k, n, method, 'SC(np=50)+Leiden-MOD', 'Recall', recall]
                acc_df.loc[len(acc_df.index)] = [k, n, method, 'SC(np=50)+Leiden-MOD', 'F1-Score', f1_score]
                #cluster_sizes = partition_statistics(graph, partition)
                #for i in range(len(cluster_sizes)):
               #     df.loc[len(df.index)] = [k, n, method, 'SC(np=50)+Leiden-MOD', i, cluster_sizes[i]]

                #partition = strict_consensus(ring, method, n_p=100, res_val=0.01)
                #cluster_sizes = partition_statistics(ring, partition)
                #for i in range(len(cluster_sizes)):
                #    df.loc[len(df.index)] = [k, n, method, 'SC(np=100)+Leiden-CPM(r=0.01)', i, cluster_sizes[i]]
    #df.to_csv('res_limit_exps_leiden_cpm_tree.csv')
    acc_df.to_csv('res_limit_exps_leiden_mod_ring_acc.csv')


if __name__ == "__main__":
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
    #sc = simple_consensus(net, algorithm='leiden-mod', n_p=args.partitions, thresh=args.threshold, max_iter=3)
    #sc = threshold_consensus(net, algorithm='leiden-mod', n_p=args.partitions, tr=args.threshold)
    if args.algorithm == 'leiden-cpm':
        #sc = simple_consensus(net, algorithm='leiden-cpm', n_p=args.partitions, thresh=args.threshold, max_iter=3, res_value=args.resolution)
        sc = threshold_consensus(net, algorithm='leiden-mod', n_p=args.partitions, tr=args.threshold, res_value=args.resolution)
    elif args.algorithm == 'leiden-mod':
        sc = threshold_consensus(net, algorithm='leiden-mod', n_p=args.partitions, tr=args.threshold)
        #sc = simple_consensus(net, algorithm='leiden-mod', n_p=args.partitions, thresh=args.threshold, max_iter=3)

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
    print("Precision, Recall, F1-score:", precision, recall, f1_score)




