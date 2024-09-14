import argparse
import os
import random
import community as cm
import igraph as ig
import leidenalg
import networkx as nx
import numpy as np
from networkx.algorithms.community import modularity, partition_quality
from sklearn.metrics.cluster import adjusted_mutual_info_score
from itertools import combinations
import matplotlib.pyplot as plt


def get_membership_list(membership_dict):
    memberships = []
    for i in range(len(membership_dict)):
        memberships.append(membership_dict[i])
    return memberships


def convert_leiden_output_to_dict(leiden_membership):
    membership_dict = {}
    for i in range(len(leiden_membership)):
        membership_dict[i] = leiden_membership[i][0]
    return membership_dict


def check_convergence(G, n_p, delta):
    count = 0
    for wt in nx.get_edge_attributes(G, 'weight').values():
        if wt != 0 and wt != n_p:
            count += 1
    if count > delta * G.number_of_edges():
        return False
    return True


def group_to_partition(partition):
    '''
    Takes in a partition, dictionary in the format {node: community_membership}
    Returns a nested list of communities [[comm1], [comm2], ...... [comm_n]]
    '''
    part_dict = {}
    for index, value in partition.items():
        if value in part_dict:
            part_dict[value].append(index)
        else:
            part_dict[value] = [index]
    return part_dict.values()


def check_arguments(args):
    if args.d > 1:
        print('delta is too high. Allowed values are between 0 and 1')
        return False
    if args.d < 0:
        print('delta is too low. Allowed values are between 0 and 1')
        return False
    if args.alg not in ('louvain', 'lpm', 'cnm', 'infomap', 'leiden'):
        print('Incorrect algorithm entered. run with -h for help')
        return False
    if args.t > 1 or args.t < 0:
        print('Incorrect tau. run with -h for help')
        return False
    return True


def communities_to_dict(communities):
    """
    Creates a [node] -> [community] lookup
    """
    result = {}
    community_index = 0
    for c in communities:
        community_mapping = ({node: community_index for index, node in enumerate(c)})
        result = {**result, **community_mapping}
        community_index += 1
    return result


def do_leiden_community_detection(data, gamma):
    networkx_graph, seed = data
    '''return leidenalg.find_partition(ig.Graph.from_networkx(networkx_graph),
                                        leidenalg.CPMVertexPartition,
                                        resolution_parameter=gamma,
                                        weights='weight',
                                        seed=seed, n_iterations=1).as_cover()'''
    return leidenalg.find_partition(ig.Graph.from_networkx(networkx_graph),
                                    leidenalg.ModularityVertexPartition,
                                    weights='weight',
                                    seed=seed, n_iterations=1).as_cover()


def get_graph_and_seed(graph, times):
    for seed in range(times):
        yield graph, seed


def thresholding(graph, n_p, thresh):
    """remove edges with weight less than thresh*n_p from the graph"""
    remove_edges = []
    for u, v in graph.edges():
        if graph[u][v]['weight'] < thresh * n_p:
            remove_edges.append((u, v))
    graph.remove_edges_from(remove_edges)
    return graph


def initialize(graph, value):
    """initialize all edges weights in graph with given constant value"""
    for u, v in graph.edges():
        graph[u][v]['weight'] = value
    return graph


def connect_singletons(graph, nextgraph):
    """keep the graph connected by adding back singletons with their maximum weight edges"""
    for node in nx.isolates(nextgraph):
        nbr, weight = sorted(graph[node].items(), key=lambda edge: edge[1]['weight'])[0]
        nextgraph.add_edge(node, nbr, weight=weight['weight'])
    return nextgraph


def triadic_closure(graph, L, n_p, communities, algorithm, mapping=None):
    for _ in range(L):
        node = np.random.choice(graph.nodes())
        neighbors = [a[1] for a in graph.edges(node)]
        if len(neighbors) >= 2:
            a, b = random.sample(set(neighbors), 2)
            if not graph.has_edge(a, b):
                graph.add_edge(a, b, weight=0)
                for i in range(n_p):
                    if algorithm == 'louvain':
                        c = communities[i]
                        if c[a] == c[b]:
                            graph[a][b]['weight'] += 1
                    elif algorithm == 'leiden':
                        node_community_lookup = communities_to_dict(communities[i])
                        if a in node_community_lookup and b in node_community_lookup and \
                                node_community_lookup[a] == \
                                node_community_lookup[b]:
                            graph[a][b]['weight'] += 1
                    elif algorithm in ('infomap', 'lpm'):
                        for c in communities[i]:
                            if a in c and b in c:
                                graph[a][b]['weight'] += 1
                    elif algorithm == 'cnm':
                        for c in communities[i]:
                            if mapping[i][a] in c and mapping[i][b] in c:
                                graph[a][b]['weight'] += 1
    return graph


def get_communities(graph, n_p, N, algorithm, gamma):
    if algorithm == 'louvain':
        return [cm.partition_at_level(cm.generate_dendrogram(graph, randomize=True, weight='weight'), 0)
                   for _ in range(n_p)]
    elif algorithm == 'leiden':
        return [do_leiden_community_detection(data, gamma) for data in get_graph_and_seed(graph, n_p)]
    elif algorithm == 'infomap':
        return [{frozenset(c) for c in ig.Graph.from_networkx(graph).community_infomap().as_cover()} for _ in
                    range(n_p)]
    elif algorithm == 'lpm':
        return [{frozenset(c) for c in ig.Graph.from_networkx(graph).community_label_propagation().as_cover()} for _
                    in range(n_p)]
    elif algorithm == 'cnm':
        return do_cnm_community_detection(graph, n_p, N)


def do_cnm_community_detection(graph, n_p, N):
    communities = []
    mapping = []
    inv_map = []

    for _ in range(n_p):
        order = list(range(N))
        random.shuffle(order)
        maps = dict(zip(range(N), order))

        mapping.append(maps)
        inv_map.append({v: k for k, v in maps.items()})
        G_c = nx.relabel_nodes(graph, mapping=maps, copy=True)
        G_igraph = ig.Graph.from_networkx(G_c)

        communities.append(G_igraph.community_fastgreedy(weights='weight').as_clustering())

    return communities, mapping, inv_map


def fast_consensus(G, algorithm='louvain', n_p=20, thresh=0.2, delta=0.02, gamma=0.01, max_iter=2):
    graph = G.copy()
    graph = initialize(graph, 1.0)
    L = G.number_of_edges()
    N = G.number_of_nodes()
    iter_count = 0

    while True:
        iter_count += 1
        if iter_count > max_iter:
            break
        print("iter ", iter_count)

        nextgraph = graph.copy()
        nextgraph = initialize(nextgraph, 0.0)

        if algorithm == 'louvain':

            communities = get_communities(graph, n_p, N, algorithm, gamma)
            #communities_additional = get_communities(graph, n_p, N, 'leiden', gamma)
            #for i in range(len(communities_additional)):
            #    communities_additional[i] = convert_leiden_output_to_dict(communities_additional[i].membership)

            #print(communities_additional[0])
            #communities = communities + communities_additional
            #for p1, p2 in combinations(communities, 2):
            #    print(adjusted_mutual_info_score(get_membership_list(p1), get_membership_list(p2)))


            for node, nbr in graph.edges():
                if graph[node][nbr]['weight'] not in (0, n_p):
                    for i in range(n_p):
                        c = communities[i]
                        if c[node] == c[nbr]:
                            nextgraph[node][nbr]['weight'] += 1
                else:
                    nextgraph[node][nbr]['weight'] = graph[node][nbr]['weight']

            nextgraph = thresholding(nextgraph, n_p, thresh)
            #print(nx.adjacency_matrix(nextgraph, weight='weight'))
            if check_convergence(nextgraph, n_p=n_p, delta=delta):
                break

            #nextgraph = triadic_closure(nextgraph, L, n_p, algorithm, None)
            nextgraph = connect_singletons(graph, nextgraph)
            graph = nextgraph.copy()
            if check_convergence(nextgraph, n_p=n_p, delta=delta):
                break

        elif algorithm == 'leiden':

            communities = get_communities(graph, n_p, N, algorithm, gamma)

            for i in range(n_p):
                node_community_lookup = communities_to_dict(communities[i])
                for community_index, _ in enumerate(communities[i]):
                    for node, nbr in graph.edges():
                        if graph[node][nbr]['weight'] not in (0, n_p):
                            if node in node_community_lookup and nbr in node_community_lookup and \
                                    node_community_lookup[node] == node_community_lookup[nbr]:
                                if node_community_lookup[node] != community_index:  # only count each community once
                                    continue
                                nextgraph[node][nbr]['weight'] += 1
                        else:
                            nextgraph[node][nbr]['weight'] = graph[node][nbr]['weight']

            #print(nx.adjacency_matrix(nextgraph, weight='weight'))
            nextgraph = thresholding(nextgraph, n_p, thresh)
            #print(nx.adjacency_matrix(nextgraph, weight='weight'))
            if check_convergence(nextgraph, n_p=n_p, delta=delta):
                break

            #nextgraph = triadic_closure(nextgraph, L, n_p, algorithm, None)
            nextgraph = connect_singletons(graph, nextgraph)
            graph = nextgraph.copy()
            if check_convergence(nextgraph, n_p=n_p, delta=delta):
                break

        elif algorithm in ('infomap', 'lpm'):

            communities = get_communities(graph, n_p, N, algorithm, gamma)

            for node, nbr in graph.edges():
                for i in range(n_p):
                    for c in communities[i]:
                        if node in c and nbr in c:
                            if not nextgraph.has_edge(node, nbr):
                                nextgraph.add_edge(node, nbr, weight=0)
                            nextgraph[node][nbr]['weight'] += 1

            nextgraph = thresholding(nextgraph, n_p, thresh)
            nextgraph = triadic_closure(nextgraph, L, n_p, algorithm, None)

            graph = nextgraph.copy()
            if check_convergence(nextgraph, n_p=n_p, delta=delta):
                break

        elif algorithm == 'cnm':

            communities, mapping, inv_map = get_communities(graph, n_p, N, algorithm, gamma)

            for i in range(n_p):
                edge_list = [(mapping[i][j], mapping[i][k]) for j, k in graph.edges()]
                for node, nbr in edge_list:
                    a, b = inv_map[i][node], inv_map[i][nbr]
                    if graph[a][b] not in (0, n_p):
                        for c in communities[i]:
                            if node in c and nbr in c:
                                nextgraph[a][b]['weight'] += 1
                    else:
                        nextgraph[a][b]['weight'] = graph[a][b]['weight']

            nextgraph = thresholding(nextgraph, n_p, thresh)
            #nextgraph = triadic_closure(nextgraph, L, n_p, algorithm, mapping)
            if check_convergence(nextgraph, n_p, delta):
                break
        else:
            break

    return get_communities(graph, n_p, N, algorithm, gamma)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Consensus clustering.')
    parser.add_argument('-f', metavar='filename', type=str, nargs='?', help='file with edgelist')
    parser.add_argument('-np', metavar='n_p', type=int, nargs='?', default=20,
                        help='number of input partitions (default: 20)', required=False)
    parser.add_argument('-t', metavar='tau', type=float, nargs='?', help='threshold for filtering weak edges')
    parser.add_argument('-d', metavar='del', type=float, nargs='?', default=0.02, required=False,
                        help='convergence parameter (default = 0.02). Converges when less than delta proportion of the edges are with wt = 1, 0')
    parser.add_argument('-g', metavar='gamma', type=float, nargs='?', default=0.01, help='resolution parameter',
                        required=False)
    parser.add_argument('--alg', metavar='alg', type=str, nargs='?', default='louvain',
                        help='choose from \'louvain\' , \'cnm\' , \'lpm\' , \'infomap\' ', required=False)

    args = parser.parse_args()

    default_tau = {'louvain': 0.2, 'cnm': 0.7, 'infomap': 0.6, 'lpm': 0.8}
    if args.t is None:
        args.t = default_tau.get(args.alg, 0.2)

    if not check_arguments(args):
        quit()

    G = nx.read_edgelist(args.f, nodetype=int)
    #G = nx.relabel_nodes(G, lambda x: x-1)
    print('#nodes, #edges, #singletons', G.number_of_nodes(), G.number_of_edges(), len(list(nx.isolates(G))))

    # relabel nodes
    mapping = dict(zip(G, range(0, G.number_of_nodes())))
    G = nx.relabel_nodes(G, mapping)

    '''fig = plt.figure(figsize=(20, 20), dpi=80)
    nx.draw(G, with_labels=True, font_size=8, node_color='lightblue')
    fig.savefig("email.pdf")'''

    output = [[]]
    #output[0] = cm.partition_at_level(cm.generate_dendrogram(G, randomize=True), 0)
    '''output[0] = leidenalg.find_partition(ig.Graph.from_networkx(G),
                                    leidenalg.ModularityVertexPartition,
                                    seed=0, n_iterations=1).as_cover()'''
    output = fast_consensus(G, algorithm=args.alg, n_p=args.np, thresh=args.t, delta=args.d, gamma=args.g)

    out_partitions_path = 'out_partitions_t' + str(args.t) + '_d' + str(args.d) + '_np' + str(args.np) + '_' + args.alg
    membership_path = 'memberships_t' + str(args.t) + '_d' + str(args.d) + '_np' + str(args.np) + '_' + args.alg

    if args.alg == 'leiden':
        out_partitions_path += '_r' + str(args.g)
        membership_path += '_r' + str(args.g)

    if not os.path.exists(out_partitions_path):
        os.makedirs(out_partitions_path)

    if not os.path.exists(membership_path):
        os.makedirs(membership_path)

    if args.alg == 'cnm':
        output = output[0]

    if args.alg == 'louvain':
        for i in range(len(output)):
            output[i] = group_to_partition(output[i])

    print('modularity, coverage, performance', modularity(G, output[0]), partition_quality(G, output[0]))
    print('number of clusters', len(output[0]))

    #if args.alg == 'leiden':
    #    print(output[0].membership)

    for i in range(len(output)):
        output_dict = communities_to_dict(output[i])
        with open(membership_path + '/' + str(i), 'w') as f:
            for k, v in sorted(output_dict.items()):
                f.write(str(k) + "\t" + str(v) + '\n')

    for i in range(len(output)):
        with open(out_partitions_path + '/' + str(i), 'w') as f:
            for community in output[i]:
                print(*community, file=f)
