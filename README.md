# Datasets for FastEnsemble paper

This repository contains the datasets and scripts used in the following paper:

Y. Tabatabaee, E. Wedell, M. Park, T. Warnow (2025). *FastEnsemble: Scalable ensemble clustering on large networks*. PLOS Complex Systems 2(10): e0000069 [preliminary version appeared at International Conference on Complex Networks and their Applications (CNA) 2024] DOI: [10.1371/journal.pcsy.0000069](https://journals.plos.org/complexsystems/article?id=10.1371/journal.pcsy.0000069)

For experiments in this study, we generated a collection of artifical networks such as ring of cliques, Erdos-Renyi (ER) graphs, and combinations of ER graphs with LFR graphs and ring of cliques. All these datasets were generated using [NetworkX](https://networkx.org). Additionally, we used a collection of 27 synthetic LFR graphs from [Park et. al. (2024)](https://link.springer.com/chapter/10.1007/978-3-031-53499-7_1), that were generated based on the properties of a collection of real networks and their Leiden clusterings with different resolutions. These datasets are available at [Illinois Data Bank](https://databank.illinois.edu/datasets/IDB-6271968). 

This repository includes the new datasets generated for this study, as well as the output of different clustering methods in each experiment.

## Dataset Description

Each directory includes the datasets and results for one set of networks used in this study. Each subdirectory is a model condition for that dataset. Below is a description of each directory:
- [lfr_training/](https://github.com/ytabatabaee/ensemble-clustering-data/tree/main/lfr_training): LFR algorithm design datasets with mixing parameters varying between `0.1` and `0.9`. The `default` model condition has 10,000 nodes and average degree of 10. The model conditions named as `d_[AVG-DEG]` have 10,000 nodes with average degree of `[AVG-DEG]` (5 or 20) and the model conditions named as `n_[NUM-NODES]` have average degree of 10 with `[NUM-NODES]` nodes (1000 or 100,000).
- [erdos_renyi/](https://github.com/ytabatabaee/ensemble-clustering-data/tree/main/erdos_renyi): Erdos-Renyi networks with 1000 nodes and density (`p`) varying between `0.001` and `0.1`.
- [erdos_renyi_lfr/](https://github.com/ytabatabaee/ensemble-clustering-data/tree/main/erdos_renyi_lfr): Erdos-Renyi network of size 1000 with density (`p`) varying between `0.001` and `0.1` connected to an LFR graph of size 1000.
- [erdos_renyi_ring/](https://github.com/ytabatabaee/ensemble-clustering-data/tree/main/erdos_renyi_ring): Erdos-Renyi network of size 1000 with density (`p`) varying between `0.001` and `0.1` connected to a Ring-of-Cliques network with 100 cliques of size 10.
- [tandon_et_al/](https://github.com/ytabatabaee/ensemble-clustering-data/tree/main/tandon_et_al): Reproduction of the 10,000 node LFR datasets from Tandon et al. (2019).
- [tree_mod/](https://github.com/ytabatabaee/ensemble-clustering-data/tree/main/tree_mod): Tree-of-Cliques networks with number of nodes (`n`) varying between 90 and 5000 and cliques of size 10 used in modularity experiments.
- [ring_mod/](https://github.com/ytabatabaee/ensemble-clustering-data/tree/main/ring_mod): Ring-of-Cliques networks with number of nodes (`n`) varying between 90 and 10,000 and cliques of size 10 used in modularity experiments.
- [ring_cpm/](https://github.com/ytabatabaee/ensemble-clustering-data/tree/main/ring_cpm) and [ring_cpm_res/](https://github.com/ytabatabaee/ensemble-clustering-data/tree/main/ring_cpm_res): Ring-of-Cliques networks with number of nodes (`n`) varying between 90 and 10,000 and cliques of size 10 used in the CPM experiments.

Below is the description of files in each directory:
- `network.dat`: network edge-list
- `community.dat`: ground-truth community structure in the form of `node`:`membership`
- `ecg.dat`: result of clustering with ECG
- `mu_dist.csv`: distribution of mixing parameter values (`mu`) for the network/ground-truth communities
- `leiden_mod.dat`: result of clustering with Leiden-modularity
- `fec_leiden_[PARAM].dat`: result of clustering with FastEnsemble with Leiden-modularity or Leiden-CPM
- `fec_nw_leiden_[PARAM].dat`: result of clustering with FastEnsemble with Leiden-modularity or Leiden-CPM without weighting
- `original_fastconsensus_louvain.clustering`: result of clustering with FastConsensus
- `strict_np[NUM-E]_leiden_mod.dat`: result of clustering with strict consensus with [NUM-E] ensembles
- `louvain.dat`: result of clustering with the Louvain algorithm
