# Datasets for FastEnsemble paper

This repository contains the datasets and scripts used in the following paper:

Y. Tabatabaee, E. Wedell, M. Park, T. Warnow (2025). *FastEnsemble: Scalable ensemble clustering on large networks*. PLOS Complex Systems 2(10): e0000069 [preliminary version appeared at International Conference on Complex Networks and their Applications (CNA) 2024] DOI: [10.1371/journal.pcsy.0000069](https://journals.plos.org/complexsystems/article?id=10.1371/journal.pcsy.0000069)

For experiments in this study, we generated a collection of artifical networks such as ring of cliques, Erdos-Renyi (ER) graphs, and combinations of ER graphs with LFR graphs and ring of cliques. All these datasets were generated using [NetworkX](https://networkx.org). Additionally, we used a collection of 27 synthetic LFR graphs from [Park et. al. (2024)](https://link.springer.com/chapter/10.1007/978-3-031-53499-7_1), that were generated based on the properties of a collection of real networks and their Leiden clusterings with different resolutions. These datasets are available at [Illinois Data Bank](https://databank.illinois.edu/datasets/IDB-6271968). 

This repository includes the new datasets generated for this study, as well as the output of different clustering methods in each experiment.

## Dataset Description

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
