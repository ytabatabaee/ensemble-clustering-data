#!/bin/bash

RES=0.5
THRESH=1.0
NP=100
NET_TYPE=lfr_training/default
CLUSTERING=fec_leiden_mod_${THRESH} #fastconsensus_t0.2_d0.02_np10_louvain.tsv #fec_leiden_mod_0.5


#for n in 90 100 500 1000 5000 10000; 
for n in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9;
#for p in 0.1 0.01 0.001 0.02 0.002 0.05 0.005;
do
   #python3 clustering_accuracy_nmi_ari.py -gt ../${NET_TYPE}/p_${p}/community.dat -p ../${NET_TYPE}/p_${p}/${CLUSTERING}.dat > ../${NET_TYPE}/p_${p}/${CLUSTERING}.acc
   python3 clustering_accuracy_nmi_ari.py -gt ../${NET_TYPE}/mu_${n}/community.dat -p ../${NET_TYPE}/mu_${n}/${CLUSTERING}.dat > ../${NET_TYPE}/mu_${n}/${CLUSTERING}.acc
   #python3 clustering_accuracy_nmi_ari.py -gt ../${NET_TYPE}/n_${n}/community.dat -p ../${NET_TYPE}/n_${n}/${CLUSTERING}.dat > ../${NET_TYPE}/n_${n}/${CLUSTERING}.acc   
   #python3 fast_ensemble.py -n ../${NET_TYPE}/mu_${n}/network.dat -t ${THRESH} -alg leiden-mod -o ../${NET_TYPE}/mu_${n}/${CLUSTERING}.dat
done
