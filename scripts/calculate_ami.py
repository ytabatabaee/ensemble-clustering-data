import multiprocessing as mp
import subprocess
import sys
import os


def run(comb):
    net, res, clustering = comb

    cmd = 'python3 clustering_accuracy_nmi_ari.py -gt ' + 'CM_paper_LFR_data/'+net+'_'+res+'_lfr/community.dat' +\
          '-p ' + 'CM_paper_LFR_data/'+net+'_'+res+'_lfr/'+clustering+'.dat' + '> ' + 'CM_paper_LFR_data/'+net+'_'+res+'_lfr/'+clustering+'.acc'

    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    _, _ = p.communicate()


if __name__ == '__main__':
    combinations = [("wiki_topcats", "0.1", "fec_0.8"),
                    ("wiki_topcats", "0.1", "leiden_0.1"),
                    ("cit_patents", "0.01", "fec_0.8"),
                    ("cit_patents", "0.1", "fec_0.8"),
                    ("cit_patents", "0.5", "fec_0.8"),
                    ("cit_patents", "0.001", "leiden_0.001"),
                    ("cit_patents", "0.01", "leiden_0.01"),
                    ("cit_patents", "0.1", "leiden_0.1"),
                    ("cit_patents", "0.5", "leiden_0.5"),
                    ("oc", "0.0001", "fec_0.8"),
                    ("oc", "0.001", "fec_0.8"),
                    ("cen", "0.01", "leiden_0.01"),
                    ("cen", "0.001", "fec_0.8"),
                    ("cen", "0.01", "fec_0.8")]
    with mp.Pool(mp.cpu_count()) as p:
        p.map(run, combinations)
