import random
import sys

import click
import numpy as np

@click.command()
@click.option("--fast-consensus-file", required=True, type=click.Path(exists=True), help="Input line per cluster format file")
@click.option("--output-prefix", required=True, type=click.Path(), help="Output file prefix")
def convert_fast_consensus_to_tsv(fast_consensus_file, output_prefix):
    """This script converts mcl label format to a two column tsv format
    """
    with open(output_prefix, "w") as f:
        with open(fast_consensus_file, "r") as fr:
            cluster_id = 0
            for line in fr:
                cluster_member_arr = line.strip().split()
                for cluster_member in cluster_member_arr:
                    f.write(f"{cluster_member}\t{cluster_id}\n")
                cluster_id += 1

if __name__ == "__main__":
    convert_fast_consensus_to_tsv()
