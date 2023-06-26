# This script is for getting full length cluster centers after using
# starcode to cluster partial sequences. Starcode must be run with the --seq-id
# option to generate an output file with truncated cluster centers and a list
# of sequence positions in the input file for each read in the cluster. This script
# takes the sequence positions, gets the full-length sequences from the starcode input, and
# calculates a consensus full-length sequence after ensuring all full-length sequences in the
# cluster are within some threshold of similarity. Optionally it also collapses sequences within a
# cluster that have the same UMI. It then outputs a list of full-length
# sequences and their read counts as a tab-separated text file.

import argparse
import csv
import os
import sys
from tqdm import tqdm
from collections import Counter
from itertools import zip_longest
from umi_tools import UMIClusterer
from utils import read_fastq
import Levenshtein

def most_common(l):
    '''Returns most common element in a list. In the case of a tie, returns the
       most common element that occurs first in the list.
    '''
    return Counter(l).most_common(1)[0][0]


def get_consensus(seqs):
    '''Get consensus sequence by calculating most common base at each position.
       In the case of a tie, the base that occurred first in the list is chosen.
       Sequences of different lengths can be present; the most frequent length is
       used for the final consensus sequence. In the case of a tie, the length of the
       first sequence in the list with one of the most frequent lengths is used.
    '''

    seq_len = most_common([len(seq) for seq in seqs])

    consensus = []
    for bases in zip_longest(*seqs, fillvalue ='N'):
        # If all bases are the same, don't try to find the most common base
        # (assuming most positions have bases that all match, this gives a massive speedup)
        if len(set(bases)) == 1:
            consensus.append(bases[0])
        else:
            consensus.append(most_common(bases))

    return "".join(consensus)[:seq_len]


def get_consensus_and_count(seqs, max_dist):
    '''Get full length consensus sequence and count from list of sequences.
       Any sequences further than Levenshtein distance of max_dist from the consensus
       are not counted and are returned as a list.

       Arguments:
           seqs (list of str): list of sequences to count
           max_dist (int): Maximum distance between sequence and the cluster consensus that is allowed;
                           sequences with larger distances will be discarded.

        Returns:
            Tuple with 3 elements: the consensus sequence, the count, and a set of sequence indices that
            did not match the consensus within the specified max_dist.
    '''
    consensus = get_consensus(seqs)
    not_matching = []
    count = 0
    for i, seq in enumerate(seqs):
        if Levenshtein.distance(consensus, seq) > max_dist:
            not_matching.append(i)
        else:
            count += 1
    return (consensus, count, set(not_matching))


def extract_umi(header, start_distance, umi_len):
    umi = header[-start_distance:-start_distance + umi_len]
    if not set(umi).issubset({'A', 'C', 'G', 'T', 'N'}):
        raise ValueError(f"UMI contains non-ACTGN character: {umi}")
    return umi


def cluster_umis(umis, threshold, clust_method="directional"):
    '''Cluster a list of UMIs using UMI-tools.

    Arguments:

        umis (list of str): List of UMI sequences to be clustered.
        threshold (int): Max distance for clustering.
        clust_method (str): Clustering algorithm for UMI-tools. Default "directional".
                          Can also be "cluster" or "adjacency".
    
    Returns:
        A list of lists with each list containing the UMIs in one cluster.
    '''

    if clust_method not in {'directional', 'cluster', 'adjacency'}:
        raise ValueError(f"clust_method must be 'directional', 'cluster', or 'adjacency'. You provided {clust_method}.")
    
    umi_counter = dict(Counter([umi.encode() for umi in umis]))
    clusterer = UMIClusterer(cluster_method=clust_method)
    clusters = clusterer(umi_counter, threshold=threshold)

    return clusters
    

def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fastq', dest='input_fastq_fp', help="Fastq file used as starcode input", required=True)
    parser.add_argument('-c', '--clustered', dest='clustered_fp', help="Starcode output file (must use --seq-id option with starcode)", required=True)
    parser.add_argument('-o', '--output', dest='output_fp', help="Output file path", required=True)
    parser.add_argument('-l', '--log', dest='log_fp', help="Log file path", required=True)
    parser.add_argument('-d', '--max_dist', dest='max_dist', type=int, default=10, help="Maximum distance between sequence and cluster consensus")
    parser.add_argument('--umi', dest='collapse_umi', action='store_true', default=False, help="Collapse sequences using UMIs")
    parser.add_argument('--umi-threshold', default=1, type=int, help="Max Levenshtein distance for collapsing UMIs")
    parser.add_argument('--umi-start', type=int, help="Distance of UMI start index from the end of the header line (len(header) - index of UMI start)", required='--umi' in sys.argv)
    parser.add_argument('--umi-len', type=int, help="Length of UMI", required='--umi' in sys.argv)
    return parser.parse_args(args)


def main(unparsed_args):
    args = parse_args(unparsed_args)

    with open(args.input_fastq_fp, 'r') as f:
        fq_reader = read_fastq(f)
        if args.collapse_umi:
            input_seqs = []
            input_umis = []
            for header, seq, qual in fq_reader:
                input_seqs.append(seq)
                input_umis.append(extract_umi(header, args.umi_start, args.umi_len))
        else:
            input_seqs = [seq for header, seq, qual in fq_reader]

    with (open(args.clustered_fp, 'r') as clustered_f,
          open(args.output_fp, 'w') as output_f,
          open(args.log_fp, 'w') as log_f):
        
        csv_reader = csv.reader(clustered_f, delimiter='\t')

        for row in tqdm(csv_reader):
            seq_indices = [int(i) for i in row[3].split(',')]
            full_seqs = [input_seqs[i - 1] for i in seq_indices] # starcode indices are 1-based
            consensus, count, not_matching = get_consensus_and_count(full_seqs, max_dist=args.max_dist)

            if not args.collapse_umi:
                output_f.write(f"{consensus}\t{count}\n")               
            else:
                umis = [input_umis[i - 1] for i in seq_indices] # starcode indices are 1-based
                umis_passed = [umi for i, umi in enumerate(umis) if i not in not_matching]
                if count > 0:
                    umis_clustered = cluster_umis(umis_passed, args.umi_threshold)
                    umi_count = len(umis_clustered)
                    output_f.write(f"{consensus}\t{count}\t{umi_count}\n")

            if len(not_matching) or count == 0:
                log_f.write(f"Consensus: {consensus}\n")
                for i in not_matching:
                    log_f.write(f"           {full_seqs[i]}\n")


if __name__ == "__main__":
    main(sys.argv[1:])