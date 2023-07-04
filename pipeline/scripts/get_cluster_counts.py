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
import multiprocessing 
import os
import random
import sys
from collections import Counter
from functools import partial
from itertools import zip_longest
from umi_tools import UMIClusterer
from utils import read_fastq
import Levenshtein

csv.field_size_limit(sys.maxsize)

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


def get_consensus_counter(seqs_count_dict):

    seq_len_dict = {}
    for seq, count in seqs_count_dict.items():
        seq_len = len(seq)
        if seq_len in seq_len_dict:
            seq_len_dict[seq_len] += count
        else:
            seq_len_dict[seq_len] = count

    seq_len = most_common(seq_len_dict)

    consensus = []

    seq_counts = [i for i in seqs_count_dict.values()]

    for bases in zip_longest(*seqs_count_dict, fillvalue ='N'):
        base_count_dict = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}
        for base, count in zip(bases, seq_counts):
            base_count_dict[base] += count
        consensus.append(most_common(base_count_dict))

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


def get_consensus_and_count_counter(seqs_count_dict, max_dist):
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
    consensus = get_consensus_counter(seqs_count_dict)
    not_matching = []
    total_count = 0
    for seq, count in seqs_count_dict.items():
        if Levenshtein.distance(consensus, seq) > max_dist:
            not_matching.append(seq)
        else:
            total_count += count
    return (consensus, total_count, set(not_matching))


def extract_umi(header, start_distance, umi_len):
    umi = header[-start_distance:-start_distance + umi_len]
    if not set(umi).issubset({'A', 'C', 'G', 'T', 'N'}):
        raise ValueError(f"UMI contains non-ACTGN character: {umi}")
    return umi


def cluster_umis(umi_count_dict, threshold, clust_method="directional"):
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
    
    umi_count_dict = {umi.encode(): count for umi, count in umi_count_dict.items()}
    clusterer = UMIClusterer(cluster_method=clust_method)
    clusters = clusterer(umi_count_dict, threshold=threshold)

    return clusters


def add_seq_to_dict(seq, cluster_index, seq_dict):
    if seq in seq_dict[cluster_index]:
        seq_dict[cluster_index][seq] += 1
    else:
        seq_dict[cluster_index][seq] = 1


def add_seq_and_umi_to_dict(seq, umi, cluster_index, seq_dict):
    if seq in seq_dict[cluster_index]:
        if umi in seq_dict[cluster_index][seq]:
            seq_dict[cluster_index][seq][umi] += 1
        else:
            seq_dict[cluster_index][seq][umi] = 1
    else:
        seq_dict[cluster_index][seq] = {umi: 1}


def get_count_string(seq_counts, max_dist):
    count_string = None
    log_string = None

    consensus, count, not_matching = get_consensus_and_count_counter(seq_counts, max_dist)
    count_string = f"{consensus}\t{count}\n"

    if len(not_matching) or count == 0:
        log_string = [f"Consensus: {consensus}\n"]
        for seq in not_matching:
            log_string.append(f"           {seq}\t{seq_counts[seq]}\n")
        log_string = "".join(log_string)

    return (count_string, log_string)


def get_count_string_umi(seq_umi_counts, max_dist, umi_threshold):
    count_string = None
    log_string = None
    seq_counts = {}
    for seq, umi_counts in seq_umi_counts.items():
        seq_counts[seq] = sum(i for i in umi_counts.values())

    consensus, count, not_matching = get_consensus_and_count_counter(seq_counts, max_dist)

    umis_passed = {}
    for seq, umi_counts_dict in seq_umi_counts.items():
        if seq not in not_matching:
            for umi, umi_count in umi_counts_dict.items():
                if umi in umis_passed:
                    umis_passed[umi] += umi_count
                else:
                    umis_passed[umi] = umi_count
    if count > 0:
        umis_clustered = cluster_umis(umis_passed, umi_threshold)
        umi_count = len(umis_clustered)
        count_string = f"{consensus}\t{count}\t{umi_count}\n"

    if len(not_matching) or count == 0:
        log_string = [f"Consensus: {consensus}\n"]
        for seq in not_matching:
            log_string.append(f"           {seq}\t{seq_counts[seq]}\n")
        log_string = "".join(log_string)

    return (count_string, log_string)
    

def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fastq', dest='input_fastq_fp', help="Fastq file used as starcode input", required=True)
    parser.add_argument('-c', '--clustered', dest='clustered_fp', help="Starcode output file (must use --seq-id option with starcode)", required=True)
    parser.add_argument('-o', '--output', dest='output_fp', help="Output file path", required=True)
    parser.add_argument('-l', '--log', dest='log_fp', help="Log file path", required=True)
    parser.add_argument('-d', '--max_dist', dest='max_dist', type=int, default=10, help="Maximum distance between sequence and cluster consensus")
    parser.add_argument('-t', '--cpus', type=int, dest='num_cpus', default=1, help="Number of cpus for multiprocessing")
    parser.add_argument('--umi', dest='collapse_umi', action='store_true', default=False, help="Collapse sequences using UMIs")
    parser.add_argument('--umi-threshold', default=1, type=int, help="Max Levenshtein distance for collapsing UMIs")
    parser.add_argument('--umi-start', type=int, help="Distance of UMI start index from the end of the header line (len(header) - index of UMI start)", required='--umi' in sys.argv)
    parser.add_argument('--umi-len', type=int, help="Length of UMI", required='--umi' in sys.argv)
    return parser.parse_args(args)


def main(unparsed_args):
    args = parse_args(unparsed_args)

    with open(args.clustered_fp, 'r') as f:
        csv_reader = csv.reader(f, delimiter='\t')
        index_cluster_mapping = {}
        cluster_sizes = []
        for i, row in enumerate(csv_reader):
            indices = [int(i) for i in row[3].split(',')]
            cluster_sizes.append(len(indices))
            for j in indices:
                index_cluster_mapping[j - 1] = i

    with open(args.input_fastq_fp, 'r') as f:
        fq_reader = read_fastq(f)
        input_seqs = {i: {} for i in range(len(cluster_sizes))}
        for i, (header, seq, qual) in enumerate(fq_reader):
            cluster_index = index_cluster_mapping[i]
            if args.collapse_umi:
                umi = extract_umi(header, args.umi_start, args.umi_len)
                add_seq_and_umi_to_dict(seq, umi, cluster_index, input_seqs)
            else:
                add_seq_to_dict(seq, cluster_index, input_seqs)

    with (open(args.output_fp, 'w') as output_f,
          open(args.log_fp, 'w') as log_f):
        
        p = multiprocessing.Pool(args.num_cpus)
        
        if args.collapse_umi:
            # shuffle input_seqs dict so that the largest clusters aren't all sent to the same core
            cluster_indices = list(input_seqs.keys())
            random.shuffle(cluster_indices)
            input_seqs_shuffled = [input_seqs[i] for i in cluster_indices]
            for count_string, log_string in p.imap(partial(get_count_string_umi, max_dist=args.max_dist, umi_threshold=args.umi_threshold), input_seqs_shuffled, chunksize=1000):
                if count_string:
                    output_f.write(count_string)
                if log_string:
                    log_f.write(log_string)

        else:
            for count_string, log_string in p.imap(partial(get_count_string, max_dist=args.max_dist), input_seqs.values(), chunksize=1000):
                if count_string:
                    output_f.write(count_string)
                if log_string:
                    log_f.write(log_string)


if __name__ == "__main__":
    main(sys.argv[1:])