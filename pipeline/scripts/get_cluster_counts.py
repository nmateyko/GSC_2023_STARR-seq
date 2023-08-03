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
    '''Get most common element in a list. In the case of a tie, returns the
       most common element that occurs first in the list.
    '''
    return Counter(l).most_common(1)[0][0]


def get_consensus(seq_counts):
    '''Get consensus sequence by calculating most common base at each position.
       In the case of a tie, the base that occurred first in the list is chosen.
       Sequences of different lengths can be present; the most frequent length is
       used for the final consensus sequence. In the case of a tie, the length of the
       first sequence in the list with one of the most frequent lengths is used.

       Arguments:
            seq_counts (dict str:int): a dictionary with DNA
            sequences as keys and the number of occurrences of the sequence as values.
        Returns:
            A string which is the consensus of the keys of the input dictionary,
            weighted by their counts.
    '''

    seq_len_dict = {}
    for seq, count in seq_counts.items():
        seq_len = len(seq)
        if seq_len in seq_len_dict:
            seq_len_dict[seq_len] += count
        else:
            seq_len_dict[seq_len] = count

    seq_len = most_common(seq_len_dict)

    consensus = []

    counts = [i for i in seq_counts.values()]

    for bases in zip_longest(*seq_counts, fillvalue ='N'):
        base_count_dict = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}
        for base, count in zip(bases, counts):
            base_count_dict[base] += count
        consensus.append(most_common(base_count_dict))

    return "".join(consensus)[:seq_len]


def get_consensus_and_count(seq_counts, max_dist):
    '''Get full length consensus sequence and count from dictionary of sequences and their counts.
       Any sequences further than Levenshtein distance of max_dist from the consensus
       are not counted and are returned as a list.

       Arguments:
           seq_counts (dict str:int): Dictionary with sequences as keys and their counts as values.
           max_dist (int): Maximum distance between sequence and the consensus that is allowed;
                           sequences with larger distances will be discarded.

        Returns:
            Tuple with 3 elements: the consensus sequence, the count, and a set of sequences that
            did not match the consensus within the specified max_dist.
    '''
    consensus = get_consensus(seq_counts)
    not_matching = []
    total_count = 0
    for seq, count in seq_counts.items():
        if Levenshtein.distance(consensus, seq) > max_dist:
            not_matching.append(seq)
        else:
            total_count += count
    return (consensus, total_count, set(not_matching))


def extract_umi(header, start_distance, umi_len):
    '''Extract a UMI sequence from a fastq header string.

       Arguments:
           header (str): Fastq read header containing the UMI.
           start distance (int): Number of positions from the end of the header to the start of the UMI;
                                 equal to len(header) - (index of first UMI character).
           umi_len (int): Length of the UMI.

        Returns:
            The umi sequence (str).
    '''
    umi = header[-start_distance:-start_distance + umi_len]
    if not set(umi).issubset({'A', 'C', 'G', 'T', 'N'}):
        raise ValueError(f"UMI contains non-ACTGN character: {umi}")
    return umi


def cluster_umis(umi_count_dict, threshold, clust_method="directional"):
    '''Cluster a dict of UMIs and counts using UMI-tools.

    Arguments:

        umis (dict str:int): Dict of UMI sequences and their counts to be clustered.
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
    '''Add a sequence to a dictionary of sequence clusters. The dict keys are cluster IDs,
    and the values are dicts that contain sequences as keys and counts as values. When the
    sequence is added to the cluster dict, its count value is incremented by 1.
    '''
    if seq in seq_dict[cluster_index]:
        seq_dict[cluster_index][seq] += 1
    else:
        seq_dict[cluster_index][seq] = 1


def add_seq_and_umi_to_dict(seq, umi, cluster_index, seq_dict):
    '''Add a sequence and its UMI to a dictionary of sequence clusters. The dict keys are cluster IDs,
    and the values are dicts that contain sequences as keys and dicts of UMIs and their counts as values.
    When the sequence is added to the cluster dict, its UMI count value is incremented by 1.
    '''
    if seq in seq_dict[cluster_index]:
        if umi in seq_dict[cluster_index][seq]:
            seq_dict[cluster_index][seq][umi] += 1
        else:
            seq_dict[cluster_index][seq][umi] = 1
    else:
        seq_dict[cluster_index][seq] = {umi: 1}


def get_count_string(seq_counts, max_dist):
    ''' Takes a dict of sequence keys and count values and calculates a consensus
    sequence and total count. Sequences further than max_dist Levenshtein distance
    from the consensus are discarded. Returns a count_string (sequence and count)
    and a log_string (discarded sequences and their counts) as a tuple.
    '''
    count_string = None
    log_string = None

    consensus, count, not_matching = get_consensus_and_count(seq_counts, max_dist)
    count_string = f"{consensus}\t{count}\n"

    if len(not_matching) or count == 0:
        log_string = [f"Consensus: {consensus}\n"]
        for seq in not_matching:
            log_string.append(f"           {seq}\t{seq_counts[seq]}\n")
        log_string = "".join(log_string)

    return (count_string, log_string)


def get_count_string_umi(seq_umi_counts, max_dist, umi_threshold):
    ''' Takes a dict of sequence keys and umi count dict values and calculates a consensus
    sequence and total count before and after umi collapsing. Sequences further than max_dist
    Levenshtein distance from the consensus are discarded, and umis are clustered
    using UMI-tools with umi_threshold. Returns a count_string (sequence, count, and collapsed count)
    and a log_string (discarded sequences and their counts) as a tuple.
    '''
    count_string = None
    log_string = None
    seq_counts = {}
    for seq, umi_counts in seq_umi_counts.items():
        seq_counts[seq] = sum(i for i in umi_counts.values())

    consensus, count, not_matching = get_consensus_and_count(seq_counts, max_dist)

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

    # Read in the starcode output and map each read index to
    # its cluster ID (line number in the starcode output file).
    with open(args.clustered_fp, 'r') as f:
        csv_reader = csv.reader(f, delimiter='\t')
        index_cluster_mapping = {}
        n_clusters = 0
        for i, row in enumerate(csv_reader):
            indices = [int(i) for i in row[3].split(',')]
            n_clusters += 1
            for j in indices:
                index_cluster_mapping[j - 1] = i

    # Read in the original (unclustered) fastq and add each
    # sequence to its corresponding cluster dictionary, or increment
    # the sequence count in the dictionary if it's already present.
    with open(args.input_fastq_fp, 'r') as f:
        fq_reader = read_fastq(f)
        seq_clust_counts = {i: {} for i in range(n_clusters)}
        for i, (header, seq, qual) in enumerate(fq_reader):
            cluster_index = index_cluster_mapping[i]
            if args.collapse_umi:
                umi = extract_umi(header, args.umi_start, args.umi_len)
                add_seq_and_umi_to_dict(seq, umi, cluster_index, seq_clust_counts)
            else:
                add_seq_to_dict(seq, cluster_index, seq_clust_counts)

    with (open(args.output_fp, 'w') as output_f,
          open(args.log_fp, 'w') as log_f):
        
        p = multiprocessing.Pool(args.num_cpus)
        
        if args.collapse_umi:
            # shuffle seq_clust_counts dict so that the largest clusters aren't all sent to the same core
            cluster_indices = list(seq_clust_counts.keys())
            random.shuffle(cluster_indices)
            seq_clust_counts_shuffled = [seq_clust_counts[i] for i in cluster_indices]

            # For each cluster in seq_clust_counts, get the consensus, count,
            # and collapsed count from the sequence count dictionary for that cluster.
            for count_string, log_string in p.imap(partial(get_count_string_umi, max_dist=args.max_dist, umi_threshold=args.umi_threshold), seq_clust_counts_shuffled, chunksize=1000):
                if count_string:
                    output_f.write(count_string)
                if log_string:
                    log_f.write(log_string)

        else:
            # For each cluster in seq_clust_counts, get the consensus and count from the sequence count dictionary for that cluster.
            for count_string, log_string in p.imap(partial(get_count_string, max_dist=args.max_dist), seq_clust_counts.values(), chunksize=1000):
                if count_string:
                    output_f.write(count_string)
                if log_string:
                    log_f.write(log_string)


if __name__ == "__main__":
    main(sys.argv[1:])