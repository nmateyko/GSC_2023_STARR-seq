# This script is for getting full length cluster centers after using
# starcode to cluster partial sequences. Starcode must be run with the --seq-id
# option to generate an output file with truncated cluster centers and a list
# of sequence positions in the input file for each read in the cluster. This script
# takes the sequence positions, gets the full-length sequences from the starcode input, and
# calculates a consensus full-length sequence after ensuring all full-length sequences in the
# cluster are within some threshold of similarity. It then outputs a list of full-length
# sequences and their read counts as a tsv.

import argparse
import csv
from collections import Counter
from itertools import zip_longest
from utils import read_fastq
import Levenshtein

def most_common(l):
    '''Returns most common element in a list.
    In the case of a tie, the element that occurs first in the list is returned.
    '''
    return Counter(l).most_common(1)[0][0]


def get_consensus(seqs):
    '''Get consensus sequence by calculating most common base at each position.
       In the case of a tie, the base that occurred first in the list is chosen.
       Sequences of different lengths can be present; the most frequent length is
       used for the final consensus sequence. In the case of a tie, the length of the
       first sequence in the list with one of the most frequent lengths is used.
    '''

    # Find most common sequence length.
    seq_len = most_common([len(seq) for seq in seqs])

    consensus = []
    for bases in zip_longest(*seqs, fillvalue ='N'):
        # If all bases are the same, don't try to find the most common base
        # (assuming most positions have bases that all match, this gives a massive speedup)
        if len(set(bases)) == 1:
            consensus.append(bases[0])
        # If not all bases match, then find the most common base
        else:
            consensus.append(most_common(bases))

    return "".join(consensus)[:seq_len]


def get_consensus_and_count(seqs, max_dist):
    '''Get full length consensus sequence and count from list of sequences.
       Any sequences further than Levenshtein distance of max_dist from the consensus
       are not counted and are returned as a list.

       Parameters:
           seqs (list of str): list of sequences to count
           max_dist (int): Maximum distance between sequence and the cluster consensus that is allowed;
                           sequences with larger distances will be discarded.

        Returns:
            Tuple with 3 elements: the consensus sequence, the count, and a list of sequences that
            did not match the consensus within the specified max_dist.
    '''
    consensus = get_consensus(seqs)
    not_matching = []
    count = 0
    for seq in seqs:
        if Levenshtein.distance(consensus, seq) > max_dist:
            not_matching.append(seq)
        else:
            count += 1
    return (consensus, count, not_matching)


def full_seq_generator(starcode_reader, input_seqs):
    '''Generator that takes in a starcode output file reader and list of full length input seqs
       and returns the full length sequences corresponding to the indexes in the starcode file.
    '''
    for row in starcode_reader:
        seq_indices = [int(i) for i in row[3].split(',')]
        full_seqs = [input_seqs[i - 1] for i in seq_indices]
        yield full_seqs


def get_full_seqs_from_starcode_clusters(input_fastq_fp, clustered_fp, output_fp, log_fp, max_dist):
    '''Converts the output of starcode (run with --seq-id) run on truncated
       sequences to a tab-separated file of full length cluster centers and counts.

       Parameters:

           input_fastq_fp (str): Path to the exact fastq file that was used
                                 to generate the starcode clustered output.
           clustered_fp (str): Path to starcode output file; expected starcode output
                               format is <seq>\t<count>\t<seq>\t<index1,index2,...,indexn>\n.
           output_fp (str): Path for output file.
           log_fp (str): Path for log file.
           max_dist (int): Maximum distance between sequence and the cluster consensus that is allowed;
                           sequences with larger distances will be discarded.   
    '''
   
    # Read in fastq and save to list. May need to index the fastq for larger files in the future
    # so that the whole fastq doesn't have to be held in memory (check out pyfastx).
    with open(input_fastq_fp, 'r') as f:
        fq_reader = read_fastq(f)
        input_seqs = [seq for header, seq, qual in fq_reader]

    with (open(clustered_fp, 'r') as clustered_f,
          open(output_fp, 'w') as output_f,
          open(log_fp, 'w') as log_f):
        
        csv_reader = csv.reader(clustered_f, delimiter='\t')

        for row in csv_reader:
            seq_indices = [int(i) for i in row[3].split(',')]
            full_seqs = [input_seqs[i - 1] for i in seq_indices]
            consensus, count, not_matching = get_consensus_and_count(full_seqs, max_dist=max_dist)
            output_f.write(f"{consensus}\t{count}\n")
            if len(not_matching):
                log_f.write(f"Consensus: {consensus}\n")
                for seq in not_matching:
                    log_f.write(f"           {seq}\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', dest='input_fastq_fp', help="Fastq file used as starcode input", required=True)
    parser.add_argument('-c', dest='clustered_fp', help="Starcode output file (must use --seq-id option with starcode)", required=True)
    parser.add_argument('-o', dest='output_fp', help="Output file", required=True)
    parser.add_argument('-l', dest='log_fp', help="Log file", required=True)
    parser.add_argument('-d', dest='max_dist', type=int, help="Maximum distance between sequence and cluster consensus", required=True)
    args = parser.parse_args()

    get_full_seqs_from_starcode_clusters(args.input_fastq_fp, args.clustered_fp, args.output_fp, args.log_fp, args.max_dist)


if __name__ == "__main__":
    main()

# get_full_seqs_from_starcode_clusters("test_files/Sahu_DNA_rep1_sample.fastq", "test_files/Sahu_DNA_rep1_sample_clustered.txt", "test_out.txt", "test_log.txt", 20)