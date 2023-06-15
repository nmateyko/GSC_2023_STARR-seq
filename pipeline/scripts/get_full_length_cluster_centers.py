# This script is for getting full length cluster centers after using
# starcode to cluster partial sequences. Starcode must be run with the --seq-id
# option to generate an output file with truncated cluster centers and a list
# of sequence positions in the input file for each read in the cluster. This script
# takes the sequence positions, gets the full-length sequences from the starcode input, and
# calculates a consensus full-length sequence after ensuring all full-length sequences in the
# cluster are within some threshold of similarity. It then outputs a list of full-length
# sequences and their read counts as a tsv.

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
       In the case of a tie, the base that occured first in the list is chosen.
       Sequences of different lengths can be present; the most frequent length is
       used for the final consensus sequence. In the case of a tie, the length of the
       first sequence in the list with one of the most frequent lengths is used.
    '''

    # Find most common sequence length.
    seq_len = most_common([len(seq) for seq in seqs])

    consensus = "".join([most_common(bases) for bases in zip_longest(*seqs, fillvalue ='N')])
    return consensus[:seq_len]

def get_full_seqs_from_starcode_clusters(input_fastq_fp, clustered_fp, output_fp, log_fp, dist_fun, max_dist):
    '''Converts the output of starcode (run with --seq-id) run on truncated
       sequences to a tab-separated file of full length cluster centers and counts.

       Parameters:

           input_fastq_fp (str): Path to the exact fastq file that was used
                                 to generate the starcode clustered output.
           clustered_fp (str): Path to starcode output file; expected starcode output
                               format is <seq>\t<count>\t<seq>\t<index1,index2,...,indexn>\n.
           output_fp (str): Path for output file.
           log_fp (str): Path for log file.
           dist_fun (returns int): Distance function for comparing similarity of full length sequences.
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
        
        csv_reader = csv.reader(f, delimiter='\t')
        for row in csv_reader:
            # Get full length sequences of cluster components
            seq_indices = [int(i) for i in row[3].split(',')]
            full_seqs = [input_seqs[i] for i in seq_indices]
            # Compare each full sequence to the consensus
            consensus = get_consensus(full_seqs)
            not_matching = []
            count = 0
            for seq in full_seqs:
                if dist_fun(consensus, seq) > max_dist:
                    not_matching.append(seq)
                else:
                    count += 1
            output_f.write(f"{consensus}\t{count}")
            if len(not_matching):
                log_f.write(f"Consensus: {consensus}")
                for seq in not_matching:
                    log_f.write(f"           {seq}")
            


get_full_seqs_from_starcode_clusters("/scratch/st-cdeboer-1/najmeh/GSC_2023_STARR-seq/pipeline/output/paired/Sahu_DNA_rep1_assembled.fastq", "/scratch/st-cdeboer-1/najmeh/GSC_2023_STARR-seq/pipeline/output/clustered/Sahu_DNA_rep1_clustered.txt", "test_out.txt", "test_log.txt", Levenshtein.distance, 10)