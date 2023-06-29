# Takes two fastq files (read 1 and read 2) and combines reads into one sequence with strandedness of read 1.
# Assumes the fastqs have been trimmed of adapter sequences and that read 2 is aligned with read 1
# (i.e. they fully overlap), but there may be some substitution errors.
# Mostly taken from https://github.com/Carldeboer/CisRegModels/blob/master/alignFastqsIntoSeqs.py

import argparse
from utils import read_fastq, revcomp
import multiprocessing 
from functools import partial


def get_alignment_score(seq1, seq2):
    '''Calculate an alignment score between two DNA sequences.
    
       Takes two sequences and calculates an alignment score between 0 and 1,
       which is the sum of matching positions divided by the total length
       of sequence compared. Sequences must be the same length.
    '''
    if len(seq1) != len(seq2):
        raise ValueError("seq1 length not equal to seq2 length")
    # don't count as a match if both bases are N
    num_matches = sum((base1 == base2) and (base1 != 'N' or base2 != 'N') for base1, base2 in zip(seq1, seq2))
    return num_matches/(len(seq1)) 


def get_consensus(read1, read2):
    '''Get consensus sequence from two reads, where reads are a (header, sequence, quality) tuple.
       Uses quality score to chose which base to keep. Returns a (header, sequence, quality) tuple.
       Returned header is from read 1, and quality string is the max quality at each position.
    '''
    consensus = []
    quality = []
    if len(read1[1]) != len(read1[2]) or len(read2[1]) != len(read2[2]):
        raise ValueError(f"Length of sequence does not match length of quality string for read {read1[0]} or {read2[0]}")
    for base1, qual1, base2, qual2 in zip(read1[1], read1[2], read2[1], read2[2]):
        if qual1 >= qual2: #converts to ASCII, corresponds to quality
            consensus.append(base1)
            quality.append(qual1)
        else:
            consensus.append(base2)
            quality.append(qual2)

    return (read1[0], "".join(consensus), "".join(quality))
           

def get_align_result_string(reads, align_threshold, seq_only):
    ''' Returns alignment result string for alignment of two reads.
    Returns a tuple: first element is the string to write if alignment was
    successful. Second element is the string to write to log file if
    alignment was not successful
    '''
    r1 = reads[0]
    r2 = reads[1]
    r1_seq = r1[1]
    r2_seq_revcomp = revcomp(r2[1])
    r2_revcomp = (r2[0], r2_seq_revcomp, r2[2][::-1])

    align_string = None
    log_string = None

    try:
        align_score = get_alignment_score(r1_seq, r2_seq_revcomp)
    except ValueError:
        align_score = -1 # seq lengths were not equal; write to log file
    if align_score > align_threshold:
        consensus = get_consensus(r1, r2_revcomp)
        if seq_only:
            align_string = consensus[1] + "\n"
        else:
            read = "\n".join([consensus[0], consensus[1], "+", consensus[2]])
            align_string = f"{read}\n"
    else:
        log_string = (f"Skipping {r1[0]}/{r2[0]} because they don't align within given "
                        f"parameters\n              {r1[2]}\nr1:           {r1_seq}\nr2 (revcomp): {r2_seq_revcomp}\n              {r2[2]}\n")
    
    return (align_string, log_string)


def pair_reads_and_save(r1_fastq_fp, r2_fastq_fp, out_fp, log_fp, align_threshold, seq_only=True):
    '''Combine paired end reads from read1 and read2 fastq files and save as a single fastq.
       Assumes read1 and read2 are the same length and are reverse complements with the potential
       for some mismatches. Saves read pairs that do not align to a log file'''
    with open(r1_fastq_fp, 'rt') as r1_file, \
         open(r2_fastq_fp, 'rt') as r2_file, \
         open(out_fp, 'wt') as out_file, \
         open(log_fp, 'wt') as log_file:

        r1_reader = read_fastq(r1_file)
        r2_reader = read_fastq(r2_file)

        for reads in zip(r1_reader, r2_reader):
            align_string, log_string = get_align_result_string(reads, align_threshold, seq_only)
            if align_string:
                out_file.write(align_string)
            else:
                log_file.write(log_string)


def pair_reads_and_save_mp(r1_fastq_fp, r2_fastq_fp, out_fp, log_fp, align_threshold, cpus, seq_only=True):
    '''Combine paired end reads from read1 and read2 fastq files and save as a single fastq.
       Uses multiprocessing for the alignment step.
       Assumes read1 and read2 are the same length and are reverse complements with the potential
       for some mismatches. Saves read pairs that do not align to a log file'''
    with open(r1_fastq_fp, 'rt') as r1_file, \
         open(r2_fastq_fp, 'rt') as r2_file, \
         open(out_fp, 'wt') as out_file, \
         open(log_fp, 'wt') as log_file:

        r1_reader = read_fastq(r1_file)
        r2_reader = read_fastq(r2_file)

        p = multiprocessing.Pool(cpus)
        
        for align_string, log_string in p.imap(partial(get_align_result_string, align_threshold=align_threshold, seq_only=seq_only), zip(r1_reader, r2_reader), chunksize=1000):
            if align_string:
                out_file.write(align_string)
            else:
                log_file.write(log_string)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--r1', dest='r1_fastq_fp', help="Input fastq read 1 file", required=True)
    parser.add_argument('--r2', dest='r2_fastq_fp', help="Input fastq read 2 file", required=True)
    parser.add_argument('--out', dest='out_fp', help="Output file", required=True)
    parser.add_argument('--log', dest='log_fp', help="Log file", required=True)
    parser.add_argument('-t', '--cpus', type=int, dest='num_cpus', default=1, help="Number of cpus for multiprocessing")
    parser.add_argument('--threshold', type=float, default=0.8, help="Fraction of aligned read1/2 that match must be greater than this value")
    parser.add_argument('--output-fastq', action='store_true', default=False, help="Output as a fastq file with header and quality (default output is sequence only)")
    args = parser.parse_args()

    save_seq_only = not args.output_fastq
    if args.num_cpus == 1:
        pair_reads_and_save(args.r1_fastq_fp, args.r2_fastq_fp, args.out_fp, args.log_fp, align_threshold=args.threshold, seq_only=save_seq_only)
    else:
        pair_reads_and_save_mp(args.r1_fastq_fp, args.r2_fastq_fp, args.out_fp, args.log_fp, align_threshold=args.threshold, cpus=args.num_cpus, seq_only=save_seq_only)


if __name__ == "__main__":
    main()