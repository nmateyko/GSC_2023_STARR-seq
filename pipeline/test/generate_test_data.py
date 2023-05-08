# This script generates two test fastq.gz files (read 1 and read 2) to run through the Snakemake pipeline

import gzip
import random
from datetime import datetime


LENGTHS = (83, 100, 132) # lengths of random sequences in library

def mutate_seq(seq, p):
    '''Mutates DNA sequence (string) with a probability of p at each position and
    returns the mutated sequence'''
    mutated_seq = []
    for base in seq:
        if random.random() <= p:
            mutated_seq.append(random.choice(tuple({'A', 'C', 'G', 'T'} - {base})))
        else:
            mutated_seq.append(base)
    return "".join(mutated_seq)

def rev_comp(seq):
    comp = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    return ''.join([comp[base] for base in seq[::-1]])

def generate_fastq(seqs, filename, read):
    with gzip.open(filename, 'wt') as f:
        for i, seq in enumerate(seqs):
            f.write(f"@NOVASEQ1:432:HG25JDSX5:1:1101:{i}:1000 {read}:N:0:NNNNNNNN+NNNNNNNN\n")
            f.write(f"{seq}\n")
            f.write(f"+\n")
            f.write(f"{'F' * len(seq)}\n")

# make random sequences of different lengths
seqs = []
for length in LENGTHS:
    for i in range(5):
        seqs.append("".join(random.choices(["A", "C", "T", "G"], k=length)))

# "amplify" each seq to get different read counts
counts = [i for i in range(2, len(seqs) + 2)] # number of reads for each sequence in seqs
seqs_amplified = []
for seq, count in zip(seqs, counts):
    seqs_amplified += [seq for i in range(count)]

# mutate seqs
mutated_seqs = [mutate_seq(seq, 0.05) for seq in seqs_amplified]

# add original seqs so that true sequences are present
full_seqs = mutated_seqs + seqs

# save seqs and expected counts; +1 because original seqs added back to mutated seqs
seq_counts = "\n".join([f"{seq}\t{count + 1}" for seq, count in zip(seqs, counts)])
date_string =  datetime.now().strftime("%Y_%m_%d_%H_%M")
with open(f"test_seq_true_counts_{date_string}.txt", 'w') as f:
    f.write(seq_counts)

# shuffle seqs
random.shuffle(full_seqs)

seqs_revcomp = [rev_comp(seq) for seq in full_seqs]

# add truseq adapter downstream of r1 and r2 seqs
r1_downstream = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG"
r2_downstream = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT"

r1_reads = [seq + r1_downstream for seq in full_seqs]
r2_reads = [seq + r2_downstream for seq in seqs_revcomp]

generate_fastq(r1_reads, f"test_r1_{date_string}.fastq.gz", "1")
generate_fastq(r2_reads, f"test_r2_{date_string}.fastq.gz", "2")