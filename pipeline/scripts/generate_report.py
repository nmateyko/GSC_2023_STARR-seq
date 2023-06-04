# This script will generate a report based on the log files from each
# step of the Snakemake pipeline to make it easier to understand how
# many reads are lost at each step of the pipeline and why.

import gzip
from utils import read_fastq


def get_paired_fastq_lengths(r1_fastq_fp, r2_fastq_fp, header_check_length, min=True):
    '''Get a list of either min or max read length in each read pair in r1/r2
       paired-end fastq files. Order of read IDs must be the same in both files.
    '''
    lengths = []
    with gzip.open(r1_fastq_fp, 'rt') as f_r1, gzip.open(r2_fastq_fp, 'rt') as f_r2:
        r1_reader = read_fastq(r1)
        r2_reader = read_fastq(r2)
        for r1, r2 in zip(r1_reader, r2_reader):
            if r1[0][:header_check_length] != r2[0][:header_check_length]:
                raise ValueError(f"Header IDs do not match: {r1[0]} and {r2[0]}")
            if min:
                lengths.append(min(len(r1[1], r2[1])))
            else:
                lengths.append(max(len(r1[1], r2[1])))
    return lengths


def generate_cutadapt_report(log, trimmed, short_1, short_2, long_1, long_2):
    