# Takes two fastq files (read 1 and read 2) and combines reads into one sequence with strandedness of read 1.
# Assumes the fastqs have been trimmed of adapter sequences and that read 2 is aligned with read 1
# (i.e. they fully overlap), but there may be some substitution errors.
# Mostly taken from https://github.com/Carldeboer/CisRegModels/blob/master/alignFastqsIntoSeqs.py

import argparse
from itertools import islice

# stolen from https://gist.github.com/jakebiesinger/759018/1b7d6bd6967780a8bbae743760c37885bdf86467
def read_fastq(fastqfile):
    '''Parse a fastq-formatted file, yielding a (header, sequence, quality) tuple'''
    fastqiter = (l.strip('\n') for l in fastqfile)  # strip trailing newlines 
    fastqiter = filter(lambda l: l, fastqiter)  # skip blank lines
    while True:
        fqlines = list(islice(fastqiter, 4))
        if len(fqlines) == 4:
            header1, seq, header2, qual = fqlines
        elif len(fqlines) == 0:
            return
        else:
            raise EOFError("Failed to parse four lines from fastq file!")

        if header1.startswith('@') and header2.startswith('+'):
            yield header1, seq, qual
        else:
            raise ValueError("Invalid header lines: %s and %s for seq %s" % (header1, header2, seq))


def revcomp(seq):
    '''Return the reverse complement of a DNA sequence'''
    if not set(seq).issubset({'A', 'C', 'G', 'T', 'N'}):
        raise ValueError(f"Sequence ({seq}) must only contain ACTGN")
    comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join([comp[i] for i in seq[::-1]])


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
        if qual1 > qual2: #converts to ASCII, corresponds to quality; if equal, doesn't matter which is better since base is the same
            consensus.append(base1)
            quality.append(qual1)
        else:
            consensus.append(base2)
            quality.append(qual2)

    return (read1[0], "".join(consensus), "".join(quality))

def pair_reads_and_save(r1_fastq_fp, r2_fastq_fp, out_fp, log_fp, align_threshold):
    '''Combine paired end reads from read1 and read2 fastq files and save as a single fastq.
    Assumes read1 and read2 are the same length and are reverse complements with the potential
    for some mismatches. Saves read pairs that do not align to a log file'''
    with open(r1_fastq_fp, 'rt') as r1, \
         open(r2_fastq_fp, 'rt') as r2, \
         open(out_fp, 'wt') as out_file, \
         open(log_fp, 'wt') as log_file:

        r1_reader = read_fastq(r1)
        r2_reader = read_fastq(r2)
        
        for r1_read, r2_read in zip(r1_reader, r2_reader):
            r1_seq = r1_read[1]
            r2_seq_revcomp = revcomp(r2_read[1])

            r2_read_revcomp = (r2_read[0], r2_seq_revcomp, r2_read[2][::-1])
            if get_alignment_score(r1_seq, r2_seq_revcomp) > align_threshold:
                consensus = get_consensus(r1_read, r2_read_revcomp)
                read = "\n".join([consensus[0], consensus[1], "+", consensus[2]])
                out_file.write(f"{read}\n")
            else:
                log_file.write(f"Skipping {r1_read[0]}/{r2_read[0]} because they don't align within given parameters\nr1:           {r1_seq}\nr2 (revcomp): {r2_seq_revcomp}\n")



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--r1', dest='r1_fastq_fp', help="Input fastq read 1 file", required=True)
    parser.add_argument('--r2', dest='r2_fastq_fp', help="Input fastq read 2 file", required=True)
    parser.add_argument('--out', dest='out_fp', help="Output file", required=True)
    parser.add_argument('--log', dest='log_fp', help="Log file", required=True)
    args = parser.parse_args()

    pair_reads_and_save(args.r1_fastq_fp, args.r2_fastq_fp, args.out_fp, args.log_fp, align_threshold=0.8)

    

if __name__ == "__main__":
    main()