import os
import pytest
import random
import re
from pair_reads import read_fastq, revcomp, get_alignment_score, get_consensus, pair_reads_and_save

def test_read_fastq():
    with open("test_files/test.fastq", 'r') as f:
        fastq_reader = read_fastq(f)
        reads = [read for read in fastq_reader]
    expected = [("@NOVASEQ1:432:HG25JDSX5:1:1101:2844:1000 2:N:0:ACTATTGC+TCCGCCTA", "GCACAAACCGTT", "AFFFFFFFFFFF"),
                ("@NOVASEQ1:432:HG25JDSX5:1:1101:3314:1000 2:N:0:ACTATTGC+TCCGCCTA", "TGATAAGTCACA", "FFFFFFFFFFFF"),
                ("@NOVASEQ1:432:HG25JDSX5:1:1101:4255:1000 2:N:0:ACTATTGC+TCCGCCTA", "AACACACCACAT", "FFFFFFFFFFF,"),
                ("@NOVASEQ1:432:HG25JDSX5:1:1101:4634:1000 2:N:0:ACTATTGC+TCCGCCTA", "CCCTGGAGGCAA", "FFFFFFFFFFFF"),
                ("@NOVASEQ1:432:HG25JDSX5:1:1101:7021:1000 2:N:0:ACTATTGC+TCCGCCTA", "ACGCTACATCTT", "FFFFFFFFFFFF")
    ]
    assert reads == expected

def test_read_fastq_EOFError():
    # file has last line deleted
    with open("test_files/test_EOFError.fastq", 'r') as f:
        with pytest.raises(EOFError, match="Failed to parse four lines from fastq file!"):
            reads = [read for read in read_fastq(f)]

def test_read_fastq_ValueError():
    # file has the '+' on one line replaced with 'A'
    with open("test_files/test_ValueError.fastq", 'r') as f:
        with pytest.raises(ValueError, match=re.escape("Invalid header lines: @NOVASEQ1:432:HG25JDSX5:1:1101:4634:1000 2:N:0:ACTATTGC+TCCGCCTA and A for seq CCCTGGAGGCAA")):
            reads = [read for read in read_fastq(f)]

@pytest.mark.parametrize("test_input, expected", [("AAAAA", "TTTTT"), ("NCTGACTG", "CAGTCAGN"), ("AAATTT", "AAATTT")])
def test_revcomp(test_input, expected):
    assert revcomp(test_input) == expected

def test_revcomp_ValueError():
    with pytest.raises(ValueError, match=r"Sequence \(ACTGH\) must only contain ACTGN"):
        revcomp("ACTGH")

@pytest.mark.parametrize("seq1, seq2, expected", [("AAAAA", "AAAAA", 1), ("ACTGACTG", "ACTGACTA", 0.875), ("ANNNN", "ANNNN", 0.2), ("GGGGG", "CCCCC", 0)])
def test_get_alignment_score(seq1, seq2, expected):
    assert get_alignment_score(seq1, seq2) == expected

def test_get_alignment_score_ValueError():
    with pytest.raises(ValueError, match=re.escape("seq1 length not equal to seq2 length")):
        get_alignment_score("AAAAA", "AAAAAA")

def test_get_consensus_same_seq():
    r1 = ("@header", "AAAAA", "FFFFF")
    assert get_consensus(r1, r1) == r1

def test_get_consensus_first_seq_better():
    r1 = ("@header", "AAAAA", "FFFFF")
    r2 = ("@header", "TAAAA", ",FFFF")
    assert get_consensus(r1, r2) == ("@header", "AAAAA", "FFFFF")

def test_get_consensus_second_seq_better():
    r1 = ("@header", "AAAAA", ",FFFF")
    r2 = ("@header", "TAAAA", "FFFFF")
    assert get_consensus(r1, r2) == ("@header", "TAAAA", "FFFFF")

def test_get_consensus_mixed():
    r1 = ("@header", "ATAAG", ",DFFE")
    r2 = ("@header", "TAAAA", "FAFF,")
    assert get_consensus(r1, r2) == ("@header", "TTAAG", "FDFFE")

def get_consensus_ValueError():
    r1 = ("@header", "ATAAGA", ",DFFE")
    r2 = ("@header", "TAAAA", "FAFF,")
    with pytest.raises(ValueError, match=re.escape("Length of sequence does not match length of quality string for read @header or @header")):
        get_consensus(r1, r2)

def test_pair_reads_and_save():
    id = random.randint(0, 100000000)
    outfile = f"test_files/test_pair_out_{id}.fastq"
    logfile = f"test_files/test_pair_out_{id}.log"
    pair_reads_and_save("test_files/pair_test_r1.fastq", "test_files/pair_test_r2.fastq", outfile, logfile, align_threshold=0.8)
    with open(outfile, 'r') as f:
        pair_reads_output = f.readlines()
    with open("test_files/paired_expected.fastq", 'r') as f:
        pair_reads_expected = f.readlines()
    os.remove(outfile)
    os.remove(logfile)
    assert pair_reads_output == pair_reads_expected

def test_pair_reads_and_save_log():
    id = random.randint(0, 100000000)
    outfile = f"test_files/test_pair_out_{id}.fastq"
    logfile = f"test_files/test_pair_out_{id}.log"
    pair_reads_and_save("test_files/pair_test_r1.fastq", "test_files/pair_test_r2.fastq", outfile, logfile, align_threshold=0.8)
    with open(logfile, 'r') as f:
        pair_reads_log = f.readlines()
    with open("test_files/paired_expected.log", 'r') as f:
        pair_reads_expected = f.readlines()
    os.remove(outfile)
    os.remove(logfile)
    assert pair_reads_log == pair_reads_expected