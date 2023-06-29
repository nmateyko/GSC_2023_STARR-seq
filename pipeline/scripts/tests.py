import os
import pytest
import datetime
import re
from utils import read_fastq, revcomp
from pair_reads import get_alignment_score, get_consensus, pair_reads_and_save, pair_reads_and_save_mp
from get_cluster_counts import most_common, get_consensus_and_count, extract_umi, cluster_umis
from get_cluster_counts import get_consensus as get_consensus_cluster
from get_cluster_counts import main as get_cluster_counts_main

############################## Test utils.py ##############################

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
    with pytest.raises(ValueError, match=re.escape("Sequence (ACTGH) must only contain ACTGN")):
        revcomp("ACTGH")


############################## Test pair_reads.py ##############################

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

def test_get_consensus_ValueError():
    r1 = ("@header", "ATAAGA", ",DFFE")
    r2 = ("@header", "TAAAA", "FAFF,")
    with pytest.raises(ValueError, match=re.escape("Length of sequence does not match length of quality string for read @header or @header")):
        get_consensus(r1, r2)

def test_pair_reads_and_save():
    try:
        file_id = datetime.datetime.now().strftime("%Y%m%d%H%M%S%f")
        outfile = f"test_files/test_pair_out_{file_id}.fastq"
        logfile = f"test_files/test_pair_out_{file_id}.log"
        pair_reads_and_save("test_files/pair_test_r1.fastq", "test_files/pair_test_r2.fastq", outfile, logfile, align_threshold=0.8, seq_only=False)
        with open(outfile, 'r') as f:
            pair_reads_output = f.readlines()
        with open("test_files/paired_expected.fastq", 'r') as f:
            pair_reads_expected = f.readlines()
        assert pair_reads_output == pair_reads_expected
    finally:
        os.remove(outfile)
        os.remove(logfile)

def test_pair_reads_and_save_mp():
    try:
        file_id = datetime.datetime.now().strftime("%Y%m%d%H%M%S%f")
        outfile = f"test_files/test_pair_out_{file_id}.fastq"
        logfile = f"test_files/test_pair_out_{file_id}.log"
        pair_reads_and_save_mp("test_files/pair_test_r1.fastq", "test_files/pair_test_r2.fastq", outfile, logfile, align_threshold=0.8, cpus=8, seq_only=False)
        with open(outfile, 'r') as f:
            pair_reads_output = f.readlines()
        with open("test_files/paired_expected.fastq", 'r') as f:
            pair_reads_expected = f.readlines()
        assert pair_reads_output == pair_reads_expected
    finally:
        os.remove(outfile)
        os.remove(logfile)
        

def test_pair_reads_and_save_seq_only():
    try:
        file_id = datetime.datetime.now().strftime("%Y%m%d%H%M%S%f")
        outfile = f"test_files/test_pair_out_{file_id}.fastq"
        logfile = f"test_files/test_pair_out_{file_id}.log"
        pair_reads_and_save("test_files/pair_test_r1.fastq", "test_files/pair_test_r2.fastq", outfile, logfile, align_threshold=0.8, seq_only=True)
        with open(outfile, 'r') as f:
            pair_reads_output = f.readlines()
        with open("test_files/paired_expected_seq_only.fastq", 'r') as f:
            pair_reads_expected = f.readlines()
        assert pair_reads_output == pair_reads_expected
    finally:
        os.remove(outfile)
        os.remove(logfile)


def test_pair_reads_and_save_log():
    try:
        file_id = datetime.datetime.now().strftime("%Y%m%d%H%M%S%f")
        outfile = f"test_files/test_pair_out_{file_id}.fastq"
        logfile = f"test_files/test_pair_out_{file_id}.log"
        pair_reads_and_save("test_files/pair_test_r1.fastq", "test_files/pair_test_r2.fastq", outfile, logfile, align_threshold=0.8, seq_only=False)
        with open(logfile, 'r') as f:
            pair_reads_log = f.readlines()
        with open("test_files/paired_expected.log", 'r') as f:
            pair_reads_expected = f.readlines()
        assert pair_reads_log == pair_reads_expected
    finally:
        os.remove(outfile)
        os.remove(logfile)



############################# Test get_cluster_counts.py #############################

@pytest.mark.parametrize("l, expected", [
    ([1], 1),
    ([1, 2, 2, 3, 4, 5, 6, 6, 6, 8],  6),
    ([1, 2, 2, 2, 3, 4, 5, 6, 6, 6, 8],  2),
    (["A", "B", "A", "C"],  "A"),
])
def test_most_common(l, expected):
    assert most_common(l) == expected

@pytest.mark.parametrize("seqs, expected", [
    (["AAAAAAAA"],  "AAAAAAAA"),
    (["AAAAAAAA", "AAAAAAAA", "TTTTTTTT"],  "AAAAAAAA"),
    (["AAAAAAAA", "TAAAAAAT", "TTTTTTTT"],  "TAAAAAAT"),
    (["AAAAAAAA", "AAAAAAAA", "TTTTTTTTT", "TTTTTTTTT"],  "AAAAAAAA"),
    (["ACGTACGTACGT", "ACGTACGTACGT", "ACCTACGTACGTAA", "ACGTACTTACGA", "TACGTACGT"],  "ACGTACGTACGT"),
])
def test_get_consensus_cluster(seqs, expected):
    assert get_consensus_cluster(seqs) == expected

@pytest.mark.parametrize("seqs, max_dist, expected", [
    (["AAAAAAAA", "AAAAATAA", "AGAAAAAA", "AAAAAAAN", "AAAAAAA"],  3, ("AAAAAAAA", 5, set())),
    (["AAAAAAAA", "AAAAAAAA", "TTTTTTTT"], 3, ("AAAAAAAA", 2, {2})),
    (["AAAAAAAA", "TAAAAAAT", "TTTTTTTT"], 2,  ("TAAAAAAT", 2, {2})),
])
def test_get_consensus_and_count(seqs, max_dist, expected):
    assert get_consensus_and_count(seqs, max_dist) == expected

@pytest.mark.parametrize("header, dist, l, expected", [
    ("@headerACGTACGT+GGGGGGGG", 17, 8, "ACGTACGT"),
    ("@headerACGTAC+GGGGGG", 13, 6, "ACGTAC")
])
def test_extract_umi(header, dist, l, expected):
    assert extract_umi(header, dist, l) == expected

def test_extract_umi_ValueError():
    header = "@headerACGTACGT+GGGGGGGG"
    with pytest.raises(ValueError, match=re.escape("UMI contains non-ACTGN character: CGTACGT+")):
        extract_umi(header, 16, 8)

@pytest.mark.parametrize("umis, t, expected", [
    (["AAAAAAAA", "AAAAAAAN", "AAAAAAAA"],  2, [["AAAAAAAA".encode(), "AAAAAAAN".encode()]]),
    (["AAAAAAAA", "AAAAAAAA", "TAAAAAAA", "AAAAAATT", "TTTTTTTT"], 1, [["AAAAAAAA".encode(), "TAAAAAAA".encode()], ["AAAAAATT".encode()],  ["TTTTTTTT".encode()]])
])
def test_cluster_umis(umis, t, expected):
    assert cluster_umis(umis, t, clust_method="directional") == expected

def test_cluster_umis_ValueError():
    with pytest.raises(ValueError, match=re.escape("clust_method must be 'directional', 'cluster', or 'adjacency'. You provided test.")):
        cluster_umis(['AAAAAAAA'], 1, clust_method="test")

def test_get_cluster_counts_main_no_umi():
    try:
        fq = "test_files/get_cluster_counts_test.fastq"
        clustered = "test_files/get_cluster_counts_test.txt"
        file_id = datetime.datetime.now().strftime("%Y%m%d%H%M%S%f")
        outfile = f"test_files/test_cluster_out_{file_id}.txt"
        logfile = f"test_files/test_cluster_{file_id}.log"
        expected_out_no_umi = f"test_files/test_cluster_expected_out_no_umi.txt"
        test_args = ['-f', fq, '-c', clustered, '-o', outfile, '-l', logfile, '-d', '3']
        get_cluster_counts_main(test_args)
        with open(outfile, 'r') as f:
            output_no_umi = f.readlines()
        with open(expected_out_no_umi, 'r') as f:
            output_expected_no_umi = f.readlines()
        assert output_no_umi == output_expected_no_umi
    finally:
        os.remove(outfile)
        os.remove(logfile)


def test_get_cluster_counts_main_umi():
    try:
        fq = "test_files/get_cluster_counts_test.fastq"
        clustered = "test_files/get_cluster_counts_test.txt"
        file_id = datetime.datetime.now().strftime("%Y%m%d%H%M%S%f")
        outfile = f"test_files/test_cluster_out_{file_id}.txt"
        logfile = f"test_files/test_cluster_{file_id}.log"
        expected_out_umi = f"test_files/test_cluster_expected_out_umi.txt"
        test_args = ['-f', fq, '-c', clustered, '-o', outfile, '-l', logfile, '-d', '3', '--umi', '--umi-threshold', '1', '--umi-start', '17', '--umi-len', '8']
        get_cluster_counts_main(test_args)
        with open(outfile, 'r') as f:
            output_umi = f.readlines()
        with open(expected_out_umi, 'r') as f:
            output_expected_umi = f.readlines()
        assert output_umi == output_expected_umi
    finally:
        os.remove(outfile)
        os.remove(logfile)


def test_get_cluster_counts_main_no_umi_in_memory():
    try:
        fq = "test_files/get_cluster_counts_test.fastq"
        clustered = "test_files/get_cluster_counts_test.txt"
        file_id = datetime.datetime.now().strftime("%Y%m%d%H%M%S%f")
        outfile = f"test_files/test_cluster_out_{file_id}.txt"
        logfile = f"test_files/test_cluster_{file_id}.log"
        expected_out_no_umi = f"test_files/test_cluster_expected_out_no_umi.txt"
        test_args = ['-f', fq, '-c', clustered, '-o', outfile, '-l', logfile, '-d', '3', '--in-memory']
        get_cluster_counts_main(test_args)
        with open(outfile, 'r') as f:
            output_no_umi = f.readlines()
        with open(expected_out_no_umi, 'r') as f:
            output_expected_no_umi = f.readlines()
        assert output_no_umi == output_expected_no_umi
    finally:
        os.remove(outfile)
        os.remove(logfile)


def test_get_cluster_counts_main_umi_in_memory():
    try:
        fq = "test_files/get_cluster_counts_test.fastq"
        clustered = "test_files/get_cluster_counts_test.txt"
        file_id = datetime.datetime.now().strftime("%Y%m%d%H%M%S%f")
        outfile = f"test_files/test_cluster_out_{file_id}.txt"
        logfile = f"test_files/test_cluster_{file_id}.log"
        expected_out_umi = f"test_files/test_cluster_expected_out_umi.txt"
        test_args = ['-f', fq, '-c', clustered, '-o', outfile, '-l', logfile, '-d', '3', '--umi', '--umi-threshold', '1', '--umi-start', '17', '--umi-len', '8', '--in-memory']
        get_cluster_counts_main(test_args)
        with open(outfile, 'r') as f:
            output_umi = f.readlines()
        with open(expected_out_umi, 'r') as f:
            output_expected_umi = f.readlines()
        assert output_umi == output_expected_umi
    finally:
        os.remove(outfile)
        os.remove(logfile)
