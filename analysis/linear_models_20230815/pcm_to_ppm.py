# This script converts a motif PCM to a PPM. The PCM contains counts of each
# base for each position, while the PPM has probabilities of each base for
# each position (sum to 1). The input file is expected to have a header line
# that starts with a >, then each tab-separated row contains the A, C, G, and T
# counts for each position. The output format has rows and columns swapped.
# The first column contains the bases, then subsequent columns represent each
# position in the sequence.

import csv
import argparse
from math import isclose

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--pcm', dest='pcm_fp', help="Input pcm file", required=True)
parser.add_argument('-o', '--output_path', dest='output_path', help="Output path; file is named <header>.ppm.", required=True)
args = parser.parse_args()

with open(args.pcm_fp, 'rt') as f:
    reader = csv.reader(f, delimiter='\t')
    header = next(reader)[0]
    if header[0] != '>':
        raise ValueError("Header does not start with '>'.")
    pcm_mat = []
    for line in reader:
        pcm_mat.append([float(i) for i in line])

output_fp = args.output_path + '/' + header[1:] + '.ppm'

ppm_mat = []

# Divide each row by its sum.
for row in pcm_mat:
    row_sum = sum(row)
    ppm_row = [i / row_sum for i in row]
    if not isclose(sum(ppm_row), 1):
        raise ValueError(f"PPM row does not add to 1: {ppm_row}. Sum = {sum(ppm_row)}")
    ppm_mat.append([str(i) for i in ppm_row])

# Rotate PPM so that cols are positions and rows are bases.
ppm_mat_rotated = [[i for i in col] for col in zip(['A', 'C', 'G', 'T'], *ppm_mat)]

with open(output_fp, 'wt') as f:
    for row in ppm_mat_rotated:
        f.write('\t'.join(row) + '\n')