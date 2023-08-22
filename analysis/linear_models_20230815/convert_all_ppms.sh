#!/usr/bin/bash
in_path=C:/Users/nmateyko/drylab/GSC_2023_STARR-seq/analysis/data/motifs/HOCOMOCOv11_core_pcm_HUMAN_mono/pcm/*.pcm
out_path=C:/Users/nmateyko/drylab/GSC_2023_STARR-seq/analysis/data/motifs/HOCOMOCOv11_core_pcm_HUMAN_mono/ppm/
# for filename in "${motif_path}/pcm/*.pcm"; do
for filename in $in_path; do
    python pcm_to_ppm.py -i $filename -o $out_path
done