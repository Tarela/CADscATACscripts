#!/bin/bash

task=$1
outdir=$2
ref_fasta=$3

for split in `seq 0 9`
do
    for dataset in train test
    do
    python get_gc_content.py \
           --input_bed $outdir/$task/svm.peaks.$task.$dataset.$split.bed \
           --ref_fasta $ref_fasta \
           --out_prefix $outdir/$task/svm.peaks.$task.$dataset.$split.gc.seq \
           --center_summit \
           --flank_size 500 \
           --store_seq
    done
done

