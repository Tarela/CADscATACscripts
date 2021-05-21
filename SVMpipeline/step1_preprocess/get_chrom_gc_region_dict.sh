#!/bin/bash

task=$1
outdir=$2

python get_chrom_gc_region_dict.py \
       --input_bed $outdir/$task/$task.candidate.negatives.tsv \
       --outf $outdir/$task/$task.candidate.negatives.gc.p

