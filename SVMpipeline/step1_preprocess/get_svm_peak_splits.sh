#!/bin/bash

task=$1
peaks=$2
outdir=$3
genome=$4
ntrain=$5

python get_svm_peak_splits.py \
       --narrowPeak $peaks \
       --ntrain $ntrain \
       --out_prefix $outdir/$task/svm.peaks.$task \
       --genome $genome

