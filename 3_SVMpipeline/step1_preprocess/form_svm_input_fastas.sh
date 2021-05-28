#!/bin/bash

task=$1
outdir=$2
ref_fasta=$3

python form_svm_input_fastas.py \
        --outf $outdir/$task/svm.inputs.$task.test.0 $outdir/$task/svm.inputs.$task.test.1 $outdir/$task/svm.inputs.$task.test.2 $outdir/$task/svm.inputs.$task.test.3 $outdir/$task/svm.inputs.$task.test.4 $outdir/$task/svm.inputs.$task.test.5 $outdir/$task/svm.inputs.$task.test.6 $outdir/$task/svm.inputs.$task.test.7 $outdir/$task/svm.inputs.$task.test.8 $outdir/$task/svm.inputs.$task.test.9 $outdir/$task/svm.inputs.$task.train.0 $outdir/$task/svm.inputs.$task.train.1 $outdir/$task/svm.inputs.$task.train.2 $outdir/$task/svm.inputs.$task.train.3 $outdir/$task/svm.inputs.$task.train.4 $outdir/$task/svm.inputs.$task.train.5 $outdir/$task/svm.inputs.$task.train.6 $outdir/$task/svm.inputs.$task.train.7 $outdir/$task/svm.inputs.$task.train.8 $outdir/$task/svm.inputs.$task.train.9 \
        --neg_pickle $outdir/$task/$task.candidate.negatives.gc.p \
        --overwrite_outf \
        --ref_fasta $ref_fasta \
        --peaks $outdir/$task/svm.peaks.$task.test.0.gc.seq $outdir/$task/svm.peaks.$task.test.1.gc.seq $outdir/$task/svm.peaks.$task.test.2.gc.seq $outdir/$task/svm.peaks.$task.test.3.gc.seq $outdir/$task/svm.peaks.$task.test.4.gc.seq $outdir/$task/svm.peaks.$task.test.5.gc.seq $outdir/$task/svm.peaks.$task.test.6.gc.seq $outdir/$task/svm.peaks.$task.test.7.gc.seq $outdir/$task/svm.peaks.$task.test.8.gc.seq $outdir/$task/svm.peaks.$task.test.9.gc.seq $outdir/$task/svm.peaks.$task.train.0.gc.seq $outdir/$task/svm.peaks.$task.train.1.gc.seq $outdir/$task/svm.peaks.$task.train.2.gc.seq $outdir/$task/svm.peaks.$task.train.3.gc.seq $outdir/$task/svm.peaks.$task.train.4.gc.seq $outdir/$task/svm.peaks.$task.train.5.gc.seq $outdir/$task/svm.peaks.$task.train.6.gc.seq $outdir/$task/svm.peaks.$task.train.7.gc.seq $outdir/$task/svm.peaks.$task.train.8.gc.seq $outdir/$task/svm.peaks.$task.train.9.gc.seq

