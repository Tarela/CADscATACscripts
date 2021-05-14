#!/bin/bash

C=$1
task=${C}
maindir=$2

outdir=${maindir}/step1_preprocess/
annotation=${maindir}/annotation/

peaks=${outdir}/celltype_peaks/${C}_m2PE_peaks.narrowPeak
ref_fasta=${annotation}/hg38.fa
chrom_sizes=${annotation}/hg38.chrom.sizes
genomewide_gc=${outdir}/step0_background/gc_hg38_nosmooth.tsv

ntrain=60000
genome=hg38

[[ -d $outdir/$task ]] || mkdir $outdir/$task

echo "starting $task $peaks"

bash get_svm_peak_splits.sh $task $peaks $outdir $genome $ntrain

echo "got svm peak splits"

bash get_gc_positives.sh $task $outdir $ref_fasta

echo "got gc content of the positive sequences"

bash get_all_negatives.sh $task $peaks $outdir $genomewide_gc

echo "got candidate negative set"

bash get_chrom_gc_region_dict.sh $task $outdir

echo "created python pickle for candidate negatives"

bash form_svm_input_fastas.sh $task $outdir $ref_fasta

echo "finished creating SVM inputs"

