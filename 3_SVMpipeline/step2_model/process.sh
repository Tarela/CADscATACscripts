#!/bin/bash
C=$1
task=${C}
maindir=$2
lsgkm_installed_folder=$3

indir=${maindir}/step1_preprocess/
outdir=${maindir}/step2_alliteration/
lsgkm_dir=${lsgkm_installed_folder}/src/

[[ -d $outdir/$task ]] || mkdir $outdir/$task
[[ -d $outdir/$task/models ]] || mkdir $outdir/$task/models
[[ -d $outdir/$task/predictions ]] || mkdir $outdir/$task/predictions

for fold in {0..9}
do
$lsgkm_dir/gkmtrain -m 10000 -v 2 -T 16 $indir/$task/svm.inputs.$task.train.$fold.positives $indir/$task/svm.inputs.$task.train.$fold.negatives $outdir/$task/models/$task.$fold
$lsgkm_dir/gkmpredict -v 2 -T 16 $indir/$task/svm.inputs.$task.test.$fold.positives $outdir/$task/models/$task.$fold.model.txt $outdir/$task/predictions/$task.$fold.positives
$lsgkm_dir/gkmpredict -v 2 -T 16 $indir/$task/svm.inputs.$task.test.$fold.negatives $outdir/$task/models/$task.$fold.model.txt $outdir/$task/predictions/$task.$fold.negatives
done
