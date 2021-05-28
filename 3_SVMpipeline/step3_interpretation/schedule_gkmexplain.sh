#!/bin/bash

cluster=$1
maindir=$2
lsgkm_dir=$3

models=${maindir}/step2_alliteration/
orig_indir=${maindir}/step3_interpretation/explain_inputs
split_indir=${maindir}/step3_interpretation/explain_split_inputs
outdir=${maindir}/step3_interpretation/explain_scores
prefix=


[[ -d $split_indir ]] || mkdir $split_indir
for task in $cluster
do
    [[ -d $split_indir/$prefix$task ]] || mkdir $split_indir/$prefix$task
    cd $split_indir/$prefix$task
    cp $orig_indir/$prefix$task.* $split_indir/$prefix$task/
    split -d -l 100 -a 3 $prefix$task.effect.fasta $prefix$task.effect.fasta.split
    split -d -l 100 -a 3 $prefix$task.noneffect.fasta $prefix$task.noneffect.fasta.split
    rm $split_indir/$prefix$task/$prefix$task.effect.fasta
    rm $split_indir/$prefix$task/$prefix$task.noneffect.fasta
done

#cd /users/soumyak/alzheimers_parkinsons/svm_snp_scoring

[[ -d $outdir/logs ]] || mkdir $outdir/logs

for task in $cluster
do
    [[ -d $outdir/$prefix$task ]] || mkdir $outdir/$prefix$task
    [[ -d $outdir/$prefix$task/split_scores ]] || mkdir $outdir/$prefix$task/split_scores
    for fold in {0..9}
    do
        for allele in effect noneffect
        do
            for split in `ls $split_indir/$prefix$task/$prefix$task.$allele.fasta.split*`
            do
                echo ${lsgkm_dir}/gkmexplain  $split $models/$prefix$task/models/$task.$fold.model.txt $outdir/$prefix$task/split_scores/fold$fold.$allele.scores.split${split:(-3)} >> ${maindir}/step3_interpretation/run_gkmexplain_tmp/${fold}_${allele}_gkmexplain_cmd.sh
            done
        done
        wait
    done
done
