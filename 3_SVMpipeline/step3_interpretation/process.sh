#!/bin/bash

C=$1
maindir=$2
lsgkm_installed_folder=$3
lsgkm_dir=${lsgkm_installed_folder}/src/
i=${C}

python run_ism.py ${C} 20 ${maindir} ${lsgkm_dir}
for j in {0..9}
    do
        paste ${maindir}/step3_interpretation/ism_scores/$i/fold$j.effect.scores ${maindir}/step3_interpretation/ism_scores/$i/fold$j.noneffect.scores | awk 'BEGIN {OFS="\t"}{print $1, $2-$4}' > ${maindir}/step3_interpretation/ism_scores/$i/fold$j.ism.scores
    done

bash schedule_gkmexplain.sh $i ${maindir}   # Run GkmExplain
bash concat_explain.sh $i          # Concatenate GkmExplain Scores

python score_kmers.py ${C} 20 ${maindir} ${lsgkm_dir}
python run_deltasvm_celltype.py ${C} 20 ${maindir}

# explain and plot
python snp_interpretation.py ${C} ${maindir}