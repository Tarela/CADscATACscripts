#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=10
#SBATCH --mem=40000
#SBATCH --time=2:00:00
#SBATCH --partition=parallel
#SBATCH --account=cphg-millerlab
for i in {1..22}
do
python make_annot.py \
--bed-file annot/annot_bed/Fibroblast_peaks_hg19.bed \
--bimfile LDSCORE/1000G_EUR_Phase3_plink/1000G.EUR.QC.${i}.bim \
--annot-file Fibroblast.chr${i}.annot.gz
done
