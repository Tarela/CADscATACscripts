#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=40000
#SBATCH --time=2:00:00
#SBATCH --partition=standard
#SBATCH --account=cphg-millerlab
for i in {1..22}
do
python make_annot.py \
--bed-file annot/annot_bed/SMC_peaks_hg19.bed \
--bimfile LDSCORE/1000G_EUR_Phase3_plink/1000G.EUR.QC.${i}.bim \
--annot-file SMC.chr${i}.annot.gz
done
