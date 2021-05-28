#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=10
#SBATCH --mem=40000
#SBATCH --time=12:00:00
#SBATCH --partition=parallel
#SBATCH --account=cphg-millerlab
for i in {1..22}
do
python ldsc.py --l2 --bfile LDSCORE/1000G_EUR_Phase3_plink/1000G.EUR.QC.${i} --ld-wind-cm 1 --annot annot/Endothelial.chr${i}.annot.gz --thin-annot --out Endothelial.chr${i} --print-snps LDSCORE/hapmap3_snps/hm.${i}.snp
done
