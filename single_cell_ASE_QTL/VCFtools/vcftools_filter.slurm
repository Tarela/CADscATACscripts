#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=40000
#SBATCH --time=8:00:00
#SBATCH --partition=standard
#SBATCH --account=cphg-millerlab
module load vcftools
module load bcftools
vcftools --gzvcf UVA_hg38_coronary_scATAC.vcf.gz  --maf 0.05  --hwe 0.000001  --remove-indels  --recode  --out UVA_hg38_coronary_scATAC_filtered
