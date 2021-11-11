#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=60GB
#SBATCH --time=4-00:00:00
#SBATCH --partition=standard
#SBATCH --account=cphg-millerlab
module load gcc/9.2.0
module load tabix
module load samtools
module load rasqual
bash ./createASVCF.sh paired_end bam.list.txt UVA_hg38_coronary_scATAC_filtered.vcf.gz UVA_hg38_coronary_scATAC_SMC_AS.vcf.gz atac
