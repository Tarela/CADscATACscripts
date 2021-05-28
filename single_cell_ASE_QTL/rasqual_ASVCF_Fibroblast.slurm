#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=10
#SBATCH --mem=60GB
#SBATCH --time=14:00:00
#SBATCH --partition=parallel
#SBATCH --account=cphg-millerlab
module load gcc/9.2.0
module load tabix
module load samtools
module load rasqual
bash ./createASVCF.sh paired_end bam.list.Fibroblast.txt UVA_hg38_coronary_scATAC_fibroblast.vcf.gz UVA_hg38_coronary_scATAC_fibroblast_AS.vcf.gz atac
