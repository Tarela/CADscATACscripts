#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=60000
#SBATCH --time=3-0:00:00
#SBATCH --partition=standard
#SBATCH --account=cphg-millerlab
module load gcc/9.2.0
module load samtools
module load tabix
module load rasqual
bash ./createASVCF.sh paired_end bam.list.bulk_2.txt UVA_hg38_coronary_bulk_ATAC_2_filtered.vcf.gz bulk_coronary_ATAC_AS_subset35.vcf.gz atac
