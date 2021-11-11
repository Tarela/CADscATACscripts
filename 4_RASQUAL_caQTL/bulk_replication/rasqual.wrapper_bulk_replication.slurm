#!/bin/bash
#SBATCH -N 1
#SBATCH --cpus-per-task=10
#SBATCH --mem=80000
#SBATCH --time=1-0:00:00
#SBATCH --partition=standard
#SBATCH --account=cphg-millerlab-VIP
module load gcc/9.2.0
module load gsl/2.4
module load rasqual
module load tabix
module load gparallel/20170822 
n_lines=$(wc -l bulk_35_input/rasqual.Bulk.35.input.txt | cut -d" " -f1)
parallel -j10 bash rasqual_bulk_replication.sh \
bulk_35_input/rasqual.Bulk.35.input.txt {} \
bulk_35_input/Y.bin \
bulk_35_input/K.bin \
bulk_35_input/X.bin \
asvcf/bulk_coronary_ATAC_AS_subset35.vcf.gz \
bulk_35_output ::: `seq $n_lines`
