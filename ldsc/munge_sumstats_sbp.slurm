#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=10
#SBATCH --mem=40000
#SBATCH --time=2:00:00
#SBATCH --partition=parallel
#SBATCH --account=cphg-millerlab
python munge_sumstats.py \
--sumstats sumstats/sbp_MVP_summary_for_ldsc.txt.gz \
--merge-alleles LDSCORE/w_hm3.snplist \
--out sbp_MVP_GWAS \
--N 300000 \
--chunksize 500000
