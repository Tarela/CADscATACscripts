#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=10
#SBATCH --mem=40000
#SBATCH --time=12:00:00
#SBATCH --partition=parallel
#SBATCH --account=cphg-millerlab
python munge_sumstats.py \
--sumstats sumstats/IMT_summary_for_ldsc.txt.gz \
--merge-alleles LDSCORE/w_hm3.snplist \
--out IMT_GWAS \
--N 60000 \
--chunksize 500000
