#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=10
#SBATCH --mem=40000
#SBATCH --time=12:00:00
#SBATCH --partition=parallel
#SBATCH --account=cphg-millerlab
python munge_sumstats.py \
--sumstats UKBB/blood_WHITE_COUNT.sumstats.gz \
--merge-alleles LDSCORE/w_hm3.snplist \
--out UKBB_White_Count \
--N 444502 \
--chunksize 500000
