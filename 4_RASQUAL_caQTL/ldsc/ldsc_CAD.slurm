#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=10
#SBATCH --mem=40000
#SBATCH --time=12:00:00
#SBATCH --partition=parallel
#SBATCH --account=cphg-millerlab
python ldsc.py \
--h2-cts CAD_GWAS.sumstats.gz \
--ref-ld-chr LDSCORE/1000G_EUR_Phase3_baseline/baseline. \
--out scATAC_CAD \
--ref-ld-chr-cts sc_atac.ldcts \
--w-ld-chr LDSCORE/weights_hm3_no_hla/weights.
