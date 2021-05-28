#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=10
#SBATCH --mem=40000
#SBATCH --time=12:00:00
#SBATCH --partition=parallel
#SBATCH --account=cphg-millerlab
module load gcc/9.2.0
module load gsl/2.4
module load rasqual
module load tabix
module load gparallel/20170822 
n_lines=$(wc -l SMC_input_3/rasqual.SMC.input.txt | cut -d" " -f1)
parallel -j10 bash rasqual_SMC.sh \
SMC_input_3/rasqual.SMC.input.txt {} \
SMC_input_3/Y.bin \
SMC_input_3/K.bin \
SMC_input_3/X.bin \
asvcf/UVA_hg38_coronary_scATAC_AS.vcf.gz \
SMC_output_3 ::: `seq $n_lines`
