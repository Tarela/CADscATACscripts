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
n_lines=$(wc -l Endothelial_input_3/rasqual.Endothelial.input.txt | cut -d" " -f1)
parallel -j10 bash rasqual.perm_Endothelial.sh \
Endothelial_input_3/rasqual.Endothelial.input.txt {} \
Endothelial_input_3/Y.bin \
Endothelial_input_3/K.bin \
Endothelial_input_3/X.bin \
asvcf/UVA_hg38_coronary_scATAC_endothelial_AS.vcf.gz \
Endothelial_output_3_perm1 ::: `seq $n_lines`
