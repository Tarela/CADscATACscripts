#!/bin/bash
maindir=$1
annotation=${maindir}/annotation/
ref_fasta=${annotation}/hg38.fa
chrom_sizes=${annotation}/hg38.chrom.sizes

python get_genomewide_gc_bins.py --ref_fasta $ref_fasta --out_prefix gc_hg38_nosmooth --output_format tsv --chrom_sizes $chrom_sizes

