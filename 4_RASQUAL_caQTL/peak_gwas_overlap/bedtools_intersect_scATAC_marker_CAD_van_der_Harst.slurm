#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=20000
#SBATCH --time=1:00:00
#SBATCH --partition=standard
#SBATCH --account=cphg-millerlab-VIP
module load gcc/7.1.0
module load bedtools
bedtools intersect -a GWAS_bed_files/van_der_Harst_CAD_SNPs_hg38_100bp.bed -b markerPeak_celltype/Endothelial_markerPeak.bed -wa > cad_overlap_snplist_Endothelial_marker.bed
bedtools intersect -a GWAS_bed_files/van_der_Harst_CAD_SNPs_hg38_100bp.bed -b markerPeak_celltype/Endothelial_markerPeak.bed -wb > cad_overlap_peaklist_Endothelial_marker.bed
bedtools intersect -a GWAS_bed_files/van_der_Harst_CAD_SNPs_hg38_100bp.bed -b markerPeak_celltype/Fibroblast_markerPeak.bed -wa > cad_overlap_snplist_Fibroblast_marker.bed
bedtools intersect -a GWAS_bed_files/van_der_Harst_CAD_SNPs_hg38_100bp.bed -b markerPeak_celltype/Fibroblast_markerPeak.bed -wb > cad_overlap_peaklist_Fibroblast_marker.bed
bedtools intersect -a GWAS_bed_files/van_der_Harst_CAD_SNPs_hg38_100bp.bed -b markerPeak_celltype/Macrophage_markerPeak.bed -wa > cad_overlap_snplist_Macrophage_marker.bed
bedtools intersect -a GWAS_bed_files/van_der_Harst_CAD_SNPs_hg38_100bp.bed -b markerPeak_celltype/Macrophage_markerPeak.bed -wb > cad_overlap_peaklist_Macrophage_marker.bed
bedtools intersect -a GWAS_bed_files/van_der_Harst_CAD_SNPs_hg38_100bp.bed -b markerPeak_celltype/Mast_markerPeak.bed -wa > cad_overlap_snplist_Mast_marker.bed
bedtools intersect -a GWAS_bed_files/van_der_Harst_CAD_SNPs_hg38_100bp.bed -b markerPeak_celltype/Mast_markerPeak.bed -wb > cad_overlap_peaklist_Mast_marker.bed
bedtools intersect -a GWAS_bed_files/van_der_Harst_CAD_SNPs_hg38_100bp.bed -b markerPeak_celltype/Pericyte_markerPeak.bed -wa > cad_overlap_snplist_Pericyte_marker.bed
bedtools intersect -a GWAS_bed_files/van_der_Harst_CAD_SNPs_hg38_100bp.bed -b markerPeak_celltype/Pericyte_markerPeak.bed -wb > cad_overlap_peaklist_Pericyte_marker.bed
bedtools intersect -a GWAS_bed_files/van_der_Harst_CAD_SNPs_hg38_100bp.bed -b markerPeak_celltype/Plasma_markerPeak.bed -wa > cad_overlap_snplist_Plasma_marker.bed
bedtools intersect -a GWAS_bed_files/van_der_Harst_CAD_SNPs_hg38_100bp.bed -b markerPeak_celltype/Plasma_markerPeak.bed -wb > cad_overlap_peaklist_Plasma_marker.bed
bedtools intersect -a GWAS_bed_files/van_der_Harst_CAD_SNPs_hg38_100bp.bed -b markerPeak_celltype/SMC_markerPeak.bed -wa > cad_overlap_snplist_SMC_marker.bed
bedtools intersect -a GWAS_bed_files/van_der_Harst_CAD_SNPs_hg38_100bp.bed -b markerPeak_celltype/SMC_markerPeak.bed -wb > cad_overlap_peaklist_SMC_marker.bed
bedtools intersect -a GWAS_bed_files/van_der_Harst_CAD_SNPs_hg38_100bp.bed -b markerPeak_celltype/T_markerPeak.bed -wa > cad_overlap_snplist_T_marker.bed
bedtools intersect -a GWAS_bed_files/van_der_Harst_CAD_SNPs_hg38_100bp.bed -b markerPeak_celltype/T_markerPeak.bed -wb > cad_overlap_peaklist_T_marker.bed
bedtools intersect -a GWAS_bed_files/van_der_Harst_CAD_SNPs_hg38_100bp.bed -b markerPeak_celltype/unknown_markerPeak.bed -wa > cad_overlap_snplist_unknown_marker.bed
bedtools intersect -a GWAS_bed_files/van_der_Harst_CAD_SNPs_hg38_100bp.bed -b markerPeak_celltype/unknown_markerPeak.bed -wb > cad_overlap_peaklist_unknown_marker.bed
