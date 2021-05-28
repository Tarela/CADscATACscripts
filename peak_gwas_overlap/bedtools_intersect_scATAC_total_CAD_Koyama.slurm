#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=20000
#SBATCH --time=1:00:00
#SBATCH --partition=standard
#SBATCH --account=cphg-millerlab
module load gcc/7.1.0
module load bedtools
bedtools intersect -a GWAS_bed_files/Koyama_CAD_SNPs_hg38_100bp.bed -b allpeak_celltype/Endothelial_peaks.bed -wa > cad_overlap_Koyama_snplist_Endothelial_total.bed
bedtools intersect -a GWAS_bed_files/Koyama_CAD_SNPs_hg38_100bp.bed -b allpeak_celltype/Endothelial_peaks.bed -wb > cad_overlap_Koyama_peaklist_Endothelial_total.bed
bedtools intersect -a GWAS_bed_files/Koyama_CAD_SNPs_hg38_100bp.bed -b allpeak_celltype/Fibroblast_peaks.bed -wa > cad_overlap_Koyama_snplist_Fibroblast_total.bed
bedtools intersect -a GWAS_bed_files/Koyama_CAD_SNPs_hg38_100bp.bed -b allpeak_celltype/Fibroblast_peaks.bed -wb > cad_overlap_Koyama_peaklist_Fibroblast_total.bed
bedtools intersect -a GWAS_bed_files/Koyama_CAD_SNPs_hg38_100bp.bed -b allpeak_celltype/Macrophage_peaks.bed -wa > cad_overlap_Koyama_snplist_Macrophage_total.bed
bedtools intersect -a GWAS_bed_files/Koyama_CAD_SNPs_hg38_100bp.bed -b allpeak_celltype/Macrophage_peaks.bed -wb > cad_overlap_Koyama_peaklist_Macrophage_total.bed
bedtools intersect -a GWAS_bed_files/Koyama_CAD_SNPs_hg38_100bp.bed -b allpeak_celltype/Mast_peaks.bed -wa > cad_overlap_Koyama_snplist_Mast_total.bed
bedtools intersect -a GWAS_bed_files/Koyama_CAD_SNPs_hg38_100bp.bed -b allpeak_celltype/Mast_peaks.bed -wb > cad_overlap_Koyama_peaklist_Mast_total.bed
bedtools intersect -a GWAS_bed_files/Koyama_CAD_SNPs_hg38_100bp.bed -b allpeak_celltype/Pericyte_peaks.bed -wa > cad_overlap_Koyama_snplist_Pericyte_total.bed
bedtools intersect -a GWAS_bed_files/Koyama_CAD_SNPs_hg38_100bp.bed -b allpeak_celltype/Pericyte_peaks.bed -wb > cad_overlap_Koyama_peaklist_Pericyte_total.bed
bedtools intersect -a GWAS_bed_files/Koyama_CAD_SNPs_hg38_100bp.bed -b allpeak_celltype/Plasma_peaks.bed -wa > cad_overlap_Koyama_snplist_Plasma_total.bed
bedtools intersect -a GWAS_bed_files/Koyama_CAD_SNPs_hg38_100bp.bed -b allpeak_celltype/Plasma_peaks.bed -wb > cad_overlap_Koyama_peaklist_Plasma_total.bed
bedtools intersect -a GWAS_bed_files/Koyama_CAD_SNPs_hg38_100bp.bed -b allpeak_celltype/SMC_peaks.bed -wa > cad_overlap_Koyama_snplist_SMC_total.bed
bedtools intersect -a GWAS_bed_files/Koyama_CAD_SNPs_hg38_100bp.bed -b allpeak_celltype/SMC_peaks.bed -wb > cad_overlap_Koyama_peaklist_SMC_total.bed
bedtools intersect -a GWAS_bed_files/Koyama_CAD_SNPs_hg38_100bp.bed -b allpeak_celltype/T_peaks.bed -wa > cad_overlap_Koyama_snplist_T_total.bed
bedtools intersect -a GWAS_bed_files/Koyama_CAD_SNPs_hg38_100bp.bed -b allpeak_celltype/T_peaks.bed -wb > cad_overlap_Koyama_peaklist_T_total.bed
bedtools intersect -a GWAS_bed_files/Koyama_CAD_SNPs_hg38_100bp.bed -b allpeak_celltype/unknown_peaks.bed -wa > cad_overlap_Koyama_snplist_unknown_total.bed
bedtools intersect -a GWAS_bed_files/Koyama_CAD_SNPs_hg38_100bp.bed -b allpeak_celltype/unknown_peaks.bed -wb > cad_overlap_Koyama_peaklist_unknown_total.bed
