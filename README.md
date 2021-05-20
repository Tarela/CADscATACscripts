# Scripts for CAD scATAC data analysis
# Single cell chromatin accessibility landscape of human coronary arteries reveals cell type specific regulatory programs and annotates disease risk variants

This repo contains all the data analysis scripts in the manuscript "Single cell chromatin accessibility landscape of human coronary arteries reveals cell type specific regulatory programs and annotates disease risk variants". In the repo, scATAC_analysis.r is for single cell ATAC-seq data analysis, scRNA_analysis.r is for single cell RNA-seq data analysis, and SVMpipeline/ is for the variant effect predictions using ATAC-seq peaks. See the details in the sections below. 

## 0. Requirements
The following environments/packages are required to run the scripts. 

Single cell data analysis (scATAC_analysis.r and scRNA_analysis.r)
- R v4.0
- ArchR v1.0.1
- Seurat v4.0.0
- dplyr v1.0.4
- patchwork v1.1.1
- colorRamps v2.3

Variant scoring (SVMpipeline)
- [lsgkm](https://github.com/kundajelab/lsgkm)
- python v3.7
- perl v5.16
- numpy v1.20
- pandas v1.1
- scipy v1.4
- IPython v7.21
- matplotlib v3.2
- plotnine v0.7
- two annotation data [hg38.fa](https://www.dropbox.com/s/el46t6onl67ejhh/hg38_clean.fa?dl=0) and [hg38.chrom.sizes](https://www.dropbox.com/s/6pf8473bot9wxlv/hg38.chrom.sizes?dl=0)

## 1. single cell ATAC-seq (scATAC) data analysis
All the scripts for the scATAC data analysis were included in the scATAC_analysis.r

The scATAC data were generated for this study and pre-processed with 10x cellranger-atac pipeline (v1.2.0, with default parameters). The output files from the cellranger-atac pipeline were used as the input for the script. The script contains following parts.
- basicQC
- cell filtering
- dimensional reduction
- cell clustering
- UMAP visualization
- batch effect detection
- scRNA-seq integration
- cell type assignment
- peak calling
- detection of celltype specific genes/peaks
- motif annotation
- chromVAR
- trajectory analysis
- co-accessibility 
- footprint analysis
- generating other plots for the manuscript

## 2. single cell RNA-seq (scRNA) data analysis
All the scripts for the scRNA data analysis were included in the scRNA_analysis.r. 

The scRNA data were generated in a previous study in the same system as scATAC (Wirka et al., Nat Med. 2019, GSE131778). The main propose for the scRNA data analysis is to generate a scRNA reference database for the integration analysis and label transferring analysis of the scATAC data. The script contains following parts.
- cell and gene filtering
- dimensional reduction 
- cell clustering
- UMAP visualization
- comparing the marker genes (Wirka et al., Nat Med. 2019, GSE131778) versus cell clustering 
- cell type assignment based on marker genes
- check knowledge based marker genes versus cell type assignment


## 3. Variant effect predictions using ATAC-seq peaks
All the scripts for the variant scoring analysis were included in the SVMpipeline/ folder. 

The pipeline and the lsgkm packge were inherited from Corces et al., Nat Genet 2020, Ghandi et al., PLoS Comput Biol 2014, Lee et al., Bioinformatics 2016, and Shrikumar et al., bioinformatics 2019. We acknowledge the helps from Dr. Anshul Kundaje and Soumya Kundu at Standford University

\#To run the variant scoring analysis. Users could download the folder to any local machine/server. Next, users should download the annotation data for the reference genome ([hg38.fa](https://www.dropbox.com/s/el46t6onl67ejhh/hg38_clean.fa?dl=0) and [hg38.chrom.sizes](https://www.dropbox.com/s/6pf8473bot9wxlv/hg38.chrom.sizes?dl=0)) and put them into the annotation/ folder. Then get in to the folders for each of the 4 steps in order (step0~step3). In each step, run the process.sh file with specific parameters. See the example below
```sh
$ cd /abspath/yourfolder/forSVMpipeline/step0_background
$ ./process.sh /abspath/yourfolder/forSVMpipeline/
$ cd /abspath/yourfolder/forSVMpipeline/step1_preprocess
$ ./process.sh SMC /abspath/yourfolder/forSVMpipeline/
$ cd /abspath/yourfolder/forSVMpipeline/step2_model/
$ ./process.sh SMC /abspath/yourfolder/forSVMpipeline/ /abspath/yourfolder/forLSGKM/
$ cd /abspath/yourfolder/forSVMpipeline/step3_interpretation/
$ ./process.sh SMC /abspath/yourfolder/forSVMpipeline/ /abspath/yourfolder/forLSGKM/
```

note that the ATAC-seq peaks were located in the celltype_peaks/ folder. Users could change the parameter SMC in the above code to any other cell type listed, or customarize a peak files in the similar format for users' own data and run the pipeline in a similar way. 
