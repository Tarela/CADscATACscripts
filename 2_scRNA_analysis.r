
# scRNAseq data processing, mainly with Seurat package, version4.0.0
library(dplyr)
library(Seurat)
library(patchwork)
set.seed(1)

# the input data GSE131778_human_coronary_scRNAseq_wirka_et_al_GEO.txt was downloaded from GEO, GSE131778, https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131778/suppl/GSE131778%5Fhuman%5Fcoronary%5FscRNAseq%5Fwirka%5Fet%5Fal%5FGEO%2Etxt%2Egz
expmat <- read.table("GSE131778_human_coronary_scRNAseq_wirka_et_al_GEO.txt",row.names=1,header=T)
data_raw <- CreateSeuratObject(counts = expmat, project = "CAscRNA", min.cells = 3, min.features = 200)

data_raw[["percent.mt"]] <- PercentageFeatureSet(data_raw, pattern = "^MT-")
pdf(file="scRNA_basicQC.pdf",width=12,height=6)
VlnPlot(data_raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(data_raw, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data_raw, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()

# QC cell filtering
genes_expCellNum <- apply(data_raw@assays$RNA@counts,1,function(x) return( length(which(x>0)) ))
data_filter_tmp <- data_raw[which(genes_expCellNum >=5),]
data_filter <- subset(data_filter_tmp, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 7.5)

# normalization
data_normalize <- NormalizeData(data_filter, normalization.method = "LogNormalize", scale.factor = 10000)

# highvar genes selection
data_highVar <- FindVariableFeatures(data_normalize, selection.method = "vst", nfeatures = 2000)
pdf(file="highVar_genes_selection_scatter.pdf",width=12,height=6)
plot1 <- VariableFeaturePlot(data_highVar)
plot2 <- LabelPoints(plot = plot1, points = top10_highVar_genes, repel = TRUE)
plot1 + plot2
dev.off()

# scaling data
all.genes <- rownames(data_highVar)
data_scale <- ScaleData(data_highVar, features = all.genes)

# dimensional reduction
data_PCA <- RunPCA(data_scale, features = VariableFeatures(object = data_scale))
print(data_PCA[["pca"]], dims = 1:5, nfeatures = 5)
pdf(file="PCA_loading_scatter.pdf",width=12,height=6)
VizDimLoadings(data_PCA, dims = 1:2, reduction = "pca")
dev.off()

pdf(file="PCA_reduction_scatter.pdf")
DimPlot(data_PCA, reduction = "pca")
dev.off()

pdf(file="PCA_loading_heatmap.pdf",width=12,height=6)
DimHeatmap(data_PCA, dims = 1:2, cells = 500, balanced = TRUE)
dev.off()

# determine the dim
data_PCA_tmp <- JackStraw(data_PCA, num.replicate = 100,dims=1:25)
data_PCA_ex <- ScoreJackStraw(data_PCA_tmp, dims = 1:25)
pdf(file="PCA_DiffDim_cmp_v1.pdf",width=12,height=6)
p1 <- JackStrawPlot(data_PCA_ex, dims = 1:15)
p2 <- ElbowPlot(data_PCA_ex)
p1+p2
dev.off()

proportion_PCA_std <- data_PCA@reductions$pca@stdev/sum(data_PCA@reductions$pca@stdev)
cdf_PCA_std <- c()
cdf_this <- 0
for(i in 1:length(proportion_PCA_std)){
    cdf_this <- cdf_this+proportion_PCA_std[i]
    cdf_PCA_std <- c(cdf_PCA_std, cdf_this)
}

pdf(file="PCA_DiffDim_cmp_v2.pdf",width=12,height=6)
par(mfrow=c(1,2),mar=c(4,4,2,2))
plot(data_PCA@reductions$pca@stdev,pch=16,ylab="std",xlab="PC")
plot(cdf_PCA_std,pch=16,ylab="CDF proportion of std",xlab="PC")
abline(h=0.5,col="red")
dev.off()

#  clustering
data_KNN_K10 <- FindNeighbors(data_PCA, dims = 1:10)
data_cluster_K10 <- FindClusters(data_KNN_K10, resolution = 0.5, random.seed=1)

# UMAP/tsne embbeding
data_UMAP_K10 <- RunUMAP(data_cluster_K10, dims = 1:10, seed.use=1)
all_marker_genes_K10 <- FindAllMarkers(data_UMAP_K10, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, random.seed=1)
MKG_K10 <- all_marker_genes_K10[which(all_marker_genes_K10[,"avg_logFC"]>=log(2) & all_marker_genes_K10[,"p_val_adj"]<1e-5),]
data_UMAPtsne_K10 <- RunTSNE(data_UMAP_K10, dims = 1:10, seed.use=1)
all_marker_genes_K10 <- FindAllMarkers(data_UMAP_K10, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, random.seed=1)
MKG_K10 <- all_marker_genes_K10[which(all_marker_genes_K10[,"avg_logFC"]>=log(2) & all_marker_genes_K10[,"p_val_adj"]<1e-5),]

# read in marker gene list from downloaded from Wirka et al 2019, GSE131778
Bcell_MKG <- as.vector(read.table("scRNA_GSE131778_makergene_list/Bcell.txt")[,1])
Endothelial_MKG <- as.vector(read.table("scRNA_GSE131778_makergene_list/Endothelial.txt")[,1])
Fibroblast_MKG <- as.vector(read.table("scRNA_GSE131778_makergene_list/Fibroblast.txt")[,1])
Fibroblast2_MKG <- as.vector(read.table("scRNA_GSE131778_makergene_list/Fibroblast2.txt")[,1])
Fibromyocyte_MKG <- as.vector(read.table("scRNA_GSE131778_makergene_list/Fibromyocyte.txt")[,1])
Macrophage_MKG <- as.vector(read.table("scRNA_GSE131778_makergene_list/Macrophage.txt")[,1])
Mast_MKG <- as.vector(read.table("scRNA_GSE131778_makergene_list/Mast.txt")[,1])
NKcell_MKG <- as.vector(read.table("scRNA_GSE131778_makergene_list/NKcell.txt")[,1])
Neuron_MKG <- as.vector(read.table("scRNA_GSE131778_makergene_list/Neuron.txt")[,1])
Pericyte1_MKG <- as.vector(read.table("scRNA_GSE131778_makergene_list/Pericyte1.txt")[,1])
Pericyte2_MKG <- as.vector(read.table("scRNA_GSE131778_makergene_list/Pericyte2.txt")[,1])
Plasma1_MKG <- as.vector(read.table("scRNA_GSE131778_makergene_list/Plasma1.txt")[,1])
Plasma2_MKG <- as.vector(read.table("scRNA_GSE131778_makergene_list/Plasma2.txt")[,1])
SMC_MKG <- as.vector(read.table("scRNA_GSE131778_makergene_list/SMC.txt")[,1])
Tcell_MKG <- as.vector(read.table("scRNA_GSE131778_makergene_list/Tcell.txt")[,1])
c10_MKG <- as.vector(read.table("scRNA_GSE131778_makergene_list/c10.txt")[,1])
c11_MKG <- as.vector(read.table("scRNA_GSE131778_makergene_list/c11.txt")[,1])
c12_MKG <- as.vector(read.table("scRNA_GSE131778_makergene_list/c12.txt")[,1])
c14_MKG <- as.vector(read.table("scRNA_GSE131778_makergene_list/c14.txt")[,1])
c18_MKG <- as.vector(read.table("scRNA_GSE131778_makergene_list/c18.txt")[,1])



## compare clusters with marker genes, generate confusion matrix  

compare_celltype <- function(usedata, outname, MKGdata){
    ### compare with DEG
    DEGovMAT <- matrix(rep(0, 21*length(table(usedata$seurat_clusters))),ncol=21)
    rownames(DEGovMAT) <- names(table(usedata$seurat_clusters))
    colnames(DEGovMAT) <- c("Bcell","Endothelial","Fibroblast","Fibroblast2","Fibromyocyte","Macrophage","Mast","NKcell","Neuron","Pericyte1","Pericyte2","Plasma1","Plasma2","SMC","Tcell","c10","c11","c12","c14","c18","total")
    for(i in rownames(DEGovMAT)){
        MKG_thisCluster <- MKGdata[which(MKGdata[,"cluster"]==i),"gene"]
        DEGovMAT[i,"Bcell"] = length(intersect(MKG_thisCluster, Bcell_MKG))
        DEGovMAT[i,"Endothelial"] = length(intersect(MKG_thisCluster, Endothelial_MKG))
        DEGovMAT[i,"Fibroblast"] = length(intersect(MKG_thisCluster, Fibroblast_MKG))
        DEGovMAT[i,"Fibroblast2"] = length(intersect(MKG_thisCluster, Fibroblast2_MKG))
        DEGovMAT[i,"Fibromyocyte"] = length(intersect(MKG_thisCluster, Fibromyocyte_MKG))
        DEGovMAT[i,"Macrophage"] = length(intersect(MKG_thisCluster, Macrophage_MKG))
        DEGovMAT[i,"Mast"] = length(intersect(MKG_thisCluster, Mast_MKG))
        DEGovMAT[i,"NKcell"] = length(intersect(MKG_thisCluster, NKcell_MKG))
        DEGovMAT[i,"Neuron"] = length(intersect(MKG_thisCluster, Neuron_MKG))
        DEGovMAT[i,"Pericyte1"] = length(intersect(MKG_thisCluster, Pericyte1_MKG))
        DEGovMAT[i,"Pericyte2"] = length(intersect(MKG_thisCluster, Pericyte2_MKG))
        DEGovMAT[i,"Plasma1"] = length(intersect(MKG_thisCluster, Plasma1_MKG))
        DEGovMAT[i,"Plasma2"] = length(intersect(MKG_thisCluster, Plasma2_MKG))
        DEGovMAT[i,"SMC"] = length(intersect(MKG_thisCluster, SMC_MKG))
        DEGovMAT[i,"Tcell"] = length(intersect(MKG_thisCluster, Tcell_MKG))
        DEGovMAT[i,"c10"] = length(intersect(MKG_thisCluster, c10_MKG))
        DEGovMAT[i,"c11"] = length(intersect(MKG_thisCluster, c11_MKG))
        DEGovMAT[i,"c12"] = length(intersect(MKG_thisCluster, c12_MKG))
        DEGovMAT[i,"c14"] = length(intersect(MKG_thisCluster, c14_MKG))
        DEGovMAT[i,"c18"] = length(intersect(MKG_thisCluster, c18_MKG))
        DEGovMAT[i,"total"] = length(MKG_thisCluster)
    }
    write.table(DEGovMAT, file=paste0(outname,"_ownClusterDEG_vs_paperDEG.txt"),row.names=T,col.names=T,sep="\t", quote=F)   
}

compare_celltype(data_UMAP_K10, "PC10", MKG_K10)

# assign cell types to clusters based on confusion matrix comparing marker gene assignment. 
own_celltype <- c("Fibroblast","Endothelial","Macrophage","Fibro","T/NK","SMC","Pericyte1","unknown1","Pericyte2","B","Plasma","unknown2","Neuron","unknown3","Mast")


# check with the knowledge based marker genes for cell  types
DSMCmarker <- c("MYOCD","SRF","TEAD3","TEAD4","ACTA2","MYH11","TAGLN","LMOD1","CNN1","TPM2","MYL9")
MSMCmarker <- c("TCF21","KLF4","FN1","LUM","TNFRSF11B","BGN")
Emarker <- c("KLF2","PECAM1","CLDN5","PLVAP","ACKR1","EGFL7", "NFKB1","NFKB2","VCAM1","SELE")

pdf(file="MKG_boxplot.pdf",width=24,height=24)
VlnPlot(data_UMAP_K10, features = c(DSMCmarker,MSMCmarker,Emarker ))
dev.off()
pdf(file="MKG_scatter.pdf",width=40,height=40)
FeaturePlot(data_UMAP_K10, features = c(DSMCmarker,MSMCmarker,Emarker ))
dev.off()
top10marker <- MKG_K10 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
pdf(file="top10marker_exp_heatmap.pdf")
DoHeatmap(data_UMAP_K10, features = top10marker$gene) + NoLegend()
dev.off()



# output processed scRNA data for integration analysis with scATAC
own_cluster <- as.vector(data_UMAP_K10$seurat_clusters)
names(own_cluster) <- names(data_UMAP_K10$seurat_clusters)

out_celltype <- matrix(rep(0, 2*length(own_cluster)),ncol=2)
rownames(out_celltype) <- names(own_cluster)
colnames(out_celltype)<-c("cluster","celltype")
idx <- rownames(out_celltype)
out_celltype[idx,"cluster"] <- own_cluster[idx]
#out_celltype[idx,"celltype"] <- 
for(i in names(table(out_celltype[,"cluster"]))){
    out_celltype[which(out_celltype[,"cluster"]==i),"celltype"] <- own_celltype[as.numeric(i)+1]
}

write.table(out_celltype,file="PC10_celltype_assignment.txt",row.names=T,col.names=T,sep="\t",quote=F)
saveRDS(data_UMAPtsne_K10, file="scRNA_PC10.rds")

