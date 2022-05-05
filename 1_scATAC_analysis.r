# scATAC data analysis pipeline, mainly with ArchR, v1.0.1
library(ArchR) # version 
set.seed(1)
addArchRGenome("hg38")

# input data from 10x cellranger-atac output
inputFiles <- c(
"scA/fragments.tsv.gz","scB/fragments.tsv.gz","scC/fragments.tsv.gz","scD/fragments.tsv.gz","scE/fragments.tsv.gz","scF/fragments.tsv.gz","scG/fragments.tsv.gz","scH/fragments.tsv.gz","scI/fragments.tsv.gz","scJ/fragments.tsv.gz","scK/fragments.tsv.gz","scL/fragments.tsv.gz","scM/fragments.tsv.gz","scN/fragments.tsv.gz","scO/fragments.tsv.gz","scP/fragments.tsv.gz","scQ/fragments.tsv.gz","scR/fragments.tsv.gz","scS/fragments.tsv.gz","scT/fragments.tsv.gz","scU/fragments.tsv.gz","scV/fragments.tsv.gz","scW/fragments.tsv.gz","scX/fragments.tsv.gz","scY/fragments.tsv.gz","scZ/fragments.tsv.gz","scAA/fragments.tsv.gz","scAB/fragments.tsv.gz","scAC/fragments.tsv.gz","scAD/fragments.tsv.gz","scAE/fragments.tsv.gz","scAF/fragments.tsv.gz","scAG/fragments.tsv.gz","scAH/fragments.tsv.gz","scAI/fragments.tsv.gz","scAJ/fragments.tsv.gz","scAK/fragments.tsv.gz","scAL/fragments.tsv.gz","scAM/fragments.tsv.gz","scAN/fragments.tsv.gz","scAO/fragments.tsv.gz","scAP/fragments.tsv.gz","scAQ/fragments.tsv.gz","scAR/fragments.tsv.gz")
names(inputFiles)<-c("scA","scB","scC","scD","scE","scF","scG","scH","scI","scJ","scK","scL","scM","scN","scO","scP","scQ","scR","scS","scT","scU","scV","scW","scX","scY","scZ","scAA","scAB","scAC","scAD","scAE","scAF","scAG","scAH","scAI","scAJ","scAK","scAL","scAM","scAN","scAO","scAP","scAQ","scAR")

# create ArchR object
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  filterTSS = 4, #Dont set this too high because you can always increase later
  filterFrags = 1000, 
  addTileMat = TRUE,
  force=FALSE,
  addGeneScoreMat = TRUE
)
projCAD1 <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "CAD",
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

# add doublet score
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1,
  force=TRUE
)

# basic QC 
proj_CAD_1 <- projCAD1
p <- ggPoint(
    x = df[,1], 
    y = df[,2], 
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
    ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")
plotPDF(p, name = "TSS-vs-Frags.pdf", ArchRProj = proj_CAD_1, addDOC = FALSE)

p1 <- plotGroups(
    ArchRProj = proj_CAD_1, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "ridges"
   )
p2 <- plotGroups(
    ArchRProj = proj_CAD_1, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )
p3 <- plotGroups(
    ArchRProj = proj_CAD_1, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "log10(nFrags)",
    plotAs = "ridges"
   )
p4 <- plotGroups(
    ArchRProj = proj_CAD_1, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "log10(nFrags)",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )
plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics.pdf", ArchRProj = proj_CAD_1, addDOC = FALSE, width = 4, height = 4)

p1 <- plotFragmentSizes(ArchRProj = proj_CAD_1)
p2 <- plotTSSEnrichment(ArchRProj = proj_CAD_1)
plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = proj_CAD_1, addDOC = FALSE, width = 5, height = 5)


# filter cells
proj_CAD_2 <- filterDoublets(proj_CAD_1,filterRatio=1.5)
idxPass <- which(proj_CAD_2$TSSEnrichment >= 7 & proj_CAD_2$nFrags >= 10000)
df2 <- getCellColData(proj_CAD_2,select = c("log10(nFrags)", "TSSEnrichment"))
cellsPass <- proj_CAD_2$cellNames[idxPass]
proj_CAD_2 <- proj_CAD_2[cellsPass, ]

p <- ggPoint(
    x = df2[,"log10(nFrags)"], 
    y = df2[,"TSSEnrichment"], 
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(3, quantile(df2[,"log10(nFrags)"], probs = 0.99)),
    ylim = c(4, quantile(df2[,"TSSEnrichment"], probs = 0.99))
) + geom_hline(yintercept = 7, lty = "dashed") + geom_vline(xintercept = 4, lty = "dashed")
plotPDF(p, name = "TSS-vs-Frags_cutoff.pdf", ArchRProj = proj_CAD_2, addDOC = FALSE)

# dimensional reduction
proj_CAD_2 <- addIterativeLSI(
    ArchRProj = proj_CAD_2,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30,
    seed=1,force=T
)

# basic clustering 
proj_CAD_2 <- addClusters(
    input = proj_CAD_2,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.8,
    force=T,seed=1
)


# UMAP embedding
proj_CAD_2 <- addUMAP(
    ArchRProj = proj_CAD_2, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine",force=T
)

p1 <- plotEmbedding(ArchRProj = proj_CAD_2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj_CAD_2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters_LSI.pdf", ArchRProj = proj_CAD_2, addDOC = FALSE, width = 5, height = 5)

# QC score projected on UMAP
p1 <- plotEmbedding(ArchRProj = proj_CAD_2, colorBy = "cellColData", name = "DoubletEnrichment", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj_CAD_2, colorBy = "cellColData", name = "DoubletScore", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = proj_CAD_2, colorBy = "cellColData", name = "PromoterRatio", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = proj_CAD_2, colorBy = "cellColData", name = "NucleosomeRatio", embedding = "UMAP")
p5 <- plotEmbedding(ArchRProj = proj_CAD_2, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP")
p6 <- plotEmbedding(ArchRProj = proj_CAD_2, colorBy = "cellColData", name = "nFrags", embedding = "UMAP")
plotPDF(p1,p2,p3,p4,p5,p6, name = "Plot-UMAP-Sample-QC_LSI.pdf", ArchRProj = proj_CAD_2, addDOC = FALSE, width = 5, height = 5)


# check batcheffect
proj_CAD_2 <- addHarmony(
    ArchRProj = proj_CAD_2,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample",force=T
)
proj_CAD_2_HAR <- addClusters(
    input = proj_CAD_2,
    reducedDims = "Harmony",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.8,
    force=T,seed=1
)

cM <- confusionMatrix(paste0(proj_CAD_2$Clusters), paste0(proj_CAD_2$Sample))
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
    mat = as.matrix(cM), 
    color = paletteContinuous("whiteBlue"), 
    border_color = "black"
)
plotPDF(p, name = "confusionMap_heatmap_LSI.pdf", ArchRProj = proj_CAD_2, addDOC = FALSE)

cM <- confusionMatrix(paste0(proj_CAD_2_HAR$Clusters), paste0(proj_CAD_2_HAR$Sample))
cM <- cM / Matrix::rowSums(cM)
p2 <- pheatmap::pheatmap(
    mat = as.matrix(cM), 
    color = paletteContinuous("whiteBlue"), 
    border_color = "black"
)

plotPDF(p, name = "confusionMap_heatmap_LSI.pdf", ArchRProj = proj_CAD_2, addDOC = FALSE)
plotPDF(p2, name = "confusionMap_heatmap_HAR.pdf", ArchRProj = proj_CAD_2, addDOC = FALSE)

p1 <- plotEmbedding(ArchRProj = proj_CAD_2_HAR, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj_CAD_2_HAR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters_HAR.pdf", ArchRProj = proj_CAD_2, addDOC = FALSE, width = 5, height = 5)


# gene score with ArchR default method
proj_CAD_2<-addGeneScoreMatrix(proj_CAD_2,force=TRUE)
proj_CAD_2 <- addImputeWeights(proj_CAD_2,seed=1)

markersGS <- getMarkerFeatures(
    ArchRProj = proj_CAD_2, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

# cluster specific genes
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1")

# knowledge based marker genes for different cell types
DSMCmarker <- c("MYOCD","SRF","TEAD3","TEAD4","ACTA2","MYH11","TAGLN","LMOD1","CNN1","TPM2","MYL9")
MSMCmarker <- c("TCF21","KLF4","FN1","LUM","TNFRSF11B","BGN")
Emarker <- c("KLF2","PECAM1","CLDN5","PLVAP","ACKR1","EGFL7", "NFKB1","NFKB2","VCAM1","SELE")
Tmarker <- c("CD8A","TCF7","RUNX3","TBX21","PRDM1")
Macrophage <- c("CD14","CD36","CD68","CD86","CSF1R","NR1H3","NR1H2","RXRA","RXRB","RXRG","IL1B","CX3CR1")
PericyteMarker <- c("NOTCH3","PDGFRB","RGS5","CSPG4")
Osteochondrogenic <- c("SOX9","RUNX2","BMP2","ALPL")
PotentialSMC <- c("CARMN","PRDM16")
driverSMC <- c("LGALS3","FN1","TNFRSF11B","COL1A1","COL4A1","COL4A2","COL6A3")
useMKG <- intersect(c(DSMCmarker,MSMCmarker,Emarker,Tmarker,Macrophage,PericyteMarker,Osteochondrogenic,PotentialSMC,driverSMC), uniq_symbol)
markerGenes <- useMKG#c(DSMCmarker,MSMCmarker,Emarker,Tmarker)


# integration with scRNAseq data
seRNA <- readRDS("scRNA_PC10.rds");
celltype_meta <- read.table("PC10_celltype_assignment.txt",row.names=1,header=T)

CT1 <- as.vector(celltype_meta[,"celltype"])
CT1[which(celltype_meta[,"celltype"]=="T/NK")] <- "T_NK"
seRNA$celltype <- CT1#celltype_meta[,"celltype"]
names(seRNA$celltype) <- rownames(celltype_meta)
seRNA$celltype <- as.factor(seRNA$celltype)

pal_RNAcelltype <- paletteDiscrete(values = seRNA$celltype)
pal_RNAcelltype[c("Fibroblast","Endothelial","Macrophage","Fibro","T_NK","SMC","Pericyte1","unknown1",
                  "Pericyte2","B","Plasma","unknown2","Neuron","unknown3","Mast")] <- rainbow(15)
p1 <- plotEmbedding(
    proj_CAD_2, 
    colorBy = "cellColData", 
    name = "predictedGroup_Un", 
    pal = pal_RNAcelltype
)
plotPDF(p1, name = "UMAP_RNAIntegration_final.pdf", ArchRProj = proj_CAD_2, addDOC = FALSE, width = 5, height = 5)

proj_CAD_3 <- addGeneIntegrationMatrix(
    ArchRProj = proj_CAD_2, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "Harmony",
    seRNA = seRNA,
    addToArrow = TRUE,
    force= TRUE,
    #groupList = groupList,
    groupRNA = "celltype",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un"
)

# RNA label transfering based celltype assignment
cM <- confusionMatrix(proj_CAD_3$Clusters, proj_CAD_3$predictedGroup_Un)
labelOld <- rownames(cM)
labelNew <- colnames(cM)[apply(cM, 1, which.max)]
labelNew[which(labelNew=="T/NK")] <- "T"
labelNew[which(labelNew=="T_NK")] <- "T"
labelNew[which(labelNew=="Pericyte1")] <- "Pericyte"
proj_CAD_3$Clusters2 <- mapLabels(proj_CAD_3$Clusters, newLabels = labelNew, oldLabels = labelOld)
label_cmp <- cbind(labelOld,labelNew)
rownames(label_cmp) <- label_cmp[,1]
write.table(label_cmp[paste0("C",seq(14)),],file="cluster_celltype_pairName_vsRNA.txt",row.names=F,col.names=F,sep="\t",quote=F)
write.table(as.matrix(cM)[paste0("C",seq(14)),],file="confusionMatrix_scRNA_scATAC_cluster.txt",row.names=T,col.names=T,sep="\t",quote=F)

p1 <- plotEmbedding(proj_CAD_3, colorBy = "cellColData", name = "Clusters2")
plotPDF(p1, name = "UMAP_vsRNA_final.pdf", ArchRProj = proj_CAD_3, addDOC = FALSE, width = 5, height = 5)


# final clustering used in manuscript based on the integration of marker genes and RNA label transfering
labelOld <- paste0("C",seq(14))
labelNew <- c("Endothelial","Fibroblast","Pericyte","SMC","SMC","SMC","SMC","unknown","Macrophage","Macrophage","Plasma","Mast","T","T")
proj_CAD_3$Clusters3 <- mapLabels(proj_CAD_3$Clusters, newLabels = labelNew, oldLabels = labelOld)
label_cmp <- cbind(labelOld,labelNew)
rownames(label_cmp) <- label_cmp[,1]
write.table(label_cmp[paste0("C",seq(14)),],file="cluster_celltype_pairName_final.txt",row.names=F,col.names=F,sep="\t",quote=F)

cM <- confusionMatrix(proj_CAD_3$Clusters3, proj_CAD_3$predictedGroup_Un)
write.table(as.matrix(cM),file="confusionMatrix_scRNA_scATAC_final.txt",row.names=T,col.names=T,sep="\t",quote=F)

p1 <- plotEmbedding(proj_CAD_3, colorBy = "cellColData", name = "Clusters3")
plotPDF(p1, name = "UMAP_celltype_final.pdf", ArchRProj = proj_CAD_3, addDOC = FALSE, width = 5, height = 5)


# peak calling with macs2 v2.1.2
pathToMacs2 <- findMacs2()
proj_CAD_4 <- addGroupCoverages(ArchRProj = proj_CAD_3, groupBy = "Clusters3")

proj_CAD_4 <- addReproduciblePeakSet(
    ArchRProj = proj_CAD_4, 
    groupBy = "Clusters3", 
    pathToMacs2 = pathToMacs2,
    extsize=100,
    cutOff=0.01,
    shift=0,
    extendSummits=200,
    promoterRegion=c(2000,2000),
    genomeSize=3e9,
    reproducibility = "(n+1)/2",
    threads = getArchRThreads()
)
proj_CAD_4 <- addPeakMatrix(proj_CAD_4)
allpeaks <- getPeakSet(proj_CAD_4)
proj_CAD_4 <- addImputeWeights(proj_CAD_4,seed=1)


# celltype specific genes (for manuscript)
markersGS <- getMarkerFeatures(
    ArchRProj = proj_CAD_4, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters3",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1")
printMKG <- function(MKGdata){
  outdata <- matrix(rep(0,nrow(MKGdata)*ncol(MKGdata)),nrow=nrow(MKGdata))
  colnames(outdata) <- colnames(MKGdata)
  for(i in 1:ncol(outdata)){
    outdata[,i] <- as.vector(MKGdata[,i])
  }
  return(outdata)
}
for(Group in names(markerList)){
  write.table(printMKG(markerList[Group][[1]]), file=paste0("markerGene_celltype_final/",Group,"_markerGene.txt"),row.names=F,col.names=T,sep="\t",quote=F)
}

p1heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1", 
  labelMarkers = c(DSMCmarker),
  transpose = TRUE
)
p1 <- ComplexHeatmap::draw(p1heatmapGS, column_title = "DSMC",heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(p1, name = "markerGene_celltype_GSheatmap_final.pdf", ArchRProj = proj_CAD_4, addDOC = FALSE, width = 5, height = 5)


# celltype specific peaks (for manuscript)
dir.create("markerPeak_celltype_final")
markersPeaks <- getMarkerFeatures(
    ArchRProj = proj_CAD_4, 
    useMatrix = "PeakMatrix", 
    groupBy = "Clusters3",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
saveRDS(markersPeaks, file="markerPeak_celltype_final/markerPeak_celltype_final.rds")
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)

output_diffpeak <- function(this_data,cellname){
    out_data <- cbind(as.vector(this_data@seqnames),
                      this_data@ranges@start,
                      this_data@ranges@start + this_data@ranges@width,
                      paste0(cellname,"_",seq(length(this_data@seqnames))),
                      this_data$Log2FC,
                      this_data$FDR)
    write.table(out_data,file=paste0("markerPeak_celltype_final/",cellname,"_markerPeak.bed"),quote=F,row.names=F,col.names=F,sep="\t")
}
output_diffpeak(markerList$Endothelial,"Endothelial")
output_diffpeak(markerList$Fibroblast,"Fibroblast")
output_diffpeak(markerList$Mast,"Mast")
output_diffpeak(markerList$Pericyte,"Pericyte")
output_diffpeak(markerList$Plasma,"Plasma")
output_diffpeak(markerList$SMC,"SMC")
output_diffpeak(markerList$T,"T")
output_diffpeak(markerList$Macrophage,"Macrophage")
output_diffpeak(markerList$unknown,"unknown")

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1",
  transpose = TRUE
)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPeaks, name = "markerPeak_celltype_SIGheatmap_final.pdf", width = 8, height = 6, ArchRProj = proj_CAD_4, addDOC = FALSE)

celltype <- "SMC"
pma <- plotMarkers(seMarker = markersPeaks, name =celltype, cutOff = "FDR <= 0.01 & Log2FC >= 1", plotAs = "MA")
pv <- plotMarkers(seMarker = markersPeaks, name =celltype, cutOff = "FDR <= 0.01 & Log2FC >= 1", plotAs = "Volcano")
plotPDF(pma, pv, name = paste0("markerPeak_scatterMaVolcano_",celltype,"_final.pdf"), width = 5, height = 5, ArchRProj = proj_CAD_4, addDOC = FALSE)


# motif annotation
proj_CAD_5 <- addMotifAnnotations(ArchRProj = proj_CAD_4, motifSet = "homer", name = "Motif",force = TRUE)
enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = proj_CAD_5,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEM, name = "enrichMotifs_homer_celltype_Heatmap_final.pdf", width = 8, height = 6, ArchRProj = proj_CAD_5, addDOC = FALSE)
usemotif <- unlist(lapply(colnames(heatmapEM@matrix), function(x) return(unlist(strsplit(x," ")[[1]][1]))))#[,1]
raw_padj_mat <- enrichMotifs@assays$data$mlog10Padj[usemotif,rownames(heatmapEM@matrix)]
write.table(enrichMotifs@assays$data$mlog10Padj,file="markerPeak_celltype_final/enrichMotifs_homer_padj_celltype_final.txt", row.names=T,col.names=T,sep="\t",quote=F)


# chromVAR
proj_CAD_6 <- addBgdPeaks(proj_CAD_5)

proj_CAD_6 <- addDeviationsMatrix(
  ArchRProj = proj_CAD_6, 
  peakAnnotation = "Motif",
  force = TRUE
)
plotVarDev <- getVarDeviations(proj_CAD_6, name = "MotifMatrix", plot = TRUE)
plotPDF(plotVarDev, name = "chromVAR_VariableMotifDeviationScores_final.pdf", width = 5, height = 5, ArchRProj = proj_CAD_6, addDOC = FALSE)

# plot key motifs/TFs as example
motifs <- c("TEAD3", "TEAD4", "TCF21", "MYOCD", "KLF4", "CARG", "SMAD3","ATF3","SPIB")

markerMotifs <- getFeatures(proj_CAD_6, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")

p <- plotGroups(ArchRProj = proj_CAD_6, 
  groupBy = "Clusters3", 
  colorBy = "MotifMatrix", 
  name = markerMotifs,
  imputeWeights = getImputeWeights(proj_CAD_6)
)

p2 <- lapply(seq_along(p), function(x){
  if(x != 1){
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6) +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
    theme(
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
    ) + ylab("")
  }else{
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6) +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
    theme(
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
    ) + ylab("")
  }
})
do.call(cowplot::plot_grid, c(list(nrow = 1, rel_widths = c(2, rep(1, length(p2) - 1))),p2))

plotPDF(p, name = "chromVAR_GroupsDeviationsImputation_final.pdf", width = 5, height = 5, ArchRProj = proj_CAD_6, addDOC = FALSE)

p <- plotEmbedding(
    ArchRProj = proj_CAD_6, 
    colorBy = "MotifMatrix", 
    name = sort(markerMotifs), 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj_CAD_6)
)

p2 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
    )
})
useP <- do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
plotPDF(useP, name = "chromVAR_UMAP_DeviationsImputation_final.pdf", width = 5, height = 5, ArchRProj = proj_CAD_6, addDOC = FALSE)


# trajectory
trajectory <- c("Endothelial", "SMC", "Fibroblast")
proj_CAD_7 <- addTrajectory(
    ArchRProj = proj_CAD_6, 
    name = "Trajectory", 
    groupBy = "Clusters3",
    trajectory = trajectory, 
    embedding = "UMAP", 
    force = TRUE
)
p <- plotTrajectory(proj_CAD_7, trajectory = "Trajectory", colorBy = "cellColData", name = "Trajectory")
plotPDF(p, name = "Trajectory_celltypeESF_UMAP_final.pdf", width = 5, height = 5, ArchRProj = proj_CAD_7, addDOC = FALSE)

trajMM  <- getTrajectory(ArchRProj = proj_CAD_7, name = "Trajectory", useMatrix = "MotifMatrix", log2Norm = FALSE)
trajGSM <- getTrajectory(ArchRProj = proj_CAD_7, name = "Trajectory", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
trajGIM <- getTrajectory(ArchRProj = proj_CAD_7, name = "Trajectory", useMatrix = "GeneIntegrationMatrix", log2Norm = FALSE)

p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"))
p2 <- trajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "horizonExtra"))
p3 <- plotTrajectoryHeatmap(trajGIM,  pal = paletteContinuous(set = "blueYellow"))

plotPDF(p1, p2, p3, name = "Trajectory_celltypeESF_Heatmaps_final.pdf", ArchRProj = proj_CAD_7, addDOC = FALSE, width = 6, height = 6)


# coAccessibility/cicero
proj_CAD_8 <- addCoAccessibility(
    ArchRProj = proj_CAD_7,
    reducedDims = "IterativeLSI"
)
cA <- getCoAccessibility(
    ArchRProj = proj_CAD_8,
    corCutOff = 0.5,
    resolution = 10000,
    returnLoops = TRUE
)
proj_CAD_8 <- addPeak2GeneLinks(
    ArchRProj = proj_CAD_8,
    reducedDims = "IterativeLSI"
)
p2g <- getPeak2GeneLinks(
    ArchRProj = proj_CAD_8,
    corCutOff = 0.45,
    resolution = 10000,
    returnLoops = TRUE
)

# check the loop calling with examples
useregion <- GRanges(seqnames = "chr3", strand = c("+"),
              ranges = IRanges(start = c(138313800), width = 100000))
p1 <- plotBrowserTrack(
    ArchRProj = proj_CAD_8, 
    groupBy = "Clusters3", 
    region=useregion,
    upstream = 0,
    downstream = 0,
    loops = cA,
    plotSummary=c("geneTrack","loopTrack","bulkTrack"),size=c(4,4,10)
)
p2 <- plotBrowserTrack(
    ArchRProj = proj_CAD_8, 
    groupBy = "Clusters3", 
    region=useregion,
    upstream = 0,
    downstream = 0,
    loops = p2g,
    plotSummary=c("geneTrack","loopTrack","bulkTrack"),size=c(4,4,10)
)
plotPDF(p1, name = "region1_coA_IGV.pdf", ArchRProj = proj_CAD_8, addDOC = FALSE, width = 6, height = 6)
plotPDF(p2, name = "region1_p2g_IGV.pdf", ArchRProj = proj_CAD_8, addDOC = FALSE, width = 6, height = 6)

p <- plotPeak2GeneHeatmap(ArchRProj = proj_CAD_8, groupBy = "Clusters3")
plotPDF(p, name = "p2g_linkHeatmaps_final.pdf", ArchRProj = proj_CAD_8, addDOC = FALSE, width = 6, height = 6)

# output the peak2gene links for plots
write.table(p2g$Peak2GeneLinks,file="peak2geneLinks_final.txt",quote=F,sep="\t",row.names=F)
write.table(cA$CoAccessibility,file="CoAccessibility_final.txt",quote=F,sep="\t",row.names=F)


# footprint analysis
motifPositions <- getPositions(proj_CAD_8)
motifs <- c("TCF21","CARG","TEAD4", "SMAD3","ATF3", "KLF4","GATA1", "CEBPA", "EBF1", "IRF4", "TBX21", "PAX5")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE,ignore.case=TRUE)))

proj_CAD_8 <- addGroupCoverages(ArchRProj = proj_CAD_8, groupBy = "Clusters3")

seFoot <- getFootprints(
  ArchRProj = proj_CAD_8, 
  positions = motifPositions[markerMotifs], 
  groupBy = "Clusters3"
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = proj_CAD_8, 
  normMethod = "Subtract",
  plotName = "Footprints-Subtract-Bias",
  addDOC = FALSE,
  smoothWindow = 5
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = proj_CAD_8, 
  normMethod = "Divide",
  plotName = "Footprints-Divide-Bias",
  addDOC = FALSE,
  smoothWindow = 5
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = proj_CAD_8, 
  normMethod = "None",
  plotName = "Footprints-No-Normalization",
  addDOC = FALSE,
  smoothWindow = 5
)

# output gene score matrix

proj <- proj_CAD_8
GeneScore<-getMatrixFromProject(proj,"GeneScoreMatrix")
GSmat<-GeneScore@assays$data@listData$GeneScoreMatrix
features<-getFeatures(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  select = NULL,
  ignoreCase = TRUE
)

rownames(GSmat)<-features
GSmat_tmp <- as.matrix(GSmat)
imputedGSmat <- imputeMatrix(mat = GSmat_tmp,imputeWeights = getImputeWeights(proj))

saveRDS(GSmat, file="geneScore_matrix_final.rds")
saveRDS(imputedGSmat, file="imputedGeneScore_matrix_final.rds")

GeneInte<-getMatrixFromProject(proj,"GeneIntegrationMatrix")
GImat<-GeneInte@assays$data@listData$GeneIntegrationMatrix
features<-getFeatures(
  ArchRProj = proj,
  useMatrix = "GeneIntegrationMatrix",
  select = NULL,
  ignoreCase = TRUE
)

rownames(GImat)<-features
saveRDS(GImat, file="geneIntegration_matrix_final.rds")

# gene score boxplot
imputGSmat <- imputedGSmat
RImat <- GImat
GSmat <- log10(as.matrix(GSmat)+0.01)
imputGSmat2 <- log10(as.matrix(imputGSmat) + 0.01)
RImat2 <- log10(as.matrix(RImat)+0.001)

GSmat_box <- function(Gname){
    pdf(file=paste0("MKG_GS_boxplot/",Gname,".pdf"),width=12,height=12)
    par(mfcol=c(2,2),mar=c(10,4,2,2))
    boxplot(GSmat[Gname, proj_CAD_8$cellNames[which(proj_CAD_8$Clusters3 == "Endothelial")] ],
            GSmat[Gname, proj_CAD_8$cellNames[which(proj_CAD_8$Clusters3 == "SMC")] ],
            GSmat[Gname, proj_CAD_8$cellNames[which(proj_CAD_8$Clusters3 == "Fibroblast")] ],
            GSmat[Gname, proj_CAD_8$cellNames[which(proj_CAD_8$Clusters3 == "Pericyte")] ],
            GSmat[Gname, proj_CAD_8$cellNames[which(proj_CAD_8$Clusters3 == "Macrophage")] ],
            GSmat[Gname, proj_CAD_8$cellNames[which(proj_CAD_8$Clusters3 == "Mast")] ],
            GSmat[Gname, proj_CAD_8$cellNames[which(proj_CAD_8$Clusters3 == "Plasma")] ],
            GSmat[Gname, proj_CAD_8$cellNames[which(proj_CAD_8$Clusters3 == "T")] ],
            names=c("Endothelial","SMC","Fibroblast","Pericyte","Macrophage","Mast","Plasma","T"),
            col=c("#D51F26","#899ED1","#272E6A","#F47D2C","#1A7738","#89288E","#FFE503","#C06DAC"),
            ylab="geneScore",main=Gname,las=2,cex.axis=2)
    
    boxplot(imputGSmat[Gname, proj_CAD_8$cellNames[which(proj_CAD_8$Clusters3 == "Endothelial")] ],
            imputGSmat[Gname, proj_CAD_8$cellNames[which(proj_CAD_8$Clusters3 == "SMC")] ],
            imputGSmat[Gname, proj_CAD_8$cellNames[which(proj_CAD_8$Clusters3 == "Fibroblast")] ],
            imputGSmat[Gname, proj_CAD_8$cellNames[which(proj_CAD_8$Clusters3 == "Pericyte")] ],
            imputGSmat[Gname, proj_CAD_8$cellNames[which(proj_CAD_8$Clusters3 == "Macrophage")] ],
            imputGSmat[Gname, proj_CAD_8$cellNames[which(proj_CAD_8$Clusters3 == "Mast")] ],
            imputGSmat[Gname, proj_CAD_8$cellNames[which(proj_CAD_8$Clusters3 == "Plasma")] ],
            imputGSmat[Gname, proj_CAD_8$cellNames[which(proj_CAD_8$Clusters3 == "T")] ],
            names=c("Endothelial","SMC","Fibroblast","Pericyte","Macrophage","Mast","Plasma","T"),
            col=c("#D51F26","#899ED1","#272E6A","#F47D2C","#1A7738","#89288E","#FFE503","#C06DAC"),
            ylab="imputted geneScore",main=Gname,las=2,cex.axis=2)
    
    boxplot(GSmat[Gname, proj_CAD_8$cellNames[which(proj_CAD_8$Clusters == "C4")] ],
            GSmat[Gname, proj_CAD_8$cellNames[which(proj_CAD_8$Clusters == "C5")] ],
            GSmat[Gname, proj_CAD_8$cellNames[which(proj_CAD_8$Clusters == "C6")] ],
            GSmat[Gname, proj_CAD_8$cellNames[which(proj_CAD_8$Clusters == "C7")] ],
            names=c("C4","C5","C6","C7"),cex.axis=2,
            ylab="geneScore",main=Gname,las=2)
    
    boxplot(imputGSmat[Gname, proj_CAD_8$cellNames[which(proj_CAD_8$Clusters == "C4")] ],
            imputGSmat[Gname, proj_CAD_8$cellNames[which(proj_CAD_8$Clusters == "C5")] ],
            imputGSmat[Gname, proj_CAD_8$cellNames[which(proj_CAD_8$Clusters == "C6")] ],
            imputGSmat[Gname, proj_CAD_8$cellNames[which(proj_CAD_8$Clusters == "C7")] ],
            names=c("C4","C5","C6","C7"),cex.axis=2,
            ylab="imputted geneScore",main=Gname,las=2)
    
    dev.off()
}
GSmat_box("TCF21")
GSmat_box("MYOCD")
GSmat_box("TEAD4")
GSmat_box("SMAD3")
GSmat_box("ATF3")
GSmat_box("KLF4")

iGSmat_box <- function(Gname){
    pdf(file=paste0("MKG_GS_boxplot/",Gname,"_imputedGS_SMC4clusters.pdf"))
    boxplot(imputGSmat[Gname, proj_CAD_8$cellNames[which(proj_CAD_8$Clusters == "C4")] ],
            imputGSmat[Gname, proj_CAD_8$cellNames[which(proj_CAD_8$Clusters == "C5")] ],
            imputGSmat[Gname, proj_CAD_8$cellNames[which(proj_CAD_8$Clusters == "C6")] ],
            imputGSmat[Gname, proj_CAD_8$cellNames[which(proj_CAD_8$Clusters == "C7")] ],
            names=c("C4","C5","C6","C7"),cex.axis=2,col=c("#9A4364","#E96E50","#EBB422","#98A5D5"),
            ylab="imputted geneScore",main=Gname,las=2)
    dev.off()
}
iGSmat_box("TCF21")
iGSmat_box("MYOCD")
iGSmat_box("MYH11")
iGSmat_box("TAGLN")
iGSmat_box("TNFRSF11B")
iGSmat_box("FN1")
iGSmat_box("CNN1")

RImat_box <- function(Gname){
    pdf(file=paste0("MKG_GS_boxplot/",Gname,"_RNAintegration_SMC4clusters.pdf"))
    boxplot(RImat2[Gname, proj_CAD_8$cellNames[which(proj_CAD_8$Clusters == "C4")] ],
            RImat2[Gname, proj_CAD_8$cellNames[which(proj_CAD_8$Clusters == "C5")] ],
            RImat2[Gname, proj_CAD_8$cellNames[which(proj_CAD_8$Clusters == "C6")] ],
            RImat2[Gname, proj_CAD_8$cellNames[which(proj_CAD_8$Clusters == "C7")] ],
            names=c("C4","C5","C6","C7"),cex.axis=2,col=c("#9A4364","#E96E50","#EBB422","#98A5D5"),
            ylab="log10 RIscore",main=Gname,las=2)
    dev.off()
}
RImat_box("TCF21")
RImat_box("MYOCD")
RImat_box("MYH11")
RImat_box("TAGLN")
RImat_box("TNFRSF11B")
RImat_box("FN1")
RImat_box("CNN1")

Gname <- "TCF21"
a <- (GSmat[Gname, proj_CAD_8$cellNames[which(proj_CAD_8$Clusters3 == "SMC")] ])
pdf(file="TCF21_GS_in_SMC_hist.pdf")
hist(log10(a+0.01),n=200,xlab="log10(TCF21genescore + 0.01)",main="TCF21 GeneScore in SMC")
dev.off()

a4 <- log10((GSmat[Gname, proj_CAD_8$cellNames[which(proj_CAD_8$Clusters == "C4")] ])+0.01)
a5 <- log10((GSmat[Gname, proj_CAD_8$cellNames[which(proj_CAD_8$Clusters == "C5")] ])+0.01)
a6 <- log10((GSmat[Gname, proj_CAD_8$cellNames[which(proj_CAD_8$Clusters == "C6")] ])+0.01)
a7 <- log10((GSmat[Gname, proj_CAD_8$cellNames[which(proj_CAD_8$Clusters == "C7")] ])+0.01)
pdf(file="TCF21_GS_in_C4567_boxplot.pdf")
boxplot(a4,a5,a6,a7,names=c("C4","C5","C6","C7"),ylab="log10(TCF21genescore + 0.01)")
dev.off()

a <- (imputedGSmat[Gname, proj_CAD_8$cellNames[which(proj_CAD_8$Clusters3 == "SMC")] ])
pdf(file="TCF21_iGS_in_SMC_hist.pdf")
hist(log10(a+0.01),n=200,xlab="log10(TCF21genescore + 0.01)",main="TCF21 GeneScore in SMC")
dev.off()

a4 <- log10((imputedGSmat[Gname, proj_CAD_8$cellNames[which(proj_CAD_8$Clusters == "C4")] ])+0.01)
a5 <- log10((imputedGSmat[Gname, proj_CAD_8$cellNames[which(proj_CAD_8$Clusters == "C5")] ])+0.01)
a6 <- log10((imputedGSmat[Gname, proj_CAD_8$cellNames[which(proj_CAD_8$Clusters == "C6")] ])+0.01)
a7 <- log10((imputedGSmat[Gname, proj_CAD_8$cellNames[which(proj_CAD_8$Clusters == "C7")] ])+0.01)
pdf(file="TCF21_iGS_in_C4567_boxplot.pdf")
boxplot(a4,a5,a6,a7,names=c("C4","C5","C6","C7"),ylab="log10(TCF21genescore + 0.01)")
dev.off()



# UMAP colored by sample meta data
newmeta <- read.table("sample_meta.csv",row.names=1,header=T,sep=",")
samplename <- as.vector(proj_CAD_8@cellColData[,"Sample"])

Segment <- rep("NA",length(samplename))
Disease_Category <- rep("NA",length(samplename))
Adventitia <- rep("NA",length(samplename))
Sex <- rep("NA",length(samplename))
Age_Group <- rep("NA",length(samplename))
Type <- rep("NA",length(samplename))

for(S in unique(samplename)){
	Segment[ which(samplename == S) ] <- newmeta[S,"Segment"]
	Disease_Category[ which(samplename == S) ] <- newmeta[S,"Disease_Category"]
	Adventitia[ which(samplename == S) ] <- newmeta[S,"Adventitia"]
	Sex[ which(samplename == S) ] <- newmeta[S,"Sex"]
	Age_Group[ which(samplename == S) ] <- newmeta[S,"Age_Group"]
	Type[ which(samplename == S) ] <- newmeta[S,"Type"]
}

proj_CAD_8$Segment <- Segment
proj_CAD_8$Disease_Category <- Disease_Category
proj_CAD_8$Adventitia <- Adventitia
proj_CAD_8$Sex <- Sex
proj_CAD_8$Age_Group <- Age_Group
proj_CAD_8$Type <- Type

for(metaname in c("Segment","Disease_Category","Adventitia","Sex","Age_Group","Type")){
p1 <- plotEmbedding(proj_CAD_8, colorBy = "cellColData", name = metaname)
plotPDF(p1, name = paste0("UMAP_meta_final_",metaname,".pdf"), ArchRProj = proj_CAD_8, addDOC = FALSE, width = 5, height = 5)	
}


# RNAlabel based diff analysis, Fibromyocyte v.s. SMC
dir.create("markerGenePeak_clusterCMP")
p1 <- plotEmbedding(proj_CAD_8, colorBy = "cellColData", name = "predictedGroup_Un")
plotPDF(p1, name = "test.pdf", ArchRProj = proj_CAD_8, addDOC = FALSE, width = 6, height = 6)
proj_CAD_8$ClustersRNA <- proj_CAD_8@cellColData[,"predictedGroup_Un"]
markersGS_FIBROvSMC <- getMarkerFeatures(
    ArchRProj = proj_CAD_8, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "predictedGroup_Un",
    useGroups = "Fibro",
    bgdGroups = "SMC",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
markersPeaks_FIBROvSMC <- getMarkerFeatures(
    ArchRProj = proj_CAD_8, 
    useMatrix = "PeakMatrix", 
    groupBy = "predictedGroup_Un",
    useGroups = "Fibro",
    bgdGroups = "SMC",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
enrichMotifs_FIBROvSMC <- peakAnnoEnrichment(
    seMarker = markersPeaks_FIBROvSMC,
    ArchRProj = proj_CAD_8,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.01 & Log2FC >= 1"
  )

markerListGene_FIBROvSMC <- getMarkers(markersGS_FIBROvSMC, cutOff = "FDR <= 0.1 & Log2FC >= 1")
markerListPeak_FIBROvSMC <- getMarkers(markersPeaks_FIBROvSMC, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)

this_data <- proj_CAD_8@peakSet
SMCvFIBRO_peak_diffscore <- cbind(as.vector(this_data@seqnames),
                                  this_data@ranges@start,
                                  this_data@ranges@start + this_data@ranges@width,
                                  paste0("peak",seq(length(this_data@ranges@start))),
                                  markersPeaks_FIBROvSMC@assays@data$FDR[[1]],
                                  markersPeaks_FIBROvSMC@assays@data$Log2FC[[1]])
colnames(SMCvFIBRO_peak_diffscore) <- c("chrm","start","end","name","FIBROvSMC_FDR","FIBROvSMC_LFC")
write.table(SMCvFIBRO_peak_diffscore,file="markerGenePeak_clusterCMP/RNAlabel_FIBROvsSMC_peak_diffscore.txt",row.names=F,col.names=T,sep="\t",quote=F)

pv <- markerPlot(seMarker = markersGS_FIBROvSMC, name = "Fibro", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "Volcano")
plotPDF(pv, name = "markersGS_FIBROvSMC_volcano.pdf", ArchRProj = proj_CAD_8, addDOC = FALSE, width = 6, height = 6)
pv <- markerPlot(seMarker = markersPeaks_FIBROvSMC, name = "Fibro", cutOff = "FDR <= 0.01 & Log2FC >= 1", plotAs = "Volcano")
plotPDF(pv, name = "markersPeaks_FIBROvSMC_volcano.pdf", ArchRProj = proj_CAD_8, addDOC = FALSE, width = 6, height = 6)


# SMC related subcluster markerGene, markerPeak
output_diffpeak_clusterCMP <- function(this_data,cellname){
    out_data <- cbind(as.vector(this_data@seqnames),
                      this_data@ranges@start,
                      this_data@ranges@start + this_data@ranges@width,
                      paste0(cellname,"_",seq(length(this_data@seqnames))),
                      this_data$Log2FC,
                      this_data$FDR)
    write.table(out_data,file=paste0("markerGenePeak_clusterCMP/",cellname,"_markerPeak.bed"),quote=F,row.names=F,col.names=F,sep="\t")
}
output_markerGene_clusterCMP <- function(this_data,outname){
    out_data <- cbind(as.vector(this_data$name),
                      this_data$Log2FC,
                      this_data$FDR)
    colnames(out_data) <- c("symbol","log2FC","FDR")
    write.table(out_data,file=paste0("markerGenePeak_clusterCMP/",outname,"_markerGene.txt"),quote=F,row.names=F,col.names=T,sep="\t")
}

markerGenePeak_clusterCMP <- function(fgCluster,bgCluster,outname){
  ############ markerGene
  markersGS_clusterCMP <- getMarkerFeatures(
      ArchRProj = proj_CAD_8, 
      useMatrix = "GeneScoreMatrix", 
      useGroups = c(fgCluster),
      bgdGroups = c(bgCluster),
      groupBy = "Clusters",
      bias = c("TSSEnrichment", "log10(nFrags)"),
      testMethod = "wilcoxon"
  )
  markerListGS_clusterCMP <- getMarkers(markersGS_clusterCMP, cutOff = "FDR <= 0.01 & Log2FC >= 1")
  output_markerGene_clusterCMP(markerListGS_clusterCMP[[1]],outname)

markerGenePeak_clusterCMP("C4","C5","C4vsC5")
markerGenePeak_clusterCMP("C4","C6","C4vsC6")
markerGenePeak_clusterCMP("C4","C7","C4vsC7")
markerGenePeak_clusterCMP("C5","C4","C5vsC4")
markerGenePeak_clusterCMP("C5","C6","C5vsC6")
markerGenePeak_clusterCMP("C5","C7","C5vsC7")
markerGenePeak_clusterCMP("C6","C4","C6vsC4")
markerGenePeak_clusterCMP("C6","C5","C6vsC5")
markerGenePeak_clusterCMP("C6","C7","C6vsC7")
markerGenePeak_clusterCMP("C7","C4","C7vsC4")
markerGenePeak_clusterCMP("C7","C5","C7vsC5")
markerGenePeak_clusterCMP("C7","C6","C7vsC6")


proj_CAD_8$tmpClusters <- rep("C0",length(proj_CAD_8$Clusters))
proj_CAD_8$tmpClusters[which(proj_CAD_8$Clusters == "C4")] <- "C47"
proj_CAD_8$tmpClusters[which(proj_CAD_8$Clusters == "C7")] <- "C47"
proj_CAD_8$tmpClusters[which(proj_CAD_8$Clusters == "C5")] <- "C56"
proj_CAD_8$tmpClusters[which(proj_CAD_8$Clusters == "C6")] <- "C56"

markerGenePeak_clusterCMP2 <- function(fgCluster,bgCluster,outname){

  markersGS_clusterCMP <- getMarkerFeatures(
      ArchRProj = proj_CAD_8, 
      useMatrix = "GeneScoreMatrix", 
      useGroups = c(fgCluster),
      bgdGroups = c(bgCluster),
      groupBy = "tmpClusters",
      bias = c("TSSEnrichment", "log10(nFrags)"),
      testMethod = "wilcoxon"
  )
  markerListGS_clusterCMP <- getMarkers(markersGS_clusterCMP, cutOff = "FDR <= 0.01 & Log2FC >= 1")
  output_markerGene_clusterCMP(markerListGS_clusterCMP[[1]],outname)

  ############ markerPeak
  markersPeaks_clusterCMP <- getMarkerFeatures(
    ArchRProj = proj_CAD_8, 
    useMatrix = "PeakMatrix", 
    useGroups = c(fgCluster),
    bgdGroups = c(bgCluster),
    groupBy = "tmpClusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )
  markerListPeak_clusterCMP <- getMarkers(markersPeaks_clusterCMP, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
  output_diffpeak_clusterCMP(markerListPeak_clusterCMP[[1]],outname)

  ############ motif
  enrichMotifs_clusterCMP <- peakAnnoEnrichment(
      seMarker = markersPeaks_clusterCMP,
      ArchRProj = proj_CAD_8,
      peakAnnotation = "Motif",
      cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )
  write.table(enrichMotifs_clusterCMP@assays$data$mlog10Padj,file=paste0("markerGenePeak_clusterCMP/",outname,"_homerMotif.txt"), row.names=T,col.names=F,sep="\t",quote=F)

}
markerGenePeak_clusterCMP2("C47","C56","C47vsC56")
markerGenePeak_clusterCMP2("C56","C47","C56vsC47")


# disease samples v.s. normal samples
feature_data <- (getCellColData(proj_CAD_8))
Ischemic_samples <- c("scA","scB","scC","scD","scG","scI","scJ","scK","scL","scM","scO","scP","scZ","scAH")
Normal_samples <- c("scQ","scR","scS","scT","scU","scV","scW","scX","scY","scAA","scAB","scAD","scAE","scAF","scAM")

dir.create("markerGenePeak_disease_final")
celltype <- "Endothelial"
diseaseMKG <- function(celltype){
  NDcluster <- rep(0, nrow(feature_data))
  NDcluster[which(feature_data[,"Sample"] %in% Ischemic_samples & feature_data[,"Clusters3"] == celltype)] <- "Ischemic"
  NDcluster[which(feature_data[,"Sample"] %in% Normal_samples & feature_data[,"Clusters3"] == celltype)] <- "Normal"
  NDcluster[which(NDcluster == 0)] <- "Other"
  proj_CAD_8$tmpCluster <- NDcluster
  print(celltype)
  print(table(NDcluster))
  #############marker gene
  markersGS_tmp <- getMarkerFeatures(
    ArchRProj = proj_CAD_8, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "tmpCluster",
    useGroups="Ischemic",
    bgdGroups="Normal",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )
  markerList_tmp <- getMarkers(markersGS_tmp, cutOff = "FDR <= 0.01 & Log2FC >= 0.58",)
  tmpdata <- as.matrix(markerList_tmp$Ischemic)[,c(5,7,8)]
  outdata <- tmpdata[order(as.numeric(tmpdata[,2])),]
  write.table(outdata, file=paste0("markerGenePeak_disease_final/markersGene_disease_",celltype,".txt"),row.names=F,col.names=T,sep="\t",quote=F)

  ############# marker peak
  markersPeaks_tmp <- getMarkerFeatures(
      ArchRProj = proj_CAD_8, 
      useMatrix = "PeakMatrix", 
      groupBy = "tmpCluster",
    useGroups="Ischemic",
    bgdGroups="Normal",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )
  markerList_peak_tmp <- getMarkers(markersPeaks_tmp, cutOff = "FDR <= 0.01 & Log2FC >= 0.58", returnGR = TRUE)
  this_data <- markerList_peak_tmp$Ischemic
  cellname <- celltype
  out_data <- cbind(as.vector(this_data@seqnames),
                    this_data@ranges@start,
                    this_data@ranges@start + this_data@ranges@width,
                    paste0(cellname,"_",seq(length(this_data@seqnames))),
                    this_data$Log2FC,
                    this_data$FDR)
  write.table(out_data,file=paste0("markerGenePeak_disease_final/markerPeak_diseaseFC15_",cellname,".bed"),quote=F,row.names=F,col.names=F,sep="\t")

  markerList_peak_tmp <- getMarkers(markersPeaks_tmp, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
  this_data <- markerList_peak_tmp$Ischemic
  cellname <- celltype
  out_data <- cbind(as.vector(this_data@seqnames),
                    this_data@ranges@start,
                    this_data@ranges@start + this_data@ranges@width,
                    paste0(cellname,"_",seq(length(this_data@seqnames))),
                    this_data$Log2FC,
                    this_data$FDR)
  write.table(out_data,file=paste0("markerGenePeak_disease_final/markerPeak_diseaseFC2_",cellname,".bed"),quote=F,row.names=F,col.names=F,sep="\t")

  markerList_peak_tmp <- getMarkers(markersPeaks_tmp, cutOff = "FDR <= 0.01 ", returnGR = TRUE)
  this_data <- markerList_peak_tmp$Ischemic
  cellname <- celltype
  out_data <- cbind(as.vector(this_data@seqnames),
                    this_data@ranges@start,
                    this_data@ranges@start + this_data@ranges@width,
                    paste0(cellname,"_",seq(length(this_data@seqnames))),
                    this_data$Log2FC,
                    this_data$FDR)
  write.table(out_data,file=paste0("markerGenePeak_disease_final/markerPeak_diseaseALL_",cellname,".bed"),quote=F,row.names=F,col.names=F,sep="\t")

}

diseaseMKG("Endothelial")
diseaseMKG("Fibroblast")
diseaseMKG("Macrophage")
diseaseMKG("Mast")
diseaseMKG("Pericyte")
diseaseMKG("Plasma")
diseaseMKG("SMC")
diseaseMKG("T")


################ generate other plots for the manuscript
DSMCmarker <- c("MYOCD","SRF","TEAD3","TEAD4","ACTA2","MYH11","TAGLN","LMOD1","CNN1","TPM2","MYL9")
MSMCmarker <- c("TCF21","KLF4","FN1","LUM","TNFRSF11B","BGN")
Emarker <- c("KLF2","PECAM1","CLDN5","PLVAP","ACKR1","EGFL7", "NFKB1","NFKB2","VCAM1","SELE")
Tmarker <- c("CD8A","TCF7","RUNX3","TBX21","PRDM1")
Macrophage <- c("CD14","CD36","CD68","CD86","CSF1R","NR1H3","NR1H2","RXRA","RXRB","RXRG","IL1B","CX3CR1")
PericyteMarker <- c("NOTCH3","PDGFRB","RGS5","CSPG4")
Osteochondrogenic <- c("SOX9","RUNX2","BMP2","ALPL")
PotentialSMC <- c("CARMN","PRDM16")
driverSMC <- c("LGALS3","FN1","TNFRSF11B","COL1A1","COL4A1","COL4A2","COL6A3")
NEWmarker <- c("SMAD2","SMAD3","SMAD4","SMAD7","ATF3","JUN","JUND","FOXC1","FOXL1","SUSD2","SUSD5","FHL3","FHL5","PALLD","C3")

thisMKG <- c("MYOCD",driverSMC)

# scatter browser view
p <- plotBrowserTrack(
  ArchRProj = proj_CAD_8,
    plotSummary = c("bulkTrack", "scTrack", "geneTrack"),
    sizes = c(8, 8, 4),
  groupBy = "Clusters3", 
  geneSymbol = DSMCmarker, 
  useMatrix = NULL,
  upstream = 100000, 
  scCellsMax = 100,
  downstream = 100000
)
plotPDF(p,name="scatter_browserTrack_DSMCmarker.pdf", width = 8, height = 8, ArchRProj = proj_CAD_8)

# gene score UMAP

geneset_scatter <- function(thisMKG,outname){
  pImp <- plotEmbedding(
      ArchRProj = proj_CAD_8, 
      colorBy = "GeneScoreMatrix", 
      name = thisMKG, 
      embedding = "UMAP",
      imputeWeights = getImputeWeights(proj_CAD_8)
  )
  
  p2Imp <- lapply(pImp, function(x){
      x + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 6.5) +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
      theme(
          axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(), 
          axis.text.y=element_blank(), 
          axis.ticks.y=element_blank()
      )
  })
  plotPDF(p2Imp, name = paste0("UMAP_geneScoreImpute_MKG_",outname,".pdf"), width = 5, height = 5, ArchRProj = proj_CAD_8, addDOC = FALSE)

}

geneset_scatter(c("ATF3","KLF4","MYOCD","SMAD3","TCF21","TEAD4"),"newUSE")
geneset_scatter(driverSMC,"driverSMC")
geneset_scatter(DSMCmarker,"differentialSMC")
geneset_scatter(MSMCmarker,"modulatedSMC")
geneset_scatter(Emarker,"Endothelial")
geneset_scatter(Tmarker,"Immune")
geneset_scatter(Macrophage,"Macrophage")
geneset_scatter(PericyteMarker,"Pericyte")
geneset_scatter(Osteochondrogenic,"Osteochondrogenic")
geneset_scatter(PotentialSMC,"PotentialSMC")
geneset_scatter(NEWmarker,"NEWmarker")

geneset_scatter_RNAintegration <- function(thisMKG,outname){
  pImp <- plotEmbedding(
      ArchRProj = proj_CAD_8, 
      colorBy = "GeneIntegrationMatrix", 
      name = thisMKG, 
      embedding = "UMAP",
      imputeWeights = getImputeWeights(proj_CAD_8)
  )
  
  p2Imp <- lapply(pImp, function(x){
      x + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 6.5) +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
      theme(
          axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(), 
          axis.text.y=element_blank(), 
          axis.ticks.y=element_blank()
      )
  })
  plotPDF(p2Imp, name = paste0("UMAP_GeneIntegration_MKG_",outname,".pdf"), width = 5, height = 5, ArchRProj = proj_CAD_8, addDOC = FALSE)
}

geneset_scatter_RNAintegration(c("ATF3","KLF4","MYOCD","SMAD3","TCF21","TEAD4"),"newUSE")
geneset_scatter_RNAintegration(driverSMC,"driverSMC")
geneset_scatter_RNAintegration(DSMCmarker,"differentialSMC")
geneset_scatter_RNAintegration(MSMCmarker,"modulatedSMC")
geneset_scatter_RNAintegration(Emarker,"Endothelial")
geneset_scatter_RNAintegration(Tmarker,"Immune")
geneset_scatter_RNAintegration(Macrophage,"Macrophage")
geneset_scatter_RNAintegration(PericyteMarker,"Pericyte")
geneset_scatter_RNAintegration(Osteochondrogenic,"Osteochondrogenic")
geneset_scatter_RNAintegration(PotentialSMC,"PotentialSMC")
geneset_scatter_RNAintegration(NEWmarker,"NEWmarker")


# gene score heatmap
markersGS <- readRDS(file="markerGene_celltype_final/markerGene_celltype_final.rds")
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1")

MKGall <- unique(c(markerList$Endothelial[,"name"],
            markerList$SMC[,"name"],
            markerList$Fibroblast[,"name"],
            markerList$Pericyte[,"name"],
            markerList$Macrophage[,"name"],
            markerList$Mast[,"name"],
            markerList$unknown[,"name"],
            markerList$Plasma[,"name"],
            markerList$T[,"name"]
            ))

celltypeIDX <- c("Endothelial","SMC","Fibroblast","Pericyte","Macrophage","Mast","unknown","Plasma","T")
MKG_GS <- celltype_GS[MKGall,celltypeIDX]
MKG_GI <- celltype_GI[intersect(MKGall,rownames(celltype_GI)),celltypeIDX]

a <- log10(MKG_GI+0.001)

bi_heatmap_fix<-function(data0,usecolor,zmin,zmax,M){
  ColorRamp <- colorRampPalette(usecolor, bias=1)(10000)   #color list
  ColorLevels <- seq(to=zmax,from=zmin, length=10000)   #number sequence
  data0[data0<zmin] <- zmin
  data0[data0>zmax] <- zmax
  ColorRamp_ex <- ColorRamp[round( (min(data0)-zmin)*10000/(zmax-zmin) ) : round( (max(data0)-zmin)*10000/(zmax-zmin) )]
  image(1:ncol(data0), 1:nrow(data0), t(data0), axes=F, col=ColorRamp_ex, xlab="", ylab="",main=M,useRaster=T,cex.main=1.2)
  box()
}

bi_heatmap_fix_heatonly<-function(data0,usecolor,zmin,zmax,M){
  ColorRamp <- colorRampPalette(usecolor, bias=1)(10000)   #color list
  ColorLevels <- seq(to=zmax,from=zmin, length=10000)   #number sequence
  data0[data0<zmin] <- zmin
  data0[data0>zmax] <- zmax
  ColorRamp_ex <- ColorRamp[round( (min(data0)-zmin)*10000/(zmax-zmin) ) : round( (max(data0)-zmin)*10000/(zmax-zmin) )]
  image(1:ncol(data0), 1:nrow(data0), t(data0), axes=F, col=ColorRamp_ex, xlab="", ylab="",main="",useRaster=T,cex.main=1.2)
  box()
}

library(colorRamps)
Z_MKG_GS <- as.matrix(scale(t(MKG_GS)))
Z_MKG_GI <- as.matrix(scale(t(MKG_GI),scale=T))
Z_MKG_GI[is.na(Z_MKG_GI)] <- -2

pdf(file="MKG_GSGI_heatmap_custom.pdf",height=4,width=12)
par(mar=c(2,2,4,6),mfrow=c(1,2))
bi_heatmap_fix(Z_MKG_GS[9:1,],colorRamps::blue2yellow(100),-2,2,"MKG geneScore")
axis(side=4,at=1:9,labels=celltypeIDX[9:1],las=2)
bi_heatmap_fix(Z_MKG_GI[9:1,],colorRamps::blue2yellow(100),-2,2,"MKG GeneIntegration")
axis(side=4,at=1:9,labels=celltypeIDX[9:1],las=2)
image(1:100,1,matrix(data=1:100, nrow=100,ncol=1),col=colorRamps::blue2yellow(100), xlab="",ylab="",cex.axis=1,xaxt="n",yaxt="n",useRaster=T)
axis(side=1,at=c(1,100),labels=c("-2","2"))
dev.off()

celltype_GS <- read.table("geneScore_matrix_celltypeLevel.txt",row.names=1,header=T)
pdf(file="MKG_GS_heatmap_custom.pdf",height=4,width=12)
par(mar=c(2,2,4,6),mfrow=c(1,2))
bi_heatmap_fix(Z_MKG_GS[9:1,],colorRamps::blue2yellow(100),-2,2,"MKG geneScore")
axis(side=4,at=1:9,labels=celltypeIDX[9:1],las=2)
image(1:100,1,matrix(data=1:100, nrow=100,ncol=1),col=colorRamps::blue2yellow(100), xlab="",ylab="",cex.axis=1,xaxt="n",yaxt="n",useRaster=T)
axis(side=1,at=c(1,100),labels=c("-2","2"))
dev.off()



# SMC subcluster trajectory

cluster_trajectory <- function(trajectory, outname){
  tmpProj <- addTrajectory(
    ArchRProj = proj_CAD_8, 
    name = "Trajectory", 
    groupBy = "Clusters",
    trajectory = trajectory, 
    embedding = "UMAP", 
    force = TRUE
)
  p <- plotTrajectory(tmpProj, trajectory = "Trajectory", colorBy = "cellColData", name = "Trajectory")
  plotPDF(p, name = paste0("Trajectory_",outname,"_UMAP_final.pdf"), width = 5, height = 5, ArchRProj = tmpProj, addDOC = FALSE)
  
  trajMM  <- getTrajectory(ArchRProj = tmpProj, name = "Trajectory", useMatrix = "MotifMatrix", log2Norm = FALSE)
  trajGSM <- getTrajectory(ArchRProj = tmpProj, name = "Trajectory", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
  trajGIM <- getTrajectory(ArchRProj = tmpProj, name = "Trajectory", useMatrix = "GeneIntegrationMatrix", log2Norm = FALSE)
  
  p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"))
  p2 <- trajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "horizonExtra"))
  p3 <- plotTrajectoryHeatmap(trajGIM,  pal = paletteContinuous(set = "blueYellow"))
  
  plotPDF(p1, p2, p3, name = paste0("Trajectory_",outname,"_Heatmaps_final.pdf"), ArchRProj = proj_CAD_8, addDOC = FALSE, width = 6, height = 6)

}
cluster_trajectory(c("C6", "C5", "C4"),"cluster654")
cluster_trajectory(c("C6", "C5", "C7"),"cluster657")
cluster_trajectory(c("C6", "C5", "C3"),"cluster653")
cluster_trajectory(c("C5", "C6", "C7"),"cluster567")




