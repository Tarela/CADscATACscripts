#Load packages
library(rasqualTools)
library(tidyverse)
library(data.table)

#Import counts matrix
macrophage_counts <- fread("Macrophage_rasqual_counts_matrix_3")
macrophage_counts <- as_tibble(macrophage_counts)


#Import gc content of each peak
gc_peaks <- fread("all_peaks_gc.bed")


#Add gc content column
macrophage_counts <- cbind(macrophage_counts, gc_peaks$`6_pct_gc`)


#Assign column names to counts matrix
colnames(macrophage_counts) <- c("Geneid", "Chr", "Start", "End", "Strand", "Length",
                                 "UV005", "UVA007", "UVA008", "UVA009", "UVA014", "UVA015",
                                 "UVA016", "UVA017", "UVA019", "UVA020", "UVA023", "UVA024",
                                 "UVA027", "UVA030", "UVA036", "UVA039", "UVA040", "UVA043",
                                 "UVA045", "UVA046", "UVA049", "UVA050", "UVA053", "UVA054",
                                 "UVA057", "UVA058", "UVA061", "UVA062", "UVA064", "UVA066",
                                 "UVA067", "UVA068", "UVA069", "UVA079", "UVA080", "UVA083",
                                 "UVA094", "UVA108", "UVA112", "UVA124", "UVA153", "GC")


#Remove UVA027 and UVA062
macrophage_counts <- select(macrophage_counts, -UVA027, -UVA062)


#Filter counts matrix - remove peaks with low read counts
macrophage_counts <- macrophage_counts %>% mutate(mean = rowMeans(macrophage_counts[,7:45], na.rm=TRUE)) %>% filter(mean >= 5)


#Assign column names to counts matrix
colnames(macrophage_counts) <- c("Geneid", "Chr", "Start", "End", "Strand", "Length",
                                 "UVA005", "UVA007", "UVA008", "UVA009", "UVA014", "UVA015",
                                 "UVA016", "UVA017", "UVA019", "UVA020", "UVA023", "UVA024",
                                 "UVA030", "UVA036", "UVA039", "UVA040", "UVA043", "UVA045", 
                                 "UVA046", "UVA049", "UVA050", "UVA053", "UVA054", "UVA057",
                                 "UVA058", "UVA061", "UVA064", "UVA066", "UVA067", "UVA068",
                                 "UVA069", "UVA079", "UVA080", "UVA083", "UVA094", "UVA108",
                                 "UVA112", "UVA124", "UVA153", "GC", "Mean")


#Prepare gene_data (Geneid is the peak name)
gene_data = dplyr::select(macrophage_counts, Geneid, Chr, Strand, Start, End, GC)
gene_data_offset <- gene_data %>% select(Geneid, GC)
colnames(gene_data_offset) <- c("gene_id", "percentage_gc_content")
print(gene_data_offset)


#Format counts matrix
colnames(macrophage_counts)
peaknames <- macrophage_counts$Geneid
macrophage_counts <- macrophage_counts[, 7:45]
rownames(macrophage_counts) <- peaknames

#Save matrix
saveRasqualMatrices(list(Macrophage = macrophage_counts), "/scratch/amt2ug/rasqualTools", file_suffix = "counts")

#Calculate size factors
size_factors = rasqualCalculateSampleOffsets(macrophage_counts, gene_data_offset, gc_correct = TRUE)
saveRasqualMatrices(list(Macrophage = size_factors), "/scratch/amt2ug/rasqualTools", file_suffix = "size_factors")


#Calculate the number of SNPs overlapping each peak
snp_coords <- read.table("snp_list_2.txt", header = TRUE, sep = "\t")
colnames(snp_coords) <- c("chr", "pos", "snp_id")

gene_metadata <- gene_data[ , 1:5]
colnames(gene_metadata) <- c("gene_id", "chr", "strand", "start", "end")

snp_counts = countSnpsOverlapingPeaks(gene_metadata, snp_coords, cis_window = 1e4)


#Format columns for RASQUAL input
snp_counts <- snp_counts %>% unite("z", chromosome_name, range_start, sep = ":", remove = FALSE) %>% unite("region", z, range_end, sep = "-", remove = FALSE)
snp_counts <- snp_counts %>% mutate(gene_name = gene_id)

#Change order of columns
#Exon start and ends correspond to peak start and end positions
snp_counts_2 <- snp_counts %>% select(gene_id, gene_name, region, cis_snp_count, feature_snp_count, exon_starts, exon_ends)

#Export table
write.table(snp_counts_2, file = "rasqual.Macrophage.input.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)