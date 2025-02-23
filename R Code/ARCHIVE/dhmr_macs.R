library(Rsamtools)
library(GenomicRanges)
library(GenomicAlignments)
library(rtracklayer)
library(ChIPseeker)
library(annotatr)
library(GenomicRanges)
library(dplyr)

# Load the TxDb object for the reference genome (e.g., hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene


# Set the output directory
output_dir <- "~/peritoneal_processing/macs_3_24_2024"

# Create the output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Download and read the blacklist file
blacklist_file <- file.path("~/reference_genomes/wgEncodeDacMapabilityConsensusExcludable.bed")
blacklist <- import(blacklist_file, format = "BED")

excluded_samples <- c("KT026", "KT027", "KT357", "KT13")

# Set the directory containing the BAM files
bam_dir <- "~/peritoneal_processing/trimmed_data_bam_9_18_2023"

# Get the list of BAM files
bam_files <- list.files(path = bam_dir, pattern = "^.*_bowtie2\\.bam$", full.names = TRUE)
bam_files <- bam_files[!sapply(excluded_samples, function(x) any(startsWith(bam_files, x)))]

# Initialize an empty list to store the peak regions
peak_list <- list()

# Process each BAM file
for (bam_file in bam_files) {
  filtered_bam <- bam_file
  
  # Check if the peak file already exists
  peak_file <- file.path(output_dir, paste0(basename(filtered_bam), "_macs2_peaks.narrowPeak"))
  if (file.exists(peak_file)) {
    message(paste("Peak file already exists for", bam_file, ". Skipping processing."))
    peaks <- read.table(peak_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                        col.names = c("chr", "start", "end", "name", "score", "strand",
                                      "signalValue", "pValue", "qValue", "peak"))
    peaks <- makeGRangesFromDataFrame(peaks, keep.extra.columns = TRUE)
    peak_list[[basename(bam_file)]] <- peaks
    next
  }
  
  # Check if the BigWig file already exists
  bigwig_file <- file.path(output_dir, paste0(basename(filtered_bam), ".bw"))
  bedgraph_file <- file.path(output_dir, paste0(basename(filtered_bam), ".bedgraph"))
  sorted_bedgraph_file <- paste0(filtered_bam, "_sorted.bedgraph")
  if (!file.exists(bigwig_file)) {
    # Convert BedGraph to BigWig format
    system(paste("/usr/local/bin/bedtools genomecov -ibam", filtered_bam, "-bg -scale",
                 "$(/tools/bin/samtools view -c", filtered_bam, ") > ", bedgraph_file))
    system(paste("LC_COLLATE=C sort -k1,1 -k2,2n", bedgraph_file, ">", sorted_bedgraph_file))
    system(paste("/tools/bin/bedGraphToBigWig", sorted_bedgraph_file, "~/5hmC-Pathway-Analysis/fastq_processing/Reference_genomes/hg19.chrom.sizes", bigwig_file))
  } else {
    message(paste("BigWig file already exists for", bam_file, ". Skipping conversion."))
  }
  
  macs_output <- file.path(output_dir, paste0(basename(filtered_bam), "_macs2"))
  # Identify potential hMRs using MACS2
  system(paste("/usr/local/bin/macs2 callpeak -t", filtered_bam, "-n", macs_output, "-p 1e-5 -f BAM -g hs"))
  
  # Read the peak regions and store them in the list
  peaks <- read.table(peak_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                      col.names = c("chr", "start", "end", "name", "score", "strand",
                                    "signalValue", "pValue", "qValue", "peak"))
  peaks <- makeGRangesFromDataFrame(peaks, keep.extra.columns = TRUE)
  peak_list[[basename(bam_file)]] <- peaks
  
  # Delete the BedGraph file
  file.remove(bedgraph_file)
  file.remove(sorted_bedgraph_file)
}

# Check if the merged peak file already exists
merged_peak_file <- file.path(output_dir, "merged_peaks_filtered.bed")
if (!file.exists(merged_peak_file)) {
  # Convert the peak_list to a GRangesList object
  peak_list <- GRangesList(peak_list)
  
  # Merge the peak regions across all samples
  merged_peaks <- Reduce(union, peak_list)
  
  # Filter peaks based on occurrence in more than 10 samples and length less than 1000 bp
  peak_counts <- countOverlaps(merged_peaks, peak_list)
  filtered_peaks <- merged_peaks[peak_counts > floor(0.5*length(bam_files))]# & width(merged_peaks) < 1000]
  
  # Exclude peaks that overlap with blacklisted regions by at least 1 base pair
  final_peaks <- subsetByOverlaps(filtered_peaks, blacklist, invert = TRUE, minoverlap = 1)
  
  # Filter out non-canonical chromosomes
  canonical_chromosomes <- paste0("chr", c(1:22))#, "X", "Y"))
  final_peaks <- final_peaks[seqnames(final_peaks) %in% canonical_chromosomes]
  
  # Export the final peak regions to a BED file
  export(final_peaks, merged_peak_file)
} else {
  message("Merged peak file already exists. Skipping merging and filtering.")
  final_peaks <- import(merged_peak_file, format = "BED")
}

# Create a second merged peak file in the desired format
merged_peak_file_saf <- file.path(output_dir, "merged_peaks_filtered_pe5_040324_gene.saf")
merged_peak_file_saf_anno <- file.path(output_dir, "merged_peaks_filtered_pe5_040324_gene_anno.saf")

if(FALSE){
  # OLD ANNOTATION APPROACH
  if (!file.exists(merged_peak_file_saf)) {
    # Annotate the peaks with the closest gene using ChIPseeker
    peak_anno <- annotatePeak(final_peaks, tssRegion = c(-1000, 1000), TxDb = txdb, annoDb = "org.Hs.eg.db")
    
    # Extract the closest gene ID for each peak
    closest_genes <- peak_anno@anno$geneId
    
    # Debug statement: Check the value of closest_genes
    message("closest_genes:")
    print(str(closest_genes))
    
    # Extract the genomic feature for each peak
    genomic_features <- peak_anno@anno$annotation
    
    # Debug statement: Check the value of genomic_features
    message("genomic_features:")
    print(str(genomic_features))
    
    # Map the genomic features to the desired categories
    feature_categories <- c(
      "Promoter" = "promoter",
      "5' UTR" = "5' UTR",
      "3' UTR" = "3' UTR",
      "Exon" = "exon",
      "Intron" = "intron",
      "Downstream" = "downstream",
      "Intergenic" = "intergenic",
      "Distal Intergenic" = "intergenic"
    )
    
    # Assign the mapped categories to the genomic features
    genomic_features <- feature_categories[genomic_features]
    
    # Replace NA values in genomic_features with "other"
    genomic_features[is.na(genomic_features)] <- "other"
    
    # Debug statement: Check the updated value of genomic_features
    message("Updated genomic_features:")
    print(str(genomic_features))
    
    # Debug statement: Check the value of final_peaks
    message("final_peaks:")
    print(str(final_peaks))
    
    # Debug statement: Check the lengths of closest_genes, genomic_features, and final_peaks
    message("Length of closest_genes:", length(closest_genes))
    message("Length of genomic_features:", length(genomic_features))
    message("Length of final_peaks:", length(final_peaks))
    
    # Convert gene IDs to gene symbols
    gene_symbols <- mapIds(org.Hs.eg.db, keys = closest_genes, keytype = "ENTREZID", column = "SYMBOL")
    
    # Replace missing gene symbols with "NA"
    gene_symbols[is.na(gene_symbols)] <- "NA"
    
    # Create a data frame with the desired format
    peak_df_anno <- data.frame(
      GeneID = gene_symbols,
      Chr = seqnames(final_peaks),
      Start = start(final_peaks),
      End = end(final_peaks),
      Strand = ".",
      genomic_feature = genomic_features,
      stringsAsFactors = FALSE
    )
    
    # Create a data frame with the desired format
    peak_df <- data.frame(
      GeneID = peak_df_anno$GeneID,
      Chr = peak_df_anno$Chr,
      Start = peak_df_anno$Start,
      End = peak_df_anno$End,
      Strand = ".",
      stringsAsFactors = FALSE
    )
    
    # Replace missing gene symbols with the concatenation of chr, start, and end
    peak_df_anno$GeneID[is.na(peak_df_anno$GeneID)] <- paste0(
      peak_df_anno$Chr[is.na(peak_df_anno$GeneID)], "_",
      peak_df_anno$Start[is.na(peak_df_anno$GeneID)], "_",
      peak_df_anno$End[is.na(peak_df_anno$GeneID)]
    )
    
    # Replace missing gene symbols with the concatenation of chr, start, and end
    peak_df$GeneID[is.na(peak_df$GeneID)] <- paste0(
      peak_df$Chr[is.na(peak_df$GeneID)], "_",
      peak_df$Start[is.na(peak_df$GeneID)], "_",
      peak_df$End[is.na(peak_df$GeneID)]
    )
    
    # Export the data frame to a SAF file
    write.table(peak_df, file = merged_peak_file_saf, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(peak_df_anno, file = merged_peak_file_saf_anno, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    
  } else {
    message("Merged peak file in SAF format already exists. Skipping formatting.")
  }
}

if(TRUE){
  # NEW ANNOTATION APPROACH
  # Integration with the existing script
  if (!file.exists(merged_peak_file_saf)) {
    source("~/5hmC-Pathway-Analysis/R Code/repeated_glmnet_funcs_v1.R")
    
    # Prepare the data frame with features from final_peaks
    feature_data <- data.frame(feature = paste(seqnames(final_peaks), start(final_peaks), end(final_peaks), sep = "_"))
    
    # Annotate the features using the annotate_features2 function
    annotated_df <- annotate_features2(feature_data)
    
    # Extract the closest gene ID and genomic feature from annotated_df
    closest_genes <- annotated_df$GeneSymbol
    genomic_features <- annotated_df$locationType
    
    # Debug statement: Check the value of closest_genes
    message("closest_genes:")
    print(str(closest_genes))
    
    # Debug statement: Check the value of genomic_features
    message("genomic_features:")
    print(str(genomic_features))
    
    # Debug statement: Check the value of final_peaks
    message("final_peaks:")
    print(str(final_peaks))
    
    # Debug statement: Check the lengths of closest_genes, genomic_features, and final_peaks
    message("Length of closest_genes:", length(closest_genes))
    message("Length of genomic_features:", length(genomic_features))
    message("Length of final_peaks:", length(final_peaks))
    
    # Convert gene IDs to gene symbols (if needed)
    gene_symbols <- mapIds(org.Hs.eg.db, keys = closest_genes, keytype = "ENTREZID", column = "SYMBOL")
    
    # Replace missing gene symbols with "NA"
    gene_symbols[is.na(gene_symbols)] <- "NA"
    
    # Create a data frame with the desired format
    peak_df_anno <- data.frame(
      GeneID = gene_symbols,
      Chr = seqnames(final_peaks),
      Start = start(final_peaks),
      End = end(final_peaks),
      Strand = ".",
      genomic_feature = genomic_features,
      stringsAsFactors = FALSE
    )
    
    # Create a data frame with the desired format
    peak_df <- data.frame(
      GeneID = peak_df_anno$GeneID,
      Chr = peak_df_anno$Chr,
      Start = peak_df_anno$Start,
      End = peak_df_anno$End,
      Strand = ".",
      stringsAsFactors = FALSE
    )
    
    # Replace missing gene symbols with the concatenation of chr, start, and end
    peak_df_anno$GeneID[is.na(peak_df_anno$GeneID)] <- paste0(
      peak_df_anno$Chr[is.na(peak_df_anno$GeneID)], "_",
      peak_df_anno$Start[is.na(peak_df_anno$GeneID)], "_",
      peak_df_anno$End[is.na(peak_df_anno$GeneID)]
    )
    
    # Replace missing gene symbols with the concatenation of chr, start, and end
    peak_df$GeneID[is.na(peak_df$GeneID)] <- paste0(
      peak_df$Chr[is.na(peak_df$GeneID)], "_",
      peak_df$Start[is.na(peak_df$GeneID)], "_",
      peak_df$End[is.na(peak_df$GeneID)]
    )
    
    # Export the data frame to a SAF file
    write.table(peak_df, file = merged_peak_file_saf, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(peak_df_anno, file = merged_peak_file_saf_anno, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    
  } else {
    message("Merged peak file in SAF format already exists. Skipping formatting.")
  }
}
# Generate hMRs for each patient and export to separate BED files
# for (i in seq_along(peak_list)) {
#   patient_id <- names(peak_list)[i]
#   patient_hmr_file <- file.path(output_dir, paste0(patient_id, "_hMRs.bed"))
#   
#   # Check if the patient hMR file already exists
#   if (!file.exists(patient_hmr_file)) {
#     patient_peaks <- peak_list[[i]]
#     overlapping_peaks <- subsetByOverlaps(patient_peaks, final_peaks)
#     
#     # Extract the original columns from the overlapping peaks
#     overlapping_peak_info <- mcols(overlapping_peaks)
#     
#     # Create a data frame with the overlapping peak information
#     overlapping_peak_df <- data.frame(
#       seqnames = seqnames(overlapping_peaks),
#       start = start(overlapping_peaks),
#       end = end(overlapping_peaks),
#       overlapping_peak_info,
#       stringsAsFactors = FALSE
#     )
#     
#     # Write the overlapping peak information to the patient hMR file
#     write.table(
#       overlapping_peak_df,
#       file = patient_hmr_file,
#       sep = "\t",
#       row.names = FALSE,
#       col.names = FALSE,
#       quote = FALSE
#     )
#   } else {
#     message(paste("Patient hMR file already exists for", patient_id, ". Skipping intersection."))
#   }
# }