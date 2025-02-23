library(Rsamtools)
library(GenomicRanges)
library(GenomicAlignments)
library(rtracklayer)
library(ChIPseeker)
library(annotatr)
library(GenomicRanges)
library(dplyr)
library(parallel)
library(DiffBind)

# For raw data processing:
library(tibble)
library(DESeq2)
library(sva)
library(ggplot2)

#BiocManager::install("consensusSeekeR")
#library(consensusSeekeR)
# Load the TxDb object for the reference genome (e.g., hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# print version number
packageVersion("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(org.Hs.eg.db)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

### Set configurations ----
# Set the output directory
comparison_set <- "PM" # "AM" for autoimmune, "PM" for peritoneal metastasis
if(comparison_set == "AM"){
  output_dir <- "/home/turagapm/peritoneal_processing/macs_6_30_2024_p_1e5/AMpos_AMneg_excl"
  selected_conditions <- c("AMpos","AMneg")
  saf_fn <- "merged_peaks_10perc_pe5_072024_gene.saf" # for readcounts
  ## For peaks
  raw_counts_file_1 <- "~/peritoneal_processing/macs_6_30_2024_p_1e5/AMpos_AMneg_excl/ReadsCount_10percent_macs3_pe5_072024.txt" 
  file_version <- "90per15_10402024_hmr_nocombat_ampos_amneg_noX"
  ## for genes
  # raw_counts_file_1 <- "~/peritoneal_processing/macs_6_30_2024_p_1e5/AMpos_AMneg_excl/ReadsCount_10percent_macs3_pe5_092824_zhang_genes.txt" 
  # file_version <- "0_gene_bodies_10302024_hmr_combat_ampos_amneg_noX"
  combat_seq_normalize = TRUE
  ninety_count_thresh <- 15
} else if(comparison_set == "PM"){
  output_dir <- "/home/turagapm/peritoneal_processing/macs_6_30_2024_p_1e5/CRC_ADE_HEALTHY"
  selected_conditions <- c("PM_positive", "PM_negative", "HEALTHY")
  # counts file expected to be in featureCounts default export format
  raw_counts_file_1 <- "~/peritoneal_processing/macs_6_30_2024_p_1e5/CRC_ADE_HEALTHY/ReadsCount_10percent_macs3_pe5_092824_peaks.txt" 
  # ReadsCount_10percent_macs3_pe5_092824_peaks.txt # for all 5hmC genomic regions (>200k features)
  # ReadsCount_10percent_macs3_pe5_092824_zhang_genes.txt" # for gene bodies (~20k features)
  saf_fn <- "merged_peaks_10perc_pe5_10162024.saf" # for filtered consensus 5hmC consensus regions (>71k features)
  file_version <- "10per_11232024_hmr_combat_healthy_noX" # text to append to file name after raw data processing
  combat_seq_normalize = TRUE
  ninety_count_thresh <- 15 # 15 for consensus region filtering
} else{
  stop("Invalid comparison_set. Must be 'AM' or 'PM'.")
}

# excluded chromosomes
excluded_chr <- c("chrX","chrY","chrM")
# Exclude any files with < 10M unique reads in any single paired end fastq file
excluded_samples <- c("KT027","KT026","KT211","KT094","KT029","KT207","KT331","KT035","KT188","KT342")
# Set metadata fileall_comp
metadata_file <- "~/5hmC-Pathway-Analysis/Raw Input/Working Inputs/all_comparisons_11_23_2024.csv" 
# Set the directory containing the BAM files
bam_dir <- "~/peritoneal_processing/trimmed_data_bam_9_18_2023"
#7_2_2024.csv" for v1 analysos or 9_03_2024 for v2 analysis
# Download and read the blacklist file
blacklist_file <- file.path("~/reference_genomes/wgEncodeDacMapabilityConsensusExcludable.bed")
blacklist <- import(blacklist_file, format = "BED")
### End of: Set configurations

### Load up files and create directories ----
# Create the output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# load metadata file
metadata <- read.csv(metadata_file, stringsAsFactors = TRUE)

# Keep only rows where condition in selected conditons
metadata <- metadata[metadata$condition %in% selected_conditions,]
# rename "X" column to "ID"
colnames(metadata)[colnames(metadata) == "X"] <- "ID"

# drop excluded samples from metadata where ID not in excluded_samples
metadata <- metadata[!metadata$ID %in% excluded_samples,]
included_samples <- metadata$ID[!metadata$ID %in% excluded_samples]


# Get the list of BAM files
bam_files <- list.files(path = bam_dir, pattern = "^.*_bowtie2\\.bam$", full.names = TRUE)
#bam_files <- bam_files[!sapply(excluded_samples, function(x) any(startsWith(bam_files, x)))]
# Drop the samples that are do not contain any of the included_samples
# Filter BAM files to include only those related to the included_samples
bam_files <- bam_files[sapply(bam_files, function(file) {
  any(sapply(included_samples, function(sample) grepl(sample, basename(file))))
})]

# Add "bamReads" column to metadata which contains the path to the BAM file for each sample ID
metadata$bamReads <- NA
for (i in seq_len(nrow(metadata))) {
  sample_id <- metadata$ID[i]
  bam_file <- bam_files[grep(sample_id, bam_files)]
  if (length(bam_file) > 0) {
    metadata$bamReads[i] <- bam_file
  }
}
### End of: Load up files and create directories

### Create metadata figure table on clinical characteristics ----
# Load necessary packages
library(dplyr)
library(tableone)
library(magrittr)

if(TRUE){
  # Ensure the output directory exists; if not, create it
  if(!dir.exists(output_dir)){
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Read in your data
  #meta <- read.csv("~/5hmC-Pathway-Analysis/Raw Input/Working Inputs/clinical_metadata_table_fig_11_23_24.csv", stringsAsFactors = FALSE)
  meta <- read.csv("~/5hmC-Pathway-Analysis/Raw Input/Working Inputs/clinical_metadata_table_fig_10_26_24.csv", stringsAsFactors = FALSE)
  # rename "ID" column to "X"
  colnames(meta)[colnames(meta) == "ID"] <- "X"
  
  # Clean and verify the 'CEA' variable
  # Identify all unique values in 'CEA' to check for non-numeric entries
  print("Unique values in 'CEA':")
  print(unique(meta$CEA))
  
  # Replace any non-numeric values with NA (e.g., "Unknown")
  meta$CEA_value <- suppressWarnings(as.numeric(meta$CEA))
  
  # Create 'CEA_status' variable indicating if CEA is known or unknown
  meta <- meta %>%
    mutate(
      CEA_status = ifelse(is.na(CEA_value), "Unknown", "Known"),
      CEA_status = factor(CEA_status, levels = c("Known", "Unknown"))
    )
  
  # Apply your specified mappings to the variables
  meta <- meta %>%
    mutate(
      # Define the groups
      Group = case_when(
        condition == "HEALTHY" ~ "Healthy controls (HC, n=73)",
        condition == "PM_positive" ~ "Cancerous, Secondary Peritoneal Carcinomatosis (PM, n=71)",
        condition == "PM_negative" ~ "Cancerous, No Secondary Peritoneal Carcinomatosis (No PM, n=41)"
      ),
      Group = factor(Group, levels = c("Healthy controls (HC, n=73)", 
                                       "Cancerous, Secondary Peritoneal Carcinomatosis (PM, n=71)", 
                                       "Cancerous, No Secondary Peritoneal Carcinomatosis (No PM, n=41)")),
      
      # Map 'sex' variable
      sex = factor(sex, levels = c("female", "male"), labels = c("Female", "Male")),
      
      # Map 'race' variable according to your specified labels
      race = factor(race,
                    levels = c("White", "Black_African_American", "Asian_Mideast_Indian", "Unknown_Other"),
                    labels = c("White", "Black, African American", "Asian, Mideast, Indian", "Unknown, Other")),
      
      # Map 'batch' variable
      batch = factor(batch, levels = c("batch_1", "batch_2"), labels = c("Batch 1", "Batch 2")),
      
      # Map 'peritoneal_mets' variable
      peritoneal_mets = factor(peritoneal_mets, levels = c("pm_present", "pm_absent"), 
                               labels = c("Peritoneal metastasis", "No peritoneal metastasis")),
      
      # Map 'other_mets' variable
      other_mets = factor(other_mets, levels = c("other_mets_present", "other_mets_absent"), 
                          labels = c("Other metastases (e.g. liver, lung...)", "No other metastases")),
      
      # Map 'primary_site' variable
      primary_site = factor(primary_site, levels = c("CRC", "ADE"), 
                            labels = c("Colorectal cancer", "Appendix adenocarcinoma")),
      
      # Map 'primary_site' variable
      primary_present = factor(primary_present, levels = c("primary_present", "primary_absent"), 
                            labels = c("Primary tumor present", "Primary tumor absent")),
      
      # Map 'treatment_status' variable
      treatment_status = factor(treatment_status,
                                levels = c("None", "Surgery only", "Chemo only", 
                                           "Chemo + radiation", "Chemo + Surgery", 
                                           "Chemo + Radiation + Surgery"),
                                labels = c("None", "Surgery only", "Chemo only", 
                                           "Chemo and radiation", "Chemo and surgery", 
                                           "Chemo, radiation and surgery")),
      
      # Map 'Path Grade' variable
      `Path Grade` = factor(Path.Grade, levels = c("1", "2", "3", "Unknown"),
                            labels = c("G1 (low)", "G2 (intermediate)", "G3 (high)", "Unknown"))
    )
  
  # Convert 'age' to numeric if it's not already
  meta$age <- as.numeric(meta$age)
  
  # Subset data excluding healthy controls where appropriate
  meta_cancer <- meta %>% filter(condition != "HEALTHY")
  
  # Reset factor levels in 'Group' after subsetting
  meta_cancer$Group <- droplevels(meta_cancer$Group)
  
  # Define variables for the overall table
  vars_overall <- c("sex", "age", "race", "batch")  # 'batch' is already included here
  catVars_overall <- c("sex", "race", "batch")
  
  # Define variables for the cancer-specific table (including 'batch')
  vars_cancer <- c("peritoneal_mets", "other_mets", "primary_site", 
                   "treatment_status", "primary_present","Path Grade", "CEA_status", "CEA_value", "batch")
  catVars_cancer <- c("peritoneal_mets", "other_mets", "primary_site", 
                      "treatment_status", "primary_present","Path Grade", "CEA_status", "batch")
  
  # Create the table including healthy controls
  table_overall <- CreateTableOne(vars = vars_overall, data = meta, factorVars = catVars_overall, 
                                  strata = "Group", includeNA = FALSE, test = TRUE)
  
  # Create the cancer-specific table including 'batch'
  table_cancer <- CreateTableOne(vars = vars_cancer, data = meta_cancer, 
                                 factorVars = catVars_cancer, strata = "Group", 
                                 includeNA = TRUE, test = TRUE)
  
  # Optional: Print the tables with p-values
  print(table_overall, showAllLevels = TRUE, quote = FALSE, noSpaces = TRUE)
  print(table_cancer, showAllLevels = TRUE, quote = FALSE, noSpaces = TRUE)
  
  ### Corrected Code Starts Here ###
  
  # Convert TableOne objects to data frames for CSV export
  # For table_overall
  table_overall_print <- print(table_overall, showAllLevels = TRUE, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
  table_overall_df <- as.data.frame(table_overall_print)
  table_overall_df <- rownames_to_column(table_overall_df, var = "Variable")
  
  # For table_cancer
  table_cancer_print <- print(table_cancer, showAllLevels = TRUE, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
  table_cancer_df <- as.data.frame(table_cancer_print)
  table_cancer_df <- rownames_to_column(table_cancer_df, var = "Variable")
  
  # Define the mapping of variables to test types
  # This mapping is based on variable characteristics:
  # - Binary categorical: Fisher's Exact Test
  # - Categorical with >2 levels: Chi-squared Test
  # - Continuous: Wilcoxon Rank Sum Test for cancer-specific (2 groups), Kruskal-Wallis Test for overall (3 groups)
  
  # Create mapping for overall table
  test_mapping_overall <- data.frame(
    Variable = vars_overall,
    Test = sapply(vars_overall, function(var){
      if(var %in% catVars_overall){
        # Categorical variables
        if(length(levels(meta[[var]])) == 2){
          return("Fisher's Exact Test")
        } else {
          return("Chi-squared Test")
        }
      } else {
        # Continuous variables (3 groups)
        return("Kruskal-Wallis Test")
      }
    }),
    stringsAsFactors = FALSE
  )
  
  # Create mapping for cancer-specific table
  test_mapping_cancer <- data.frame(
    Variable = vars_cancer,
    Test = sapply(vars_cancer, function(var){
      if(var %in% catVars_cancer){
        # Categorical variables
        if(length(levels(meta_cancer[[var]])) == 2){
          return("Fisher's Exact Test")
        } else {
          return("Chi-squared Test")
        }
      } else {
        # Continuous variables (2 groups)
        return("Wilcoxon Rank Sum Test")
      }
    }),
    stringsAsFactors = FALSE
  )
  
  # Combine both mappings into one dataframe with a 'Table' column
  test_types_df <- bind_rows(
    test_mapping_overall %>% mutate(Table = "Overall"),
    test_mapping_cancer %>% mutate(Table = "Cancer_Specific")
  )
  
  # Save the test types mapping to a separate CSV
  write.csv(test_types_df, file.path(output_dir, "clinical_metadata_test_types_v2.csv"), row.names = FALSE)
  
  # Save the TableOne tables to CSV in output_dir
  write.csv(table_overall_df, file.path(output_dir, "clinical_metadata_table_overall_v2.csv"), row.names = FALSE)
  write.csv(table_cancer_df, file.path(output_dir, "clinical_metadata_table_cancer_specific_v2.csv"), row.names = FALSE)
  
  ### Corrected Code Ends Here ###
  
  # Optional: Print the test types dataframe for verification
  print("Test Types Mapping:")
  print(test_types_df)
}

### Run MACS3 on all bam files in bam_dir ----
process_bam_file <- function(bam_file, output_dir) {
  filtered_bam <- bam_file
  peak_list <- list()
  
  # Check if the peak file already exists
  peak_file <- file.path(output_dir, paste0(basename(filtered_bam), "_macs3_peaks.xls"))
  if (file.exists(peak_file)) {
    message(paste("Peak file already exists for", bam_file, ". Skipping processing."))
    # peaks <- read.table(peak_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
    #                     col.names = c("chr", "start", "end", "name", "score", "strand",
    #                                   "signalValue", "pValue", "qValue", "peak"))
    # peaks <- makeGRangesFromDataFrame(peaks, keep.extra.columns = TRUE)
    # peaks <- readNarrowPeakFile(peak_file, extractRegions = TRUE, extractPeaks = TRUE)
  } else{
    macs_output <- file.path(output_dir, paste0(basename(filtered_bam), "_macs3"))
    # Identify potential hMRs using MACS3
    system(paste(
      "/usr/local/bin/macs3 callpeak",
      "-t", filtered_bam,          # Input BAM file (treatment)
      "-n", macs_output,           # Output prefix for peak files
      "--nolambda",                # Disables local lambda model (appropriate without control)
      "-f BAMPE",                  # Specifies paired-end BAM format
      "-p 1e-5",                   # Sets stringent p-value cutoff
      "-g hs"                      # Specifies human genome size
    ))    
    # Read the peak regions and store them in the list
    # peaks <- read.table(peak_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
    #                     col.names = c("chr", "start", "end", "name", "score", "strand",
    #                                   "signalValue", "pValue", "qValue", "peak"))
    # peaks <- makeGRangesFromDataFrame(peaks, keep.extra.columns = TRUE)
    # Add name to peaks$peak and peaks$narrowPeak based on basename and length of peaks
    # peaks <- readNarrowPeakFile(peak_file, extractRegions = TRUE, extractPeaks = TRUE)
  }
  
  # # Replace any occurrence of "." with "_" in the basename
  # base_name <- gsub("\\.", "_", basename(filtered_bam))
  # 
  # # Add name to peaks$peak and peaks$narrowPeak based on basename and length of peaks
  # peaks[["peak"]]$name <- paste0(base_name, "_peak_", seq_along(peaks[["peak"]]))
  # names(peaks[["peak"]]) <- rep(base_name, length(peaks[["peak"]]))
  # peaks[["narrowPeak"]]$name <- paste0(base_name, "_peak_", seq_along(peaks[["narrowPeak"]]))
  # names(peaks[["narrowPeak"]]) <- rep(base_name, length(peaks[["narrowPeak"]]))
  
  # return(peaks)
  return(peak_file)
}

# Number of cores to use
num_cores <- 50

# Use mclapply to run in parallel
final_peaks_list <- mclapply(bam_files, process_bam_file, output_dir = output_dir, mc.cores = num_cores)
# flatten list
final_peaks_list <- unlist(final_peaks_list, recursive = FALSE)

### Find consensus peaks using DiffBind ----
# Add final_peaks_list to the metadata as the "Peaks" column by finding the one that contains the row's "ID"
metadata$Peaks <- NA
for (i in seq_len(nrow(metadata))) {
  sample_id <- metadata$ID[i]
  peaks <- final_peaks_list[grep(sample_id, final_peaks_list)]
  if (length(peaks) > 0) {
    metadata$Peaks[i] <- peaks
  }
}

# Add a column called "PeakCaller" that contains "macs" the length of the dataframe
metadata$PeakCaller <- rep("macs", nrow(metadata))

# Read in sample sheet:
healthy_diseased <- dba(
  sampleSheet = metadata,
  #minOverlap = 0.1,
  bRemoveM = TRUE,
  config = data.frame(
    AnalysisMethod = "DBA_DESEQ2",
    DataType = "DBA_DATA_FRAME",
    RunParallel = TRUE,
    #fragmentSize = 180,
    ReportInit = "Healthy_Diseased",
    bUsePval = TRUE,
    doBlacklist = DBA_BLACKLIST_HG19
  )
)

# 1. All peaks
all_peaks <- dba.peakset(healthy_diseased, bRetrieve=TRUE)
all_peaks_df <- as.data.frame(all_peaks)

#overlap_rate <- dba.overlap(healthy_diseased, mode=DBA_OLAP_RATE)
# plot overlap_rate
#plot(overlap_rate)

#consensus_peaks_05 <- dba.peakset(healthy_diseased, 1:185 , consensus=FALSE, minOverlap=0.05,bRetrieve=TRUE)
#consensus_peaks_05_df <- as.data.frame(consensus_peaks_05)

#determine number of samples in healthy diseased
consensus_peaks_10 <- dba.peakset(healthy_diseased, 1:nrow(metadata) , consensus=FALSE, minOverlap=0.1,bRetrieve=TRUE)
consensus_peaks_10_df <- as.data.frame(consensus_peaks_10)

# Exclude unwanted chromosomes
# final_peaks <- seqnames, start and end columns from consensus_peaks_05_df
final_peaks <- consensus_peaks_10_df[, c("seqnames", "start", "end")]
# Drop rows where seqnames do not match chr1, chr2, chr3, etc... chr22 exactly
#final_peaks <- final_peaks[grepl("chr[1-9]|chr1[0-9]|chr2[0-2]", final_peaks$seqnames),]
# drop rows where seqnames contain excluded_chr
final_peaks <- final_peaks[!grepl(paste(excluded_chr, collapse = "|"), final_peaks$seqnames),]
# drop rows where seqnames contain "random"
final_peaks <- final_peaks[!grepl("random", final_peaks$seqnames),]

### End of: Find consensus peaks using DiffBind

### Plot peaks by width to determine optimal width cutoff -------------------
# # Create dataframe with percent of all peaks with width > 1000, 2000, 3000, 4000, 5000, 6000, 7000,8000,9000,10000
# peak_width_counts <- data.frame(width = c(1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000),
#                                 count = sapply(c(1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000), function(width) {
#                                   sum(consensus_peaks_10_df$width > width)
#                                 }))
# # add percent column
# peak_width_counts$percent <- peak_width_counts$count / nrow(consensus_peaks_10_df) * 100
# 
# # plot peak_width_counts
# plot(peak_width_counts$width, peak_width_counts$count, type = "b", xlab = "Peak Width", ylab = "Number of Peaks",
#      main = "Number of Peaks with Width Greater Than Threshold")
# 
# # Plot distribution of consensus_peaks_10_df$width as a boxplot
# boxplot(all_peaks_df$width, main = "Distribution of Peak Widths in 10% Consensus Peaks", ylab = "Width",
#         #set y axis range to 0-5000
#         ylim = c(0, 5000)
#         )
#
### end of peak width plot 


# # 3. Peaks occurring in at least 10% of samples with summits = 500
# First, re-count with the summit parameter
#healthy_diseased_summit <- dba.count(healthy_diseased, 1:109, minOverlap = 0.1, summits=500,bRemoveDuplicates=TRUE,bParallel=FALSE)
# Step 2: Extract Peak Count Data
#peak_counts <- dba.peakset(healthy_diseased_summit, bRetrieve = TRUE)

# Step 3: Convert to Data Frame
#consensus_peaks_summit_df <- as.data.frame(peak_counts)

# Print information about each dataframe
#print(paste("Number of all peaks:", nrow(all_peaks_df)))
#print(paste("Number of peaks in 5% consensus:", nrow(consensus_peaks_05_df)))
#print(paste("Number of peaks in summit-based 10% consensus:", nrow(consensus_peaks_summit_df)))

### Generate annotated .saf files ----
# Create a second merged peak file in the desired format
merged_peak_file_saf <- file.path(output_dir, saf_fn) # 072224_gene.saf for v1 analysis # 090324_gene.saf for v2 analysis
merged_peak_file_saf_anno <- file.path(output_dir, saf_fn) # 072224_gene.saf for v1 analysis # 090324_gene.saf for v2 analysis

if(TRUE){
  # NEW ANNOTATION APPROACH
  # Integration with the existing script
  if (!file.exists(merged_peak_file_saf)) {
    source("~/5hmC-Pathway-Analysis/R Code/repeated_glmnet_funcs_v1.R")
    
    # Prepare the data frame with features from final_peaks
    feature_data <- data.frame(feature = paste(final_peaks$seqnames, final_peaks$start, final_peaks$end, sep = "_"))
    #print row 96799 of feature_data
    #print(feature_data[96799,])
    
    # Annotate the features using the annotate_features2 function
    annotated_df <- annotate_features2(feature_data)
    
    #print row 96799 of annotated_df
    #print(annotated_df[96799,])
    # Extract the closest gene ID and genomic feature from annotated_df
    closest_genes <- annotated_df$annot.symbol
    genomic_features <- annotated_df$locationType
    
    # Debug statement: Check the value of closest_genes
    message("closest_genes:")
    print(str(closest_genes))
    # print number of NA
    message("Number of NA in closest_genes:", sum(is.na(closest_genes)))
    
    
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
    
    # if any closest_genes are NA:
    if(sum(is.na(closest_genes)) > 0){
      print("There are missing gene symbols. Replacing with chr_start_end.")
      # Convert gene IDs to gene symbols (if needed)
      gene_symbols <- mapIds(org.Hs.eg.db, keys = closest_genes, keytype = "ENTREZID", column = "SYMBOL")
      
      # Replace missing gene symbols with "NA"
      gene_symbols[is.na(gene_symbols)] <- "NA"
    } else {
      print("No missing gene symbols.")
      gene_symbols <- closest_genes
    }
    
    # Create a data frame with the desired format
    peak_df_anno <- data.frame(
      GeneID = gene_symbols,
      Chr = annotated_df$chr,
      Start = annotated_df$start,
      End = annotated_df$end,
      Strand = rep(".", nrow(annotated_df)),
      genomic_feature = genomic_features,
      stringsAsFactors = FALSE
    )
    
    # Create a data frame with the desired format
    peak_df <- data.frame(
      GeneID = paste0(
        peak_df_anno$Chr, "_",
        peak_df_anno$Start, "_",
        peak_df_anno$End),
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
    
    #print row 96799 of peak_df
    #print(peak_df[96799,])
    
    # # Replace missing gene symbols with the concatenation of chr, start, and end
    # peak_df$GeneID[is.na(peak_df$GeneID)] <- paste0(
    #   peak_df$Chr[is.na(peak_df$GeneID)], "_",
    #   peak_df$Start[is.na(peak_df$GeneID)], "_",
    #   peak_df$End[is.na(peak_df$GeneID)]
    # )
    
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
### End of: Generate annotated .saf files


### Run feature counts on specific .saf file with consensus peaks ----
# Set paths
output_count <- file.path(output_dir, "ReadsCount_10percent_macs3_pe5_092824_peaks.txt") 
# hmr: ReadsCount_10percent_macs3_pe5_092824_peaks.txt
# genebodies: ReadsCount_10percent_macs3_pe5_092824_zhang_genes.txt
# 092824 for v3 analysis # 090324 for v2 analysis # 072224 for v1 analysis
#saf_file <- file.path(output_dir, "merged_peaks_10perc_pe5_092824.saf") # 072224_gene.saf for v1 analysis # 090324_gene.saf for v2 analysis
saf_file <- "~/reference_genomes/hg19_genebody_fromZhangLab.saf"
# hmr: merged_peaks_10perc_pe5_092824.saf
# genebodies: hg19_genebody_fromZhangLab.saf
bam_dir <- "/home/turagapm/peritoneal_processing/trimmed_data_bam_9_18_2023" 
log_file <- file.path(output_dir, "featureCounts.log")

# Define the command
command <- paste(
  "/tools/bin/featureCounts -T 64",
  paste0("-a ", saf_file), #
  "-F SAF -p --countReadPairs -Q 20 -C --ignoreDup",
  paste0("-o ", output_count),
  file.path(bam_dir, "*_bowtie2.bam"),
  paste0("2>> ", log_file)
)

# Check if the output file exists
if (!file.exists(output_count)) {
  # Execute the command
  system(command)
} else {
  message("Output file already exists. Skipping featureCounts.")
}

### End of: Run feature counts on specific .saf file with consensus peaks

### rawcounts_processing.R with additional filtering: configuration ----
# Purpose of this script is process raw featureCounts file
# and reduce it to only the required sample in the correct order.

### CONFIGURATION
# set working directory
setwd("~/5hmC-Pathway-Analysis/")

# settings
counts1_sample_pathname <- "../../peritoneal_processing/trimmed_data_bam_9_18_2023/"
counts1_sample_extension <- "_bowtie2.bam"

### Load Sample File
sample_file <- metadata_file
selected_conds <- selected_conditions
sample_data <- read.csv(
  sample_file,
  strip.white = TRUE,
)

### Select relevant samples
# Keep only rows in sample_data where condition is in selected_conditions
sample_data <- sample_data[sample_data$condition %in% selected_conditions, ]

# set X to row names
rownames(sample_data) <- sample_data$X
# drop column called X
sample_data <- sample_data[!names(sample_data) %in% "X"]

# Drop rows where row names is in excluded samples
sample_data <- sample_data[!rownames(sample_data) %in% excluded_samples,]
# drop rows where X is in excluded samples
#sample_data <- sample_data[!sample_data$X %in% excluded_samples, ]

# Read in data
raw_counts_data_1 <- read.table(
  raw_counts_file_1,
  sep = "\t",
  header = TRUE,
  skip = 1,
  check.names = F
)

# drop rows where chromosome is in excluded_chr
raw_counts_data_1 <- raw_counts_data_1[!grepl(paste(excluded_chr, collapse = "|"), raw_counts_data_1$Chr),]

### rawcounts_processing.R: Plot read counts statistics on length ----
# Plot number of rows by length with percentages and cutoff at 10,000
lengths <- raw_counts_data_1$Length

# Create a data frame for lengths
lengths_df <- data.frame(Length = lengths)

# Cap lengths at 5000 and create a new column
lengths_df$Length_Capped <- ifelse(lengths_df$Length > 5000, 5001, lengths_df$Length)

# Define breaks from 0 to 5000 with 500 as bin width
breaks <- seq(0, 5000, by=500)

# Include an extra break for lengths greater than 5000
breaks <- c(breaks, Inf)

# Create length groups
lengths_df$Length_Group <- cut(lengths_df$Length_Capped, breaks=breaks, include.lowest=TRUE, right=FALSE)

# Replace the last bin label with '>=5000'
levels(lengths_df$Length_Group)[length(levels(lengths_df$Length_Group))] <- ">=5000"

# Calculate percentages
lengths_pct <- lengths_df %>%
  group_by(Length_Group) %>%
  summarise(Count = n()) %>%
  mutate(Percentage = (Count / sum(Count)) * 100)

# Convert Length_Group levels to more readable format, up to 5000
lengths_pct$Length_Group <- factor(lengths_pct$Length_Group, 
                                   levels = levels(lengths_pct$Length_Group), 
                                   labels = c("[0,500)", "[500,1000)", "[1000,1500)", 
                                              "[1500,2000)", "[2000,2500)", "[2500,3000)", 
                                              "[3000,3500)", "[3500,4000)", "[4000,4500)", 
                                              "[4500,5000)", ">=5000"))

# Calculate the total number of rows
total_rows <- sum(lengths_pct$Count)

# Plot the histogram with n= in the title and % in the y-axis labels
length_plot <- ggplot(lengths_pct, aes(x=Length_Group, y=Percentage)) +
  geom_bar(stat='identity', fill="grey") +
  geom_text(aes(label=sprintf("%.1f%%", Percentage)), 
            vjust=-0.5, size=3.5) +  # Adding percentage labels on top of bars
  xlab('Length') +
  ylab('Percentage') +
  ggtitle(paste('Histogram of Lengths (n=', total_rows, ')', sep = '')) +  # Adding total n in title
  theme_minimal() +  # Remove grey background
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +  # Add % to y-axis labels
  geom_vline(xintercept=10.5, linetype="dashed", color="#1F78B4", size=1) +  # Vertical dashed line before the last bar
  annotate("text", x=10, y=max(lengths_pct$Percentage) + 2, label="length <= 5,000bp", color="#1F78B4", angle=90, vjust=1,hjust=1)  # Label for the vertical line

### rawcounts_processing.R: Filter out rows with length > 5000 ----
# if zhang_genes not in raw_counts_file_1, then
if(!grepl("zhang_genes", raw_counts_file_1)){
  ## drop length and format
  # Drop rows where "Length" is > 5000
  raw_counts_data_1 <- raw_counts_data_1[raw_counts_data_1$Length <= 5000, ]
} else{
  print("Keeping all features...")
}
# print outline, first couple row names and col names
print(dim(raw_counts_data_1))

# Drop preamble
colnames(raw_counts_data_1) <- gsub("/home/turagapm/peritoneal_processing/trimmed_data_bam_9_18_2023/", "", colnames(raw_counts_data_1))

### PROCESS FILE TO DEFINED OUTPUT
counts_trimmed <- subset(raw_counts_data_1,select=-c(Chr,Start,End,Strand,Length)) # was raw_counts
counts_final <- subset(raw_counts_data_1,select=-c(Chr,Start,End,Strand,Length)) # was raw_counts

# Simplify sample names in count file to remove pathname and extensions
colnames(counts_final) <- gsub(counts1_sample_extension, "", colnames(counts_trimmed))
colnames(counts_final) <- gsub(counts1_sample_pathname, "", colnames(counts_final))

# Eliminate unnecessary samples in count file and sample data file
sample_names <- rownames(sample_data)#dplyr::select(sample_data, c(1)))
counts_final_sample <- counts_final[c('Geneid', sample_names)]
#counts_final_sample <- counts_final_sample[, !names(counts_final_sample) %in% excluded_samples]

# print outline, first couple row names and col names
print(dim(counts_final_sample))
# PRINT IF ANY ROWS HAVE na IN GENEID
print(paste("Number of rows with NA in Geneid:", sum(is.na(counts_final_sample$Geneid))))

### rawcounts_processing.R: PLOT FREQUENCY OF COUNTS ----

# Update read thresholds to have a step of 5
read_thresholds <- seq(0, 50, by = 5)
num_rows <- numeric(length(read_thresholds))

total_rows <- nrow(counts_final_sample)

for (i in seq_along(read_thresholds)) {
  threshold <- read_thresholds[i]
  keep <- rowSums(counts_final_sample >= threshold) >= 0.9 * ncol(counts_final_sample)
  num_rows[i] <- sum(keep)
}

# Calculate percentages
percent_rows <- (num_rows / total_rows) * 100

# Create a data frame for plotting
plot_data <- data.frame(ReadThreshold = read_thresholds, Percentage = percent_rows)

# Plot using ggplot2 with vertical line at x=20 and label for 90% samples with >= 20 reads
library(ggplot2)

# Calculate the dynamic x-intercept based on ninety_count_thresh
# Ensure that ReadThreshold is numeric for proper sorting and comparison
read_thresholds <- sort(unique(plot_data$ReadThreshold))

# Check if ninety_count_threshold exists in ReadThreshold
if (!(ninety_count_thresh %in% read_thresholds)) {
  stop("ninety_count_thresh is not present in ReadThreshold values.")
}

# Find the position of ninety_count_thresh in the sorted ReadThresholds
position <- which(read_thresholds == ninety_count_thresh) - 1

# Calculate xintercept to place the line between the current and next bar
xintercept <- position + 0.5

freq_plot <- ggplot(plot_data, aes(x = factor(ReadThreshold), y = Percentage)) +
  geom_bar(stat = 'identity', fill = "grey") +
  xlab('Read Count Threshold') +
  ylab('Percentage of Rows') +
  ggtitle(paste('Percentage of Rows with >=90% Samples Above Read Count Threshold (n=', total_rows, ')', sep = '')) +  # Adding n=total_rows in title
  theme_bw() +  # Remove grey background
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.grid.major = element_blank(),  # Optionally remove gridlines
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank()
  ) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 105)) +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")), 
            vjust = -0.5, hjust = 0.5, size = 3) +
  geom_vline(xintercept = xintercept, linetype = "dashed", color = "#1F78B4", size = 1) +  # Vertical dashed line at x=20 (position 4 when step size is 5)
  annotate("text", x = xintercept+0.3, y = max(plot_data$Percentage) + 2, 
           label = paste0("90% samples with >= ",ninety_count_thresh," reads"), color = "#1F78B4", angle = 90, vjust = 0, hjust = 0.9) +
  scale_fill_identity()  # Use the custom color scale


# print length of counts_final_sample
print(paste("Length of counts_final_sample:", nrow(counts_final_sample)))
print(paste("Number of rows with NA in Geneid:", sum(is.na(counts_final_sample$Geneid))))


### rawcounts_processing.R: FILTER ROWS WITH LOW READ COUNTS in >90% of SAMPLES ----
# if "zhang_genes" in raw_counts_file_1, then
if(grepl("zhang_genes", raw_counts_file_1)){
  # Drop rows from counts_data that don't have at least 90% of columns with 10 reads or more
  keep <- rowSums(counts_final_sample[, sample_names] >= ninety_count_thresh) >= 0.9 * ncol(counts_final_sample[, sample_names])
  # Print kept versus total
  print(paste("Kept", sum(keep), "of", nrow(counts_final_sample), "rows"))
  counts_final_sample <- counts_final_sample[keep,]
} else{
  # Drop rows from counts_data that don't have at least 90% of columns with 10 reads or more
  keep <- rowSums(counts_final_sample[, sample_names] >= ninety_count_thresh) >= 0.9 * ncol(counts_final_sample[, sample_names])
  # Print kept versus total
  print(paste("Kept", sum(keep), "of", nrow(counts_final_sample), "rows"))
  counts_final_sample <- counts_final_sample[keep,]
}

# print length of counts_final_sample
print(paste("Length of counts_final_sample:", nrow(counts_final_sample)))
print(paste("Number of rows with NA in Geneid:", sum(is.na(counts_final_sample$Geneid))))

### rawcounts_processing.R: Normalize and prep for export ----
sample_names <- rownames(sample_data)#dplyr::select(sample_data, c(1)))
#sample_names <- sample_names[sample_names != 'X']

if(all(sample_names %in% colnames(counts_final_sample))==FALSE){
  print("Error: Some samples in your conditions file do not exist in your counts file. Please correct this and re-run.")
  stop()
}

### RE-ORDER SAMPLE DATA FILE BY SAMPLE NAME AND BY CLASS IN ASCENDING ORDER
sample_data_ordered <- sample_data[order(sample_data$condition, rownames(sample_data)), ]
#sample_data_ordered <- sample_data_ordered[!rownames(sample_data_ordered) %in% excluded_samples, ]

# Reorder counts_final_sample columns
counts_final_ordered <- counts_final_sample[, rownames(sample_data_ordered)]

### RE-ORDER SAMPLES IN COUNTS FILE BY SAMPLE NAME AND BY CLASS IN ASCENDING ORDER
class_vector <- counts_final_ordered[1, ]
class_vector[1, ] <- names(counts_final_ordered)
class_vector[] <- sample_data_ordered$condition[match(unlist(class_vector), rownames(sample_data_ordered))]
class_vector[1,1] <- "sample_ORDER"

# Sort by sample name first
class_vector_ordered <- class_vector#[, c(1, order(colnames(class_vector[, 2:ncol(class_vector)])) + 1)]

# Then sort by sample class
#class_vector_ordered <- class_vector_ordered[, c(1, order(unlist(class_vector_ordered[1, 2:ncol(class_vector_ordered)], use.names=FALSE)) + 1)]

#counts_final_ordered <- counts_final_sample[, colnames(class_vector_ordered)]


#colnames(sample_data_ordered)[1] <- ""

### Combat-Seq normalization if enabled
if (combat_seq_normalize) {
  counts_final_sample_combat <- counts_final_sample
  counts_final_matrix <- as.matrix(counts_final_ordered)
  batch_vector <- as.numeric(factor(sample_data_ordered$batch))
  
  combat_counts <- ComBat_seq(counts_final_matrix, 
                              batch = sample_data_ordered$batch, 
                              group = sample_data_ordered$condition)
  counts_final_sample_combat[-1] <- as.data.frame(combat_counts)
  # reorder counts_final_sample columns to match order of sample_data_ordered
  counts_final_sample_combat <- counts_final_sample_combat[, rownames(sample_data_ordered)]
  counts_final_sample <- cbind(counts_final_sample$Geneid,counts_final_sample_combat)
} else{
  counts_final_sample_combat <- counts_final_sample[,c("Geneid",rownames(sample_data_ordered))]
  counts_final_sample <- counts_final_sample[,c("Geneid",rownames(sample_data_ordered))]
}

# print length of counts_final_sample
print(paste("Length of counts_final_sample:", nrow(counts_final_sample)))
print(paste("Number of rows with NA in Geneid:", sum(is.na(counts_final_sample$Geneid))))

# rename "counts_final[keep,]$Geneid" column to "Geneid"
colnames(counts_final_sample)[1] <- "Geneid"



### rawcounts_processing.R: CREATE GSEA FILES ----
sample_count <- nrow(sample_data_ordered)
class1_name <- sample_data_ordered[2, 1]
class2_name <- sample_data_ordered[nrow(sample_data_ordered), 1]
class1_count <- sum(sample_data_ordered$condition == class1_name)
class2_count <- sum(sample_data_ordered$condition == class2_name)

ssGSEA_cls_list <- list(row1=c(sample_count , "2" , "1"), row2=c("#" , class1_name , class2_name), row3=c(rep(0, class1_count), rep(1, class2_count)))
ssGSEA_cls <- as.data.frame(do.call(rbind, ssGSEA_cls_list))
rownames(ssGSEA_cls) <- c()
names(ssGSEA_cls) <- NULL
ssGSEA_cls[1, 4:ncol(ssGSEA_cls)] <- c(rep("", ncol(ssGSEA_cls) - 3))
ssGSEA_cls[2, 4:ncol(ssGSEA_cls)] <- c(rep("", ncol(ssGSEA_cls) - 3))

# ssGSEA needs a gene length file for TPM normalization
gene_length <- raw_counts_data_1[, c('Geneid', 'Length')]

### CREATE GSEA CLASS FILES
GSEA_cls_list <- list(row1=c(sample_count , "2" , "1"), row2=c("#" , class1_name , class2_name), row3=c(rep(class1_name, class1_count), rep(class2_name, class2_count)))
GSEA_cls <- as.data.frame(do.call(rbind, GSEA_cls_list))
rownames(GSEA_cls) <- c()
names(GSEA_cls) <- NULL
GSEA_cls[1, 4:ncol(GSEA_cls)] <- c(rep("", ncol(GSEA_cls) - 3))
GSEA_cls[2, 4:ncol(GSEA_cls)] <- c(rep("", ncol(GSEA_cls) - 3))

### CREATE TPM FILE FOR ssGSEA
# counts_final_ordered_rownames <- counts_final_sample
# rownames(counts_final_ordered_rownames) <- counts_final_sample[, 1]
# 
# gene_length_rownames <- gene_length
# rownames(gene_length_rownames) <- gene_length[, 1]
# gene_length_rownames <- within(gene_length_rownames, rm(Geneid))
# 
# gene_number <- nrow(counts_final_ordered_rownames)
# 
# counts_to_tpm <- function (counts_final_ordered_rownames, gene_length_rownames) {
#   x <- counts_final_ordered_rownames / gene_length_rownames
#   return (t(t(x) * 1e6 / colSums(x)))
# }
# 
# tpm <- counts_to_tpm(counts_final_ordered_rownames, gene_length_rownames[, 1])
# tpm_dataframe <- as.data.frame(tpm)
# 
# description = matrix(c(rep("na", gene_number)), gene_number, 1)
# tpm_genepattern <- add_column(tpm_dataframe, description, .before=1)
# tpm_genepattern <- rownames_to_column(tpm_genepattern, var = "NAME")



### rawcounts_processing.R: OUTPUT RESUTLING FILES ----
# Create root output folder if it doesn't exist
if (!file.exists("./Output/")) {
  dir.create("./Output/")
}

if (!file.exists("./Output/Raw Data Processing/")) {
  dir.create("./Output/Raw Data Processing/")
}

dir.create(paste("./Output/Raw Data Processing/", class1_name, "_", class2_name, "_", file_version, sep=""))

# Output counts .csv file for DESeq2 normalization followed by GSEA -or- counts to TPM conversion for ssGSEA
write.csv(counts_final_sample, file = paste("./Output/Raw Data Processing/", class1_name, "_", class2_name, "_", file_version, "/", class1_name, "_", class2_name, "_DESeq2_rawcounts.csv", sep = ""), row.names = FALSE, quote = FALSE)

# Output conditions .csv file for DESeq2 normalization
write.csv(sample_data_ordered, file = paste("./Output/Raw Data Processing/", class1_name, "_", class2_name, "_", file_version, "/", class1_name, "_", class2_name, "_DESeq2_conditions.csv", sep = ""), row.names = TRUE, quote = FALSE)

# Output phenotype .cls file for GSEA 
write.table(GSEA_cls, file = paste("./Output/Raw Data Processing/", class1_name, "_", class2_name, "_", file_version, "/", class1_name, "_", class2_name, "_GSEA_phenotype.cls", sep = ""), quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)

# Output normalized counts .txt file for GSEA
#write.table(tpm_genepattern, file = paste("./Output/Raw Data Processing/", class1_name, "_", class2_name, "_", file_version, "/", class1_name, "_", class2_name, "_ssGSEA_tpm.tsv", sep = ""), quote = FALSE, sep = '\t', row.names = FALSE)
# print columns in counts_final_sample[-1]
print(colnames(counts_final_sample[-1]))
print(rownames(sample_data_ordered))
# DESeq2 normalization without scaling
dds <- DESeqDataSetFromMatrix(countData = counts_final_sample[-1], colData = sample_data_ordered, design = ~ 1)
dds <- DESeq(dds)
normalized_counts <- counts(dds, normalized = TRUE)
normalized_counts_df <- as.data.frame(normalized_counts)
normalized_counts_df <- cbind(Geneid = counts_final_sample$Geneid, normalized_counts_df)

# Output DESeq2 normalized counts .csv file
write.csv(normalized_counts_df, file = paste("./Output/Raw Data Processing/", class1_name, "_", class2_name, "_", file_version, "/", class1_name, "_", class2_name, "_DESeq2_normcounts.csv", sep = ""), row.names = FALSE, quote = FALSE)

# output freq_plot and length_plot to file
ggsave(paste("./Output/Raw Data Processing/", class1_name, "_", class2_name, "_", file_version, "/", "frequency_plot_noX.png", sep = ""), freq_plot, width = 8, height = 6, units = "in")
ggsave(paste("./Output/Raw Data Processing/", class1_name, "_", class2_name, "_", file_version, "/", "length_plot_noX.png", sep = ""), length_plot, width = 8, height = 6, units = "in")
# also save as svg
ggsave(paste("./Output/Raw Data Processing/", class1_name, "_", class2_name, "_", file_version, "/", "frequency_plot_noX.svg", sep = ""), freq_plot, width = 8, height = 6, units = "in")
ggsave(paste("./Output/Raw Data Processing/", class1_name, "_", class2_name, "_", file_version, "/", "length_plot_noX.svg", sep = ""), length_plot, width = 8, height = 6, units = "in")


### GENERATE PEAKS PER GENOMIC FEATURE TYPE BOXPLOT ----
if(TRUE){
  source("~/5hmC-Pathway-Analysis/R Code/repeated_glmnet_funcs_v1.R")
  # define path of output
  if(comparison_set == "AM"){
    output_dir_name <- "~/peritoneal_processing/macs_6_30_2024_p_1e5/AMpos_AMneg_excl"
    total_read_perbam <- "total_reads_11_23_2024.csv"
    merged_peak_file_saf <- file.path(output_dir_name, saf_fn)
  }else if(comparison_set =="PM"){
    output_dir_name <- "~/peritoneal_processing/macs_6_30_2024_p_1e5/CRC_ADE_HEALTHY/peaks_per_featuretype/"
    total_read_perbam <- "total_reads_11_23_2024.csv"
    merged_peak_file_saf <- file.path(output_dir, saf_fn)
  }
  
  # define total_read_perbam path
  total_read_perbam_path <- file.path(output_dir_name, total_read_perbam)
  
  ## Load Sample File
  sample_file <- metadata_file
  #excluded_samples <- c("KT026", "KT26","KT027","KT27")
  selected_conds <- selected_conditions
  sample_data <- read.csv(
    sample_file,
    strip.white = TRUE,
  )
  
  ## Code to copy .xls files
  # # Extract unique elements from "X" column
  # x_values <- unique(sample_data$X)
  # 
  # # Define the parent directory of output_dir_name
  # parent_dir <- dirname(output_dir_name)
  # 
  # # Get list of .xls files in the parent directory
  # xls_files <- dir(parent_dir, pattern = "\\.xls$", full.names = TRUE)
  # library(fs)
  # # Loop through each .xls file and check if it matches any element in x_values
  # for (file in xls_files) {
  #   file_name <- basename(file)
  #   match_found <- any(sapply(x_values, function(x) grepl(x, file_name)))
  #   
  #   # If the file name contains an element from x_values, copy it to output_dir_name
  #   if (match_found) {
  #     file_copy(file, output_dir_name)
  #   }
  # }
  
  ## Select relevant samples
  # Keep only rows in sample_data where condition is in selected_conds
  sample_data <- sample_data[sample_data$condition %in% selected_conds, ]
  
  # make column X the rownames and drop col X
  rownames(sample_data) <- sample_data$X
  sample_data <- sample_data[, -1]
  # print number of rows
  # Drop rows where row names is in excluded samples
  sample_data <- sample_data[!rownames(sample_data) %in% excluded_samples,]
  
  # print number of samples in each condition
  print(table(sample_data$condition))
  
  # get file names ending in .xls and containing one of the sample names in folder: ~/peritoneal_processing/macs_6_30_2024_p_1e5
  peak_files <- list.files(output_dir, pattern = ".xls", full.names = TRUE)
  #peak_files <- peak_files[grepl(paste(rownames(sample_data), collapse = "|"), peak_files)]
  
  ## Required if intersect_filter is TRUE, otherwise optional
  # load consensus peak .saf file
  
  merged_peak_file_df <- read.table(merged_peak_file_saf, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  # Assign correct column names
  colnames(merged_peak_file_df) <- c("GeneID", "Chr", "Start", "End", "Strand")
  
  # Optional if the .saf file was filtered
  # drop rows where GeneID not in counts_final_sample$Geneid
  merged_peak_file_df <- merged_peak_file_df[merged_peak_file_df$GeneID %in% counts_final_sample$Geneid,]
  
  # Create GRanges object
  merged_peak_file_gr <- makeGRangesFromDataFrame(merged_peak_file_df, 
                                                  keep.extra.columns = TRUE,
                                                  seqnames.field = "Chr",
                                                  start.field = "Start",
                                                  end.field = "End",
                                                  strand.field = "Strand",
                                                  ignore.strand = FALSE)
  
  # Function to process a single peak file
  process_peak_file <- function(peak_file, merged_peak_file_gr, intersect_filter = TRUE) {
    output_fn = file.path(output_dir_name, gsub(".xls", ".csv", basename(peak_file)))
    # if folder does not exist, create it
    if (!dir.exists(dirname(output_fn))) {
      dir.create(dirname(output_fn))
    }
    
    message("Processing file:", peak_file)
    
    if(!file.exists(output_fn)){
      message("File does not exist:", output_fn)
      
      peak_data <- read.table(peak_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      print(head(peak_data))
      peak_data <- peak_data[peak_data$length <= 5000, ]
      
      # Convert peak data to GRanges object
      message("Converting peak data to GRanges object...")
      tryCatch({
        peak_data_gr <- makeGRangesFromDataFrame(peak_data, keep.extra.columns = TRUE)
        message("GRanges object successfully created.")
      }, error = function(e) {
        stop("Error in creating GRanges object: ", e$message)
      })
      
      # Intersect with merged_peak_file_gr if intersect_filter is TRUE
      if(intersect_filter){
        message("Performing intersection with merged_peak_file_gr...")
        overlaps <- findOverlaps(peak_data_gr, merged_peak_file_gr)
        peak_data <- peak_data[queryHits(overlaps),]
      } else {
        message("Skipping intersection as intersect_filter is FALSE.")
      }
      
      #peak_data <- peak_data[grepl("chr[1-9]|chr1[0-9]|chr2[0-2]", peak_data$chr),]
      # drop rows where chr is in excluded_chr
      # print number of rows
      # print("Before filtering chrs:")
      # print(head(peak_data))
      # print(nrow(peak_data))
      # # print head of rows with NA in start or end
      # print(head(peak_data[is.na(peak_data$start) | is.na(peak_data$end),]))
      peak_data <- peak_data[!peak_data$chr %in% excluded_chr,]
      peak_data <- peak_data[!grepl("random", peak_data$chr),]
      peak_data <- peak_data[!grepl("chrUn", peak_data$chr),]
      # print("After filtering chrs:")
      # print(head(peak_data))
      # print(nrow(peak_data))
      # print(head(peak_data[is.na(peak_data$start) | is.na(peak_data$end),]))
      # # print rows where chr is NA or empty
      # print("Rows where chr column is not of format chr1-22:")
      # print(head(peak_data[!grepl("chr[1-9]|chr1[0-9]|chr2[0-2]", peak_data$chr),]))
      
      feature_data <- data.frame(feature = paste(peak_data$chr, peak_data$start, peak_data$end, sep = "_"))
      # print any na rows in feature_data
      # print("NA rows in feature_data:")
      # print(head(feature_data[is.na(feature_data$feature),]))
      
      annotated_df <- annotate_features2(feature_data, annots = c('hg19_basicgenes'))
      
      closest_genes <- annotated_df$annot.symbol
      genomic_features <- annotated_df$locationType
      
      if(sum(is.na(closest_genes)) > 0){
        message("There are missing gene symbols. Replacing with 'NA'.")
        gene_symbols <- mapIds(org.Hs.eg.db, keys = closest_genes, keytype = "ENTREZID", column = "SYMBOL")
        gene_symbols[is.na(gene_symbols)] <- "NA"
      } else {
        message("No missing gene symbols.")
        gene_symbols <- closest_genes
      }
      print(head(annotated_df[is.na(annotated_df$start),]))
      peak_df_anno <- data.frame(
        geneID = gene_symbols,
        chr = annotated_df$chr,
        start = annotated_df$start,
        end = annotated_df$end,
        strand = rep(".", nrow(annotated_df)),
        genomic_feature = genomic_features,
        closest_genes = closest_genes,
        stringsAsFactors = FALSE
      )
      
      write.csv(peak_df_anno, file = output_fn, row.names = FALSE)
      message("Processed file: ", peak_file)
    } else {
      message("Annotated peak file already exists. Skipping annotation:", output_fn)
    }
  }
  
  # Use mclapply to process files in parallel
  mclapply(peak_files, function(peak_file) {
    process_peak_file(peak_file, merged_peak_file_gr, intersect_filter = FALSE)  # Set to FALSE to skip intersection
  }, mc.cores = 50)
  
  
  ## for each peak file, load into a dataframe and count the number of peaks per genomic feature type into a new dataframe
  
  # Create list of annotated peak files
  annotated_peak_files <- list.files(output_dir_name, pattern = "peaks.csv", full.names = TRUE)
  
  # Initialize an empty dataframe
  peak_counts <- data.frame()
  
  # For each annotated_peak_file, load into a dataframe and count the number of peaks per genomic feature type
  for (i in seq_along(annotated_peak_files)) {
    peak_file <- annotated_peak_files[i]  # Changed from [1] to [i] to iterate through all files
    peak_data <- read.csv(peak_file, header = TRUE, stringsAsFactors = FALSE)
    
    # Extract sample name from the filename
    sample_name <- gsub(".csv", "", basename(peak_file))
    sample_name <- gsub("annotated_peaks_", "", sample_name)
    sample_name <- gsub("_bowtie2.bam_macs3_peaks", "", sample_name)
    
    # Count the number of peaks per genomic feature type
    feature_counts <- table(peak_data$genomic_feature)
    
    # If peak_counts is empty, initialize it with correct column names
    if (nrow(peak_counts) == 0) {
      col_names <- c("sample", names(feature_counts))
      peak_counts <- data.frame(matrix(ncol = length(col_names), nrow = 0))
      colnames(peak_counts) <- col_names
    }
    
    # Create a new row with the sample name and feature counts
    new_row <- c(sample_name, as.vector(feature_counts))
    
    # Add the new row to peak_counts, using sample_name as row name
    peak_counts[sample_name, ] <- new_row
  }
  
  # set sample_data on rownames
  rownames(peak_counts) <- peak_counts$sample
  peak_counts <- peak_counts[, -1]
  
  # sort columns in this order: 1-5kb Upstream, Promoter, 5' UTR, Exon, Intron, 3' UTR, lncRNA, Enhancer, Intergenic
  peak_counts <- peak_counts[, c("1-5kb Upstream", "Promoter", "5' UTR", "Exon", "Intron", "3' UTR", "Intergenic")] #"lncRNA", "Enhancer", 
  
  # merge peak_counts with sample_data on rownames, keeping only rows in sample_data
  peak_counts_merged <- merge(data.frame(
    condition = sample_data[,"condition"],
    row.names = rownames(sample_data)
  ), peak_counts, by = "row.names", all = TRUE)
  
  rownames(peak_counts_merged) <- peak_counts_merged$Row.names
  peak_counts_merged$Row.names <- NULL
  # drop column specifically called Row.names
  
  peak_counts_merged[, 2:ncol(peak_counts_merged)] <- lapply(peak_counts_merged[, 2:ncol(peak_counts_merged)], as.numeric)
  
  # add total column
  peak_counts_merged$total <- rowSums(peak_counts_merged[,2:ncol(peak_counts_merged)])
  
  #print nrows
  ### Normalize
  library(parallel)
  
  if (file.exists(total_read_perbam_path)) {
    total_reads <- read.csv(total_read_perbam_path, row.names = 1)
  } else {
    # Function to get total reads from BAM file
    get_total_reads <- function(bam_file) {
      # Define flags
      # Exclude unmapped (4), secondary alignments (256), duplicates (1024), supplementary alignments (2048)
      exclude_flags <- 4 + 256 + 1024 + 2048  # Sum to 3332
      # Include properly paired reads (2) and first read in pair (64)
      include_flags <- 2 + 64  # Sum to 66
      # Construct command
      cmd <- paste("/tools/bin/samtools view -c -f", include_flags, "-q 20 -F", exclude_flags, bam_file)
      total <- as.numeric(system(cmd, intern = TRUE))
      return(total / 1e6)  # Convert to millions of reads
    }
    
    # Get BAM files and sample names
    bam_files <- metadata$bamReads
    
    # Determine the number of cores to use (leave one core free)
    num_cores <- 50
    
    total_reads <- unlist(mclapply(bam_files, get_total_reads, mc.cores = num_cores))
    write.csv(total_reads, file = total_read_perbam_path)
  }
  
  # Extract KT### numbers using a regular expression
  kt_numbers <- sub(".*/(KT[0-9]+)_.*", "\\1", bam_files)
  
  # Display the result
  print(kt_numbers)
  
  # Extract the 'x' column from total_reads as a named vector
  total_reads_vector <- as.vector(total_reads)#$x)
  names(total_reads_vector) <- kt_numbers #rownames(total_reads)
  total_reads <- total_reads_vector
  
  # drop excluded samples from total_reads
  total_reads <- total_reads[!names(total_reads) %in% excluded_samples]
  # drop excluded samples from peak_counts_merged
  peak_counts_merged <- peak_counts_merged[rownames(peak_counts_merged) %in% names(total_reads),]
  
  peak_counts_merged_norm <- peak_counts_merged
  
  # Keep only "Total" column in peak_counts_merged_norm
  #peak_counts_merged_norm <- peak_counts_merged_norm[, c("condition", "Exon")]
  
  # save 
  # Check if all rownames in peak_counts_merged_norm are in total_reads
  missing_samples <- setdiff(rownames(peak_counts_merged_norm), names(total_reads))
  if (length(missing_samples) > 0) {
    stop("The following samples are missing in total_reads: ", paste(missing_samples, collapse = ", "))
  }
  
  # Normalize peak counts to peaks per million reads
  for (col in 2:ncol(peak_counts_merged_norm)) {
    sample_names <- rownames(peak_counts_merged_norm)
    peak_counts_merged_norm[, col] <- peak_counts_merged_norm[, col] / total_reads[sample_names]
  }
  
  ### END NORM
  library(reshape2)
  library(ggpubr) # for stat_pvalue_manual function
  ## plot boxplot using plotly of peak counts per genomic feature type colored by condition
  peak_counts_melted <- melt(peak_counts_merged_norm, id.vars = c("condition"))
  #peak_counts_melted <- melt(peak_counts_merged, id.vars = c("condition"))
  
  peak_counts_melted$condition <- as.factor(peak_counts_melted$condition)
  peak_counts_melted$variable <- as.factor(peak_counts_melted$variable)
  
  # set value to integer
  peak_counts_melted$value <- as.numeric(peak_counts_melted$value)
  
  # Determine y-axis range programmatically, handling NA and zero values
  valid_values <- peak_counts_melted$value[!is.na(peak_counts_melted$value) & peak_counts_melted$value > 0]
  
  if(length(valid_values) > 0) {
    y_min <- floor(log10(min(valid_values)))
    y_max <- ceiling(log10(max(valid_values)))
  } else {
    # Fallback if there are no valid values
    y_min <- 0
    y_max <- 1
    warning("No valid positive values found for y-axis scaling. Using default range.")
  }
  
  # Add some padding to the y-axis limits
  y_min <- max(y_min, -10)  # Limit the lower bound to avoid going too low
  y_max <- min(y_max, 10)   # Limit the upper bound to avoid going too high
  
  library(ggplot2)
  library(scales)  # For log_breaks()
  
  # Define custom colors with more muted blue and red
  custom_colors <- c("PM_negative" = "#1F78B4A6", "PM_positive" = "#E31A1CA6","HEALTHY" = "grey",
                     "AMpos" = "#E31A1CA6","AMneg"= "#1F78B4A6")  # Cornflower Blue and Tomato
  
  
  # Custom function to generate breaks and labels
  log10_breaks <- function(x) {
    minx <- floor(log10(min(x)))
    maxx <- ceiling(log10(max(x)))
    
    major_breaks <- 10^(minx:maxx)
    major_breaks <- as.integer(major_breaks)
    minor_breaks <- rep(1:9, maxx - minx + 1) * rep(major_breaks, each = 9) / 10
    
    list(
      breaks = major_breaks,
      minor_breaks = minor_breaks[minor_breaks > min(x) & minor_breaks < max(x)],
      labels = scales::label_number(accuracy = 1)(major_breaks)
    )
  }
  
  # Create the plot
  p <- ggplot(peak_counts_melted, aes(x = variable, y = value, fill = condition)) +
    geom_jitter(
      aes(color = condition),
      position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
      size = 1,
      alpha = 0.5
    ) +
    geom_boxplot(
      outlier.shape = NA,
      position = position_dodge(width = 0.8),
      width = 0.7
    ) +
    scale_fill_manual(values = custom_colors) +
    scale_color_manual(values = custom_colors) +
    scale_y_log10(
      breaks = log10_breaks(c(20,6000))$breaks,
      minor_breaks = log10_breaks(c(20,6000))$minor_breaks
    ) +
    # # set y range
    # coord_cartesian(ylim = c(30, 6000)) +
    # scale_y_continuous(
    #   limits = c(50, 500),
    #   breaks = seq(50, 500, by = 50)  # Adjust the breaks as needed
    # ) +
    theme_minimal() +
    labs(
      title = "5hmC Peaks by Genomic Feature Type",
      x = "Genomic Feature Type",
      y = "Peaks per million reads (log10)",
      fill = "Condition",
      color = "Condition"
    ) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8),
      axis.text.y = element_text(size = 8),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      legend.position = "right",
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_line(color = "grey95"),
      panel.border = element_blank(),
      plot.margin = margin(l = 50, r = 10, b = 10, t = 10, unit = "pt")
    )
  
  # Display the plot
  print(p)
  
  ## Statistical Testing
  # Summary statistics
  summary_stats <- peak_counts_melted %>%
    group_by(condition, variable) %>%
    summarise(
      samples = n(),
      mean = mean(value),
      sd = sd(value),
      median = median(value),
      IQR = IQR(value),
      .groups = 'drop'
    )
  
  # Define the statistical test to be used
  test_used <- "t-test"  # Change to "Wilcoxon test" if needed
  
  # Perform pairwise t-tests for each variable
  test_results_list <- list()
  
  for (var in unique(peak_counts_melted$variable)) {
    data_var <- subset(peak_counts_melted, variable == var)
    
    # Pairwise t-test with Holm correction
    pairwise_test <- pairwise.t.test(data_var$value, data_var$condition,
                                     p.adjust.method = "holm", pool.sd = TRUE)
    
    # Convert to data frame
    test_results <- as.data.frame(pairwise_test$p.value)
    test_results$group1 <- rownames(test_results)
    
    # Reshape to long format
    test_results_long <- test_results %>%
      pivot_longer(cols = -group1, names_to = "group2", values_to = "p") %>%
      mutate(variable = var) %>%
      filter(!is.na(p))
    
    # Append to list
    test_results_list[[var]] <- test_results_long
  }
  
  # Combine all results
  test_results_df <- bind_rows(test_results_list)
  
  # Merge with summary statistics
  test_results_combined <- test_results_df %>%
    left_join(summary_stats, by = c("variable", "group1" = "condition")) %>%
    rename(samples1 = samples, mean1 = mean, sd1 = sd, median1 = median, IQR1 = IQR) %>%
    left_join(summary_stats, by = c("variable", "group2" = "condition")) %>%
    rename(samples2 = samples, mean2 = mean, sd2 = sd, median2 = median, IQR2 = IQR)
  
  # Add the test used
  test_results_combined$test <- test_used
  
  # Compute y.position for plotting
  max_y_values <- peak_counts_melted %>%
    group_by(variable) %>%
    summarize(max_value = max(value, na.rm = TRUE), .groups = 'drop')
  
  # Merge and compute y.position
  test_results_combined <- test_results_combined %>%
    left_join(max_y_values, by = "variable") %>%
    group_by(variable) %>%
    arrange(variable, p) %>%
    mutate(
      y.position = max_value + (row_number()) * (max_value * 0.05)
    ) %>%
    ungroup()
  
  # Reorder columns
  test_results_combined <- test_results_combined %>%
    dplyr::select(variable, group1, group2,
           samples1, mean1, sd1, median1, IQR1,
           samples2, mean2, sd2, median2, IQR2,
           p, test, y.position) %>%
    arrange(variable, group1, group2)
  
  # Debug: Check the structure and head of test_results_combined
  print("Structure of test_results_combined:")
  str(test_results_combined)
  print("Head of test_results_combined:")
  print(head(test_results_combined))

  
  ### Individual Region Plots with p-values from CSV
  # Create the plot with facets
  p2 <- ggplot(peak_counts_melted, aes(x = condition, y = value, fill = condition)) +
    geom_boxplot(outlier.shape = NA, width = 0.7, alpha = 0.8) +
    geom_jitter(aes(color = condition), width = 0.2, size = 1, alpha = 0.5) +
    scale_fill_manual(values = custom_colors) +
    scale_color_manual(values = custom_colors) +
    facet_wrap(~ variable, scales = "free_y", ncol = 3) +
    theme_minimal() +
    labs(
      title = "5hmC Peaks by Genomic Feature Type",
      x = "Condition",
      y = "Peaks per million reads (log10)",
      fill = "Condition",
      color = "Condition"
    ) +
    theme(
      axis.text.y = element_text(size = 8),
      axis.title.y = element_text(size = 10),
      legend.position = "right",
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_line(color = "grey95"),
      panel.border = element_blank(),
      strip.background = element_rect(fill = "grey95"),
      strip.text = element_text(size = 10, face = "bold"),
      plot.margin = margin(l = 50, r = 10, b = 10, t = 10, unit = "pt")
    )
  
  # Debug: Check if 'test_results_combined' has the necessary columns
  required_columns <- c("variable", "group1", "group2", "p", "y.position")
  missing_columns <- setdiff(required_columns, colnames(test_results_combined))
  if(length(missing_columns) > 0){
    stop("The following required columns are missing in test_results_combined: ", 
         paste(missing_columns, collapse = ", "))
  }
  
  # Debug: Print a sample of test_results_combined to verify
  print("Sample of test_results_combined for plotting:")
  print(head(test_results_combined))
  
  # Add p-values annotation with appropriate labels
  test_results_combined <- test_results_combined %>%
    mutate(p_label = case_when(
      p < 0.001 ~ "***",
      p < 0.01 ~ "**",
      p < 0.05 ~ "*",
      TRUE ~ "NS"
    ))
  
  # Add p-values from CSV file using stat_pvalue_manual
  p2 <- p2 + stat_pvalue_manual(
    data = test_results_combined,
    label = "p_label",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    step.increase = 0.1,
    size = 3,
    inherit.aes = FALSE  # Important to prevent inheriting main plot aesthetics
  )
  
  # Display the individual region plot
  print(p2)
  
  ## Save test results to CSV
  write.csv(test_results_combined, 
            file = "~/5hmC-Pathway-Analysis/Output/feature_type/peak_counts_per_genomic_feature_type_PMposPMneg_HEALTHY_10_29_24.csv", 
            row.names = FALSE)
  
  # Save the plots
  ggsave("~/5hmC-Pathway-Analysis/Output/feature_type/peak_counts_per_genomic_feature_type_PMposPMneg_HEALTHY_10_29_24.png", 
         plot = p, width = 8.5, height = 3.5, units = "in", dpi = 300)
  ggsave("~/5hmC-Pathway-Analysis/Output/feature_type/peak_counts_per_genomic_feature_type_PMposPMneg_HEALTHY_10_29_24.svg", 
         plot = p, width = 8.5, height = 3.5, units = "in")
  
  # Save the individual region plot
  ggsave("~/5hmC-Pathway-Analysis/Output/feature_type/PM_peak_counts_per_genomic_feature_type_subplots_with_ggsignif_adjusted_10_29_24.png", 
         plot = p2, width = 15, height = 12, units = "in", dpi = 300)
  ggsave("~/5hmC-Pathway-Analysis/Output/feature_type/PM_peak_counts_per_genomic_feature_type_subplots_with_ggsignif_adjusted_10_29_24.svg", 
         plot = p2, width = 15, height = 12, units = "in")
}

### Generate piecharts with "(n= )" in the title ----
source("~/5hmC-Pathway-Analysis/R Code/repeated_glmnet_funcs_v1.R")
# Load consensus peak .saf file
#merged_peak_file_saf <- file.path(output_dir, "merged_peaks_10perc_pe5_092824.saf")
merged_peak_file_saf <- file.path(output_dir, saf_fn)
merged_peak_file_df <- read.table(merged_peak_file_saf, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# Only keep rows where V1 is in counts_final_sample$Geneid
merged_peak_file_df <- merged_peak_file_df[merged_peak_file_df$V1 %in% counts_final_sample$Geneid,]
# Print number of rows
print(nrow(merged_peak_file_df))
print(nrow(counts_final_sample))

# Print number of rows where counts_final_sample$Geneid is NA
print(sum(is.na(counts_final_sample$Geneid)))

# Convert V1 to a single column dataframe
merged_peak_file_df_plot_enh <- data.frame(feature = merged_peak_file_df$V1)
# Annotate features
merged_peak_file_df_plot_enh <- annotate_features2(merged_peak_file_df_plot_enh, annots = c('hg19_custom_genehancer'))

# Create a dataframe with 'feature' where locationType == Enhancer
merged_peak_file_df_plot_enh_only <- merged_peak_file_df_plot_enh[merged_peak_file_df_plot_enh$locationType == "Enhancer",]
merged_peak_file_df_plot_enh_only <- data.frame(feature = merged_peak_file_df_plot_enh_only$feature)
merged_peak_file_df_plot_enh_only <- annotate_features2(merged_peak_file_df_plot_enh_only, annots = c('hg19_basicgenes'))
# Add column called "enhancerType" that is equal to "Enhancer"
merged_peak_file_df_plot_enh_only$enhancerType <- "Enhancer"

merged_peak_file_df_plot_inter_only <- merged_peak_file_df_plot_enh[merged_peak_file_df_plot_enh$locationType == "Intergenic",]
merged_peak_file_df_plot_inter_only <- data.frame(feature = merged_peak_file_df_plot_inter_only$feature)
merged_peak_file_df_plot_inter_only <- annotate_features2(merged_peak_file_df_plot_inter_only, annots = c('hg19_basicgenes'))
# Add column called "enhancerType" that is equal to "Non-Enhancer"
merged_peak_file_df_plot_inter_only$enhancerType <- "Non-enhancer"

# Concatenate two dataframes into one called merged_peak_file_df_plot
merged_peak_file_df_plot <- rbind(merged_peak_file_df_plot_enh_only, merged_peak_file_df_plot_inter_only)
# Keep only specified columns
merged_peak_file_df_plot <- merged_peak_file_df_plot[, c("feature", "chr", "start", "end", "all_symbols", "locationType", "enhancerType")]

# Sort locationType in the specified order
merged_peak_file_df_plot$locationType <- factor(merged_peak_file_df_plot$locationType, levels = c("1-5kb Upstream", "Promoter", "5' UTR", "Exon", "Intron", "3' UTR", "Intergenic"))

# Load required libraries
library(ggplot2)
library(dplyr)
library(scales)

# Create stacked bar plot for enhancerType
enhancer_counts <- merged_peak_file_df_plot %>%
  count(enhancerType) %>%
  mutate(percentage = n / sum(n) * 100)

# Get number of rows for the title
n_total <- nrow(merged_peak_file_df_plot)

enh_plot <- ggplot(enhancer_counts, aes(x = "", y = n, fill = enhancerType)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_text(aes(label = paste0(enhancerType, "\n", sprintf("%.1f%%", percentage))),
            position = position_stack(vjust = 0.5),
            color = "white", size = 4) +
  scale_fill_manual(values = c("Enhancer" = "#696969", "Non-enhancer" = "#A9A9A9")) +
  scale_y_continuous(labels = comma_format(big.mark = ",")) +
  labs(title = paste("Enh Types (n =", n_total, ")"),
       y = "Peak Count",
       x = NULL) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 12),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

# Create pie chart for locationType
library(ggrepel)

source("~/5hmC-Pathway-Analysis/R Code/repeated_glmnet_funcs_v1.R")

# **Retrieve the color tables**
color_tables <- create_color_tables()

location_counts <- merged_peak_file_df_plot %>%
  count(locationType) %>%
  mutate(percentage = n / sum(n) * 100,
         label = paste0(sprintf("%.1f%%", percentage)))

pie_plot <- ggplot(location_counts, aes(x = "", y = n, fill = locationType)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  geom_label_repel(aes(label = label),
                   position = position_stack(vjust = 0.5),
                   size = 5, show.legend = FALSE) +
  scale_fill_manual(values = color_tables$location_type_table) +  # **Apply the new color table**
  labs(title = "Distribution of Location Types",
       fill = "Location Type") +
  theme_void() +
  theme(legend.position = "right",
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        # increase font of legend
        legend.text = element_text(size = 12))

# Save plots as PNG and SVG files
ggsave("~/5hmC-Pathway-Analysis/Output/feature_type/Enhancer_Type_Distribution_PM90cp15_11_9_24.png", enh_plot, width = 4, height = 6, units = "in", dpi = 300)
ggsave("~/5hmC-Pathway-Analysis/Output/feature_type/Enhancer_Type_Distribution_PM90cp15_11_9_24.svg", enh_plot, width = 4, height = 6, units = "in")
ggsave("~/5hmC-Pathway-Analysis/Output/feature_type/Location_Type_Distribution_PM90cp15_11_9_24.png", pie_plot, width = 6, height = 4, units = "in", dpi = 300)
ggsave("~/5hmC-Pathway-Analysis/Output/feature_type/Location_Type_Distribution_PM90cp15_11_9_24.svg", pie_plot, width = 6, height = 4, units = "in")


### GENERATE METAPLOTS -----
# Load required libraries
library(data.table)
library(GenomicFeatures)
library(dplyr)
library(rtracklayer)
library(ggplot2)
library(foreach)
library(doParallel)
# Load the TxDb object for the reference genome (e.g., hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# Ensure data.table is attached
if (!is.element("package:data.table", search())) attachNamespace('data.table')

# Ensure data.table is attached
if (!is.element("package:data.table", search())) attachNamespace('data.table')

# Set global variables
ROOT_FOLDER <- "~/5hmC-Pathway-Analysis/Output/Raw Data Processing/HEALTHY_PM_positive_10per_10172024_hmr_combat_pm_healthy/"
FILENAME_TRUNK <- "HEALTHY_PM_positive_DESeq2"
OUTPUT_FOLDER_METAGENE <- "~/5hmC-Pathway-Analysis/Output/metagene"
BIGWIG_FOLDER <- "~/peritoneal_processing/trimmed_data_bw_8_18_2024_bigwig_RPGC/" #"~/peritoneal_processing/trimmed_data_bw_9_18_2023_bigwig/"
peak_file_prefix <- paste0(OUTPUT_FOLDER_METAGENE, "/sampled_peaks")

# Debug function
debug_print <- function(message, data = NULL) {
  cat("\nDEBUG:", message, "\n")
  if (!is.null(data)) {
    print(head(data))
    cat("Dimensions:", dim(data), "\n")
  }
}

# Function to read and prepare metadata
read_metadata <- function(meta_name) {
  meta_tracks <- fread(meta_name)
  
  if (!"sample_id" %in% names(meta_tracks)) {
    meta_tracks[, sample_id := as.character(meta_tracks[[1]])]
  } else {
    meta_tracks[, sample_id := as.character(sample_id)]
  }
  
  meta_tracks[, `:=`(
    scaling_factor = 1,
    bigwig_fn = paste0(BIGWIG_FOLDER, sample_id, "_bowtie2.bw")
  )]
  
  setcolorder(meta_tracks, c("sample_id", setdiff(names(meta_tracks), "sample_id")))
  
  debug_print("Processed metadata:", meta_tracks)
  meta_tracks
}

prepare_gene_bodies <- function(n_samples = 0, gene_file = "all_genes", separate_strands = FALSE) {
  message("Starting prepare_gene_bodies function.")
  
  if (gene_file == "all_genes") {
    message("Using all genes from TxDb.Hsapiens.UCSC.hg19.knownGene.")
    all_gene_bodies <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
    
    # drop genes on chrX, chrY and chrM
    all_gene_bodies <- all_gene_bodies[!seqnames(all_gene_bodies) %in% excluded_chr]
    
    
    # Convert GRanges to data.frame and rename columns
    all_gene_bodies_df <- as.data.frame(all_gene_bodies)
    colnames(all_gene_bodies_df) <- c("V1", "V2", "V3", "V4", "V5")
    
    all_gene_bodies_df$V1 <- as.character(seqnames(all_gene_bodies))
    all_gene_bodies_df$V2 <- start(all_gene_bodies)
    all_gene_bodies_df$V3 <- end(all_gene_bodies)
    all_gene_bodies_df$V4 <- all_gene_bodies$gene_id
    all_gene_bodies_df$V5 <- as.character(strand(all_gene_bodies))
    
    
    set.seed(123)
    if (n_samples > 0 && n_samples <= 1) {
      n_samples_count <- ceiling(n_samples * nrow(all_gene_bodies_df))
      message(paste("Sampling", n_samples_count, "genes from", nrow(all_gene_bodies_df), "total genes."))
      all_gene_bodies_df <- all_gene_bodies_df[sample(1:nrow(all_gene_bodies_df), n_samples_count), ]
    }
  } else {
    message(paste("Reading gene data from file:", gene_file))
    all_gene_bodies_df <- fread(gene_file, header = FALSE)
    if(n_samples > 0 && n_samples <= 1){
      n_samples_count <- ceiling(n_samples * nrow(all_gene_bodies_df))
      message(paste("Sampling", n_samples_count, "genes from", nrow(all_gene_bodies_df), "total genes in file."))
      all_gene_bodies_df <- all_gene_bodies_df[sample(1:nrow(all_gene_bodies_df), n_samples_count), ]
    }
  }
  
  if (separate_strands) {
    message("Separating strands.")
    message("All gene bodies:")
    print(head(all_gene_bodies_df))
    
    if ("V5" %in% colnames(all_gene_bodies_df)) {
      if ("." %in% all_gene_bodies_df$V5) {
        message("Replacing '.' with '+' in strand information.")
        all_gene_bodies_df$V5[all_gene_bodies_df$V5 == "."] <- "+"
      }
      
      positive_strand <- all_gene_bodies_df[all_gene_bodies_df$V5 == "+", ]
      message(paste("Number of genes on positive strand:", nrow(positive_strand)))
      
      if (nrow(positive_strand) == 0) {
        stop("No positive strand data found.")
      }
      
      if (any(all_gene_bodies_df$V5 == "-")) {
        negative_strand <- all_gene_bodies_df[all_gene_bodies_df$V5 == "-", ]
        message(paste("Number of genes on negative strand:", nrow(negative_strand)))
      } else {
        message("No genes found on negative strand. Creating dummy negative strand data.")
        negative_strand <- positive_strand[1, , drop = FALSE]
        negative_strand$V5 <- "-"
        negative_strand$V2 <- 1
        negative_strand$V3 <- 2
      }
      return(list(positive = positive_strand, negative = negative_strand))
    } else {
      stop("Strand information (V5) not found in the gene data.")
    }
  } else {
    message("Not separating strands. Returning all gene bodies.")
    return(list(all = all_gene_bodies_df))
  }
}

metagene_analysis <- function(
    separate_strands = TRUE,
    gene_file = "all_genes",
    n_genes = 100,
    regenerate_matrices = FALSE,
    before = 500,
    after = 500,
    region_body_length = 500,
    binsize=10
) {
  
  print("Preparing gene bodies...")
  
  # Function to export gene bodies to BED file
  export_gene_bodies <- function(gene_bodies, output_file_prefix) {
    for (strand in names(gene_bodies)) {
      fwrite(
        as.data.table(gene_bodies[[strand]]),
        file = paste0(output_file_prefix, "_", strand, ".bed"),
        sep = "\t",
        col.names = FALSE
      )
      debug_print(paste(strand, "gene bodies exported to:"), 
                  paste0(output_file_prefix, "_", strand, ".bed"))
    }
  }
  
  # Function to parse metagene data
  parse_metagene <- function(matrix_file) {
    df <- fread(cmd = paste("zcat", matrix_file, "| grep -v '^#'"), 
                header = FALSE)
    
    numeric_cols <- sapply(df, is.numeric)
    numeric_cols[1:5] <- FALSE
    
    if(sum(numeric_cols) == 0) {
      stop("No numeric columns found in the matrix file.")
    }
    
    df_numeric <- df[, ..numeric_cols]
    
    col_means <- colMeans(df_numeric, na.rm = TRUE)
    
    as.data.table(t(col_means))
  }
  
  # Function to process all samples
  process_samples <- function(meta_tracks, output_metagene_fn) {
    result <- rbindlist(lapply(seq_len(nrow(meta_tracks)), function(i) {
      sample <- meta_tracks$sample_id[i]
      result_list <- lapply(names(gene_bodies), function(strand) {
        matrix_file <- paste0(output_metagene_fn, "_", sample, "_", strand, "_matrix.gz")
        data <- parse_metagene(matrix_file)
        data[, `:=`(sample_id = sample, strand = strand)]
        data
      })
      rbindlist(result_list)
    }))
    
    setkey(result, sample_id, strand)
    result
  }
  
  # Function to prepare data for plotting
  prepare_plot_data <- function(df_metagene, meta_tracks) {
    debug_print("Preparing plot data. Input data:", df_metagene)
    # Convert to data.table if not already
    
    # Convert to data.table if not already
    if (!is.data.table(df_metagene)) {
      df_metagene <- as.data.table(df_metagene)
    }
    
    long_data <- melt(df_metagene, id.vars = c("sample_id", "strand"), variable.name = "Pos", value.name = "Score")
    
    # Ensure long_data is a data.table
    if (!is.data.table(long_data)) {
      long_data <- as.data.table(long_data)
    }
    
    long_data[, Pos := as.numeric(gsub("V", "", Pos))]
    
    # Flip the x-axis for negative strand if separate_strands is TRUE
    if (separate_strands) {
      max_pos <- max(long_data$Pos)
      long_data[strand == "negative", Pos := max_pos - Pos + 1]
    }
    
    long_data <- merge(long_data, meta_tracks, by = "sample_id", all.x = TRUE)
    
    if (nrow(long_data) == 0) {
      stop("No data after merging. Check sample_id consistency between metagene data and metadata.")
    }
    
    long_data[, Score := Score / mean(Score), by = .(sample_id, strand)]
    
    # Average scores across strands
    aggregated_data <- long_data[, .(
      score_median = median(Score),
      score_mean = mean(Score),
      score_sd = sd(Score)
    ), by = .(Pos, condition)]
    
    setkey(aggregated_data, Pos, condition)
    aggregated_data
  }
  
  # Function to create the combined metagene plot with NO SMOOTHING ----
  # create_combined_metagene_plot_nosmooth <- function(plot_data) {
  #   debug_print("Creating combined plot with data:", plot_data)
  #   if (nrow(plot_data) == 0) {
  #     stop("No data available for plotting. Check previous steps.")
  #   }
  #   
  #   custom_colors <- c("#696969", "#3e75b0", "#d2321f")
  #   
  #   #construct labels. If gene_file == "all_genes"
  #   if (gene_file == "all_genes") {
  #     labels_x <- c(paste0("-",before), "TSS", "TES", paste0("+",after))
  #   } else {
  #     labels_x <- c(paste0("-",before), "Peak start", "Peak end", paste0("+",after))
  #   }
  #   
  #   ggplot(plot_data, aes(x = Pos, color = condition, fill = condition)) +
  #     geom_ribbon(aes(ymin = score_mean - score_sd, ymax = score_mean + score_sd), alpha = 0.3, color = NA) +
  #     geom_line(aes(y = score_mean), size = 1) +
  #     scale_x_continuous(breaks = c(6,
  #                                   6 + before/binsize,
  #                                   6 + before/binsize + region_body_length/binsize,
  #                                   6 + before/binsize + region_body_length/binsize + after/binsize),
  #                        labels = labels_x) +
  #     geom_vline(xintercept = c(6 + before/binsize, 6 + before/binsize +region_body_length/binsize), linetype = "dashed", color = "grey", alpha = 0.5) +
  #     scale_color_manual(values = custom_colors) +
  #     scale_fill_manual(values = custom_colors) +
  #     labs(x = "Gene regions", 
  #          y = "Average Score", 
  #          title = "Combined Metagene Plot",
  #          color = "Condition",
  #          fill = "Condition") +
  #     theme_classic() +
  #     theme(legend.position = "right",
  #           plot.title = element_text(hjust = 0.5, face = "bold"),
  #           aspect.ratio = 0.5)
  # }
  # 
  # Function to smooth the data by condition ----
  smooth_data_by_condition <- function(data, span = 0.3) {
    # Apply LOESS smoothing to the mean and standard deviation for each condition
    data <- data %>%
      group_by(condition) %>%
      mutate(
        smoothed_mean = predict(loess(score_mean ~ Pos, span = span)),
        smoothed_sd = predict(loess(score_sd ~ Pos, span = span))
      ) %>%
      ungroup()
    
    return(data)
  }
  
  # Function to create the combined metagene plot
  create_combined_metagene_plot <- function(plot_data, smoothing_span = 0.025) {
    debug_print("Creating combined plot with data:", plot_data)
    if (nrow(plot_data) == 0) {
      stop("No data available for plotting. Check previous steps.")
    }
    
    custom_colors <- c("darkgrey", "#1F78B4A6", "#E31A1CA6")
    
    # Construct labels. If gene_file == "all_genes"
    if (gene_file == "all_genes") {
      labels_x <- c(paste0("-", before), "TSS", "TES", paste0("+", after))
    } else {
      labels_x <- c(paste0("-", before), "Peak start", "Peak end", paste0("+", after))
    }
    
    # Smooth the plot data by condition
    smoothed_data <- smooth_data_by_condition(plot_data, span = smoothing_span)
    
    # Plot the smoothed data
    ggplot(smoothed_data, aes(x = Pos, color = condition, fill = condition)) +
      # Smoothed ribbon
      geom_ribbon(aes(ymin = smoothed_mean - smoothed_sd, ymax = smoothed_mean + smoothed_sd), 
                  alpha = 0.3, color = NA) +
      # Smoothed line
      geom_line(aes(y = smoothed_mean), size = 0.7) +
      scale_x_continuous(breaks = c(6,
                                    6 + before / binsize,
                                    6 + before / binsize + region_body_length / binsize,
                                    6 + before / binsize + region_body_length / binsize + after / binsize),
                         labels = labels_x) +
      geom_vline(xintercept = c(6 + before / binsize, 6 + before / binsize + region_body_length / binsize), 
                 linetype = "dashed", color = "grey", alpha = 0.5) +
      scale_color_manual(values = custom_colors) +
      scale_fill_manual(values = custom_colors) +
      labs(x = "Gene regions", 
           y = "Average Score", 
           title = "Combined Metagene Plot",
           color = "Condition",
           fill = "Condition") +
      theme_classic() +
      theme(legend.position = "right",
            plot.title = element_text(hjust = 0.5, face = "bold"),
            aspect.ratio = 0.5)
  }
  
  # Main execution
  meta_name <- paste0(ROOT_FOLDER, FILENAME_TRUNK, "_conditions.csv")
  meta_tracks <- read_metadata(meta_name)
  
  gene_bodies <- prepare_gene_bodies(n_samples = n_genes, gene_file = gene_file, separate_strands = separate_strands)
  
  bed_file_prefix <- paste0(OUTPUT_FOLDER_METAGENE, "/sampled_peaks")
  export_gene_bodies(gene_bodies, bed_file_prefix)
  
  output_metagene_fn <- paste0(OUTPUT_FOLDER_METAGENE, "/", FILENAME_TRUNK, "_metagene_peaks")
  
  if (regenerate_matrices) {
    print("Running deepTools computeMatrix. This may take some time.")
    
    num_cores <- 15
    registerDoParallel(cores = num_cores)
    
    results <- foreach(i = 1:nrow(meta_tracks), .combine = c, .packages = c("data.table")) %dopar% {
      tryCatch({
        for (strand in names(gene_bodies)) {
          bigwig_file <- meta_tracks$bigwig_fn[i]
          regions_bed <- paste0(bed_file_prefix, "_", strand, ".bed")
          output_prefix <- paste0(output_metagene_fn, "_", meta_tracks$sample_id[i])
          
          cmd <- paste(
            "computeMatrix scale-regions",
            "-S", bigwig_file,
            "-R", regions_bed,
            paste0("--quiet -b ", before, " -a ", after, " --binSize ", binsize, " -p 3 --regionBodyLength ", region_body_length),
            "-o", paste0(output_prefix, "_", strand, "_matrix.gz")
          )
          
          cat("\nDEBUG: Running deepTools command:", cmd, "\n")
          system(cmd)
        }
        NULL  # Return NULL if successful
      }, error = function(e) {
        list(error = TRUE, message = as.character(e))
      })
    }
    
    errors <- Filter(function(x) !is.null(x) && inherits(x, "list") && x$error, results)
    if (length(errors) > 0) {
      print("Errors occurred during parallel execution:")
      for (error in errors) {
        print(error$message)
      }
    }
    
    stopImplicitCluster()
  }
  
  df_metagene <- process_samples(meta_tracks, output_metagene_fn)
  
  plot_data <- prepare_plot_data(df_metagene, meta_tracks)
  
  if (nrow(plot_data) > 0) {
    plot <- create_combined_metagene_plot(plot_data)
    ggsave(paste0(OUTPUT_FOLDER_METAGENE, "/", FILENAME_TRUNK, "_metagene_plot.png"), plot, width = 10, height = 4)
    print("Metagene analysis and plotting completed. Check the output files.")
  } else {
    print("Error: No data available for plotting. Check the debug output for more information.")
  }
  
  # Save the sample IDs per condition to CSV for QC purposes
  processed_samples <- unique(df_metagene$sample_id)
  samples_processed <- meta_tracks[sample_id %in% processed_samples, .(sample_id, condition)]
  fwrite(samples_processed, file = paste0(OUTPUT_FOLDER_METAGENE, "/", FILENAME_TRUNK, "_samples_per_condition.csv"))
  
  # Return main dataframes and plot
  list(
    meta_tracks = meta_tracks,
    df_metagene = df_metagene,
    plot_data = plot_data,
    plot = plot
  )
}


# load consensus peak .saf file
merged_peak_file_saf <- file.path("~/5hmC-Pathway-Analysis/Output/DESeq2/Results/significant_results_annotated_enhancer_PM_negativePM_positive_ntop 1000 _ nocombat_10per0702024_seed126.txt")
#merged_peak_file_saf <- file.path(output_dir, "merged_peaks_10perc_pe5_072224_gene.saf")
merged_peak_file_df <- read.table(merged_peak_file_saf, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# if V1 row 1 contains "chr"
if(grepl("chr", merged_peak_file_df$V1[1])){
  # drop first row
  merged_peak_file_df <- merged_peak_file_df[-1,]
  # keep only rows where V5 > 0
  merged_peak_file_df <- merged_peak_file_df[merged_peak_file_df$V5 < 0,]
  # keep only columns V20, V1, V2, V3, V16, V23 in that order
  merged_peak_file_df <- merged_peak_file_df[, c(20, 1, 2, 3, 16, 21)]
  # rename columns to V1, V2, V3, V4, V5
  colnames(merged_peak_file_df) <- c("V1", "V2", "V3", "V4", "V5", "V6")
}

# add a column filled with the 2nd to last column minus the 3rd to last column
# add it second to last
merged_peak_file_df$length <- 0
# output to .bed file with no row names, no column names keeping only columns in this order:
# V2 V3 V4 length V5
write.table(merged_peak_file_df[, c(2, 3, 4, 6, 5)], file = paste0("~/5hmC-Pathway-Analysis/Output/metagene","/sampled_peaks", ".bed"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


results <- metagene_analysis(
  gene_file = "all_genes", # paste0("~/5hmC-Pathway-Analysis/Output/metagene","/sampled_peaks", ".bed"), #"all_genes",
  separate_strands = TRUE,
  regenerate_matrices = FALSE,
  binsize=10,
  n_genes=0,
  before = 3000,
  after = 3000,
  region_body_length = 6000
)

meta_tracks <- results$meta_tracks
df_metagene <- results$df_metagene
plot_data <- results$plot_data
plot <- results$plot
plot

# Save plot to metagene output folder
ggsave(paste0("~/5hmC-Pathway-Analysis/Output/metagene/", "HEALTHY_PM_positive_DESeq2_metagene_plot_10_19_24_3.png"), plot, width = 6, height = 2, units = "in", dpi = 300)
# save plot to svg
ggsave(paste0("~/5hmC-Pathway-Analysis/Output/metagene/", "HEALTHY_PM_positive_DESeq2_metagene_plot_10_19_24_3.svg"), plot, width = 6, height = 2, units = "in")


### Plot enrichment boxplots by region type ----
# set to true if you want to include healthy in heatmap
# define contrast groups

ver <- "8_19_2024"
comparisons <- list(
  c("PM_negative", "PM_positive"),
  c("PM_negative", "HEALTHY"),
  c("PM_positive", "HEALTHY")
)

if(TRUE){
  ### For plotting healthy with PMpos and PM neg in heatmap
  meta_name_all <- "~/5hmC-Pathway-Analysis/Output/Raw Data Processing/HEALTHY_PM_positive_10per_080424_hmr_nocombat_pm_healthy_X/HEALTHY_PM_positive_DESeq2_conditions.csv"
  counts_name_all <- "~/5hmC-Pathway-Analysis/Output/Raw Data Processing/HEALTHY_PM_positive_10per_080424_hmr_nocombat_pm_healthy_X/HEALTHY_PM_positive_DESeq2_rawcounts.csv"
  counts_data_all <- read.csv(counts_name_all,row.names = 1)
  meta_all <-  read.csv(meta_name_all,row.names = 1)
  
  # drop all columns in meta where condition is not equal to healthy or PM_positive
  meta_all <- meta_all[meta_all$condition %in% c("PM_positive", "PM_negative","HEALTHY"),]
  
  # drop columns from counts_data_all that are not in meta_all
  counts_data_all <- counts_data_all[, colnames(counts_data_all) %in% rownames(meta_all)]
  
  #contrast_groups <- unique conditions
  groups_all <- unique(meta_all$condition)
  # contrast_groups <- c()
  # for (i in 1:nrow(groups_all)) {
  #   contrast_groups <- c(contrast_groups, groups_all[i,1])
  # }
  # # drop NA elements
  # contrast_groups <- contrast_groups[!is.na(contrast_groups)]
  
  # Create DESeq2 object without specifying a design
  dds_all <- DESeqDataSetFromMatrix(countData = counts_data_all, 
                                    colData = meta_all, 
                                    design = ~ 1)  # Using '~ 1' as a placeholder design
  
  # Estimate size factors and normalize
  #dds_all <- estimateSizeFactors(dds_all)
  #normalized_counts_all <- counts(dds_all, normalized=TRUE)
  normalized_counts_all <- counts(dds_all, normalized = FALSE) / (colSums(counts(dds_all, normalized = FALSE)) / 1e6)
  
  normalized_counts_all_tb <- data.frame(counts_data_all) %>% 
    rownames_to_column(var="gene") %>% 
    as_tibble()
  
  ### If plotting boxplots for enrichment in different types
  # read in all_results_annotated_enhancer file
  all_results_annotated_enhancer <- read.csv("~/5hmC-Pathway-Analysis/Output/DESeq2/Results/all_results_annotated_PM_positivePM_negative_ntop 1000 _ nocombat_10per081424_seed126.txt",
                                             header = TRUE, sep = "\t")
  # keep only columns called "feature" and "locationType"
  all_results_annotated_enhancer <- all_results_annotated_enhancer[,c("feature","locationType")]
  
  norm_sig_all <- normalized_counts_all_tb %>%
    filter(gene %in% all_results_annotated_enhancer$feature) %>%
    data.frame() %>%
    column_to_rownames(var = "gene")
  
  nsamples_all <- ncol(normalized_counts_all_tb)
} else {
  meta_all <- meta
  norm_sig_all <- norm_sig
}

if(TRUE){
  
  annotation <- meta_all %>%
    dplyr::select(condition, primary_present, ovr_histopath, chemo_6weeks, primary_site)
  
  log_counts_data_mat <- data.matrix(norm_sig_all, rownames.force = NA)
  log_counts_data_mat <- t(scale(t(log_counts_data_mat), scale = FALSE, center = FALSE))
  
  # Prepare data for boxplot
  boxplot_data <- as.data.frame(log_counts_data_mat)
  boxplot_data$gene <- rownames(boxplot_data)
  boxplot_data_long <- tidyr::pivot_longer(boxplot_data, -gene, names_to = "sample", values_to = "normalized_count")
  
  # Add condition information
  boxplot_data_long$condition <- annotation$condition[match(boxplot_data_long$sample, rownames(annotation))]
  
  # Add locationType information
  boxplot_data_long$locationType <- all_results_annotated_enhancer$locationType[match(boxplot_data_long$gene, all_results_annotated_enhancer$feature)]
  
  # Define colors for conditions
  condition_colors <- c("PM_positive" = "#E31A1CA6", "PM_negative" = "#1F78B4A6", "HEALTHY" = "darkgrey")
  
  # Calculate the number of regions in each locationType
  n_regions <- boxplot_data_long %>%
    group_by(locationType) %>%
    summarise(n_regions = n_distinct(gene), .groups = 'drop')
  
  # Function to calculate upper limit (Q3 + 1.5 * IQR)
  calculate_upper_limit <- function(x) {
    q3 <- quantile(x, 0.75)
    iqr <- IQR(x)
    return(q3 + 1.5 * iqr)
  }
  
  # Create a list to store individual plots
  plots <- list()
  
  # Create a plot for each locationType
  for (loc in unique(boxplot_data_long$locationType)) {
    data_subset <- boxplot_data_long %>% filter(locationType == loc)
    
    n <- n_regions %>%
      filter(locationType == loc) %>%
      pull(n_regions)
    
    # Calculate upper limit for y-axis
    y_upper_limit <- calculate_upper_limit(data_subset$normalized_count)
    
    p <- ggplot(data_subset, aes(x = condition, y = normalized_count, fill = condition)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      theme_minimal(base_size = 12) +
      labs(title = paste(loc, "(n =", n, ")"),
           x = "Condition",
           y = "Normalized Count") +
      scale_fill_manual(values = condition_colors) +
      coord_cartesian(ylim = c(0, y_upper_limit)) +  # Set y-axis limit
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
      )
    
    plots[[loc]] <- p
  }
  
  # Combine all plots
  combined_plot <- gridExtra::grid.arrange(grobs = plots, ncol = 3)
  
  # Display the combined plot
  print(combined_plot)
  
  # Save combined plot
  ggsave(paste("./Output/DESeq2/enrichment_boxplots/normalized_count_boxplot_subplots_", paste(contrast_groups, collapse = "_"), "_", ver, ".png", sep = ""), 
         combined_plot, width = 20, height = 12)
  
  ggsave(paste("./Output/DESeq2/enrichment_boxplots/normalized_count_boxplot_subplots_", paste(contrast_groups, collapse = "_"), "_", ver, ".svg", sep = ""), 
         combined_plot, width = 20, height = 12)
  
}