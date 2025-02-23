#Purpose of this script is to visualize 5hmc by genetic feature

### INSTALL LIBRARIES
# Setup, uncomment follow
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("CGPfunctions")
#install.packages("DEGreport") 
#install.packages("lasso2")
#BiocManager::install("DEGreport")
# install.packages("svglite")
# BiocManager::install("CGPfunctions", force = TRUE)
#remove.packages(lib="DEGreport")
library(BiocManager)
library(GenomicAlignments)
library(Rsamtools)
library(GenomicFeatures)
library(dplyr)
library(ggplot2)

###CONFIGURATION
#set working directory, select where you extracted folder
setwd("~/5hmC-Pathway-Analysis/")

##whole dataset balanced
#counts_name <- "./Output/Raw Data Processing/1-2-3-4PMpos_8-9-11PMneg_vFINAL/1-2-3-4PMpos_8-9-11PMneg_DESeq2_rawcounts.csv"
#meta_name <- "./Output/Raw Data Processing/1-2-3-4PMpos_8-9-11PMneg_vFINAL/1-2-3-4PMpos_8-9-11PMneg_DESeq2_conditions.csv"
#validation
#counts_name <- "./Output/Randomization/1-2-3-4PMpos_8-9-11PMneg_DESeq2_whole_balance/1-2-3-4PMpos_8-9-11PMneg_validation_rawcounts.csv"
#meta_name <- "./Output/Randomization/1-2-3-4PMpos_8-9-11PMneg_DESeq2_whole_balance/1-2-3-4PMpos_8-9-11PMneg_validation_conditions.csv"
#training
#counts_name <- "./Output/Randomization/1-2-3-4PMpos_8-9-11PMneg_DESeq2_whole_balance/1-2-3-4PMpos_8-9-11PMneg_training_rawcounts.csv"
#meta_name <- "./Output/Randomization/1-2-3-4PMpos_8-9-11PMneg_DESeq2_whole_balance/1-2-3-4PMpos_8-9-11PMneg_training_conditions.csv"

##CRC-HGA vs nonCANCERcontrols
meta_name <- "./Output/Raw Data Processing/CRC_HGA_NONcancerCONTROLS_CRC-HGA_nonCANCER/CRC_HGA_NONcancerCONTROLS_DESeq2_conditions.csv"
counts_name <- "./Output/Raw Data Processing/CRC_HGA_NONcancerCONTROLS_CRC-HGA_nonCANCER/CRC_HGA_NONcancerCONTROLS_DESeq2_rawcounts.csv"

#read in data, define what counts & conditions files
counts_data <- read.csv(counts_name,row.names = 1)
meta <-  read.csv(meta_name,row.names = 1)

# Assuming you have a GTF file for your genome
txdb <- makeTxDbFromGFF(file="~/reference_genomes/gencode.v43.annotation.gtf.gz")

# Define the regions of interest
promoters <- promoters(txdb, upstream=2000, downstream=500)
three_utr <- threeUTRsByTranscript(txdb, use.names=TRUE)
five_utr <- fiveUTRsByTranscript(txdb, use.names=TRUE)
exons <- exonsBy(txdb, by="tx", use.names=TRUE)
introns <- intronsByTranscript(txdb, use.names=TRUE)

# Get list of bam files
bam_files <- list.files(path = "your_directory", pattern = "*.bam", full.names = TRUE)

# Initialize dataframe for feature counts
feature_counts <- data.frame()

# Iterate over each bam file, extracting counts for each feature
for(file in bam_files) {
  sample_name <- gsub(".bam", "", basename(file))
  
  # Load Bam file
  param <- ScanBamParam(which=promoters)
  reads <- readGAlignments(file, param=param)
  
  # Extract counts for each region and bind them into a dataframe
  count_data <- data.frame(
    sample = sample_name,
    promoter_count = sum(countOverlaps(promoters, reads)),
    three_utr_count = sum(countOverlaps(three_utr, reads)),
    five_utr_count = sum(countOverlaps(five_utr, reads)),
    exon_count = sum(countOverlaps(exons, reads)),
    intron_count = sum(countOverlaps(introns, reads))
  )
  
  feature_counts <- rbind(feature_counts, count_data)
}

# Read in metadata and merge it with the counts data
metadata <- read.csv("conditions.csv")
data <- left_join(feature_counts, metadata, by = "sample")

# Melt the data for easier plotting
data_melt <- reshape2::melt(data, id.vars = c("sample", "condition"))

# Plot the data
ggplot(data_melt, aes(x=variable, y=value, fill=condition)) +
  geom_boxplot() +
  labs(x="Feature", y="Count") +
  theme_bw()
