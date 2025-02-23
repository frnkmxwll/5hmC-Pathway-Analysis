### rawcounts_processing.R ###
# Purpose of this script is process raw featureCounts file
# and reduce it to only the required sample in the correct order.

### INSTALL LIBRARIES
# Setup, uncomment follow
# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("DESeq2")
# install.packages("dplyr")
# BiocManager::install("sva")

library(dplyr)
library(tibble)
library(DESeq2)
library(sva)

### CONFIGURATION
# set working directory
setwd("~/5hmC-Pathway-Analysis/")

# settings
combine_files = FALSE
combat_seq_normalize = FALSE

# counts file expected to be in featureCounts default export format
raw_counts_file_1 <- "~/peritoneal_processing/DESeq2/ReadsCount_AutosomeGenes_macs2_pe5_040124.txt"
raw_counts_file_2 <- "./Raw Input/Working Inputs/ReadsCount_AutosomeGenes_ProteinCoding_Pilot_Dataset_added0.txt"

counts1_sample_pathname <- "../../peritoneal_processing/trimmed_data_bam_9_18_2023/"
counts1_sample_extension <- "_bowtie2.bam"
counts2_sample_pathname <- "/media/CLab3b/xiaolong/cfPeri/Bam/"
counts2_sample_extension <- ".bam"

sample_file <- "~/5hmC-Pathway-Analysis/Raw Input/Working Inputs/all_comparisons_6_16_2024.csv"
excluded_samples <- c("KT026","KT027", "KT156")

file_version <- "06162024_hmr_nocombat"
selected_conds <- list("ADE_AMneg","ADE_AMpos")

# Read in data
raw_counts_data_1 <- read.table(
  raw_counts_file_1,
  sep = "\t",
  header = TRUE,
  skip = 1,
  check.names = F
)

# Drop preamble
colnames(raw_counts_data_1) <- gsub("/home/turagapm/peritoneal_processing/trimmed_data_bam_9_18_2023/", "", colnames(raw_counts_data_1))

raw_counts_data_2 <- read.table(
  raw_counts_file_2,
  sep = "\t",
  header = TRUE,
  skip = 1,
  check.names = F
)

sample_data <- read.csv(
  sample_file,
  strip.white = TRUE,
)

### Select relevant samples
# Keep only rows in sample_data where condition is in selected_conds
sample_data <- sample_data[sample_data$condition %in% selected_conds, ]

### COMBINE COUNT FILES (if applicable)
if(combine_files==TRUE){
  if (all(raw_counts_data_1[,1] %in% raw_counts_data_2[,1])==0){
    print("Error: List of genes in your two input counts files are different. Please correct and re-run.")
    stop()
  }
  raw_counts <- inner_join(raw_counts_data_1, raw_counts_data_2)
} else {
  raw_counts <- raw_counts_data_1
  print("Proceeding without combining counts files.")
}

### PROCESS FILE TO DEFINED OUTPUT
counts_trimmed <- subset(raw_counts,select=-c(Chr,Start,End,Strand,Length))
counts_final <- subset(raw_counts,select=-c(Chr,Start,End,Strand,Length))

# Simplify sample names in count file to remove pathname and extensions
colnames(counts_final) <- gsub(counts1_sample_extension, "", colnames(counts_trimmed))
colnames(counts_final) <- gsub(counts1_sample_pathname, "", colnames(counts_final))

if(combine_files==TRUE && counts2_sample_extension != counts1_sample_extension){
  colnames(counts_final) <- gsub(counts2_sample_extension, "", colnames(counts_final))
}

if(combine_files==TRUE && counts2_sample_pathname != counts1_sample_pathname){
  colnames(counts_final) <- gsub(counts2_sample_pathname, "", colnames(counts_final))
}

### VALIDATION
sample_data <- sample_data[! sample_data$X %in% excluded_samples, ]
sample_names <- dplyr::select(sample_data, c(1))
sample_names <- sample_names[sample_names != 'X']

if(all(sample_names %in% colnames(counts_final))==FALSE){
  print("Error: Some samples in your conditions file do not exist in your counts file. Please correct this and re-run.")
  stop()
}

# Eliminate unnecessary samples in count file and sample data file
counts_final_sample <- counts_final[c('Geneid', sample_names)]
counts_final_sample <- counts_final_sample[, !colnames(counts_final_sample) %in% excluded_samples]

### RE-ORDER SAMPLES IN COUNTS FILE BY SAMPLE NAME AND BY CLASS IN ASCENDING ORDER
class_vector <- counts_final_sample[1, ]
class_vector[1, ] <- names(counts_final_sample)
class_vector[] <- sample_data$condition[match(unlist(class_vector), sample_data$X)]
class_vector[1,1] <- "sample_ORDER"

# Sort by sample name first
class_vector_ordered <- class_vector[, c(1, order(colnames(class_vector[, 2:ncol(class_vector)])) + 1)]

# Then sort by sample class
class_vector_ordered <- class_vector_ordered[, c(1, order(unlist(class_vector_ordered[1, 2:ncol(class_vector_ordered)], use.names=FALSE)) + 1)]

counts_final_ordered <- counts_final_sample[, colnames(class_vector_ordered)]

### RE-ORDER SAMPLE DATA FILE BY SAMPLE NAME AND BY CLASS IN ASCENDING ORDER
sample_data_ordered <- sample_data[order(sample_data$condition, sample_data$X), ]
sample_data_ordered <- sample_data_ordered[!sample_data_ordered$X %in% excluded_samples, ]

colnames(sample_data_ordered)[1] <- ""

### Combat-Seq normalization if enabled
if (combat_seq_normalize) {
  counts_final_matrix <- as.matrix(counts_final_sample[-1])
  batch_vector <- as.numeric(factor(sample_data_ordered$batch))
  combat_counts <- ComBat_seq(counts_final_matrix, batch = meta$batch, group = meta$condition)
  counts_final_sample[-1] <- as.data.frame(combat_counts)
}

### CREATE ssGSEA CLASS FILES
sample_count <- nrow(sample_data_ordered)
class1_name <- sample_data_ordered[2, 2]
class2_name <- sample_data_ordered[nrow(sample_data_ordered), 2]
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
counts_final_ordered_rownames <- counts_final_ordered[, -1]
rownames(counts_final_ordered_rownames) <- counts_final_ordered[, 1]

gene_length_rownames <- gene_length
rownames(gene_length_rownames) <- gene_length[, 1]
gene_length_rownames <- within(gene_length_rownames, rm(Geneid))

gene_number <- nrow(counts_final_ordered_rownames)

counts_to_tpm <- function (counts_final_ordered_rownames, gene_length_rownames) {
  x <- counts_final_ordered_rownames / gene_length_rownames
  return (t(t(x) * 1e6 / colSums(x)))
}

tpm <- counts_to_tpm(counts_final_ordered_rownames, gene_length_rownames[, 1])
tpm_dataframe <- as.data.frame(tpm)

description = matrix(c(rep("na", gene_number)), gene_number, 1)
tpm_genepattern <- add_column(tpm_dataframe, description, .before=1)
tpm_genepattern <- rownames_to_column(tpm_genepattern, var = "NAME")

# Create root output folder if it doesn't exist
if (!file.exists("./Output/")) {
  dir.create("./Output/")
}

if (!file.exists("./Output/Raw Data Processing/")) {
  dir.create("./Output/Raw Data Processing/")
}

dir.create(paste("./Output/Raw Data Processing/", class1_name, "_", class2_name, "_", file_version, sep=""))

### OUTPUT RESUTLING FILES

# Output counts .csv file for DESeq2 normalization followed by GSEA -or- counts to TPM conversion for ssGSEA
write.csv(counts_final_ordered, file = paste("./Output/Raw Data Processing/", class1_name, "_", class2_name, "_", file_version, "/", class1_name, "_", class2_name, "_DESeq2_rawcounts.csv", sep = ""), row.names = FALSE, quote = FALSE)

# Output conditions .csv file for DESeq2 normalization
write.csv(sample_data_ordered, file = paste("./Output/Raw Data Processing/", class1_name, "_", class2_name, "_", file_version, "/", class1_name, "_", class2_name, "_DESeq2_conditions.csv", sep = ""), row.names = FALSE, quote = FALSE)

# Output phenotype .cls file for GSEA 
write.table(GSEA_cls, file = paste("./Output/Raw Data Processing/", class1_name, "_", class2_name, "_", file_version, "/", class1_name, "_", class2_name, "_GSEA_phenotype.cls", sep = ""), quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)

# Output normalized counts .txt file for GSEA
write.table(tpm_genepattern, file = paste("./Output/Raw Data Processing/", class1_name, "_", class2_name, "_", file_version, "/", class1_name, "_", class2_name, "_ssGSEA_tpm.tsv", sep = ""), quote = FALSE, sep = '\t', row.names = FALSE)

# DESeq2 normalization without scaling
dds <- DESeqDataSetFromMatrix(countData = counts_final_sample[-1], colData = sample_data_ordered, design = ~ 1)
dds <- DESeq(dds)
normalized_counts <- counts(dds, normalized = TRUE)
normalized_counts_df <- as.data.frame(normalized_counts)
normalized_counts_df <- cbind(Geneid = counts_final_sample$Geneid, normalized_counts_df)

# Output DESeq2 normalized counts .csv file
write.csv(normalized_counts_df, file = paste("./Output/Raw Data Processing/", class1_name, "_", class2_name, "_", file_version, "/", class1_name, "_", class2_name, "_DESeq2_normcounts.csv", sep = ""), row.names = FALSE, quote = FALSE)
