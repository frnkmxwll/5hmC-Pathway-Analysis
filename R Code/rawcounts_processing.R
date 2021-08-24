# Purpose of this script is process raw featureCounts file
# and reduce it to only the required sample in the correct order.

### INSTALL LIBRARIES
# Setup, uncomment follow
# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("DESeq2")
# install.packages("dplyr")

library(dplyr) # note: dply is necessary as data frames are too large for R's native merge function, need dplyr's inner join function.
library(DESeq2)

### CONFIGURATION
# set working directory
setwd("G:/My Drive/Graduate School/Chuan Lab/Peritoneal Disease/Raw data preparation")

# settings
combine_files = TRUE #set to TRUE to combine multiple counts files, FALSE if working with a single counts file.

# counts file expected to be in featureCounts default export format
raw_counts_file_1 <- "./Input/ReadsCount_AutosomeGenes_ProteinCoding_wholeGB.txt"
raw_counts_file_2 <- "./Input/ReadsCount_AutosomeGenes_ProteinCoding_wholeGB.txt" #if applicable

# The following variables will simplify the sample names to eliminate any leading pathname and trailing extensions.
# e.g. the feature counts columns sometimes look like "/media/CLab3b/xiaolong/cfPeri/Bam/KT19.bam"
# define these variables for script to trim this text, and simplify column names to "KT19" in above example
counts1_sample_pathname <- "X.media.CLab3b.xiaolong.cfPeri.Bam."
counts1_sample_extension <- ".bam"

counts2_sample_pathname <- "X.media.CLab3b.xiaolong.cfPeri.Bam."
counts2_sample_extension <- ".bam"

# sample file expected to be in here: shorturl.at/qCU19
sample_file <- "./Input/SAMPLE_FILE.csv"

#read in data
raw_counts_data_1 <- read.table(
  raw_counts_file_1,
  sep = "\t",
  header = TRUE,
  skip = 1
  )

raw_counts_data_2 <- read.table(
  raw_counts_file_2,
  sep = "\t",
  header = TRUE,
  skip = 1
  )

sample_data <- read.csv(
  sample_file, 
  strip.white = TRUE,
  )

### COMBINE COUNT FILES (if applicable)


if(combine_files==TRUE){
  if (all(raw_counts_data_1[,1] %in% raw_counts_data_2[,1])==0){
    print("Error: List of genes in your two input counts files are different. Please correct an re-run.")
    stop()
  }
  raw_counts <- inner_join(raw_counts_data_1, raw_counts_data_2)
}

if(combine_files==FALSE){
  raw_counts <- raw_counts_data_1
  print("Proceeding without combining counts files.")
}

### PROCESS FILE TO DEFINED OUTPUT
# Eliminate necessary / unnecessary metadata columns in count file
counts_trimmed <- subset(raw_counts,select=-c(Chr,Start,End,Strand,Length))
counts_final <- subset(raw_counts,select=-c(Chr,Start,End,Strand,Length))

# Simplify sample names in count file

colnames(counts_final) <- gsub(counts1_sample_extension, "", colnames(counts_trimmed))

colnames(counts_final) <- gsub(counts1_sample_pathname, "", colnames(counts_final))

#this runs only if you are combining 2 files with different sample extensions
if(combine_files==TRUE && counts2_sample_extension != counts1_sample_extension){
  colnames(counts_final) <- gsub(counts2_sample_extension, "", colnames(counts_final))
} 

#this runs only if you are combining 2 files with different sample pathnames
if(combine_files==TRUE && counts2_sample_pathname != counts1_sample_pathname){
  colnames(counts_final) <- gsub(counts2_sample_pathname, "", colnames(counts_final))
} 

### VALIDATION
# confirm all sample in sample / conditions table, exist in counts file.
sample_names <- select(sample_data,c(1))
sample_names <- sample_names[sample_names != 'X']


if(all(sample_names %in% colnames(counts_final))==FALSE){
  print("Error: Some sample in your conditions file, do not exist in your counts file. Please correct this and re-run.")
  stop()
}

# Eliminate unnecessary samples in count file
counts_final_sample <- counts_final[c('Geneid',sample_names)]

### RE-ORDER SAMPLE BASED ON CLASS
# Create lookup vector that includes all sample and their class from the sample_data input file
class_vector <- counts_final_sample[1,]
class_vector[1,] <- names(counts_final_sample)
class_vector[] <- sample_data$condition[match(unlist(class_vector), sample_data$X)]
class_vector[1,1] <- "sample_ORDER"

# sort by sample name first; note: the +1 after the order function is to make sure we don't include the "Geneid" column in our sort.
class_vector_ordered <- class_vector[,c(1,order(colnames(class_vector[,2:ncol(class_vector)]))+1)]

# then sort by sample class; note: the +1 after the order function is to make sure we don't include the "Geneid" column in our sort.
class_vector_ordered <- class_vector_ordered [,c(1,order(unlist(class_vector_ordered[1,2:ncol(class_vector_ordered)], use.names=FALSE))+1)]

# then select columns in this order and feed into new dataframe
counts_final_ordered <- counts_final_sample[,colnames(class_vector_ordered)]

### RE-ORDER SAMPLE DATA FILE
sample_data_ordered <- sample_data[order(sample_data$condition,sample_data$X),]

#when important sample_data, the blank header over the sample names got filled in with an "X" this needs to be removed for GSEA.
colnames(sample_data_ordered)[1] <- ""

### CREATE ssGSEA CLASS FILES
sample_count <- nrow(sample_data_ordered)
class1_name <- sample_data_ordered[2,2]
class2_name <- sample_data_ordered[nrow(sample_data_ordered),2]
class1_count <- sum(sample_data_ordered$condition == class1_name)
class2_count <- sum(sample_data_ordered$condition == class2_name)

ssGSEA_cls_list <- list(row1=c(sample_count , "2" , "1"),row2=c("#" , class1_name , class2_name), row3=c(rep(0,class1_count),rep(1,class2_count)))

ssGSEA_cls <- as.data.frame(do.call(rbind, ssGSEA_cls_list))
rownames(ssGSEA_cls) <- c()
names(ssGSEA_cls) <- NULL
ssGSEA_cls [1,4:ncol(ssGSEA_cls)] <- c(rep("",ncol(ssGSEA_cls)-3))
ssGSEA_cls [2,4:ncol(ssGSEA_cls)] <- c(rep("",ncol(ssGSEA_cls)-3))

# ssGSEA needs a gene length file for TPM normalization
gene_length <- raw_counts_data_1[,c('Geneid','Length')]

### CREATE GSEA CLASS FILES
GSEA_cls_list <- list(row1=c(sample_count , "2" , "1"),row2=c("#" , class1_name , class2_name), row3=c(rep(class1_name,class1_count),rep(class2_name,class2_count)))

GSEA_cls <- as.data.frame(do.call(rbind, GSEA_cls_list))
rownames(GSEA_cls) <- c()
names(GSEA_cls) <- NULL
GSEA_cls [1,4:ncol(GSEA_cls)] <- c(rep("",ncol(GSEA_cls)-3))
GSEA_cls [2,4:ncol(GSEA_cls)] <- c(rep("",ncol(GSEA_cls)-3))

### CREATE TPM FILE FOR ssGSEA
# Need tables with row names for tpm conversion, clunky but gets job done
counts_final_ordered_rownames <- counts_final_ordered[,-1]
rownames(counts_final_ordered_rownames) <- counts_final_ordered[,1]

gene_length_rownames <- gene_length
rownames(gene_length_rownames) <- gene_length[,1]
gene_length_rownames <- within(gene_length_rownames, rm(Geneid))

#define number of genes
gene_number <- nrow(counts_final_ordered_rownames)

#calculation to convert counts to tpm

counts_to_tpm <- function (counts_final_ordered_rownames,gene_length_rownames) {
  x <- counts_final_ordered_rownames / gene_length_rownames
  return (t(t(x)*1e6/colSums(x)))
}

tpm <- counts_to_tpm (counts_final_ordered_rownames, gene_length_rownames[,1])

tpm_dataframe <- as.data.frame(tpm)

# add description column after Geneid populated with 'na' required format by ssGSEA
description = matrix(c(rep("na",gene_number)),gene_number,1)
tpm_genepattern <- add_column(tpm_dataframe,description, .before=1)
tpm_genepattern <- rownames_to_column(tpm_genepattern, var = "NAME")

### CREATE GSEA NORMALIZED COUNTS FILE

#load data into DESeq2 object dds
dds <- DESeqDataSetFromMatrix(countData = counts_final_ordered_rownames, colData = sample_data_ordered, design = ~ condition)

#load up size factors into dds object in order to normalize using median of ratios method
dds <- DESeq(dds)

#use counts 'normalized=true' function to pull out normalized counts
normalized_counts <- counts(dds,normalized=TRUE)

##save normalized counts to file. Un-comment this line if you need a normalized counts file to be used in GSEA

normalized_counts_GSEA <- as.data.frame(normalized_counts)

description = matrix(c(rep("na",gene_number)),gene_number,1)

normalized_counts_GSEA <- add_column(normalized_counts_GSEA,description, .before=1)

normalized_counts_GSEA <- rownames_to_column(normalized_counts_GSEA, var = "NAME")

### OUTPUT RESUTLING FILES
# output counts .csv file for DESeq2 normalization followed by GSEA -or- counts to TPM conversion for ssGSEA
write.csv(counts_final_ordered,file = paste("./Output/",class1_name,"_",class2_name,"_DESeq2_rawcounts.csv",sep = ""),row.names = FALSE,quote=FALSE)

# output conditions .csv file for DESeq2 normalization
write.csv(sample_data_ordered,file = paste("./Output/",class1_name,"_",class2_name,"_DESeq2_conditions.csv",sep = ""),row.names = FALSE,quote=FALSE)

# output phenotype .cls file for GSEA 
write.table(GSEA_cls, file=paste("./Output/",class1_name,"_",class2_name,"_GSEA_phenotype.cls",sep = ""), quote=FALSE, sep='\t', row.names = FALSE, col.names = FALSE)

# output normalized counts .txt file for GSEA
write.table(normalized_counts_GSEA, file=paste("./Output/",class1_name,"_",class2_name,"_GSEA_normcounts.txt", sep = ""), sep="\t", quote=F, row.names=FALSE)

# output phenotype .cls file for ssGSEA 
write.table(ssGSEA_cls, file=paste("./Output/",class1_name,"_",class2_name,"_ssGSEA_phenotype.cls",sep = ""), quote=FALSE, sep='\t', row.names = FALSE, col.names = FALSE)

# output gene length .csv file for ssGSEA counts to TPM normalization, no longer required as tpm conversion has been built in here.
# write.csv(gene_length,file = "./Output/gene_length.csv",row.names = FALSE,quote=FALSE)

#output tpm tsv file
write.table(tpm_genepattern,file=paste("./Output/",class1_name,"_",class2_name,"_ssGSEA_tpm.tsv",sep = ""), quote=FALSE, sep='\t', row.names=FALSE)
