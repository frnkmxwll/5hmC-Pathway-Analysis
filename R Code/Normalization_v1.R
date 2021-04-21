#PURPOSE OF THIS FILE IS TO SAVE NORMALIZED READCOUNTS TO A .txt FILE
#make sure DESeq2 package is installed available
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("DESeq2", quietly = TRUE))
  BiocManager::install("DESeq2")

library(DESeq2)

#set working directory
setwd("C:/Users/ymali/Google Drive/Personal Documents/Chuan Lab/Peritoneal Disease/Data Analysis/Normalization")

#read in data
counts_data <- read.csv("./Input/counts_LGA_v1.csv",row.names = 1)

conditions <-  read.csv("./Input/conditions_LGA_v1.csv",row.names = 1)

#check columns are equal
all(colnames(counts_data) %in% rownames(conditions))
all(colnames(counts_data) == rownames(conditions))

#load data into DESeq2 object dds
dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = conditions, design = ~ condition)

#load up size factors into dds object in order to normalize using median of ratios method
dds <- estimateSizeFactors(dds)

#use counts 'normalized=true' function to pull out normalized counts
normalized_counts_data <- counts(dds,normalized=TRUE)

#save to file
write.table(normalized_counts_data, file="./Output/normalized_counts_LGA_v1.txt", sep="\t", quote=F, col.names=NA)

