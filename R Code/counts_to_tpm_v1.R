
setwd("C:/Users/ymali/Google Drive/Personal Documents/Chuan Lab/Peritoneal Disease/Data Analysis/ssGSEA/Inputs/Normalization")


counts_name <- "./genebodies_pm_raw_counts_v1.csv"
length_name <- "./genebodies_length.csv"
#read in data, define what counts & conditions files

counts_data <- read.csv(counts_name,row.names = 1)
length_data <- read.csv(length_name,row.names = 1)

counts_to_tpm <- function (counts_data,length_data) {
  x <- counts_data / length_data
  return (t(t(x)*1e6/colSums(x)))
}

tpm <- counts_to_tpm (counts_data, length_data[,1])
tpm

counts_data

write.csv(tpm,'genebodies_pm_TPM_v1.csv')