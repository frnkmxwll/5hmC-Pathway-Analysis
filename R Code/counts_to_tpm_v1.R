#Per this documentation, ssGSEA requires TPM and not counts.
#https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Using_RNA-seq_Datasets_with_GSEA
#This blog post proposed solution for converting counts to TPM equivalent.
#https://www.biostars.org/p/335187/

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

write.table(tpm,file='genebodies_pm_TPM_v2.tsv', quote=FALSE, sep='\t', col.names = NA)