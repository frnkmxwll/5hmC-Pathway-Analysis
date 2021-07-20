#Per this documentation, ssGSEA requires TPM and not counts.
#https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Using_RNA-seq_Datasets_with_GSEA
#This blog post proposed solution for converting counts to TPM equivalent.
#https://www.biostars.org/p/335187/

library(tibble)

setwd("C:/Users/ymali/Google Drive/Personal Documents/Chuan Lab/Peritoneal Disease/Data Analysis/ssGSEA/Inputs/Normalization/survival jul 15/")


counts_name <- "./genebodies_PROG_raw_counts_v2.csv"
length_name <- "./genebodies_length.csv"
gene_number <- 19100

#read in data, define what counts & conditions files

counts_data <- read.csv(counts_name,row.names = 1)
length_data <- read.csv(length_name,row.names = 1)

counts_to_tpm <- function (counts_data,length_data) {
  x <- counts_data / length_data
  return (t(t(x)*1e6/colSums(x)))
}

tpm <- counts_to_tpm (counts_data, length_data[,1])
tpm

tpm_dataframe <- as.data.frame(tpm)

description = matrix(c(rep("na",gene_number)),gene_number,1)
description

tpm_genepattern <- add_column(tpm_dataframe,description, .before=1)
tpm_genepattern

tpm_genepattern <-rownames_to_column(tpm_genepattern, var = "NAME")

#counts_data

write.table(tpm_genepattern,file='genebodies_PROG_TPM_v2.tsv', quote=FALSE, sep='\t', row.names=FALSE)