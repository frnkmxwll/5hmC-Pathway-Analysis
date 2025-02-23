#Purpose of this script is to visualize differential gene expression heatmaps,
#PCA plot and other charts between two comparison groups.

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
library(DESeq2) #for DESeq2 analysis
library(pROC)
library(tibble)
library(dplyr)
library(tidyverse)
library(ggplot2) #for plotting
library(ggrepel)
#library(DEGreport)
library(RColorBrewer)
library(pheatmap) #for heatmaps
library(viridis)
library(colorspace)
#library(M3C) #for UMAP
library(ashr)
library(CGPfunctions) #needed for plotxtabs2 
library(ggpubr) # needed for whisker plot charting
library(aod) #for logistic regression
library(broom)
library(leaps) #for logistic regression optimizaton
library(devtools) #used for complex heatmaps
library(ComplexHeatmap)#used for complex heatmaps
library(circlize)

###CONFIGURATION
#set working directory, select where you extracted folder
setwd("~/5hmC-Pathway-Analysis/")

##whole dataset
counts_name <- "./Output/Raw Data Processing/CRCpmonlyPOS_HEALTHY_whole_combatseq/CRCpmonlyPOS_HEALTHY_DESeq2_rawcounts.csv"
meta_name <- "./Output/Raw Data Processing/CRCpmonlyPOS_HEALTHY_whole_combatseq/CRCpmonlyPOS_HEALTHY_DESeq2_conditions.csv"

#define padj cutoff, you may need to run with several padj values until you have an appropriate number of significant results.
#used to select significant genes for results tables, PCA plots, heatmaps and UMAP plots.
cutoff_type = 1 # 0=padj cutoff, default; 1=lfc & pvalue cutoff
padj.cutoff = 0.1 # 0.1 default
pvalue.cutoff = 0.05
lfc.cutoff = 0.137504 # 0.137504 ~ 10%, 0.263034 ~ 20% change, 0.32 ~ 25%, 0.584963 ~ 50% change, 

#Select version for all output files (e.g. 1, 2, 3, ...)
ver <- "CRC_Healthy"

#read in data, define what counts & conditions files
counts_data <- read.csv(counts_name,row.names = 1)
meta <-  read.csv(meta_name,row.names = 1)
gene_number <- nrow(counts_data)

# define contrast groups
groups <- unique(meta[c("condition")])
contrast_groups <- c("condition",groups[1,1], groups[2,1])

###VALIDATION
#check columns are equal
all(colnames(counts_data) %in% rownames(meta))
all(colnames(counts_data) == rownames(meta))

###CREATE DESeq2 OBJECT
#load data into DESeq2 object dds
design_formula <- ~ condition
#+ sex + age + race + batch
#design_formula <- ~ condition + batch
dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = meta, design = design_formula)

#load up size factors into dds object in order to normalize using median of ratios method
dds <- DESeq(dds)

keep <- rowSums( counts(dds) >= 1 ) >= 2
#dds <- dds[keep,]

#use counts 'normalized=true' function to pull out normalized counts
#normalized_counts <- counts(dds,normalized=TRUE)
#normalize_counts <- rlog(dds, blind = TRUE)
vst_data <- vst(dds, blind = TRUE)
normalized_counts <- assay(vst_data)

normalized_counts <- data.frame(t(normalized_counts))

# Suppose your coefficients are stored in a named vector like this
# coefs <- c(gene1 = 0.5, gene2 = -0.2, ...)
coefs <- c(
  PHLDA3 = 0.58834405,
  SNTG2 = 0.079772,
  FBXL7 = 0.52887915,
  FST = 1.37235249,
  PERP = 4.11269898,
  SULF1 = 2.62835328,
  RUNX1T1 = 1.44801754,
  ADAM20 = 0.07389831,
  NME3 = -0.24957724,
  SPAG4 = 0.51709632
)

# subset df_counts
normalized_counts_subset <- normalized_counts[, names(coefs)]

# check if all the required genes are present
if (!all(names(coefs) %in% colnames(normalized_counts_subset))) {
  stop("Not all required genes are present in df_counts")
}

# multiply counts by coefficients and sum
prediction <- rowSums(sweep(normalized_counts_subset, 2, coefs, "*"))

# add the prediction to df_meta
meta$prediction <- prediction
prediction

