#PURPOSE OF THIS SCRIPT IS TO VISUALIZE HEATMAP
#tutorial taken from here: https://github.com/hbctraining/DGE_workshop/tree/master/lessons

library(DESeq2)
library(magrittr)
library(tibble)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(DEGreport)
library(RColorBrewer)
library(pheatmap)
library(viridis)
library(colorspace)

#set working directory
setwd("C:/Users/ymali/Google Drive/Personal Documents/Chuan Lab/Peritoneal Disease/Data Analysis")

#read in data
counts_data <- read.csv("./Normalization/Input/counts_LGA_v1.csv",row.names = 1)
conditions <-  read.csv("./Normalization/Input/conditions_LGA_v1.csv",row.names = 1)

#check columns are equal
all(colnames(counts_data) %in% rownames(conditions))
all(colnames(counts_data) == rownames(conditions))

#load data into DESeq2 object dds
dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = conditions, design = ~ condition)

#load up size factors into dds object in order to normalize using median of ratios method
dds <- DESeq(dds)

#use counts 'normalized=true' function to pull out normalized counts
normalized_counts <- counts(dds,normalized=TRUE)

#save to file
#write.table(normalized_counts_data, file="normalized_counts_PM_v1.txt", sep="\t", quote=F, col.names=NA)

#define padj cutoff
padj.cutoff <- 0.5

groups <- unique(conditions[c("condition")])

## Define contrasts, extract results table, and shrink the log2 fold changes
contrast <- c("condition", groups[1,1], groups[2,1])
res_table_unshrunken <- results(dds, contrast=contrast, alpha = padj.cutoff)
res_table <- lfcShrink(dds, coef=2, res=res_table_unshrunken)

#display res_table
#res_table %>% data.frame() %>% View()

#Display summary results
# summary(res_table)

#convert the results table into a tibble:
res_table_tb <- res_table %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

#subset table to only keep significant genes using cutoff
sig <- res_table_tb %>% 
  filter(padj < padj.cutoff)

normalized_counts <- normalized_counts %>% 
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

nsamples <- ncol(normalized_counts)

### Extract normalized expression for significant genes from the OE and control samples (4:9), and set the gene column (1) to row names
norm_sig <- normalized_counts[,c(1,2:nsamples)] %>% 
  filter(gene %in% sig$gene) %>% 
  data.frame() %>%
  column_to_rownames(var = "gene") 

### Annotate our heatmap (optional)
annotation <- conditions %>% 
  select(condition)

### save results files
write.table(res_table, file=paste("./DESeq2_Visualization/Results/all_results_",contrast[2],contrast[3],"_v2.txt"), sep="\t", quote=F, col.names=NA)
write.table(sig, file=paste("./DESeq2_Visualization/Results/significant_results_",contrast[2],contrast[3],"_v2.txt"), sep="\t", quote=F, col.names=NA)

### save heatmap

chart_title <- paste(contrast[2],"/",contrast[3],"padj<",padj.cutoff)
chart_title
png(paste("./DESeq2_Visualization/Heatmaps/sig_heatmap_",contrast[2],contrast[3],"_v2.png"), width = 900, height = 1200)
pheatmap(norm_sig, 
         main = chart_title,
         color = diverging_hcl(15,"Blue-Red2"), 
         cluster_rows = T,
         cluster_cols = T,
         show_rownames = T,
         annotation = annotation, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20,
)
dev.off()

### save PCA plots
rld <- vst(dds, blind=TRUE)
sig_genes <-  sig$gene
sig_genes

png(paste("./DESeq2_Visualization/PCA/sig_PCA_",contrast[2],contrast[3],"_v2.png"), width = 900, height = 1200)
plotPCA(
  rld[sig_genes,], 
  intgroup = "condition", 
  )
dev.off()
