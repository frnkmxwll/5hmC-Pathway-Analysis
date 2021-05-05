#PURPOSE OF THIS SCRIPT IS TO VISUALIZE DGE HEATMAPS & PCA PLOTS BETWEEN TWO COMPARISON GROUPS
#tutorial taken from here: https://github.com/hbctraining/DGE_workshop/tree/master/lessons

###LIBRARIES INSTALL ALL IF FIRST TIME RUNNING###
#install.packages("package_name")
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
library(M3C)

###CONFIGURATION
#set working directory, select where you extracted folder
setwd("C:/Users/ymali/Google Drive/Personal Documents/Chuan Lab/Peritoneal Disease/Data Analysis")

#read in data, define what counts & conditions files
counts_data <- read.csv("./Input/counts_pm_v1.csv",row.names = 1)
conditions <-  read.csv("./Input/conditions_pm_v1.csv",row.names = 1)

#define padj cutoff, you may need to run with several padj values until you have an appropriate number of significant results.
padj.cutoff <- 0.05

#Select version for all output files (e.g. 1, 2, 3, ...)
ver <- 4

###VALIDATION
#check columns are equal
all(colnames(counts_data) %in% rownames(conditions))
all(colnames(counts_data) == rownames(conditions))

###CREATE DESeq2 OBJECT
#load data into DESeq2 object dds
dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = conditions, design = ~ condition)

#load up size factors into dds object in order to normalize using median of ratios method
dds <- DESeq(dds)

#use counts 'normalized=true' function to pull out normalized counts
normalized_counts <- counts(dds,normalized=TRUE)

groups <- unique(conditions[c("condition")])

#Define contrasts, extract results table, and shrink the log2 fold changes
contrast <- c("condition", groups[1,1], groups[2,1])
res_table_unshrunken <- results(dds, contrast=contrast, alpha = padj.cutoff)
res_table <- lfcShrink(dds, coef=2, res=res_table_unshrunken)

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

#determine number of samples
nsamples <- ncol(normalized_counts)

#Extract normalized expression for significant genes from the two comparison groups and set the gene column (1) to row names
norm_sig <- normalized_counts[,c(1,2:nsamples)] %>% 
  filter(gene %in% sig$gene) %>% 
  data.frame() %>%
  column_to_rownames(var = "gene") 


### DISPLAY KEY TABLES (uncomment specific rows for debugging)
#norm_sig
#sig
#res_table
#res_table %>% data.frame() %>% View()
#summary(res_table)


### SAVE RESULTS TABLES TO TEXT FILES ###
write.table(res_table, file=paste("./DESeq2_Visualization/Results/all_results_",contrast[2],contrast[3],"_v",ver,".txt", sep = ""), sep="\t", quote=F, col.names=NA)
write.table(sig, file=paste("./DESeq2_Visualization/Results/significant_results_",contrast[2],contrast[3],"_v",ver,".txt", sep = ""), sep="\t", quote=F, col.names=NA)
### save normalized counts to file. Un-comment this line if you need a normalized counts file to be used in GSEA
#write.table(normalized_counts_data, file=paste("normalized_counts_PM_",ver,".txt", sep = ""), sep="\t", quote=F, col.names=NA)

### GENERATE HEATMAP ###
#Annotate our heatmap (optional)
annotation <- conditions %>% 
  select(condition)

#Save heatmap to png
chart_title <- paste(contrast[2],"/",contrast[3],"padj<",padj.cutoff)
chart_title
png(paste("./DESeq2_Visualization/Heatmaps/sig_heatmap_",contrast[2],contrast[3],"_v",ver,".png", sep = ""), width = 900, height = 1200)
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

###GENERATE PCA PLOT

#See details for below operation in lesson 3 of DGE workshop
rld <- vst(dds, blind=TRUE)
res_genes <- row.names(res_table)
sig_genes <-  sig$gene

#save PCA plot to png
#In the below replace sig_genes with res_genes if you want to perform PCA analysis on all genes rather than just on significant genes.
png(paste("./DESeq2_Visualization/PCA/sig_PCA_",contrast[2],contrast[3],"_v",ver,".png", sep = ""), width = 900, height = 1200)
plotPCA(
  rld[sig_genes,], 
  intgroup = "condition"
  )
dev.off()


### GENERATE UMAP PLOT
#save UMAP plot to png
png(paste("./DESeq2_Visualization/UMAP/sig_UMAP_",contrast[2],contrast[3],"_v",ver,".png", sep = ""), width = 900, height = 1200)
umap(norm_sig, labels=as.factor(conditions$condition),printres = FALSE, seed = FALSE,
     axistextsize = 18, legendtextsize = 18, dotsize = 5,
     textlabelsize = 4, legendtitle = "Group", controlscale = FALSE,
     scale = 1,     printwidth = 22, text = FALSE)
dev.off()
