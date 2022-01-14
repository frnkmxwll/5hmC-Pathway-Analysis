#Purpose of this script is to visualize differential gene expression heatmaps,
#PCA plot and other charts between two comparison groups.

### INSTALL LIBRARIES
# Setup, uncomment follow
# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("DESeq2")
# install.packages("dplyr")
# BiocManager::install("sva")

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
library(ashr)

###CONFIGURATION
#set working directory, select where you extracted folder
setwd("~/5hmC-Pathway-Analysis/")

counts_name <- "./Output/Raw Data Processing/1or3pdONLYpos_8or11metNEGtumorNEG_whole_combatseq/1or3pdONLYpos_8or11metNEGtumorNEG_DESeq2_rawcounts.csv"
meta_name <- "./Output/Raw Data Processing/1or3pdONLYpos_8or11metNEGtumorNEG_whole_combatseq/1or3pdONLYpos_8or11metNEGtumorNEG_DESeq2_conditions_bin.csv"

#read in data, define what counts & conditions files
counts_data <- read.csv(counts_name,row.names = 1)
meta <-  read.csv(meta_name,row.names = 1)

#define padj cutoff, you may need to run with several padj values until you have an appropriate number of significant results.
#used to select significant genes for results tables, PCA plots, heatmaps and UMAP plots.
cutoff_type = 1 # 0=padj cutoff, default; 1=lfc & pvalue cutoff
padj.cutoff = 0.1 # 0.1 default
pvalue.cutoff = 0.01
lfc.cutoff = 0.15 #0.3 ~ 20% change, 0.15 ~ 10%

#Select version for all output files (e.g. 1, 2, 3, ...)

ver <- "1o3vs8o11_p0p01lfc0p1batchonly"
gene_number <- nrow(counts_data)

###VALIDATION
#check columns are equal
all(colnames(counts_data) %in% rownames(meta))
all(colnames(counts_data) == rownames(meta))

###CREATE DESeq2 OBJECT
#load data into DESeq2 object dds
design_formula <- ~ condition + batch# + sex + age + race
design_formula

dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = meta, design = design_formula)

#load up size factors into dds object in order to normalize using median of ratios method
dds <- DESeq(dds)

#use counts 'normalized=true' function to pull out normalized counts
normalized_counts <- counts(dds,normalized=TRUE)

groups <- unique(meta[c("condition")])

#Define contrasts, extract results table, and shrink the log2 fold changes
contrast_groups <- c("condition",groups[1,1], groups[2,1])

res_table <- results(dds, contrast=contrast_groups, alpha = padj.cutoff)

#Consider replacing above line with shrunken values for fold change below:
#res_table_unshrunken <- results(dds, contrast=contrast_groups, alpha = padj.cutoff)
#res_table <- lfcShrink(dds, coef = 2, res=res_table_unshrunken)

metadata(res_table)$filterThreshold

plot(metadata(res_table)$filterNumRej, 
     type="b", ylab="number of rejections",
     xlab="quantiles of filter")
lines(metadata(res_table)$lo.fit, col="red")
abline(v=metadata(res_table)$filterTheta)
dev.off()

#convert the results table into a tibble:
res_table_tb <- res_table %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

#subset table to only keep significant genes using cutoff
if(cutoff_type == 0){
  sig <- res_table_tb %>% 
    filter(padj < padj.cutoff)
}

if(cutoff_type == 1){
  sig <- res_table_tb %>% 
    filter(abs(log2FoldChange) > lfc.cutoff, pvalue < pvalue.cutoff)
}

normalized_counts_tb <- normalized_counts %>% 
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

#determine number of samples and significant genes
nsamples <- ncol(normalized_counts_tb)
nsig <- nrow(sig)

#Extract normalized expression for significant genes from the two comparison groups and set the gene column (1) to row names
norm_sig <- normalized_counts_tb[,c(1,2:nsamples)] %>% 
  filter(gene %in% sig$gene) %>% 
  data.frame() %>%
  column_to_rownames(var = "gene") 

###DISPLAY KEY TABLES (uncomment specific rows for debugging)
#normalized_counts_tb
#norm_sig
#sig
#res_table
#res_table %>% data.frame() %>% View()
#summary(res_table)

### CREATE ROOT OUTPUT FOLDERS IF THEY DON'T EXIST
if (file.exists("./Output/")) {
  cat("The folder already exists")
} else {
  dir.create("./Output/")
}

if (file.exists("./Output/DESeq2/")) {
  cat("The folder already exists")
} else {
  dir.create("./Output/DESeq2/")
}

###SAVE RESULTS TABLES TO TEXT FILES
if (file.exists("./Output/DESeq2/Results/")) {
  cat("The folder already exists")
} else {
  dir.create("./Output/DESeq2/Results/")
}

write.table(res_table, file=paste("./Output/DESeq2/Results/all_results_",contrast_groups[2],contrast_groups[3],"_",ver,".txt", sep = ""), sep="\t", quote=F, col.names=NA)
write.table(sig, file=paste("./Output/DESeq2/Results/significant_results_",contrast_groups[2],contrast_groups[3],"_",ver,".txt", sep = ""), sep="\t", quote=F, col.names=NA)


###GENERATE VOLCANO PLOTS
if (file.exists("./Output/DESeq2/Volcano/")) {
  cat("The folder already exists")
} else {
  dir.create("./Output/DESeq2/Volcano/")
}

#Obtain logical vector where TRUE values denote padj values < cutoff and fold change > 1.5 in either direction

#working with log2 fold changes so this translates to an actual fold change of 2^lfc.cutoff. used in volcano plots
res_table_tb_volcano <- res_table_tb %>% 
  mutate(threshold_OE = padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff)

res_table_tb_volcano <- res_table_tb_volcano %>% arrange(padj) %>% mutate(genelabels = "")

res_table_tb_volcano$genelabels[1:nsig] <- res_table_tb_volcano$gene[1:nsig]

png(paste("./Output/DESeq2/Volcano/volcano_",contrast_groups[2],contrast_groups[3],"_",ver,".png", sep = ""), width = 900, height = 900)
ggplot(res_table_tb_volcano, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = threshold_OE)) +
  ggtitle(paste(contrast_groups[2],"/",contrast_groups[3],"Enrichment")) +
  geom_text_repel(aes(label = genelabels,size=20)) +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) +
  xlim(-0.6, 0.6)
dev.off()

###GENERATE HEATMAP
if (file.exists("./Output/DESeq2/Heatmaps/")) {
  cat("The folder already exists")
} else {
  dir.create("./Output/DESeq2/Heatmaps/")
}

#Annotate our heatmap (optional)
annotation <- meta %>% 
  select(condition)

#Save heatmap to png
if(cutoff_type == 0){
heatmap_title <- paste(contrast_groups[2],"/",contrast_groups[3],"padj <",padj.cutoff)
}

if(cutoff_type == 1){
  heatmap_title <- paste(contrast_groups[2],"/",contrast_groups[3],"pvalue <",pvalue.cutoff,"|lfc|>",lfc.cutoff)
}

png(paste("./Output/DESeq2/Heatmaps/sig_heatmap_",contrast_groups[2],contrast_groups[3],"_",ver,".png", sep = ""), width = 900, height = 1200)
pheatmap(norm_sig, 
         main = heatmap_title,
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
if (file.exists("./Output/DESeq2/PCA/")) {
  cat("The folder already exists")
} else {
  dir.create("./Output/PCA/Results/")
}

#See details for below operation in lesson 3 of DGE workshop
rld <- vst(dds, blind=TRUE)
res_genes <- row.names(res_table)
sig_genes <-  sig$gene

#save PCA plot to png
#In the below replace sig_genes with res_genes if you want to perform PCA analysis on all genes rather than just on significant genes.
png(paste("./Output/DESeq2/PCA/sig_PCA_",contrast_groups[2],contrast_groups[3],"_",ver,".png", sep = ""), width = 900, height = 1200)

plotPCA_labels <- plotPCA(
  rld[sig_genes,],
  #rownames(c(row.names(meta)))
  )

#reverse order of labels to match heatmap labelling
plotPCA_labels$data$group<-factor(plotPCA_labels$data$group,levels = rev(levels(plotPCA_labels$data$group)))

plotPCA_labels + geom_text(aes(label = name),position=position_nudge(y = 0.07),) + ggtitle(heatmap_title)

dev.off()


###GENERATE UMAP PLOT
if (file.exists("./Output/DESeq2/UMAP/")) {
  cat("The folder already exists")
} else {
  dir.create("./Output/DESeq2/UMAP/")
}

#save UMAP plot to png
png(paste("./Output/DESeq2/UMAP/sig_UMAP_",contrast_groups[2],contrast_groups[3],"_",ver,".png", sep = ""), width = 900, height = 1200)
umap(norm_sig, labels=as.factor(meta$condition),printres = FALSE, seed = FALSE,
     axistextsize = 18, legendtextsize = 18, dotsize = 5,
     textlabelsize = 4, legendtitle = "Group", controlscale = FALSE,
     scale = 1,     printwidth = 22, text = FALSE)
dev.off()


###SAVE CONFIG TABLES TO TEXT FILES
if (file.exists("./Output/DESeq2/Config/")) {
  cat("The folder already exists")
} else {
  dir.create("./Output/DESeq2/Config/")
}

config <- c(
  paste("counts file name:", counts_name), 
  paste("conditions file name:", meta_name), 
  paste("padj cut off",padj.cutoff),
  paste("output file name:", ver),
  paste("volcano lfc cutoff:", lfc.cutoff),
  paste("cutoff type (0=padj cutoff, default; 1=lfc & pvalue cutoff):", cutoff_type),
  paste("padj cutoff (if type=0):", padj.cutoff),
  paste("pvalue cutoff (if type=1):", pvalue.cutoff),
  paste("lfc cutoff (if type=1):", lfc.cutoff),
  paste("design formula:", c(design_formula))
  )

config_frame <- config
write.table(config, file=paste("./Output/DESeq2/Config/config_",contrast_groups[2],contrast_groups[3],"_",ver,".txt", sep = ""), sep="\t", quote=F, col.names=NA)
