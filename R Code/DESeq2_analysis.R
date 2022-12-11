#Purpose of this script is to visualize differential gene expression heatmaps,
#PCA plot and other charts between two comparison groups.

### INSTALL LIBRARIES
# Setup, uncomment follow
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("CGPfunctions")
#install.packages("DEGreport") 
#install.packages("lasso2")
#BiocManager::install("DEGreport")
 install.packages("svglite")
# BiocManager::install("CGPfunctions", force = TRUE)
#remove.packages(lib="DEGreport")
library(BiocManager)
library(DESeq2) #for DESeq2 analysis
library(magrittr)
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

### Used for distance clustering
robust_dist = function(x, y) {
  qx = quantile(x, c(0.1, 0.9))
  qy = quantile(y, c(0.1, 0.9))
  l = x > qx[1] & x < qx[2] & y > qy[1] & y < qy[2]
  x = x[l]
  y = y[l]
  sqrt(sum((x - y)^2))
}


###CONFIGURATION
#set working directory, select where you extracted folder
setwd("~/5hmC-Pathway-Analysis/")
#training
#counts_name <- "./Output/Randomization/1-2-3-4PMpos_8-9-11PMneg_DESeq2__final/1-2-3-4PMpos_8-9-11PMneg_training_rawcounts.csv"
#meta_name <- "./Output/Randomization/1-2-3-4PMpos_8-9-11PMneg_DESeq2__final/1-2-3-4PMpos_8-9-11PMneg_training_conditions.csv"
#validation
#counts_name <- "./Output/Randomization/1o3_pdONLYpos_8o11_metNEGtumNEG_DESeq2_v3/1o3_pdONLYpos_8o11_metNEGtumNEG_training_rawcounts.csv"
#meta_name <- "./Output/Randomization/1o3_pdONLYpos_8o11_metNEGtumNEG_DESeq2_v3/1o3_pdONLYpos_8o11_metNEGtumNEG_training_conditions.csv"
#whole dataset
counts_name <- "./Output/Raw Data Processing/1-2-3-4PMpos_8-9-11PMneg_vFINAL/1-2-3-4PMpos_8-9-11PMneg_DESeq2_rawcounts.csv"
meta_name <- "./Output/Raw Data Processing/1-2-3-4PMpos_8-9-11PMneg_vFINAL/1-2-3-4PMpos_8-9-11PMneg_DESeq2_conditions.csv"
#read in data, define what counts & conditions files
counts_data <- read.csv(counts_name,row.names = 1)
meta <-  read.csv(meta_name,row.names = 1)

#include to exclude specific sample
#counts_data <- subset(counts_data,select=-c(KT126))
#meta <- meta[!(row.names(meta) %in% c("KT126")),]

#define padj cutoff, you may need to run with several padj values until you have an appropriate number of significant results.
#used to select significant genes for results tables, PCA plots, heatmaps and UMAP plots.
cutoff_type = 1 # 0=padj cutoff, default; 1=lfc & pvalue cutoff
padj.cutoff = 0.05 # 0.1 default
pvalue.cutoff = 0.001
lfc.cutoff = 0.137504 # 0.137504 ~ 10%, 0.263034 ~ 20% change, 0.584963 ~ 50% change, 
min_samples = 0.8
min_count = 5

#Select version for all output files (e.g. 1, 2, 3, ...)

ver <- "pvalue0p001_lfc0p213_FINAL_wF"
gene_number <- nrow(counts_data)

#Set desired outputs:
output_results_tables = 1
output_volcano = 1
output_heatmap = 1
output_PCA = 1
output_UMAP = 0
output_config = 1
output_xtabs = 1

# define contrast groups
groups <- unique(meta[c("condition")])
contrast_groups <- c("condition",groups[1,1], groups[2,1])

###VALIDATION
#check columns are equal
all(colnames(counts_data) %in% rownames(meta))
all(colnames(counts_data) == rownames(meta))

###CREATE DESeq2 OBJECT
#load data into DESeq2 object dds
design_formula <- ~ condition #+ sex + age + race + batch
#design_formula <- ~ condition + batch
dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = meta, design = design_formula)

#load up size factors into dds object in order to normalize using median of ratios method
dds <- DESeq(dds)

keep <- rowSums( counts(dds) >= 1 ) >= 2
#dds <- dds[keep,]

#use counts 'normalized=true' function to pull out normalized counts
normalized_counts <- counts(dds,normalized=TRUE)

#drop genes with insufficient counts
normalized_counts <- normalized_counts[rowSums(normalized_counts >= min_count) >= min_samples, ]

#Save normalized counts (uncommment if desired, will save to results folder.)
write.csv(normalized_counts, file=paste("./Output/DESeq2/Results/normalized_counts_",contrast_groups[2],contrast_groups[3],"_",ver,".csv", sep = ""))


# extract results table, and shrink the log2 fold changes
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

if (file.exists("./Output/DESeq2/Results/")) {
  cat("The folder already exists")
} else {
  dir.create("./Output/DESeq2/Results/")
}

if (file.exists("./Output/DESeq2/Volcano/")) {
  cat("The folder already exists")
} else {
  dir.create("./Output/DESeq2/Volcano/")
}

if (file.exists("./Output/DESeq2/Heatmaps/")) {
  cat("The folder already exists")
} else {
  dir.create("./Output/DESeq2/Heatmaps/")
}

if (file.exists("./Output/DESeq2/PCA/")) {
  cat("The folder already exists")
} else {
  dir.create("./Output/DESeq2/PCA/")
}

if (file.exists("./Output/DESeq2/Xtabs/")) {
  cat("The folder already exists")
} else {
  dir.create("./Output/DESeq2/Xtabs/")
}

if (file.exists("./Output/DESeq2/Config/")) {
  cat("The folder already exists")
} else {
  dir.create("./Output/DESeq2/Config/")
}

###SAVE RESULTS TABLES TO TEXT FILES
if(output_results_tables == 1){
  write.table(res_table, file=paste("./Output/DESeq2/Results/all_results_",contrast_groups[2],contrast_groups[3],"_",ver,".txt", sep = ""), sep="\t", quote=F, col.names=NA)
  write.table(sig, file=paste("./Output/DESeq2/Results/significant_results_",contrast_groups[2],contrast_groups[3],"_",ver,".txt", sep = ""), sep="\t", quote=F, col.names=NA)
}

###GENERATE VOLCANO PLOTS

if(output_volcano == 1){
  #Obtain logical vector where TRUE values denote padj values < cutoff and fold change > 1.5 in either direction
  #working with log2 fold changes so this translates to an actual fold change of 2^lfc.cutoff. used in volcano plots
  if (cutoff_type == 0){
  res_table_tb_volcano <- res_table_tb %>% 
    mutate(threshold_OE = padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff)
  
  res_table_tb_volcano <- res_table_tb_volcano %>% arrange(padj) %>% mutate(genelabels = "")
  
  res_table_tb_volcano$genelabels[1:nsig] <- res_table_tb_volcano$gene[1:nsig]
  
  png(paste("./Output/DESeq2/Volcano/volcano_",contrast_groups[2],contrast_groups[3],"_",ver,".png", sep = ""), width = 900, height = 900)
    print({
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
      })
  dev.off()
  }
  
  if (cutoff_type == 1){
    res_table_tb_volcano <- res_table_tb %>% 
      mutate(threshold_OE = pvalue < pvalue.cutoff & abs(log2FoldChange) >= lfc.cutoff)
    
    res_table_tb_volcano <- res_table_tb_volcano %>% arrange(desc(threshold_OE),pvalue) %>% mutate(genelabels = "")
    
    res_table_tb_volcano$genelabels[1:nsig] <- res_table_tb_volcano$gene[1:nsig]
    
    lfc_max = max(res_table_tb_volcano$log2FoldChange)
    lfc_min = -min(res_table_tb_volcano$log2FoldChange)
    lfc_max = max(lfc_max,lfc_min)
    
    png(paste("./Output/DESeq2/Volcano/volcano_",contrast_groups[2],contrast_groups[3],"_",ver,".png", sep = ""), width = 900, height = 900)
    print({
      ggplot(res_table_tb_volcano, aes(x = log2FoldChange, y = -log10(pvalue))) +
      geom_point(aes(colour = threshold_OE)) +
      ggtitle(paste(contrast_groups[2],"/",contrast_groups[3],"Enrichment")) +
      geom_text_repel(aes(label = genelabels,size=20)) +
      xlab("log2 fold change") + 
      ylab("-log10 p-value") +
      theme(legend.position = "none",
            plot.title = element_text(size = rel(1.5), hjust = 0.5),
            axis.title = element_text(size = rel(1.25))) +
      xlim(-lfc_max, lfc_min)
      })
    dev.off()
  }
}

###GENERATE HEATMAP
if(output_heatmap == 1){
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
  
  log_counts_data_mat = data.matrix(norm_sig, rownames.force = NA)
  log_counts_data_mat = t(scale(t(log_counts_data_mat)))
  ha = HeatmapAnnotation(condition=annotation$condition,col = list(condition = c("8-9-11PMneg" = "#ff9289", "1-2-3-4PMpos" = "#00dae0")))
  
  svg(paste("./Output/DESeq2/Heatmaps/sig_heatmap_",contrast_groups[2],contrast_groups[3],"_",ver,".svg", sep = ""), width = 1500, height = 1200)
  Heatmap(log_counts_data_mat, 
          top_annotation = ha,
          #col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
          col = colorRamp2(c(-4, -2,-1,0,1, 2,4), c("#4676b5", "#82b1d3","#dbeff6","#fefebd","#fee395", "#fc9961","#d73027")),
          clustering_distance_rows = "pearson",
          clustering_distance_columns = robust_dist,
#  pheatmap(norm_sig, 
#           main = heatmap_title,
# #         color = diverging_hcl(15,"Blue-Red2"), 
#           cluster_rows = T,
#           clustering_distance_cols = "euclidean",
##           clustering_distance_rows = "manhattan",
#           cluster_cols = T,
#           show_rownames = T,
#           show_colnames = F,
#           annotation = annotation, 
#           border_color = NA, 
#           fontsize = 10, 
#           scale = "row", 
#           fontsize_row = 10, 
#           height = 20,
  )
  dev.off()
}

###GENERATE PCA PLOT
if(output_PCA == 1) {
  #See details for below operation in lesson 3 of DGE workshop
  rld <- vst(dds, blind=TRUE)
  res_genes <- row.names(res_table)
  sig_genes <-  sig$gene
  
  for(each in c(cols,"condition")){
    #In the below replace sig_genes with res_genes if you want to perform PCA analysis on all genes rather than just on significant genes.
    plotPCA_labels <- plotPCA(rld[sig_genes,], intgroup = c(each)) #rownames(c(row.names(meta)))
    #reverse order of labels to match heatmap labelling
    #plotPCA_labels$data$group<-factor(plotPCA_labels$data$group,levels = rev(levels(plotPCA_labels$data$group)))
    #plotPCA_labels$data$condition<-factor(plotPCA_labels$data$condition,levels = rev(levels(plotPCA_labels$data$condition)))
    
    
    svg(paste("./Output/DESeq2/PCA/sig_PCA_",contrast_groups[2],contrast_groups[3],"_",ver,each,".svg", sep = ""))
      #plotPCA_labels + geom_text(aes(label = name),position=position_nudge(y = 0.07),) + ggtitle(heatmap_title) 
    print({plotPCA_labels + ggtitle(heatmap_title) + stat_ellipse(level=0.9)}, cex.main=3)
    dev.off()
  }
}

###GENERATE UMAP PLOT (Curent issues with M3C library)
#if (file.exists("./Output/DESeq2/UMAP/")) {
#  cat("The folder already exists")
#} else {
#  dir.create("./Output/DESeq2/UMAP/")
#}

#save UMAP plot to png
#png(paste("./Output/DESeq2/UMAP/sig_UMAP_",contrast_groups[2],contrast_groups[3],"_",ver,".png", sep = ""), width = 900, height = 1200)
#umap(norm_sig, labels=as.factor(meta$condition),printres = FALSE, seed = FALSE,
#     axistextsize = 18, legendtextsize = 18, dotsize = 5,
#     textlabelsize = 4, legendtitle = "Group", controlscale = FALSE,
#     scale = 1,     printwidth = 22, text = FALSE)
#dev.off()

### SAVE CROSS TABULATION CHARTS
if(output_xtabs == 1){
  # Convert meta columns to factors (except for condition and age)
  cols <- c(colnames(meta[,!(names(meta) %in% c("condition","age"))]))
  meta_factor <- meta
  meta_factor[cols] <- lapply(meta_factor[cols],factor)
  summary(meta_factor)
  
  #Save PlotXTabs2 https://www.rdocumentation.org/packages/CGPfunctions/versions/0.6.3/topics/PlotXTabs2 to a grid
  
  jpeg(paste("./Output/DESeq2/Xtabs/",contrast_groups[2],contrast_groups[3],"_",ver,".png", sep = ""),
       width=2600, height=1000, res = 144)
    plot_sex <- PlotXTabs2(meta_factor,condition,sex,title="Sex")
    plot_race <- PlotXTabs2(meta_factor,condition,race, title="Race")
    plot_batch <- PlotXTabs2(meta_factor,condition,batch, title = "Batch")
    #plot_primary <- PlotXTabs2(meta_factor,condition,primary_present, title = "Primary Present")
    plot_site <- PlotXTabs2(meta_factor,condition,primary_site, title = "Primary Site")
    plot_chemo <- PlotXTabs2(meta_factor,condition,chemo_6weeks, title = "Chemo <6weeks")
    plot_age <- ggboxplot(meta_factor, x = "condition", y = "age", 
              color = "condition", palette = c("#00AFBB", "#E7B800"),
              ylab = "Age", xlab = "Condition", title = "Age") + stat_compare_means(method = "t.test")
    #  include other mets plot if applicable, and add "plot_other" to ggarrange.
    #  plot_other <- PlotXTabs2(meta_factor,condition,other_mets)
  
    ggarrange(plot_sex, plot_race, plot_batch, plot_chemo, plot_site, plot_age, ncol = 3, nrow = 2)
  #            labels = c("Sex", "Race", "Batch", "Primary Present","Primary Site", "Chemo 6 Weeks", "Age", "Other Mets (n/a"),
  #            vjust	= 0.5,
          
          
    dev.off()
}

###SAVE CONFIG TABLES TO TEXT FILES
if(output_config == 1){
  # Create folders
  
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
  
  write.table(config, file=paste("./Output/DESeq2/Config/config_",contrast_groups[2],contrast_groups[3],"_",ver,".txt", sep = ""), sep="\t", quote=F, col.names=NA)
}

