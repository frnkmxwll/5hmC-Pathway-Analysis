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
#library(BiocManager)
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
#library(CGPfunctions) #needed for plotxtabs2 
library(ggpubr) # needed for whisker plot charting
#library(aod) #for logistic regression
library(broom)
library(leaps) #for logistic regression optimizaton
library(devtools) #used for complex heatmaps
library(ComplexHeatmap)#used for complex heatmaps
library(circlize)
library(BiocParallel)
# Source the file containing annotate_features2 function
source("~/5hmC-Pathway-Analysis/R Code/repeated_glmnet_funcs_v1.R")

###CONFIGURATION
#set working directory, select where you extracted folder
setwd("~/5hmC-Pathway-Analysis/")

comp_set <- "PM" # AM or PM
if(comp_set == "AM"){
  # # # Peak analysis ----
  # subset_conditions <- list(list("condition", list("AMpos","AMneg")))
  # counts_name <- "./Output/Raw Data Processing/AMneg_AMpos_90per15_10292024_hmr_combat_ampos_amneg_noX/AMneg_AMpos_DESeq2_rawcounts.csv"
  # meta_name <- "./Output/Raw Data Processing/AMneg_AMpos_90per15_10292024_hmr_combat_ampos_amneg_noX/AMneg_AMpos_DESeq2_conditions.csv"
  # ver.text <- "15co_11_09_2024_peaks_combat_AM_noPM"
  # design_formula <- ~ condition + age_norm + race + sex + ovr_histopath + peritoneal_mets # primary_site # primary_present # peritoneal_mets
  # region_type= 0 # 0=macs_peak_locations 1=gene symbols
  # cutoff_type = 1 # 0=padj cutoff, default; 1=lfc & pvalue cutoff; 2 = top stat genes
  # padj.cutoff = 0.001 # 0.1 default
  # pvalue.cutoff = 0.001 #0.001 for heatmap; 0.01 for ORA; lfc 0.26 for heatmap; 0.2 for ORA
  # num_top_genes = 2000
  # lfc.cutoff = 0.2 # 0.0704 ~ 5%, 0.137504 ~ 10%, 0.2016 ~ 15% change 0.263034 ~ 20% change, 0.32 ~ 25%, 0.415 ~ 33%, 0.58 ~ 50% change,
  # genebodies_type = FALSE
  # randomize_conds = FALSE
  # if(randomize_conds == TRUE){
  #   ver.text <- paste0(ver.text,"_RAND")
  # }
  
  #Genebod analysis ----
  subset_conditions <- list(list("condition", list("AMpos","AMneg")))
  counts_name <- "./Output/Raw Data Processing/AMneg_AMpos_0_gene_bodies_10302024_hmr_combat_ampos_amneg_noX/AMneg_AMpos_DESeq2_rawcounts.csv"
  meta_name <- "./Output/Raw Data Processing/AMneg_AMpos_0_gene_bodies_10302024_hmr_combat_ampos_amneg_noX/AMneg_AMpos_DESeq2_conditions.csv"
  ver.text <- "0_genebodies_11_9_2024_combat_AM"
  design_formula <- ~ condition + age_norm + race + sex + ovr_histopath + peritoneal_mets  # primary_site # primary_present # peritoneal_mets
  region_type= 0 # 0=macs_peak_locations 1=gene symbols
  cutoff_type = 1 # 0=padj cutoff, default; 1=lfc & pvalue cutoff; 2 = top stat genes
  padj.cutoff = 0.001 # 0.1 default
  pvalue.cutoff = 0.01 #0.001 for heatmap; 0.01 for ORA; lfc 0.26 for heatmap; 0.2 for ORA
  num_top_genes = 2000
  lfc.cutoff = 0.263 # 0.0704 ~ 5%, 0.137504 ~ 10%, 0.2016 ~ 15% change 0.263034 ~ 20% change, 0.32 ~ 25%, 0.415 ~ 33%, 0.58 ~ 50% change,
  genebodies_type = TRUE
  randomize_conds = FALSE
  # if randomize_conds == TRUE  add "_RAND" to ver.text:
  if(randomize_conds == TRUE){
    ver.text <- paste0(ver.text,"_RAND")
  }
  
  ## PM Settings ----
} else if(comp_set == "PM"){
  subset_conditions <- list(list("condition", list("PM_negative","PM_positive")))
  counts_name <- "./Output/Raw Data Processing/HEALTHY_PM_positive_10per_11232024_hmr_combat_healthy_noX/HEALTHY_PM_positive_DESeq2_rawcounts.csv"
  meta_name <- "./Output/Raw Data Processing/HEALTHY_PM_positive_10per_11232024_hmr_combat_healthy_noX/HEALTHY_PM_positive_DESeq2_conditions.csv"
  ver.text <- "12_20_2024_peaks_combat_PM"
  design_formula <- ~ condition + batch + age_norm + race + sex + primary_site #primary_present
  #+ primary_site #primary_site #+ peritoneal_mets+ ovr_histopath  + primary_site + peritoneal_mets
  region_type= 0 # 0=macs_peak_locations 1=gene symbols
  cutoff_type = 2 # 0=padj cutoff, default; 1=lfc & pvalue cutoff; 2 = top stat genes
  padj.cutoff = 0.05 # 0.1 default
  pvalue.cutoff = 0.001
  num_top_genes = 2000
  lfc.cutoff = 0.2 # 0.0704 ~ 5%, 0.137504 ~ 10%, 0.2016 ~ 15% change 0.263034 ~ 20% change, 0.32 ~ 25%, 0.415 ~ 33%, 0.58 ~ 50% change, 
  genebodies_type = FALSE
  randomize_conds = FALSE
  # if randomize_conds == TRUE  add "_RAND" to ver.text:
  if(randomize_conds == TRUE){
    ver.text <- paste0(ver.text,"_RAND")
  }
}

#include to exclude specific sample
#counts_data <- subset(counts_data,select=-c(KT126))
#meta <- meta[!(row.names(meta) %in% c("KT126")),]

if(cutoff_type == 0){
  ver <- paste0("padj",padj.cutoff,"_lfc",lfc.cutoff,"_",ver.text)
}

if(cutoff_type == 1){
  ver <- paste0("pval",pvalue.cutoff,"_lfc",lfc.cutoff,"_",ver.text)
}

if(cutoff_type == 2){
  ver <- paste0("ntop",num_top_genes,"_",ver.text)
}

#Set desired outputs:
output_results_tables = 1
output_volcano = 0
output_heatmap = 1
output_PCA = 0
output_UMAP = 0
output_config = 1
output_xtabs = 0
output_norm_counts = 0 #(normcounts, vsg and rlog normalized files)


### LOAD UP DATA ----
#read in data, define what counts & conditions files
counts_data <- read.csv(counts_name,row.names = 1)
meta_complete <-  read.csv(meta_name,row.names = 1)
gene_number <- nrow(counts_data)

# randomize condition labels
set.seed(123)
if(randomize_conds == TRUE){
  meta_complete$condition <- sample(meta_complete$condition)
}

if(genebodies_type == TRUE){
  # # drop rows from counts_data that don't have at least 10% of columns with 10 reads or more 
  keep <- rowSums(counts_data >= 10) >= 0.1 * ncol(counts_data) # for genes

  # # print rows kept
  print(paste("Rows kept: ", sum(keep), " out of ", nrow(counts_data), " (", sum(keep) / nrow(counts_data) * 100, "%)"))
  counts_data <- counts_data[keep,]
}

#subset_conditions <- NULL


source("~/5hmC-Pathway-Analysis/R Code/repeated_glmnet_funcs_v1.R")
color_tables <- create_color_tables()

meta_complete <- meta_complete %>%
  mutate(
    peritoneal_mets = ifelse(condition == "HEALTHY", "healthy", peritoneal_mets),
    other_mets = ifelse(condition == "HEALTHY", "healthy", other_mets),
    primary_present = ifelse(condition == "HEALTHY", "healthy", primary_present),
    chemo_6weeks = ifelse(condition == "HEALTHY", "healthy", chemo_6weeks),
    primary_site = ifelse(condition == "HEALTHY", "healthy", primary_site)
  )

# replace forbidden characters "/" or "-" with "_" in values
meta_complete$race <- gsub("/", "_", meta_complete$race)
meta_complete$race <- gsub("-", "_", meta_complete$race)

# Add a column called age_norm that corresponds to centered and scaled age column
meta_complete$age_norm <- as.numeric(scale(meta_complete$age, center = TRUE, scale = TRUE))

# Add a column called "age_cat" that corresponds to age categories low or high (below or above median)
meta_complete$age_cat <- ifelse(meta_complete$age < median(meta_complete$age), "low", "high")

# if AM, then load file from ~/5hmC-Pathway-Analysis/Output/Raw Data Processing/AMneg_AMpos_10_90per_10282024_hmr_combat_ampos_amneg_noX/
if(comp_set == "AM"){
  # read in csv
  updated_race_df = read.csv("~/5hmC-Pathway-Analysis/Output/Raw Data Processing/AMneg_AMpos_10_90per_10282024_hmr_combat_ampos_amneg_noX/updated_race.csv")
  # replace "0" with "zero", "1" with "one", etc... up to 4 in race column
  # Replace values in the 'race' column with corresponding strings
  updated_race_df$race <- as.character(updated_race_df$race)
  updated_race_df$race[updated_race_df$race == "0"] <- "zero"
  updated_race_df$race[updated_race_df$race == "1"] <- "one"
  updated_race_df$race[updated_race_df$race == "2"] <- "two"
  updated_race_df$race[updated_race_df$race == "3"] <- "three"
  updated_race_df$race[updated_race_df$race == "4"] <- "four"
  
  # merge updated race column from updated_race_df with meta_complete, where ID are row names. df looks like:
  # ID,race
  # KT001,0
  # Ensure 'ID' column is set as row names in updated_race_df for compatibility
  rownames(updated_race_df) <- updated_race_df$ID
  
  # Assuming `meta_complete` is already defined and has matching IDs in its row names
  # Merge updated race column from updated_race_df into meta_complete by matching row names
  meta_complete$race <- updated_race_df[rownames(meta_complete), "race"]
  
  # Verify the merge result
  print(head(meta_complete))
}

# set race_cat to "white" or "other"
meta_complete$race_cat <- ifelse(meta_complete$race == "White", "White", "Other")


### Correct metadata fields ----
# Subset meta based on the conditions provided
if (!is.null(subset_conditions)) {
  for (condition in subset_conditions) {
    column_name <- condition[[1]]
    value <- unlist(condition[[2]])
    if (column_name %in% colnames(meta_complete)) {
      # keep only rows where column value in column_name are in value
      meta <- meta_complete[meta_complete[, column_name] %in% value, , drop = FALSE]
    } else {
      warning(paste("Column", column_name, "not found in meta data. Skipping this condition."))
    }
  }
} else{
  meta <- meta_complete
}


# Convert categorical variables to factors
# meta$condition <- factor(meta$condition)
# meta$race <- factor(meta$race)
# meta$sex <- factor(meta$sex)
# meta$ovr_histopath <- factor(meta$ovr_histopath)
# meta$peritoneal_mets <- factor(meta$peritoneal_mets)
# meta$primary_site <- factor(meta$primary_site)
# meta$primary_present <- factor(meta$primary_present)
# meta$chemo_6weeks <- factor(meta$chemo_6weeks)
# meta$other_mets <- factor(meta$other_mets)
# meta$primary_site <- factor(meta$primary_site)
# meta$age_cat <- factor(meta$age_cat)
# meta$race_cat <- factor(meta$race_cat)

counts_data <- counts_data[, colnames(counts_data) %in% rownames(meta)]


# In condition column, replace any value not equal to HEALTHY with "DISEASED"
#meta$condition <- ifelse(meta$condition != "HEALTHY", "DISEASED", "HEALTHY")

# define contrast groups
groups <- unique(meta[c("condition")])
contrast_groups <- c("condition",groups[2,1],groups[1,1])

#check columns are equal
all(colnames(counts_data) %in% rownames(meta))
all(colnames(counts_data) == rownames(meta))

###CREATE DESeq2 OBJECT ----
#load data into DESeq2 object dds
BiocParallel::register(MulticoreParam(workers=20), default = TRUE)

dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = meta, design = design_formula)

###

#load up size factors into dds object in order to normalize using median of ratios method
dds <- DESeq(dds,parallel=TRUE)


#use counts 'normalized=true' function to pull out normalized counts
normalized_counts <- counts(dds,normalized=TRUE)

# extract results table, and shrink the log2 fold changes
if(cutoff_type == 0){
  res_table <- results(dds, contrast=contrast_groups, alpha = padj.cutoff)
}
if(cutoff_type == 1){
  res_table <- results(dds, contrast=contrast_groups)
}
if(cutoff_type == 2){
  res_table <- results(dds, contrast=contrast_groups)
}

#Consider replacing above line with shrunken values for fold change below:
#res_table_unshrunken <- results(dds, contrast=contrast_groups, alpha = padj.cutoff)
#res_table <- lfcShrink(dds, coef = 2, res=res_table_unshrunken)

metadata(res_table)$filterThreshold

# plot(metadata(res_table)$filterNumRej, 
#      type="b", ylab="number of rejections",
#      xlab="quantiles of filter")
# lines(metadata(res_table)$lo.fit, col="red")
# abline(v=metadata(res_table)$filterTheta)
# dev.off()

#convert the results table into a tibble:
res_table_tb <- res_table %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

#subset table to only keep significant genes using cutoff
if(cutoff_type == 0){
  sig <- res_table_tb %>% 
    filter(padj <= padj.cutoff, abs(log2FoldChange) >= lfc.cutoff)
}

if(cutoff_type == 1){
  sig <- res_table_tb %>% 
    filter(abs(log2FoldChange) >= lfc.cutoff, pvalue <= pvalue.cutoff)
}

if(cutoff_type == 2){
  # Sort by stat in descending order
  res_table_sorted <- res_table_tb %>% 
    arrange(desc(abs(stat)))
  
  # Filter res_table_tb down to num_top_genes first rows only
  sig <- res_table_sorted %>% 
    slice_head(n = num_top_genes)
}

normalized_counts_tb <- data.frame(normalized_counts) %>% 
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

### OUPUT NORM  COUNTS ----
if (output_norm_counts==1){
  vst_data <- vst(dds, blind = TRUE)  # perform variance stabilizing transformation
  vst_data <- assay(vst_data)  # extract transformed data for glmnet
  
  rlog_data <- rlog(dds, blind = TRUE)  # perform variance stabilizing transformation
  rlog_data <- assay(rlog_data)  # extract transformed data for glmnet
  
  #Save normalized counts (uncommment if desired, will save to results folder.)
  write.csv(normalized_counts, file=paste("./Output/DESeq2/Results/normalized_counts_",contrast_groups[2],contrast_groups[3],"_",ver,".csv", sep = ""))
  write.csv(vst_data, file=paste("./Output/DESeq2/Results/vst_counts_",contrast_groups[2],contrast_groups[3],"_",ver,".csv", sep = ""))
  write.csv(rlog_data, file=paste("./Output/DESeq2/Results/rlog_counts_",contrast_groups[2],contrast_groups[3],"_",ver,".csv", sep = ""))
  
  vst_counts_tb <- vst_data %>% 
    data.frame() %>%
    rownames_to_column(var="gene") %>% 
    as_tibble()
  
  vst_sig <- vst_counts_tb[,c(1,2:nsamples)] %>% 
    filter(gene %in% sig$gene) %>% 
    data.frame() %>%
    column_to_rownames(var = "gene")
}

### CREATE ROOT OUTPUT FOLDERS IF THEY DON'T EXIST ----
if(TRUE){
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
}

###GENERATE VOLCANO PLOTS----
if(output_volcano == 1){
  if (cutoff_type == 0){
    heatmap_title <- paste(contrast_groups[2], "/", contrast_groups[3], "padj <", padj.cutoff)
  }
  
  if (cutoff_type == 1){
    heatmap_title <- paste(contrast_groups[2], "/", contrast_groups[3], "pvalue <", pvalue.cutoff, "|lfc|>", lfc.cutoff)
  }
  
  if (cutoff_type == 2){
    heatmap_title <- paste(contrast_groups[2], "/", contrast_groups[3], "top", num_top_genes, "genes")
  }
  
  res_table_tb_volcano <- res_table_tb %>% 
    arrange(desc(abs(stat))) %>%
    mutate(rank = row_number())
  
  if (cutoff_type == 2) {
    sig_genes <- res_table_tb_volcano %>% 
      filter(rank <= num_top_genes)
  } else {
    sig_genes <- res_table_tb_volcano %>% 
      filter((cutoff_type == 0 & padj < padj.cutoff | cutoff_type == 1 & pvalue < pvalue.cutoff) & 
               abs(log2FoldChange) > lfc.cutoff)
  }
  
  res_table_tb_volcano <- res_table_tb_volcano %>%
    mutate(threshold_OE = case_when(
      gene %in% sig_genes$gene & log2FoldChange >= 0 ~ "Upregulated",
      gene %in% sig_genes$gene & log2FoldChange < 0 ~ "Downregulated",
      TRUE ~ "Not significant"
    ))
  
  lfc_max <- max(abs(sig_genes$log2FoldChange)) * 1.1
  y_max <- if(cutoff_type == 0) max(-log10(sig_genes$padj)) * 1.1 else max(-log10(sig_genes$pvalue)) * 1.1
  
  # Round up y_max to the next whole number
  y_max <- ceiling(y_max)
  
  # round lfc_max to nearest 0.1
  lfc_max <- round(lfc_max, 1)
  # Define the increment
  increment <- 0.2
  
  # Determine the smallest and largest multiples of increment within the range
  min_val <- -lfc_max
  max_val <- lfc_max
  
  # Calculate the sequence from the smallest multiple to the largest multiple
  sequence <- seq(from = floor(min_val / increment) * increment,
                  to = ceiling(max_val / increment) * increment,
                  by = increment)
  
  # Ensure 0 is included if it falls within the range
  if (0 < min(sequence) || 0 > max(sequence)) {
    sequence <- c(sequence, 0)
  }
  
  # Remove duplicates and sort the sequence
  sequence <- sort(unique(round(sequence, 1)))
  
  plot_volcano <- function(data, y_col, cutoff) {
    p <- ggplot(data, aes(x = log2FoldChange, y = -log10(!!sym(y_col)))) +
      theme_minimal() +
      geom_point(aes(colour = threshold_OE)) +
      scale_color_manual(values = c("Upregulated" = "#E31A1CB6", "Downregulated" = "#1F78B4B6", "Not significant" = "#CCCCCC"),
                         labels = c("Downregulated", "Not significant", "Upregulated"),
                         name = "Gene regulation") +
      ggtitle(heatmap_title) +
      xlab("log2 fold change") + 
      ylab(paste("-log10", if(y_col == "padj") "adjusted" else "", "p-value")) +
      theme(legend.position = "right",
            legend.title = element_text(size = rel(1)),
            legend.text = element_text(size = rel(0.8)),
            plot.title = element_text(size = rel(1.5), hjust = 0.5),
            axis.title = element_text(size = rel(1.25)),
            axis.text.x = element_text(size = rel(1.2)),
            axis.text.y = element_text(size = rel(1.2)),
            panel.grid.major = element_line(color = "grey90", size = 0.2),
            panel.grid.minor = element_line(color = "grey95", size = 0.1)) +
      coord_cartesian(xlim = c(-lfc_max, lfc_max), ylim = c(0, y_max)) +
      scale_x_continuous(breaks = sequence) +
      scale_y_continuous(breaks = seq(0, y_max, 2))
    
    if (cutoff_type != 2) {
      p <- p +
        geom_vline(xintercept = c(-lfc.cutoff, lfc.cutoff), linetype = "dashed", color = "darkgrey") +
        geom_hline(yintercept = -log10(cutoff), linetype = "dashed", color = "darkgrey") +
        annotate("text", x = lfc.cutoff, y = y_max, label = paste0("+", lfc.cutoff), hjust = -0.2, color = "darkgrey") +
        annotate("text", x = -lfc.cutoff, y = y_max, label = paste0("-", lfc.cutoff), hjust = 1.2, color = "darkgrey") +
        annotate("text", x = 0, y = -log10(cutoff), label = paste0("p-value cutoff: ", cutoff), vjust = -1, hjust = 1, color = "grey")
    }
    
    return(p)
  }
  
  if (cutoff_type == 0) {
    p <- plot_volcano(res_table_tb_volcano, "padj", padj.cutoff)
  } else if (cutoff_type == 1) {
    p <- plot_volcano(res_table_tb_volcano, "pvalue", pvalue.cutoff)
  } else {
    p <- plot_volcano(res_table_tb_volcano, "pvalue", NA)
  }
  
  ggsave(paste0("./Output/DESeq2/Volcano/volcano_", contrast_groups[2], contrast_groups[3], "_", ver, ".png"), 
         plot = p, width = 6, height = 4, dpi = 300)
  
  ggsave(paste0("./Output/DESeq2/Volcano/volcano_", contrast_groups[2], contrast_groups[3], "_", ver, ".svg"), 
         plot = p, width = 6, height = 4)
}

###GENERATE HEATMAP ----
# set to true if you want to include healthy in heatmap
if(FALSE){
  ### For plotting healthy with PMpos and PM neg in heatmap
  counts_name_all <- counts_name
  
  ### Gene Bodies
  #meta_name_all <- meta_name #"~/5hmC-Pathway-Analysis/Output/Raw Data Processing/HEALTHY_PM_positive_genebodies_082624_nocombat_pm_healthy_X/HEALTHY_PM_positive_DESeq2_conditions.csv"
  #counts_name_all <- counts_name # "~/5hmC-Pathway-Analysis/Output/Raw Data Processing/HEALTHY_PM_positive_genebodies_082624_nocombat_pm_healthy_X/HEALTHY_PM_positive_DESeq2_rawcounts.csv"
  
  counts_data_all <- read.csv(counts_name_all,row.names = 1)
  meta_all <- meta_complete
  
  # keep only rows where batch = batch_2
  meta_all <- meta_all[meta_all$batch == "batch_2",]
  
  # drop rows where condition == PM_negative
  #meta_all <- meta_all[meta_all$condition != "HEALTHY",]
  
  # drop columns from counts_data_all that are not in meta_all
  counts_data_all <- counts_data_all[, colnames(counts_data_all) %in% rownames(meta_all)]
  
  groups_all <- unique(meta_all[c("condition")])
  contrast_groups <- NULL
  # for each element in groups_all add it to contrast_groups
  for (i in 1:nrow(groups_all)) {
    contrast_groups <- c(contrast_groups, groups_all[i,1])
  }
  
  # drop NA elements
  contrast_groups <- contrast_groups[!is.na(contrast_groups)]
  
  # Create DESeq2 object without specifying a design
  dds_all <- DESeqDataSetFromMatrix(countData = counts_data_all, 
                                    colData = meta_all, 
                                    design = ~ 1)  # Using '~ 1' as a placeholder design
  #dds_all <- dds_all[keep,]
  
  # Estimate size factors and normalize
  dds_all <- estimateSizeFactors(dds_all)
  normalized_counts_all <- counts(dds_all, normalized=TRUE)
  normalized_counts_all_tb <- data.frame(normalized_counts_all) %>% 
    rownames_to_column(var="gene") %>% 
    as_tibble()
  
  ### --- If plotting boxplots for enrichment in different types ---
  # read in all_results_annotated_enhancer file
  #all_results_annotated_enhancer <- read.csv("~/5hmC-Pathway-Analysis/Output/DESeq2/Results/all_results_annotated_enhancer_PM_positivePM_negative_pval0.001_lfc0.2_10_13_2024_hmr_combat.txt",
  #                                           header = TRUE, sep = "\t")
  
  # keep only columns called "feature" and "locationType"
  #all_results_annotated_enhancer <- all_results_annotated_enhancer[,c("feature","locationType")]
  
  # # keep only rows where locationType is "1-5kb Upstream"
  # all_results_annotated_enhancer <- all_results_annotated_enhancer[
  #   all_results_annotated_enhancer$locationType == "1-5kb Upstream",]
  # norm_sig_all <- normalized_counts_all_tb %>%  #[,c(1,2:nsamples_all)]
  #   filter(gene %in% all_results_annotated_enhancer$feature) %>%
  #   data.frame() %>%
  #   column_to_rownames(var = "gene")
  
  log_counts_dataframe <- normalized_counts_all_tb
  
  nsamples_all <- ncol(normalized_counts_tb)
  ## Plot only sig genes
  #Extract normalized expression for significant genes from the two comparison groups and set the gene column (1) to row names
  
  norm_sig_all <- normalized_counts_all_tb %>%  #[,c(1,2:nsamples_all)]
    filter(gene %in% sig$gene) %>%
    data.frame() %>%
    column_to_rownames(var = "gene")
  
  # add "_healthy" to ver.text
  ver <- paste0(ver,"_healthy")
  
} else {
  meta_all <- meta
  norm_sig_all <- norm_sig
  log_counts_dataframe <- normalized_counts_tb
  dds_all <- dds
}

if(output_heatmap == 1){
  source("~/5hmC-Pathway-Analysis/R Code/repeated_glmnet_funcs_v1.R")
  color_tables <- create_color_tables()
  
  # annotation <- meta_all %>%
  #   dplyr::select(condition, primary_present,  ovr_histopath, chemo_6weeks, primary_site,peritoneal_mets)
  annotation <- meta_all 
  
  if(cutoff_type == 0){
    heatmap_title <- paste(contrast_groups[2],"/",contrast_groups[3],"padj <",padj.cutoff)
  }
  
  if(cutoff_type == 1){
    heatmap_title <- paste(contrast_groups[2],"/",contrast_groups[3],"pvalue <",pvalue.cutoff,"|lfc|>",lfc.cutoff)
  }
  
  if(cutoff_type == 2){
    heatmap_title <- paste(contrast_groups[2],"/",contrast_groups[3],"top ",num_top_genes," genes")
  }
  
  
  # convert to dataframe
  log_counts_dataframe <- as.data.frame(log_counts_dataframe)
  
  # convert "gene" column to row names
  rownames(log_counts_dataframe) <- log_counts_dataframe$gene
  # drop gene column
  log_counts_dataframe <- log_counts_dataframe[,-1]
  
  #log_counts_data_mat <- data.matrix(norm_sig_all, rownames.force = NA)
  log_counts_data_mat <- data.matrix(log_counts_dataframe, rownames.force = NA)
  
  # Calculate z-scores for all rows
  zscore_matrix <- t(scale(t(log_counts_data_mat), center = TRUE, scale = TRUE))
  
  # Ensure rownames are preserved
  rownames(zscore_matrix) <- rownames(log_counts_dataframe)
  
  # Filter the z-score matrix to keep only rows present in norm_sig_all
  zscore_matrix_filtered <- zscore_matrix[rownames(zscore_matrix) %in% rownames(norm_sig_all), ]
  
  zscore_matrix_filtered <- data.matrix(zscore_matrix_filtered, rownames.force = NA)
  
  
  # group1 <- paste(contrast_groups[[3]])
  # group2 <- paste(contrast_groups[[2]])
  # # if [[4]] is not null then group3 <- paste(contrast_groups[[4]])
  # if (length(contrast_groups) > 3) {
  #   group3 <- paste(contrast_groups[[4]])
  # }
  # 
  # # Keep the condition_table colors unchanged
  # if (length(contrast_groups) > 3) {
  #   condition_table <- c("#1F78B4A6", "#E31A1CA6", "#CCCCCC")
  #   names(condition_table) <- c(group1, group2, group3)
  # } else {
  #   condition_table <- c("#1F78B4A6", "#E31A1CA6")
  #   names(condition_table) <- c(group1, group2)
  # }
  # 
  # # Updated colors using colorblind-safe palette
  # # Primary Site Colors
  # primary_site_colors <- c("#009E73", "lightgreen", "#CC79A7", "#CCCCCC")
  # primary_site_table <- setNames(primary_site_colors, c("CRC", "ADE", "appendiceal_mucocele", "healthy"))
  # 
  # # Pathology Colors
  # path_colors <- c("#56B4E9", "#009E73", "#CCCCCC")
  # path_table <- setNames(path_colors, c("AMN", "Adenocarcinoma", "unknown_or_na"))
  # 
  # # Primary Present Colors
  # primary_present_colors <- c("#F0E442", "#E69F00", "#CCCCCC")
  # primary_present_table <- setNames(primary_present_colors, c("primary_absent", "primary_present", "healthy"))
  # 
  # # Sex Colors
  # sex_colors <- c("lightpink", "#CC79A7")
  # sex_table <- setNames(sex_colors, c("male", "female"))
  # 
  # # Other Mets Colors
  # other_mets_colors <- c("#E69F00", "#009E73", "#CCCCCC")
  # other_mets_table <- setNames(other_mets_colors, c("other_mets_present", "other_mets_absent", "healthy"))
  # 
  # # Chemo Colors
  # chemo_colors <- c("#999999", "#56B4E9")
  # chemo_table <- setNames(chemo_colors, c("No", "Yes"))
  # 
  # # Batch Colors
  # batch_colors <- c("#a98de8", "#daccee")
  # batch_table <- setNames(batch_colors, c("batch_1", "batch_2"))
  # 
  # # **Age Category Colors**
  # age_cat_colors <- c("#F5B5AD", "#EF7F70")  # Blue for "low", Vermillion for "high"
  # age_cat_table <- setNames(age_cat_colors, c("low", "high"))
  # 
  # # **Race Colors**
  # race_colors <- c("#9ED4F2", "#56B4E9", "#CDBDF1", "#a98de8")  # Yellow, Bluish Green, Sky Blue, Reddish Purple
  # race_table <- setNames(race_colors, c("White", "Black_African_American", "Asian_Mideast_Indian", "Unknown_Other"))
  
  ha <- HeatmapAnnotation(
    condition = annotation$condition,
    primary_present = annotation$primary_present,
    primary_site = annotation$primary_site,
    sex = annotation$sex,
    batch = annotation$batch,
    age_cat = annotation$age_cat,
    race = annotation$race,
    col = list(
      condition = color_tables$condition_table,
      primary_site = color_tables$primary_site_table,
      sex = color_tables$sex_table,
      primary_present = color_tables$primary_present_table,
      batch = color_tables$batch_table,
      age_cat = color_tables$age_cat_table,
      race = color_tables$race_table
    )
  )
  
  if(comp_set =="AM"){
    ha <- HeatmapAnnotation(
      condition = annotation$condition,
      primary_present = annotation$primary_present,
      path = annotation$ovr_histopath,
      peritoneal_mets = annotation$peritoneal_mets,
      sex = annotation$sex,
      batch = annotation$batch,
      age_cat = annotation$age_cat,
      race = annotation$race,
      col = list(
        condition = color_tables$condition_table,
        primary_present = color_tables$primary_present_table,
        path = color_tables$path_table,
        peritoneal_mets = color_tables$pm_mets_table_am,
        sex = color_tables$sex_table,
        batch = color_tables$batch_table,
        age_cat = color_tables$age_cat_table,
        race = color_tables$race_table
      ),
    annotation_height = unit(c(10,5,5,5,5,5,5,5), "mm"),
    simple_anno_size = unit(4, "mm")
    )
  }

  # Define robust distance function and set robust_val
  robust_val <- 0.05
  use_robust_dist <- TRUE  # Set to TRUE if you want to use robust_dist
  
  robust_dist <- function(x, y) {
    qx = quantile(x, c(robust_val, 1 - robust_val))
    qy = quantile(y, c(robust_val, 1 - robust_val))
    l = x >= qx[1] & x <= qx[2] & y >= qy[1] & y <= qy[2]
    x = x[l]
    y = y[l]
    sqrt(sum((x - y)^2))
  }
  
  # Choose clustering distance based on use_robust_dist
  if(use_robust_dist){
    clustering_distance <- robust_dist
  } else {
    clustering_distance <- "euclidean"
  }
  
  heatmap_out <- Heatmap(
    zscore_matrix_filtered,
    top_annotation = ha,
    col = colorRamp2(c(-3, -2, -1, 0, 1, 2, 3), c("#4676b5", "#82b1d3", "#dbeff6", "white", "#fee395", "#fc9961", "#d73027")),
    clustering_distance_rows = clustering_distance,
    clustering_distance_columns = clustering_distance,
    clustering_method_columns = "ward.D2",
    clustering_method_rows = "ward.D2",
    column_title = heatmap_title,
    show_row_names = FALSE,
    show_column_names = FALSE,
    show_row_dend = FALSE,
    name = "z-score",
    heatmap_legend_param = list(ncol = 2)
  )
  
  svg(paste("./Output/DESeq2/Heatmaps/ward2_sig_heatmap_",contrast_groups[2],contrast_groups[3],"_",ver,"PMnegPMpos.svg", sep = ""),
      width = 6, height = 8) # short: 10, 4.5 or long: 5, 6
  heatmap_out
  dev.off()
  
  png(paste("./Output/DESeq2/Heatmaps/ward2_sig_heatmap_",contrast_groups[2],contrast_groups[3],"_",ver,"PMnegPMpos.png", sep = ""),
      # set width and height
      width = 6, height = 8, units = 'in', res = 300)
  heatmap_out
  dev.off()
  
  heatmap_out
}

### GENERATE BOXPLOTS ENRICHED IN SPECIFIC GENES----
if(FALSE){
  # Load necessary libraries (if not already loaded)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggpubr)   # For stat_compare_means and stat_pvalue_manual
  library(rstatix)  # For statistical tests and integration with ggplot2
  
  # Create color tables
  color_tables <- create_color_tables()
  
  # Define the list of genes you want to plot
  #selected_genes <- c("UGT2B15", "SLC3A1", "CEBPA", "HLA-DMB", "HES5", "FLT4")  # Replace with your genes
  # selected_genes <- c("chr17_2308117_2309225",
  #                     "chr21_39878796_39879742",
  #                     "chr6_112178879_112181580",
  #                     "chr10_76703755_76704867")
  # individual_titles <- c("MNT\nGH17J002389 (1-5kb Upstr.)",
  #                        "ERG\nGH21J038503 (Intron)",
  #                        "FYN\nGH06J111857-59 (Intron)",
  #                        "KAT6B\nGH10J074944 (Intron)")
  selected_genes <- c("HES5","CDH5","THBD","RRAS")
  individual_titles <- c("HES5\n gene body",
                         "CDH5\n gene body",
                         "THBD\n gene body",
                         "RRAS\n gene body")
  
  # Create a named vector for labelling facets
  gene_labels <- setNames(individual_titles, selected_genes)
  
  # Filter the normalized counts for the selected genes
  selected_counts <- normalized_counts_tb %>%
    filter(gene %in% selected_genes) %>%
    pivot_longer(cols = -gene, names_to = "Sample", values_to = "Count")
  
  # Merge with meta data to get condition information
  selected_counts <- selected_counts %>%
    left_join(meta %>% rownames_to_column(var = "Sample"), by = "Sample")
  
  # Perform Wilcoxon test for each gene and save p-values
  stat_test <- selected_counts %>%
    group_by(gene) %>%
    wilcox_test(Count ~ condition) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj")  # Adds a 'p.adj.signif' column based on adjusted p-values
  
  # Add a 'test' column indicating the statistical test used
  stat_test <- stat_test %>%
    mutate(test = "Wilcoxon")  # All tests are Wilcoxon in this script
  
  # Add y.position for label placement
  stat_test <- stat_test %>%
    add_y_position(fun = "max", scales = "free_y", step.increase = 0.1)
  
  # Remove 'groups' column before saving to CSV
  stat_test_savable <- stat_test %>%
    select(gene, p, p.adj, p.adj.signif, y.position, test)  # Select relevant columns
  
  # Save p-values and test information to CSV
  write.csv(stat_test_savable, paste("./Output/DESeq2/specific_gene_box/", 
                                     contrast_groups[[2]][1], contrast_groups[[2]][2], 
                                     "_selected_genes_p_values.csv", sep = "_"), 
            row.names = FALSE)
  
  # Compute median values for each gene and condition
  medians_df <- selected_counts %>%
    group_by(gene, condition) %>%
    summarize(median_count = median(Count), .groups = 'drop')
  
  # Create boxplots for each gene, colored by condition, with statistical test results and median annotations
  p3 <- ggplot(selected_counts, aes(x = condition, y = Count)) +
    geom_boxplot(aes(fill = condition), outlier.shape = NA, alpha = 0.6) +  # Map fill within geom_boxplot
    geom_jitter(aes(color = condition), width = 0.2, alpha = 0.8, show.legend = FALSE) +  # Map color within geom_jitter and hide its legend
    facet_wrap(~ gene, scales = "free_y", labeller = labeller(gene = gene_labels)) +
    theme_bw() +
    xlab(NULL) +  # Remove x-axis title
    ylab("Normalized Counts") +
    ggtitle("Expression of Selected Genes") +
    theme(
      legend.position = "right",             # Position the legend on the right
      axis.text.x = element_blank(),         # Remove x-axis tick labels
      axis.ticks.x = element_blank()         # Remove x-axis ticks
    ) +
    # Add median annotations
    geom_text(data = medians_df, 
              aes(x = condition, y = median_count, label = round(median_count, 1)), 
              color = "white", 
              size = 3, 
              vjust = 0.5) +
    # Add statistical test annotations
    stat_pvalue_manual(
      stat_test,
      label = "p.adj.signif",    # Correct label variable
      tip.length = 0,
      label.x = 1.5              # Centers the annotation between the two conditions
    ) +
    # Define consistent colors for fill and jitter using color_tables
    scale_fill_manual(values = color_tables$condition_table[c("AMpos", "AMneg")]) +
    scale_color_manual(values = color_tables$condition_table[c("AMpos", "AMneg")]) +
    guides(fill = guide_legend(title = "Condition"))  # Add a single legend for fill
  
  # Display the plot
  print(p3)
  
  # Ensure the output directory exists
  if(!dir.exists("./Output/DESeq2/specific_gene_box/")) {
    dir.create("./Output/DESeq2/specific_gene_box/", recursive = TRUE)
  }
  
  # Save the plot to SVG
  svg_filename <- paste("./Output/DESeq2/specific_gene_box/", 
                        contrast_groups[[2]][1], contrast_groups[[2]][2], 
                        "_selected_genes_boxplots_with_stats.svg", sep = "_")
  svg(svg_filename, width = 6, height = 5)
  print(p3)
  dev.off()
  
  # Save the plot to PNG
  png_filename <- paste("./Output/DESeq2/specific_gene_box/", 
                        contrast_groups[[2]][1], contrast_groups[[2]][2], 
                        "_selected_genes_boxplots_with_stats.png", sep = "_")
  png(png_filename, width = 6, height = 5, units = 'in', res = 300)
  print(p3)
  dev.off()
}

### GENERATE BOXPLOTS ENRICHED IN SIG REGIONS:----
if(FALSE){
  # Load required libraries
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggsignif)
  library(gridExtra)
  library(broom)  # For tidy output of test results
  
  # Read the CSV file
  all_results <- read.csv("~/5hmC-Pathway-Analysis/Output/DESeq2/Results/all_results_annotated_PM_negativePM_positive_pval0.001_lfc0.2_11_23_2024_peaks_nocombat_PM_healthy.txt", sep = "\t")
  
  # Add column called "region_type" to all_results
  all_results <- all_results %>%
    mutate(region_type = ifelse(log2FoldChange > 0, "PM", "No_PM"))
  
  # Create annotation dataframe
  annotation <- meta_all %>%
    dplyr::select(condition, primary_present, ovr_histopath, chemo_6weeks, primary_site)
  
  # Prepare data for boxplot
  boxplot_data <- as.data.frame(norm_sig_all)
  boxplot_data$region <- rownames(boxplot_data)
  boxplot_data_long <- tidyr::pivot_longer(boxplot_data, -region, names_to = "sample", values_to = "normalized_count")
  
  # Add region_type column to boxplot_data_long by merging on region and feature
  boxplot_data_long <- boxplot_data_long %>%
    left_join(all_results %>% dplyr::select(feature, region_type), by = c("region" = "feature"))
  
  # Add condition information
  boxplot_data_long$condition <- annotation$condition[match(boxplot_data_long$sample, rownames(annotation))]
  
  # Calculate fold change values
  boxplot_data_long <- boxplot_data_long %>%
    group_by(region, region_type) %>%
    mutate(
      baseline_mean = mean(normalized_count[condition == "HEALTHY"]),
      fold_change = (normalized_count / baseline_mean - 1) * 100  # Convert to percentage
    ) %>%
    ungroup()
  
  # Define colors for conditions
  condition_colors <- color_tables$condition_table
  
  # Calculate the number of regions in each region type
  n_regions <- boxplot_data_long %>%
    group_by(region_type) %>%
    summarise(n_regions = n_distinct(region), .groups = 'drop')
  
  # Create a named vector for easier lookup
  n_regions_vector <- setNames(n_regions$n_regions, n_regions$region_type)
  
  # Create a custom labeller function to include (n=n_regions) in subtitles
  region_type_labeller <- function(value) {
    n <- n_regions_vector[value]
    return(paste(value, " (n=", n, " regions)", sep = ""))
  }
  
  # Function to calculate the whiskers for each group
  get_whiskers <- function(y) {
    q1 <- quantile(y, 0.25, na.rm = TRUE)
    q3 <- quantile(y, 0.75, na.rm = TRUE)
    iqr <- q3 - q1
    lower_whisker <- max(min(y, na.rm = TRUE), q1 - 1.5 * iqr)
    upper_whisker <- min(max(y, na.rm = TRUE), q3 + 1.5 * iqr)
    return(c(lower_whisker, upper_whisker))
  }
  
  # Function to calculate boxplot statistics
  calculate_boxplot_stats <- function(data) {
    stats <- data %>%
      summarise(
        mean = mean(fold_change, na.rm = TRUE),
        median = median(fold_change, na.rm = TRUE),
        q1 = quantile(fold_change, 0.25, na.rm = TRUE),
        q3 = quantile(fold_change, 0.75, na.rm = TRUE),
        lower_whisker = min(fold_change, na.rm = TRUE),
        upper_whisker = max(fold_change, na.rm = TRUE)
      )
    return(stats)
  }
  
  # Initialize list to store results
  results <- list()
  
  # Create separate plots for each region_type and calculate statistics
  plot_list <- lapply(unique(boxplot_data_long$region_type), function(rt) {
    data_subset <- subset(boxplot_data_long, region_type == rt)
    
    # Calculate boxplot statistics
    boxplot_stats <- data_subset %>%
      group_by(condition) %>%
      do(calculate_boxplot_stats(.))
    
    # Perform one-way ANOVA
    anova_result <- aov(fold_change ~ condition, data = data_subset)
    anova_summary <- summary(anova_result)
    
    # Perform Tukey's HSD post-hoc test
    tukey_result <- TukeyHSD(anova_result)
    tukey_df <- as.data.frame(tukey_result$condition)
    tukey_df$comparison <- rownames(tukey_df)
    tukey_df$region_type <- rt
    
    # Split 'comparison' into 'group1' and 'group2'
    tukey_df <- tukey_df %>%
      separate(comparison, into = c("group1", "group2"), sep = "-")
    
    # Store results in a list
    results[[rt]] <- list(
      boxplot_stats = boxplot_stats,
      anova_summary = anova_summary,
      tukey_results = tukey_df
    )
    
    # Create pairwise comparisons for plot
    comparisons_from_tukey <- mapply(c, tukey_df$group1, tukey_df$group2, SIMPLIFY = FALSE)
    
    # Extract adjusted p-values
    pvalues <- tukey_df$`p adj`
    
    # Calculate the whiskers for each condition
    whiskers <- do.call(rbind, tapply(data_subset$fold_change, data_subset$condition, get_whiskers))
    lower_whisker <- min(whiskers[,1])
    upper_whisker <- max(whiskers[,2])
    whisker_range <- upper_whisker - lower_whisker
    
    # Calculate y positions for significance bars
    y_start <- upper_whisker - 0.25 * whisker_range/2
    
    # Calculate y-axis limits
    y_min <- -60
    y_max <- 90
    
    # Calculate breaks at 50% intervals
    break_interval <- 50
    min_break <- floor(y_min / break_interval) * break_interval
    max_break <- ceiling(y_max / break_interval) * break_interval
    breaks <- seq(min_break, max_break, by = break_interval)
    
    plot <- ggplot(data_subset, aes(x = condition, y = fold_change, fill = condition)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      theme_minimal(base_size = 20) +
      labs(title = region_type_labeller(rt),
           x = "Condition",
           y = "Fold Change (%)") +
      geom_signif(
        comparisons = comparisons_from_tukey,
        p_values = pvalues,
        step_increase = 0.02,
        map_signif_level = TRUE,
        y_position = rep(y_start, length(comparisons_from_tukey)),
        tip_length = 0
      ) +
      scale_fill_manual(values = condition_colors) +
      theme(
        legend.position = "right",
        strip.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        legend.background = element_rect(fill = "white", color = NA),
        axis.text.x = element_blank()  # This removes x-axis labels
      ) +
      scale_y_continuous(
        limits = c(y_min, y_max),
        breaks = breaks,
        labels = function(x) paste0(x, "%"),
        expand = c(0, 0)  # This removes padding
      ) +
      # **Add Median Labels in White**
      geom_text(
        data = boxplot_stats,
        aes(x = condition, y = median, label = sprintf("%.1f", median)),
        vjust = -0.5,            # Adjust vertical position above the median line
        size = 2.5,                # Adjust text size as needed
        color = "white",         # Set text color to white
        fontface = "bold"        # Optional: make text bold for better visibility
      )
    
    return(list(plot = plot, stats = results[[rt]]))
  })
  
  # Combine all results
  all_boxplot_stats <- bind_rows(lapply(plot_list, function(x) x$stats$boxplot_stats))
  all_anova_results <- bind_rows(lapply(plot_list, function(x) broom::tidy(x$stats$anova_summary[[1]])))
  all_tukey_results <- bind_rows(lapply(plot_list, function(x) x$stats$tukey_results))
  
  # Save results to CSV
  write.csv(all_boxplot_stats, "./Output/DESeq2/enrichment_boxplots/boxplot_statistics.csv", row.names = FALSE)
  write.csv(all_anova_results, "./Output/DESeq2/enrichment_boxplots/anova_results.csv", row.names = FALSE)
  write.csv(all_tukey_results, "./Output/DESeq2/enrichment_boxplots/tukey_results.csv", row.names = FALSE)
  
  # Extract plots from plot_list
  plots <- lapply(plot_list, function(x) x$plot)
  
  # Arrange the plots in a grid
  boxplot <- grid.arrange(grobs = plots, ncol = 2)
  
  # Save boxplot
  ggsave("./Output/DESeq2/enrichment_boxplots/enrichment_boxplot_by_region_type.png", 
         boxplot, width = 9, height = 5, dpi = 300)
  
  ggsave("./Output/DESeq2/enrichment_boxplots/enrichment_boxplot_by_region_type.svg", 
         boxplot, width = 9, height = 5)
}

### GENERATE BOXPLOTS BY REGION TYPE:----
if(FALSE){
  # Read the annotated results file
  annotated_results <- read.csv("~/5hmC-Pathway-Analysis/Output/DESeq2/Results/all_results_annotated_PM_positivePM_negative_pval0.001_lfc0.2_10_13_2024_hmr_combat.txt", sep="\t")
  
  # Create a lookup table for region types
  region_type_lookup <- annotated_results %>%
    select(feature, annot.type) %>%
    distinct()
  
  # Add region type to normalized_counts_all
  normalized_counts_df <- as.data.frame(normalized_counts_all) %>%
    rownames_to_column("feature") %>%
    left_join(region_type_lookup, by = "feature")
  
  # Melt the dataframe for easier plotting
  melted_counts <- normalized_counts_df %>%
    pivot_longer(cols = -c(feature, annot.type), names_to = "sample", values_to = "count")
  
  # Join with metadata
  melted_counts_with_meta <- melted_counts %>%
    left_join(meta_all, by = c("sample" = "row.names"))
  
  # Calculate mean counts per region type and condition
  mean_counts <- melted_counts_with_meta %>%
    group_by(annot.type, condition) %>%
    summarize(mean_count = mean(count, na.rm = TRUE))
  
  # Create the plot
  ggplot(mean_counts, aes(x = annot.type, y = mean_count, fill = condition)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Region Type", y = "Mean Normalized Count", 
         title = "Enrichment by Region Type",
         fill = "Condition") +
    scale_fill_brewer(palette = "Set1")
  
  # If you want to add error bars, you can calculate standard error and add them to the plot
  counts_with_se <- melted_counts_with_meta %>%
    group_by(annot.type, condition) %>%
    summarize(
      mean_count = mean(count, na.rm = TRUE),
      se = sd(count, na.rm = TRUE) / sqrt(n())
    )
  
  ggplot(counts_with_se, aes(x = annot.type, y = mean_count, fill = condition)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_errorbar(aes(ymin = mean_count - se, ymax = mean_count + se), 
                  position = position_dodge(width = 0.9), width = 0.25) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Region Type", y = "Mean Normalized Count", 
         title = "Enrichment by Region Type with Standard Error",
         fill = "Condition") +
    scale_fill_brewer(palette = "Set1")
  
}

## Generate barplots----
if(FALSE){
  # Load required libraries
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(rstatix)
  
  # Read the CSV file
  all_results <- read.csv("~/5hmC-Pathway-Analysis/Output/DESeq2/Results/all_results_annotated_PM_positivePM_negative_ntop 2000 _ nocombat_10per082024_seed126.txt", sep = "\t")
  
  # Prepare the annotation data
  annotation <- meta_all %>%
    dplyr::select(condition, primary_present, ovr_histopath, chemo_6weeks, primary_site)
  
  # Filter out 'lncrna_gencode' and 'NA' categories from all_results
  all_results <- all_results %>%
    filter(!annot.type %in% c("hg19_lncrna_gencode", "NA")) %>%
    filter(!is.na(annot.type))
  
  # Prepare data for analysis
  barplot_data <- as.data.frame(normalized_counts_all)
  barplot_data$region <- rownames(barplot_data)
  barplot_data_long <- tidyr::pivot_longer(barplot_data, -region, names_to = "sample", values_to = "enrichment_score")
  
  # Add condition information
  barplot_data_long$condition <- annotation$condition[match(barplot_data_long$sample, rownames(annotation))]
  
  # Add region type information
  region_types <- all_results %>%
    dplyr::select(feature, annot.type) %>%
    dplyr::rename(region = feature, region_type = annot.type)
  
  barplot_data_long <- barplot_data_long %>%
    left_join(region_types, by = "region") %>%
    filter(!is.na(region_type))  # Remove any rows where region_type is NA
  
  
  # Update region type names
  region_type_mapping <- c(
    "hg19_genes_1to5kb" = "1to5kb",
    "hg19_genes_promoters" = "promoters",
    "hg19_genes_5UTRs" = "5UTRs",
    "hg19_genes_introns" = "introns",
    "hg19_genes_exons" = "exons",
    "hg19_genes_3UTRs" = "3UTRs"
  )
  
  barplot_data_long <- barplot_data_long %>%
    mutate(region_type = region_type_mapping[region_type]) %>%
    filter(!is.na(region_type))  # Remove any rows where region_type is NA after mapping
  
  # Set the desired order for region types and conditions
  region_type_order <- c("1to5kb", "promoters", "5UTRs", "introns", "exons", "3UTRs")
  
  # Sort barplot_data_long by region_type in order of region_type_order
  barplot_data_long <- barplot_data_long %>%
    mutate(region_type = factor(region_type, levels = region_type_order))
  
  # Create a summary of the data for plotting
  summary_data <- barplot_data_long %>%
    group_by(region_type, condition) %>%
    summarise(
      mean_enrichment = mean(enrichment_score, na.rm = TRUE),
      sd_enrichment = sd(enrichment_score, na.rm = TRUE),
      n = n(),
      se_enrichment = sd_enrichment / sqrt(n),
      ci_lower = mean_enrichment - qt(0.975, df = n - 1) * se_enrichment,
      ci_upper = mean_enrichment + qt(0.975, df = n - 1) * se_enrichment,
      .groups = 'drop'
    )
  
  # Count the number of regions for each region type
  region_counts <- barplot_data_long %>%
    group_by(region_type) %>%
    summarise(n_regions = n_distinct(region), .groups = 'drop')
  
  # Add region counts to the x-axis labels
  region_labels <- paste0(region_type_order, "\n(n=", region_counts$n_regions, ")")
  
  condition_order <- c("HEALTHY", "PM_negative", "PM_positive")
  
  # Perform pairwise comparisons using the full dataset
  stat.test <- barplot_data_long %>%
    group_by(region_type) %>%
    pairwise_t_test(
      enrichment_score ~ condition, 
      p.adjust.method = "bonferroni",
      comparisons = list(
        c("HEALTHY", "PM_negative"),
        c("HEALTHY", "PM_positive"),
        c("PM_negative", "PM_positive")
      )
    ) %>%
    add_xy_position(x = "region_type", fun="mean", dodge = 0.8)
  
  # Define colors for conditions
  condition_colors <- color_tables$condition_table
  
  # install ggbreak package
  #install.packages("ggbreak")
  library(ggbreak)
  
  
  # Create the plot using summary_data
  ggplot <- ggplot(summary_data, aes(x = region_type, y = mean_enrichment, fill = condition)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                  position = position_dodge(width = 0.8), width = 0.25) +
    scale_fill_manual(values = condition_colors) +
    labs(title = "5hmC Enrichment Scores by Region Type",
         x = "Region Type",
         y = "Mean 5hmC Enrichment Score",
         fill = "Condition") +
    theme_pubr() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top",
      panel.grid.major.y = element_line(color = "grey90")
    ) +
    scale_x_discrete(limits = region_type_order, labels = region_labels) +
    scale_y_continuous(
      breaks = c(0, 5, seq(25, 55, by = 5)),  # Specify the breaks you want to show
      labels = c("0", "5", seq(25, 55, by = 5)),  # Specify the labels for these breaks
      expand = expansion(mult = c(0, 0.1))
    )+
    scale_y_break(c(5, 25))
  
  # Add significance levels using stat_pvalue_manual
  ggplot <- ggplot +
    stat_pvalue_manual(
      stat.test,
      label = "p.adj.signif",
      xmin="xmin",
      xmax="xmax",
      inherit.aes = FALSE
    )
  
  # Print the plot
  print(ggplot)
  
  
  # Create an empty plot for the custom significance legend
  sig_legend_plot <- ggplot(data = data.frame(x = 1, y = 1), aes(x = x, y = y)) + 
    theme_void() +  # Remove all grid, axes, and background
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 10, face = "bold"),  # Customize title
      legend.text = element_text(size = 8)  # Customize text size
    ) +
    scale_color_manual(
      name = "Significance Levels",  # Title for the legend
      values = c("ns" = "black", "*" = "black", "**" = "black", "***" = "black"),  # Colors
      labels = c("ns (p ≥ 0.05)", "* (p < 0.05)", "** (p < 0.01)", "*** (p < 0.001)"),  # Labels
      guide = guide_legend(override.aes = list(size = 2))  # Adjust point size in legend
    ) +
    # Dummy points to create legend entries
    geom_point(aes(color = "ns")) +
    geom_point(aes(color = "*")) +
    geom_point(aes(color = "**")) +
    geom_point(aes(color = "***"))
  
  # Print the empty plot with the custom legend
  print(sig_legend_plot)
  
  # save significance legend plot
  ggsave("./Output/DESeq2/enrichment_boxplots/5hmc_enrichment_grouped_barplot_with_significance_legend.png", 
         sig_legend_plot, width = 5, height = 1)
  # save as svg
  ggsave("./Output/DESeq2/enrichment_boxplots/5hmc_enrichment_grouped_barplot_with_significance_legend.svg", 
         sig_legend_plot, width = 5, height = 1)
  
  # Save barplot
  ggsave("./Output/DESeq2/enrichment_boxplots/5hmc_enrichment_grouped_barplot_with_significance.png", 
         ggplot, width = 10, height = 5, dpi = 300)
  
  ggsave("./Output/DESeq2/enrichment_boxplots/5hmc_enrichment_grouped_barplot_with_significance.svg", 
         ggplot, width = 10, height = 5)
  
  # Identify list columns
  list_cols <- sapply(stat.test, is.list)
  
  # Flatten list columns
  stat.test_flat <- stat.test %>%
    mutate(across(where(is.list), ~ sapply(.x, paste, collapse = ";"))) %>%
    mutate(across(where(is.list), as.character))
  
  # Write the flattened tibble to CSV
  write.csv(stat.test_flat, "./Output/DESeq2/enrichment_boxplots/5hmc_enrichment_grouped_barplot_with_significance_stat_test.csv", row.names = FALSE)
  
  # save summary_data to csv
  write.csv(summary_data, "./Output/DESeq2/enrichment_boxplots/5hmc_enrichment_grouped_barplot_with_significance_summary_data.csv", row.names = FALSE)
  
}

###GENERATE PCA PLOT----
if(output_PCA == 1) {
  #See details for below operation in lesson 3 of DGE workshop
  rld <- vst(dds_all, blind=TRUE) #dds_all for healthy, pmpos, pm neg, or dds pmpos pmneg
  assay(rld) <- t(scale(t(assay(rld)), center = TRUE, scale = TRUE))
  
  res_genes <- row.names(res_table)
  sig_genes <-  sig$gene
  cols <- c(colnames(meta_all[,!(names(meta_all) %in% c("condition","age"))]))
  
  # Define custom colors
  # Updated custom colors using the new color scheme
  
  condition_colors <- c(
    color_tables$condition_table,
    color_tables$pm_mets_table_am, # pm_mets_table_am or pm_mets_table
    color_tables$batch_table,
    color_tables$sex_table,
    color_tables$primary_site_table,
    color_tables$race_table,
    color_tables$primary_present_table,
    color_tables$path_table,
    color_tables$age_cat_table
  )
  
  # Decide whether to use sig_genes or res_genes for PCA
  considered_factors <- "sig" # or sig for significant genes, res for all genes
  if(considered_factors == "sig"){
    pca_genes <- sig_genes
  } else {
    pca_genes <- res_genes
  }

  pca_genes_type <- if(identical(pca_genes, sig_genes)) "sig_genes" else "res_genes"
  
  for(each in c("condition","age_cat","peritoneal_mets","batch","race","sex","primary_present","primary_site","ovr_histopath")){ #(each in c(cols,"condition")){
    # Construct title
    if(cutoff_type == 0){
      PCA_title <- paste(each," ,", contrast_groups[2],"/",contrast_groups[3],"padj <",padj.cutoff)
    }
    
    if(cutoff_type == 1){
      PCA_title <- paste(each," ,", contrast_groups[2],"/",contrast_groups[3],"pvalue <",pvalue.cutoff,"|lfc|>",lfc.cutoff)
    }
    
    #In the below replace sig_genes with res_genes if you want to perform PCA analysis on all genes rather than just on significant genes.
    plotPCA_labels <- plotPCA(rld[pca_genes,], intgroup = c(each), ntop = 500) #rownames(c(row.names(meta)))
    #reverse order of labels to match heatmap labelling
    # Extract sample names for labeling
    plotPCA_labels$data$sample <- rownames(plotPCA_labels$data)
    
    # Modify the plot to have a white background with no borders and no axis lines
    plotPCA_labels <- plotPCA_labels + 
      theme_minimal() +  # Use theme_minimal() as a base
      theme(
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.background = element_rect(fill = "white", colour = NA),
        legend.box.background = element_rect(fill = "white", colour = NA)
      ) + 
      stat_ellipse(level=0.9) +
      scale_color_manual(values = condition_colors)  # Apply custom colors
    
    png(paste("./Output/DESeq2/PCA/",considered_factors,"_PCA_",contrast_groups[2],contrast_groups[3],"_",ver,each,".png", sep = ""))
    print({plotPCA_labels + 
        ggtitle(PCA_title) + 
        stat_ellipse(level=0.9)
    }, cex.main=3)
    Sys.sleep(2)
    dev.off()
    
    svg(paste("./Output/DESeq2/PCA/",considered_factors,"_PCA_",contrast_groups[2],contrast_groups[3],"_",ver,each,".svg", sep = ""))
    print({plotPCA_labels + 
        ggtitle(PCA_title) + 
        stat_ellipse(level=0.9)
    }, cex.main=3)
    dev.off()
  }
}

### SAVE CONFIG TABLES TO TEXT FILES----
if(output_config == 1){
  # Create folders
  
  analysis_date <- Sys.time()
  
  subset_conditions_str <- paste(capture.output(dput(subset_conditions)), collapse="")
  design_formula_str <- paste(deparse(design_formula), collapse="")
  
  # General Settings
  config_general <- list(
    "Analysis date and time" = as.character(analysis_date),
    "Comparison set" = as.character(comp_set),
    "Subset conditions" = subset_conditions_str,
    "Counts file name" = as.character(counts_name),
    "Conditions file name" = as.character(meta_name),
    "Version text" = as.character(ver.text),
    "Output version" = as.character(ver),
    "Design formula" = design_formula_str,
    "Region type" = as.character(region_type),
    "Number of samples" = as.character(ncol(counts_data)),
    # if gene filtering was applied
    "Contrast groups" = paste(contrast_groups, collapse=", ")
  )
  
  if(genebodies_type == TRUE){
    config_general[["Gene filtering"]] <- "Gene bodies only"
    config_general[["Number of genes before filtering"]] <- as.character(nrow(counts_data) + sum(!keep))
    config_general[["Number of genes after filtering"]] <- as.character(nrow(counts_data))
    config_general[["Prefiltering settings"]] <- paste("Rows kept: ", nrow(counts_data), " out of ", nrow(counts_data) + sum(!keep), " (", round(nrow(counts_data) / (nrow(counts_data) + sum(!keep)) * 100, 2), "%)")
  }
  
  # Volcano Plot Parameters
  config_volcano <- list(
    "Output volcano plot" = as.character(output_volcano),
    "Cutoff type" = as.character(cutoff_type),
    "padj cutoff" = as.character(padj.cutoff),
    "pvalue cutoff" = as.character(pvalue.cutoff),
    "log2 fold change cutoff" = as.character(lfc.cutoff)
  )
  
  # Heatmap Parameters
  config_heatmap <- list(
    "Output heatmap" = as.character(output_heatmap),
    "Heatmap title" = as.character(heatmap_title),
    "Number of significant genes" = as.character(nsig),
    "Heatmap clustering method" = "ward.D2",
    "Heatmap distance metric" = if(use_robust_dist) paste0("robust_dist (robust_val = ", robust_val, ")") else "euclidean",
    "Robust distance used" = as.character(use_robust_dist),
    "Robust_val" = as.character(robust_val)
  )
  
  if(output_PCA == 1){
    # PCA Plot Parameters
    config_pca <- list(
      "Output PCA" = as.character(output_PCA),
      "PCA ntop" = "500",
      "PCA transformation" = "VST with z-score normalization",
      "PCA groups" = paste(cols, collapse=", "),
      "PCA genes used" = as.character(pca_genes_type)
    )
  } else {
    config_pca <- list(
      "Output PCA" = as.character(output_PCA)
    )
  }
  
  
  # UMAP Parameters
  config_umap <- list(
    "Output UMAP" = as.character(output_UMAP)
    # Add UMAP-specific parameters here if applicable
  )
  
  # Other Outputs
  config_other_outputs <- list(
    "Output results tables" = as.character(output_results_tables),
    "Output config" = as.character(output_config),
    "Output cross-tabs" = as.character(output_xtabs),
    "Output normalized counts" = as.character(output_norm_counts)
  )
  
  # Combine all configurations into one list with sections
  config_sections <- list(
    "General Settings" = config_general,
    "Volcano Plot Parameters" = config_volcano,
    "Heatmap Parameters" = config_heatmap,
    "PCA Plot Parameters" = config_pca,
    "UMAP Parameters" = config_umap,
    "Other Outputs" = config_other_outputs
  )
  
  # Function to write each section to the config file
  config_file_name <- paste0("./Output/DESeq2/Config/config_",contrast_groups[2],contrast_groups[3],"_",ver,".txt")
  con <- file(config_file_name, "w")
  for(section_name in names(config_sections)){
    writeLines(paste0("### ", section_name, " ###"), con)
    section <- config_sections[[section_name]]
    section_lines <- sapply(names(section), function(x) paste(x, ":", as.character(section[[x]])))
    writeLines(section_lines, con)
    writeLines("", con)  # Add a blank line between sections
  }
  close(con)
  
  # Save session info
  session_info <- capture.output(sessionInfo())
  session_info_file_name <- paste0("./Output/DESeq2/Config/session_info_",contrast_groups[2],contrast_groups[3],"_",ver,".txt")
  writeLines(session_info, con=session_info_file_name)
}

### GREAT Enrichment Analysis ----
if(TRUE){
  # Load necessary package
  library(GenomicRanges)
  
  # Split the 'gene' column into 'chr', 'start', and 'end'
  res_table_tb <- transform(
    res_table_tb,
    chr = sapply(strsplit(gene, "_"), `[`, 1),
    start = as.numeric(sapply(strsplit(gene, "_"), `[`, 2)),
    end = as.numeric(sapply(strsplit(gene, "_"), `[`, 3))
  )
  
  # Create background GRanges object with all rows
  background_grange <- GRanges(
    seqnames = res_table_tb$chr,
    ranges = IRanges(start = res_table_tb$start, end = res_table_tb$end),
    baseMean = res_table_tb$baseMean,
    log2FoldChange = res_table_tb$log2FoldChange,
    lfcSE = res_table_tb$lfcSE,
    stat = res_table_tb$stat,
    pvalue = res_table_tb$pvalue,
    padj = res_table_tb$padj
  )
  
  # Define the configurable number of rows to include
  num_rows <- 1000 # Adjust as needed
  
  # Create ampos_grange with rows having the greatest 'stat' values
  ampos_df <- res_table_tb[order(-res_table_tb$stat), ][1:num_rows, ]
  ampos_grange <- GRanges(
    seqnames = ampos_df$chr,
    ranges = IRanges(start = ampos_df$start, end = ampos_df$end),
    baseMean = ampos_df$baseMean,
    log2FoldChange = ampos_df$log2FoldChange,
    lfcSE = ampos_df$lfcSE,
    stat = ampos_df$stat,
    pvalue = ampos_df$pvalue,
    padj = ampos_df$padj
  )
  
  # Create amneg_grange with rows having the smallest 'stat' values
  amneg_df <- res_table_tb[order(res_table_tb$stat), ][1:num_rows, ]
  amneg_grange <- GRanges(
    seqnames = amneg_df$chr,
    ranges = IRanges(start = amneg_df$start, end = amneg_df$end),
    baseMean = amneg_df$baseMean,
    log2FoldChange = amneg_df$log2FoldChange,
    lfcSE = amneg_df$lfcSE,
    stat = amneg_df$stat,
    pvalue = amneg_df$pvalue,
    padj = amneg_df$padj
  )
  
  amall_df <- res_table_tb[order(abs(res_table_tb$stat), decreasing = TRUE), ][1:(2*num_rows), ]
  amall_grange <- GRanges(
    seqnames = amall_df$chr,
    ranges = IRanges(start = amall_df$start, end = amall_df$end),
    baseMean = amall_df$baseMean,
    log2FoldChange = amall_df$log2FoldChange,
    lfcSE = amall_df$lfcSE,
    stat = amall_df$stat,
    pvalue = amall_df$pvalue,
    padj = amall_df$padj
  )
  
  # Load the required library
  # install rGREAT from bioconductor
  #BiocManager::install("rGREAT")
  #library(devtools)
  #install_github("jokergoo/rGREAT")
  
  library(rGREAT)
  # print what version i am using
  packageVersion("rGREAT")
  
  # Set the analysis configuration
  gene_set_collection <- "msigdb:C2:CP:KEGG"  # or GO:BP or msigdb:C2:CP:REACTOME or msigdb:C2:CP:KEGG
  extend_from_setting <- "gene"  # Extend from the full gene
  mode_setting = "basalPlusExt" # or 'basalPlusExt', 'twoClosest' and 'oneClosest'.
  
  # Perform enrichment analysis on the ampos_grange
  res_ampos <- great(
    ampos_grange,
    gene_set_collection,
    min_gene_set_size = 15,
    background = background_grange,
    "hg19"#,  # Adjust genome version to match your data (e.g., "hg38" or "hg19")
    #mode = mode_setting,
    #extend_from = extend_from_setting
  )
  
  # Retrieve the enrichment table for ampos_grange
  enrichment_ampos <- getEnrichmentTable(res_ampos)
  # sort by p_adjust_hyper smallest to larget
  enrichment_ampos <- enrichment_ampos[order(enrichment_ampos$p_adjust_hyper),]
  
  # Perform enrichment analysis on the amneg_grange
  res_amneg <- great(
    amneg_grange,
    gene_set_collection,
    min_gene_set_size = 15,
    "hg19",  # Adjust genome version to match your data
    background = background_grange#,
    #mode = mode_setting,
    #extend_from = extend_from_setting
  )
  
  # Retrieve the enrichment table for amneg_grange
  enrichment_amneg <- getEnrichmentTable(res_amneg)
  # sort by p_adjust_hyper smallest to larget
  enrichment_amneg <- enrichment_amneg[order(enrichment_amneg$p_adjust_hyper),]

  
  # Perform enrichment analysis on the amall_grange
  res_amall <- great(
    amall_grange,
    gene_set_collection,
    min_gene_set_size = 15,
    "hg19",  # Adjust genome version to match your data
    background = background_grange
    #mode = "basalPlusExt"#,
    #extend_from = extend_from_setting
  )
  
  enrichment_amall <- getEnrichmentTable(res_amall)
  # sort by p_adjust_hyper smallest to larget
  enrichment_amall <- enrichment_amall[order(enrichment_amall$p_adjust_hyper),]
  print("Enrichment results for amall_grange:")
  print(head(enrichment_amall))
  
  print("Enrichment results for ampos_grange:")
  print(head(enrichment_ampos))
  
  print("Enrichment results for amneg_grange:")
  print(head(enrichment_amneg))
  
  ## SAVE ALL CSVs to Output/DESeq2/GSEA_regions
  write.csv(enrichment_ampos, paste0("./Output/DESeq2/GSEA_regions/enrichment_ampos_",contrast_groups[2],contrast_groups[3],"_",ver,".csv"), row.names = FALSE)
  write.csv(enrichment_amneg, paste0("./Output/DESeq2/GSEA_regions/enrichment_amneg_",contrast_groups[2],contrast_groups[3],"_",ver,".csv"), row.names = FALSE)
  write.csv(enrichment_amall, paste0("./Output/DESeq2/GSEA_regions/enrichment_amall_",contrast_groups[2],contrast_groups[3],"_",ver,".csv"), row.names = FALSE)
}

### Visualize rGREAT rresult ----
if(FALSE){
  # Load necessary libraries
  library(rGREAT)
  library(msigdbr)
  library(clusterProfiler)
  library(KEGGREST)
  library(dplyr)
  library(tidyr)
  
  # Assuming 'enrichment_amneg' is your enrichment table from previous steps
  
  # Step 1: Get MSigDB KEGG pathways
  msigdb_kegg <- msigdbr(
    species = "Homo sapiens",
    category = "C2",
    subcategory = "CP:KEGG"
  )
  
  # Extract all genes involved in KEGG pathways
  all_kegg_genes <- unique(msigdb_kegg$human_gene_symbol)
  
  # Step 2: Map 'id's to KEGG IDs and Descriptions
  
  # Create a mapping between gs_name and KEGG IDs from gs_exact_source
  pathway_mapping <- msigdb_kegg %>%
    dplyr::select(gs_name, gs_exact_source) %>%
    distinct() %>%
    mutate(KEGG_ID = sub("path:", "", gs_exact_source))
  
  # Merge the enrichment table with the pathway mapping
  enrichment_amneg_mapped <- enrichment_ampos %>%
    left_join(pathway_mapping, by = c("id" = "gs_name"))
  
  # Get KEGG pathway descriptions
  kegg_descriptions <- keggList("pathway", "hsa")
  kegg_df <- data.frame(
    KEGG_ID = sub("path:", "", names(kegg_descriptions)),
    Description = sub(" - Homo sapiens \\(human\\)", "", kegg_descriptions),
    stringsAsFactors = FALSE
  )
  
  # Merge descriptions into the enrichment table
  enrichment_amneg_mapped <- enrichment_amneg_mapped %>%
    left_join(kegg_df, by = "KEGG_ID")
  
  # Step 3: Get genes associated with each pathway
  # Group msigdb_kegg to get genes per pathway
  pathway_genes <- msigdb_kegg %>%
    group_by(gs_name) %>%
    summarize(genes_in_pathway = list(unique(human_gene_symbol)))
  
  # Get the genes associated with your input regions
  region_gene_assoc <- getRegionGeneAssociations(res_amneg)
  input_genes <- unique(unlist(region_gene_assoc$annotated_genes))
  
  # Filter input genes to those present in KEGG pathways
  input_genes_in_kegg <- intersect(input_genes, all_kegg_genes)
  total_input_genes_in_kegg <- length(input_genes_in_kegg)
  
  # Step 4: Compute GeneRatio and BgRatio
  # Total number of genes in the input set
  total_input_genes <- length(input_genes)
  
  # For the background, get associated genes
  region_gene_assoc_bg <- getRegionGeneAssociations(
    great(
      background_grange,
      gene_set_collection,
      "hg19",
      min_gene_set_size = 15,
      mode = mode_setting,
      extend_from = extend_from
    )
  )
  background_genes <- unique(unlist(region_gene_assoc_bg$annotated_genes))
  background_genes_in_kegg <- intersect(background_genes, all_kegg_genes)
  total_background_genes_in_kegg <- length(background_genes_in_kegg)
  
  # For each pathway, find overlapping genes with input genes
  enrichment_amneg_mapped <- enrichment_amneg_mapped %>%
    left_join(pathway_genes, by = c("id" = "gs_name")) %>%
    rowwise() %>%
    mutate(
      overlapping_genes = list(intersect(genes_in_pathway, input_genes)),
      geneID = paste(overlapping_genes, collapse = "/"),
      Count = length(overlapping_genes),
      # Compute GeneRatio and BgRatio with correct denominators
      GeneRatio = paste(Count, total_input_genes_in_kegg, sep = "/"),
      BgRatio = paste(length(genes_in_pathway[[1]]), total_background_genes_in_kegg, sep = "/"),
      pvalue = p_value_hyper,
      p.adjust = p_adjust_hyper,
      qvalue = p_adjust_hyper  # Assuming adjusted p-values and q-values are the same
    ) %>%
    ungroup()
  
  # if enrichment_amneg_mapped$Description is NA, replace with $id value, but replace "KEGG_" with "" and set to all lower case except the first letter
  enrichment_amneg_mapped$Description[is.na(enrichment_amneg_mapped$Description)] <- enrichment_amneg_mapped$id[is.na(enrichment_amneg_mapped$Description)]
  enrichment_amneg_mapped$Description <- tolower(sub("KEGG_", "", enrichment_amneg_mapped$Description))
  
  # Step 5: Prepare the final data frame
  final_df <- enrichment_amneg_mapped %>%
    mutate(
      category = NA,      # KEGG categories can be added if available
      subcategory = NA    # KEGG subcategories can be added if available
    ) %>%
    dplyr::select(
      category,
      subcategory,
      ID = KEGG_ID,
      Description,
      GeneRatio,
      BgRatio,
      pvalue,
      p.adjust,
      qvalue,
      geneID,
      Count
    )
  
  # Filter significant terms if needed (e.g., p.adjust < 0.2)
  final_df <- final_df %>%
    filter(p.adjust < 0.05)  # Adjust cutoff as per your requirement
  
  
  
  # Display the final data frame
  print("Final enrichment results:")
  print(final_df)
  
  
  # Load the DOSE package for enrichResult class
  library(DOSE)
  
  # Prepare geneSets for the enrichResult object
  geneSets <- msigdb_kegg %>%
    group_by(gs_id) %>%
    summarize(genes = list(unique(human_gene_symbol))) %>%
    deframe()
  
  # Create the enrichResult object
  enrich_res <- new(
    "enrichResult",
    result = final_df,
    pvalueCutoff = 0.2,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    gene = input_genes,
    universe = background_genes,
    geneSets = geneSets,
    organism = "hsa",
    keytype = "SYMBOL",
    ontology = "KEGG",
    readable = FALSE
  )
  
  # Display the enrichResult object
  print(enrich_res)
  dotplot(enrich_res, showCategory=40)
}

## Verision two: visualize rGREAT results as network plot: ----
if(TRUE){
  # Load necessary libraries
  library(tidyverse)
  library(tidygraph)
  library(ggraph)
  # print version number
  packageVersion("ggraph")
  library(GenomicRanges)
  library(org.Hs.eg.db)  # For gene annotation
  library(scatterpie)     # For pie chart nodes
  library(ggrepel)
  
  # Set the pie_plot_setting variable
  pie_plot_setting <- FALSE  # Set to TRUE for pie charts, FALSE for dark grey pathways
  
  selected_plot <- "all" # or pos or neg or all
  if(selected_plot == "pos"){
    selected_res <- res_ampos
  } else if(selected_plot == "neg"){
    selected_res <- res_amneg
  } else {
    selected_res <- res_amall
  }
  
  # Extract the enrichment table
  great_table <- getEnrichmentTable(selected_res)
  
  # Sort by p_adjust_hyper from smallest to largest
  # great_table <- great_table[order(great_table$p_adjust_hyper),]
  
  # Filter for significant pathways (adjust p-value threshold as needed)
  significant_pathways_og <- great_table %>%
    filter(p_adjust_hyper < 0.05)
  significant_pathways_og
  
  # sort by p_adjust_hyper
  significant_pathways_og <- significant_pathways_og[order(significant_pathways_og$fold_enrichment_hyper),]
  significant_pathways_og <- significant_pathways_og[1:15,]
  
  # kepp only rows with fold_enrichment_hyper > 2
  #significant_pathways_og <- significant_pathways_og[significant_pathways_og$fold_enrichment_hyper > 2,]
  
  # Keep only top 20 with smallest p_adjust_hyper
  #significant_pathways <- significant_pathways_og[1:20,]
  significant_pathways <- significant_pathways_og
  
  # Extract the pathway IDs
  pathway_ids <- significant_pathways$id
  
  # Extract the gene sets for significant pathways
  gene_sets <- selected_res@gene_sets[names(selected_res@gene_sets) %in% pathway_ids]
  
  # Create a GRanges object from res_table_tb
  deseq2_gr <- GRanges(
    seqnames = res_table_tb$chr,
    ranges = IRanges(start = res_table_tb$start, end = res_table_tb$end),
    strand = "*",
    gene = res_table_tb$gene,
    baseMean = res_table_tb$baseMean,
    log2FoldChange = res_table_tb$log2FoldChange,
    lfcSE = res_table_tb$lfcSE,
    stat = res_table_tb$stat,
    pvalue = res_table_tb$pvalue,
    padj = res_table_tb$padj
  )
  
  # Filter to keep only regions with p-value < 0.005
  #deseq2_gr <- deseq2_gr[deseq2_gr$pvalue < 0.01]
  
  # Filter to keep only regions with log2FoldChange > 0.13
  #deseq2_gr <- deseq2_gr[deseq2_gr$log2FoldChange > 0.2]
  
  # if selected_plot == "all" sort by abs, if "pos" sort by stat in decreasing, or "neg" sort by stat in increasing
  if(selected_plot == "all"){
    deseq2_gr <- deseq2_gr[order(abs(deseq2_gr$stat), decreasing = TRUE)]
  } else if(selected_plot == "pos"){
    deseq2_gr <- deseq2_gr[order(deseq2_gr$stat, decreasing = TRUE)]
  } else {
    deseq2_gr <- deseq2_gr[order(deseq2_gr$stat)]
  }
  
  # filter to keep only top 200 rows with highest stat value
  deseq2_gr <- deseq2_gr[1:200]
  
  # Extract extended TSS regions with gene IDs from rGREAT output
  extended_tss <- res_ampos@extended_tss  # Contains gene_id in metadata
  
  # Perform overlaps between DESeq2 regions and extended TSS regions
  overlaps <- findOverlaps(deseq2_gr, extended_tss)
  
  # Create a data frame mapping regions to gene IDs
  region_gene_df <- data.frame(
    region_index = queryHits(overlaps),
    gene_id = extended_tss$gene_id[subjectHits(overlaps)]
  )
  
  # Add region information and fold changes to region_gene_df
  region_gene_df <- region_gene_df %>%
    mutate(
      region_name = deseq2_gr$gene[region_index],
      log2FoldChange = deseq2_gr$log2FoldChange[region_index],
      pvalue = deseq2_gr$pvalue[region_index],     # Include p-value
      padj = deseq2_gr$padj[region_index]
    )
  
  # Map Entrez Gene IDs to Gene Symbols
  gene_ids <- unique(region_gene_df$gene_id)
  
  # Use mapIds to get gene symbols
  gene_symbols <- mapIds(org.Hs.eg.db, keys = as.character(gene_ids),
                         column = "SYMBOL", keytype = "ENTREZID",
                         multiVals = "first")
  
  # Create a data frame with gene_id and gene_symbol
  gene_id_symbol_df <- data.frame(
    gene_id = gene_ids,
    gene_symbol = gene_symbols,
    stringsAsFactors = FALSE
  )
  
  # Merge gene symbols into region_gene_df
  region_gene_df <- left_join(region_gene_df, gene_id_symbol_df, by = "gene_id")
  
  # Add a count of unique region_name per gene_symbol
  region_gene_df <- region_gene_df %>%
    group_by(gene_symbol) %>%
    mutate(
      unique_regions = n_distinct(region_name)
    ) %>%
    ungroup() # Ungroup to avoid issues in summarizing
  
  # save dataframe as csv
  write.csv(region_gene_df, paste("./Output/DESeq2/GSEA_regions/region_gene_df_",contrast_groups[2],contrast_groups[3],"_",ver,".csv", sep = ""), row.names = FALSE)
  
  
  # Now use reframe with conditional formatting
  region_gene_df <- region_gene_df %>%
    group_by(gene_symbol) %>%
    reframe(
      region_index = dplyr::first(region_index), # Retain the first region_index or customize as needed
      gene_id = dplyr::first(gene_id),           # Retain the first gene_id or customize as needed
      region_name = paste(unique(region_name), collapse = ", "), # Concatenate unique region_names
      log2FoldChange = mean(log2FoldChange), # Average log2FoldChange or customize as needed
      pvalue = min(pvalue),               # Take minimum pvalue or customize as needed
      padj = min(padj),                   # Take minimum adjusted p-value or customize as needed
      gene_symbol = if (unique(unique_regions) > 1) {
        paste(gene_symbol, "(", unique(unique_regions), ")", sep = "")
      } else {
        gene_symbol
      }
    )
  
  # convert log2FoldChange to PercentChange
  region_gene_df$FoldChange <- (2^region_gene_df$log2FoldChange - 1) * 100
  
  
  # save dataframe as csv
  write.csv(region_gene_df, paste("./Output/DESeq2/GSEA_regions/region_gene_df_collapsed",contrast_groups[2],contrast_groups[3],"_",ver,".csv", sep = ""), row.names = FALSE)
  
  # drop duplicate rows
  region_gene_df <- unique(region_gene_df)
  
  # Create a data frame mapping gene IDs to pathways
  gene_pathway_df <- enframe(gene_sets, name = "pathway", value = "gene_ids") %>%
    unnest(cols = c(gene_ids)) %>%
    rename(gene_id = gene_ids)
  
  # Update pathway names: remove "KEGG_" and convert to lowercase
  gene_pathway_df$pathway <- tolower(sub("KEGG_", "", gene_pathway_df$pathway))
  
  # Merge to get regions connected to pathways via genes
  region_gene_pathway_df <- inner_join(region_gene_df, gene_pathway_df, by = "gene_id")
  
  
  
  # Compute counts of positive and negative FoldChange regions for each pathway
  pathway_counts <- region_gene_pathway_df %>%
    group_by(pathway) %>%
    summarise(
      num_regions = n(),
      num_positive = sum(FoldChange > 0),
      num_negative = sum(FoldChange < 0)
    )
  
  # Prepare edges data frame (same as before)
  edges <- region_gene_pathway_df %>%
    dplyr::select(region_name, pathway) %>%
    distinct()
  
  # Nodes for regions (include FoldChange, pvalue, and gene_symbol)
  region_nodes <- region_gene_pathway_df %>%
    dplyr::select(name = region_name, FoldChange, pvalue, gene_symbol) %>%
    distinct() %>%
    mutate(type = "Region")
  
  # Nodes for pathways
  pathway_nodes <- edges %>%
    dplyr::select(name = pathway) %>%
    distinct() %>%
    mutate(type = "Pathway", FoldChange = NA, pvalue = NA, gene_symbol = NA)
  
  # Merge counts into pathway_nodes
  pathway_nodes <- pathway_nodes %>%
    left_join(pathway_counts, by = c("name" = "pathway"))
  
  # Combine nodes
  nodes <- bind_rows(region_nodes, pathway_nodes) %>%
    mutate(id = row_number())
  
  # Create a lookup table for node IDs
  node_lookup <- nodes %>% dplyr::select(name, id)
  
  # Map node IDs in edges
  edges <- edges %>%
    left_join(node_lookup, by = c("region_name" = "name")) %>%
    rename(from = id) %>%
    left_join(node_lookup, by = c("pathway" = "name")) %>%
    rename(to = id)
  
  # Ensure edges have necessary columns
  edges <- edges %>%
    dplyr::select(from, to)
  
  # Create a tidygraph object
  graph <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE)
  
  # Compute node degree and filter nodes with at least 1 edge
  filtered_graph <- graph %>%
    activate(nodes) %>%
    mutate(node_degree = centrality_degree()) %>%  # Add degree as a node attribute
    filter(node_degree >= 1)                      # Filter nodes with at least 1 edge
  
  # Compute size for region nodes based on -log10(pvalue)
  filtered_graph <- filtered_graph %>%
    mutate(size = ifelse(type == "Region", -log10(pvalue + 1e-16), NA))
  
  # Optionally rescale sizes to a reasonable range (e.g., 2 to 5)
  filtered_graph <- filtered_graph %>%
    mutate(size = ifelse(type == "Region",
                         scales::rescale(size, to = c(2, 4)),
                         5))  # Set pathway node sizes to 5
  
  # Compute the layout and extract node positions
  layout <- create_layout(filtered_graph, layout = "kk") #fr
  
  # Extract node data including x and y positions
  node_data <- as_tibble(layout)
  
  # Prepare data for pathway nodes (pie charts)
  pathway_pie_data <- node_data %>%
    filter(type == "Pathway") %>%
    dplyr::select(name, x, y, num_positive, num_negative)
  
  # Handle any NA values (e.g., if num_positive or num_negative is NA)
  pathway_pie_data <- pathway_pie_data %>%
    mutate(
      num_positive = ifelse(is.na(num_positive), 0, num_positive),
      num_negative = ifelse(is.na(num_negative), 0, num_negative)
    )
  
  
  
  # Separate region and pathway data
  region_node_data <- node_data %>%
    filter(type == "Region")
  
  pathway_pie_data <- node_data %>%
    filter(type == "Pathway")
  
  # Base plot with edges and region nodes
  p_net_pos <- ggraph(layout) +
    # Draw edges
    geom_edge_link(color = "grey", alpha = 1) +
    # Draw region nodes, color by FoldChange and size by -log10(pvalue)
    geom_node_point(data = region_node_data, aes(x = x, y = y, color = FoldChange, size = size))
  
  # Conditionally add pathway nodes as pie charts or dark grey points
  if (pie_plot_setting) {
    p_net_pos <- p_net_pos +
      # Draw pathway nodes as pie charts
      geom_scatterpie(
        data = pathway_pie_data,
        aes(x = x, y = y),
        cols = c("num_positive", "num_negative"),
        pie_scale = 1.5,
        color = NA
      ) +
      # Customize fill scale for pie chart segments
      scale_fill_manual(
        values = c("num_positive" = "#B71C1C", "num_negative" = "#0D47A1"),
        labels = c("Positive", "Negative"),
        name = "Region Direction"
      )
  } else {
    p_net_pos <- p_net_pos +
      # Draw pathway nodes as dark grey points
      geom_node_point(
        data = pathway_pie_data,
        aes(x = x, y = y),
        color = "#7f7f7f",
        alpha=1,
        size = 5  # Adjust size as needed
      )
  }
  
  # Add labels to pathway nodes using ggrepel
  p_net_pos <- p_net_pos +
    geom_text_repel(
      data = pathway_pie_data,
      aes(x = x, y = y, label = name),
      size = 3,
      fontface = "bold",
      box.padding = 0.5,
      point.padding = 0.3,
      force = 2,
      segment.color = "grey50",
      max.overlaps = Inf
    ) +
    # Add labels to region nodes (gene symbols) using ggrepel
    geom_text_repel(
      data = region_node_data,
      aes(x = x, y = y, label = gene_symbol),
      size = 3,
      box.padding = 0.5,
      point.padding = 0.3,
      force = 2,
      segment.color = "grey50",
      color = "darkgrey",  # Added this line to make region labels dark grey
      max.overlaps = Inf
    )
  
  # Check if all FoldChange values are positive, negative, or mixed
  if (all(region_node_data$FoldChange > 0)) {
    # All positive: light red to red gradient
    p_net_pos <- p_net_pos + scale_color_gradient(low = "#FFEBEE", high = "#B71C1C", name = "Percent Change")
  } else if (all(region_node_data$FoldChange < 0)) {
    # All negative: blue to light blue gradient
    p_net_pos <- p_net_pos + scale_color_gradient(low = "#0D47A1", high = "#CCDFFF", name = "Percent Change")
  } else {
    # Mixed: standard gradient with white midpoint
    p_net_pos <- p_net_pos + scale_color_gradient2(
      low = "#0D47A1", 
      mid = "white", 
      high = "#B71C1C", 
      midpoint = 0,
      name = "Log2 Fold Change"
    )
  }
  
  # Add theme and labels
  p_net_pos <- p_net_pos +
    # Theme adjustments
    theme_minimal() +
    labs(title = "Network of KEGG Pathway Enrichments Driven by Genomic Regions") +
    theme(legend.position = "bottom")
  
  # If pie charts are enabled, ensure the legend includes the pie chart information
  if (pie_plot_setting) {
    p_net_pos <- p_net_pos + guides(fill = guide_legend(override.aes = list(alpha = 1)))
  }
  
  # drop grid lines
  p_net_pos <- p_net_pos + theme(panel.grid = element_blank())
  
  # drop x and y axis labels and titles
  p_net_pos <- p_net_pos + theme(
    axis.text.x = element_blank(), 
    axis.text.y = element_blank(),
    axis.title.x = element_blank(), 
    axis.title.y = element_blank()
  )
  
  # change "size" label to "p-value"
  p_net_pos <- p_net_pos + labs(size = "p-value")
  
  # move legends to be vertical on right side
  p_net_pos <- p_net_pos + theme(legend.position = "right")

  
  # Print the plot
  p_net_pos
  
  # creat folder if it doesn't exist
  dir.create("./Output/DESeq2/GSEA_regions", showWarnings = FALSE)
  
  # Save as an editable SVG using ggsave
  ggsave(
    filename = paste("./Output/DESeq2/GSEA_regions/network_pos_", contrast_groups[2], contrast_groups[3], "_", ver, "PMnegPMpos.svg", sep = ""),
    plot = p_net_pos,
    width = 10, height = 7,
    device = "svg"
  )
  
  # Save as PNG for additional compatibility
  ggsave(
    filename = paste("./Output/DESeq2/GSEA_regions/network_pos_", contrast_groups[2], contrast_groups[3], "_", ver, "PMnegPMpos.png", sep = ""),
    plot = p_net_pos,
    width = 10, height = 7,
    dpi = 300,
    units = "in"
  )
  
  # save significant_pathways_og as csv
  write.csv(significant_pathways_og, paste("./Output/DESeq2/GSEA_regions/significant_pathways_",contrast_groups[2],contrast_groups[3],"_",ver,".csv", sep = ""), row.names = FALSE)

  
  # Load the ggplot2 library
  library(ggplot2)
  
  # Filter the data based on your adjusted p-value threshold
  significant_pathways_og <- great_table %>%
    filter(p_adjust_hyper < 0.05)  # Adjust the threshold based on your analysis
  
  # Replace "KEGG_" in the 'id' column with "" and convert to lowercase
  significant_pathways_og$id <- tolower(gsub("KEGG_", "", significant_pathways_og$id))
  
  # Replace "_" with " " in the 'id' column
  significant_pathways_og$id <- gsub("_", " ", significant_pathways_og$id)
  
  # Sort the data frame by 'fold_enrichment_hyper' in decreasing order
  significant_pathways_og <- significant_pathways_og[order(significant_pathways_og$fold_enrichment_hyper, decreasing = TRUE),]
  
  # Convert 'id' into a factor with levels in the order of 'id' in the sorted data frame
  significant_pathways_og$id <- factor(significant_pathways_og$id, levels = significant_pathways_og$id)
  
  significant_pathways_og_plot <- ggplot(significant_pathways_og, aes(x = fold_enrichment_hyper, y = reorder(id, fold_enrichment_hyper))) +
    geom_point(aes(color = p_adjust_hyper, size = genome_fraction)) +
    scale_color_gradient(low = "#2e2e2e", high = "#CCCCCC") +
    labs(x = "Hypergeometric Fold Enrichment",
         y = "Pathway",
         color = "Adjusted p-value",
         size = "Genome Fraction") +
    scale_x_continuous(expand = expansion(mult = c(0.15, 0.15))) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8),
          panel.grid.minor.x = element_blank())
  
  # Save the plot as an editable SVG using ggsave
  ggsave(
    filename = paste("./Output/DESeq2/GSEA_regions/significant_pathways_", contrast_groups[2], contrast_groups[3], "_", ver, "_gsea_kegg_dot.svg", sep = ""),
    plot = significant_pathways_og_plot,
    width = 6, height = 4,
    device = "svg"
  )
  
  # Save the plot as a PNG for compatibility
  ggsave(
    filename = paste("./Output/DESeq2/GSEA_regions/significant_pathways_", contrast_groups[2], contrast_groups[3], "_", ver, "_gsea_kegg_dot.png", sep = ""),
    plot = significant_pathways_og_plot,
    width = 6, height = 4,
    dpi = 300,
    units = "in"
  )
  
  if(FALSE){
    # Load necessary libraries
    library(tidyverse)
    library(vegan)
    library(ggplot2)
    library(msigdbr)
    library(ggdendro)
    
    # Step 1: Filter significant pathways from your enrichment analysis
    # Assuming 'enrichment_amall' is your enrichment analysis result from GREAT
    significant_pathways_og <- enrichment_amall %>%
      filter(p_adjust_hyper < 0.05)  # Adjust the threshold based on your analysis
    
    # Step 2: Prepare pathway names for matching
    significant_pathways_og$gs_name <- significant_pathways_og$id
    
    # Step 3: Use msigdbr to get KEGG pathways and their gene symbols
    msigdbr_df <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
    
    # Convert pathway names to uppercase for consistent matching
    msigdbr_df$gs_name_upper <- toupper(msigdbr_df$gs_name)
    
    # Step 4: Filter msigdbr to include only your significant pathways
    filtered_msigdb <- msigdbr_df %>%
      filter(gs_name_upper %in% significant_pathways_og$gs_name)
    
    # Step 5: Create a list of genes for each pathway
    pathway_genes_list <- filtered_msigdb %>%
      group_by(gs_name_upper) %>%
      summarize(genes = list(unique(gene_symbol))) %>%
      deframe()
    
    # Verify that all pathways have been matched
    if(length(pathway_genes_list) < nrow(significant_pathways_og)) {
      warning("Some pathways were not matched with gene sets.")
    }
    # Step 6: Create a binary matrix of pathways vs genes
    all_genes <- unique(unlist(pathway_genes_list))
    
    binary_matrix <- matrix(0, nrow=length(pathway_genes_list), ncol=length(all_genes))
    rownames(binary_matrix) <- names(pathway_genes_list)
    colnames(binary_matrix) <- all_genes
    
    for(i in seq_along(pathway_genes_list)){
      pathway_name <- names(pathway_genes_list)[i]
      genes_in_pathway <- pathway_genes_list[[i]]
      # Ensure genes are in the column names
      genes_in_pathway <- genes_in_pathway[genes_in_pathway %in% all_genes]
      binary_matrix[pathway_name, genes_in_pathway] <- 1
    }
    
    # Step 7: Compute the Jaccard distance matrix
    dist_matrix <- vegdist(binary_matrix, method="jaccard")
    
    # Step 8: Perform hierarchical clustering
    hc <- hclust(dist_matrix, method="average")
    
    # **NEW**: Set a distance threshold to define clusters
    distance_threshold <- 0.80  # Adjust this value as needed
    
    # Cut the dendrogram at the specified distance threshold to get clusters
    clusters <- cutree(hc, h = distance_threshold)
    
    # Convert hclust object to a dendrogram object
    dend <- as.dendrogram(hc)
    
    # Step 10: Modify labels
    # Replace "KEGG_" with "" and "_" with " " in the labels
    old_labels <- labels(dend)
    new_labels <- old_labels %>%
      str_replace_all("KEGG_", "") %>%
      str_replace_all("_", " ")
    labels(dend) <- new_labels
    
    # Use dendextend to color branches according to clusters
    dend_colored <- color_branches(dend, h = distance_threshold)
    
    par(mar = c(5, 1, 1, 30) + 0.1)  # Adjust margins to fit labels
    # Plot the dendrogram with colored clusters
    plot(dend_colored, 
         main = paste("Dendrogram of Pathways (Threshold =", distance_threshold, ")"), 
         xlab = "Distance",
         horiz = TRUE,  # This parameter makes the plot vertical
         cex.lab = 0.8, cex.axis = 0.8, cex.main = 1,
         las = 1)  # Makes y-axis labels horizontal
    
    # Add light grey boxes around clusters meeting the cutoff
    rect.dendrogram(dend_colored, 
                    h = distance_threshold, 
                    border = "lightgrey", 
                    lty = 1, 
                    lwd = 2, 
                    horiz = TRUE)
    
    
    
  }
}

# print ggplot version
packageVersion("ggplot2")

### SAVE RESULTS FILES----
if(TRUE){
  # Create named list for filenames
  # Load required libraries
  library(data.table)
  library(dplyr)
  library(tibble)
  library(GenomicRanges)
  library(IRanges)
  library(annotatr)
  library(tidyr)
  
  # Create named list for filenames
  create_filename_list <- function(contrast_groups, ver) {
    list(
      res_table = paste0("./Output/DESeq2/Results/all_results_", contrast_groups[2], contrast_groups[3], "_", ver, ".txt"),
      res_table_annotated = paste0("./Output/DESeq2/Results/all_results_annotated_", contrast_groups[2], contrast_groups[3], "_", ver, ".txt"),
      res_table_annotated_enhancer = paste0("./Output/DESeq2/Results/all_results_annotated_enhancer_", contrast_groups[2], contrast_groups[3], "_", ver, ".txt"),
      sig = paste0("./Output/DESeq2/Results/significant_results_", contrast_groups[2], contrast_groups[3], "_", ver, ".txt"),
      sig_annotated = paste0("./Output/DESeq2/Results/significant_results_annotated_", contrast_groups[2], contrast_groups[3], "_", ver, ".txt"),
      sig_annotated_enhancer = paste0("./Output/DESeq2/Results/significant_results_annotated_enhancer_", contrast_groups[2], contrast_groups[3], "_", ver, ".txt"),
      GSEA = paste0("./Output/DESeq2/Results/GSEA_", contrast_groups[2], contrast_groups[3], "_", ver, ".rnk"),
      GSEA_annotated = paste0("./Output/DESeq2/Results/GSEA_annotated_", contrast_groups[2], contrast_groups[3], "_", ver, ".rnk"),
      GSEA_annotated_enhancer = paste0("./Output/DESeq2/Results/GSEA_annotated_enhancer_", contrast_groups[2], contrast_groups[3], "_", ver, ".rnk"),
      ORA_all = paste0("./Output/DESeq2/Results/ORA_all_", contrast_groups[2], contrast_groups[3], "_", ver, ".txt"),
      ORA_plus = paste0("./Output/DESeq2/Results/ORA_plus_", contrast_groups[2], contrast_groups[3], "_", ver, ".txt"),
      ORA_minus = paste0("./Output/DESeq2/Results/ORA_minus_", contrast_groups[2], contrast_groups[3], "_", ver, ".txt"),
      ORA_all_annotated = paste0("./Output/DESeq2/Results/ORA_all_annotated_", contrast_groups[2], contrast_groups[3], "_", ver, ".txt"),
      ORA_plus_annotated = paste0("./Output/DESeq2/Results/ORA_plus_annotated_", contrast_groups[2], contrast_groups[3], "_", ver, ".txt"),
      ORA_minus_annotated = paste0("./Output/DESeq2/Results/ORA_minus_annotated_", contrast_groups[2], contrast_groups[3], "_", ver, ".txt"),
      ORA_all_annotated_enhancer = paste0("./Output/DESeq2/Results/ORA_all_annotated_enhancer_", contrast_groups[2], contrast_groups[3], "_", ver, ".txt"),
      ORA_plus_annotated_enhancer = paste0("./Output/DESeq2/Results/ORA_plus_annotated_enhancer_", contrast_groups[2], contrast_groups[3], "_", ver, ".txt"),
      ORA_minus_annotated_enhancer = paste0("./Output/DESeq2/Results/ORA_minus_annotated_enhancer_", contrast_groups[2], contrast_groups[3], "_", ver, ".txt")
    )
  }
  
  # Main function to process and save results
  process_and_save_results <- function(output_results_tables, res_table, sig, res_table_tb, contrast_groups, ver, region_type, genebodies = FALSE) {
    fn_list <- create_filename_list(contrast_groups, ver)
    
    if (output_results_tables == 1) {
      cat("Processing and saving main results...\n")
      
      # Define the files to process based on genebodies parameter
      files_to_process <- if (genebodies) {
        c("res_table", "sig", "GSEA", "ORA_all", "ORA_plus", "ORA_minus")
      } else {
        names(fn_list)
      }
      
      for (file_type in files_to_process) {
        if (!file.exists(fn_list[[file_type]])) {
          switch(file_type,
                 "res_table" = {
                   write.table(res_table, file = fn_list$res_table, sep="\t", quote=F, col.names=NA)
                   cat("Saved all results to:", fn_list$res_table, "\n")
                 },
                 "sig" = {
                   write.table(sig, file = fn_list$sig, sep="\t", quote=F, col.names=NA)
                   cat("Saved significant results to:", fn_list$sig, "\n")
                 },
                 "GSEA" = {
                   gsea_data <- res_table_tb[, c("gene", "log2FoldChange")]
                   fwrite(gsea_data, file = fn_list$GSEA, sep = "\t", quote = FALSE, col.names = FALSE)
                   cat("Saved GSEA results to:", fn_list$GSEA, "\n")
                 },
                 "ORA_all" = {
                   ora_all <- sig$gene
                   fwrite(list(gene = ora_all), file = fn_list$ORA_all, sep = "\t", quote = FALSE, col.names = FALSE)
                   cat("Saved ORA all results to:", fn_list$ORA_all, "\n")
                 },
                 "ORA_plus" = {
                   ora_plus <- sig$gene[sig$log2FoldChange > 0]
                   fwrite(list(gene = ora_plus), file = fn_list$ORA_plus, sep = "\t", quote = FALSE, col.names = FALSE)
                   cat("Saved ORA plus results to:", fn_list$ORA_plus, "\n")
                 },
                 "ORA_minus" = {
                   ora_minus <- sig$gene[sig$log2FoldChange < 0]
                   fwrite(list(gene = ora_minus), file = fn_list$ORA_minus, sep = "\t", quote = FALSE, col.names = FALSE)
                   cat("Saved ORA minus results to:", fn_list$ORA_minus, "\n")
                 }
          )
        } else {
          cat("File already exists:", fn_list[[file_type]], "\n")
        }
      }
      
      if (region_type == 0 && !genebodies) {
        process_region_type_zero(fn_list, contrast_groups, ver, res_table, sig, res_table_tb)
      }
    }
    
    return(fn_list)
  }
  
  process_region_type_zero <- function(fn_list, contrast_groups, ver, res_table, sig, res_table_tb) {
    cat("Processing region type zero...\n")
    
    # Load enhancer lookup file
    enhancer_lookup_df <- fread("~/reference_genomes/GeneHancer_AnnotSV_gene_association_scores_v5.20.txt")
    enhancer_lookup_df <- enhancer_lookup_df[is_elite == 1]
    
    # Create a lookup table for GH symbols
    gh_lookup <- enhancer_lookup_df[startsWith(GHid, "GH") & nchar(GHid) >= 8, .(GHid, symbol)]
    setkey(gh_lookup, GHid)
    
    # Source the file containing annotate_features2 function
    source("~/5hmC-Pathway-Analysis/R Code/repeated_glmnet_funcs_v1.R")
    
    # Process res_table
    process_file("res_table", res_table, fn_list, contrast_groups, ver, gh_lookup)
    
    # Process sig
    process_file("sig", sig, fn_list, contrast_groups, ver, gh_lookup)
    
    # Process GSEA
    process_file("GSEA", res_table_tb[, c("gene", "log2FoldChange")], fn_list, contrast_groups, ver, gh_lookup)
    
    # Process ORA
    process_file("ORA_all", data.frame(feature = sig$gene), fn_list, contrast_groups, ver, gh_lookup)
    process_file("ORA_plus", data.frame(feature = sig$gene[sig$log2FoldChange > 0]), fn_list, contrast_groups, ver, gh_lookup)
    process_file("ORA_minus", data.frame(feature = sig$gene[sig$log2FoldChange < 0]), fn_list, contrast_groups, ver, gh_lookup)
  }
  
  process_file <- function(file_type, data, fn_list, contrast_groups, ver, gh_lookup) {
    cat("Processing file:", file_type, "\n")
    
    # Check if annotated and enhancer annotated files already exist
    if (file.exists(fn_list[[paste0(file_type, "_annotated")]]) && 
        file.exists(fn_list[[paste0(file_type, "_annotated_enhancer")]])) {
      cat("Annotated files already exist for:", file_type, "\n")
      return()
    }
    
    # Convert data to data.frame if it's not already
    if (!is.data.frame(data)) {
      data <- as.data.frame(data)
    }
    
    # Handle different data structures
    if (file_type == "res_table") {
      # For res_table, we need to add a feature column
      data$feature <- rownames(data)
      data <- data[, c("feature", setdiff(names(data), "feature"))]
    } else if (file_type == "sig") {
      # For sig, rename the "gene" column to "feature"
      if ("gene" %in% colnames(data)) {
        colnames(data)[colnames(data) == "gene"] <- "feature"
      } else {
        cat("Warning: 'gene' column not found in sig data. Current column names:", colnames(data), "\n")
      }
    } else if (file_type %in% c("GSEA", "ORA_all", "ORA_plus", "ORA_minus")) {
      # For GSEA and ORA files, ensure the first column is named 'feature'
      colnames(data)[1] <- "feature"
    }
    
    # Ensure the first column is named 'feature'
    if (ncol(data) > 0 && colnames(data)[1] != "feature") {
      colnames(data)[1] <- "feature"
    } else if (ncol(data) == 0) {
      cat("Error: Data is empty for file type:", file_type, "\n")
      return()
    }
    
    # Remove rows with NA in feature column
    na_count <- sum(is.na(data$feature))
    if (na_count > 0) {
      cat("Found", na_count, "rows with NA values in the feature column in", file_type, "\n")
      data <- data[!is.na(data$feature), ]
    }
    
    # Apply the annotate_features2 function
    tryCatch({
      annotated_df <- annotate_features2(data)
      
      # Save annotated version for all file types
      if (!file.exists(fn_list[[paste0(file_type, "_annotated")]])) {
        if (file_type == "GSEA") {
          fwrite(annotated_df[, c("annot.symbol", "log2FoldChange")], 
                 file = fn_list[[paste0(file_type, "_annotated")]], 
                 sep = "\t", quote = FALSE, col.names = FALSE)
        } else if (file_type %in% c("ORA_all", "ORA_plus", "ORA_minus")) {
          annotated_symbols <- data.frame(annot.symbol = annotated_df[, "annot.symbol"])
          fwrite(annotated_symbols, 
                 file = fn_list[[paste0(file_type, "_annotated")]], 
                 sep = "\t", quote = FALSE, col.names = FALSE)
        } else {
          fwrite(annotated_df, file = fn_list[[paste0(file_type, "_annotated")]], sep = "\t", quote = FALSE)
        }
        cat("Saved annotated file:", fn_list[[paste0(file_type, "_annotated")]], "\n")
      }
      
      # Process enhanced file
      process_enhanced_file(file_type, annotated_df, fn_list, gh_lookup)
    }, error = function(e) {
      cat("Error in annotate_features2 for", file_type, ":", conditionMessage(e), "\n")
      print(str(data))  # Print structure of data for debugging
    })
  }
  
  process_enhanced_file <- function(file_type, annotated_df, fn_list, gh_lookup) {
    cat("Processing enhanced file:", file_type, "\n")
    
    setDT(annotated_df)
    
    annotated_df_enh <- annotated_df[, .(
      symbols_gh_full = list(process_symbols(all_symbols, gh_lookup))
    ), by = names(annotated_df)]
    
    annotated_df_enh <- annotated_df_enh[, .(
      feature = rep(feature, sapply(symbols_gh_full, length)),
      symbols_gh_full = unlist(symbols_gh_full)
    ), by = setdiff(names(annotated_df_enh), c("feature", "symbols_gh_full"))]
    
    if (!file.exists(fn_list[[paste0(file_type, "_annotated_enhancer")]])) {
      if (file_type %in% c("GSEA", "ORA_all", "ORA_plus", "ORA_minus")) {
        # For GSEA and ORA files, only write specific columns
        if (file_type == "GSEA") {
          fwrite(annotated_df_enh[, .(symbols_gh_full, log2FoldChange)], 
                 file = fn_list[[paste0(file_type, "_annotated_enhancer")]], 
                 sep = "\t", quote = FALSE, col.names = FALSE)
        } else {
          fwrite(annotated_df_enh[, .(symbols_gh_full)], 
                 file = fn_list[[paste0(file_type, "_annotated_enhancer")]], 
                 sep = "\t", quote = FALSE, col.names = FALSE)
        }
      } else {
        fwrite(annotated_df_enh, file = fn_list[[paste0(file_type, "_annotated_enhancer")]], sep = "\t", quote = FALSE)
      }
      cat("Saved enhanced annotated file:", fn_list[[paste0(file_type, "_annotated_enhancer")]], "\n")
    } else {
      cat("Enhanced annotated file already exists:", fn_list[[paste0(file_type, "_annotated_enhancer")]], "\n")
    }
  }
  
  # Function to process symbols efficiently
  process_symbols <- function(symbols, gh_lookup) {
    symbols_list <- strsplit(symbols, ", ")[[1]]
    gh_symbols <- symbols_list[startsWith(symbols_list, "GH") & nchar(symbols_list) >= 8]
    non_gh_symbols <- symbols_list[!(symbols_list %in% gh_symbols)]
    
    if (length(gh_symbols) > 0) {
      matched_symbols <- gh_lookup[gh_symbols, on = "GHid"]$symbol
      return(c(matched_symbols, non_gh_symbols))
    } else {
      return(symbols_list)
    }
  }
  
  # Main execution
  fn_list <- process_and_save_results(output_results_tables, 
                                      res_table, sig, 
                                      res_table_tb, 
                                      contrast_groups, ver, 
                                      region_type, genebodies = genebodies_type)
  
  # New function to generate BED files
  generate_bed_files <- function(sig, contrast_groups, ver) {
    # Separate into up-regulated and down-regulated
    sig_up <- sig[sig$log2FoldChange > 0, ]
    sig_down <- sig[sig$log2FoldChange < 0, ]
    
    # Function to create bed entry from gene
    create_bed_entry <- function(row) {
      gene_name <- as.character(row["gene"])  
      split_gene_name <- strsplit(gene_name, "_")[[1]] #Split into columns with "_"
      
      chr <- split_gene_name[1]
      start <- as.numeric(split_gene_name[2])
      end <- as.numeric(split_gene_name[3])
      
      data.frame(
        chr = chr,
        start = start - 1, # subtract 1 to match the provided format (0-based)
        end = end,
        score = end-start, 
        strand = "+"
      )
    }
    
    # Create bed files if they don't exist, one per strand
    for (direction in c("positive", "negative")) {
      for (strand in c("+", "-")) {
        filename <- paste0("./Output/DESeq2/Results/", contrast_groups[2], "_", direction, "_", strand,"_",ver, ".bed")
        if (!file.exists(filename)) {
          if (direction == "positive") {
            bed_data <- sig_up %>% filter(grepl(strand, gene)) %>% rowwise() %>% do(create_bed_entry(.))
          } else {
            bed_data <- sig_down %>% filter(grepl(strand, gene)) %>% rowwise() %>% do(create_bed_entry(.))
          }
          
          # Filter empty dataframes
          if (nrow(bed_data) > 0) {
            fwrite(bed_data, file = filename, sep = "\t", col.names = FALSE)
            cat("Saved BED file:", filename, "\n")
          }
        } else {
          cat("BED file already exists:", filename, "\n")
        }
      }
    }
  }
  
  # Generate BED files
  generate_bed_files(sig, contrast_groups, ver)
}


### GSEA and ORA ANALYSIS prep ----
if(TRUE){
  fn_list <- create_filename_list(contrast_groups, ver)
  # load sig file into sig
  sig <- fread(fn_list$sig, header = TRUE)
  
  if(genebodies_type == TRUE){
    # load res_table into res_table
    res_table <- fread(fn_list$res_table, header = TRUE)
    # Load GSEA and ORA files into a dataframe
    ora_all_enh_annotated <- fread(fn_list$ORA_all, header = FALSE)
    ora_plus_enh_annotated <- fread(fn_list$ORA_plus, header = FALSE)
    ora_neg_enh_annotated <- fread(fn_list$ORA_minus, header = FALSE)
    gsea_annotated_enhancer <- fread(fn_list$GSEA, header = FALSE)
  } else{
    # load res_table into res_table
    res_table <- fread(fn_list$res_table_annotated_enhancer, header = TRUE)
    # Load GSEA_annotated into a dataframe
    ora_all_enh_annotated <- fread(fn_list$ORA_all_annotated_enhancer, header = FALSE)
    ora_plus_enh_annotated <- fread(fn_list$ORA_plus_annotated_enhancer, header = FALSE)
    ora_neg_enh_annotated <- fread(fn_list$ORA_minus_annotated_enhancer, header = FALSE)
    gsea_annotated_enhancer <- fread(fn_list$GSEA_annotated_enhancer, header = FALSE)
  }
  
  
  # drop blank, "NA" and na rows
  ora_all_enh_annotated <- ora_all_enh_annotated[!is.na(V1) & V1 != "NA" & V1 != "",]
  ora_plus_enh_annotated <- ora_plus_enh_annotated[!is.na(V1) & V1 != "NA" & V1 != "",]
  ora_neg_enh_annotated <- ora_neg_enh_annotated[!is.na(V1) & V1 != "NA" & V1 != "",]
  gsea_annotated_enhancer <- gsea_annotated_enhancer[!is.na(V1) & V1 != "NA" & V1 != "",]
  
  ### Keep rows with largest absolute value of V2
  # Ensure gsea_annotated_enhancer is a data.table
  setDT(gsea_annotated_enhancer)
  # Create new column V3 as absolute value of V2
  gsea_annotated_enhancer[, V3 := abs(V2)]
  # Keep row with largest V3 value for each unique V1
  gsea_annotated_enhancer <- gsea_annotated_enhancer[order(-V3), .SD[1], by = V1]
  # Drop V3 column
  gsea_annotated_enhancer[, V3 := NULL]
  
  # Load required libraries
  library(biomaRt)
  library(data.table)
  library(pbapply)  # For progress bars
  
  # Function to convert gene IDs to Ensembl IDs or EntrezIDs
  convert_to_ensembl <- function(genes) {
    # Connect to the Ensembl database
    cat("Connecting to Ensembl database...\n")
    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    
    # Prepare a vector of unique gene IDs
    unique_genes <- unique(genes)
    cat("Number of unique genes:", length(unique_genes), "\n")
    
    # Create a list of attribute queries for different ID types
    queries <- list(
      list(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "ensembl_gene_id")),
      list(filters = "hgnc_symbol", attributes = c("hgnc_symbol", "ensembl_gene_id")),
      list(filters = "entrezgene_id", attributes = c("entrezgene_id", "ensembl_gene_id")),
      list(filters = "refseq_mrna", attributes = c("refseq_mrna", "ensembl_gene_id"))
    )
    
    # Perform queries and combine results
    cat("Querying Ensembl database...\n")
    results <- lapply(queries, function(q) {
      result <- getBM(filters = q$filters, attributes = q$attributes, values = unique_genes, mart = ensembl)
      setnames(result, c("input_id", "ensembl_gene_id"))
      return(result)
    })
    
    # Combine all results
    cat("Combining results...\n")
    all_results <- rbindlist(results, use.names = TRUE, fill = TRUE)
    
    # Remove duplicates and NA values
    all_results <- unique(all_results[!is.na(ensembl_gene_id)])
    
    # Create a named vector for easy mapping
    ensembl_map <- setNames(all_results$ensembl_gene_id, all_results$input_id)
    
    # Map the original genes to Ensembl IDs
    cat("Mapping genes to Ensembl IDs...\n")
    ensembl_ids <- ensembl_map[genes]
    
    # Return the mapped Ensembl IDs, keeping original IDs where no match was found
    ifelse(is.na(ensembl_ids), genes, ensembl_ids)
  }
  
  convert_to_entrez <- function(genes) {
    # Connect to the Ensembl database
    cat("Connecting to Ensembl database...\n")
    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    
    # Prepare a vector of unique gene IDs
    unique_genes <- unique(genes)
    cat("Number of unique genes:", length(unique_genes), "\n")
    
    # Create a list of attribute queries for different ID types
    queries <- list(
      list(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "entrezgene_id")),
      list(filters = "hgnc_symbol", attributes = c("hgnc_symbol", "entrezgene_id")),
      list(filters = "entrezgene_id", attributes = c("entrezgene_id", "entrezgene_id")),
      list(filters = "refseq_mrna", attributes = c("refseq_mrna", "entrezgene_id"))
    )
    
    # Perform queries and combine results
    cat("Querying Ensembl database...\n")
    results <- lapply(queries, function(q) {
      result <- getBM(filters = q$filters, attributes = q$attributes, values = unique_genes, mart = ensembl)
      setnames(result, c("input_id", "entrezgene_id"))
      return(result)
    })
    
    # Combine all results
    cat("Combining results...\n")
    all_results <- rbindlist(results, use.names = TRUE, fill = TRUE)
    
    # Remove duplicates and NA values
    all_results <- unique(all_results[!is.na(entrezgene_id)])
    
    # Create a named vector for easy mapping
    entrez_map <- setNames(all_results$entrezgene_id, all_results$input_id)
    
    # Map the original genes to Entrez IDs
    cat("Mapping genes to Entrez IDs...\n")
    entrez_ids <- entrez_map[genes]
    
    # Return the mapped Entrez IDs, keeping original IDs where no match was found
    ifelse(is.na(entrez_ids), genes, entrez_ids)
  }
  
  # Function to process and convert a dataframe
  process_dataframe <- function(df, name) {
    cat("Processing", name, "...\n")
    
    # Convert the first column to Ensembl IDs
    df[, V1 := convert_to_entrez(V1)]
    
    # Print the first few rows
    cat("First few rows of", name, "after conversion:\n")
    print(head(df))
    
    # Return the converted dataframe
    df
  }
  
  # Main processing function
  main_processing <- function() {
    # List of dataframes to process
    dfs <- list(
      ora_plus_enh_annotated = ora_plus_enh_annotated,
      ora_neg_enh_annotated = ora_neg_enh_annotated,
      ora_all_enh_annotated = ora_all_enh_annotated,
      gsea_annotated_enhancer = gsea_annotated_enhancer
    )
    
    # Process each dataframe
    processed_dfs <- lapply(names(dfs), function(name) {
      df <- dfs[[name]]
      cat("\nProcessing", name, "...\n")
      processed_df <- process_dataframe(df, name)
      
      # Optionally, save the updated dataframe
      output_file <- paste0(name, "_ensembl.txt")
      fwrite(processed_df, file = output_file, sep = "\t", quote = FALSE, col.names = FALSE)
      cat("Saved processed data to:", output_file, "\n")
      
      processed_df
    })
    
    # Name the processed dataframes
    names(processed_dfs) <- names(dfs)
    
    # Print summary of conversion
    cat("\nConversion summary:\n")
    for (name in names(processed_dfs)) {
      cat(name, ": ", nrow(processed_dfs[[name]]), "rows\n")
    }
    
    # Return the processed dataframes
    processed_dfs
  }
  
  # Add this function at the beginning of your script
  save_results_to_csv <- function(result, filename) {
    # Convert the result to a data frame
    result_df <- as.data.frame(result)
    
    # Write the data frame to a CSV file
    write.csv(result_df, file = paste0("./Output/GSEA/", filename, ".csv"), row.names = FALSE)
  }
  
  # Run the main processing function
  result <- main_processing()
  
  # Assign the results back to the original variable names
  ora_all_enh_annotated <- result$ora_all_enh_annotated
  ora_plus_enh_annotated <- result$ora_plus_enh_annotated
  ora_neg_enh_annotated <- result$ora_neg_enh_annotated
  gsea_annotated_enhancer <- result$gsea_annotated_enhancer
  
  # Drop any rows where first column is not an integer
  ora_all_enh_annotated <- ora_all_enh_annotated[grepl("^[0-9]+$", ora_all_enh_annotated$V1),]
  ora_plus_enh_annotated <- ora_plus_enh_annotated[grepl("^[0-9]+$", ora_plus_enh_annotated$V1),]
  ora_neg_enh_annotated <- ora_neg_enh_annotated[grepl("^[0-9]+$", ora_neg_enh_annotated$V1),]
  gsea_annotated_enhancer <- gsea_annotated_enhancer[grepl("^[0-9]+$", gsea_annotated_enhancer$V1),]
  #ora_plus_enh_annotated <- ora_plus_enh_annotated[grepl("^ENSG", ora_plus_enh_annotated$V1),]
  #ora_neg_enh_annotated <- ora_neg_enh_annotated[grepl("^ENSG", ora_neg_enh_annotated$V1),]
  #gsea_annotated_enhancer <- gsea_annotated_enhancer[grepl("^ENSG", gsea_annotated_enhancer$V1),]
  # save gsea_annotated_enhancer as tsv
  write.table(gsea_annotated_enhancer, file = "./Output/DESeq2/Results/gsea_annotated_enhancer.tsv", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  
  # Ensure the second column is numeric
  #gsea_annotated_enhancer$V2 <- as.numeric(gsea_annotated_enhancer$V2)
  
  # If there are duplicate V1s, aggregate the V2s by averaging
  #gsea_annotated_enhancer <- gsea_annotated_enhancer[, .(V2 = mean(V2)), by = V1]
  
  # Convert to a named vector (values from V2, names from V1)
  gsea_annotated_enhancer_vec <- setNames(gsea_annotated_enhancer$V2, gsea_annotated_enhancer$V1)
  
  # Sort the named vector in descending order by values
  sorted_gsea_annotated_enhancer_vec <- gsea_annotated_enhancer_vec[order(-gsea_annotated_enhancer_vec)]
  
  # ### REPEAT FOR GOBP
  # library(clusterProfiler)
  # library(DOSE)
  # library(patchwork)
  # library(ggplot2)
  # library(enrichplot)
  # library(data.table)
  # 
  # Helper functions
  safe_numeric <- function(x) {
    as.numeric(sub("/.*", "", x)) / as.numeric(sub(".*/", "", x))
  }
  
  create_combined_barplot <- function(pos_data, neg_data, top_pos = 50, top_neg = 50, x_label = "Percent of genes DHM") {
    
    total_lower = length(neg_data$p.adjust)
    total_upper = length(pos_data$p.adjust)
    total_rows = total_lower + total_upper
    percent_lower = total_lower / total_rows
    percent_upper = 1-percent_lower
    
    # Get full ranges for consistent scaling
    p_adjust_range <- range(c(pos_data$p.adjust, neg_data$p.adjust[1:top_neg]), na.rm = TRUE)
    gene_ratio_range <- range(c(safe_numeric(pos_data$GeneRatio), safe_numeric(neg_data$GeneRatio[1:top_neg])), na.rm = TRUE)
    
    # Positive bar plot
    p1 <- barplot(pos_data, showCategory = top_pos, x = "GeneRatio") +
      ggtitle("Barplot for ORA Positive") +
      xlab(NULL) +  
      scale_fill_gradient(low = "red", high = "blue", limits = p_adjust_range) +
      coord_cartesian(xlim = gene_ratio_range) +  
      theme(legend.position = "none", axis.text.y = element_text(size = rel(1.2)), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            panel.grid.minor = element_blank())
    
    # Negative bar plot
    p2 <- barplot(neg_data, showCategory = top_neg, x = "GeneRatio") +
      ggtitle("Barplot for ORA Negative") +
      xlab(x_label) +  
      scale_fill_gradient(low = "red", high = "blue", limits = p_adjust_range) +
      coord_cartesian(xlim = gene_ratio_range) +  
      theme(legend.position = "none", axis.text.y = element_text(size = rel(1.2)),
            panel.grid.minor = element_blank())
    
    # Combine plots & add legend
    combined_plot <- p1 / p2 + plot_layout(heights = c(percent_upper, percent_lower))
    legend <- get_legend(barplot(pos_data, showCategory = top_pos, x = "GeneRatio") + scale_fill_gradient(low = "red", high = "blue", limits = p_adjust_range) + guides(fill = guide_colorbar(title = "p.adjust")))
    final_plot <- combined_plot | legend
    
    return(final_plot)
  }

  save_plot <- function(plot, filename, width = 8, height = 6, dpi = 300) {
    # Save as PNG
    png(paste0("./Output/GSEA/", filename, ".png"), width = width, height = height, units = "in", res = dpi)
    print(plot)
    dev.off()
    
    # Save as SVG
    svg(paste0("./Output/GSEA/", filename, ".svg"), width = width, height = height)
    print(plot)
    dev.off()
  }
}

#### GSEA and ORA ANALYSIS and VIZ for Regions ----
if(FALSE){
  library(clusterProfiler)
  library(enrichplot)
  pval_cut = 0.1
  ora_pos_kegg <- enrichKEGG(gene = ora_plus_enh_annotated$V1, organism = 'hsa', keyType = "kegg", pAdjustMethod = "BH", minGSSize = 15, maxGSSize = 500, pvalueCutoff = pval_cut, qvalueCutoff = pval_cut)
  ora_neg_kegg <- enrichKEGG(gene = ora_neg_enh_annotated$V1, organism = 'hsa', keyType = "kegg", pAdjustMethod = "BH", minGSSize = 15, maxGSSize = 500, pvalueCutoff = pval_cut, qvalueCutoff = pval_cut)
  
  pval_cut = 0.2
  gsea1_kegg <- gseKEGG(geneList = sorted_gsea_annotated_enhancer_vec, organism = 'hsa', keyType = "kegg",minGSSize = 15, maxGSSize = 500, pvalueCutoff = pval_cut, pAdjustMethod = "BH",verbose = TRUE)
  
  # # GOBP Analysis
  # ora_pos_gobp <- enrichGO(gene = ora_plus_enh_annotated$V1, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = pval_cut, qvalueCutoff = pval_cut)
  # ora_neg_gobp <- enrichGO(gene = ora_neg_enh_annotated$V1, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = pval_cut, qvalueCutoff = pval_cut)
  # ora_all_gobp <- enrichGO(gene = ora_all_enh_annotated$V1, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = pval_cut, qvalueCutoff = pval_cut)
  # gsea1_gobp <- gseGO(geneList = sorted_gsea_annotated_enhancer_vec, OrgDb = org.Hs.eg.db, ont = "BP", keyType = "ENTREZID", minGSSize = 5, maxGSSize = 2000, pvalueCutoff = pval_cut, pAdjustMethod = "BH",verbose = TRUE)

  # ----- ORA Plots -----
  
  # GOBP Visualization
  # combined_gobp_plot <- create_combined_barplot(ora_pos_gobp, ora_neg_gobp, x_label = "Percent of genes DHM", top_neg = 10,top_pos = 10)
  # ora_pos_gobp_term <- pairwise_termsim(ora_pos_gobp)
  # ora_neg_gobp_term <- pairwise_termsim(ora_neg_gobp)
  # ora_pos_gobp_net <- emapplot(ora_pos_gobp_term,showCategory=5,layout = "kk",node_label_size=8,min_edge = 0.1) 
  # ora_neg_gobp_net <- emapplot(ora_neg_gobp_term,showCategory=5,layout = "kk",node_label_size=8,min_edge = 0.1)
  # 
  # KEGG Visualization
  combined_kegg_plot <- create_combined_barplot(ora_pos_kegg, ora_neg_kegg)
  
  ora_pos_kegg_term <- pairwise_termsim(ora_pos_kegg)
  ora_neg_kegg_term <- pairwise_termsim(ora_neg_kegg)
  ora_pos_kegg_net <- emapplot(ora_pos_kegg_term,showCategory=35,layout = "kk",node_label_size=8,min_edge = 0.25) 
  ora_neg_kegg_net <- emapplot(ora_neg_kegg_term,showCategory=10,layout = "kk",node_label_size=8,min_edge = 0.1)
  
  gsea1_kegg_read <- setReadable(gsea1_kegg, 'org.Hs.eg.db', 'ENTREZID')
  gsea1_kegg_term_all <- pairwise_termsim(gsea1_kegg)
  gsea1_kegg_term_all_net <- emapplot(gsea1_kegg_term_all,showCategory=45,layout = "kk",node_label_size=8,min_edge = 0.25) 
  
  # Separate the positive and negative enrichment results
  # Subset the gseaResult object for positive and negative NES values
  # Filter the GSEA results for positive and negative NES
  gsea1_kegg_pos <- gsea1_kegg %>% filter(NES > 0) %>% 
    slice(1:40)
  gsea1_kegg_neg <- gsea1_kegg %>% filter(NES < 0) %>% 
    slice(1:40)
  
  # Calculate pairwise term similarities for each subset
  gsea1_kegg_term_pos <- pairwise_termsim(gsea1_kegg_pos)
  gsea1_kegg_term_neg <- pairwise_termsim(gsea1_kegg_neg)
  gsea1_kegg_term_all <- pairwise_termsim(gsea1_kegg_read)
  
  
  gsea1_kegg_term_pos_net <- emapplot(gsea1_kegg_term_pos,showCategory=35,layout = "kk",node_label_size=8,min_edge = 0.2) 
  gsea1_kegg_term_neg_net <- emapplot(gsea1_kegg_term_neg,showCategory=15,layout = "kk",node_label_size=8,min_edge = 0.2) 
  gsea1_kegg_term_all_net <- emapplot(gsea1_kegg_term_all,showCategory=50,layout = "kk",node_label_size=8,min_edge = 0.2) 
  
  
  # Generate tree plots for positive and negative terms
  gsea_tree_pos <- treeplot(gsea1_kegg_term_pos, nCluster = 8,
                            cluster.params = list(method = "complete",color = NULL, label_words_n = 4,
                                                  label_format = 20),
                            showCategory = 50,hilight = FALSE)
  gsea_tree_neg <- treeplot(gsea1_kegg_term_neg,nCluster = 4,
                            cluster.params = list(method = "ward.D2",color = NULL, label_words_n = 5,
                                                  label_format = 30),showCategory = 50,hilight = FALSE)
  gsea_tree_all <- treeplot(gsea1_kegg_term_all, nCluster = 14,
                            cluster.params = list(method = "complete",color = NULL, label_words_n = 4,
                                                  label_format = 20),
                            showCategory = 60,hilight = FALSE)
  
  
  # SAVE FILES
  save_plot(gsea_tree_pos, width=10,height=8,paste0(contrast_groups[2], contrast_groups[3], "_", ver, "_GSEA_tree_pos"))
  save_plot(gsea_tree_neg, width=10,height=8,paste0(contrast_groups[2], contrast_groups[3], "_", ver, "_GSEA_tree_neg"))
  }

# ---- GSEA Plots  for regions---- 
if(FALSE){
  # gsea1_gobp_dotplot <- dotplot(gsea1_gobp, x = "NES", size = "GeneRatio", showCategory = 30) + 
  #   ggtitle("Dotplot for GSEA GOBP") + 
  #   geom_vline(xintercept = 0, linetype = "dashed", color = "black")
  # gsea1_kegg_dotplot <- dotplot(gsea1_kegg, x = "NES", size = "GeneRatio", showCategory = 30) + 
  #   ggtitle("Dotplot for GSEA KEGG") + 
  #   geom_vline(xintercept = 0, linetype = "dashed", color = "black")
  # returns n top results as tibble from result file
  return_top_n <- function(gsea1_kegg,n_neg,n_pos, padj_cutoff){ 
    # Convert the list to a tibble
    gsea1_kegg_df <- as_tibble(gsea1_kegg@result)
    
    # filter out rows where p.adjust > padj_cutoff
    gsea1_kegg_df <- gsea1_kegg_df %>% filter(p.adjust < padj_cutoff)
    
    # Filter and select top 10 rows where NES < 0
    negative_nes <- gsea1_kegg_df %>%
      filter(NES < 0) %>%
      arrange(p.adjust) %>%
      slice_head(n = n_neg)
    
    # Filter and select top 10 rows where NES > 0
    positive_nes <- gsea1_kegg_df %>%
      filter(NES > 0) %>%
      arrange(p.adjust) %>%
      slice_head(n = n_pos)
    
    # Combine the results
    top_ten_kegg <- bind_rows(negative_nes, positive_nes)
    
    # Sort the final result by p.adjust
    top_ten_kegg <- top_ten_kegg %>% arrange(p.adjust)
    
    # add geneRatio column calculated by count of genes in core_enrichment column (separated by /) divided by setSize
    top_ten_kegg$GeneRatio <- sapply(seq_len(nrow(top_ten_kegg)), function(i) {
      core_genes <- strsplit(top_ten_kegg$core_enrichment[i], "/")[[1]]
      length(core_genes) / top_ten_kegg$setSize[i]
    })
    
    return(top_ten_kegg)
  }
  
  top_kegg <- return_top_n(gsea1_kegg,15,15,0.025)
  # top_gobp <- return_top_n(gsea1_gobp,15,15,0.025)
  
  gsea_kegg_dot <- ggplot(top_kegg, aes(x = NES, y = reorder(Description, NES))) +
    geom_point(aes(color = p.adjust, size = GeneRatio)) +
    scale_color_gradient(low = "red", high = "blue") +
    labs(x = "Normalized Enrichment Score (NES)",
         y = "Pathway",
         color = "Adjusted p-value",
         size = "Gene Ratio") +
    scale_x_continuous(expand = expansion(mult = c(0.15, 0.15)))+
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8),
          #remove minor break lines on x axis
          panel.grid.minor.x = element_blank())
  
  # # Create the plot
  # gsea_gobp_dot <- ggplot(top_gobp, aes(x = NES, y = reorder(Description, NES))) +
  #   geom_point(aes(color = p.adjust, size = GeneRatio)) +
  #   scale_color_gradient(low = "red", high = "blue") +
  #   labs(x = "Normalized Enrichment Score (NES)",
  #        y = "Pathway",
  #        color = "Adjusted p-value",
  #        size = "Gene Ratio") +
  #   scale_x_continuous(expand = expansion(mult = c(0.15, 0.15)))+
  #   theme_minimal() +
  #   theme(axis.text.y = element_text(size = 8),
  #         panel.grid.minor = element_blank())
  
  gsea1_kegg_read <- setReadable(gsea1_kegg, 'org.Hs.eg.db', 'ENTREZID')
  gsea1_kegg_term <- pairwise_termsim(gsea1_kegg)
  gsea_kegg_net <- emapplot(gsea1_kegg_term,showCategory=30)
  
  # gsea1_gobp_read <- setReadable(gsea1_gobp, 'org.Hs.eg.db', 'ENTREZID')
  # gsea1_gobp_term <- pairwise_termsim(gsea1_gobp)
  # gsea_gobp_net <- emapplot(gsea1_gobp_term,showCategory=30)
}

# ----- Save region Plots -----
if(FALSE){
  
  
  # # GOBP Plots
  # save_plot(combined_gobp_plot, width=7,height=8,paste0(contrast_groups[2], contrast_groups[3], "_", ver, "_ora_pos_neg_gobp"))
  # 
  # # GSEA GOBP Dot Plot
  # save_plot(gsea_gobp_dot, width=6.5,height=6,paste0(contrast_groups[2], contrast_groups[3], "_", ver, "_gsea_gobp_dot"))
  # 
  # # GSEA GOBP Network Plot
  # save_plot(gsea_gobp_net, width=12,height=9,paste0(contrast_groups[2], contrast_groups[3], "_", ver, "_gsea_gobp_net"))
  
  # KEGG Plots
  save_plot(combined_kegg_plot, width=6,height=6,paste0(contrast_groups[2], contrast_groups[3], "_", ver, "_ora_pos_neg_kegg"))
  
  # GSEA KEGG Dot Plot
  save_plot(gsea_kegg_dot, width=6.5,height=5,paste0(contrast_groups[2], contrast_groups[3], "_", ver, "_gsea_kegg_dot"))
  
  # GSEA KEGG Network Plot
  save_plot(gsea_kegg_net, width=12,height=9,paste0(contrast_groups[2], contrast_groups[3], "_", ver, "_gsea_kegg_net"))
  
  # ORA KEGG Network Plots
  save_plot(ora_pos_kegg_net, paste0(contrast_groups[2], contrast_groups[3], "_", ver, "_ora_pos_kegg_net"))
  save_plot(ora_neg_kegg_net, paste0(contrast_groups[2], contrast_groups[3], "_", ver, "_ora_neg_kegg_net"))
  # 
  # # ORA GOBP Network Plots
  # save_plot(ora_pos_gobp_net, paste0(contrast_groups[2], contrast_groups[3], "_", ver, "_ora_pos_gobp_net"))
  # save_plot(ora_neg_gobp_net, paste0(contrast_groups[2], contrast_groups[3], "_", ver, "_ora_neg_gobp_net"))
  
  # Add this function at the beginning of your script
  save_results_to_csv <- function(result, filename) {
    # Convert the result to a data frame
    result_df <- as.data.frame(result)
    
    # Write the data frame to a CSV file
    write.csv(result_df, file = paste0("./Output/GSEA/", filename, ".csv"), row.names = FALSE)
  }
  
  # After your analysis sections, add these lines to export the results:
  
  # # Export GOBP results
  # save_results_to_csv(ora_pos_gobp, paste0(contrast_groups[2], contrast_groups[3], "_", ver, "_ora_pos_gobp"))
  # save_results_to_csv(ora_neg_gobp, paste0(contrast_groups[2], contrast_groups[3], "_", ver, "_ora_neg_gobp"))
  # save_results_to_csv(ora_all_gobp, paste0(contrast_groups[2], contrast_groups[3], "_", ver, "_ora_all_gobp"))
  # save_results_to_csv(gsea1_gobp, paste0(contrast_groups[2], contrast_groups[3], "_", ver, "_gsea_gobp"))
  
  # Export KEGG results
  save_results_to_csv(ora_pos_kegg, paste0(contrast_groups[2], contrast_groups[3], "_", ver, "_ora_pos_kegg"))
  save_results_to_csv(ora_neg_kegg, paste0(contrast_groups[2], contrast_groups[3], "_", ver, "_ora_neg_kegg"))
  save_results_to_csv(ora_all_kegg, paste0(contrast_groups[2], contrast_groups[3], "_", ver, "_ora_all_kegg"))
  save_results_to_csv(gsea1_kegg, paste0(contrast_groups[2], contrast_groups[3], "_", ver, "_gsea_kegg"))
  
  
}

### KEGG ORA BAR PLOTS AND GENE NETWORK PLOTS ----
if(FALSE){
  pval_cut = 0.1
  # KEGG Analysis
  ora_pos_kegg <- enrichKEGG(gene = ora_plus_enh_annotated$V1, organism = 'hsa', keyType = "kegg", pAdjustMethod = "BH", minGSSize = 15, maxGSSize = 500, pvalueCutoff = pval_cut, qvalueCutoff = pval_cut)
  ora_neg_kegg <- enrichKEGG(gene = ora_neg_enh_annotated$V1, organism = 'hsa', keyType = "kegg", pAdjustMethod = "BH", minGSSize = 15, maxGSSize = 500, pvalueCutoff = pval_cut, qvalueCutoff = pval_cut)
  ora_all_kegg <- enrichKEGG(gene = ora_all_enh_annotated$V1, organism = 'hsa', keyType = "kegg", pAdjustMethod = "BH", minGSSize = 15, maxGSSize = 500, pvalueCutoff = pval_cut, qvalueCutoff = pval_cut)
  # gsea1_kegg <- gseKEGG(geneList = sorted_gsea_annotated_enhancer_vec, organism = 'hsa', keyType = "kegg",minGSSize = 15, maxGSSize = 500, pvalueCutoff = pval_cut, pAdjustMethod = "BH",verbose = TRUE)
  
  # KEGG Visualization
  combined_kegg_plot <- create_combined_barplot(ora_pos_kegg, ora_neg_kegg)
  
  
  if(genebodies_type == TRUE){
    sig_annotated <- fread(fn_list$res_table, header = TRUE)
    # Convert gene symbols to Entrez IDs
    sig_annotated$entrez <- convert_to_entrez(sig_annotated$V1)
    # Filter rows based on pvalue < 0.01 and |log2FoldChange| > 0.2
    sig_annotated <- sig_annotated[pvalue < 0.05 & abs(log2FoldChange) > 0.13]
    log2fc_vector <- setNames(sig_annotated$log2FoldChange, sig_annotated$V1)
  } else {
    sig_annotated <- fread(fn_list$res_table_annotated_enhancer, header = TRUE)
    # Convert gene symbols to Entrez IDs
    sig_annotated$entrez <- convert_to_entrez(sig_annotated$symbols_gh_full)
    # Filter rows based on pvalue < 0.01 and |log2FoldChange| > 0.2
    sig_annotated <- sig_annotated[pvalue < 0.05 & abs(log2FoldChange) > 0.13]
    log2fc_vector <- setNames(sig_annotated$log2FoldChange, sig_annotated$symbols_gh_full)
  }
  
  
  # print names
  
  # make ora_pos and ora_neg geneID readable
  ora_pos_kegg_x <- setReadable(ora_pos_kegg, 'org.Hs.eg.db', 'ENTREZID')
  ora_neg_kegg_x <- setReadable(ora_neg_kegg, 'org.Hs.eg.db', 'ENTREZID')
  ora_all_kegg_x <- setReadable(ora_all_kegg, 'org.Hs.eg.db', 'ENTREZID')
  # 
  # # make ora_pos and ora_neg geneID readable for both GOBP and KEGG
  # ora_pos_gobp_x <- setReadable(ora_pos_gobp, 'org.Hs.eg.db', 'ENTREZID')
  # ora_neg_gobp_x <- setReadable(ora_neg_gobp, 'org.Hs.eg.db', 'ENTREZID')
  # ora_all_gobp_x <- setReadable(ora_all_gobp, 'org.Hs.eg.db', 'ENTREZID')
  
  # Load gene ids from ORA results by ora_pos_kegg_x$geneID, but splitting "/" into new rows
  ora_pos_kegg_x_gene_list <- ora_pos_kegg_x$geneID %>% strsplit("/") %>% unlist()
  ora_neg_kegg_x_gene_list <- ora_neg_kegg_x$geneID %>% strsplit("/") %>% unlist()
  ora_all_kegg_x_gene_list <- ora_all_kegg_x$geneID %>% strsplit("/") %>% unlist()
  
  # # Load gene ids from ORA results for both GOBP and KEGG
  # ora_pos_gobp_x_gene_list <- ora_pos_gobp_x$geneID %>% strsplit("/") %>% unlist()
  # ora_neg_gobp_x_gene_list <- ora_neg_gobp_x$geneID %>% strsplit("/") %>% unlist()
  # ora_all_gobp_x_gene_list <- ora_all_gobp_x$geneID %>% strsplit("/") %>% unlist()
  # 
  
  # drop duplicates
  ora_pos_kegg_x_gene_list <- unique(ora_pos_kegg_x_gene_list)
  ora_neg_kegg_x_gene_list <- unique(ora_neg_kegg_x_gene_list)
  ora_all_kegg_x_gene_list <- unique(ora_all_kegg_x_gene_list)
  
 
  
  # # drop duplicates
  # ora_pos_gobp_x_gene_list <- unique(ora_pos_gobp_x_gene_list)
  # ora_neg_gobp_x_gene_list <- unique(ora_neg_gobp_x_gene_list)
  # ora_all_gobp_x_gene_list <- unique(ora_all_gobp_x_gene_list)
  
  create_fc_vector <- function(gene_list, log2fc_vector) {
    # Create a named vector to store results
    result <- numeric(length(gene_list))
    names(result) <- gene_list
    
    # Create a lookup table for faster searching
    lookup_table <- log2fc_vector[!duplicated(names(log2fc_vector))]
    
    # Iterate through the gene list
    for (i in seq_along(gene_list)) {
      gene <- gene_list[i]
      
      # Look up the gene in the lookup table
      fc <- lookup_table[gene]
      
      # If not found, assign 0
      if (is.na(fc)) {
        fc <- 0
      }
      
      result[i] <- fc
    }
    
    return(result)
  }
  
  ### ORA GENE PLOTS 
  # Function to generate the cnetplot with customizations
  # Function to generate the cnetplot with customizations
  # Function to generate the cnetplot with customizations
  generate_cnetplot <- function(ora_x, description, log2fc_vector, num_path = 5, ora_type = "all", min_genes_per_set = 0, circ_param = TRUE,pval=0.25) {
    # Filter log2fc_vector based on genes present in ora_x
    genes_in_ora <- unique(unlist(strsplit(ora_x@result$geneID, "/")))
    log2fc_vector <- log2fc_vector[names(log2fc_vector) %in% genes_in_ora]
    
    # Drop duplicates
    log2fc_vector <- log2fc_vector[!duplicated(names(log2fc_vector))]
    
    # Convert log2fc to fold change based on ora_type
    if (ora_type == "pos") {
      ora_fc <- 2^log2fc_vector
      enrichment_type <- "Positive"
    } else if (ora_type == "neg") {
      ora_fc <- (2^log2fc_vector) #1 / (2^log2fc_vector)
      enrichment_type <- "Negative"
    } else {
      ora_fc <- ifelse(log2fc_vector >= 0, 2^log2fc_vector, 1 / (2^abs(log2fc_vector)))
      enrichment_type <- "All"
    }
    
    # Convert fold change to percentage
    ora_fc <- (ora_fc - 1) * 100
    
    # Apply scaling to p.adjust values for categories only
    scale_size <- function(x, min_size=0.5, max_size=4) {
      scaled_x <- (x - min(x)) / (max(x) - min(x))
      return(scaled_x * (max_size - min_size) + min_size)
    }
    
    # Filter the ora_x object based on the gene count
    gene_counts <- sapply(strsplit(ora_x@result$geneID, "/"), length)
    ora_x@result <- ora_x@result[gene_counts >= min_genes_per_set, ]
    
    category_sizes <- ora_x@result$p.adjust[ora_x@result$p.adjust < pval]
    actual_p_adjust <- ora_x@result$p.adjust[ora_x@result$p.adjust < pval]
    
    if (length(category_sizes) > num_path) {
      category_sizes <- category_sizes[1:num_path]
      actual_p_adjust <- actual_p_adjust[1:num_path]
    }
    
    category_sizes <- scale_size(category_sizes)
    
    legend_breaks <- category_sizes
    legend_labels <- actual_p_adjust
    legend_labels <- rev(legend_labels)
    
    # Determine color scheme based on fold change values
    if (all(ora_fc >= 0)) {
      low_color <- "#FFEBEE"
      high_color <- "#B71C1C"
    } else if (any(ora_fc >= 0) && any(ora_fc < 0)) {
      low_color <- "#0D47A1"
      high_color <- "#B71C1C"
    } else {
      low_color <- "#0D47A1"
      high_color <- "#E3F2FD"
    }
    
    # Modify the cnetplot call with updated parameter usage
    p <- cnetplot(
      ora_x,
      foldChange=ora_fc,
      showCategory = length(category_sizes),
      layout = "kk",
      circular = circ_param,
      colorEdge = TRUE,
      node_label = "all",
      cex_gene = 2,
      cex_label_gene = 0,
      cex_label_category = 1,
      color_category = 'grey30'
    )
    
    adjusted_breaks <- legend_breaks[1:length(category_sizes)]
    adjusted_labels <- legend_labels[1:length(category_sizes)]
    
    adjusted_labels <- round(adjusted_labels, 3)
    
    adjusted_breaks <- as.numeric(adjusted_breaks)
    adjusted_labels <- as.character(round(as.numeric(adjusted_labels), 3))
    
    size_scale <- adjusted_breaks
    names(size_scale) <- adjusted_labels
    
    # Create the plot title
    plot_title <- paste(enrichment_type, "Enriched Gene Set Network")
    
    p_ggraph <- p +
      scale_color_gradient(low = low_color, high = high_color, name = "% Fold Change") +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(color = "black"),
        legend.title = element_text(color = "black"),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16)  # Center the title and make it bold
      ) +
      ggtitle(plot_title)  # Add the title to the plot
    
    p_adjust_values <- setNames(ora_x@result$p.adjust, ora_x@result$Description)
    
    p_ggraph$data$p_adjust <- p_adjust_values[p_ggraph$data$name]
    
    p_ggraph <- p_ggraph + 
      scale_size_continuous(
        name = "DHGs in Pathway"
      )
    
    p_ggraph$layers[[1]]$aes_params$edge_width <- 0.25
    p_ggraph$layers[[1]]$aes_params$edge_alpha <- 0.5
    
    plot_data <- ggplot_build(p_ggraph)
    
    p_ggraph <- p_ggraph + 
      geom_text_repel(
        data = plot_data$data[[4]],
        aes(x = plot_data$data[[4]]$x, y = plot_data$data[[4]]$y, label = plot_data$data[[4]]$label),
        size = 4,
        max.overlaps = 10,
        fontface = "bold"
      )
    
    return(p_ggraph)
  }
  #ora_fig <- generate_cnetplot(ora_all_x, ora_all_x_gene_list, log2fc_vector)
  
  ora_pos_kegg_fig <- generate_cnetplot(ora_pos_kegg_x, ora_pos_kegg_x_gene_list, log2fc_vector,num_path=10, ora_type = "pos", min_genes_per_set = 0, circ_param = TRUE,pval=0.1)
  ora_neg_kegg_fig <- generate_cnetplot(ora_neg_kegg_x, ora_neg_kegg_x_gene_list, log2fc_vector,num_path=10, ora_type = "neg", min_genes_per_set = 0, circ_param = TRUE,pval=0.1)
  
  ora_pos_kegg_fig_net <- generate_cnetplot(ora_pos_kegg_x, ora_pos_kegg_x_gene_list, log2fc_vector,num_path=10, ora_type = "pos", min_genes_per_set = 0, circ_param = FALSE)
  ora_neg_kegg_fig_net <- generate_cnetplot(ora_neg_kegg_x, ora_neg_kegg_x_gene_list, log2fc_vector,num_path=10, ora_type = "neg", min_genes_per_set = 0, circ_param = FALSE)
  
  
  # ORA KEGG Network Plots
  save_plot(ora_pos_kegg_fig_net, 
            paste0(contrast_groups[2], contrast_groups[3], "_", ver, "_ora_pos_kegg_net"),width = 8, height = 6, dpi = 300)
  save_plot(ora_neg_kegg_fig_net, 
            paste0(contrast_groups[2], contrast_groups[3], "_", ver, "_ora_neg_kegg_net"),width = 10, height = 10, dpi = 300)
  save_plot(combined_kegg_plot,
            paste0(contrast_groups[2], contrast_groups[3], "_", ver, "_ora_DotPlot_kegg"),width = 7, height = 7, dpi = 300)
  save_results_to_csv(ora_pos_kegg, paste0(contrast_groups[2], contrast_groups[3], "_", ver, "_ORApos_kegg"))
  save_results_to_csv(ora_neg_kegg, paste0(contrast_groups[2], contrast_groups[3], "_", ver, "_ORAneg_kegg"))
  
  
  # # Generate GOBP plots
  # ora_pos_gobp_fig <- generate_cnetplot(ora_pos_gobp_x, ora_pos_gobp_x_gene_list, log2fc_vector, num_path=15, ora_type = "pos", min_genes_per_set = 0, circ_param = TRUE)
  # ora_neg_gobp_fig <- generate_cnetplot(ora_neg_gobp_x, ora_neg_gobp_x_gene_list, log2fc_vector, num_path=5, ora_type = "neg", min_genes_per_set = 0, circ_param = TRUE)
  # 
  # ora_pos_gobp_fig_net <- generate_cnetplot(ora_pos_gobp_x, ora_pos_gobp_x_gene_list, log2fc_vector, num_path=15, ora_type = "pos", min_genes_per_set = 0, circ_param = FALSE)
  # ora_neg_gobp_fig_net <- generate_cnetplot(ora_neg_gobp_x, ora_neg_gobp_x_gene_list, log2fc_vector, num_path=4, ora_type = "neg", min_genes_per_set = 0, circ_param = FALSE)
}

### GSEA GENE PLOTS VERSION 1 ----
if(FALSE){
  
  # read in sig_annotated file
  sig_annotated <- fread(fn_list$res_table_annotated_enhancer, header = TRUE)
  
  # Filter rows based on pvalue < 0.01 and |log2FoldChange| > 0.2
  sig_annotated <- sig_annotated[pvalue < 0.01 & abs(log2FoldChange) > 0.2]
  
  # Convert gene symbols to Entrez IDs
  sig_annotated$entrez <- convert_to_entrez(sig_annotated$symbols_gh_full)
  
  # Create a named vector of log2FoldChange values
  log2fc_vector <- setNames(sig_annotated$log2FoldChange, sig_annotated$symbols_gh_full)
  
  # Load necessary libraries
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(stringr)
  
  # Define your p-value cutoff
  pval_cut <- 0.1  # Adjust as needed
  
  # Make GSEA results readable
  gsea1_kegg_x <- setReadable(gsea1_kegg, 'org.Hs.eg.db', 'ENTREZID')
  gsea1_gobp_x <- setReadable(gsea1_gobp, 'org.Hs.eg.db', 'ENTREZID')
  
  # Function to extract gene lists from enrichResult or gseaResult
  extract_gene_list <- function(ora_x) {
    if (inherits(ora_x, "gseaResult")) {
      # For GSEA, extract core_enrichment
      gene_lists <- ora_x@result$core_enrichment
    } else {
      # For ORA, extract geneID
      gene_lists <- ora_x@result$geneID
    }
    # Split by "/" and unlist
    gene_list <- unique(unlist(strsplit(gene_lists, "/")))
    return(gene_list)
  }
  
  # Extract gene lists for GSEA
  gsea1_kegg_x_gene_list <- extract_gene_list(gsea1_kegg_x)
  gsea1_gobp_x_gene_list <- extract_gene_list(gsea1_gobp_x)
  
  # Function to create a fold change vector
  create_fc_vector <- function(gene_list, log2fc_vector) {
    # Initialize a numeric vector with zeros
    result <- numeric(length(gene_list))
    names(result) <- gene_list
    
    # Assign fold change values where available
    matched_genes <- intersect(names(log2fc_vector), gene_list)
    result[matched_genes] <- log2fc_vector[matched_genes]
    
    return(result)
  }
  
  gsea_kegg_fc_vector <- create_fc_vector(gsea1_kegg_x_gene_list, log2fc_vector)
  gsea_gobp_fc_vector <- create_fc_vector(gsea1_gobp_x_gene_list, log2fc_vector)
  
  # drop  items == to 0 from gsea_kegg_fc_vector
  gsea_kegg_fc_vector <- gsea_kegg_fc_vector[gsea_kegg_fc_vector != 0]
  
  # Function to generate the cnetplot with customizations
  generate_cnetplot <- function(ora_x, log2fc_vector, num_path = 5, ora_type = "all", min_genes_per_set = 0, circ_param = TRUE, pval = 0.25) {
    # Determine if the input is GSEA or ORA
    is_gsea <- inherits(ora_x, "gseaResult")
    
    # Extract gene lists based on the type
    if (is_gsea) {
      gene_column <- "core_enrichment"
      gene_list <- extract_gene_list(ora_x)
    } else {
      gene_column <- "geneID"
      gene_list <- unique(unlist(strsplit(ora_x$geneID, "/")))
    }
    
    # Filter log2fc_vector based on genes present in ora_x
    log2fc_filtered <- log2fc_vector[names(log2fc_vector) %in% gene_list]
    
    # Drop duplicates
    log2fc_filtered <- log2fc_filtered[!duplicated(names(log2fc_filtered))]
    
    # Convert log2fc to fold change based on ora_type
    if (ora_type == "pos") {
      ora_fc <- 2^log2fc_filtered
      enrichment_type <- "Positive"
    } else if (ora_type == "neg") {
      ora_fc <- 1 / (2^log2fc_filtered)
      enrichment_type <- "Negative"
    } else {
      ora_fc <- ifelse(log2fc_filtered >= 0, 2^log2fc_filtered, 1 / (2^abs(log2fc_filtered)))
      enrichment_type <- "All"
    }
    
    # Convert fold change to percentage
    ora_fc <- (ora_fc - 1) * 100
    
    # Apply scaling to p.adjust values for categories only
    scale_size <- function(x, min_size = 0.5, max_size = 4) {
      scaled_x <- (x - min(x)) / (max(x) - min(x))
      return(scaled_x * (max_size - min_size) + min_size)
    }
    
    # Filter the ora_x object based on the gene count
    if (is_gsea) {
      gene_counts <- sapply(strsplit(ora_x@result$core_enrichment, "/"), length)
    } else {
      gene_counts <- sapply(strsplit(ora_x$geneID, "/"), length)
    }
    ora_x@result <- ora_x@result[gene_counts >= min_genes_per_set, ]
    
    # Filter based on p-value cutoff
    ora_x@result <- ora_x@result[ora_x@result$p.adjust < pval, ]
    
    # Limit to the top pathways
    if (nrow(ora_x@result) > num_path) {
      ora_x@result <- ora_x@result[1:num_path, ]
    }
    
    # Prepare category sizes for scaling
    category_sizes <- ora_x@result$p.adjust
    category_sizes <- scale_size(category_sizes)
    
    # Determine color scheme based on fold change values
    if (all(ora_fc >= 0)) {
      low_color <- "#FFEBEE"
      high_color <- "#B71C1C"
    } else if (any(ora_fc >= 0) && any(ora_fc < 0)) {
      low_color <- "#0D47A1"
      mid_color <- "white"
      high_color <- "#B71C1C"
    } else {
      low_color <- "#0D47A1"
      high_color <- "#E3F2FD"
    }
    
    # Generate the cnetplot
    p <- cnetplot(
      ora_x,
      foldChange = ora_fc,
      showCategory = num_path,
      layout = "kk",
      circular = circ_param,
      colorEdge = TRUE,
      node_label = "all",
      cex_gene = 2,
      cex_label_gene = 0,
      cex_label_category = 1,
      color_category = 'grey30'
    )
    
    # Create the plot title
    plot_title <- paste(enrichment_type, "Enriched Gene Set Network")
    
    #  Customize the plot
    p_ggraph <- p +
      scale_color_gradient2(low = low_color, mid = mid_color, high = high_color, name = "% Fold Change") +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(color = "black"),
        legend.title = element_text(color = "black"),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
      ) +
      ggtitle(plot_title)
    
    # Add p.adjust values to the nodes
    p_adjust_values <- setNames(ora_x@result$p.adjust, ora_x@result$Description)
    p_ggraph$data$p_adjust <- p_adjust_values[p_ggraph$data$name]
    
    # Adjust edge aesthetics if present
    if(length(p_ggraph$layers) >= 1){
      p_ggraph$layers[[1]]$aes_params$edge_width <- 0.25
      p_ggraph$layers[[1]]$aes_params$edge_alpha <- 0.5
    }
    
    # Check if 'x' and 'y' columns exist
    if(!all(c("x", "y") %in% colnames(p_ggraph$data))) {
      stop("The plot data does not contain 'x' and 'y' columns required for geom_text_repel.")
    }
    
    # Add labels using geom_text_repel for all nodes
    p_ggraph <- p_ggraph + 
      geom_text_repel(
        data = subset(p_ggraph$data, !is.na(name)),  # Ensure that we have valid names for labeling
        aes(x = x, y = y, label = name),
        size = 4,
        max.overlaps = 10,
        fontface = "bold"
      )
    
    return(p_ggraph)
  }
  
  # Load necessary libraries
  library(dplyr)
  library(clusterProfiler)
  
  # Extract the names of the genes in the filtered fold change vector
  filtered_gene_names <- names(gsea_kegg_fc_vector)
  
  # Extract the result data frame from the S4 object
  gsea_data_df <- as.data.frame(gsea1_kegg_x)
  
  # Function to filter core enrichment genes based on the filtered genes
  filter_core_enrichment <- function(gsea_data, filtered_genes) {
    # Loop through each row in the data frame and filter core enrichment
    gsea_data$core_enrichment <- sapply(gsea_data$core_enrichment, function(genes_str) {
      # Split the core enrichment genes
      genes <- unlist(strsplit(genes_str, "/"))
      # Keep only the genes that are in the filtered gene list
      filtered_core_genes <- genes[genes %in% filtered_genes]
      # Collapse the genes back into a single string
      paste(filtered_core_genes, collapse = "/")
    })
    return(gsea_data)
  }
  
  # Apply the filtering function to the data frame
  gsea_data_df_filtered <- filter_core_enrichment(gsea_data_df, filtered_gene_names)
  
  # Update the original S4 object with the filtered data frame
  # We need to use the @ operator to access and modify the slot properly
  gsea1_kegg_x@result <- gsea_data_df_filtered
  
  # Generate GSEA KEGG Plot with the modified object
  gsea_kegg_fig <- generate_cnetplot(
    ora_x = gsea1_kegg_x,
    log2fc_vector = gsea_kegg_fc_vector,
    num_path = 20,           # Adjust the number of pathways as needed
    ora_type = "all",        # Can be "pos", "neg", or "all"
    min_genes_per_set = 0,   # Minimum number of genes per set
    circ_param = FALSE,       # Circular layout
    pval = 0.2               # P-value cutoff
  )
  
  # # Generate GSEA GOBP Plot
  # gsea_gobp_fig <- generate_cnetplot(
  #   ora_x = gsea1_gobp_x,
  #   log2fc_vector = gsea_gobp_fc_vector,
  #   num_path = 10,           # Adjust the number of pathways as needed
  #   ora_type = "all",        # Can be "pos", "neg", or "all"
  #   min_genes_per_set = 0,   # Minimum number of genes per set
  #   circ_param = TRUE,       # Circular layout
  #   pval = 0.2                # P-value cutoff
  # )
  
  print(gsea_kegg_fig)
  #print(gsea_gobp_fig)
}

### GSEA GENE PLOT VERSION 2 ----
if(TRUE){
  library(clusterProfiler)
  library(ggraph)
  
  pval_cut = 0.00001
  gsea1_kegg_2 <- gseKEGG(geneList = sorted_gsea_annotated_enhancer_vec, organism = 'hsa', keyType = "kegg", minGSSize = 15, maxGSSize = 500, pvalueCutoff = pval_cut, pAdjustMethod = "BH",verbose = TRUE)
  # Make GSEA results readable
  gsea1_kegg_2 <- setReadable(gsea1_kegg_2, 'org.Hs.eg.db', 'ENTREZID')
  
  return_top_n <- function(gsea1_kegg,n_neg,n_pos, padj_cutoff){ 
    # Convert the list to a tibble
    gsea1_kegg_df <- as_tibble(gsea1_kegg@result)
    
    # filter out rows where p.adjust > padj_cutoff
    gsea1_kegg_df <- gsea1_kegg_df %>% filter(p.adjust < padj_cutoff)
    
    # Filter and select top 10 rows where NES < 0
    negative_nes <- gsea1_kegg_df %>%
      filter(NES < 0) %>%
      arrange(p.adjust) %>%
      slice_head(n = n_neg)
    
    # Filter and select top 10 rows where NES > 0
    positive_nes <- gsea1_kegg_df %>%
      filter(NES > 0) %>%
      arrange(p.adjust) %>%
      slice_head(n = n_pos)
    
    # Combine the results
    top_ten_kegg <- bind_rows(negative_nes, positive_nes)
    
    # Sort the final result by p.adjust
    top_ten_kegg <- top_ten_kegg %>% arrange(p.adjust)
    
    # add geneRatio column calculated by count of genes in core_enrichment column (separated by /) divided by setSize
    top_ten_kegg$GeneRatio <- sapply(seq_len(nrow(top_ten_kegg)), function(i) {
      core_genes <- strsplit(top_ten_kegg$core_enrichment[i], "/")[[1]]
      length(core_genes) / top_ten_kegg$setSize[i]
    })
    
    return(top_ten_kegg)
  }
  
  if(genebodies_type == TRUE){
    # read in sig_annotated file
    sig_annotated <- fread(fn_list$res_table, header = TRUE)
    # Remove rows where gene symbol is NA
    sig_annotated <- sig_annotated[!is.na(V1) & V1 != "NA"]
    # Filter rows based on pvalue < 0.01 and |log2FoldChange| > 0.2
    #sig_annotated <- sig_annotated[pvalue < 0.05 & abs(log2FoldChange) > 0.136]
    # select 400 rows with highest abs(stat) value
    sig_annotated <- sig_annotated[order(abs(stat), decreasing = TRUE),]
    sig_annotated <- sig_annotated[1:1000,]
    # Convert gene symbols to Entrez IDs
    sig_annotated$entrez <- convert_to_entrez(sig_annotated$V1)
    # Create a named vector of log2FoldChange values
    log2fc_vector <- setNames(sig_annotated$log2FoldChange, sig_annotated$V1)
    pval_vector <- setNames(sig_annotated$pvalue, sig_annotated$V1)
  } else{
    # read in sig_annotated file
    sig_annotated <- fread(fn_list$res_table_annotated_enhancer, header = TRUE)
    # Remove rows where gene symbol is NA
    sig_annotated <- sig_annotated[!is.na(symbols_gh_full) & symbols_gh_full != "NA"]
    # Filter rows based on pvalue < 0.01 and |log2FoldChange| > 0.2
    sig_annotated <- sig_annotated[pvalue < 0.05 & abs(log2FoldChange) > 0.136]
    # Convert gene symbols to Entrez IDs
    sig_annotated$entrez <- convert_to_entrez(sig_annotated$symbols_gh_full)
    # Create a named vector of log2FoldChange values
    log2fc_vector <- setNames(sig_annotated$log2FoldChange, sig_annotated$V1)
    pval_vector <- setNames(sig_annotated$pvalue, sig_annotated$V1)
  }
  
  # Function to extract gene lists from enrichResult or gseaResult and keep pathway info
  extract_gene_list_with_pathway <- function(ora_x) {
    gene_list_pathway <- list()
    if (inherits(ora_x, "gseaResult")) {
      gene_column <- "core_enrichment"
    } else {
      gene_column <- "geneID"
    }
    
    for (i in 1:nrow(ora_x@result)) {
      pathway_name <- ora_x@result$Description[i]
      genes <- unlist(strsplit(ora_x@result[[gene_column]][i], "/"))
      gene_list_pathway[[pathway_name]] <- genes
    }
    
    return(gene_list_pathway)
  }
  
  # Extract gene lists with pathways for GSEA results
  gsea1_kegg_x_gene_list_with_pathway <- extract_gene_list_with_pathway(gsea1_kegg_2)
  
  # Function to filter top 5 genes with highest |log2FoldChange| for each pathway
  filter_top_genes_per_pathway <- function(gene_list_pathway, log2fc_vector, top_n = 5) {
    filtered_gene_list <- list()
    
    for (pathway in names(gene_list_pathway)) {
      genes <- gene_list_pathway[[pathway]]
      
      # Get the log2 fold change values for these genes
      gene_fc_values <- log2fc_vector[genes]
      gene_fc_values <- gene_fc_values[!is.na(gene_fc_values)]  # Remove NAs
      
      # Sort by absolute value of log2 fold change and select the top N genes
      top_genes <- names(sort(abs(gene_fc_values), decreasing = TRUE))[1:min(top_n, length(gene_fc_values))]
      
      # Store the filtered gene list for the pathway
      filtered_gene_list[[pathway]] <- top_genes
    }
    
    return(filtered_gene_list)
  }
  
  # Apply the filtering to keep top 5 genes per pathway
  filtered_gene_list_pathway <- filter_top_genes_per_pathway(gsea1_kegg_x_gene_list_with_pathway, 
                                                             log2fc_vector, 
                                                             top_n = 5)
  
  # Function to update core enrichment in gseaResult or enrichResult
  update_core_enrichment <- function(ora_x, filtered_gene_list_pathway) {
    if (inherits(ora_x, "gseaResult")) {
      gene_column <- "core_enrichment"
    } else {
      gene_column <- "geneID"
    }
    
    for (i in 1:nrow(ora_x@result)) {
      pathway_name <- ora_x@result$Description[i]
      if (pathway_name %in% names(filtered_gene_list_pathway)) {
        ora_x@result[[gene_column]][i] <- paste(filtered_gene_list_pathway[[pathway_name]], collapse = "/")
      }
    }
    
    return(ora_x)
  }
  
  # Load necessary libraries
  library(clusterProfiler)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  
  # Function to generate the cnetplot with customizations
  generate_cnetplot <- function(ora_x, log2fc_vector, pval_vector=NULL, num_path = 5, ora_type = "all", min_genes_per_set = 0, circ_param = TRUE, pval = 0.25) {
    # Determine if the input is GSEA or ORA
    is_gsea <- inherits(ora_x, "gseaResult")
    
    # Filter pathways based on NES score if available and requested
    if (is_gsea && ora_type %in% c("pos", "neg")) {
      if ("NES" %in% colnames(ora_x@result)) {
        if (ora_type == "pos") {
          ora_x@result <- ora_x@result[ora_x@result$NES > 0, ]
          enrichment_type <- "Positive"
        } else if (ora_type == "neg") {
          ora_x@result <- ora_x@result[ora_x@result$NES < 0, ]
          enrichment_type <- "Negative"
        }
      } else {
        stop("NES score not found in GSEA result. Unable to filter pathways by NES.")
      }
    } else {
      enrichment_type <- "All"
    }
    
    # Extract gene lists based on the type
    if (is_gsea) {
      gene_column <- "core_enrichment"
      gene_list <- unique(unlist(strsplit(ora_x@result$core_enrichment, "/")))
    } else {
      gene_column <- "geneID"
      gene_list <- unique(unlist(strsplit(ora_x$geneID, "/")))
    }
    
    # Filter log2fc_vector based on genes present in ora_x
    log2fc_filtered <- log2fc_vector[names(log2fc_vector) %in% gene_list]
    
    # Drop duplicates
    log2fc_filtered <- log2fc_filtered[!duplicated(names(log2fc_filtered))]
    
    # Convert log2fc to fold change for plotting
    ora_fc <- ifelse(log2fc_filtered >= 0, 2^log2fc_filtered, 1 / (2^abs(log2fc_filtered)))
    
    # Convert fold change to percentage
    ora_fc <- (ora_fc - 1) * 100
    
    # Apply scaling to p.adjust values for categories only
    scale_size_category <- function(x, min_size = 0.5, max_size = 4) {
      scaled_x <- (x - min(x)) / (max(x) - min(x))
      return(scaled_x * (max_size - min_size) + min_size)
    }
    
    # Filter the ora_x object based on the gene count
    if (is_gsea) {
      gene_counts <- sapply(strsplit(ora_x@result$core_enrichment, "/"), length)
    } else {
      gene_counts <- sapply(strsplit(ora_x$geneID, "/"), length)
    }
    ora_x@result <- ora_x@result[gene_counts >= min_genes_per_set, ]
    
    # Filter based on p-value cutoff
    ora_x@result <- ora_x@result[ora_x@result$p.adjust < pval, ]
    
    # Limit to the top pathways
    if (nrow(ora_x@result) > num_path) {
      ora_x@result <- ora_x@result[1:num_path, ]
    }
    
    # Prepare category sizes for scaling (optional, not used in size scaling for gene nodes)
    category_sizes <- ora_x@result$p.adjust
    category_sizes <- scale_size_category(category_sizes)
    
    # Determine color scheme based on fold change values
    if (all(ora_fc >= 0)) {
      low_color <- "white"
      mid_color <- "#FFEBEE"
      high_color <- "#B71C1C"
    } else if (any(ora_fc >= 0) && any(ora_fc < 0)) {
      low_color <- "#0D47A1"
      mid_color <- "white"
      high_color <- "#B71C1C"
    } else {
      low_color <- "#0D47A1"
      mid_color <- "#E3F2FD"
      high_color <- "white"
    }
    
    # Generate the cnetplot
    p <- cnetplot(
      ora_x,
      foldChange = ora_fc,
      showCategory = num_path,
      layout = "kk",
      circular = circ_param,
      colorEdge = FALSE,  # Set colorEdge to FALSE to make all edges the same color
      edge_color = "grey",  # Set edges to grey
      node_label = "all",
      cex_gene = 1,
      cex_category = 1,
      cex_label_gene = 0,
      cex_label_category = 0
    )
    
    # Get the pathway names
    if ("Description" %in% colnames(ora_x@result)) {
      pathway_names <- ora_x@result$Description
    } else {
      pathway_names <- ora_x@result$ID
    }
    
    # Get gene names
    gene_names <- names(log2fc_filtered)
    
    # If pval_vector is provided, filter it
    if (!is.null(pval_vector)) {
      pval_filtered <- pval_vector[names(pval_vector) %in% gene_names]
      pval_filtered <- pval_filtered[!duplicated(names(pval_filtered))]
    } else {
      # If pval_vector is not provided, set default p-value
      pval_filtered <- rep(0.05, length(gene_names))
      names(pval_filtered) <- gene_names
    }
    
    # Compute -log10(pval) for gene nodes and cap between 1 and 4
    gene_logpvals <- -log10(pval_filtered)
    gene_logpvals <- pmin(pmax(gene_logpvals, 0), 5)  # Cap between 1 and 4
    
    # Set pathway node size (constant)
    pathway_node_size <- 6
    
    # Assign node types
    p$data$node_type <- ifelse(p$data$name %in% pathway_names, "pathway",
                               ifelse(p$data$name %in% gene_names, "gene", NA))
    
    # Initialize pval and logpval columns
    p$data$pval <- NA
    p$data$logpval <- NA
    
    # Assign pval and logpval to gene nodes
    gene_indices <- which(p$data$node_type == "gene")
    gene_names_in_pdata <- p$data$name[gene_indices]
    
    p$data$pval[gene_indices] <- pval_filtered[gene_names_in_pdata]
    p$data$logpval[gene_indices] <- gene_logpvals[gene_names_in_pdata]
    
    # Assign size variable: genes use logpval, pathways use NA
    p$data$size_var <- NA
    p$data$size_var[gene_indices] <- p$data$logpval[gene_indices]
    # Pathway nodes will be plotted separately with fixed size
    
    # Assign label properties
    p$data$label_color <- NA
    p$data$label_fontface <- NA
    p$data$label_color[p$data$node_type == "gene"] <- "grey50"
    p$data$label_fontface[p$data$node_type == "gene"] <- "plain"
    p$data$label_color[p$data$node_type == "pathway"] <- "black"
    p$data$label_fontface[p$data$node_type == "pathway"] <- "bold"
    
    # Create the plot title
    plot_title <- paste(enrichment_type, "Enriched Gene Set Network")
    
    # Customize the plot
    p_ggraph <- p +
      scale_color_gradient2(low = low_color, mid = mid_color, high = high_color, name = "% Fold Change") +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(color = "black"),
        legend.title = element_text(color = "black"),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
      ) +
      ggtitle(plot_title)
    
    # Subset data for gene and pathway nodes
    gene_nodes <- subset(p$data, node_type == "gene")
    pathway_nodes <- subset(p$data, node_type == "pathway")
    
    # Plot pathway nodes with fixed size
    p_ggraph <- p_ggraph +
      geom_node_point(aes(x = x, y = y, color = color),
                      data = pathway_nodes,
                      size = pathway_node_size,
                      show.legend = FALSE)  # Hide pathway nodes from the legend
    
    # Plot gene nodes with size mapped to logpval
    p_ggraph <- p_ggraph +
      geom_node_point(aes(x = x, y = y, size = logpval, color = color),
                      data = gene_nodes) +
      scale_size_continuous(
        name = "-log10(p-value)",
        breaks = 1:4,
        labels = 1:4,
        limits = c(1, 4),
        guide = guide_legend(order = 1)
      )
    
    # Add gene labels in grey
    p_ggraph <- p_ggraph +
      geom_text_repel(
        data = gene_nodes,
        aes(x = x, y = y, label = name),
        size = 3,
        color = "grey50",
        max.overlaps = 20
      )
    
    # Add pathway labels in black and bold
    p_ggraph <- p_ggraph +
      geom_text_repel(
        data = pathway_nodes,
        aes(x = x, y = y, label = name),
        size = 3,
        fontface = "bold",
        color = "black",
        max.overlaps = 20
      )
    
    return(p_ggraph)
  }
  
  
  
  # Function to generate the cnetplot with customizations
  generate_cnetplot_spec <- function(ora_x, log2fc_vector, num_path = 5, ora_type = "all", min_genes_per_set = 0, circ_param = TRUE, pval = 0.25, selected_pathways = NULL) {
    # Determine if the input is GSEA or ORA
    is_gsea <- inherits(ora_x, "gseaResult")
    
    # Filter pathways based on NES score if available and requested
    if (is_gsea && ora_type %in% c("pos", "neg")) {
      if ("NES" %in% colnames(ora_x@result)) {
        if (ora_type == "pos") {
          ora_x@result <- ora_x@result[ora_x@result$NES > 0, ]
          enrichment_type <- "Positive"
        } else if (ora_type == "neg") {
          ora_x@result <- ora_x@result[ora_x@result$NES < 0, ]
          enrichment_type <- "Negative"
        }
      } else {
        stop("NES score not found in GSEA result. Unable to filter pathways by NES.")
      }
    } else {
      enrichment_type <- "All"
    }
    
    # Extract gene lists based on the type
    if (is_gsea) {
      gene_column <- "core_enrichment"
      gene_list <- extract_gene_list(ora_x)
    } else {
      gene_column <- "geneID"
      gene_list <- unique(unlist(strsplit(ora_x$geneID, "/")))
    }
    
    # Filter log2fc_vector based on genes present in ora_x
    log2fc_filtered <- log2fc_vector[names(log2fc_vector) %in% gene_list]
    
    # Drop duplicates
    log2fc_filtered <- log2fc_filtered[!duplicated(names(log2fc_filtered))]
    
    # Convert log2fc to fold change for plotting
    ora_fc <- ifelse(log2fc_filtered >= 0, 2^log2fc_filtered, 1 / (2^abs(log2fc_filtered)))
    
    # Convert fold change to percentage
    ora_fc <- (ora_fc - 1) * 100
    
    # Apply scaling to p.adjust values for categories only
    scale_size <- function(x, min_size = 0.5, max_size = 4) {
      scaled_x <- (x - min(x)) / (max(x) - min(x))
      return(scaled_x * (max_size - min_size) + min_size)
    }
    
    # Filter the ora_x object based on the gene count
    if (is_gsea) {
      gene_counts <- sapply(strsplit(ora_x@result$core_enrichment, "/"), length)
    } else {
      gene_counts <- sapply(strsplit(ora_x$geneID, "/"), length)
    }
    ora_x@result <- ora_x@result[gene_counts >= min_genes_per_set, ]
    
    # Filter based on p-value cutoff
    ora_x@result <- ora_x@result[ora_x@result$p.adjust < pval, ]
    
    # If selected_pathways is provided, filter to include only those pathways
    if (!is.null(selected_pathways)) {
      ora_x@result <- ora_x@result[ora_x@result$Description %in% selected_pathways, ]
    } else {
      # Limit to the top pathways
      if (nrow(ora_x@result) > num_path) {
        ora_x@result <- ora_x@result[1:num_path, ]
      }
    }
    
    # Prepare category sizes for scaling
    category_sizes <- ora_x@result$p.adjust
    category_sizes <- scale_size(category_sizes)
    
    # Determine color scheme based on fold change values
    if (all(ora_fc >= 0)) {
      low_color <- "white"
      mid_color <- "#FFEBEE"
      high_color <- "#B71C1C"
    } else if (any(ora_fc >= 0) && any(ora_fc < 0)) {
      low_color <- "#0D47A1"
      mid_color <- "white"
      high_color <- "#B71C1C"
    } else {
      low_color <- "#0D47A1"
      mid_color <- "#E3F2FD"
      high_color <- "white"
    }
    
    print("head(ora_x)")
    print(head(ora_x))
    
    # Generate the cnetplot
    p <- cnetplot(
      ora_x,
      foldChange = ora_fc,
      showCategory = nrow(ora_x@result),
      layout = "kk",
      circular = circ_param,
      colorEdge = TRUE,
      node_label = "all",
      cex_gene = 2,
      cex_label_gene = 0,
      cex_label_category = 0,
      color_category = 'grey30'
    )
    
    # Create the plot title
    plot_title <- paste(enrichment_type, "Enriched Gene Set Network")
    
    # Customize the plot
    p_ggraph <- p +
      scale_color_gradient2(low = low_color, mid = mid_color, high = high_color, name = "% Fold Change") +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(color = "black"),
        legend.title = element_text(color = "black"),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
      ) +
      ggtitle(plot_title)
    
    # Add labels using geom_text_repel for all nodes
    p_ggraph <- p_ggraph + 
      geom_text_repel(
        data = subset(p_ggraph$data, !is.na(name)),  # Ensure that we have valid names for labeling
        aes(x = x, y = y, label = name),
        size = 3,
        max.overlaps = 20
      )
    
    return(p_ggraph)
  }
  
  
  # Update the gsea1_kegg_x object with the filtered gene list
  gsea1_kegg_x <- update_core_enrichment(gsea1_kegg_2, filtered_gene_list_pathway)
  
  # Generate the plots
  
  # Function to extract gene lists from enrichResult or gseaResult
  extract_gene_list <- function(ora_x) {
    if (inherits(ora_x, "gseaResult")) {
      # For GSEA, extract core_enrichment
      gene_lists <- ora_x@result$core_enrichment
    } else {
      # For ORA, extract geneID
      gene_lists <- ora_x@result$geneID
    }
    # Split by "/" and unlist
    gene_list <- unique(unlist(strsplit(gene_lists, "/")))
    return(gene_list)
  }
  
  top_kegg <- return_top_n(gsea1_kegg_2,10,25,pval_cut)
  
  gsea_kegg_dot <- ggplot(top_kegg, aes(x = NES, y = reorder(Description, NES))) +
    geom_point(aes(color = p.adjust, size = GeneRatio)) +
    scale_color_gradient(low = "#2e2e2e", high = "#CCCCCC") +
    labs(x = "Normalized Enrichment Score (NES)",
         y = "Pathway",
         color = "Adjusted p-value",
         size = "Gene Ratio") +
    scale_x_continuous(expand = expansion(mult = c(0.15, 0.15)))+
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8),
          #remove minor break lines on x axis
          panel.grid.minor.x = element_blank())
  
  
  gsea1_kegg_x@result <- gsea1_kegg_x@result[!is.na(gsea1_kegg_x@result$core_enrichment) & gsea1_kegg_x@result$core_enrichment != "NA", ]
  
  # Extract gene lists for GSEA
  gsea1_kegg_x_gene_list <- extract_gene_list(gsea1_kegg_x)
  
  # Function to create a fold change vector
  create_fc_vector <- function(gene_list, log2fc_vector) {
    # Initialize a numeric vector with zeros
    result <- numeric(length(gene_list))
    names(result) <- gene_list
    
    # Assign fold change values where available
    matched_genes <- intersect(names(log2fc_vector), gene_list)
    result[matched_genes] <- log2fc_vector[matched_genes]
    
    return(result)
  }
  
  gsea_kegg_fc_vector <- create_fc_vector(gsea1_kegg_x_gene_list, log2fc_vector)
  gsea_kegg_pval_vector <- create_fc_vector(gsea1_kegg_x_gene_list, pval_vector)
  
  # drop  items == to 0 from gsea_kegg_fc_vector
  gsea_kegg_fc_vector <- gsea_kegg_fc_vector[gsea_kegg_fc_vector != 0]
  gsea_kegg_pval_vector <- gsea_kegg_pval_vector[gsea_kegg_pval_vector != 0]
  
  # 1. Plot for all pathways
  gsea_kegg_fig_all <- generate_cnetplot(
    ora_x = gsea1_kegg_x,
    log2fc_vector = gsea_kegg_fc_vector,
    pval_vector = gsea_kegg_pval_vector,
    num_path = 20,           # Adjust the number of pathways as needed
    ora_type = "all",        # Plot all enriched pathways
    min_genes_per_set = 1,   # Minimum number of genes per set
    circ_param = FALSE,      # Non-circular layout
    pval = 0.00001              # P-value cutoff
  )
  
  # 2. Plot for positively enriched pathways (circular)
  gsea_kegg_fig_pos <- generate_cnetplot_nocluster(
    ora_x = gsea1_kegg_x,
    log2fc_vector = gsea_kegg_fc_vector,
    pval_vector = gsea_kegg_pval_vector,
    num_path = 40,           # Adjust the number of pathways as needed
    ora_type = "pos",        # Plot positively enriched pathways (NES > 0)
    min_genes_per_set = 3,   # Minimum number of genes per set
    circ_param = FALSE,       # Circular layout
    pval = 0.0001               # P-value cutoff
  )
  
  ### Clustering, not currently working
  # generate_cnetplot_cluster(
  #   ora_x = gsea1_kegg_x,
  #   log2fc_vector = gsea_kegg_fc_vector,
  #   gsea_tree_data = gsea_tree_all$data, # Provide your gsea_tree data here
  #   num_path = 10,
  #   ora_type = "all",
  #   min_genes_per_set = 3,
  #   circ_param = FALSE,
  #   pval = 0.001
  # )
  
  
  # 3. Plot for negatively enriched pathways (circular)
  gsea_kegg_fig_neg <- generate_cnetplot(
    ora_x = gsea1_kegg_x,
    log2fc_vector = gsea_kegg_fc_vector,
    pval_vector = gsea_kegg_pval_vector,
    num_path = 50,           # Adjust the number of pathways as needed
    ora_type = "neg",        # Plot negatively enriched pathways (NES < 0)
    min_genes_per_set = 0,   # Minimum number of genes per set
    circ_param = FALSE,       # Circular layout
    pval = 0.0001              # P-value cutoff
  )
  
  my_pathways <- c("Drug metabolism - cytochrome P450",
                   "Retinol metabolism",
                   "Chemical carcinogenesis - DNA adducts",
                   "Bile secretion",
                   "Ascorbate and aldarate metabolism",
                   "Complement and coagulation cascades",
                   "Platelet activation",
                   "Natural killer cell mediated cytotoxicity",
                   #"Graft-versus-host disease",
                   "Autoimmune thyroid disease",
                   "Antigen processing and presentation",
                   "Th1 and Th2 cell differentiation",
                   "Neutrophil extracellular trap formation",
                   "Pathways in cancer","PI3K-Akt signaling pathway",
                   "B cell receptor signaling pathway",
                   "TNF signaling pathway","JAK-STAT signaling pathway",
                   #"Inflammatory bowel disease",
                   "Fc epsilon RI signaling pathway",
                   "Rap1 signaling pathway","VEGF signaling pathway",
                   "Focal adhesion","Cell adhesion molecules",
                   "T cell receptor signaling pathway",
                   "Toll-like receptor signaling pathway",
                   #"Shigellosis",
                   "Leukocyte transendothelial migration")

  # 1. Plot for all pathways
  gsea_kegg_fig_all_spec <- generate_cnetplot_spec(
    ora_x = gsea1_kegg_x,
    log2fc_vector = gsea_kegg_fc_vector,
    selected_pathways = my_pathways,
    num_path = 40,           # Adjust the number of pathways as needed
    ora_type = "all",        # Plot all enriched pathways
    min_genes_per_set = 2,   # Minimum number of genes per set
    circ_param = FALSE,      # Non-circular layout
    pval = 0.001              # P-value cutoff
  )
  
  # Display the plots
  # gsea_kegg_dot
  # gsea_kegg_fig_all
  # gsea_kegg_fig_pos
  # gsea_kegg_fig_neg
  
  
  # Save SVGs
  svg(paste("./Output/GSEA/",contrast_groups[2],contrast_groups[3],"_",ver,"_gsea_kegg_fig_pos.svg", sep = ""), width = 11, height = 8)
  gsea_kegg_fig_pos
  dev.off()
  svg(paste("./Output/GSEA/",contrast_groups[2],contrast_groups[3],"_",ver,"_gsea_kegg_fig_neg.svg", sep = ""), width = 11, height = 8)
  gsea_kegg_fig_neg
  dev.off()
  svg(paste("./Output/GSEA/",contrast_groups[2],contrast_groups[3],"_",ver,"_gsea_kegg_dot.svg", sep = ""), width = 7, height = 8)
  gsea_kegg_dot
  dev.off()
  
  # Save PNGs
  png(paste("./Output/GSEA/",contrast_groups[2],contrast_groups[3],"_",ver,"_gsea_kegg_fig_pos.png", sep = ""), width = 11, height = 8, units = 'in', res = 300)
  gsea_kegg_fig_pos
  dev.off()
  png(paste("./Output/GSEA/",contrast_groups[2],contrast_groups[3],"_",ver,"_gsea_kegg_fig_neg.png", sep = ""), width = 11, height = 8, units = 'in', res = 300)
  gsea_kegg_fig_neg
  dev.off()
  png(paste("./Output/GSEA/",contrast_groups[2],contrast_groups[3],"_",ver,"_gsea_kegg_dot.png", sep = ""), width = 7, height = 8, units = 'in', res = 300)
  gsea_kegg_dot
  dev.off()
  # save gsea_kegg_fig_all
  save_plot(gsea_kegg_fig_all, width=10,height=8,paste0(contrast_groups[2], contrast_groups[3], "_", ver, "_gsea_kegg_all_genes"))
  
  # save gsea_kegg_fig_all
  #save_plot(gsea_kegg_fig_all, width=15,height=10,paste0(contrast_groups[2], contrast_groups[3], "_", ver, "_gsea_kegg_all_genes"))
  
  save_results_to_csv(gsea1_kegg_2, paste0(contrast_groups[2], contrast_groups[3], "_", ver, "_gsea_kegg"))
  
}



### OTHER GSEA / KEGG PLOTS DEPRECTATE??? ----
if(FALSE){
  #install.packages("remotes")
  #remotes::install_github("GuangchuangYu/enrichplot")
  gsea1_pair_kegg <- enrichplot::pairwise_termsim(gsea1_kegg)
  gsea_fig_kegg <- enrichplot::treeplot(gsea1_pair_kegg,
                                        hilight = TRUE, 
                                        offset = 10,               # Adjust distance between bar and tree
                                        offset_tiplab = 0.5,        # Adjust distance between nodes and branches
                                        label_format = 20,          # Ensure labels wrap lines at 30 characters
                                        xlim = c(0, 10),
                                        fontsize=3,
                                        cluster.params = list (label_words_n = 6,
                                                               method = "average"
                                        )
  ) 
  
  # gsea_fig2 <- barplot(
  #   gsea1@result,
  #   x = "NES",
  #   color = "p.adjust",
  #   showCategory = 12) + 
  #   # increase font size
  #   theme(axis.text.x = element_text(size = 20),
  #         axis.text.y = element_text(size = 20),
  #         axis.title = element_text(size = 20),
  #         plot.title = element_text(size = 20),
  #         legend.text = element_text(size = 20),
  #         legend.title = element_text(size = 20)) +
  #   ggtitle("dotplot for GSEA") + geom_vline(xintercept = 0, linetype = "dashed", color = "black")
  # 
  
  # Load required libraries
  library(ggplot2)
  
  # Custom function to create logarithmic breaks
  log_breaks <- function(data, n = 5) {
    min_val <- min(data)
    max_val <- max(data)
    breaks <- 10^seq(log10(min_val), log10(max_val), length.out = n)
    return(breaks)
  }
  
  # Assuming your data frame is named 'gsea1@result'
  # Select and arrange the data
  plot_data <- gsea1_kegg@result
  plot_data <- plot_data[, c("Description", "NES", "p.adjust")]
  plot_data <- plot_data[order(plot_data$p.adjust, decreasing = TRUE), ]
  # select rows where p.adjust is < 0.1
  plot_data <- plot_data[plot_data$p.adjust < 0.005,]
  rownames(plot_data)
  # Determine breaks for p.adjust
  p_adjust_breaks <- log_breaks(plot_data$p.adjust, n = 4)
  
  # Determine x-axis limits and breaks
  x_min <- floor(min(plot_data$NES))
  x_max <- ceiling(max(plot_data$NES))
  x_breaks <- seq(x_min, x_max, by = 1)
  
  # Create the plot
  # Create the plot
  gsea_fig2 <- ggplot(plot_data, aes(x = NES, y = reorder(Description, NES), fill = p.adjust)) +
    geom_bar(stat = "identity") +
    scale_fill_gradient(low = "red", high = "blue", 
                        name = "p.adjust",
                        trans = "log10",
                        breaks = p_adjust_breaks,
                        labels = scales::scientific_format()(p_adjust_breaks)) +
    scale_x_continuous(limits = c(x_min, x_max), breaks = x_breaks) +
    labs(x = "NES", y = NULL) +
    theme_minimal() +
    theme(
      panel.grid.major.x = element_line(color = "grey90", size = 0.2),
      panel.grid.major.y = element_line(color = "grey90", size = 0.2),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15, color = "black"),
      axis.title.x = element_text(size = 15),
      plot.title = element_text(hjust = 0, size = 12, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 12),
      legend.key.size = unit(0.5, "cm"),
      plot.margin = margin(10, 10, 10, 10)
    ) +
    ggtitle("B") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.5)
  
  
  for (i in 1:length(all_genesets)) {
    for (j in 1:length(all_genesets)) {
      similarity_matrix[i, j] <- calc_jaccard(all_genesets[[i]], all_genesets[[j]])
    }
  }
  
  # Create network graph
  graph <- graph_from_adjacency_matrix(similarity_matrix, weighted = TRUE, mode = "undirected", diag = FALSE)
  
  # Set a threshold for edge visibility (e.g., 0.1 for 10% similarity)
  E(graph)$weight[E(graph)$weight < 0.1] <- 0
  graph <- delete_edges(graph, E(graph)[E(graph)$weight == 0])
  
  # Set node colors based on source
  V(graph)$color <- ifelse(V(graph)$name %in% names(gsea_genes), "red", "blue")
  
  # Plot network
  plot(graph, 
       vertex.size = 10,
       vertex.label = V(graph)$name,
       vertex.label.cex = 0.6,
       vertex.label.color = "black",
       edge.width = E(graph)$weight * 5,
       layout = layout_with_fr(graph),
       main = "Geneset Similarity Network")
  
  # Add legend for node colors
  legend("bottomright", 
         legend = c("GSEA", "ORA"), 
         fill = c("red", "blue"), 
         title = "Source")
  
  # Save all files to the output folder /GSEA
  png(paste("./Output/GSEA/",contrast_groups[2],contrast_groups[3],"_",ver,"_ora_pos.png", sep = ""), width = 800, height = 500)
  ora_pos_kegg_fig
  dev.off()
  png(paste("./Output/GSEA/",contrast_groups[2],contrast_groups[3],"_",ver,"_ora_neg.png", sep = ""), width = 800, height = 500)
  ora_neg_kegg_fig
  dev.off()
  png(paste("./Output/GSEA/",contrast_groups[2],contrast_groups[3],"_",ver,"_ora_pos_net.png", sep = ""), width = 800, height = 500)
  ora_pos_kegg_fig_net
  dev.off()
  png(paste("./Output/GSEA/",contrast_groups[2],contrast_groups[3],"_",ver,"_ora_neg_net.png", sep = ""), width = 800, height = 500)
  ora_neg_kegg_fig_net
  dev.off()
  png(paste("./Output/GSEA/",contrast_groups[2],contrast_groups[3],"_",ver,"_gsea_gobp.png", sep = ""), width = 800, height = 400)
  gsea_fig2
  dev.off()
  # Save SVGs
  svg(paste("./Output/GSEA/",contrast_groups[2],contrast_groups[3],"_",ver,"_ora_pos.svg", sep = ""), width = 800, height = 500)
  ora_pos_kegg_fig
  dev.off()
  svg(paste("./Output/GSEA/",contrast_groups[2],contrast_groups[3],"_",ver,"_ora_neg.svg", sep = ""), width = 800, height = 500)
  ora_neg_kegg_fig
  dev.off()
  svg(paste("./Output/GSEA/",contrast_groups[2],contrast_groups[3],"_",ver,"_ora_pos_net.svg", sep = ""), width = 800, height = 500)
  ora_pos_kegg_fig_net
  dev.off()
  svg(paste("./Output/GSEA/",contrast_groups[2],contrast_groups[3],"_",ver,"_ora_neg_net.svg", sep = ""), width = 800, height = 500)
  ora_neg_kegg_fig_net
  dev.off()
  svg(paste("./Output/GSEA/",contrast_groups[2],contrast_groups[3],"_",ver,"_gsea_gobp.svg", sep = ""), width = 800, height = 500)
  gsea_fig2
  dev.off()
  
  # Save ora_all_x@result to a file in Output/GSEA folder:
  write.table(ora_pos_kegg_x@result, file = paste("./Output/GSEA/",contrast_groups[2],contrast_groups[3],"_",ver,"_ora_pos_kegg_x.tsv", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(ora_neg_kegg_x@result, file = paste("./Output/GSEA/",contrast_groups[2],contrast_groups[3],"_",ver,"_ora_neg_kegg_x.tsv", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
}

### SAVE CROSS TABULATION CHARTS ----
if(output_xtabs == 1){
  # Convert meta columns to factors (except for condition and age)
  meta_factor <- meta
  #meta_factor[cols] <- lapply(meta_factor[cols],factor)
  summary(meta_factor)
  
  #Save PlotXTabs2 https://www.rdocumentation.org/packages/CGPfunctions/versions/0.6.3/topics/PlotXTabs2 to a grid
  
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
  
  png(paste("./Output/DESeq2/Xtabs/",contrast_groups[2],contrast_groups[3],"_",ver,".png", sep = ""),
      width=2600, height=1000, res = 144)
  ggarrange(plot_sex, plot_race, plot_batch, plot_chemo, plot_site, plot_age, ncol = 3, nrow = 2)
  #            labels = c("Sex", "Race", "Batch", "Primary Present","Primary Site", "Chemo 6 Weeks", "Age", "Other Mets (n/a"),
  #            vjust	= 0.5,
  dev.off()
  
  svg(paste("./Output/DESeq2/Xtabs/",contrast_groups[2],contrast_groups[3],"_",ver,".svg", sep = ""),
      width=2600, height=1000)
  ggarrange(plot_sex, plot_race, plot_batch, plot_chemo, plot_site, plot_age, ncol = 3, nrow = 2)
  dev.off()
  
}

### CORRELATION PLOTS ----
if(FALSE){
  # Correlation analysis between read counts
  correlation_matrix <- cor(norm_sig, method = "pearson")
  
  # Create annotation data frame
  annotation_data <- data.frame(
    #condition = meta$condition,
    primary_site = meta$primary_site,
    primary_present = meta$primary_present,
    peritoneal_mets = meta$peritoneal_mets
  )
  rownames(annotation_data) <- colnames(norm_sig)
  
  # Create annotation colors
  annotation_colors <- list(
    #condition = condition_table,
    primary_site = primary_site_table,
    primary_present = primary_present_table,
    peritoneal_mets = peritoneal_mets_table
  )
  
  # Create sample labels for x-axis and y-axis
  sample_labels <- colnames(norm_sig)
  
  png(paste("./Output/DESeq2/Correlation/correlation_heatmap_",contrast_groups[2],contrast_groups[3],"_",ver,".png", sep = ""), width = 900, height = 900)
  pheatmap(correlation_matrix, 
           color = colorRampPalette(c("#4575b4", "#ffffbf", "#d73027"))(100),
           main = "Correlation Heatmap",
           fontsize = 10,
           annotation_col = annotation_data,
           annotation_colors = annotation_colors,
           annotation_names_col = TRUE,
           annotation_legend = TRUE,
           clustering_distance_cols = "euclidean",
           clustering_distance_rows = "euclidean",
           show_rownames = FALSE,
           show_colnames = FALSE,
           labels_row = NULL,
           labels_col = NULL
  )
  dev.off()
  
  # Ensure that the sample names in norm_sig and meta are aligned
  common_samples <- intersect(colnames(norm_sig), rownames(meta))
  norm_sig_filtered <- norm_sig[, common_samples]
  meta_filtered <- meta[common_samples, , drop = FALSE]
  
  # Print the dimensions of norm_sig_filtered and meta_filtered
  cat("Dimensions of norm_sig_filtered:", dim(norm_sig_filtered), "\n")
  cat("Dimensions of meta_filtered:", dim(meta_filtered), "\n")
  
  # Statistical test to compare peak intensities between the two groups
  peak_intensities <- data.frame(Intensity = colMeans(norm_sig_filtered), 
                                 Group = meta_filtered$condition)
  wilcox_test <- wilcox.test(Intensity ~ Group, data = peak_intensities)
  cat("Wilcoxon Rank Sum Test:\n")
  print(wilcox_test)
  
}

### CEllular Deconvolution Plot ----
if(TRUE){
  # Load the source script and necessary libraries
  source("~/5hmC-Pathway-Analysis/R Code/repeated_glmnet_funcs_v1.R")
  library(tidyr)
  library(dplyr)
  library(ggplot2)
  library(scales)    # For percent_format()
  library(rstatix)   # For statistical tests and tidy output
  library(ggforce)   # For advanced faceting (if needed)
  library(broom)     # For tidying statistical test outputs
  
  # Define the file path
  file_path <- "~/5hmC-Pathway-Analysis/Raw Input/Working Inputs/tissue_of_origin_11_24_2024.csv"
  
  # Read the data without headers and specify NA strings
  data_raw <- read.csv(file_path, header = FALSE, stringsAsFactors = FALSE, na.strings = c("#N/A", "NA", ""))
  
  # Extract conditions and sample IDs from the first two rows
  conditions <- as.character(data_raw[1, 3:ncol(data_raw)])
  sample_ids <- as.character(data_raw[2, 3:ncol(data_raw)])
  
  # Create a dataframe mapping sample IDs to conditions
  sample_info <- data.frame(sample_id = sample_ids, condition = conditions, stringsAsFactors = FALSE)
  
  # Extract the main data starting from the third row
  data <- data_raw[3:nrow(data_raw), ]
  
  # Assign column names
  colnames(data) <- c("cell_type_system", "cell_type", sample_ids)
  
  # Reshape the data from wide to long format
  data_long <- data %>%
    pivot_longer(cols = -c(cell_type_system, cell_type), names_to = "sample_id", values_to = "percent")
  
  # Convert 'percent' to numeric, handling non-numeric entries
  data_long$percent <- as.numeric(data_long$percent)
  
  # Merge with sample_info to add the 'condition' column
  data_long <- data_long %>%
    left_join(sample_info, by = "sample_id")
  
  # Remove rows with NA in 'percent' or 'condition' if necessary
  data_long <- data_long %>%
    filter(!is.na(percent) & !is.na(condition))
  
  # Get the condition_table
  color_tables <- create_color_tables()
  condition_table <- color_tables$condition_table
  
  # Ensure that the 'condition' in data_long matches the names in condition_table
  data_long$condition <- factor(data_long$condition, levels = names(condition_table))
  
  # Sum the percent values for each cell_type_system per sample
  data_summed <- data_long %>%
    group_by(cell_type_system, sample_id, condition) %>%
    summarize(total_percent = sum(percent, na.rm = TRUE), .groups = "drop")
  
  # Calculate median total_percent for each cell_type_system and condition
  medians_summed <- data_summed %>%
    group_by(cell_type_system, condition) %>%
    summarize(median_total_percent = median(total_percent, na.rm = TRUE), .groups = "drop")
  
  # Adjust the median annotations to correct x positions
  dodge <- position_dodge(width = 0.75)
  
  # Create the first faceted boxplot for summed values by cell_type_system with adjusted y-axis limits
  p1 <- ggplot(data_summed, aes(x = condition, y = total_percent, fill = condition)) +
    geom_boxplot(outlier.shape = NA, position = dodge, coef = 1.5) +  # Disable outliers and set position dodge
    geom_text(
      data = medians_summed,
      aes(
        x = condition,
        y = median_total_percent,
        label = sprintf("%.1f%%", median_total_percent * 100),
        group = condition
      ),
      position = dodge,
      color = "white",
      size = 3.5,
      vjust = -0.5
    ) +
    theme_bw() +
    facet_wrap(~cell_type_system, scales = "free", ncol = 2) +  # Free scales for x and y axes
    labs(
      title = "Total Percent by Cell Type System (Faceted, Free X and Y Axes)",
      x = "Condition",
      y = "Total Percent"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 8)
    ) +
    scale_fill_manual(values = condition_table) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    coord_cartesian(clip = "off")  # Allow annotations outside plot area
  
  # For the second plot, facet by cell_type
  # Calculate median percent for each cell_type and condition
  medians_long <- data_long %>%
    group_by(cell_type, condition) %>%
    summarize(median_percent = median(percent, na.rm = TRUE), .groups = "drop")
  
  # Adjust the position dodge
  dodge <- position_dodge(width = 0.75)
  
  # Create the second plot, faceted by cell_type
  p2 <- ggplot(data_long, aes(x = condition, y = percent, fill = condition)) +
    geom_boxplot(outlier.shape = NA, position = dodge, coef = 1.5) +  # Disable outliers and set position dodge
    geom_text(
      data = medians_long,
      aes(
        x = condition,
        y = median_percent,
        label = sprintf("%.1f%%", median_percent * 100),
        group = condition
      ),
      position = dodge,
      color = "white",
      size = 3,
      vjust = -0.5
    ) +
    theme_bw() +
    facet_wrap(~cell_type, scales = "free", ncol = 3) +  # Facet by cell_type
    labs(
      title = "Percent by Condition (Faceted by Cell Type, Free X and Y Axes)",
      x = "Condition",
      y = "Percent"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      strip.text = element_text(size = 8)
    ) +
    scale_fill_manual(values = condition_table) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    coord_cartesian(clip = "off")  # Allow annotations outside plot area
  
  # ----------------------------
  # Updated Statistical Tests
  # ----------------------------
  
  # Create folder if it doesn't exist
  output_dir <- "./Output/DESeq2/cell_decon"
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Function to perform pairwise Wilcoxon tests safely
  perform_pairwise_wilcox <- function(data, value_col, group_col, id_col) {
    # Check for at least two conditions
    conditions_present <- unique(data[[group_col]])
    if (length(conditions_present) < 2) {
      message("Not enough conditions for ", unique(data[[id_col]]), ". Skipping.")
      return(NULL)
    }
    
    # Check for sufficient variability
    variability_check <- data %>%
      group_by_at(group_col) %>%
      summarize(n_unique = n_distinct(.data[[value_col]]), .groups = "drop") %>%
      filter(n_unique > 1)
    
    if (nrow(variability_check) < 2) {
      message("Insufficient variability in ", unique(data[[id_col]]), ". Skipping.")
      return(NULL)
    }
    
    # Perform the pairwise Wilcoxon test
    test <- tryCatch({
      pairwise.wilcox.test(
        x = data[[value_col]],
        g = data[[group_col]],
        p.adjust.method = "BH"
      )
    }, error = function(e) {
      message("Could not perform test for ", unique(data[[id_col]]), ": ", e$message)
      return(NULL)
    })
    
    if (is.null(test)) {
      return(NULL)
    }
    
    # Extract p-values into a data frame
    pvals <- as.data.frame(as.table(test$p.value))
    colnames(pvals) <- c("group1", "group2", "p")
    
    # Add id_col
    pvals[[id_col]] <- unique(data[[id_col]])
    
    # Add adjusted p-values and significance levels
    pvals <- pvals %>%
      mutate(
        p.adj = p.adjust(p, method = "BH"),
        p.adj.signif = case_when(
          p.adj <= 0.001 ~ "***",
          p.adj <= 0.01 ~ "**",
          p.adj <= 0.05 ~ "*",
          p.adj <= 0.1 ~ ".",
          TRUE ~ "ns"
        )
      )
    
    return(pvals)
  }
  
  # Perform tests for cell_type_system
  tests_system_pairwise <- data_summed %>%
    group_split(cell_type_system) %>%
    lapply(function(df) {
      perform_pairwise_wilcox(df, value_col = "total_percent", group_col = "condition", id_col = "cell_type_system")
    }) %>%
    bind_rows()
  
  # Perform tests for cell_type
  tests_cell_type_pairwise <- data_long %>%
    group_split(cell_type) %>%
    lapply(function(df) {
      perform_pairwise_wilcox(df, value_col = "percent", group_col = "condition", id_col = "cell_type")
    }) %>%
    bind_rows()
  
  # ----------------------------
  # Compute and Merge Average Values
  # ----------------------------
  
  # 1. Compute average total_percent per cell_type_system and condition
  avg_system <- data_summed %>%
    group_by(cell_type_system, condition) %>%
    summarize(avg_total_percent = mean(total_percent, na.rm = TRUE), .groups = "drop")
  
  # 2. Merge average values into tests_system_pairwise
  tests_system_pairwise <- tests_system_pairwise %>%
    # Merge for group1
    left_join(avg_system, by = c("cell_type_system", "group1" = "condition")) %>%
    rename(group1_avg = avg_total_percent) %>%
    # Merge for group2
    left_join(avg_system, by = c("cell_type_system", "group2" = "condition")) %>%
    rename(group2_avg = avg_total_percent)
  
  # 3. Compute average percent per cell_type and condition
  avg_cell_type <- data_long %>%
    group_by(cell_type, condition) %>%
    summarize(avg_percent = mean(percent, na.rm = TRUE), .groups = "drop")
  
  # 4. Merge average values into tests_cell_type_pairwise
  tests_cell_type_pairwise <- tests_cell_type_pairwise %>%
    # Merge for group1
    left_join(avg_cell_type, by = c("cell_type", "group1" = "condition")) %>%
    rename(group1_avg = avg_percent) %>%
    # Merge for group2
    left_join(avg_cell_type, by = c("cell_type", "group2" = "condition")) %>%
    rename(group2_avg = avg_percent)
  
  
  # Save the pairwise test results to CSV
  if (!is.null(tests_system_pairwise)) {
    write.csv(tests_system_pairwise, file.path(output_dir, "pairwise_tests_cell_type_system.csv"), row.names = FALSE)
  } else {
    message("No pairwise tests were performed for cell_type_system.")
  }
  
  # Save the pairwise test results to CSV
  if (!is.null(tests_cell_type_pairwise)) {
    write.csv(tests_cell_type_pairwise, file.path(output_dir, "pairwise_tests_cell_type.csv"), row.names = FALSE)
  } else {
    message("No pairwise tests were performed for cell_type.")
  }
  
  # -------------------------
  # Save Plots as PNG and SVG
  # -------------------------
  
  # Define plot file names
  plot1_png <- file.path(output_dir, "total_percent_by_cell_type_system.png")
  plot1_svg <- file.path(output_dir, "total_percent_by_cell_type_system.svg")
  plot2_png <- file.path(output_dir, "percent_by_condition_cell_type.png")
  plot2_svg <- file.path(output_dir, "percent_by_condition_cell_type.svg")
  
  # Save the first plot (p1) as PNG
  ggsave(filename = plot1_png, plot = p1, width = 10, height = 8, dpi = 300)
  
  # Save the first plot (p1) as SVG
  ggsave(filename = plot1_svg, plot = p1, width = 10, height = 8)
  
  # Save the second plot (p2) as PNG
  ggsave(filename = plot2_png, plot = p2, width = 12, height = 16, dpi = 300)
  
  # Save the second plot (p2) as SVG
  ggsave(filename = plot2_svg, plot = p2, width = 12, height =16)
  
  # Optional: Print a message indicating that plots have been saved
  message("Plots and pairwise statistical test results have been saved in the directory: ", output_dir)
  
}

### TCGA STATSTICAL TEST: ----
# Load libraries
library(readr)     # for reading TSV
library(ggplot2)   # for plotting
library(ggsignif)  # for adding significance brackets easily

# 1. Read data (adjust column names if different)
df_test <- read_tsv(
  "~/5hmC-Pathway-Analysis/Raw Input/Working Inputs/test_genes_tcga.txt"
)
df_rand <- read_tsv(
  "~/5hmC-Pathway-Analysis/Raw Input/Working Inputs/rand_genes_tcga.txt"
)

# Extract the p-value column (assuming it’s called "p-value")
test_pvals <- df_test[["p-value"]]
rand_pvals <- df_rand[["p-value"]]

# 2. Combine into one data frame for plotting
df_plot <- data.frame(
  pval  = c(test_pvals, rand_pvals),
  group = c(
    rep("Test Genes", length(test_pvals)),
    rep("Random Genes", length(rand_pvals))
  )
)

# 3. Box plot + jitter + statistical annotation
#    We'll use geom_signif from ggsignif to show the test result
ggplot(df_plot, aes(x = group, y = pval)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  geom_signif(
    comparisons = list(c("Test Genes", "Random Genes")),
    test = "wilcox.test",                       # Wilcoxon rank-sum test
    test.args = list(alternative = "less"),     # one-sided test
    map_signif_level = TRUE                     # show significance markers
  ) +
  theme_minimal() +
  labs(
    title = "Comparison of p-value Distributions",
    x     = NULL,
    y     = "p-value"
  )
