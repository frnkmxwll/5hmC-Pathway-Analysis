#PURPOSE OF THIS SCRIPT IS TO COMPARE SIMILARITY OF TWO ORDERED GENE LISTS
#Documentation from this link: https://bioconductor.org/packages/release/bioc/manuals/OrderedList/man/OrderedList.pdf
#Tutorial on OrderedList with examples: https://mirrors.nju.edu.cn/bioconductor/bioc/1.9/vignettes/OrderedList/inst/doc/tr_2006_01.pdf

### INSTALL LIBRARIES
# Setup, uncomment follow
# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("DESeq2")
# install.packages("dplyr")
# BiocManager::install("sva")
 library(OrderedList)
 library(dplyr)

 ###CONFIGURATION
 #set working directory, select where you extracted folder
 setwd("~/5hmC-Pathway-Analysis/")
 training_results_path <- "./Output/DESeq2/Results/all_results_METnegPMpos_lfc_training_v1.txt"
 training_sig_results_path <- "./Output/DESeq2/Results/significant_results_METnegPMpos_lfc_training_v1.txt"
 validation_results_path <- "./Output/DESeq2/Results/all_results_METnegPMpos_lfc_validation_v1.txt"
 
 #read in data, define what counts & conditions files
 training_results <- read.delim(training_results_path)
 validation_results <-  read.delim(validation_results_path)
 training_sig_results <- read.delim(training_sig_results_path)
 
 # eliminate rows with lfc=NA (this was due to all 0 reads likely.)
 training_results_complete <- training_results[complete.cases(training_results),]
 validation_results_complete <- validation_results[complete.cases(validation_results),]

 # eliminate genes only present in the training or validation sets.
 common <- intersect(training_results_complete$X, validation_results_complete$X)  
 training_results_complete <- training_results_complete[training_results_complete$X %in% common,] 
 validation_results_complete <- validation_results_complete[validation_results_complete$X %in% common,] 
 
 # order the results by lfc
 training_results_ordered <- training_results_complete[order(training_results_complete$log2FoldChange),]
 validation_results_ordered <- validation_results_complete[order(validation_results_complete$log2FoldChange),]
 
 #Need to renumber row names for new order 1- number of genes
 row.names(training_results_ordered) <- 1:nrow(training_results_ordered)
 row.names(validation_results_ordered) <- 1:nrow(validation_results_ordered) 

 #use the compareList function to compare the training and validation ranks
 training_list <- as.character(training_results_ordered[,1])
 validation_list <- as.character(validation_results_ordered[,1])
 training_validation_compared <- compareLists(training_list,validation_list, two.sided = TRUE)

 #note you can set the max rank for genes you want to plot and find overlap here.
  overlapping_genes <- getOverlap(training_validation_compared, max.rank=1000)

plot(overlapping_genes)
intersect(overlapping_genes$intersect, training_sig_results$gene)


training_sig_ranks <- training_results_ordered[which(training_results_ordered$X %in% training_sig_results$gene),]
validation_sig_ranks <- validation_results_ordered[which(validation_results_ordered$X %in% training_sig_results$gene),]

training_sig_ranks_pos <- subset(training_sig_ranks,training_sig_ranks$log2FoldChange > 0)
training_sig_ranks_neg <- subset(training_sig_ranks,training_sig_ranks$log2FoldChange < 0)

validation_sig_ranks_pos <- subset(validation_sig_ranks,validation_sig_ranks$log2FoldChange > 0)
validation_sig_ranks_neg <- subset(validation_sig_ranks,validation_sig_ranks$log2FoldChange < 0)

sum(as.numeric(rownames(training_sig_ranks_pos)))/nrow(validation_sig_ranks_pos)
sum(as.numeric(rownames(validation_sig_ranks_pos)))/nrow(validation_sig_ranks_pos)

sum(as.numeric(rownames(training_sig_ranks_neg)))/nrow(validation_sig_ranks_neg)
sum(as.numeric(rownames(validation_sig_ranks_neg)))/nrow(validation_sig_ranks_neg)
       