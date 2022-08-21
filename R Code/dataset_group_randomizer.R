# Purpose of this script is to take in a file that is ready for DESeq2 analysis
# and split it into a training and validation set.

### INSTALL LIBRARIES
# Setup, uncomment follow
# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("DESeq2")
# install.packages("dplyr")
# BiocManager::install("sva")
# install.packages("caTools")

#library(caTools)
#library(dplyr)


### CONFIGURATION
# set working directory
setwd("~/5hmC-Pathway-Analysis/")
counts_name <- "./Output/Raw Data Processing/1-2-3-4PMpos_8-9-11PMneg__combatseq/1-2-3-4PMpos_8-9-11PMneg_DESeq2_rawcounts.csv"
meta_name <- "./Output/Raw Data Processing/1-2-3-4PMpos_8-9-11PMneg__combatseq/1-2-3-4PMpos_8-9-11PMneg_DESeq2_conditions.csv"

file_version <- "_final"
random_seed=151

# settings
# recommend 70-30 or 80-20 split between training and validation.
training_fraction = 60
create_sig_training_counts = FALSE #SET TO TRUE ONLY AFTER HAVING RUN DESEQ2, WHEN RUNNING FOR FIRST TIME SET TO FALSE.
#When set to TRUE, outputs counts files with rows only for genes found to be significantly enriched.

#read in data, define what counts & conditions files
counts_data <- read.csv(counts_name,row.names = 1)
meta <-  read.csv(meta_name,row.names = 1)

#Used to subset genes from validation counts table based on sig results from DESeq2 on training counts
sig_counts_name <- "./Output/DESeq2/Results/significant_results_METnegPMpos_lfc_training_v1.txt"
sig_counts <- read.delim(sig_counts_name)

# Define two comparison groups
first_condition = distinct(meta,condition)[1,1]
second_condition = distinct(meta,condition)[2,1]

# Split data into two comparison groups
first_condition_meta = meta[meta$condition==first_condition,]
second_condition_meta = meta[meta$condition==second_condition,]

# Split Data into Training and Testing in R
# reference: https://www.programmingr.com/examples/neat-tricks/sample-r-function/randomly-split-data-in-r/
# UPDATE TO DEFINE SAMPLE SIZE FOR BOTH COMPARISON GROUP TABLES
first_sample_size = floor((training_fraction/100)*nrow(first_condition_meta))
second_sample_size = floor((training_fraction/100)*nrow(second_condition_meta))
set.seed(random_seed)

# randomly split data in r
# UPDATE TO SPLIT BOTH COMPARISON GROUP TABLES
first_picked = sample(seq_len(nrow(first_condition_meta)),size = first_sample_size)
second_picked = sample(seq_len(nrow(second_condition_meta)),size = second_sample_size)

first_training =first_condition_meta[first_picked,]
first_validation =first_condition_meta[-first_picked,]

second_training =second_condition_meta[second_picked,]
second_validation =second_condition_meta[-second_picked,]

#Combine both training meta tables together, and both validation meta tables together.
training_meta = rbind(first_training,second_training)
training_meta_ordered = training_meta[order(row.names(training_meta)), ]
training_meta_ordered = training_meta_ordered[order(training_meta_ordered$condition),]

validation_meta = rbind(first_validation,second_validation)
validation_meta_ordered = validation_meta[order(row.names(validation_meta)), ]
validation_meta_ordered = validation_meta_ordered[order(validation_meta_ordered$condition),]

#Use training and validation meta tables to select columns from count_data table
training_counts = counts_data[,colnames(counts_data) %in% rownames(training_meta_ordered)]
validation_counts = counts_data[,colnames(counts_data) %in% rownames(validation_meta_ordered)]


### OUTPUT RESUTLING FILES

if(create_sig_training_counts==FALSE){
  #Create root output folder if it doesn't exist.
  folder <- "./Output/Randomization/"
  if (file.exists(folder)) {
    cat("The folder already exists")
  } else {
    
    dir.create("./Output/Randomization/")
  }
  
  dir.create(paste("./Output/Randomization/",first_condition,"_",second_condition,"_","DESeq2_",file_version,sep=""))
  
  # output conditions files for DESeq2
  write.csv(
    training_meta_ordered,
    file = paste("./Output/Randomization/",first_condition,"_",second_condition,"_","DESeq2_",file_version,"/",first_condition,"_",second_condition,"_training_conditions.csv",sep = ""),
    row.names = TRUE,
    quote=FALSE)
  
  write.csv(
    validation_meta_ordered,
    file = paste("./Output/Randomization/",first_condition,"_",second_condition,"_","DESeq2_",file_version,"/",first_condition,"_",second_condition,"_validation_conditions.csv",sep = ""),
    row.names = TRUE,
    quote=FALSE)
  
  # output counts files for DESeq2
  write.csv(
    training_counts,
    file = paste("./Output/Randomization/",first_condition,"_",second_condition,"_","DESeq2_",file_version,"/",first_condition,"_",second_condition,"_training_rawcounts.csv",sep = ""),
    row.names = TRUE,
    quote=FALSE)
  
  write.csv(
    validation_counts,
    file = paste("./Output/Randomization/",first_condition,"_",second_condition,"_","DESeq2_",file_version,"/",first_condition,"_",second_condition,"_validation_rawcounts.csv",sep = ""),
    row.names = TRUE,
    quote=FALSE)
  
  #output config file
  config <- c(
    paste("input counts file name:", counts_name), 
    paste("input conditions file name:", meta_name), 
    paste("file version:",file_version),
    paste("Training fraction:", training_fraction,"%"),
    paste("Output folder:", "./Output/Randomization/",first_condition,"_",second_condition,"_","DESeq2_",file_version,sep=""),
    paste("Random seed:", random_seed)
    )
  
  write.table(
    config, 
    file=paste("./Output/Randomization/",first_condition,"_",second_condition,"_","DESeq2_",file_version,"/",first_condition,"_",second_condition,"_config",".txt", sep = ""), 
    sep="\t", 
    quote=F, 
    col.names=NA
    )
}

if(create_sig_training_counts==TRUE){
  validation_counts_sigonly = validation_counts[row.names(validation_counts) %in% sig_counts$gene,]
  training_counts_sigonly = training_counts[row.names(training_counts) %in% sig_counts$gene,]
  
  write.csv(
    validation_counts_sigonly,
    file = paste("./Output/Randomization/",first_condition,"_",second_condition,"_","DESeq2_",file_version,"/",first_condition,"_",second_condition,"_validation_sigonlycounts.csv",sep = ""),
    row.names = TRUE,
    quote=FALSE)
  
  write.csv(
    training_counts_sigonly,
    file = paste("./Output/Randomization/",first_condition,"_",second_condition,"_","DESeq2_",file_version,"/",first_condition,"_",second_condition,"_training_sigonlycounts.csv",sep = ""),
    row.names = TRUE,
    quote=FALSE)
}
