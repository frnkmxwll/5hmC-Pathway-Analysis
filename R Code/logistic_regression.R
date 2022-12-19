#Purpose of this script is to visualize differential gene expression heatmaps,
#PCA plot and other charts between two comparison groups.

### INSTALL LIBRARIES
# Setup, uncomment follow
# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("DESeq2")
# install.packages("augment")
# BiocManager::install("ggiraphExtra")
#devtools::install_github("cardiomoon/ggiraphExtra")

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
#library(M3C)
library(ashr)
library(CGPfunctions) #needed for plotxtabs2 
library(ggpubr) # needed for whisker plot charting
library(aod) #for logistic regression
library(broom)
library(leaps) #for logistic regression optimizaton
library(pROC)
library(boot)
#library(ggfortify) #used for pca clustering
library(cluster) #used for pca clustering
#library(ggiraph)
#library(ggiraphExtra)
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

### Raw counts
#counts_name_training <- "./Output/Randomization/1o3_pdONLYpos_8o11_metNEGtumNEG_DESeq2_v3/1o3_pdONLYpos_8o11_metNEGtumNEG_training_rawcounts.csv"
#counts_name_validation <- "./Output/Randomization/1o3_pdONLYpos_8o11_metNEGtumNEG_DESeq2_v3/1o3_pdONLYpos_8o11_metNEGtumNEG_validation_rawcounts.csv"

### Norm counts
counts_name_training <- "./Output/Randomization/1o3_pdONLYpos_8o11_metNEGtumNEG_DESeq2_final/1o3_pdONLYpos_8o11_metNEGtumNEG_training_rawcounts.csv"
counts_name_validation <- "./Output/Randomization/1o3_pdONLYpos_8o11_metNEGtumNEG_DESeq2_final/1o3_pdONLYpos_8o11_metNEGtumNEG_validation_rawcounts.csv"

meta_name_training <- "./Output/Randomization/1o3_pdONLYpos_8o11_metNEGtumNEG_DESeq2_final/1o3_pdONLYpos_8o11_metNEGtumNEG_training_conditions.csv"
meta_name_validation <- "./Output/Randomization/1o3_pdONLYpos_8o11_metNEGtumNEG_DESeq2_final/1o3_pdONLYpos_8o11_metNEGtumNEG_validation_conditions.csv"
counts_data_training <- read.csv(counts_name_training,row.names = 1)
counts_data_validation <- read.csv(counts_name_validation,row.names = 1)

counts_name <- "./Output/Randomization/1o3_pdONLYpos_8o11_metNEGtumNEG_DESeq2_v3/1o3_pdONLYpos_8o11_metNEGtumNEG_training_rawcounts.csv"
meta_name <- "./Output/Randomization/1o3_pdONLYpos_8o11_metNEGtumNEG_DESeq2_v3/1o3_pdONLYpos_8o11_metNEGtumNEG_training_conditions.csv"
previous_results_name <- "./Output/DESeq2/logit/1o3_pdONLYpos8o11_metNEGtumNEG_training_vFINAL/individual_predictors.csv"

#read in data, define what counts & conditions files
counts_data <- read.csv(counts_name,row.names = 1)
meta <-  read.csv(meta_name,row.names = 1)
previous_logit_results <- read.csv(previous_results_name)
previous_results_cutoff <- 0.01
previous_logit_results <- previous_logit_results[previous_logit_results$`p.value` < previous_results_cutoff, ]
previous_logit_results <- previous_logit_results[complete.cases(previous_logit_results), ]

rownames(counts_data_training) <- str_replace(c(rownames(counts_data_training)),'\\-',"_")
rownames(counts_data_training) <- str_replace(c(rownames(counts_data_training)),'\\.',"_")
rownames(counts_data_validation) <- str_replace(c(rownames(counts_data_validation)),'\\-',"_")
rownames(counts_data_validation) <- str_replace(c(rownames(counts_data_validation)),'\\.',"_")

meta_training <-  read.csv(meta_name_training,row.names=1)
meta_validation <- read.csv(meta_name_validation,row.names=1)

#create conditions vector expected as input by glmnet (0,0,0,...,1,1,1,,...)
sample_count_training <- nrow(meta_training)
class1_name <- meta_training[2,1]
class2_name <- meta_training[nrow(meta_training),1]
class1_count_training <- sum(meta_training$condition == class1_name)
class2_count_training <- sum(meta_training$condition == class2_name)
conditions_vector=c(rep(0,class1_count_training),rep(1,class2_count_training))


###CONFIGURATION
#set working directory, select where you extracted folder
setwd("~/5hmC-Pathway-Analysis/")

#include to exclude specific sample
#counts_data <- subset(counts_data,select=-c(KT126))
#meta <- meta[!(row.names(meta) %in% c("KT126")),]

#define padj cutoff, you may need to run with several padj values until you have an appropriate number of significant results.
#used to select significant genes for results tables, PCA plots, heatmaps and UMAP plots.
data_set_var = 2 # 1=full dataset, 2 = training only, 3 = validation only
cutoff_type = 2 # 0=padj cutoff, default; 1=lfc & pvalue cutoff; 2 = no cutoff (genome wide); 3=gene list; 4=previously run sig list;
padj.cutoff = 0.1 # 0.1 default
pvalue.cutoff = 0.1
lfc.cutoff = 0.137504 #0.263034 ~ 20% change, 0.137504 ~ 10%, 0.584963 ~ 50% change, 

#Select version for all output files (e.g. 1, 2, 3, ...)

ver <- "training_vFINAL"
gene_number <- nrow(counts_data)

#Set desired outputs:
output_results_tables = 1
output_heatmap = 1
output_PCA = 1
output_config = 1

# define contrast groups
groups <- unique(meta[c("condition")])
contrast_groups <- c("condition",groups[1,1], groups[2,1])

if (file.exists(paste("./Output/DESeq2/logit/",contrast_groups[2],contrast_groups[3],"_",ver,"/",sep=""))) {
  cat("The folder already exists")
} else {
  dir.create(paste("./Output/DESeq2/logit/",contrast_groups[2],contrast_groups[3],"_",ver,"/",sep=""))
}

### Define datasets
if (data_set_var == 1){
  counts_data = counts_data
  meta = meta
}

if (data_set_var == 2){
  counts_data = counts_data_training
  meta = meta_training
}

if (data_set_var == 3){
  counts_data = counts_data_validation
  meta = meta_validation
}


###VALIDATION
#check columns are equal
all(colnames(counts_data) %in% rownames(meta))
all(colnames(counts_data) == rownames(meta))

###CREATE DESeq2 OBJECT
#load data into DESeq2 object dds
design_formula <- ~ condition #sex + age + race + batch + condition 

dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = meta, design = design_formula)

#load up size factors into dds object in order to normalize using median of ratios method
dds <- DESeq(dds)

#use counts 'normalized=true' function to pull out normalized counts
normalized_counts <- counts(dds,normalized=TRUE)

# extract results table, and shrink the log2 fold changes
res_table <- results(dds, contrast=contrast_groups, alpha = padj.cutoff)

#Consider replacing above line with shrunken values for fold change below:
#res_table_unshrunken <- results(dds, contrast=contrast_groups, alpha = padj.cutoff)
#res_table <- lfcShrink(dds, coef = 2, res=res_table_unshrunken)

#metadata(res_table)$filterThreshold

#plot(metadata(res_table)$filterNumRej, 
#     type="b", ylab="number of rejections",
#     xlab="quantiles of filter")
#lines(metadata(res_table)$lo.fit, col="red")
#abline(v=metadata(res_table)$filterTheta)
#dev.off()

#convert the results table into a tibble:
res_table_tb <- res_table %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

normalized_counts_tb <- normalized_counts %>% 
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

#Extract normalized expression for significant genes from the two comparison groups and set the gene column (1) to row names
norm_sig <- normalized_counts_tb[,c(1,2:nsamples)] %>% 
  filter(gene %in% sig$gene) %>% 
  data.frame() %>%
  column_to_rownames(var = "gene") 

### subset table to only keep significant genes using cutoff
if(cutoff_type == 0){
  sig <- res_table_tb %>% 
    filter(padj < padj.cutoff)
}

if(cutoff_type == 1){
  sig <- res_table_tb %>% 
    filter(abs(log2FoldChange) > lfc.cutoff, pvalue < pvalue.cutoff)
}

###Genome wide regression
if(cutoff_type == 2){
  sig <- res_table_tb
  counts_data_training_subset = counts_data_training
  counts_data_validation_subset = counts_data_validation
}

### Specific genes only
if(cutoff_type == 3){
  ### genes in logit with pvalue <0.05
  #counts_data_training_subset = counts_data_training[c("KRT34","SCGB1D4","FERD3L","HEXIM1","CLEC4M","OR5AC2","OR1L4","KRT19","KIF4B","IL10","OR1D2","TMEM236","CD1C","WNT8A","RPL29","APOF","HLA_DMB","LAIR2","METTL7B","PGLYRP3","KRTAP19_6","RPE65","ZNF816","CTD_2287O16_3","AL359195_1","CTD_2545G14_7","SBDS","CD300E","FKBP9","FCER1G","PCYT1A","PSAP","HCST","CTSS","IQCF5","TAAR5","ZNF572","SLC35B1","LGALS13","CD300LB","KRTAP21_2","PCDHB4","PRPS1L1","CLPS","S100A11","HIST1H4D","CAMP","MARCH1","KRT31","HRH4","CWF19L1","RP11_793H13_10","MS4A6E","FKBP1B","REG3A","OR7G3","PAX3","COMTD1","C21orf62","LRRC8C","FCRL5","HES3","ZNF134","ADGRE3","NLRC4","CLEC4G","PILRA","RP4_576H24_4","RP11_444E17_6","IL23A","NAIP","CEBPB","LITAF","LILRA4","SCGB1D2","LTB4R","SPRR2D","FOS","OTOP1","GJB7","RPL3","CTC_435M10_3","KLK3","HPD","SLAMF9","SLC9A4","STX10","LA16c_431H6_6","ANGPT4","TNFSF14","SPCS1","ARHGAP30","CXCR1","ITK","C6orf226","ITGAM","MEFV","LCE2C","SGPP2","SEMG2","CYP11B2","ISG20L2","C1orf54","LAMTOR4","HRH2","SLC22A14","ARRDC5","SV2A","MS4A1","ZNF284","CRIP3","SLC46A2","OR7E24","DDIT3","MS4A13","WDR46"),]
  #counts_data_validation_subset = counts_data_validation[c("KRT34","SCGB1D4","FERD3L","HEXIM1","CLEC4M","OR5AC2","OR1L4","KRT19","KIF4B","IL10","OR1D2","TMEM236","CD1C","WNT8A","RPL29","APOF","HLA_DMB","LAIR2","METTL7B","PGLYRP3","KRTAP19_6","RPE65","ZNF816","CTD_2287O16_3","AL359195_1","CTD_2545G14_7","SBDS","CD300E","FKBP9","FCER1G","PCYT1A","PSAP","HCST","CTSS","IQCF5","TAAR5","ZNF572","SLC35B1","LGALS13","CD300LB","KRTAP21_2","PCDHB4","PRPS1L1","CLPS","S100A11","HIST1H4D","CAMP","MARCH1","KRT31","HRH4","CWF19L1","RP11_793H13_10","MS4A6E","FKBP1B","REG3A","OR7G3","PAX3","COMTD1","C21orf62","LRRC8C","FCRL5","HES3","ZNF134","ADGRE3","NLRC4","CLEC4G","PILRA","RP4_576H24_4","RP11_444E17_6","IL23A","NAIP","CEBPB","LITAF","LILRA4","SCGB1D2","LTB4R","SPRR2D","FOS","OTOP1","GJB7","RPL3","CTC_435M10_3","KLK3","HPD","SLAMF9","SLC9A4","STX10","LA16c_431H6_6","ANGPT4","TNFSF14","SPCS1","ARHGAP30","CXCR1","ITK","C6orf226","ITGAM","MEFV","LCE2C","SGPP2","SEMG2","CYP11B2","ISG20L2","C1orf54","LAMTOR4","HRH2","SLC22A14","ARRDC5","SV2A","MS4A1","ZNF284","CRIP3","SLC46A2","OR7E24","DDIT3","MS4A13","WDR46"),]
  ### genes in final model
  counts_data_training_subset = counts_data_training[c("ZNF134","RPL3"),]
  counts_data_validation_subset = counts_data_validation[c("ZNF134","RPL3"),]
  sig <- filter(res_table_tb,gene=="ZNF134"|gene=="RPL3")
} 

### Specific genes only from previous genome wide analysis
if(cutoff_type == 4){
  ### genes in logit with pvalue <0.05
  #counts_data_training_subset = counts_data_training[c("KRT34","SCGB1D4","FERD3L","HEXIM1","CLEC4M","OR5AC2","OR1L4","KRT19","KIF4B","IL10","OR1D2","TMEM236","CD1C","WNT8A","RPL29","APOF","HLA_DMB","LAIR2","METTL7B","PGLYRP3","KRTAP19_6","RPE65","ZNF816","CTD_2287O16_3","AL359195_1","CTD_2545G14_7","SBDS","CD300E","FKBP9","FCER1G","PCYT1A","PSAP","HCST","CTSS","IQCF5","TAAR5","ZNF572","SLC35B1","LGALS13","CD300LB","KRTAP21_2","PCDHB4","PRPS1L1","CLPS","S100A11","HIST1H4D","CAMP","MARCH1","KRT31","HRH4","CWF19L1","RP11_793H13_10","MS4A6E","FKBP1B","REG3A","OR7G3","PAX3","COMTD1","C21orf62","LRRC8C","FCRL5","HES3","ZNF134","ADGRE3","NLRC4","CLEC4G","PILRA","RP4_576H24_4","RP11_444E17_6","IL23A","NAIP","CEBPB","LITAF","LILRA4","SCGB1D2","LTB4R","SPRR2D","FOS","OTOP1","GJB7","RPL3","CTC_435M10_3","KLK3","HPD","SLAMF9","SLC9A4","STX10","LA16c_431H6_6","ANGPT4","TNFSF14","SPCS1","ARHGAP30","CXCR1","ITK","C6orf226","ITGAM","MEFV","LCE2C","SGPP2","SEMG2","CYP11B2","ISG20L2","C1orf54","LAMTOR4","HRH2","SLC22A14","ARRDC5","SV2A","MS4A1","ZNF284","CRIP3","SLC46A2","OR7E24","DDIT3","MS4A13","WDR46"),]
  #counts_data_validation_subset = counts_data_validation[c("KRT34","SCGB1D4","FERD3L","HEXIM1","CLEC4M","OR5AC2","OR1L4","KRT19","KIF4B","IL10","OR1D2","TMEM236","CD1C","WNT8A","RPL29","APOF","HLA_DMB","LAIR2","METTL7B","PGLYRP3","KRTAP19_6","RPE65","ZNF816","CTD_2287O16_3","AL359195_1","CTD_2545G14_7","SBDS","CD300E","FKBP9","FCER1G","PCYT1A","PSAP","HCST","CTSS","IQCF5","TAAR5","ZNF572","SLC35B1","LGALS13","CD300LB","KRTAP21_2","PCDHB4","PRPS1L1","CLPS","S100A11","HIST1H4D","CAMP","MARCH1","KRT31","HRH4","CWF19L1","RP11_793H13_10","MS4A6E","FKBP1B","REG3A","OR7G3","PAX3","COMTD1","C21orf62","LRRC8C","FCRL5","HES3","ZNF134","ADGRE3","NLRC4","CLEC4G","PILRA","RP4_576H24_4","RP11_444E17_6","IL23A","NAIP","CEBPB","LITAF","LILRA4","SCGB1D2","LTB4R","SPRR2D","FOS","OTOP1","GJB7","RPL3","CTC_435M10_3","KLK3","HPD","SLAMF9","SLC9A4","STX10","LA16c_431H6_6","ANGPT4","TNFSF14","SPCS1","ARHGAP30","CXCR1","ITK","C6orf226","ITGAM","MEFV","LCE2C","SGPP2","SEMG2","CYP11B2","ISG20L2","C1orf54","LAMTOR4","HRH2","SLC22A14","ARRDC5","SV2A","MS4A1","ZNF284","CRIP3","SLC46A2","OR7E24","DDIT3","MS4A13","WDR46"),]
  ### genes in final model
  counts_data_training_subset = counts_data_training[previous_logit_results$term,]
  counts_data_training_subset = counts_data_training_subset[rowSums(is.na(counts_data_training_subset)) != ncol(counts_data_training_subset),]

  counts_data_validation_subset = counts_data_validation[previous_logit_results$term,]
  counts_data_validation_subset = counts_data_validation_subset[rowSums(is.na(counts_data_validation_subset)) != ncol(counts_data_validation_subset),]
  
  sig <- filter(res_table_tb,gene %in% previous_logit_results$gene)
} 

#determine number of samples and significant genes
nsamples <- ncol(normalized_counts_tb)
nsig <- nrow(sig)

###LOGISTIC REGRESSION
#combine metadata and counts data

#retain only significant genes
meta_factor <- meta
meta_factor_training <- meta_training
meta_factor_validation <- meta_validation

sig_gene_names <- sig$gene
meta_sig_colnames = c(colnames(meta_factor),sig_gene_names)
predictors = merge(meta_factor,t(norm_sig),by.x=0,by.y=0,all.x=TRUE, all.y=TRUE)
predictors_training = merge(meta_factor_training,t(counts_data_training_subset),by.x=0,by.y=0,all.x=TRUE, all.y=TRUE)
predictors_validation = merge(meta_factor_validation,t(counts_data_validation_subset),by.x=0,by.y=0,all.x=TRUE, all.y=TRUE)

#predictors = predictors[,c(colnames(meta_factor),counts_data_training_subset)]
#predictors_training = predictors_training[,c(colnames(meta_factor),counts_data_training_subset)]
#predictors_validation = predictors_validation[,c(colnames(meta_factor),counts_data_training_subset)]

#convert meta columns to factors
predictors[colnames(meta_factor)] <- lapply(predictors[colnames(meta_factor)],factor)
predictors_training[colnames(meta_factor)] <- lapply(predictors_training[colnames(meta_factor)],factor)
predictors_validation[colnames(meta_factor)] <- lapply(predictors_validation[colnames(meta_factor)],factor)

#replace outcome with 1=PM pos and 0 =PM neg (confirm that contrfast groups match)
predictors$condition = gsub(contrast_groups[2],'1',predictors$condition)
predictors$condition = gsub(contrast_groups[3],'0',predictors$condition)

predictors_training$condition = gsub(contrast_groups[2],'1',predictors_training$condition)
predictors_training$condition = gsub(contrast_groups[3],'0',predictors_training$condition)

predictors_validation$condition = gsub(contrast_groups[2],'1',predictors_validation$condition)
predictors_validation$condition = gsub(contrast_groups[3],'0',predictors_validation$condition)

# remove "-" from gene names which cause problems in gene formula
colnames(predictors) <- gsub ("-","_",colnames(predictors))

#predictors[,c(sig_gene_names)] <- lapply(predictors[,c(sig_gene_names)],as.numeric)
predictors$condition <- as.numeric(predictors$condition)
predictors$age <- as.numeric(predictors$age)
predictors <- subset(predictors,select=-c(peritoneal_mets,other_mets,primary_present))

predictors_training$condition <- as.numeric(predictors_training$condition)
predictors_training$age <- as.numeric(predictors_training$age)
predictors_training <- subset(predictors_training,select=-c(peritoneal_mets,other_mets,primary_present))

predictors_validation$condition <- as.numeric(predictors_validation$condition)
predictors_validation$age <- as.numeric(predictors_validation$age)
predictors_validation <- subset(predictors_validation,select=-c(peritoneal_mets,other_mets,primary_present))

predictor_count <- ncol(predictors)

# Plot all predictors individually
for (single_predictor_index in 2:predictor_count) {
#  logit_variables = paste(colnames(predictors[-1]), collapse=" + ")
  logit_variables = colnames(predictors[single_predictor_index])
  logit_formula = as.formula(paste("condition ~ ",logit_variables, collapse=""))
  
  mylogit <- glm(logit_formula, data = predictors, family = "binomial")
  logistic_results <- summary(mylogit)
  #print(single_predictor_index)
  #print(summary(mylogit))
  if (!exists("tidy_glm")){
    tidy_glm = tidy(mylogit)
    } 
  else{
    tidy_glm <- rbind(tidy_glm,tidy(mylogit))
    #print(tidy_glm)
  }
  if(single_predictor_index == predictor_count){
    tidy_glm = tidy_glm[!(tidy_glm$term == "(Intercept)"),]
    write.csv(tidy_glm, paste("./Output/DESeq2/logit/",contrast_groups[2],contrast_groups[3],"_",ver,"/","individual_predictors.csv",sep=""))
    rm(tidy_glm)
  }
}
# Plot clinical confounding variables
logit_variables = paste("batch","chemo_6weeks", "primary_site", "age", "sex", sep = " + ")
print(logit_variables)
logit_formula = as.formula(paste("condition ~ ",logit_variables, collapse=""))

mylogit <- glm(logit_formula, data = predictors, family = "binomial")
logistic_results <- summary(mylogit)
print(summary(mylogit))
if (!exists("tidy_glm")){
  tidy_glm = tidy(mylogit)
} else{
  tidy_glm <- rbind(tidy_glm,tidy(mylogit))
  print(tidy_glm)
}
tidy_glm = tidy_glm[!(tidy_glm$term == "(Intercept)"),]
tidy_glm = tidy_glm[!duplicated(tidy_glm$term),]
write.csv(tidy_glm, paste("./Output/DESeq2/logit/",contrast_groups[2],contrast_groups[3],"_",ver,"/","confounders_batch_chemo_site_age_sex.csv",sep=""))
rm(tidy_glm)

test_prob = predict(mylogit, type = "response")

png(paste("./Output/DESeq2/logit/",contrast_groups[2],contrast_groups[3],"_",ver,"/","ROC_ClinicalOnly_",contrast_groups[2],contrast_groups[3],ver,".png",sep=""), width = 900, height = 900)
roc(predictors$condition ~ test_prob, plot = TRUE, print.auc = TRUE,main = paste("Clinical Variables Only"), cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
dev.off()


if(cutoff_type != 2){
  # Plot clinical confounding variables & 1 Gene at a time.
  for (single_predictor_index in 8:predictor_count) {
  #  logit_variables = paste("batch", "sex", "age", "chemo_6weeks", "primary_site", colnames(predictors[single_predictor_index]), sep = " + ") # all clinical variables
    logit_variables = paste("chemo_6weeks", colnames(predictors[single_predictor_index]), sep = " + ") # chemo only
    print(logit_variables)
    logit_formula = as.formula(paste("condition ~ ",logit_variables, collapse=""))
    
    mylogit <- glm(logit_formula, data = predictors, family = "binomial")
    logistic_results <- summary(mylogit)
    print(single_predictor_index)
    print(summary(mylogit))
    
    test_prob = predict(mylogit, type = "response")
    
    png(paste("./Output/DESeq2/logit/",contrast_groups[2],contrast_groups[3],"_",ver,"/","ROC_fullmodel_",single_predictor_index,"_",contrast_groups[2],contrast_groups[3],ver,".png",sep=""), width = 900, height = 900)
    roc(predictors$condition ~ test_prob, plot = TRUE, print.auc = TRUE,main = paste("Clinical Variables +",colnames(predictors[single_predictor_index])), cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2,ci=TRUE)
    dev.off()
    
    if (!exists("tidy_glm")){
      tidy_glm = tidy(mylogit)
    } else{
      tidy_glm <- rbind(tidy_glm,tidy(mylogit))
      print(tidy_glm)
    }
    if(single_predictor_index == predictor_count){
      tidy_glm = tidy_glm[!(tidy_glm$term == "(Intercept)"),]
      tidy_glm = tidy_glm[!duplicated(tidy_glm$term),]
      write.csv(tidy_glm, paste("./Output/DESeq2/logit/",contrast_groups[2],contrast_groups[3],"_",ver,"/","individual_genes_with_confounders_v6.csv",sep=""))
      rm(tidy_glm)
    }
  }
}



# for details see https://homepages.uc.edu/~lis6/Teaching/ML19Spring/Lab/lab6_VarSel.html
# Optimize model for AIC
all_variables = paste("batch", "sex", "age", "chemo_6weeks", "primary_site", paste(colnames(predictors[c(8:46)]),collapse = " + "),sep = " + ")

#bacward
#p05_variables = paste("batch", "sex", "age", "chemo_6weeks", "primary_site", "KRT34", "FERD3L", "HEXIM1","SCGB1D4","HLA_DMB","TMEM236","IL10","PSAP","HCST","RPL3","HIST1H4D","KRTAP21-2","IL10RA","ITK",sep = " + ")

#forward
p05_variables = paste("batch", "sex", "age", "chemo_6weeks", "primary_site","WNT8A","KIF4B","ZNF572","OR1L4","RPL41","KRT34","IFNL2",sep = " + ")
clinical_variables = paste("batch", "sex", "age", "chemo_6weeks", "primary_site",sep = " + ")
print(p05_variables)
print(clinical_variables)

#define full and null models
logit_formula = as.formula(paste("condition ~ ",logit_variables, collapse=""))
p05_formula= as.formula(paste("condition ~ ",p05_variables, collapse=""))
clinical_formula = as.formula(paste("condition ~ ",clinical_variables, collapse=""))

#Step forward
nullmodel<- lm(clinical_formula, data=predictors)
fullmodel<- lm(p05_formula, data=predictors)
model.step.b<- step(nullmodel,scope=list(lower=nullmodel,upper=fullmodel),direction='forward')



#Output final model
#logit_variables = paste("batch", "sex", "age", "chemo_6weeks", "primary_site", "SLC35B1",	"MANF",	"HIST1H4D",	"RPL3",	"ANGPT4",	"PLSCR4",	"FERD3L",	"BAP1",	"SCGB1D4",	"ITK",	"IL10RA",	"TMEM109",	"FAM214B",sep = " + ")
logit_variables = paste(chemo_6weeks,rownames(counts_data_training_subset)[1],rownames(counts_data_training_subset)[2],  sep = " + ")
#logit_variables = paste(colnames(predictors[,2:19108]), collapse = " + ")   
logit_formula = as.formula(paste("condition ~ ",logit_variables, collapse=""))
print(logit_formula)

logit_rpl3 <- glm(condition ~ RPL3, data = predictors_training, family = "binomial")
logit_znf134 <- glm(condition ~ ZNF134, data = predictors_training, family = "binomial")
logit_rpl3_znf134 <- glm(condition ~ ZNF134 + RPL3, data = predictors_training, family = "binomial")
logit_chemo_rpl3 <- glm(condition ~ chemo_6weeks + RPL3, data = predictors_training, family = "binomial")
logit_chemo_znf134 <- glm(condition ~ chemo_6weeks + ZNF134, data = predictors_training, family = "binomial")
logit_full <- glm(condition ~ chemo_6weeks + ZNF134 + RPL3, data = predictors_training, family = "binomial")

#logistic_results <- summary(mylogit)
#print(summary(mylogit))

#test_prob = predict(mylogit, type = "response")
pred_rpl3_val = predict(logit_rpl3, predictors_validation, type = "response")
pred_znf134_val = predict(logit_znf134, predictors_validation, type = "response")
pred_rpl3_znf134_val = predict(logit_rpl3_znf134, predictors_validation, type = "response")
pred_chemo_rpl3_val = predict(logit_chemo_rpl3, predictors_validation, type = "response")
pred_chemo_znf134_val = predict(logit_chemo_znf134, predictors_validation, type = "response")
pred_full_val = predict(logit_full, predictors_validation, type = "response")

pred_rpl3_tra = predict(logit_rpl3, predictors_training, type = "response")
pred_znf134_tra = predict(logit_znf134, predictors_training, type = "response")
pred_rpl3_znf134_tra = predict(logit_rpl3_znf134, predictors_training, type = "response")
pred_chemo_rpl3_tra = predict(logit_chemo_rpl3, predictors_training, type = "response")
pred_chemo_znf134_tra = predict(logit_chemo_znf134, predictors_training, type = "response")
pred_full_tra = predict(logit_full, predictors_training, type = "response")

png(paste("./Output/DESeq2/logit/",contrast_groups[2],contrast_groups[3],"_",ver,"/","ROC_fullmodel_valid_",contrast_groups[2],contrast_groups[3],ver,".png",sep=""), width = 900, height = 900)
  roc(predictors_validation$condition ~ test_prob, plot = TRUE, print.auc = TRUE,main = paste("full model against validation set"), cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, ci=TRUE)
dev.off()

#test_prob = predict(mylogit, type = "response")
test_prob = predict(logit_full, predictors_training, type = "response")

png(paste("./Output/DESeq2/logit/",contrast_groups[2],contrast_groups[3],"_",ver,"/","ROC_fullmodel_training_",contrast_groups[2],contrast_groups[3],ver,".png",sep=""), width = 900, height = 900)
roc(predictors_training$condition ~ test_prob, plot = TRUE, print.auc = TRUE,main = paste("full model against training set"), cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, ci=TRUE)
dev.off()

rocobj_rpl3_val <- roc(predictors_validation$condition ~ pred_rpl3_val)
rocobj_znf134_val <- roc(predictors_validation$condition ~ pred_znf134_val)
rocobj_rpl3_znf134_val <- roc(predictors_validation$condition ~ pred_rpl3_znf134_val)
rocobj_chemo_rpl3_val <- roc(predictors_validation$condition ~ pred_chemo_rpl3_val)
rocobj_chemo_znf134_val <- roc(predictors_validation$condition ~ pred_chemo_znf134_val)
rocobj_full_val <- roc(predictors_validation$condition ~ pred_full_val)

rocobj_rpl3_tra <- roc(predictors_training$condition ~ pred_rpl3_tra)
rocobj_znf134_tra <- roc(predictors_training$condition ~ pred_znf134_tra)
rocobj_rpl3_znf134_tra <- roc(predictors_training$condition ~ pred_rpl3_znf134_tra)
rocobj_chemo_rpl3_tra <- roc(predictors_training$condition ~ pred_chemo_rpl3_tra)
rocobj_chemo_znf134_tra <- roc(predictors_training$condition ~ pred_chemo_znf134_tra)
rocobj_full_tra <- roc(predictors_training$condition ~ pred_full_tra)

roc.val_list <- list(rocobj_rpl3_val,rocobj_znf134_val,rocobj_full_val)
roc.tra_list <- list(rocobj_rpl3_tra,rocobj_znf134_tra,rocobj_full_tra)

roc.tra_val_full_list <- list(rocobj_full_tra,rocobj_full_val)
roc.tra_val_rpl3_list <- list(rocobj_rpl3_tra,rocobj_rpl3_val)
roc.tra_val_znf134_list <- list(rocobj_znf134_tra,rocobj_znf134_val)

###source: https://stackoverflow.com/questions/66505014/how-to-add-auc-to-a-multiple-roc-graph-with-procs-ggroc
### Create data labels and chart for combined AUC Curves (VALIDATION)
val_labels <- data.frame()
i=1
for (item in roc.val_list){
  paste0("AUC = ",round(item$auc,digits=3)) -> val_labels[i,"auc"]
  paste0("95% CI: ",round(ci.auc(item)[1],digits=3),"-",round(ci.auc(item)[2],digits=3)) -> val_labels[i,"ci"]
  i=i+1
}

model_names <- c("RPL3","ZNF134", "Full Model")
val_labels["model_name"]<- paste(model_names)
val_labels["full_label"]<- paste(val_labels$model_name,val_labels$auc,val_labels$ci)                  

png(paste("./Output/DESeq2/logit/",contrast_groups[2],contrast_groups[3],"_",ver,"/",contrast_groups[2],contrast_groups[3],ver,"val_ROC",".png",sep=""), width = 900, height = 1200)
validation_roc <- ggroc(roc.val_list)
validation_roc + theme_minimal() + ggtitle("Validation ROC") + scale_color_discrete(labels=val_labels$full_label) + geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed") 
dev.off()

### Create data labels and chart for combined AUC Curves (TRAINING)
tra_labels <- data.frame()
i=1
for (item in roc.tra_list){
  paste0("AUC = ",round(item$auc,digits=3)) -> tra_labels[i,"auc"]
  paste0("95% CI: ",round(ci.auc(item)[1],digits=3),"-",round(ci.auc(item)[2],digits=3)) -> tra_labels[i,"ci"]
  i=i+1
}

tra_labels["model_name"]<- paste(model_names)
tra_labels["full_label"]<- paste(tra_labels$model_name,tra_labels$auc,tra_labels$ci)  

png(paste("./Output/DESeq2/logit/",contrast_groups[2],contrast_groups[3],"_",ver,"/",contrast_groups[2],contrast_groups[3],ver,"tra_ROC",".png",sep=""), width = 900, height = 1200)
training_roc <- ggroc(roc.tra_list)
training_roc + theme_minimal() + ggtitle("Training ROC") + scale_color_discrete(labels=tra_labels$full_label) + geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed") 
dev.off()

tra_val_full_roc <- ggroc(roc.tra_val_full_list)
tra_val_znf134_roc <- ggroc(roc.tra_val_znf134_list)
tra_val_rpl3_roc <- ggroc(roc.tra_val_rpl3_list)

#Combined validation and training plots
png(paste("./Output/DESeq2/logit/",contrast_groups[2],contrast_groups[3],"_",ver,"/",contrast_groups[2],contrast_groups[3],ver,"tra_val_fullROC",".png",sep=""), width = 900, height = 900)
tra_val_full_roc + theme_minimal() + ggtitle("Training and Validation Full Model ROC") + scale_color_discrete(labels=c(tra_labels$full_label[[3]],val_labels$full_label[[3]])) + geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed") 
dev.off()

png(paste("./Output/DESeq2/logit/",contrast_groups[2],contrast_groups[3],"_",ver,"/",contrast_groups[2],contrast_groups[3],ver,"tra_val_znf134ROC",".png",sep=""), width = 900, height = 900)
tra_val_znf134_roc + theme_minimal() + ggtitle("Training and Validation ZNF134 only ROC") + scale_color_discrete(labels=c(tra_labels$full_label[[3]],val_labels$full_label[[3]])) + geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed") 
dev.off()

png(paste("./Output/DESeq2/logit/",contrast_groups[2],contrast_groups[3],"_",ver,"/",contrast_groups[2],contrast_groups[3],ver,"tra_val_rpl3ROC",".png",sep=""), width = 900, height = 900)
tra_val_rpl3_roc + theme_minimal() + ggtitle("Training and Validation RPL3 only ROC") + scale_color_discrete(labels=c(tra_labels$full_label[[3]],val_labels$full_label[[3]])) + geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed") 
dev.off()


rm(logit_variables)

### Validation heatmap
#Save heatmap based on logistic regression results:
#log_counts_data<- rownames(counts_data_training_subset)
log_counts_data <- counts_data_validation_subset
rownames(log_counts_data) <- str_replace(c(rownames(log_counts_data)),'\\.',"_")
rownames(log_counts_data) <- str_replace(c(rownames(log_counts_data)),'\\-',"_")
#log_counts_data <- log_counts_data[c(rownames(counts_data_validation_subset)),c(rownames(meta_validation))]
heatmap_title <- paste(contrast_groups[2],"/",contrast_groups[3],"logistic regression sig genes p-value <0.025")

annotation <- meta_validation %>% 
  select(condition)

#transpose
#log_counts_data <- t(log_counts_data)
log_counts_data_mat = data.matrix(log_counts_data, rownames.force = NA)
log_counts_data_mat = t(scale(t(log_counts_data_mat)))
ha = HeatmapAnnotation(condition=annotation$condition,col = list(condition = c("1o3_pdONLYpos" = "#ff9289", "8o11_metNEGtumNEG" = "#00dae0")))

png(paste("./Output/DESeq2/logit/",contrast_groups[2],contrast_groups[3],"_",ver,"/",contrast_groups[2],contrast_groups[3],ver,"validation_heatmap_eucl",".png",sep=""), width = 900, height = 300)
#pheatmap(log_counts_data,
#         main = heatmap_title,
#         clustering_distance_cols = "euclidean",
#         cluster_rows = T,
#         cluster_cols = T,
#         show_rownames = T,
#         annotation = annotation,
#         #         annotation_row = annotation,  #uncomment for transpose
#         border_color = NA, 
#         fontsize = 10, 
#         scale = "row", 
#         #         scale = "column", #uncomment for transpose
#         fontsize_row = 10, 
#         height = 20,
Heatmap(log_counts_data_mat, 
        top_annotation = ha,
        #col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
        col = colorRamp2(c(-4, -2,-1,0,1, 2,4), c("#4676b5", "#82b1d3","#dbeff6","#fefebd","#fee395", "#fc9961","#d73027")),
        clustering_distance_rows = "pearson",
        clustering_distance_columns = "euclidean"
)
dev.off()

###Training heatmap
#Save heatmap based on logistic regression results:
#log_counts_data<- rownames(counts_data_training_subset)
log_counts_data <- counts_data_training_subset
rownames(log_counts_data) <- str_replace(c(rownames(log_counts_data)),'\\.',"_")
rownames(log_counts_data) <- str_replace(c(rownames(log_counts_data)),'\\-',"_")
#log_counts_data <- log_counts_data[c(rownames(counts_data_training_subset)),c(rownames(meta_training))]
heatmap_title <- paste(contrast_groups[2],"/",contrast_groups[3],"logistic regression sig genes p-value <0.025")

#transpose
#log_counts_data <- t(log_counts_data) #uncomment for transpose

annotation <- meta_training %>% 
  select(condition)

log_counts_data_mat = data.matrix(log_counts_data[rownames(log_counts_data) %in% previous_logit_results$term,], rownames.force = NA)
log_counts_data_mat = t(scale(t(log_counts_data_mat)))
ha = HeatmapAnnotation(condition=annotation$condition,col = list(condition = c("1o3_pdONLYpos" = "#ff9289", "8o11_metNEGtumNEG" = "#00dae0")))

png(paste("./Output/DESeq2/logit/",contrast_groups[2],contrast_groups[3],"_",ver,"/",contrast_groups[2],contrast_groups[3],ver,"training_heatmap_eucl",".png",sep=""), width = 900, height = 300)
#pheatmap(log_counts_data_mat,
#         main = heatmap_title,
##         clustering_distance_cols = "binary",
#         cluster_rows = T,
#         cluster_cols = T,
#         show_rownames = T,
#         top_annotation = ha,
##         annotation_row = annotation,  #uncomment for transpose
#         border_color = NA, 
#         fontsize = 10, 
#         scale = "row", 
##         scale = "column", #uncomment for transpose
#         fontsize_row = 10, 
#         height = 20
#)
Heatmap(log_counts_data_mat, 
        top_annotation = ha,
        #col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
        col = colorRamp2(c(-4, -2,-1,0,1, 2,4), c("#4676b5", "#82b1d3","#dbeff6","#fefebd","#fee395", "#fc9961","#d73027")),
        clustering_distance_rows = "pearson",
        clustering_distance_columns = "pearson",
)
dev.off()


###Whole dataset heatmap
#Save heatmap based on logistic regression results:
#log_counts_data<- rownames(counts_data_training_subset)
log_counts_data <- normalized_counts
rownames(log_counts_data) <- str_replace(c(rownames(log_counts_data)),'\\.',"_")
rownames(log_counts_data) <- str_replace(c(rownames(log_counts_data)),'\\-',"_")
log_counts_data <- log_counts_data[c(rownames(counts_data_training_subset)),]
heatmap_title <- paste(contrast_groups[2],"/",contrast_groups[3],"logistic regression sig genes p-value <0.01")

annotation <- meta %>% 
  select(condition)

#transpose
#log_counts_data <- t(log_counts_data)

png(paste("./Output/DESeq2/logit/",contrast_groups[2],contrast_groups[3],"_",ver,"/",contrast_groups[2],contrast_groups[3],ver,"wholeset100_def",".png",sep=""), width = 900, height = 1200)
pheatmap(log_counts_data,
         main = heatmap_title,
         #         clustering_distance_cols = "binary",
         cluster_rows = T,
         cluster_cols = T,
         show_rownames = T,
         annotation = annotation,
         #         annotation_row = annotation,  #uncomment for transpose
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         #         scale = "column", #uncomment for transpose 
         fontsize_row = 10, 
         height = 20,
)
dev.off()

###GENERATE PCA PLOT
#See details for below operation in lesson 3 of DGE workshop
rld <- vst(dds, blind=TRUE)
res_genes <- row.names(res_table)
sig_genes <-  sig$gene
#p < 0.05
#sig_genes <- c("KRT34","SCGB1D4","FERD3L","HEXIM1","CLEC4M","OR5AC2","OR1L4","KRT19","KIF4B","IL10","OR1D2","TMEM236","CD1C","WNT8A","RPL29","APOF","HLA_DMB","LAIR2","METTL7B","PGLYRP3","KRTAP19_6","RPE65","ZNF816","CTD_2287O16_3","AL359195_1","CTD_2545G14_7","SBDS","CD300E","FKBP9","FCER1G","PCYT1A","PSAP","HCST","CTSS","IQCF5","TAAR5","ZNF572","SLC35B1","LGALS13","CD300LB","KRTAP21_2","PCDHB4","PRPS1L1","CLPS","S100A11","HIST1H4D","CAMP","MARCH1","KRT31","HRH4","CWF19L1","RP11_793H13_10","MS4A6E","FKBP1B","REG3A","OR7G3","PAX3","COMTD1","C21orf62","LRRC8C","FCRL5","HES3","ZNF134","ADGRE3","NLRC4","CLEC4G","PILRA","RP4_576H24_4","RP11_444E17_6","IL23A","NAIP","CEBPB","LITAF","LILRA4","SCGB1D2","LTB4R","SPRR2D","FOS","OTOP1","GJB7","RPL3","CTC_435M10_3","KLK3","HPD","SLAMF9","SLC9A4","STX10","LA16c_431H6_6","ANGPT4","TNFSF14","SPCS1","ARHGAP30","CXCR1","ITK","C6orf226","ITGAM","MEFV","LCE2C","SGPP2","SEMG2","CYP11B2","ISG20L2","C1orf54","LAMTOR4","HRH2","SLC22A14","ARRDC5","SV2A","MS4A1","ZNF284","CRIP3","SLC46A2","OR7E24","DDIT3","MS4A13","WDR46")
#sig_genes <- c("KRT34","SCGB1D4","FERD3L","HEXIM1","CLEC4M","OR5AC2","OR1L4","KRT19","KIF4B","IL10","OR1D2","TMEM236","CD1C","WNT8A","RPL29","APOF","HLA_DMB","LAIR2","METTL7B","PGLYRP3","KRTAP19_6","RPE65","ZNF816","CTD_2287O16_3","AL359195_1","CTD_2545G14_7","SBDS","CD300E","FKBP9","FCER1G","PCYT1A","PSAP","HCST","CTSS","IQCF5","TAAR5","ZNF572","SLC35B1","LGALS13","CD300LB","KRTAP21_2","PCDHB4","PRPS1L1","CLPS","S100A11","HIST1H4D","CAMP","44621","KRT31","HRH4","CWF19L1","RP11_793H13_10","MS4A6E","FKBP1B","REG3A","OR7G3","PAX3","COMTD1","C21orf62","LRRC8C","FCRL5","HES3","ZNF134","ADGRE3","NLRC4","CLEC4G","PILRA","RP4_576H24_4","RP11_444E17_6","IL23A","NAIP","CEBPB","LITAF","LILRA4","SCGB1D2","LTB4R","SPRR2D","FOS","OTOP1","GJB7","RPL3","CTC_435M10_3","KLK3","HPD","SLAMF9","SLC9A4","STX10","LA16c_431H6_6","ANGPT4","TNFSF14","SPCS1","ARHGAP30","CXCR1","ITK","C6orf226","ITGAM","MEFV","LCE2C","SGPP2","SEMG2","CYP11B2","ISG20L2","C1orf54","LAMTOR4","HRH2","SLC22A14","ARRDC5","SV2A","MS4A1","ZNF284","CRIP3","SLC46A2","OR7E24","DDIT3","MS4A13","WDR46","NXT1","PILRB","S100A12","SLC7A7","IL10RA","CST11","ZNF729","RNF175","KRTAP4_5","MS4A4A","OR5B17","AURKB","GATS","OR6P1","GJB6","CASP14","PENK","C19orf45","TMEM154","KIAA0226L","PADI2","KRTAP26_1","AC107021_1","HK3","AJUBA","NPIPA3","OR10J3","ZNF320","ACP6","TIGD2","LY75_CD302","STARD10","CLEC4E","RRS1","SLC28A3","CAPG","CEACAM4","RP11_257K9_8","MS4A6A","TMEM38A","chemo_6weeksYes","MYO1F","HCRT","SARS2","SIGLECL1","SLC2A14","WDR77","IL1RN","GAS7","CALB2","MRFAP1L1","KRT32","CRCP","ZNF333","APOBEC3A","BCL2L15","CXCR4","EXD1","CEACAM8","RAB39A","ZNF397","SPINT4","MNX1","S100A9","TIFAB","S100A1","CD200R1L","C1orf194","HSH2D","HES1","PCDH17","PHKG1","MPEG1","LRRD1","FOXE3","DEF8","HIST1H2AG","OR2B2","OR7G2","SP140","TMEM50B","PDPR","CCR3","KLF6","CTSK","AC010547_9","G6PC","MSRB1","KIF15","LGALS9B","RP11_599B13_6","IL36G","DUS2","MTNR1A","EEF1A1","IFNL2","OR2T10","RHOG","ADORA2A","KRT38","NUPR1","MZF1","ZNF486","ID1","HCK","MARCO","APMAP","RFPL3","ASIP","GOLGA8O","APEX1","RP11_49K24_6","MRM1","CX3CR1","EIF5A","CHRNA6","NLRP12","TYROBP","LYPD6B","LYZ","CD53","TREM1","MS4A14","FOLR3","C19orf38","CCR8","GRPEL1","TCP10L2","NWD1","CNR2","SMPDL3B","TUBA4A","SLC31A2","ALAS1","TKT","OR10G6","TMEM180","DHRS9","FOXB2","SYTL3","SLC26A8","CD14","HTR3A","JUN","MANF","CD209","SUPT4H1","MS4A5","AK5","TMEM45B","MT4","IFNA6","RAB31","TVP23A","FCRL6","S100A8","SERPINB8","C5orf38","TAS2R5","TOP3A","CKAP4","CHGB","KRTAP5_9","TMEM130","ADRB2","SLC18A1","DNAJC5","TRAPPC5","POLE3","RPL18","OR6V1","HOXA10","SNX11","CCDC125","CLEC17A","FPR1","SMN2","ID4","KRT33A","NLRP2","HLX","GNS","LILRB2","C16orf70","ARHGAP9","ALOX15","CLDN8","COX7C","COX6B1","SPRR1B","GABRR2","CCDC140","PRR13","OCM","MRPL11","ACOXL","FAM72A","FGD3","NGB","C14orf144","MRFAP1","CD44","MPZL3","FMOD","TAGAP","DEFB107B","SERPINA4","LHX4","COMMD3_BMI1","PROKR2","OR2V1","OR51E1","IL32","NDUFS3","OR14K1","GLI1","GLTP","HTR3B","RP4_583P15_15","NLRP3","EEF2KMT","CDHR2","MED11","RDH16","GMFG","LEP","C1orf162","OR1E1","POLR2G","CPLX2","RAB11FIP1","COPS6","AC018755_1","LAIR1","KIF22","AP001024_2","CEP97","SLC39A11","CCL27","PGD","FAM105A","ARSG","FLT3","AK2","FCAR","ICAM3","MAPKAPK3","FCRLB","FGF6","IFI27L2","SLC2A5","CDA","HIST1H3B","SEC31B","E2F2","NEURL2","SELV","NSUN5","LY75","KRT5","TMEM191B","TSTD1","KRT36","CEACAM3","SLC22A15","FCGR3A","MGARP","NFKBID","TMEM52B","KAT8","MS4A18","LTA4H","RASGRP1","RAB37","FABP1","RXRB","FAM27E2","SIRPD","SIGMAR1","CSNK1A1L","ARIH2OS","SURF1","CD48","DNAJC5B","OR52I1","NCR3LG1","TNP1","DAGLB","SERPINA1","PUS3","SNX10","PLAC8","DIO2","MPPE1","SLC35F2","LGALS3","PRB2","ADPRH","OR6X1","LINC00694","C10orf12","ZAR1L","NCR3","CEACAM1","HOXC10","GLIPR2","RAB3D","LMNB1","ECH1","FAM214B","CD68","PDAP1","FCGBP","MAP10","AP4B1","TMEM81","SLC22A7","SEMA4A","KRTAP19_1","TXNL4B","TDRKH","FCGR2B","DRC3","NUDT19","A4GNT","ZNF813","VAV1","RPL41","RRP12","PLXNC1","CTD_2116N17_1","OLAH","NFKBIA","RPL26L1","POLR2J3","PKD2L1","CYSLTR2","ZSCAN26","BORCS7","AADACL3","ZNF788","PFKFB4","OR2Z1","AL669831_1","TRADD","DPAGT1","RUNX1T1","HIGD1B","MYOF","FOXM1","PATE3","RBM8A","PBX2","OR6C70","PLA2G2D","PPP1R17","HAVCR2","THOC3","BRI3BP","CFL1","TMPRSS4","PCGF1","OR10G8","S100P","TTC14","APBB1IP","CORO2A","CBR3","PTAFR","TAS2R39","ADAT1","RHOH","PRMT5","PGLYRP4","COLGALT2","AICDA","C12orf76","CTB_54O9_9","C1orf106","UHRF1","RP3_369A17_6","PSMB8","IZUMO3","GMIP","CATIP","CCDC74A")
#sig_genes <- c("RPL3","ZNF134")
rownames(rld) <- str_replace(c(rownames(rld)),'\\.',"_")
rownames(rld) <- str_replace(c(rownames(rld)),'\\-',"_")

#save PCA plot to png
#In the below replace sig_genes with res_genes if you want to perform PCA analysis on all genes rather than just on significant genes.
png(paste("./Output/DESeq2/logit/",contrast_groups[2],contrast_groups[3],"_",ver,"/","log_sig_PCA_",contrast_groups[2],contrast_groups[3],ver,".png",sep=""), width = 900, height = 1200)

plotPCA_labels <- plotPCA(
  rld[sig_genes,], intgroup = c("condition")
  #rownames(c(row.names(meta)))
) + coord_fixed(ratio=1) + stat_ellipse(level = 0.9)
 
#reverse order of labels to match heatmap labelling
plotPCA_labels$data$group<-factor(plotPCA_labels$data$group,levels = rev(levels(plotPCA_labels$data$group)))
plotPCA_labels$data$condition<-factor(plotPCA_labels$data$condition,levels = rev(levels(plotPCA_labels$data$condition)))

plotPCA_labels 
#+ geom_text(aes(label = name),position=position_nudge(y = 0.07),) + ggtitle(heatmap_title) 

dev.off()

### Additional logit statistics
### Source: https://thomaselove.github.io/432-notes/logistic-regression-and-the-resect-data.html
### Source: https://www.projectpro.io/recipes/select-best-cutoff-point-for-problem-roc-auc-curve-r#mcetoc_1g4q3v8dnk
confint(mylogit, level = 0.95)

mylogit_aug <- augment(logit_full, type.predict = "response")

#Plot each gene and where training case is
png(paste("./Output/DESeq2/logit/",contrast_groups[2],contrast_groups[3],"_",ver,"/","gene1prob_",contrast_groups[2],contrast_groups[3],ver,".png",sep=""), width = 900, height = 1200,pointsize=60)
  ggplot(mylogit_aug, aes(x = RPL3, y = condition)) +
    geom_jitter(height = 0.05) +
    geom_line(aes(x = RPL3, y = .fitted), 
              col = "blue") +
    labs(title = "Logistic Regression from Model mylogit")+
    theme(text = element_text(size = 30))
dev.off()
  
png(paste("./Output/DESeq2/logit/",contrast_groups[2],contrast_groups[3],"_",ver,"/","gene2prob_",contrast_groups[2],contrast_groups[3],ver,".png",sep=""), width = 900, height = 1200)
  ggplot(mylogit_aug, aes(x = ZNF134, y = condition)) +
    geom_jitter(height = 0.05) +
    geom_line(aes(x = ZNF134, y = .fitted), 
              col = "blue") +
    labs(title = "Logistic Regression from Model mylogit")+
    theme(text = element_text(size = 30))
dev.off()

png(paste("./Output/DESeq2/logit/",contrast_groups[2],contrast_groups[3],"_",ver,"/","gene3prob_",contrast_groups[2],contrast_groups[3],ver,".png",sep=""), width = 900, height = 1200)
  ggplot(mylogit_aug, aes(x = RPL3 + ZNF134, y = condition,color=chemo_6weeks)) +
    geom_jitter(height = 0.05) +
    geom_line(aes(x = RPL3 + ZNF134, y = .fitted)) +
    labs(title = "Logistic Regression from Model mylogit")+
    theme(text = element_text(size = 30))
dev.off()

png(paste("./Output/DESeq2/logit/",contrast_groups[2],contrast_groups[3],"_",ver,"/","full_model_",contrast_groups[2],contrast_groups[3],ver,".png",sep=""), width = 900, height = 1200)
ggplot(mylogit_aug, aes(x = c(chemo_6weeks,RPL3,ZNF134), y = condition)) +
  geom_jitter(height = 0.05) +
  geom_line(aes(x = c(chemo_6weeks,RPL3,ZNF134), y = .fitted), 
            col = "blue") +
  labs(title = "Logistic Regression from Model mylogit")+
  theme(text = element_text(size = 30))
dev.off()

#Plot to select cutoff
ROCR_pred_test <- prediction(test_prob,predictors_validation$condition)
ROCR_perf_test <- performance(ROCR_pred_test,'tpr','fpr')

png(paste("./Output/DESeq2/logit/",contrast_groups[2],contrast_groups[3],"_",ver,"/","cutoff_ROC_",contrast_groups[2],contrast_groups[3],ver,".png",sep=""), width = 900, height = 1200)
  plot(ROCR_perf_test,colorize=TRUE,print.cutoffs.at=seq(0.1,by=0.1), cex.lab=1, cex.axis=5)
dev.off()


mylogit_aug %$%
  confusionMatrix(
    data = factor(.fitted >= 0.5),
    reference = factor(condition == 1),
    positive = "TRUE"
  )

mylogit_aug %$%
  confusionMatrix(
    data = predictors_validation$condition,
    reference = ROCR_pred_test@labels
  )

plot(logit_formula, data=predictors_training, col="red4")
lines(logit_formula, newdat, col="green4", lwd=2)

