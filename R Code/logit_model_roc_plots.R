#Purpose of this script is generate ROC plots for specific genes
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

###CONFIGURATION
#set working directory, select where you extracted folder
setwd("~/5hmC-Pathway-Analysis/")
ver <- "model_genes_norm_counts_v1"

### Norm counts
counts_name_training <- "./Output/DESeq2/Results/normalized_counts_1o3_pdONLYpos8o11_metNEGtumNEG_pvalue0p1_lfc0p26_training_csarb.csv"
counts_name_validation <- "./Output/DESeq2/Results/normalized_counts_1o3_pdONLYpos8o11_metNEGtumNEG_pvalue0p005_lfc0p13_validation_csarb.csv"

meta_name_training <- "./Output/Randomization/1o3_pdONLYpos_8o11_metNEGtumNEG_DESeq2_v3/1o3_pdONLYpos_8o11_metNEGtumNEG_training_conditions.csv"
meta_name_validation <- "./Output/Randomization/1o3_pdONLYpos_8o11_metNEGtumNEG_DESeq2_v3/1o3_pdONLYpos_8o11_metNEGtumNEG_validation_conditions.csv"

previous_results_name <- "./Output/DESeq2/logit/1o3_pdONLYpos8o11_metNEGtumNEG_genome_wide_v1/individual_predictors.csv"

previous_logit_results <- read.csv(previous_results_name)
previous_results_cutoff <- 0.025
previous_logit_results <- previous_logit_results[previous_logit_results$`p.value` < previous_results_cutoff, ]

### Load up data
counts_data_training <- read.csv(counts_name_training,row.names = 1)
counts_data_validation <- read.csv(counts_name_validation,row.names = 1)

# Correct gene names for forbidden characters
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
conditions_vector_tra=c(rep(0,class1_count_training),rep(1,class2_count_training))

class1_count_validation <- sum(meta_validation$condition == class1_name)
class2_count_validation <- sum(meta_validation$condition == class2_name)
conditions_vector_val=c(rep(0,class1_count_validation),rep(1,class2_count_validation))

counts_data_training_subset = counts_data_training[c("ZNF134","RPL3"),]
counts_data_validation_subset = counts_data_validation[c("ZNF134","RPL3"),]

### Update meta files for logit
# define contrast groups
groups <- unique(meta_training[c("condition")])
contrast_groups <- c("condition",groups[1,1], groups[2,1])

meta_factor_training <- meta_training
meta_factor_validation <- meta_validation

predictors_training = merge(meta_factor_training,t(counts_data_training_subset),by.x=0,by.y=0,all.x=TRUE, all.y=TRUE)
predictors_validation = merge(meta_factor_validation,t(counts_data_validation_subset),by.x=0,by.y=0,all.x=TRUE, all.y=TRUE)

#convert meta columns to factors
predictors_training[colnames(meta_factor_training)] <- lapply(predictors_training[colnames(meta_factor_training)],factor)
predictors_validation[colnames(meta_factor_validation)] <- lapply(predictors_validation[colnames(meta_factor_validation)],factor)

#replace outcome with 1=PM pos and 0 =PM neg (confirm that contrfast groups match)
predictors_training$condition = gsub(contrast_groups[2],'1',predictors_training$condition)
predictors_training$condition = gsub(contrast_groups[3],'0',predictors_training$condition)

predictors_validation$condition = gsub(contrast_groups[2],'1',predictors_validation$condition)
predictors_validation$condition = gsub(contrast_groups[3],'0',predictors_validation$condition)

predictors_training$condition <- as.numeric(predictors_training$condition)
predictors_training$age <- as.numeric(predictors_training$age)
predictors_training <- subset(predictors_training,select=-c(peritoneal_mets,other_mets,primary_present))

predictors_validation$condition <- as.numeric(predictors_validation$condition)
predictors_validation$age <- as.numeric(predictors_validation$age)
predictors_validation <- subset(predictors_validation,select=-c(peritoneal_mets,other_mets,primary_present))

### Generate logit models
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

### Generate ROC objs for plots
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

roc.val_list <- list(rocobj_rpl3_val,rocobj_znf134_val,rocobj_full_val,rocobj_rpl3_znf134_val)
roc.tra_list <- list(rocobj_rpl3_tra,rocobj_znf134_tra,rocobj_full_tra,rocobj_rpl3_znf134_tra)

roc.tra_val_full_list <- list(rocobj_full_tra,rocobj_full_val)
roc.tra_val_rpl3_list <- list(rocobj_rpl3_tra,rocobj_rpl3_val)
roc.tra_val_znf134_list <- list(rocobj_znf134_tra,rocobj_znf134_val)
roc.tra_val_rpl3_znf134_list <- list(rocobj_rpl3_znf134_tra,rocobj_rpl3_znf134_val)

tra_val_full_roc <- ggroc(roc.tra_val_full_list)
tra_val_znf134_roc <- ggroc(roc.tra_val_znf134_list)
tra_val_rpl3_roc <- ggroc(roc.tra_val_rpl3_list)
tra_val_rpl3_znf134_roc <- ggroc(roc.tra_val_rpl3_znf134_list)

#Generate labels
val_labels <- data.frame()
i=1
for (item in roc.val_list){
  paste0("AUC = ",round(item$auc,digits=3)) -> val_labels[i,"auc"]
  paste0("95% CI: ",round(ci.auc(item)[1],digits=3),"-",round(ci.auc(item)[2],digits=3)) -> val_labels[i,"ci"]
  i=i+1
}

tra_labels <- data.frame()
i=1
for (item in roc.tra_list){
  paste0("AUC = ",round(item$auc,digits=3)) -> tra_labels[i,"auc"]
  paste0("95% CI: ",round(ci.auc(item)[1],digits=3),"-",round(ci.auc(item)[2],digits=3)) -> tra_labels[i,"ci"]
  i=i+1
}

model_names <- c("RPL3","ZNF134", "Full Model", "RPL3 + ZNF134")
val_labels["model_name"]<- paste(model_names)
val_labels["full_label"]<- paste("Validation set,",val_labels$model_name,val_labels$auc,val_labels$ci)                
tra_labels["model_name"]<- paste(model_names)
tra_labels["full_label"]<- paste("Training set,",tra_labels$model_name,tra_labels$auc,tra_labels$ci)  


#Combined validation and training plots
png(paste("./Output/DESeq2/logit/",contrast_groups[2],contrast_groups[3],"_",ver,"/",contrast_groups[2],contrast_groups[3],ver,"tra_val_fullROC",".png",sep=""), width = 900, height = 900)
tra_val_full_roc + theme_minimal() + ggtitle("Training and Validation Full Model ROC") + scale_color_discrete(labels=c(tra_labels$full_label[[3]],val_labels$full_label[[3]])) + geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed") 
dev.off()

png(paste("./Output/DESeq2/logit/",contrast_groups[2],contrast_groups[3],"_",ver,"/",contrast_groups[2],contrast_groups[3],ver,"tra_val_znf134ROC",".png",sep=""), width = 900, height = 900)
tra_val_znf134_roc + theme_minimal() + ggtitle("Training and Validation ZNF134 only ROC") + scale_color_discrete(labels=c(tra_labels$full_label[[2]],val_labels$full_label[[2]])) + geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed") 
dev.off()

png(paste("./Output/DESeq2/logit/",contrast_groups[2],contrast_groups[3],"_",ver,"/",contrast_groups[2],contrast_groups[3],ver,"tra_val_rpl3ROC",".png",sep=""), width = 900, height = 900)
tra_val_rpl3_roc + theme_minimal() + ggtitle("Training and Validation RPL3 only ROC") + scale_color_discrete(labels=c(tra_labels$full_label[[1]],val_labels$full_label[[1]])) + geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed") 
dev.off()

png(paste("./Output/DESeq2/logit/",contrast_groups[2],contrast_groups[3],"_",ver,"/",contrast_groups[2],contrast_groups[3],ver,"tra_val_rpl3_znf134ROC",".png",sep=""), width = 900, height = 900)
tra_val_rpl3_znf134_roc + theme_minimal() + ggtitle("Training and Validation RPL3 + ZNF134 ROC") + scale_color_discrete(labels=c(tra_labels$full_label[[4]],val_labels$full_label[[4]])) + geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed") 
dev.off()


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

png(paste("./Output/DESeq2/logit/",contrast_groups[2],contrast_groups[3],"_",ver,"/",contrast_groups[2],contrast_groups[3],ver,"validation_heatmap_euc",".png",sep=""), width = 900, height = 300)
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
        cluster_rows = FALSE,
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

log_counts_data_mat = data.matrix(log_counts_data, rownames.force = NA)
log_counts_data_mat = t(scale(t(log_counts_data_mat)))
ha = HeatmapAnnotation(condition=annotation$condition,col = list(condition = c("1o3_pdONLYpos" = "#ff9289", "8o11_metNEGtumNEG" = "#00dae0")))

png(paste("./Output/DESeq2/logit/",contrast_groups[2],contrast_groups[3],"_",ver,"/",contrast_groups[2],contrast_groups[3],ver,"training_heatmap_euc",".png",sep=""), width = 900, height = 300)
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
        cluster_rows = FALSE,
        clustering_distance_columns = "euclidean"
)
dev.off()


### Specific genes only from previous genome wide analysis
counts_data_training_subset = counts_data_training[previous_logit_results$term,]
counts_data_training_subset = counts_data_training_subset[rowSums(is.na(counts_data_training_subset)) != ncol(counts_data_training_subset),]

counts_data_validation_subset = counts_data_validation[previous_logit_results$term,]
counts_data_validation_subset = counts_data_validation_subset[rowSums(is.na(counts_data_validation_subset)) != ncol(counts_data_validation_subset),]

###Whole dataset heatmap
#Save heatmap based on logistic regression results:
#log_counts_data<- rownames(counts_data_training_subset)
log_counts_data <- counts_data_training_subset
rownames(log_counts_data) <- str_replace(c(rownames(log_counts_data)),'\\.',"_")
rownames(log_counts_data) <- str_replace(c(rownames(log_counts_data)),'\\-',"_")
log_counts_data <- log_counts_data[c(rownames(counts_data_training_subset)),]
heatmap_title <- paste(contrast_groups[2],"/",contrast_groups[3],"logistic regression sig genes p-value <0.025")

annotation <- meta_training %>% 
  select(condition)

log_counts_data_mat = data.matrix(log_counts_data, rownames.force = NA)
log_counts_data_mat = t(scale(t(log_counts_data_mat)))
ha = HeatmapAnnotation(condition=annotation$condition,col = list(condition = c("1o3_pdONLYpos" = "#ff9289", "8o11_metNEGtumNEG" = "#00dae0")))

png(paste("./Output/DESeq2/logit/",contrast_groups[2],contrast_groups[3],"_",ver,"/",contrast_groups[2],contrast_groups[3],ver,"wholeset100_rob",".png",sep=""), width = 900, height = 1200)

Heatmap(log_counts_data_mat, 
        top_annotation = ha,
        #col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
        col = colorRamp2(c(-4, -2,-1,0,1, 2,4), c("#4676b5", "#82b1d3","#dbeff6","#fefebd","#fee395", "#fc9961","#d73027")),
        clustering_distance_rows = "pearson",
        clustering_distance_columns = robust_dist,
)
dev.off()