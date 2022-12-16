# Purpose of this script is to summarize into an epigenetic score using a multivariable logistic regression model 
# and the elastic net regularization (or LASSO)

#tutorials found here:
#https://cran.r-project.org/web/packages/glmnet/vignettes/glmnet.pdf
#https://stats.stackexchange.com/questions/72251/an-example-lasso-regression-using-glmnet-for-binary-outcome
#https://stackoverflow.com/questions/18130338/plotting-an-roc-curve-in-glmnet
#

### INSTALL LIBRARIES
# Setup, uncomment follow
# install.packages("dplyr")
#install.packages("glmnet", repos = "https://cran.us.r-project.org")
#install.packages("svMisc")

library(glmnet)
library(dplyr)
library(ROCR)
library(svMisc) 
library (glmnetUtils)

###CONFIGURATION
#set working directory, select where you extracted folder
setwd("~/5hmC-Pathway-Analysis/")
counts_name_training <- "./Output/DESeq2/Results/normalized_counts_1o3_pdONLYpos8o11_metNEGtumNEG_pvalue0p01_lfc0p13_FINALtraining.csv"
counts_name_validation <- "./Output/DESeq2/Results/normalized_counts_1o3_pdONLYpos8o11_metNEGtumNEG_pvalue0p01_lfc0p13_FINALvalidation.csv"
meta_name_training <- "./Output/Randomization/1o3_pdONLYpos_8o11_metNEGtumNEG_DESeq2_final/1o3_pdONLYpos_8o11_metNEGtumNEG_training_conditions.csv"
meta_name_validation <- "./Output/Randomization/1o3_pdONLYpos_8o11_metNEGtumNEG_DESeq2_final/1o3_pdONLYpos_8o11_metNEGtumNEG_training_conditions.csv"
#sig_results_file <- "./Output/DESeq2/Results/significant_results_1o3_pdONLYpos8o11_metNEGtumNEG__pvalue0p05_lfc0p13_training.txt"
sig_results_file <- "./Output/DESeq2/logit/1-2-3-4PMpos8-9-11PMneg_training_genomewide_0p01_lfc0p13/individual_predictors_0p01.txt"
ver <- "whole_dataset_genomewide_logit_0p01_lfc0p13_nfold10"
alpha_value <- 1
number_of_runs = 50
nfolds_val = 10

#read in data, define what counts & conditions files
#note: glmnet expects a matrix as input, not a dataframe, 
#with genes in columns & samples in rows so we must transpose
counts_data_training <- read.csv(counts_name_training,row.names = 1)
counts_data_validation <- read.csv(counts_name_validation,row.names = 1)

rownames(counts_data_training) <- str_replace(c(rownames(counts_data_training)),'\\-',"_")
rownames(counts_data_training) <- str_replace(c(rownames(counts_data_training)),'\\.',"_")
rownames(counts_data_validation) <- str_replace(c(rownames(counts_data_validation)),'\\-',"_")
rownames(counts_data_validation) <- str_replace(c(rownames(counts_data_validation)),'\\.',"_")

counts_data_training <- t(as.matrix(counts_data_training))
counts_data_validation <- t(as.matrix(counts_data_validation))

meta_training <-  read.csv(meta_name_training,row.names=1)
meta_validation <- read.csv(meta_name_validation,row.names=1)


sig_results <- read.table(sig_results_file,header=TRUE)
sig_results_genes <- sig_results$gene

###Exclude genes where >5 samples >5 counts
#uncomment two lines to only look at sig genes from DESeq2 or other list
counts_data_training <- counts_data_training[,sig_results_genes]
counts_data_validation <- counts_data_validation[,sig_results_genes]
#flter out genes with <5 samples with >5 counts
counts_data_training_gt <- counts_data_training>5
gt_5_row <- colCounts(counts_data_training_gt, value=TRUE)
counts_data_training <- rbind(counts_data_training,gt_5_row)
counts_data_training <- counts_data_training[,counts_data_training["gt_5_row",]>5]
counts_data_training <- head(counts_data_training,-1)

#create conditions vector expected as input by glmnet (0,0,0,...,1,1,1,,...)
sample_count_training <- nrow(meta_training)
class1_name <- meta_training[2,1]
class2_name <- meta_training[nrow(meta_training),1]
class1_count_training <- sum(meta_training$condition == class1_name)
class2_count_training <- sum(meta_training$condition == class2_name)
conditions_vector=c(rep(0,class1_count_training),rep(1,class2_count_training))

if (file.exists(paste("./Output/glmnet/",class1_name,class2_name,"_",ver,"/",sep="")) == FALSE ) {
  dir.create(paste("./Output/glmnet/",class1_name,class2_name,"_",ver,"/",sep=""))
}

#create glmnet object
#fit <- glmnet(counts_data_training,conditions_vector, family = "binomial")

#run cross validation 200 times and save coefficients to an external file.
pb <- progress_bar$new(
  format = " downloading [:bar] :percent eta: :eta",
  total = number_of_runs, clear = FALSE, width= 200)

for (i in 1:number_of_runs) {
  pb$tick()
  set.seed(i)

### optimize lambda only, use alpha defined in config above.
  cvfit <- cv.glmnet(counts_data_training, conditions_vector, family = "binomial", type.measure = "auc", nfolds=5, alpha=alpha_value,  type.gaussian = "naive")
  tmp_coeffs <- coef(cvfit, s = "lambda.min")

### Optimize alpha & lambda together
#  cvfit <- cva.glmnet(counts_data_training, conditions_vector, family = "binomial", type.measure = "auc", nfolds=3,  type.gaussian = "naive")
#  alpha_v <- cvfit$alpha
#  error <- sapply(cvfit$modlist, function(mod) {min(mod$cvm)})
#  alpha_v <-  alpha_v[which.min(error)]
#  c("alpha",i,"= ",alpha_v)
#  tmp_coeffs <- coef(cvfit, s = "lambda.min", alpha = alpha_v)
  
  new_results <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)
  if(i==1){
    compiled_results = new_results
  }
  if(i>1){
    compiled_results = rbind(compiled_results,new_results)
  }
}

png(paste("./Output/glmnet/",class1_name,class2_name,"_",ver,"/","cvfit",class1_name,class2_name,"_",ver,".png", sep = ""), width = 1200, height = 900)
plot(cvfit)
dev.off()


#Save coefficients from 200 runs of cross validation models to file. 
#Ideally script would be expanded to find all predictors present in >X% of runs. 
#It does not currently do this.
write.table(
  compiled_results, 
  file=paste("./Output/glmnet/",class1_name,class2_name,"_",ver,"/",class1_name,"_",class2_name,"_",ver,"_coef",".txt", sep = ""), 
  sep="\t",
  quote = F
)
new_results <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)
