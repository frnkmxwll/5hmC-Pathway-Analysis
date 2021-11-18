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
# install.packages("svMisc")

library(glmnet)
library(dplyr)
library(ROCR)
library(svMisc) 

###CONFIGURATION
#set working directory, select where you extracted folder
setwd("~/5hmC-Pathway-Analysis/")
counts_name_training <- "./Output/Randomization/METneg_PMonlyPOS_DESeq2_whole_combatseq_50/METneg_PMonlyPOS_training_rawcounts.csv"
counts_name_validation <- "./Output/Randomization/METneg_PMonlyPOS_DESeq2_whole_combatseq_50/METneg_PMonlyPOS_validation_rawcounts.csv"
meta_name_training <- "./Output/Randomization/METneg_PMonlyPOS_DESeq2_whole_combatseq_50/METneg_PMonlyPOS_training_conditions.csv"
meta_name_validation <- "./Output/Randomization/METneg_PMonlyPOS_DESeq2_whole_combatseq_50/METneg_PMonlyPOS_validation_conditions.csv"
ver <- "v3_a05_50"
alpha_value <- 0.5
cv_folds <- 3 #number of folds to run in each cross validation
cv_runs <- 200 #number of cross validations to run to build model
run_cutoff <- 0.8 #genes must feature in this percent of cross-validation run to feature in the model

#read in data, define what counts & conditions files
#note: glmnet expects a matrix as input, not a dataframe, 
#with genes in columns & samples in rows so we must transpose
counts_data_training <- t(as.matrix(read.csv(counts_name_training,row.names = 1)))
counts_data_validation <- t(as.matrix(read.csv(counts_name_validation,row.names = 1)))

meta_training <-  read.csv(meta_name_training,row.names=1)
meta_validation <- read.csv(meta_name_validation,row.names=1)

#create conditions vector expected as input by glmnet (0,0,0,...,1,1,1,,...)
sample_count_training <- nrow(meta_training)
class1_name <- meta_training[2,1]
class2_name <- meta_training[nrow(meta_training),1]
class1_count_training <- sum(meta_training$condition == class1_name)
class2_count_training <- sum(meta_training$condition == class2_name)
conditions_vector=c(rep(0,class1_count_training),rep(1,class2_count_training))

#create glmnet object
#fit <- glmnet(counts_data_training,conditions_vector, family = "binomial")

#run cross validation the number of times defined in config cv_runs value and save coefficients to an external file.
for (i in 1:cv_runs) {
  progress(i, max.value = cv_runs)
  set.seed(i)
  cvfit <- cv.glmnet(counts_data_training, conditions_vector, family = "binomial", nfolds=cv_folds, alpha=alpha_value, type.measure = "auc")

  tmp_coeffs <- coef(cvfit, s = "lambda.min")
  new_results <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)
  if(i==1){
    compiled_results = new_results
  }
  if(i>1){
    compiled_results = rbind(compiled_results,new_results)
  }
}

#plot(cvfit)

#Save coefficients from cv_runs of cross validation models to file. if necessary for debugging
#write.table(
#  compiled_results, 
#  file=paste("./Output/glmnet/",class1_name,"_",class2_name,"_",ver,"_all_coef",".txt", sep = ""), 
#  sep="\t",
#  quote = F
#)

#Count occurences of each unique gene across all runs, and count their occurences
compiled_results <- compiled_results[!(compiled_results$name=="(Intercept)"),] #remove unnecessary "Intercept" rows
unique_gene_occurences <- table(compiled_results$name) #count number of unique occurrences of each gene
unique_gene_occurences_df <- as.data.frame(rbind(unique_gene_occurences)) #convert table to dataframe
unique_gene_occurences_df <- rbind(unique_gene_occurences_df,unique_gene_occurences_df[c("unique_gene_occurences"),]/cv_runs) #Add row calculating occurence rate of each gene
rownames(unique_gene_occurences_df) <- c("occurence_count","occurence_rate") #update occurence rate row name
subset_meeting_cutoff <- colnames(unique_gene_occurences_df[,c(unique_gene_occurences_df[2,]>= run_cutoff)]) #select only genes meeting inclusing cutoff run_cutoff
#unique_gene_occurences_df[2,]>= run_cutoff
counts_data_training_subset = counts_data_training[,subset_meeting_cutoff]
counts_data_validation_subset = counts_data_validation[,subset_meeting_cutoff]
subset_meeting_cutoff

set.seed(3)
fit <- glmnet(counts_data_training_subset,conditions_vector, family = "binomial", alpha = alpha_value)
cvfit <- cv.glmnet(counts_data_training_subset, conditions_vector, family = "binomial", type.measure = "auc", nfolds=cv_folds, alpha = alpha_value)
#plot(fit, xvar = "lambda", label = TRUE)

cvfit.predict.prob = predict(fit, counts_data_validation_subset, type = "response", s = cvfit$lambda.min)
#cvfit.predict.prob > 0.5
#predict.g [predict.g == "TRUE"] <- 1
#predict.g [predict.g == "FALSE"] <- 0
predicted_labels = as.data.frame(cvfit.predict.prob)

###PREDICT RESULTS FROM MODEL
#Set labels to 1 / 0, this is required to compare prediction results against ground truth.
actual_labels = meta_validation$condition
actual_labels [actual_labels == class2_name] <- 1
actual_labels [actual_labels == class1_name] <- 0

prediction_results=prediction(predicted_labels[,1],actual_labels)
prediction_results

###PLOT ROC CURVE
perf <- performance(prediction_results,"tpr", "fpr")

#Determine AUC
auc_ROCR <- performance(prediction_results,"auc")

#Display final model coefficients
tmp_coeffs <- coef(cvfit, s = "lambda.min")
new_results <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)

###OUTPUT ALL FILES
#ROC curve
png(paste("./Output/glmnet/",class1_name,"_",class2_name,"_",ver,"_2ROC",".png", sep = ""), width = 900, height = 900)
plot(
  perf,colorize=FALSE,
  col="black", 
  main = paste(
    class2_name," vs ",class1_name,"|",
    "Alpha=",alpha_value,"|",
    "AUC=",round(auc_ROCR@y.values[[1]],2)
    )
  )
lines(c(1,0),c(1,0),col="gray",lty=4)
dev.off()

#all coeffiecients
#write.table(
#  tmp_coeffs, 
#  file=paste("./Output/glmnet/",class1_name,"_",class2_name,"_",ver,"_final_coefs",".txt", sep = ""), 
#  sep="\t",
#  quote = F
#)

#Config file
config <- c(
  paste("training counts file name:", counts_name_training), 
  paste("training conditions file name:", meta_name_training), 
  paste("validation counts file name:", counts_name_validation), 
  paste("validation conditions file name:", meta_name_validation), 
  paste("output file name:", ver),
  paste("alpha value (1=lasso regression, 0.5=elastic net , 0=ridge regression:", alpha_value),
  paste("number of folds to run in each cross validation:", cv_folds),
  paste("number of cross validations ran to build model:", cv_runs),
  paste("genes must feature in this percent of cross-validation run to feature in the model:", run_cutoff)
)
write.table(config, file=paste("./Output/glmnet/",class1_name,"_",class2_name,"_",ver,"_config",".txt", sep = ""), sep="\t", quote=F, col.names=NA)

#final coefficients
write.table(
  new_results, 
  file=paste("./Output/glmnet/",class1_name,"_",class2_name,"_",ver,"_final_coefs",".txt", sep = ""), 
  sep="\t",
  quote = F
)

