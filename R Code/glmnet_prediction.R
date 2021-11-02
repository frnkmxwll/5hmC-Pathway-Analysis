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
counts_name_training <- "./Output/Randomization/METneg_PMpos_DESeq2_whole_combatseq/METneg_PMpos_training_rawcounts.csv"
counts_name_validation <- "./Output/Randomization/METneg_PMpos_DESeq2_whole_combatseq/METneg_PMpos_validation_rawcounts.csv"
meta_name_training <- "./Output/Randomization/METneg_PMpos_DESeq2_whole_combatseq/METneg_PMpos_training_conditions.csv"
meta_name_validation <- "./Output/Randomization/METneg_PMpos_DESeq2_whole_combatseq/METneg_PMpos_validation_conditions.csv"
ver <- "v3"
alpha_value <- 0.5

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

#run cross validation 200 times and save coefficients to an external file.
for (i in 1:50) {
  progress(i,progress.bar=TRUE)
  set.seed(i)
  cvfit <- cv.glmnet(counts_data_training, conditions_vector, family = "binomial", type.measure = "auc", nfolds=7, alpha=alpha_value)
  tmp_coeffs <- coef(cvfit, s = "lambda.min")
  new_results <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)
  if(i==1){
    compiled_results = new_results
  }
  if(i>1){
    compiled_results = rbind(compiled_results,new_results)
  }
}

#Save coefficients from 200 runs of cross validation models to file. 
#Ideally script woul be expanded to find all predictors present in >X% of runs. 
#It does not currently do this.
write.table(
  compiled_results, 
  file=paste("./Output/glmnet/",class1_name,"_",class2_name,"_",ver,"_coef",".txt", sep = ""), 
  sep="\t",
  quote = F
)

counts_data_training_subset = counts_data_training[,c("IL12A","HOXA10","RPL19","P2RX5-TAX1BP3","IQCF5","CD1C")]
counts_data_validation_subset = counts_data_validation[,c("IL12A","HOXA10","RPL19","P2RX5-TAX1BP3","IQCF5","CD1C")]

set.seed(3)
fit <- glmnet(counts_data_training_subset,conditions_vector, family = "binomial", alpha = alpha_value)
cvfit <- cv.glmnet(counts_data_training_subset, conditions_vector, family = "binomial", type.measure = "auc", nfolds=7, alpha = alpha_value)
#plot(fit, xvar = "lambda", label = TRUE)

cvfit.predict.prob = predict(fit, counts_data_validation_subset, type = "response", s = cvfit$lambda.min)
#cvfit.predict.prob > 0.5
#predict.g [predict.g == "TRUE"] <- 1
#predict.g [predict.g == "FALSE"] <- 0
predicted_labels = as.data.frame(cvfit.predict.prob)

###PREDICT RESULTS FROM MODEL
#Set labels to 1 / 0, this is required to compare prediction results against ground truth.
actual_labels = meta_validation$condition
actual_labels [actual_labels == "PMpos"] <- 1
actual_labels [actual_labels == "METneg"] <- 0

prediction_results=prediction(predicted_labels[,1],actual_labels)
prediction_results

###PLOT ROC CURVE
perf <- performance(prediction_results,"tpr", "fpr")
plot(perf,colorize=FALSE,col="black")
lines(c(1,0),c(1,0),col="gray",lty=4)

#Determine AUC
auc_ROCR <- performance(prediction_results,"auc")
auc_ROCR@y.values[[1]]

#Display final model coefficients
tmp_coeffs <- coef(cvfit, s = "lambda.min")
tmp_coeffs
new_results <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)

