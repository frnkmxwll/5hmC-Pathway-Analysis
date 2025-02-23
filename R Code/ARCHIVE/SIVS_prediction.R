# Purpose of this script is to summarize into an epigenetic score using a multivariable logistic regression model 
# and the elastic net regularization (or LASSO)

#tutorials found here:
#https://cran.r-project.org/web/packages/glmnet/vignettes/glmnet.pdf
#https://stats.stackexchange.com/questions/72251/an-example-lasso-regression-using-glmnet-for-binary-outcome
#https://stackoverflow.com/questions/18130338/plotting-an-roc-curve-in-glmnet

### INSTALL LIBRARIES
# Setup, uncomment follow
# install.packages("dplyr")
#install.packages("glmnet", repos = "https://cran.us.r-project.org")
#install.packages("svMisc")
#install.packages("sivs", repos = "https://cran.rstudio.com")

library(sivs)
library(varhandle)

library(glmnet)
library(dplyr)
library(ROCR)
library(svMisc) 
library (glmnetUtils)
library (progress)
library (caret)
library(dplyr)
library(pheatmap)
library(ComplexHeatmap)#used for complex heatmaps


###CONFIGURATION
#set working directory, select where you extracted folder
setwd("~/5hmC-Pathway-Analysis/")

## 1o3 vs 8o11
#counts_name_training <- "./Output/DESeq2/Results/vst_counts_1o3_pdONLYpos8o11_metNEGtumNEG_pval0p05_lfc0p13_training.csv"
#counts_name_validation <- "./Output/DESeq2/Results/vst_counts_1o3_pdONLYpos8o11_metNEGtumNEG_pval0p05_lfc0p13_validation.csv"
#norm_counts_name_training <- "./Output/DESeq2/Results/normalized_counts_1o3_pdONLYpos8o11_metNEGtumNEG_pval0p05_lfc0p13_training.csv"
#norm_counts_name_validation <- "./Output/DESeq2/Results/normalized_counts_1o3_pdONLYpos8o11_metNEGtumNEG_pval0p05_lfc0p13_validation.csv"
#meta_name_training <- "./Output/Randomization/1o3_pdONLYpos_8o11_metNEGtumNEG_DESeq2_final/1o3_pdONLYpos_8o11_metNEGtumNEG_training_conditions.csv"
#meta_name_validation <- "./Output/Randomization/1o3_pdONLYpos_8o11_metNEGtumNEG_DESeq2_final/1o3_pdONLYpos_8o11_metNEGtumNEG_validation_conditions.csv"
#sig_results_file <- "./Output/DESeq2/Results/all_results_1o3_pdONLYpos8o11_metNEGtumNEG_pval0p05_lfc0p13_training.txt"
#nfolds_val = 3

## all data
#counts_name_training <- "./Output/DESeq2/Results/vst_counts_1-2-3-4PMpos8-9-11PMneg_pvalue0p01_lfc0p13_training.csv"
#counts_name_validation <- "./Output/DESeq2/Results/vst_counts_1-2-3-4PMpos8-9-11PMneg_pvalue0p01_lfc0p13_validation.csv"
#counts_name_training <- "./Output/DESeq2/Results/normalized_counts_1-2-3-4PMpos8-9-11PMneg_pval0p05_lfc0p26_training.csv"
#counts_name_validation <- "./Output/DESeq2/Results/normalized_counts_1-2-3-4PMpos8-9-11PMneg_pval0p05_lfc0p26_validation.csv"
#norm_counts_name_training <- "./Output/DESeq2/Results/normalized_counts_1-2-3-4PMpos8-9-11PMneg_pval0p05_lfc0p26_training.csv"
#norm_counts_name_validation <- "./Output/DESeq2/Results/normalized_counts_1-2-3-4PMpos8-9-11PMneg_pval0p05_lfc0p26_validation.csv"
#meta_name_training <- "./Output/Randomization/1-2-3-4PMpos_8-9-11PMneg_DESeq2_1234_8911_final_dropped/1-2-3-4PMpos_8-9-11PMneg_training_conditions.csv"
#meta_name_validation <- "./Output/Randomization/1-2-3-4PMpos_8-9-11PMneg_DESeq2_1234_8911_final_dropped/1-2-3-4PMpos_8-9-11PMneg_validation_conditions.csv"
#sig_results_file <- "./Output/DESeq2/Results/all_results_1-2-3-4PMpos8-9-11PMneg_pval0p05_lfc0p26_training.txt"
#nfolds_val = 6

## all data balanced
counts_name_training <- "./Output/DESeq2/Results/vst_counts_1-2-3-4PMpos8-9-11PMneg_pval0p05_lfc0p13_balanced_training.csv"
counts_name_validation <- "./Output/DESeq2/Results/vst_counts_1-2-3-4PMpos8-9-11PMneg_pval0p05_lfc0p13_balanced_validation.csv"
norm_counts_name_training <- "./Output/DESeq2/Results/normalized_counts_1-2-3-4PMpos8-9-11PMneg_pval0p05_lfc0p13_balanced_training.csv"
norm_counts_name_validation <- "./Output/DESeq2/Results/normalized_counts_1-2-3-4PMpos8-9-11PMneg_pval0p05_lfc0p13_balanced_validation.csv"
meta_name_training <- "./Output/Randomization/1-2-3-4PMpos_8-9-11PMneg_DESeq2_whole_balance/1-2-3-4PMpos_8-9-11PMneg_training_conditions.csv"
meta_name_validation <- "./Output/Randomization/1-2-3-4PMpos_8-9-11PMneg_DESeq2_whole_balance/1-2-3-4PMpos_8-9-11PMneg_validation_conditions.csv"
sig_results_file <- "./Output/DESeq2/Results/all_results_1-2-3-4PMpos8-9-11PMneg_pval0p05_lfc0p13_balanced_training.txt"
nfolds_val = 6

## all data no chemo
#counts_name_training <- "./Output/DESeq2/Results/normalized_counts_1-2-3-4PMpos8-9-11PMneg_nochemo_pval0p05_lfc0p26_training.csv"
#counts_name_validation <- "./Output/DESeq2/Results/normalized_counts_1-2-3-4PMpos8-9-11PMneg_nochemo_pval0p05_lfc0p26_validation.csv"
#norm_counts_name_training <- "./Output/DESeq2/Results/normalized_counts_1-2-3-4PMpos8-9-11PMneg_nochemo_pval0p05_lfc0p26_training.csv"
#norm_counts_name_validation <- "./Output/DESeq2/Results/normalized_counts_1-2-3-4PMpos8-9-11PMneg_nochemo_pval0p05_lfc0p26_validation.csv"
#meta_name_training <- "./Output/Randomization/1-2-3-4PMpos_8-9-11PMneg_DESeq2_nochemo/1-2-3-4PMpos_8-9-11PMneg_training_conditions.csv"
#meta_name_validation <- "./Output/Randomization/1-2-3-4PMpos_8-9-11PMneg_DESeq2_nochemo/1-2-3-4PMpos_8-9-11PMneg_validation_conditions.csv"
#sig_results_file <- "./Output/DESeq2/Results/all_results_1-2-3-4PMpos8-9-11PMneg_nochemo_pval0p05_lfc0p26_training.txt"

sig_type = "DESEQ" # "DESEQ" or "LOGIT"
#sig_results_file <- "./Output/DESeq2/logit/1o3_pdONLYpos8o11_metNEGtumNEG_training_genome_wide_vFINAL/individual_predictors.csv"
ver <- "0p05_lfc0p13_nfold6_balance_6_7"
alpha_value <- 1
number_of_runs = 50
pvalue_cutoff = 0.05
log2FoldChange_cutoff = 0.137504 # 0.137504 ~ 10%, 0.263034 ~ 20% change, 0.584963 ~ 50% change, 

seed_val = 10

#read in data, define what counts & conditions files
#note: glmnet expects a matrix as input, not a dataframe, 
#with genes in columns & samples in rows so we must transpose
counts_data_training <- read.csv(counts_name_training,row.names = 1)
counts_data_validation <- read.csv(counts_name_validation,row.names = 1)
norm_counts_data_training <- read.csv(norm_counts_name_training,row.names = 1)
norm_counts_data_validation <- read.csv(norm_counts_name_validation,row.names = 1)

#rownames(counts_data_training) <- str_replace(c(rownames(counts_data_training)),'\\-',"_")
#rownames(counts_data_training) <- str_replace(c(rownames(counts_data_training)),'\\.',"_")
#rownames(counts_data_validation) <- str_replace(c(rownames(counts_data_validation)),'\\-',"_")
#rownames(counts_data_validation) <- str_replace(c(rownames(counts_data_validation)),'\\.',"_")

counts_data_training <- as.data.frame(t(counts_data_training))
counts_data_validation <- as.data.frame(t(counts_data_validation))
norm_counts_data_training <- t(as.matrix(norm_counts_data_training))
norm_counts_data_validation <- t(as.matrix(norm_counts_data_validation))

meta_training <- read.csv(meta_name_training,row.names=1)
meta_validation <- read.csv(meta_name_validation,row.names=1)

# Define mapping
map <- list("pilot" = 0,
            "large_cohort" = 1,
            "female" = 0,
            "male" = 1,
            "pm_present" = 1,
            "pm_absent" = 0,
            "other_mets_absent" = 0,
            "other_mets_present" = 1,
            "primary_present" = 1,
            "primary_absent" = 0,
            "CRC" = 0,
            "HGA" = 1,
            "Yes" = 1,
            "No" = 0)

# Apply mapping
meta_training[] <- lapply(meta_training, function(x) {
  x <- as.character(x)
  for(i in 1:length(map)){
    x[x == names(map)[i]] <- map[[i]]
  }
  return(x)
})

meta_validation[] <- lapply(meta_validation, function(x) {
  x <- as.character(x)
  for(i in 1:length(map)){
    x[x == names(map)[i]] <- map[[i]]
  }
  return(x)
})

if (sig_type == "LOGIT"){
  sig_results <- read.table(sig_results_file,header=TRUE,sep=",")
  sig_results <- sig_results[sig_results[,"p.value"]<pvalue_cutoff,]
  sig_results <- na.omit(sig_results)
  sig_results_genes <- sig_results$term
}
if (sig_type == "DESEQ"){
  sig_results <- read.table(sig_results_file,header=TRUE)
  sig_results <- subset(sig_results, pvalue < pvalue_cutoff & abs(log2FoldChange) > log2FoldChange_cutoff)
  sig_results <- na.omit(sig_results)
  sig_results_genes <- rownames(sig_results)
}

#sig_results_genes <- str_replace(c(sig_results_genes),'\\.',"_")
#sig_results_genes <- str_replace(c(sig_results_genes),'\\-',"_")

###Exclude genes where >5 samples <5 counts
#uncomment two lines to only look at sig genes from DESeq2 or other list
#counts_data_training <- counts_data_training[,sig_results_genes]

# get the intersection of sig_genes and column names of counts_training
existing_sig_genes <- intersect(sig_results_genes, colnames(counts_data_training))
existing_sig_genes <- intersect(existing_sig_genes, colnames(counts_data_validation))

# subset counts_training to select only the existing columns
counts_data_training <- counts_data_training[,existing_sig_genes]
counts_data_validation <- counts_data_validation[, existing_sig_genes]

# # Add relevant clinical variables.
# # Convert row names to a column in both data frames
# counts_data_validation$id <- rownames(counts_data_validation)
# meta_validation$id <- rownames(meta_validation)
# counts_data_training$id <- rownames(counts_data_training)
# meta_training$id <- rownames(meta_training)
# 
# # Perform the merge by 'id'
# counts_data_validation <- merge(counts_data_validation, meta_validation[c("id", "primary_site", "chemo_6weeks")], by = "id", all.x = TRUE)
# counts_data_training <- merge(counts_data_training, meta_training[c("id", "primary_site", "chemo_6weeks")], by = "id", all.x = TRUE)
# 
# # Convert the 'primary_site' and 'chemo_6weeks' columns to numeric
# counts_data_validation$primary_site <- as.numeric(as.character(counts_data_validation$primary_site))
# counts_data_validation$chemo_6weeks <- as.numeric(as.character(counts_data_validation$chemo_6weeks))
# counts_data_training$primary_site <- as.numeric(as.character(counts_data_training$primary_site))
# counts_data_training$chemo_6weeks <- as.numeric(as.character(counts_data_training$chemo_6weeks))
# 
# # Convert the 'id' column back to row names if needed
# rownames(counts_data_validation) <- counts_data_validation$id
# counts_data_validation$id <- NULL
# rownames(counts_data_training) <- counts_data_training$id
# counts_data_training$id <- NULL

counts_data_training <- as.matrix(counts_data_training)
counts_data_validation <- as.matrix(counts_data_validation)

#create conditions vector expected as input by glmnet (0,0,0,...,1,1,1,,...)
sample_count_training <- nrow(meta_training)
class1_name <- meta_training[2,1]
class2_name <- meta_training[nrow(meta_training),1]
class1_count_training <- sum(meta_training$condition == class1_name)
class2_count_training <- sum(meta_training$condition == class2_name)
conditions_vector_training=c(rep(1,class1_count_training),rep(-1,class2_count_training))

class1_count_validation <- sum(meta_validation$condition == class1_name)
class2_count_validation <- sum(meta_validation$condition == class2_name)
conditions_vector_validation=c(rep(1,class1_count_validation),rep(-1,class2_count_validation))

conditions_vector_training_sivs = factor(conditions_vector_training)

sivs_obj <- sivs::sivs(
  x = counts_data_training,
  y = conditions_vector_training_sivs,
  family = "binomial",
  iter.count = 100,
  nfolds=6,
  alpha=1
  )

layout(mat = matrix(c(1,2,
                      3,3),
                    nrow = 2,
                    byrow = T))
{
  plot(sivs_obj)
  layout(1)
}

vimp_df = data.frame(sivs_obj$vimp)
selected_predictors = head(rownames(vimp_df),9)

selection_freq_df = data.frame(sivs_obj$selection.freq)
colnames(selection_freq_df) <- "freq"
selected_predictors = rownames(selection_freq_df)[selection_freq_df$freq >= 95][-1]

selected_predictors=sivs::suggest(sivs_obj, strictness = 0.01)
selected_predictors
#selected_predictors=c()

# Fit model to only top predictors
counts_data_validation_top <- counts_data_validation[,selected_predictors]
counts_data_training_top <- counts_data_training[,selected_predictors]
norm_counts_data_training_top <- norm_counts_data_training[,selected_predictors]
norm_counts_data_validation_top <- norm_counts_data_validation[,selected_predictors]

cvfit <- cv.glmnet(counts_data_training_top, conditions_vector_training, family = "binomial", type.measure = "auc", alpha=0.875, nfolds = nfolds_val, type.gaussian = "naive")
tmp_coeffs <- coef(cvfit, s = "lambda.min") # or lambda.min


# Compute the probabilities for the validation set
probabilities_val <- predict(cvfit, newx = counts_data_validation_top, s = "lambda.min", type = "response")
probabilities_tra <- predict(cvfit, newx = counts_data_training_top, s = "lambda.min", type = "response")


# Create a prediction object with ROCR
pred_val <- prediction(probabilities_val, conditions_vector_validation)
pred_tra <- prediction(probabilities_tra, conditions_vector_training)

# Generate performance metrics
roc_perf_val <- performance(pred_val, "tpr", "fpr")
roc_perf_tra <- performance(pred_tra, "tpr", "fpr")

# Calculate AUC
auc_obj <- performance(pred_tra, measure = "auc")
auc <- as.numeric(auc_obj@y.values)
print(paste("Training AUC: ", auc))

auc_obj <- performance(pred_val, measure = "auc")
auc <- as.numeric(auc_obj@y.values)
print(paste("Validation AUC: ", auc))

# Plot the ROC curve
plot(roc_perf_val, colorize = TRUE, print.cutoffs.at = seq(0,1,by=0.1), text.adj = c(-0.2,1.7),print.auc = TRUE)
# Add a diagonal line representing 50% AUC
abline(a = 0, b = 1, lty = 2, col = "gray")

# Plot the ROC curve
plot(roc_perf_tra, colorize = TRUE, print.cutoffs.at = seq(0,1,by=0.1), text.adj = c(-0.2,1.7),print.auc = TRUE)
# Add a diagonal line representing 50% AUC
abline(a = 0, b = 1, lty = 2, col = "gray")


## Emprically determine optimal alpha value:
alpha_values <- seq(0,1, by = 0.025)  # Change the step size as desired

# Initialize a list to store cv.glmnet objects for each alpha
cvfits <- list()

# Loop over alpha values
for (i in seq_along(alpha_values)) {
  cvfits[[i]] <- cv.glmnet(counts_data_training, conditions_vector_training, family = "binomial", type.measure = "auc", nfolds = nfolds_val, alpha = alpha_values[i], type.gaussian = "naive")
}

# Get the CV errors
cverrors <- sapply(cvfits, function(x) min(x$cvm))

# Find the alpha with the smallest CV error
optimal_alpha <- alpha_values[which.min(cverrors)]
c("optimal alpha:", optimal_alpha)

coefficient_table <- data.frame(run_number = numeric(), predictor = character(), coefficient = numeric())  # Table to store coefficients for each predictor
frequency_table <- data.frame(predictor = character(), frequency = numeric(), average_coefficient = numeric(), min_coefficient = numeric(), max_coefficient = numeric(), variance_coefficient = numeric())  # Frequency table


if (file.exists(paste("./Output/glmnet/",class1_name,class2_name,"_",ver,"/",sep="")) == FALSE ) {
  dir.create(paste("./Output/glmnet/",class1_name,class2_name,"_",ver,"/",sep=""))
}

#create glmnet object
#fit <- glmnet(counts_data_training,conditions_vector_training, family = "binomial")

#run cross validation 200 times and save coefficients to an external file.
pb <- progress_bar$new(
  format = " downloading [:bar] :percent eta: :eta",
  total = number_of_runs, clear = FALSE, width= 200)

for (i in 1:number_of_runs) {
  pb$tick()
  set.seed(i)
  
  ### optimize lambda only, use alpha defined in config above.
  cvfit <- cv.glmnet(counts_data_training, conditions_vector_training, family = "binomial", type.measure = "auc", nfolds=nfolds_val, alpha=optimal_alpha,  type.gaussian = "naive")
  tmp_coeffs <- coef(cvfit, s = "lambda.1se")
  
  new_results <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)
  if(i==1){
    compiled_results = new_results
  }
  if(i>1){
    compiled_results = rbind(compiled_results,new_results)
  }
  
  predictor_names <- rownames(tmp_coeffs)
  coefficients <- as.numeric(tmp_coeffs)
  run_numbers <- rep(i, length(predictor_names))
  coefficient_table <- rbind(coefficient_table, data.frame(run_number = run_numbers, predictor = predictor_names, coefficient = coefficients))
}

new_table <- coefficient_table %>%
  filter(coefficient != 0) %>%
  group_by(predictor) %>%
  summarise(
    freq_coef = n(),
    mean_coef = mean(coefficient),
    min_coef = min(coefficient),
    max_coef = max(coefficient)
  )

new_table$abs_mean <- abs(new_table$mean_coef)
new_table <- new_table %>%
  arrange(desc(freq_coef), desc(abs_mean))

new_table <- new_table %>%
  filter(predictor != "(Intercept)")

print(new_table, n=300)

# Select top predictors
selected_predictors <- new_table %>%
  filter(freq_coef > 0.99*number_of_runs) %>%
  pull(predictor)

# Fit model to only top predictors
counts_data_validation_top <- counts_data_validation[,selected_predictors]
counts_data_training_top <- counts_data_training[,selected_predictors]
norm_counts_data_training_top <- norm_counts_data_training[,selected_predictors]
norm_counts_data_validation_top <- norm_counts_data_validation[,selected_predictors]

cvfit <- cv.glmnet(counts_data_training_top, conditions_vector_training, family = "binomial", type.measure = "auc", alpha=0.1, nfolds = nfolds_val, type.gaussian = "naive")
tmp_coeffs <- coef(cvfit, s = "lambda.min") # or lambda.min


# Compute the probabilities for the validation set
probabilities_val <- predict(cvfit, newx = counts_data_validation_top, s = "lambda.min", type = "response")
probabilities_tra <- predict(cvfit, newx = counts_data_training_top, s = "lambda.min", type = "response")


# Create a prediction object with ROCR
pred_val <- prediction(probabilities_val, conditions_vector_validation)
pred_tra <- prediction(probabilities_tra, conditions_vector_training)

# Generate performance metrics
roc_perf_val <- performance(pred_val, "tpr", "fpr")
roc_perf_tra <- performance(pred_tra, "tpr", "fpr")

# Calculate AUC
auc_obj <- performance(pred_val, measure = "auc")
auc <- as.numeric(auc_obj@y.values)
print(paste("Validation AUC: ", auc))

auc_obj <- performance(pred_tra, measure = "auc")
auc <- as.numeric(auc_obj@y.values)
print(paste("Training AUC: ", auc))

# Plot the ROC curve
plot(roc_perf_val, colorize = TRUE, print.cutoffs.at = seq(0,1,by=0.1), text.adj = c(-0.2,1.7),print.auc = TRUE)
# Add a diagonal line representing 50% AUC
abline(a = 0, b = 1, lty = 2, col = "gray")

# Plot the ROC curve
plot(roc_perf_tra, colorize = TRUE, print.cutoffs.at = seq(0,1,by=0.1), text.adj = c(-0.2,1.7),print.auc = TRUE)
# Add a diagonal line representing 50% AUC
abline(a = 0, b = 1, lty = 2, col = "gray")

png(paste("./Output/glmnet/",class1_name,class2_name,"_",ver,"/","cvfit",class1_name,class2_name,"_",ver,"_a_",optimal_alpha,".png", sep = ""), width = 1200, height = 900)
plot(cvfit)
dev.off()


#Save coefficients from 200 runs of cross validation models to file. 
#Ideally script would be expanded to find all predictors present in >X% of runs. 
#It does not currently do this.
write.table(
  compiled_results, 
  file=paste("./Output/glmnet/",class1_name,class2_name,"_",ver,"/",class1_name,"_",class2_name,"_",ver,"_a_",optimal_alpha,"_coef",".txt", sep = ""), 
  sep="\t",
  quote = F
)

new_results <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)

### GENERATE HEATMAPS
counts_data_training_top_t = t(scale(counts_data_training_top))

# assuming your metadata table is called metadata
condition <- meta_training$condition
names(condition) <- rownames(meta_training)

# create a data frame for the annotation
annotation_df <- data.frame(
  Condition = condition
)

# Convert 'Condition' to a factor for the annotation
annotation_df$Condition <- as.factor(annotation_df$Condition)


# Convert 'Condition' to a factor for the annotation
groups <- unique(meta_training[c("condition")])
contrast_groups <- c("condition",groups[1,1], groups[2,1])
group1 = paste(contrast_groups[[3]])
group2 = paste(contrast_groups[[2]])
condition_table = c()
condition_table <- c("#ff9289","#00dae0")
names(condition_table) = c(group2,group1)
ha = HeatmapAnnotation(condition=annotation_df$Condition,col = list(condition = condition_table))

#pheatmap(, annotation_col = annotation_df)
heatmap_tra <- Heatmap(counts_data_training_top_t, 
        top_annotation = ha,
        name = "Training",
        #col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
        #col = colorRamp2(c(-4, -2,-1,0,1, 2,4), c("#4676b5", "#82b1d3","#dbeff6","#fefebd","#fee395", "#fc9961","#d73027")),
        #clustering_distance_rows = "euclidean",
        clustering_distance_columns = "euclidean")
        #column_title = heatmap_title)

## VALIDATION
counts_data_validation_top_t = t(scale(counts_data_validation_top))

# assuming your metadata table is called metadata
condition <- meta_validation$condition
names(condition) <- rownames(meta_validation)

# create a data frame for the annotation
annotation_df <- data.frame(
  Condition = condition
)

# Convert 'Condition' to a factor for the annotation
annotation_df$Condition <- as.factor(annotation_df$Condition)
# Convert 'Condition' to a factor for the annotation
groups <- unique(meta_training[c("condition")])
contrast_groups <- c("condition",groups[1,1], groups[2,1])
group1 = paste(contrast_groups[[3]])
group2 = paste(contrast_groups[[2]])
condition_table = c()
condition_table <- c("#ff9289","#00dae0")
names(condition_table) = c(group2,group1)
ha = HeatmapAnnotation(condition=annotation_df$Condition,col = list(condition = condition_table))

#pheatmap(, annotation_col = annotation_df)
heatmap_val <- Heatmap(counts_data_validation_top_t, 
        top_annotation = ha,
        name = "Validation",
        #col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
        #col = colorRamp2(c(-4, -2,-1,0,1, 2,4), c("#4676b5", "#82b1d3","#dbeff6","#fefebd","#fee395", "#fc9961","#d73027")),
        #clustering_distance_rows = "euclidean",
        clustering_distance_columns = "euclidean")
#column_title = heatmap_title)

# Create the combined plot with shared row clustering
heatmap_tra + heatmap_val

