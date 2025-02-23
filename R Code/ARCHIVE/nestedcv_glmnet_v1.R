###dataset_group_randomizer.R###
# Purpose: Split a DESeq2-ready file into training, validation, and holdout sets.
# ### INSTALL LIBRARIES IF NECESSARY
# install.packages("caret")
# install.packages("dplyr")
# install.packages("BiocManager")
# BiocManager::install("DESeq2")
# BiocManager::install("limma")
# BiocManager::install("sva")
# install.packages("doParallel")
# install.packages("foreach")
# install.packages("progress")
# install.packages("glmnet")
# install.packages("ROCR")
# install.packages("future")
# install.packages("future.apply")
# install.packages("plotly")
# install.packages("progressr")
# install.packages("future")
#install.packages('Rfast')
#install.packages("nestedcv")
# install.packages("pbapply")
# BiocManager::install("kaleido")
# BiocManager::install("ChIPseeker")
# BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
# BiocManager::install("org.Hs.eg.db")

## Kaleidoscope for plotly exports
#install.packages('reticulate')
#reticulate::install_miniconda()
#reticulate::conda_install('r-reticulate', 'python-kaleido')
#reticulate::conda_install('r-reticulate', 'plotly', channel = 'plotly')
#reticulate::use_miniconda('r-reticulate')

### Load Libraries
library(caret)  # for glmnet model training
library(dplyr)
#library(limma)
library(DESeq2)  # for normalization and DGE analysis
#library(sva)  # for combat-seq batch normalization
library(future)
library(future.apply)

### Glmnet libraries
library(doParallel)
library(BiocParallel)
library(foreach)
#library(progress)
library(glmnet)
library(nestedcv)
library(caret)
library(ROCR)  # Add this line to load the ROCR package

library(future)
library(future.apply)
library(plotly)
library(limma)
library(sva)

# For Annotations
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

library(progressr)
library(pbapply)

source("~/5hmC-Pathway-Analysis/R Code/repeated_glmnet_funcs_v1.R")
handlers("txtprogressbar")

### Custom nestcv
#setwd("C:/Users/ymali/AppData/Local/R/win-library/4.2/nestedcv_YM/")
#library(devtools)
#devtools::install("C:/Users/ymali/AppData/Local/R/win-library/4.2/nestedcv_YM/")
#library(nestedcvYM)





#handlers(global = TRUE)

#closeAllConnections()

### CONFIGURATION
## General Settings
# Detect operating system and store in variable
os <- Sys.info()[['sysname']]

if (os == "Windows") {
  # Windows
  root_folder <- "C:/Users/ymali/PycharmProjects/peritoneal/R/data/"
  output_folder <- "C:/Users/ymali/PycharmProjects/peritoneal/R/outputs/"
} else {
  # Linux
  root_folder <- "./Output/Raw Data Processing/PMneg_PMpos_06022024_hmr_nocombat/"
  output_folder <- "./Output/repeated_holdout/"
}
# Location of raw_data_processing output

filename_trunk <- "PMneg_PMpos_DESeq2"
counts_name <- paste0(root_folder, filename_trunk, "_rawcounts.csv")
normcounts_name <- paste0(root_folder, filename_trunk, "_normcounts.csv")
meta_name <- paste0(root_folder, filename_trunk, "_conditions.csv")
file_version <- "hmr_nocombat_602020_06022024_CRCHGA_seed1"
save_final_files <- FALSE
debug_progress <- TRUE
debug_data_summary <- FALSE
random_seed <- 10
set.seed(random_seed)


## Partition Settings
training_fraction <- 70
validation_fraction <- 30
holdout_fraction <- 0
min_samples <- 0.8
min_count <- 5
create_sig_training_counts <- FALSE
combat_norm <- FALSE
interaction_params <- c("condition", "primary_present","batch")

## DESeq2 Settings
cutoff_type = 2 # 0=padj cutoff, default; 1=lfc & pvalue cutoff; 2 = top stat genes
padj_cutoff = 0.01 # 0.1 default
pvalue_cutoff = 0.001
lfc_cutoff = 0.138 # 0.137504 ~ 10%, 0.263034 ~ 20% change, 0.32 ~ 25%, 0.415 ~ 33%, 0.58 ~ 50% change,
num_top_genes = 1000
design_formula <- ~ condition + primary_present + batch

## glmnet Settings
nfolds_val = 5
rep_val = 1
feat_occ_thresholds <- seq(0.5, 1, by = 0.1) #c(0.7,0.75,0.8,0.85,0.9,0.95,0.98,1)#
alpha_values <- seq(0.5, 1, by = 0.1)
standardize_bit = FALSE # should caret center and standardize before training models?
number_of_runs <- 5
part_rand_seed = seq(100,101, by=1)

rapid = TRUE
if(rapid){
  num_raw_genes = 1000
  num_top_genes = 100
  feat_occ_thresholds <- seq(0.9, 1, by = 0.05) #c(0.7,0.75,0.8,0.85,0.9,0.95,0.98,1)#
  alpha_values <- seq(0.9, 1, by = 0.1)
  standardize_bit = FALSE # should caret center and standardize before training models?
  number_of_runs <- 10
  part_rand_seed = seq(100,102, by=1)
}

### Main Script

## Load data
# Load and preprocess data
data_prepared <- load_and_prepare_data(counts_name, normcounts_name, meta_name, debug_data_summary, debug_progress, rapid, num_raw_genes)
# Add stratum collum based on interaction parameters
data_prepared$meta$stratum <- do.call(interaction, data_prepared$meta[interaction_params])
conditions <- define_conditions(data_prepared$meta)

# for each random_seed in part_rand_seed
#for(random_seed in part_rand_seed){
# Temp set rand seed for testing
random_seed = 1
# Set seed
set.seed(random_seed)
# Split data into training and validation set, order them and normalize.
partitions <- partition_data(data_prepared$meta, interaction_params, training_fraction, validation_fraction, holdout_fraction, random_seed, debug_progress)
meta_training <- partitions$training_meta
meta_validation <- partitions$validation_meta
meta_holdout <- partitions$holdout_meta

ordered_filtered_counts <- order_and_filter_data(meta_training, meta_validation, partitions$holdout_meta, data_prepared$counts_data, min_count, min_samples, debug_progress)
normalized_counts <- normalize_counts(ordered_filtered_counts$training, ordered_filtered_counts$validation, ordered_filtered_counts$holdout, combat_norm, debug_progress)


# Prepare class vectors
class1_name <- meta_training[2, 1]
class2_name <- meta_training[nrow(meta_training), 1]
conditions_vector_training <- create_conditions_vector(meta_training, class1_name, class2_name)
conditions_vector_validation <- create_conditions_vector(meta_validation, class1_name, class2_name)
conditions_vector_training_fac <- as.factor(make.names(conditions_vector_training))
conditions_vector_validation_fac <- as.factor(make.names(conditions_vector_validation))
class_weights <- ifelse(conditions_vector_training == 1, 1 / sum(meta_training$condition == class1_name), 1 / sum(meta_training$condition == class2_name))

# number of cores available
num_cores <- parallel::detectCores()

coeff_freq_df <- data.frame()

start_time <- Sys.time()
for(random_seed in part_rand_seed){
  # Set seed
  set.seed(random_seed, "L'Ecuyer-CMRG")
  # Split data into training and validation set, order them and normalize.
  
  
  repeated_fit <- nestcv.glmnet(
    conditions_vector_training_fac,
    t(normalized_counts$training_normcounts),
    family = "binomial",
    weights = class_weights,
    alphaSet = alpha_values,
    cv.cores = floor(num_cores/2),
    n_outer_folds = 5, 
    n_inner_folds = 10,
    #penalty.factor = rep(1,ncol(training_normcounts_matrix)),
    finalCV = FALSE, #use median results
    verbose = TRUE,
    filterFUN = deseq2_filter,
    filter_options = list(
      training_data = ordered_filtered_counts$training,
      design_formula = design_formula,
      cutoff_type = cutoff_type,
      padj_cutoff = padj_cutoff,
      pvalue_cutoff = pvalue_cutoff,
      lfc_cutoff = lfc_cutoff,
      num_top_genes = num_top_genes
    )
  )
  
  temp_coeff_df <- repeated_fit$final_coef
  # convert rownames of temp_coeff_df to a column called feature
  temp_coeff_df$feature <- rownames(temp_coeff_df)
  # drop any rows where feature contains "intercept" (case insensitive)
  temp_coeff_df <- temp_coeff_df[!grepl("intercept", temp_coeff_df$feature, ignore.case = TRUE),]
  # drop rownames
  rownames(temp_coeff_df) <- NULL
  
  # if coeff_freq_df is empty dataframe, assign the first iteration's coefficients_df to it
  if (nrow(coeff_freq_df) == 0) {
    coeff_freq_df <- temp_coeff_df
  }  else {
    # append rownames to the end of coeff_freq_df (without overwriting), dropping any rows where rowname == "(Intercept)"
    coeff_freq_df <- rbind(coeff_freq_df, temp_coeff_df)
  }
}
end_time <- Sys.time()
repeatcv_op_time <- end_time - start_time
print(paste("Time taken for repeatedcv: ", repeatcv_op_time))

# Initialize thresh_perf_df
thresh_perf_df <- data.frame()
best_auc <- 0
for (threshold in feat_occ_thresholds) {
  # Filter the coefficients based on the threshold
  filtered_coefficients <- coeff_freq_df %>%
    group_by(feature) %>%
    filter(n() >= threshold) %>%
    summarise(
      freq_coef = n(),
      mean_coef = mean(coef),
      min_coef = min(coef),
      max_coef = max(coef)
    ) %>%
    arrange(desc(freq_coef), desc(abs(mean_coef)))
  
  # Create a list of the top features
  top_features <- filtered_coefficients$feature
  
  # Subset the normalized counts data using top_features
  norm_counts_data_training_top <- subset_normalized_counts(normalized_counts$training_normcounts, top_features)
  
  repeated_fit <- nestcv.glmnet(
    conditions_vector_training_fac,
    norm_counts_data_training_top,
    family = "binomial",
    weights = class_weights,
    alphaSet = seq(0, 1, 0.1),
    cv.cores = 16,
    n_outer_folds = 5, 
    n_inner_folds = 10,
    #penalty.factor = rep(1,ncol(training_normcounts_matrix)),
    verbose = TRUE,
    filterFUN = NULL, 
    filter_options = NULL
  )
  
  # Predict probabilities
  probabilities_tra <- predict(repeated_fit, newdata = as.matrix(norm_counts_data_training_top), type = "response")
  
  # Create prediction objects
  pred_tra <- prediction(probabilities_tra, conditions_vector_training_fac)
  
  # Calculate ROC performance
  roc_perf_tra <- performance(pred_tra, "tpr", "fpr")
  
  # Calculate AUC values
  auc_tra <- as.numeric(performance(pred_tra, measure = "auc")@y.values)
  
  # Use youden's J to find optimal threshold
  perf_sens <- performance(pred_tra, measure = "sens")
  perf_spec <- performance(pred_tra, measure = "spec")
  
  sens_values <- unlist(perf_sens@y.values)
  spec_values <- unlist(perf_spec@y.values)
  class_thresholds <- unlist(perf_sens@x.values)
  youden_j <- sens_values + spec_values - 1
  optimal_index <- which.max(youden_j)
  opt_class_threshold <- class_thresholds[optimal_index]
  cat("Optimal threshold (Youden's J):", opt_class_threshold, "\n")
  
  # Calculate confusion matrix
  cm_tra <- confusionMatrix(factor(ifelse(probabilities_tra > opt_class_threshold, "X1", "X0")), conditions_vector_training_fac)
  
  # Calculate metrics
  metrics_row <- data.frame(
    threshold = threshold,
    num_features = length(top_features),
    auc_tra = auc_tra,
    sensitivity_tra = cm_tra$byClass["Sensitivity"],
    specificity_tra = cm_tra$byClass["Specificity"],
    recall_tra = cm_tra$byClass["Recall"],
    precision_tra = cm_tra$byClass["Precision"],
    accuracy_tra = cm_tra$overall["Accuracy"],
    f1_tra = (2*(cm_tra$byClass["Precision"] * cm_tra$byClass["Recall"])/(cm_tra$byClass["Recall"] + cm_tra$byClass["Precision"])),
    opt_class_thresh = opt_class_threshold
  )
  
  # if thresh_perf_df is empty dataframe, assign the first iteration's metrics_row to it
  if (nrow(thresh_perf_df) == 0) {
    thresh_perf_df <- metrics_row
  }
  else {# append metrics_row to the end of thresh_perf_df
    thresh_perf_df <- rbind(thresh_perf_df, metrics_row)
  }
  
  # Append metrics_row to thresh_perf_df
  thresh_perf_df <- rbind(thresh_perf_df, metrics_row)
  
  # if auc_tra is greater than best_auc, update store the best model
  if (auc_tra > best_auc) {
    best_auc <- auc_tra
    best_threshold <- threshold
    best_repeated_fit <- repeated_fit
  }
}



predy <- predict(best_repeated_fit, newdata = t(normalized_counts$validation_normcounts), s = "lambda.min", type = "class")
predy <- as.vector(predy)
predyp <- predict(best_repeated_fit, newdata = t(normalized_counts$validation_normcounts), s = "lambda.min", type = "response")
predyp <- as.vector(predyp)
output <- data.frame(testy = make.names(conditions_vector_validation), predy = predy, predyp = predyp)

## Results on test set
## shows bias since univariate filtering was applied to whole dataset
predSummary(output)

# Outer CV ROC
plot(best_repeated_fit$roc, main = "Outer fold ROC", font.main = 1, col = 'blue')
legend("bottomright", legend = paste0("AUC = ", signif(pROC::auc(best_repeated_fit$roc), 3)), bty = 'n')

# Inner CV ROC
best_repeated_fit.inroc <- innercv_roc(best_repeated_fit)
plot(best_repeated_fit.inroc, main = "Inner fold ROC", font.main = 1, col = 'red')
legend("bottomright", legend = paste0("AUC = ", signif(pROC::auc(best_repeated_fit.inroc), 3)), bty = 'n')

# Plot ROC curve for test set
valroc <- pROC::roc(output$testy, output$predyp, direction = "<", quiet = TRUE)
testroc <- plot(best_repeated_fit$roc, col = "blue", main = "ROC", las = 1)

lines(valroc, col = 'blue')
lines(testroc, col = 'red')
legend('bottomright', legend = c("Left-out inner CV folds", 
                                 "Test partition"), 
       col = c("red", "blue"), lty = 1, lwd = 2, bty = "n")



coefficients_df <- best_repeated_fit$final_coef

# convert row names to column called feature

coefficients_df$feature <- rownames(coefficients_df)


# Drop all rows where the feature is "(Intercept)" and save to new dataframe coefficients_df_features
coefficients_df_features <- coefficients_df[coefficients_df$feature != "(Intercept)", ]

# Add gene symbols to coefficients_df_features
coefficients_df_features <- annotate_features(coefficients_df_features)