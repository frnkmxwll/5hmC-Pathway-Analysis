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
# Install required packages if not already installed
library(glmnet)
library(dplyr)
library(ROCR)
#library(svMisc) 
#library (glmnetUtils)
library (progress)
library (caret)
library(dplyr)
library(doParallel)

### For plotting
library(circlize)
library(pheatmap)
library(ComplexHeatmap)#used for complex heatmaps

### For annotations
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ChIPseeker)

###CONFIGURATION
#set working directory, select where you extracted folder
setwd("~/5hmC-Pathway-Analysis/")

## CRC only
sig_results_file <- "./Output/DESeq2/Results/all_results_PMnegPMpos_ntop 500 _ CRCHGA_nocombat_602020_06022024_cbpp_seed126.txt"
root_folder <- "./Output/Randomization/PMneg_PMpos_DESeq2_hmr_nocombat_602020_06022024_CRCHGA_cbpp_seed126/"
filename_trunk <- "PMneg_PMpos"

counts_name_training <- paste0(root_folder, filename_trunk, "_training_rawcounts.csv")
norm_counts_name_training <- paste0(root_folder, filename_trunk, "_training_normcounts.csv")
meta_name_training <- paste0(root_folder, filename_trunk, "_training_conditions.csv")

counts_name_validation <- paste0(root_folder, filename_trunk, "_validation_rawcounts.csv")
norm_counts_name_validation <- paste0(root_folder, filename_trunk, "_validation_normcounts.csv")
meta_name_validation <- paste0(root_folder, filename_trunk, "_validation_conditions.csv")

counts_name_holdout <- paste0(root_folder, filename_trunk, "_holdout_rawcounts.csv")
norm_counts_name_holdout <- paste0(root_folder, filename_trunk, "_holdout_normcounts.csv")
meta_name_holdout <- paste0(root_folder, filename_trunk, "_holdout_conditions.csv")

#counts_name_holdout <- "./Output/Raw Data Processing/LGA_PMneg_LGA_PMpos_05252024_hmr_combat/LGA_PMneg_LGA_PMpos_DESeq2_rawcounts.csv"
#norm_counts_name_holdout <- "./Output/Raw Data Processing/LGA_PMneg_LGA_PMpos_05252024_hmr_combat/LGA_PMneg_LGA_PMpos_DESeq2_normcounts.csv"
#meta_name_holdout <- "./Output/Raw Data Processing/LGA_PMneg_LGA_PMpos_05252024_hmr_combat/LGA_PMneg_LGA_PMpos_DESeq2_conditions.csv"

sig_type = "DESEQ" # "DESEQ" or "LOGIT"

ver.text <- "CRCHGA_hmr_combat_602020"
number_of_runs = 50
#pvalue_cutoff = 0.01
#log2FoldChange_cutoff =  0.26  # 0.138 ~ 10%, 0.263 ~ 20% change, 0.32 ~ 25%, 0.415 ~ 33%, 0.58 ~ 50% change, 
num_top_genes <- 1000 # Define the number of top genes to select
nfolds_val = 5
seed_val = 11
feat_occ_thresholds <- seq(0.7, 1, by = 0.05) #c(0.7,0.75,0.8,0.85,0.9,0.95,0.98,1)#
alpha_values <- seq(0.1, 1, by = 0.1)
lambda_values <- 10^seq(-4, 1, length = 20)
train_single = FALSE
standardize_bit = FALSE

#ver <- paste("pval",pvalue_cutoff,"_lfc",log2FoldChange_cutoff,"_nfold",nfolds_val,"+",ver.text)
ver <- paste("ngenes",num_top_genes,"_nfold",nfolds_val,"+",ver.text)

#read in data, define what counts & conditions files
#note: glmnet expects a matrix as input, not a dataframe, 
#with genes in columns & samples in rows so we must transpose
norm_counts_data_training <- read.csv(norm_counts_name_training,row.names = 1)
norm_counts_data_validation <- read.csv(norm_counts_name_validation,row.names = 1)
norm_counts_data_holdout <- read.csv(norm_counts_name_holdout,row.names = 1)

norm_counts_data_training <- t(as.matrix(norm_counts_data_training))
norm_counts_data_validation <- t(as.matrix(norm_counts_data_validation))
norm_counts_data_holdout <- t(as.matrix(norm_counts_data_holdout))

meta_training <- read.csv(meta_name_training,row.names=1)
meta_validation <- read.csv(meta_name_validation,row.names=1)
meta_holdout <- read.csv(meta_name_holdout,row.names=1)

# Define mapping
map <- list("pilot" = 0,
            "large_cohort" = 1,
            "female" = 0,
            "male" = 1,
            "other_mets_absent" = 0,
            "other_mets_present" = 1,
            "Yes" = 1,
            "No" = 0,
            "LGA_PMneg" = "PMneg",
            "LGA_PMpos" = "PMpos",
            "HGA_PMneg" = "PMneg",
            "HGA_PMpos" = "PMpos"
            )

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

meta_holdout[] <- lapply(meta_holdout, function(x) {
  x <- as.character(x)
  for(i in 1:length(map)){
    x[x == names(map)[i]] <- map[[i]]
  }
  return(x)
})

# Read DESeq2 results
sig_results <- read.table(sig_results_file, header = TRUE)
sig_results <- na.omit(sig_results)

# Rank genes by the stat column and select the top N genes
sig_results <- sig_results %>%
  arrange(desc(abs(stat)))

top_genes <- head(sig_results, num_top_genes)
sig_results_genes <- rownames(top_genes)

# get the intersection of sig_genes and column names of counts_training
existing_sig_genes <- intersect(sig_results_genes, colnames(norm_counts_data_training))
existing_sig_genes <- intersect(existing_sig_genes, colnames(norm_counts_data_validation))
existing_sig_genes <- intersect(existing_sig_genes, colnames(norm_counts_data_holdout))

# subset counts_training to select only the existing columns
norm_counts_data_training <- norm_counts_data_training[,existing_sig_genes]
norm_counts_data_validation <- norm_counts_data_validation[, existing_sig_genes]
norm_counts_data_holdout <- norm_counts_data_holdout[, existing_sig_genes]

#create conditions vector expected as input by glmnet (0,0,0,...,1,1,1,...)
class1_name <- meta_training[2,1]
class2_name <- meta_training[nrow(meta_training),1]
class1_count_training <- sum(meta_training$condition == class1_name)
class2_count_training <- sum(meta_training$condition == class2_name)
conditions_vector_training <- c(rep(1, class1_count_training), rep(0, class2_count_training))

class1_count_validation <- sum(meta_validation$condition == class1_name)
class2_count_validation <- sum(meta_validation$condition == class2_name)
conditions_vector_validation <- c(rep(1, class1_count_validation), rep(0, class2_count_validation))

class1_count_holdout <- sum(meta_holdout$condition == class1_name)
class2_count_holdout <- sum(meta_holdout$condition == class2_name)
conditions_vector_holdout <- c(rep(1, class1_count_holdout), rep(0, class2_count_holdout))

conditions_vector_training_fac <- as.factor(make.names(conditions_vector_training))
conditions_vector_validation_fac <- as.factor(make.names(conditions_vector_validation))
conditions_vector_holdout_fac <- as.factor(make.names(conditions_vector_holdout))

### Define class weights
class_weights <- ifelse(conditions_vector_training == 1, 
                        1 / class1_count_training, 
                        1 / class2_count_training)

set.seed(123)

### TRAIN MANY MODELS

if (file.exists(paste("./Output/glmnet/", class1_name, class2_name, "_", ver, "/", sep = "")) == FALSE) {
  dir.create(paste("./Output/glmnet/", class1_name, class2_name, "_", ver, "/", sep = ""))
}

optimal_auc_val <- -Inf
optimal_model <- NULL
optimal_alpha <- NA
optimal_threshold <- NA
best_model_metrics <- NULL
roc_data_tra_best <- list()
roc_data_val_best <- list()

model_metrics <- data.frame()

library(doParallel)
library(foreach)
library(progress)
library(dplyr)
library(glmnet)

for (alpha_value in alpha_values) {
  print("Starting on:")
  print(alpha_value)
  
  
  # Number of cores to use
  num_cores <- 50 #detectCores() - 1
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  # Progress bar setup
  pb <- progress_bar$new(
    format = " downloading [:bar] :percent eta: :eta",
    total = number_of_runs, clear = FALSE, width = 200
  )
  
  # Initialize coefficient table
  coefficient_table <- data.frame()
  
  # Parallel processing
  compiled_results <- foreach(i = 1:number_of_runs, .combine = rbind, .packages = c('glmnet', 'dplyr', 'caret')) %dopar% {
    set.seed(i)
    
    # Combine stratification variables into a single factor
    stratification_factor <- interaction(meta_training$condition, meta_training$primary_present)
    
    # Create stratified k-folds
    customIndex <- createFolds(stratification_factor, k = nfolds_val, list = TRUE)
    
    # Convert customIndex to a format usable by cv.glmnet
    foldid <- rep(NA, length(stratification_factor))
    for (fold in seq_along(customIndex)) {
      foldid[customIndex[[fold]]] <- fold
    }
    
    # Train the model using cv.glmnet with stratified k-folds
    cvfit <- cv.glmnet(
      x = norm_counts_data_training,
      y = conditions_vector_training,
      family = "binomial",
      type.measure = "auc",
      nfolds = nfolds_val,
      alpha = alpha_value,
      weights = class_weights,
      foldid = foldid,
      standardize = standardize_bit
    )
    
    tmp_coeffs <- coef(cvfit, s = "lambda.min")
    new_results <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)
    
    predictor_names <- rownames(tmp_coeffs)
    coefficients <- as.numeric(tmp_coeffs)
    run_numbers <- rep(i, length(predictor_names))
    data.frame(run_number = run_numbers, predictor = predictor_names, coefficient = coefficients)
  }
  
  # Close the parallel cluster
  stopCluster(cl)
  
  # Update progress bar in the main process
  for (i in 1:number_of_runs) {
    pb$tick()
  }
  
  # Post-processing the results
  coefficient_table <- compiled_results %>%
    filter(coefficient != 0) %>%
    group_by(predictor) %>%
    summarise(
      freq_coef = n(),
      mean_coef = mean(coefficient),
      min_coef = min(coefficient),
      max_coef = max(coefficient)
    )
  
  coefficient_table$abs_mean <- abs(coefficient_table$mean_coef)
  coefficient_table <- coefficient_table %>%
    arrange(desc(freq_coef), desc(abs_mean)) %>%
    filter(predictor != "(Intercept)")
  
  print(coefficient_table, n = 300)
  
  ### START GEMINI ###
  # Lists for ROC data across thresholds
  roc_data_tra <- list()
  roc_data_val <- list()
  
  for (threshold in feat_occ_thresholds) {
    
    # --- FEATURE SELECTION ---
    selected_predictors <- coefficient_table %>%
      filter(freq_coef >= threshold * number_of_runs) %>%
      pull(predictor)
    
    # Filter the dataset based on selected predictors
    norm_counts_data_training_top <- norm_counts_data_training[, selected_predictors]
    norm_counts_data_validation_top <- norm_counts_data_validation[, selected_predictors]
    
    # --- MODEL TRAINING & TUNING ---
    # Combine stratification variables into a single factor
    stratification_factor <- interaction(meta_training$condition, meta_training$primary_present)
    
    # Define a custom index for balanced k-folds based on the combined stratification factor
    customIndex <- createFolds(stratification_factor, k = nfolds_val, returnTrain = TRUE)
    
    cl <- makePSOCKcluster(20)
    registerDoParallel(cl)
    
    # Define training control
    train_control <- trainControl(
      method = "repeatedcv",
      number = nfolds_val,
      repeats = 2,
      summaryFunction = twoClassSummary,
      classProbs = TRUE,
      index = customIndex,
      verboseIter = TRUE
    )
    
    # Train the model using grid search
    set.seed(123)
    glmnet_model <- train(
      x = norm_counts_data_training_top, 
      y = conditions_vector_training_fac, 
      method = "glmnet", 
      trControl = train_control, 
      #tuneGrid = expand.grid(alpha = seq(0, 0.2, by = 0.05), lambda = lambda_values),  # Use alpha value and let glmnet decide lambda
      metric = "ROC",
      weights = class_weights,
      standardize=standardize_bit
    )
    
    # Stop the parallel cluster
    stopCluster(cl)
    
    # View the best parameters
    best_params <- glmnet_model$bestTune
    print(best_params)
    #print(as.matrix(coef(glmnet_model$finalModel, s = best_params$lambda)))
    
    # --- PERFORMANCE EVALUATION ---
    # Ensure to use the correct 'type' argument for probabilities
    probabilities_val <- predict(glmnet_model, newdata = as.matrix(norm_counts_data_validation_top), type = "prob")
    probabilities_tra <- predict(glmnet_model, newdata = as.matrix(norm_counts_data_training_top), type = "prob")
    
    # Ensure that 'probabilities_val' and 'probabilities_tra' are correctly formatted for the 'prediction' function
    pred_val <- prediction(probabilities_val[, 2], conditions_vector_validation)
    pred_tra <- prediction(probabilities_tra[, 2], conditions_vector_training)
    
    # Generate performance metrics
    roc_perf_val <- performance(pred_val, "tpr", "fpr")
    roc_perf_tra <- performance(pred_tra, "tpr", "fpr")
    
    # Store AUC values
    auc_tra <- as.numeric(performance(pred_tra, measure = "auc")@y.values)
    auc_val <- as.numeric(performance(pred_val, measure = "auc")@y.values)
    
    # Confusion matrix and additional metrics
    pred_class_val <- ifelse(probabilities_val[, 2] > 0.5, "X1", "X0")
    cm_val <- confusionMatrix(factor(pred_class_val), conditions_vector_validation_fac)
    sensitivity_val <- cm_val$byClass["Sensitivity"]
    specificity_val <- cm_val$byClass["Specificity"]
    recall_val <- cm_val$byClass["Recall"]
    precision_val <- cm_val$byClass["Precision"]
    accuracy_val <- cm_val$overall["Accuracy"]
    f1_val <- 2 * (precision_val * recall_val) / (precision_val + recall_val)
    
    pred_class_tra <- ifelse(probabilities_tra[, 2] > 0.5, "X1", "X0")
    cm_tra <- confusionMatrix(factor(pred_class_tra), conditions_vector_training_fac)
    sensitivity_tra <- cm_tra$byClass["Sensitivity"]
    specificity_tra <- cm_tra$byClass["Specificity"]
    recall_tra <- cm_tra$byClass["Recall"]
    precision_tra <- cm_tra$byClass["Precision"]
    accuracy_tra <- cm_tra$overall["Accuracy"]
    f1_tra <- 2 * (precision_tra * recall_tra) / (precision_tra + recall_tra)
    
    # --- STORE METRICS AND ROC DATA ---
    metrics_row <- data.frame(
      alpha = alpha_value,
      threshold = threshold,
      num_features = length(selected_predictors),
      auc_tra = auc_tra,
      auc_val = auc_val,
      sensitivity_tra = sensitivity_tra,
      specificity_tra = specificity_tra,
      recall_tra = recall_tra,
      precision_tra = precision_tra,
      accuracy_tra = accuracy_tra,
      f1_tra = f1_tra,
      sensitivity_val = sensitivity_val,
      specificity_val = specificity_val,
      recall_val = recall_val,
      precision_val = precision_val,
      accuracy_val = accuracy_val,
      f1_val = f1_val,
      number_of_runs = number_of_runs,
      #pvalue_cutoff = pvalue_cutoff,
      #log2FoldChange_cutoff = log2FoldChange_cutoff,
      num_top_genes = num_top_genes,
      nfolds_val = nfolds_val,
      root_folder = root_folder,
      sig_results_file = sig_results_file
    )
    
    model_metrics <- rbind(model_metrics, metrics_row)
    
    roc_data_tra[[as.character(threshold)]] <- list(roc_perf_tra, auc_tra)
    roc_data_val[[as.character(threshold)]] <- list(roc_perf_val, auc_val)
    
    # Retain the model with the optimal AUC on the validation set
    if (auc_val > optimal_auc_val) {
      optimal_auc_val <- auc_val
      best_model <- glmnet_model
      optimal_alpha <- alpha_value
      optimal_threshold <- threshold
      best_model_metrics <- metrics_row
      roc_data_tra_best <- roc_data_tra
      roc_data_val_best <- roc_data_val
      # Extract the best parameters and coefficients for the best model
      best_params <- glmnet_model$bestTune
      best_coefficients <- as.matrix(coef(glmnet_model$finalModel, s = best_params$lambda))
    }
  }
}

# --- PLOTTING ROC CURVES FOR THE BEST MODEL ---
plot_roc_curves <- function(roc_data, title) {
  if (length(roc_data) == 0) {
    stop("roc_data is empty. No ROC curves to plot.")
  }
  
  plot(roc_data[[1]][[1]], main = title, col = 1, lty = 1, type = "l") # First curve
  if (length(roc_data) > 1) {
    for (i in 2:length(roc_data)) {
      plot(roc_data[[i]][[1]], add = TRUE, col = i, lty = i, type = "l")
    }
  }
  
  legend_labels <- sapply(names(roc_data), function(thresh) {
    auc_value <- roc_data[[thresh]][[2]]
    paste0("Threshold: ", thresh, " (AUC: ", round(auc_value, 2), ")")
  })
  legend("bottomright", legend = legend_labels, col = 1:length(roc_data), lty = 1:length(roc_data), title = "Threshold", cex = 0.8)
  abline(a = 0, b = 1, lty = 2, col = "gray") # Add 50% AUC line
}

# After the model training and evaluation loop:
# Only plot the ROC curves for the best model

if (length(roc_data_tra_best) > 0) {
  plot_roc_curves(roc_data_tra_best, "Training ROC Curves for Best Model")
} else {
  warning("No training ROC data available for the best model.")
}

if (length(roc_data_val_best) > 0) {
  plot_roc_curves(roc_data_val_best, "Validation ROC Curves for Best Model")
} else {
  warning("No validation ROC data available for the best model.")
}

# Output the best model metrics and parameters
print(paste("Optimal Alpha: ", optimal_alpha))
print(paste("Optimal Threshold: ", optimal_threshold))
print("Best Model Metrics:")
print(best_model_metrics)

print("Best Parameters:")
print(best_params)

#print("Best Coefficients:")
#print(best_coefficients)

# Return the optimal model
#optimal_model

### END GEMINI ###


### Specific Run
if (train_single==TRUE){
  ### TRAIN MANY MODELS
  if (file.exists(paste("./Output/glmnet/",class1_name,class2_name,"_",ver,"/",sep="")) == FALSE ) {
    dir.create(paste("./Output/glmnet/",class1_name,class2_name,"_",ver,"/",sep=""))
  }
  
  pb <- progress_bar$new(
    format = " downloading [:bar] :percent eta: :eta",
    total = number_of_runs, clear = FALSE, width= 200)
  
  coefficient_table = data.frame()
  
  for (i in 1:number_of_runs) {
    pb$tick()
    set.seed(i)
    
    ###
    # Combine stratification variables into a single factor
    stratification_factor <- interaction(meta_training$condition, meta_training$primary_present)
    
    # Create stratified k-folds
    customIndex <- createFolds(stratification_factor, k = nfolds_val, list = TRUE)
    
    # Convert customIndex to a format usable by cv.glmnet
    foldid <- rep(NA, length(stratification_factor))
    for (fold in seq_along(customIndex)) {
      foldid[customIndex[[fold]]] <- fold
    }
    
    # Train the model using cv.glmnet with stratified k-folds
    cvfit <- cv.glmnet(
      x = norm_counts_data_training, 
      y = conditions_vector_training, 
      family = "binomial", 
      type.measure = "auc", 
      nfolds = nfolds_val, 
      alpha = 0.5, 
      weights = class_weights,
      foldid = foldid
    )
    ###
    #cvfit <- cv.glmnet(norm_counts_data_training, conditions_vector_training, family = "binomial", type.measure="auc", nfolds=nfolds_val, alpha=0.95, weights = class_weights)
    tmp_coeffs <- coef(cvfit, s = "lambda.min")
    
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
    filter(freq_coef >= 0.95 * number_of_runs) %>%
    pull(predictor)
  
  print(selected_predictors)
  
  norm_counts_data_training_top <- norm_counts_data_training[, selected_predictors]
  norm_counts_data_validation_top <- norm_counts_data_validation[, selected_predictors]
  norm_counts_data_holdout_top <- norm_counts_data_holdout[, selected_predictors]
  
  
  # Define the parameter grid
  tune_grid <- expand.grid(
    alpha = seq(0, 1, by = 0.05),
    lambda_values
  )
  
  
  # Combine stratification variables into a single factor
  stratification_factor <- interaction(meta_training$condition) #meta_training$primary_site
  
  # Define a custom index for balanced k-folds based on the combined stratification factor
  customIndex <- createFolds(stratification_factor, k = nfolds_val, returnTrain = TRUE)
  
  # Define training control
  train_control <- trainControl(method = "repeatedcv",
                                number=nfolds_val,
                                repeats = 10,
                                summaryFunction = twoClassSummary,
                                ## Estimate class probabilities
                                classProbs = TRUE,
                                index = customIndex,
                                verboseIter = TRUE)
  
  # Train the model using grid search
  set.seed(123)
  glmnet_model <- train(
    x = norm_counts_data_training_top, 
    y = conditions_vector_training_fac, 
    method = "glmnet", 
    trControl = train_control, 
    tuneGrid = tune_grid, 
    metric = "ROC",
    weights = class_weights
  )
  
  # View the best parameters
  best_params <- glmnet_model$bestTune
  print(best_params)
  
  # Use the final model retrained on the entire training data with the best parameters
  best_model <- glmnet_model$finalModel
  
  # Extract coefficients using the optimal lambda value
  best_coefficients <- coef(best_model, s = best_params$lambda)
}

# Convert the coefficients to a matrix
coefficients_matrix <- as.matrix(best_coefficients)

# Filter non-zero coefficients
non_zero_coefficients <- coefficients_matrix#[coefficients_matrix != 0, , drop = FALSE]

# Convert to a readable format and print
non_zero_coefficients_list <- as.data.frame(non_zero_coefficients)

selected_predictors <- rownames(non_zero_coefficients)[-1]

norm_counts_data_training_top <- norm_counts_data_training[, selected_predictors]
norm_counts_data_validation_top <- norm_counts_data_validation[, selected_predictors]
norm_counts_data_holdout_top <- norm_counts_data_holdout[, selected_predictors]

# Compute the probabilities for the validation, holdout, and training sets
probabilities_val <- predict(best_model$finalModel, newx = as.matrix(norm_counts_data_validation_top), s = best_params$lambda, type = "response")
probabilities_hol <- predict(best_model$finalModel, newx = as.matrix(norm_counts_data_holdout_top), s = best_params$lambda, type = "response")
probabilities_tra <- predict(best_model$finalModel, newx = as.matrix(norm_counts_data_training_top), s = best_params$lambda, type = "response")

# Create prediction objects with ROCR
pred_val <- prediction(probabilities_val, conditions_vector_validation)
pred_hol <- prediction(probabilities_hol, conditions_vector_holdout)
pred_tra <- prediction(probabilities_tra, conditions_vector_training)

# Generate performance metrics
roc_perf_val <- performance(pred_val, "tpr", "fpr")
roc_perf_hol <- performance(pred_hol, "tpr", "fpr")
roc_perf_tra <- performance(pred_tra, "tpr", "fpr")

# Calculate AUC for training set
auc_obj_tra <- performance(pred_tra, measure = "auc")
auc_tra <- as.numeric(auc_obj_tra@y.values)
print(paste("Training AUC: ", auc_tra))

# Calculate AUC for validation set
auc_obj_val <- performance(pred_val, measure = "auc")
auc_val <- as.numeric(auc_obj_val@y.values)
print(paste("Validation AUC: ", auc_val))

# Calculate AUC for holdout set
auc_obj_hol <- performance(pred_hol, measure = "auc")
auc_hol <- as.numeric(auc_obj_hol@y.values)
print(paste("Holdout AUC: ", auc_hol))

# Calculate sensitivity and specificity at various thresholds
perf_sens <- performance(pred_val, measure = "sens")
perf_spec <- performance(pred_val, measure = "spec")

# Extract sensitivity and specificity values
sens_values <- unlist(perf_sens@y.values)
spec_values <- unlist(perf_spec@y.values)
thresholds <- unlist(perf_sens@x.values)

# Calculate Youden's J statistic
youden_j <- sens_values + spec_values - 1

# Find the threshold that maximizes Youden's J
optimal_index <- which.max(youden_j)
opt_threshold <- thresholds[optimal_index]
print(paste("Optimal threshold (Youden's J):", opt_threshold))

# Convert probabilities to class labels using the optimal threshold
predicted_labels_tra <- ifelse(probabilities_tra >= opt_threshold, 1, 0)
predicted_labels_val <- ifelse(probabilities_val >= opt_threshold, 1, 0)
predicted_labels_hol <- ifelse(probabilities_hol >= opt_threshold, 1, 0)

# Calculate performance metrics for training set
confusionMatrix(factor(predicted_labels_tra), factor(conditions_vector_training))

# Calculate performance metrics for validation set
confusionMatrix(factor(predicted_labels_val), factor(conditions_vector_validation))

# Calculate performance metrics for holdout set
confusionMatrix(factor(predicted_labels_hol), factor(conditions_vector_holdout))

# Plot the ROC curve for training set
plot_tra <- plot(roc_perf_tra, colorize = TRUE, print.cutoffs.at = seq(0,1,by=0.1), text.adj = c(-0.2,1.7),print.auc = TRUE)
# Add a diagonal line representing 50% AUC
abline(a = 0, b = 1, lty = 2, col = "gray")

# Plot the ROC curve for validation set
plot_val <- plot(roc_perf_val, colorize = TRUE, print.cutoffs.at = seq(0,1,by=0.1), text.adj = c(-0.2,1.7),print.auc = TRUE)
# Add a diagonal line representing 50% AUC
abline(a = 0, b = 1, lty = 2, col = "gray")

# Plot the ROC curve for holdout set
plot_hol <- plot(roc_perf_hol, colorize = TRUE, print.cutoffs.at = seq(0,1,by=0.1), text.adj = c(-0.2,1.7),print.auc = TRUE)
# Add a diagonal line representing 50% AUC
abline(a = 0, b = 1, lty = 2, col = "gray")

### ANNOTATE FEATURES WITH SYMBOLS AND GENOMIC FEATURE TYPES

# Function to annotate features
annotate_features <- function(features) {
  # Split the features into chromosome, start, and end positions
  chr_positions <- strsplit(features, "_")
  chr <- sapply(chr_positions, function(x) x[1])
  start <- as.numeric(sapply(chr_positions, function(x) x[2]))
  end <- as.numeric(sapply(chr_positions, function(x) x[3]))
  
  # Create a GRanges object
  gr <- GRanges(seqnames = chr, ranges = IRanges(start = start, end = end))
  
  # Load the TxDb annotation database
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  
  # Annotate the peaks with gene information
  peakAnno <- annotatePeak(gr, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db", overlap = "all")
  
  # Extract the gene symbols and location types
  gene_symbols <- as.character(peakAnno@anno$SYMBOL)
  location_types <- as.character(peakAnno@anno$annotation)
  
  # Replace missing gene symbols with "N/A"
  gene_symbols[is.na(gene_symbols)] <- "N/A"
  
  return(data.frame(Feature = features, GeneSymbol = gene_symbols, LocationType = location_types))
}

# Drop the intercept row by name
coeff_matrix <- non_zero_coefficients[rownames(non_zero_coefficients) != "(Intercept)", , drop = FALSE]

# Convert back to sparse matrix if needed
final_model_coeffs <- as(coeff_matrix, "dgCMatrix")

selected_features_final_model <- rownames(final_model_coeffs)
selected_features_final_model_coeffs <- final_model_coeffs[, 1]

# Create a dataframe for the features and their coefficients
features_coeffs_df <- data.frame(Feature = selected_features_final_model, Coefficient = selected_features_final_model_coeffs)

# Annotate the selected features
annotated_features_df <- annotate_features(selected_features_final_model)

# Combine the coefficients with the annotations
final_features_annotated_df <- merge(features_coeffs_df, annotated_features_df, by = "Feature")

# Print the final annotated features with coefficients
print("Final annotated features with coefficients:")
print(final_features_annotated_df)

# HEATMAP FUNCTION
generate_heatmaps <- function(training_data, validation_data, holdout_data, annotated_features, meta_training, meta_validation, meta_holdout) {
  counts_data_training_top_t <- t(scale(training_data, scale = TRUE, center = TRUE))
  rownames(counts_data_training_top_t) <- annotated_features$GeneSymbol
  
  counts_data_validation_top_t <- t(scale(validation_data, scale = TRUE, center = TRUE))
  rownames(counts_data_validation_top_t) <- annotated_features$GeneSymbol
  
  counts_data_holdout_top_t <- t(scale(holdout_data, scale = TRUE, center = TRUE))
  rownames(counts_data_holdout_top_t) <- annotated_features$GeneSymbol
  
  groups <- unique(meta_training$condition)
  contrast_groups <- c("condition", groups[1], groups[2])
  group1 <- paste(contrast_groups[[3]])
  group2 <- paste(contrast_groups[[2]])
  
  # Extract relevant metadata columns, including batch
  annotation_training <- meta_training %>%
    dplyr::select(condition, primary_site, primary_present, peritoneal_mets, batch)
  
  annotation_validation <- meta_validation %>%
    dplyr::select(condition, primary_site, primary_present, peritoneal_mets, batch)
  
  annotation_holdout <- meta_holdout %>%
    dplyr::select(condition, primary_site, primary_present, peritoneal_mets, batch)
  
  # Define color schemes for annotation
  condition_table <- c("#E31A1C", "#1F78B4")
  names(condition_table) <- c(group2, group1)
  
  primary_site_colors <- c("#33A02C", "#6A3D9A", "#FF7D9A")
  primary_site_table <- setNames(primary_site_colors, c("CRC", "HGA", "LGA"))
  
  primary_present_colors <- c("#FF7F00", "#FFFF99", "#CCCCCC")
  primary_present_table <- setNames(primary_present_colors, c("primary_absent", "primary_present", "healthy"))
  
  peritoneal_mets_colors <- c("#D985B2", "#00B3B3", "#CCCCCC")
  peritoneal_mets_table <- setNames(peritoneal_mets_colors, c("pm_absent", "pm_present", "healthy"))
  
  # Define color schemes for batches
  batch_table <- c("#D3D3D3", "#696969") 
  names(batch_table) <- unique(annotation_training$batch)
  
  # Create heatmap annotation
  ha_training <- HeatmapAnnotation(
    condition = annotation_training$condition,
    primary_site = annotation_training$primary_site,
    primary_present = annotation_training$primary_present,
    peritoneal_mets = annotation_training$peritoneal_mets,
    batch = annotation_training$batch,
    col = list(
      condition = condition_table,
      primary_site = primary_site_table,
      primary_present = primary_present_table,
      peritoneal_mets = peritoneal_mets_table,
      batch = batch_table
    ),
    show_annotation_name = FALSE
  )
  
  ha_validation <- HeatmapAnnotation(
    condition = annotation_validation$condition,
    primary_site = annotation_validation$primary_site,
    primary_present = annotation_validation$primary_present,
    peritoneal_mets = annotation_validation$peritoneal_mets,
    batch = annotation_validation$batch,
    col = list(
      condition = condition_table,
      primary_site = primary_site_table,
      primary_present = primary_present_table,
      peritoneal_mets = peritoneal_mets_table,
      batch = batch_table
    ),
    show_annotation_name = FALSE
  )
  
  ha_holdout <- HeatmapAnnotation(
    condition = annotation_holdout$condition,
    primary_site = annotation_holdout$primary_site,
    primary_present = annotation_holdout$primary_present,
    peritoneal_mets = annotation_holdout$peritoneal_mets,
    batch = annotation_holdout$batch,
    col = list(
      condition = condition_table,
      primary_site = primary_site_table,
      primary_present = primary_present_table,
      peritoneal_mets = peritoneal_mets_table,
      batch = batch_table
    )
  )
  
  # Create heatmaps
  heatmap_tra <- Heatmap(counts_data_training_top_t, 
                         top_annotation = ha_training,
                         col = colorRamp2(c(-4, -2, -1, 0, 1, 2, 4), c("#4676b5", "#82b1d3", "#dbeff6", "#fefebd", "#fee395", "#fc9961", "#d73027")),
                         show_row_names = TRUE,
                         show_column_names = FALSE,
                         clustering_distance_columns = "euclidean")
  
  heatmap_val <- Heatmap(counts_data_validation_top_t, 
                         top_annotation = ha_validation,
                         col = colorRamp2(c(-4, -2, -1, 0, 1, 2, 4), c("#4676b5", "#82b1d3", "#dbeff6", "#fefebd", "#fee395", "#fc9961", "#d73027")),
                         show_row_names = TRUE,
                         show_column_names = FALSE,
                         clustering_distance_columns = "euclidean")
  
  heatmap_hol <- Heatmap(counts_data_holdout_top_t, 
                         top_annotation = ha_holdout,
                         col = colorRamp2(c(-4, -2, -1, 0, 1, 2, 4), c("#4676b5", "#82b1d3", "#dbeff6", "#fefebd", "#fee395", "#fc9961", "#d73027")),
                         show_row_names = TRUE,
                         show_column_names = FALSE,
                         clustering_distance_columns = "euclidean")
  
  return(heatmap_tra + heatmap_val + heatmap_hol)
}

### GENERATE HEATMAPS
heatmap <- generate_heatmaps(norm_counts_data_training_top, norm_counts_data_validation_top, norm_counts_data_holdout_top, annotated_features_df, meta_training, meta_validation, meta_holdout)
png(paste("~/5hmC-Pathway-Analysis/sig_heatmap_120.png", sep = ""), width = 1200, height = 900)
heatmap
dev.off()

