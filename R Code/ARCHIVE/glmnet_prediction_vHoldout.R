### LOAD LIBRARIES
library(glmnet)
library(dplyr)
library(ROCR)
library(svMisc)
library(glmnetUtils)
library(progress)
library(caret)
library(pheatmap)
library(ComplexHeatmap)
library(pROC)
library(ChIPseeker)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

### CONFIGURATION
setwd("~/5hmC-Pathway-Analysis/")

## File paths
sig_results_file <- "./Output/DESeq2/Results/all_results_PMnegPMpos_pval 0.001 _lfc 0.13 _ CRCHGA_combati_602020_052924_cp_seed122.txt"
root_folder <- "./Output/Randomization/PMneg_PMpos_DESeq2_hmr_combat_602020_052924_CRCHGA_cp_seed122/"
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

# Parameters
sig_type <- "DESEQ"
ver.text <- "CRCHGA_hmr_combati_602020_cp_seed122"
cutoff_type <- 1
number_of_runs <- 40
padj_cutoff <- 0.001
pvalue_cutoff <- 0.01
log2FoldChange_cutoff <- 0.13
nfolds_val <- 10
seed_val <- 7
alpha_values <- seq(0, 1, by = 0.1)  # Define alpha values
feature_occurrence_thresholds <- c(0.8)#0.75,0.8,0.85,0.9,0.95,0.98)
skip_rfe <- TRUE  # Set to TRUE to skip RFE, FALSE to apply RFE

### SET VERSION STRING BASED ON CUTOFF TYPE
if(cutoff_type == 0){
  ver <- paste("padj", padj_cutoff, "_lfc", log2FoldChange_cutoff, "_nfold", nfolds_val, "+", ver.text)
} else {
  ver <- paste("pval", pvalue_cutoff, "_lfc", log2FoldChange_cutoff, "_nfold", nfolds_val, "+", ver.text)
}

### READ AND PREPROCESS DATA

# Load and transpose data
norm_counts_data_training <- read.csv(norm_counts_name_training, row.names = 1)
norm_counts_data_validation <- read.csv(norm_counts_name_validation, row.names = 1)
norm_counts_data_holdout <- read.csv(norm_counts_name_holdout, row.names = 1)

norm_counts_data_training <- t(as.matrix(norm_counts_data_training))
norm_counts_data_validation <- t(as.matrix(norm_counts_data_validation))
norm_counts_data_holdout <- t(as.matrix(norm_counts_data_holdout))

# Load metadata
meta_training <- read.csv(meta_name_training, row.names = 1)
meta_validation <- read.csv(meta_name_validation, row.names = 1)
meta_holdout <- read.csv(meta_name_holdout, row.names = 1)

### APPLY MAPPING TO META DATA

# Define mapping
map <- list("pilot" = 0, "large_cohort" = 1, "female" = 0, "male" = 1, "other_mets_absent" = 0, "other_mets_present" = 1, "Yes" = 1, "No" = 0, "LGA_PMneg" = "PMneg", "LGA_PMpos" = "PMpos", "HGA_PMneg" = "PMneg", "HGA_PMpos" = "PMpos")

# Function to apply mapping
apply_mapping <- function(meta_data) {
  meta_data[] <- lapply(meta_data, function(x) {
    x <- as.character(x)
    for (i in 1:length(map)) {
      x[x == names(map)[i]] <- map[[i]]
    }
    return(x)
  })
  return(meta_data)
}

# Apply mapping to all metadata
meta_training <- apply_mapping(meta_training)
meta_validation <- apply_mapping(meta_validation)
meta_holdout <- apply_mapping(meta_holdout)

### FILTER SIGNIFICANT GENES

if (sig_type == "LOGIT") {
  sig_results <- read.table(sig_results_file, header = TRUE, sep = ",")
  sig_results <- sig_results[sig_results[,"p.value"] < pvalue_cutoff,]
  sig_results <- na.omit(sig_results)
  sig_results_genes <- sig_results$term
} else {
  sig_results <- read.table(sig_results_file, header = TRUE)
  sig_results <- subset(sig_results, pvalue < pvalue_cutoff & abs(log2FoldChange) > log2FoldChange_cutoff)
  sig_results <- na.omit(sig_results)
  sig_results_genes <- rownames(sig_results)
}

# Get the intersection of significant genes and column names of count data
existing_sig_genes <- intersect(sig_results_genes, colnames(norm_counts_data_training))
existing_sig_genes <- intersect(existing_sig_genes, colnames(norm_counts_data_validation))
existing_sig_genes <- intersect(existing_sig_genes, colnames(norm_counts_data_holdout))

# Subset count data to select only the existing significant genes
norm_counts_data_training <- norm_counts_data_training[, existing_sig_genes]
norm_counts_data_validation <- norm_counts_data_validation[, existing_sig_genes]
norm_counts_data_holdout <- norm_counts_data_holdout[, existing_sig_genes]

### PREPARE CONDITION VECTORS AND WEIGHTS

class1_name <- meta_training[2, 1]
class2_name <- meta_training[nrow(meta_training), 1]

# Function to create condition vectors
conditions_vector <- function(meta_data, class1_name, class2_name) {
  class1_count <- sum(meta_data$condition == class1_name)
  class2_count <- sum(meta_data$condition == class2_name)
  return(c(rep(1, class1_count), rep(0, class2_count)))
}

conditions_vector_training <- conditions_vector(meta_training, class1_name, class2_name)
conditions_vector_validation <- conditions_vector(meta_validation, class1_name, class2_name)
conditions_vector_holdout <- conditions_vector(meta_holdout, class1_name, class2_name)

# Convert class levels to valid R variable names
conditions_vector_training_fac <- as.factor(make.names(conditions_vector_training))
conditions_vector_validation_fac <- as.factor(make.names(conditions_vector_validation))
conditions_vector_holdout_fac <- as.factor(make.names(conditions_vector_holdout))

# Calculate class weights
class_weights <- ifelse(conditions_vector_training == 1, 1 / sum(conditions_vector_training == 1), 1 / sum(conditions_vector_training == 0))
class_weights_val <- ifelse(conditions_vector_validation == 1, 1 / sum(conditions_vector_validation == 1), 1 / sum(conditions_vector_validation == 0))


### OPTIMIZATION FUNCTIONS

# Function to optimize hyperparameters using cross-validation on training set
optimize_hyperparameters <- function(training_data, validation_data, training_labels, validation_labels, class_weights) {
  tune_grid <- expand.grid(alpha = seq(0, 1, by = 0.1), lambda = 10^seq(-3, 0, length = 30))
  train_control <- trainControl(method = "cv", number = 10, summaryFunction = twoClassSummary, classProbs = TRUE, verboseIter = TRUE)
  glmnet_model <- train(x = training_data, y = training_labels, method = "glmnet", trControl = train_control, tuneGrid = tune_grid, metric = "ROC", weights = class_weights)
  return(glmnet_model$bestTune)
}


# Function to run multiple models and collect coefficients
run_models_collect_coefficients <- function(data, labels, class_weights, alpha_values, n_runs = 100, nfolds_val = 10) {
  pb <- progress_bar$new(total = length(alpha_values) * n_runs, format = " running [:bar] :percent eta: :eta")
  coefficient_list <- list()
  for (alpha in alpha_values) {
    for (i in 1:n_runs) {
      pb$tick()
      set.seed(i)
      cvfit <- cv.glmnet(data, labels, family = "binomial", type.measure = "class", nfolds = nfolds_val, alpha = alpha, weights = class_weights)
      coeffs <- as.matrix(coef(cvfit, s = "lambda.1se"))
      coefficient_list[[length(coefficient_list) + 1]] <- coeffs
    }
  }
  return(coefficient_list)
}

# Function to aggregate coefficients and find frequently selected features
aggregate_coefficients <- function(coefficient_list, threshold) {
  all_features <- unique(unlist(lapply(coefficient_list, function(x) rownames(x)[x != 0])))
  feature_freq <- sapply(all_features, function(f) mean(sapply(coefficient_list, function(x) any(rownames(x) == f & x[rownames(x) == f, ] != 0))))
  selected_features <- names(feature_freq)[feature_freq >= threshold]
  return(selected_features)
}

# Function to ensure selected features are in the training data
ensure_features_in_data <- function(features, data) {
  return(features[features %in% colnames(data)])
}

# Function to optimize threshold using Youden's J statistic
optimize_threshold <- function(model, data, labels) {
  probabilities <- predict(model, newx = data, s = "lambda.min", type = "response")
  pred <- prediction(probabilities, labels)
  perf <- performance(pred, "sens", "spec")
  J <- perf@y.values[[1]] + perf@x.values[[1]] - 1
  opt_index <- which.max(J)
  optimal_threshold <- perf@alpha.values[[1]][opt_index]
  return(optimal_threshold)
}

### FEATURE SELECTION AND MODEL TRAINING

# Define a range of feature occurrence thresholds to optimize
best_auc <- 0
best_threshold <- NULL
best_features <- NULL

# Iterate over feature occurrence thresholds to find the best model
for (threshold in feature_occurrence_thresholds) {
  coefficient_list <- run_models_collect_coefficients(norm_counts_data_training, conditions_vector_training_fac, class_weights, alpha_values, number_of_runs, nfolds_val)
  selected_features <- aggregate_coefficients(coefficient_list, threshold)
  
  # Ensure selected features are in the training and validation data
  selected_features <- ensure_features_in_data(selected_features, norm_counts_data_training)
  selected_features <- ensure_features_in_data(selected_features, norm_counts_data_validation)
  
  if (length(selected_features) > 0) {
    norm_counts_data_training_top <- norm_counts_data_training[, selected_features]
    norm_counts_data_validation_top <- norm_counts_data_validation[, selected_features]
    
    # Optimize hyperparameters using the training set and validate on the validation set
    best_params <- optimize_hyperparameters(norm_counts_data_training_top, norm_counts_data_validation_top, conditions_vector_training_fac, conditions_vector_validation_fac, class_weights)
    
    # Train final model on the training set with optimized hyperparameters
    final_model <- glmnet(x = norm_counts_data_training_top, y = conditions_vector_training_fac, family = "binomial", alpha = best_params$alpha, lambda = best_params$lambda, weights = class_weights)
    
    # Predict on the validation set
    probabilities_validation <- predict(final_model, newx = norm_counts_data_validation_top, s = "lambda.min", type = "response")
    
    # Compute AUC on the validation set
    roc_obj <- roc(conditions_vector_validation, as.numeric(probabilities_validation))
    auc <- auc(roc_obj)
    
    if (auc > best_auc) {
      best_auc <- auc
      best_threshold <- threshold
      best_features <- selected_features
    }
  }
}



### FINAL MODEL TRAINING AND EVALUATION

print("Best feature occurrence threshold:")
print(best_threshold)
print("Selected features for final model:")
print(best_features)

# Use the best threshold and features to train the final model
norm_counts_data_training_top <- norm_counts_data_training[, best_features]
norm_counts_data_validation_top <- norm_counts_data_validation[, best_features]
norm_counts_data_holdout_top <- norm_counts_data_holdout[, best_features]

best_params <- optimize_hyperparameters(norm_counts_data_training_top, norm_counts_data_validation_top, conditions_vector_training_fac, conditions_vector_validation_fac, class_weights)

final_model <- glmnet(x = norm_counts_data_training_top, y = conditions_vector_training_fac, family = "binomial", alpha = best_params$alpha, lambda = best_params$lambda, weights = class_weights)

# Extract and print features in the final model
final_model_coeffs <- coef(final_model, s = "lambda.min")
selected_features_final_model <- rownames(final_model_coeffs)[final_model_coeffs[, 1] != 0]
print("Selected features in the final model:")
print(selected_features_final_model)

# Optimize threshold on the validation set
optimal_threshold <- optimize_threshold(final_model, norm_counts_data_validation_top, conditions_vector_validation_fac)
print("Optimal threshold:")
print(optimal_threshold)

# Predict probabilities with the final model
probabilities_training <- predict(final_model, newx = norm_counts_data_training_top, s = "lambda.min", type = "response")
probabilities_validation <- predict(final_model, newx = norm_counts_data_validation_top, s = "lambda.min", type = "response")
probabilities_holdout <- predict(final_model, newx = norm_counts_data_holdout_top, s = "lambda.min", type = "response")

# Predict with the optimized threshold
predict_with_threshold <- function(probabilities, threshold) {
  return(ifelse(probabilities >= threshold, 1, 0))
}

predicted_labels_training <- predict_with_threshold(probabilities_training, optimal_threshold)
predicted_labels_validation <- predict_with_threshold(probabilities_validation, optimal_threshold)
predicted_labels_holdout <- predict_with_threshold(probabilities_holdout, optimal_threshold)

# Compute metrics
compute_metrics <- function(true_labels, predicted_labels) {
  confusion <- confusionMatrix(factor(predicted_labels), factor(true_labels))
  roc_obj <- roc(true_labels, as.numeric(predicted_labels))
  auc <- auc(roc_obj)
  return(list(confusion = confusion, auc = auc))
}

metrics_training <- compute_metrics(conditions_vector_training, predicted_labels_training)
metrics_validation <- compute_metrics(conditions_vector_validation, predicted_labels_validation)
metrics_holdout <- compute_metrics(conditions_vector_holdout, predicted_labels_holdout)

# Print final metrics
print(paste("Training AUC: ", metrics_training$auc))
print(paste("Validation AUC: ", metrics_validation$auc))
print(paste("Holdout AUC: ", metrics_holdout$auc))
print(metrics_training$confusion)
print(metrics_validation$confusion)
print(metrics_holdout$confusion)

# Print optimized and hard-coded key variables
cat("\nOptimized and Hard-coded Key Variables:\n")
cat("Number of runs:", number_of_runs, "\n")
cat("Alpha values tested:", alpha_values, "\n")
cat("Best alpha:", best_params$alpha, "\n")
cat("Best lambda:", best_params$lambda, "\n")
cat("Best feature occurrence threshold:", best_threshold, "\n")
cat("Optimal threshold:", optimal_threshold, "\n")
cat("Number of selected features:", length(selected_features_final_model), "\n")

### APPLY RFE
if (!skip_rfe) {
  apply_rfe <- function(training_data, validation_data, training_labels, validation_labels, num_features) {
    control <- rfeControl(functions = rfFuncs, method = "cv", number = 10, saveDetails = TRUE)
    rfe_result <- rfe(x = training_data, y = training_labels, sizes = seq(1, num_features, by = 1), rfeControl = control)
    
    # Validate RFE-selected features on the validation set
    selected_features <- predictors(rfe_result)
    training_data_rfe <- training_data[, selected_features]
    validation_data_rfe <- validation_data[, selected_features]
    
    # Optimize hyperparameters again on the RFE-selected features
    best_params_rfe <- optimize_hyperparameters(training_data_rfe, validation_data_rfe, training_labels, validation_labels, class_weights)
    
    return(list(selected_features = selected_features, best_params_rfe = best_params_rfe))
  }
  
  ### FINAL FEATURE SELECTION AND MODEL TRAINING WITH RFE
  
  # Apply RFE to ensure sparsest model
  rfe_results <- apply_rfe(norm_counts_data_training_top, norm_counts_data_validation_top, conditions_vector_training_fac, conditions_vector_validation_fac, length(selected_features_final_model))
  selected_features <- rfe_results$selected_features
  best_params <- rfe_results$best_params_rfe
  
  # Use the selected features from RFE
  norm_counts_data_training_top <- norm_counts_data_training_top[, selected_features]
  norm_counts_data_validation_top <- norm_counts_data_validation_top[, selected_features]
  norm_counts_data_holdout_top <- norm_counts_data_holdout_top[, selected_features]
  
  # Train the final model with RFE-selected features
  final_model <- glmnet(x = norm_counts_data_training_top, y = conditions_vector_training_fac, family = "binomial", alpha = best_params$alpha, lambda = best_params$lambda, weights = class_weights)
  
} else {
  # Use the best threshold and features to train the final model
  norm_counts_data_training_top <- norm_counts_data_training[, best_features]
  norm_counts_data_validation_top <- norm_counts_data_validation[, best_features]
  norm_counts_data_holdout_top <- norm_counts_data_holdout[, best_features]
  
  best_params <- optimize_hyperparameters(norm_counts_data_training_top, norm_counts_data_validation_top, conditions_vector_training_fac, conditions_vector_validation_fac, class_weights)
  
  final_model <- glmnet(x = norm_counts_data_training_top, y = conditions_vector_training_fac, family = "binomial", alpha = best_params$alpha, lambda = best_params$lambda, weights = class_weights)
}

# Extract and print features in the final model
final_model_coeffs <- coef(final_model, s = "lambda.min")
selected_features_final_model <- rownames(final_model_coeffs)[final_model_coeffs[, 1] != 0]
print("Selected features in the final model:")
print(selected_features_final_model)

# Optimize threshold on the validation set
optimal_threshold <- optimize_threshold(final_model, norm_counts_data_validation_top, conditions_vector_validation_fac)
print("Optimal threshold:")
print(optimal_threshold)

# Predict probabilities with the final model
probabilities_training <- predict(final_model, newx = norm_counts_data_training_top, s = "lambda.min", type = "response")
probabilities_validation <- predict(final_model, newx = norm_counts_data_validation_top, s = "lambda.min", type = "response")
probabilities_holdout <- predict(final_model, newx = norm_counts_data_holdout_top, s = "lambda.min", type = "response")

# Predict with the optimized threshold
predict_with_threshold <- function(probabilities, threshold) {
  return(ifelse(probabilities >= threshold, 1, 0))
}

predicted_labels_training <- predict_with_threshold(probabilities_training, optimal_threshold)
predicted_labels_validation <- predict_with_threshold(probabilities_validation, optimal_threshold)
predicted_labels_holdout <- predict_with_threshold(probabilities_holdout, optimal_threshold)

# Compute metrics
compute_metrics <- function(true_labels, predicted_labels) {
  confusion <- confusionMatrix(factor(predicted_labels), factor(true_labels))
  roc_obj <- roc(true_labels, as.numeric(predicted_labels))
  auc <- auc(roc_obj)
  return(list(confusion = confusion, auc = auc))
}

metrics_training <- compute_metrics(conditions_vector_training, predicted_labels_training)
metrics_validation <- compute_metrics(conditions_vector_validation, predicted_labels_validation)
metrics_holdout <- compute_metrics(conditions_vector_holdout, predicted_labels_holdout)

# Print final metrics
print(paste("Training AUC: ", metrics_training$auc))
print(paste("Validation AUC: ", metrics_validation$auc))
print(paste("Holdout AUC: ", metrics_holdout$auc))
print(metrics_training$confusion)
print(metrics_validation$confusion)
print(metrics_holdout$confusion)

# Print optimized and hard-coded key variables
cat("\nOptimized and Hard-coded Key Variables:\n")
cat("Number of runs:", number_of_runs, "\n")
cat("Alpha values tested:", alpha_values, "\n")
cat("Best alpha:", best_params$alpha, "\n")
cat("Best lambda:", best_params$lambda, "\n")
cat("Best feature occurrence threshold:", best_threshold, "\n")
cat("Optimal threshold:", optimal_threshold, "\n")
cat("Number of selected features:", length(selected_features_final_model), "\n")

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


# Convert to a regular matrix to manipulate
coeff_matrix <- as.matrix(final_model_coeffs)

# Drop the intercept row by name
coeff_matrix <- coeff_matrix[rownames(coeff_matrix) != "(Intercept)", , drop = FALSE]

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

### HEATMAP FUNCTION
generate_heatmaps <- function(training_data, validation_data, holdout_data, annotated_features, meta_training, meta_validation, meta_holdout) {
  counts_data_training_top_t <- t(scale(training_data, scale = TRUE, center = TRUE))
  rownames(counts_data_training_top_t) <- annotated_features$GeneSymbol
  
  counts_data_validation_top_t <- t(scale(validation_data, scale = TRUE, center = TRUE))
  rownames(counts_data_validation_top_t) <- annotated_features$GeneSymbol
  
  counts_data_holdout_top_t <- t(scale(holdout_data, scale = TRUE, center = TRUE))
  rownames(counts_data_holdout_top_t) <- annotated_features$GeneSymbol
  
  groups <- unique(meta_training[c("condition")])
  contrast_groups <- c("condition", groups[1, 1], groups[2, 1])
  group1 <- paste(contrast_groups[[3]])
  group2 <- paste(contrast_groups[[2]])
  
  # Extract relevant metadata columns
  annotation_training <- meta_training %>%
    dplyr::select(condition, primary_site, primary_present, peritoneal_mets)
  
  # Only include annotation on the holdout
  annotation_holdout <- meta_holdout %>%
    dplyr::select(condition, primary_site, primary_present, peritoneal_mets)
  
  annotation_validation <- meta_validation %>%
    dplyr::select(condition, primary_site, primary_present, peritoneal_mets)
  
  # Define color schemes for annotation
  condition_table <- c("#E31A1C", "#1F78B4")
  names(condition_table) <- c(group2, group1)
  
  primary_site_colors <- c("#33A02C", "#6A3D9A", "#FF7D9A")
  primary_site_table <- setNames(primary_site_colors, c("CRC", "HGA", "LGA"))
  
  primary_present_colors <- c("#FF7F00", "#FFFF99", "#CCCCCC")
  primary_present_table <- setNames(primary_present_colors, c("primary_absent", "primary_present", "healthy"))
  
  peritoneal_mets_colors <- c("#D985B2", "#00B3B3", "#CCCCCC")
  peritoneal_mets_table <- setNames(peritoneal_mets_colors, c("pm_absent", "pm_present", "healthy"))
  
  # Create heatmap annotation
  ha_training <- HeatmapAnnotation(
    condition = annotation_training$condition,
    primary_site = annotation_training$primary_site,
    primary_present = annotation_training$primary_present,
    peritoneal_mets = annotation_training$peritoneal_mets,
    col = list(
      condition = condition_table,
      primary_site = primary_site_table,
      primary_present = primary_present_table,
      peritoneal_mets = peritoneal_mets_table
    ),
    show_annotation_name = FALSE
  )
  
  ha_validation <- HeatmapAnnotation(
    condition = annotation_validation$condition,
    primary_site = annotation_validation$primary_site,
    primary_present = annotation_validation$primary_present,
    peritoneal_mets = annotation_validation$peritoneal_mets,
    col = list(
      condition = condition_table,
      primary_site = primary_site_table,
      primary_present = primary_present_table,
      peritoneal_mets = peritoneal_mets_table
    ),
    show_annotation_name = FALSE
  )
  
  # Create holdout heatmap annotation
  ha_holdout <- HeatmapAnnotation(
    condition = annotation_holdout$condition,
    primary_site = annotation_holdout$primary_site,
    primary_present = annotation_holdout$primary_present,
    peritoneal_mets = annotation_holdout$peritoneal_mets,
    col = list(
      condition = condition_table,
      primary_site = primary_site_table,
      primary_present = primary_present_table,
      peritoneal_mets = peritoneal_mets_table
    )
  )
  
  # Create heatmaps
  heatmap_tra <- Heatmap(counts_data_training_top_t, 
                         top_annotation = ha_training,
                         col = colorRamp2(c(-4, -2, -1, 0, 1, 2, 4), c("#4676b5", "#82b1d3", "#dbeff6", "#fefebd", "#fee395", "#fc9961", "#d73027")),
                         show_row_names = TRUE,
                         show_column_names = FALSE)
  
  heatmap_val <- Heatmap(counts_data_validation_top_t, 
                         top_annotation = ha_validation,
                         col = colorRamp2(c(-4, -2, -1, 0, 1, 2, 4), c("#4676b5", "#82b1d3", "#dbeff6", "#fefebd", "#fee395", "#fc9961", "#d73027")),
                         show_row_names = TRUE,
                         show_column_names = FALSE)
  
  heatmap_hol <- Heatmap(counts_data_holdout_top_t, 
                         top_annotation = ha_holdout,
                         col = colorRamp2(c(-4, -2, -1, 0, 1, 2, 4), c("#4676b5", "#82b1d3", "#dbeff6", "#fefebd", "#fee395", "#fc9961", "#d73027")),
                         show_row_names = TRUE,
                         show_column_names = FALSE)
  
  return(heatmap_tra + heatmap_val + heatmap_hol)
}

### GENERATE HEATMAPS
generate_heatmaps(norm_counts_data_training_top, norm_counts_data_validation_top, norm_counts_data_holdout_top, annotated_features_df, meta_training, meta_validation, meta_holdout)


# Reverse the probabilities to flip the score correlation
probabilities_training_flipped <- 1 - probabilities_training
probabilities_validation_flipped <- 1 - probabilities_validation
probabilities_holdout_flipped <- 1 - probabilities_holdout

# Load required libraries
library(ggplot2)
library(patchwork)

# Define a function to create boxplots for each dataset and metadata combination
create_boxplots <- function(data, meta_data, set_name) {
  # Merge the scores with the metadata
  merged_data <- data.frame(score = as.numeric(data), meta_data)
  
  # Create individual boxplots for each metadata column
  p1 <- ggplot(merged_data, aes(x = condition, y = score)) + 
    geom_boxplot() + 
    ggtitle(paste(set_name, "by Condition"))
  
  p2 <- ggplot(merged_data, aes(x = primary_site, y = score)) + 
    geom_boxplot() + 
    ggtitle(paste(set_name, "by Primary Site"))
  
  p3 <- ggplot(merged_data, aes(x = primary_present, y = score)) + 
    geom_boxplot() + 
    ggtitle(paste(set_name, "by Primary Present"))
  
  p4 <- ggplot(merged_data, aes(x = batch, y = score)) + 
    geom_boxplot() + 
    ggtitle(paste(set_name, "by Batch"))
  
  # Create combined boxplots for condition + primary_site and condition + primary_present
  merged_data$condition_primary_site <- interaction(merged_data$condition, merged_data$primary_site)
  merged_data$condition_primary_present <- interaction(merged_data$condition, merged_data$primary_present)
  
  p5 <- ggplot(merged_data, aes(x = condition_primary_site, y = score)) + 
    geom_boxplot() + 
    ggtitle(paste(set_name, "by Condition and Primary Site"))
  
  p6 <- ggplot(merged_data, aes(x = condition_primary_present, y = score)) + 
    geom_boxplot() + 
    ggtitle(paste(set_name, "by Condition and Primary Present"))
  
  # Combine the plots into a grid
  combined_plot <- p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 2)
  
  return(combined_plot)
}

# Assuming you have the scores in 'probabilities_training', 'probabilities_validation', and 'probabilities_holdout'
# Create boxplots for training set
training_boxplots <- create_boxplots(probabilities_training_flipped, meta_training, "Training Set")

# Create boxplots for validation set
validation_boxplots <- create_boxplots(probabilities_validation_flipped, meta_validation, "Validation Set")

# Create boxplots for holdout set
holdout_boxplots <- create_boxplots(probabilities_holdout_flipped, meta_holdout, "Holdout Set")

# Print the plots
print(training_boxplots)
print(validation_boxplots)
print(holdout_boxplots)
