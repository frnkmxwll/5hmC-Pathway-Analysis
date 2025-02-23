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

library(doParallel)

### CONFIGURATION
setwd("~/5hmC-Pathway-Analysis/")

## CRC only
sig_results_file <- "./Output/DESeq2/Results/all_results_PMnegPMpos_pval 0.001 _lfc 0.138 _ CRCHGA_combati_7030_053124_cs_seed126.txt"
root_folder <- "./Output/Randomization/PMneg_PMpos_DESeq2_hmr_combat_7030_053124_CRCHGA_cs_seed126/"
filename_trunk <- "PMneg_PMpos"

counts_name_training <- paste0(root_folder, filename_trunk, "_training_rawcounts.csv")
norm_counts_name_training <- paste0(root_folder, filename_trunk, "_training_normcounts.csv")
meta_name_training <- paste0(root_folder, filename_trunk, "_training_conditions.csv")

counts_name_validation <- paste0(root_folder, filename_trunk, "_validation_rawcounts.csv")
norm_counts_name_validation <- paste0(root_folder, filename_trunk, "_validation_normcounts.csv")
meta_name_validation <- paste0(root_folder, filename_trunk, "_validation_conditions.csv")

# Parameters
ver.text <- "CRCHGA_hmr_combati_7030_cp_seed126"
number_of_runs <- 10
pvalue_cutoff <- 0.01
log2FoldChange_cutoff <- 0.13
nfolds_val <- 10
seed_val <- 7
alpha_values <- seq(0, 1, by = 0.1)  # Define alpha values
feat_occ_thresholds <- c(0.7,0.8)#0.75,0.8,0.85,0.9,0.95,0.98)
skip_rfe <- TRUE  # Set to TRUE to skip RFE, FALSE to apply RFE

ver <- paste("pval", pvalue_cutoff, "_lfc", log2FoldChange_cutoff, "_nfold", nfolds_val, "+", ver.text)

### READ AND PREPROCESS DATA
# Load and transpose data
norm_counts_data_training <- read.csv(norm_counts_name_training, row.names = 1)
norm_counts_data_validation <- read.csv(norm_counts_name_validation, row.names = 1)

norm_counts_data_training <- t(as.matrix(norm_counts_data_training))
norm_counts_data_validation <- t(as.matrix(norm_counts_data_validation))

# Load metadata
meta_training <- read.csv(meta_name_training, row.names = 1)
meta_validation <- read.csv(meta_name_validation, row.names = 1)

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

### FILTER SIGNIFICANT GENES
sig_results <- read.table(sig_results_file, header = TRUE)
sig_results <- subset(sig_results, pvalue < pvalue_cutoff & abs(log2FoldChange) > log2FoldChange_cutoff)
sig_results <- na.omit(sig_results)
sig_results_genes <- rownames(sig_results)

# Get the intersection of significant genes and column names of count data
existing_sig_genes <- intersect(sig_results_genes, colnames(norm_counts_data_training))
existing_sig_genes <- intersect(existing_sig_genes, colnames(norm_counts_data_validation))

# Subset count data to select only the existing significant genes
norm_counts_data_training <- norm_counts_data_training[, existing_sig_genes]
norm_counts_data_validation <- norm_counts_data_validation[, existing_sig_genes]

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

# Convert class levels to valid R variable names
conditions_vector_training_fac <- as.factor(make.names(conditions_vector_training))
conditions_vector_validation_fac <- as.factor(make.names(conditions_vector_validation))

# Calculate class weights
class_weights <- ifelse(conditions_vector_training == 1, 1 / sum(conditions_vector_training == 1), 1 / sum(conditions_vector_training == 0))
class_weights_val <- ifelse(conditions_vector_validation == 1, 1 / sum(conditions_vector_validation == 1), 1 / sum(conditions_vector_validation == 0))


### SIVS FEATURE SELECTION
library(sivs)
set.seed(12345)
sivs_obj <- sivs::sivs(x = norm_counts_data_training, 
                       y = factor(conditions_vector_training), 
                       family="binomial",
                       progressbar = TRUE,
                       parallel.cores = "grace",
                       nfolds=10,
                       iter.count=100)

plot(sivs_obj)

# Define the threshold
threshold <- 90

# Remove the intercept from the selection frequency list
selection_freq_no_intercept <- sivs_obj$selection.freq[names(sivs_obj$selection.freq) != "(Intercept)"]

# Extract the features with selection frequency above the threshold
selected_features <- names(selection_freq_no_intercept[selection_freq_no_intercept >= threshold])

selected_features <- sivs::suggest(sivs_obj, strictness = 0.01)




### DEFINE OPTIMIZATION FUNCTIONS
cl <- makePSOCKcluster(5)
registerDoParallel(cl)

# Combine stratification variables into a single factor
stratification_factor <- interaction(meta_training$condition, meta_training$batch, meta_training$primary_site)

# Define a custom index for balanced k-folds based on the combined stratification factor
customIndex <- createFolds(stratification_factor, k = 5, returnTrain = TRUE)

# Define training control
train_control <- trainControl(method = "repeatedcv",
                              number=10,
                              repeats = 5,
                              summaryFunction = twoClassSummary,
                              ## Estimate class probabilities
                              classProbs = TRUE,
                              index = customIndex)

# Define the tuning grid for alpha and lambda
tune_grid <- expand.grid(
  alpha = seq(0, 1, by = 0.05),
  lambda = 10^seq(-4, -1, length = 100)
)

# Perform the grid search
set.seed(123) # For reproducibility
model <- train(
  x = norm_counts_data_training[, selected_features],
  y = conditions_vector_training_fac,
  preProc = c("center"),
  method = "glmnet",
  metric="ROC",
  trControl = train_control,
  tuneGrid = tune_grid,
  weights = class_weights,
  allowParallel = TRUE
)

# Stop the parallel cluster
stopCluster(cl)

# Inspect the results
print(model$bestTune)

# Extract the best model
best_model <- model$finalModel

# Get the optimal lambda value from the bestTune
optimal_lambda <- model$bestTune$lambda

# Extract coefficients using the optimal lambda value
coefficients <- coef(best_model, s = optimal_lambda)

# Convert the coefficients to a matrix
coefficients_matrix <- as.matrix(coefficients)

# Filter non-zero coefficients
non_zero_coefficients <- coefficients_matrix[coefficients_matrix != 0, , drop = FALSE]

# Convert to a readable format and print
non_zero_coefficients_list <- as.data.frame(non_zero_coefficients)
print(non_zero_coefficients_list)

# Make predictions on the validation set
predictions <- predict(model, newdata = norm_counts_data_validation, type = "prob",preProc = c("center"))

# Extract the predicted probabilities for the positive class
predicted_probabilities <- predictions[, "X1"]

# Calculate AUC
roc_curve <- roc(conditions_vector_validation_fac, predicted_probabilities)
auc_value <- auc(roc_curve)
print(paste("AUC:", auc_value))

# Make class predictions
predicted_classes <- predict(model, newdata = norm_counts_data_validation)

# Generate confusion matrix
conf_matrix <- confusionMatrix(predicted_classes, conditions_vector_validation_fac)
print(conf_matrix)

