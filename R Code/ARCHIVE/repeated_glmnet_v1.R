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
library(foreach)
#library(progress)
library(glmnet)
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

source("~/5hmC-Pathway-Analysis/R Code/repeated_glmnet_funcs_v1.R")
handlers("txtprogressbar")
#handlers(global = TRUE)

#closeAllConnections()

### CONFIGURATION

## General Settings
## General Settings
setwd("~/5hmC-Pathway-Analysis/")
# Location of raw_data_processing output
root_folder <- "./Output/Raw Data Processing/PMneg_PMpos_06022024_hmr_nocombat/"
filename_trunk <- "PMneg_PMpos_DESeq2"
counts_name <- paste0(root_folder, filename_trunk, "_rawcounts.csv")
normcounts_name <- paste0(root_folder, filename_trunk, "_normcounts.csv")
meta_name <- paste0(root_folder, filename_trunk, "_conditions.csv")
file_version <- "hmr_nocombat_602020_06022024_CRCHGA_seed1"
output_folder <- "./Output/repeated_holdout/"
save_final_files <- FALSE
debug_progress <- TRUE
debug_data_summary <- FALSE
random_seed <- 10
set.seed(random_seed)


## Partition Settings
training_fraction <- 60
validation_fraction <- 20
holdout_fraction <- 20
min_samples <- 0.8
min_count <- 5
create_sig_training_counts <- FALSE
combat_norm <- FALSE
interaction_params <- c("condition", "primary_present","batch") 

## DESeq2 Settings
cutoff_type = 2 # 0=padj cutoff, default; 1=lfc & pvalue cutoff; 2 = top stat genes
padj.cutoff = 0.01 # 0.1 default
pvalue.cutoff = 0.001
lfc.cutoff = 0.138 # 0.137504 ~ 10%, 0.263034 ~ 20% change, 0.32 ~ 25%, 0.415 ~ 33%, 0.58 ~ 50% change,
num_top_genes = 1000
design_formula <- ~ condition + primary_present + batch

## glmnet Settings
nfolds_val = 5
rep_val = 2
feat_occ_thresholds <- seq(0.7, 1, by = 0.1) #c(0.7,0.75,0.8,0.85,0.9,0.95,0.98,1)#
alpha_values <- seq(0.5, 1, by = 0.1)
standardize_bit = FALSE # should caret center and standardize before training models?
number_of_runs <- 10
part_rand_seed = seq(100,105, by=1)

rapid = TRUE
if(rapid){
  num_raw_genes = 100
  num_top_genes = 10
  feat_occ_thresholds <- seq(0.9, 1, by = 0.05) #c(0.7,0.75,0.8,0.85,0.9,0.95,0.98,1)#
  alpha_values <- seq(0.9, 1, by = 0.1)
  standardize_bit = FALSE # should caret center and standardize before training models?
  number_of_runs <- 10
  part_rand_seed = seq(100,102, by=1)
}

### Define Functions

### Main Script
results <- initialize_results_dfs()

# Load and preprocess data
data_prepared <- load_and_prepare_data(counts_name, normcounts_name, meta_name, debug_data_summary, debug_progress, rapid, num_raw_genes)
# Add stratum collum based on interaction parameters
data_prepared$meta$stratum <- do.call(interaction, data_prepared$meta[interaction_params])
conditions <- define_conditions(data_prepared$meta)

#plan(list(multisession, sequential))#list(tweak(multisession, workers = 10), tweak(multisession, workers = 9),tweak(multisession, workers = 4)))
plan(multisession)
# print number of cores available
print(paste("Number of cores available: ", availableCores()))


start_time <- Sys.time()

with_progress({
  handlers("txtprogressbar")
  # Initialize progressor
  p_outer <- progressr::progressor(along = part_rand_seed)
  # Initialize counter
  iteration_counter <- 0
  # Use future_lapply with progress monitoring
  results_list <- future.apply::future_lapply(part_rand_seed, function(random_seed) {
    # Increment counter
    iteration_counter <<- iteration_counter + 1
    # Update progress
    p_outer(sprintf("Split %d", random_seed), class = if (iteration_counter %% 1 == 0) "sticky", amount = 0)
    
    set.seed(random_seed)
    
    # Load and preprocess data
    partitions <- partition_data(data_prepared$meta, interaction_params, training_fraction, validation_fraction, holdout_fraction, random_seed, debug_progress)
    ordered_filtered_counts <- order_and_filter_data(partitions$training_meta, partitions$validation_meta, partitions$holdout_meta, data_prepared$counts_data, min_count, min_samples, debug_progress)
    normalized_counts <- normalize_counts(ordered_filtered_counts$training, ordered_filtered_counts$validation, ordered_filtered_counts$holdout, combat_norm, debug_progress)
    
    # Save partitioned and normalized data files if necessary
    if (save_final_files) {
      save_files(ordered_filtered_counts$training, ordered_filtered_counts$validation, ordered_filtered_counts$holdout, normalized_counts$training_normcounts, normalized_counts$validation_normcounts, normalized_counts$holdout_normcounts, conditions$first_condition, conditions$second_condition, file_version, output_folder, debug_progress)
    }
    
    # DESeq2 Analysis on Training Data
    res_table_training <- perform_deseq2(ordered_filtered_counts$training, design_formula, cutoff_type, padj_cutoff, pvalue_cutoff, lfc_cutoff, num_top_genes, debug_progress, debug_data_summary)
    
    # Prepare metadata and normalized counts for glmnet
    meta_training <- partitions$training_meta
    meta_validation <- partitions$validation_meta
    meta_holdout <- partitions$holdout_meta
    
    sig_results_genes <- rownames(res_table_training)
    
    norm_counts_data_training <- subset_normalized_counts(normalized_counts$training_normcounts, sig_results_genes)
    norm_counts_data_validation <- subset_normalized_counts(normalized_counts$validation_normcounts, sig_results_genes)
    norm_counts_data_holdout <- subset_normalized_counts(normalized_counts$holdout_normcounts, sig_results_genes)
    
    class1_name <- meta_training[2, 1]
    class2_name <- meta_training[nrow(meta_training), 1]
    conditions_vector_training <- create_conditions_vector(meta_training, class1_name, class2_name)
    conditions_vector_validation <- create_conditions_vector(meta_validation, class1_name, class2_name)
    conditions_vector_holdout <- create_conditions_vector(meta_holdout, class1_name, class2_name)
    
    conditions_vector_training_fac <- as.factor(make.names(conditions_vector_training))
    conditions_vector_validation_fac <- as.factor(make.names(conditions_vector_validation))
    conditions_vector_holdout_fac <- as.factor(make.names(conditions_vector_holdout))
    
    class_weights <- ifelse(conditions_vector_training == 1, 1 / sum(meta_training$condition == class1_name), 1 / sum(meta_training$condition == class2_name))
    
    # Train glmnet model
    result <- train_glmnet_model(
      norm_counts_data_training = norm_counts_data_training,
      conditions_vector_training = conditions_vector_training,
      alpha_values = alpha_values,
      number_of_runs = number_of_runs,
      nfolds_val = nfolds_val,
      standardize_bit = standardize_bit,
      class_weights = class_weights,
      meta_training = meta_training,
      debug_progress = TRUE,  # Enable debug statements
      norm_counts_data_validation = norm_counts_data_validation,
      conditions_vector_validation = conditions_vector_validation,
      conditions_vector_training_fac = conditions_vector_training_fac,
      conditions_vector_validation_fac = conditions_vector_validation_fac,
      feat_occ_thresholds = feat_occ_thresholds,
      num_top_genes = num_top_genes,
      rep_val = rep_val,
      root_folder = root_folder
    )
    
    # Evaluate the best model
    evaluation_results <- evaluate_model(
      best_model = result$best_model,
      best_params = result$best_params,
      norm_counts_data_training = norm_counts_data_training,
      norm_counts_data_validation = norm_counts_data_validation,
      norm_counts_data_holdout = norm_counts_data_holdout,
      conditions_vector_training = conditions_vector_training,
      conditions_vector_validation = conditions_vector_validation,
      conditions_vector_holdout = conditions_vector_holdout
    )
    
    # Store metrics and coefficients
    # Add new columns to best_model_metrics directly
    result$best_model_metrics$model_id <- random_seed
    result$best_model_metrics$random_seed <- random_seed
    result$best_model_metrics$auc_hol <- evaluation_results$auc_hol
    result$best_model_metrics$sensitivity_hol <- evaluation_results$sensitivity_hol
    result$best_model_metrics$specificity_hol <- evaluation_results$specificity_hol
    result$best_model_metrics$recall_hol <- evaluation_results$recall_hol
    result$best_model_metrics$precision_hol <- evaluation_results$precision_hol
    result$best_model_metrics$accuracy_hol <- evaluation_results$accuracy_hol
    result$best_model_metrics$f1_hol <- evaluation_results$f1_hol
    
    # Assign modified best_model_metrics to metrics_row
    metrics_row <- result$best_model_metrics
    
    coefficients <- as.data.frame(as.matrix(evaluation_results$best_coefficients))
    coefficients <- coefficients[coefficients$s1 != 0, , drop = FALSE]
    coefficients_df <- data.frame(feature = rownames(coefficients), coefficient = coefficients$s1, model_id = random_seed)
    
    p_outer(message = sprintf("Processed Splits: %d/%d", iteration_counter, length(part_rand_seed)))
    return(list(metrics_row = metrics_row, coefficients = coefficients_df, roc_data_tra = evaluation_results$roc_perf_tra, roc_data_val = evaluation_results$roc_perf_val, roc_data_hol = evaluation_results$roc_perf_hol, best_model=result$best_model))
  }, future.packages = c("caret", "dplyr", "limma", "DESeq2", "sva", "doParallel", "foreach", "progress", "glmnet", "ROCR", "future.apply", "BiocParallel"), future.seed = TRUE)
})

end_time <- Sys.time()
end_time - start_time

### Prepare results for visualization
metrics_df <- do.call(rbind, lapply(results_list, function(x) x$metrics_row))

roc_data_tra_list <- lapply(results_list, function(x) x$roc_data_tra)
roc_data_val_list <- lapply(results_list, function(x) x$roc_data_val)
roc_data_hol_list <- lapply(results_list, function(x) x$roc_data_hol)

# Identify the best and median models based on validation set
best_model_idx <- which.max(metrics_df$auc_val)
median_model_idx <- which.min(abs(metrics_df$auc_val - median(metrics_df$auc_val)))

#print metrics row from best model and median model
print(metrics_df[best_model_idx, ])
print(metrics_df[median_model_idx, ])

# Pull coefficients from the median model
coefficients_df <- do.call(rbind, lapply(results_list, function(x) x$coefficients))

# Drop all rows where the feature is "(Intercept)" and save to new dataframe coefficients_df_features
coefficients_df_features <- coefficients_df[coefficients_df$feature != "(Intercept)", ]

# Add gene symbols to coefficients_df_features
coefficients_df_features <- annotate_features(coefficients_df_features)

# Calculate feature appearance percentage
feature_counts <- coefficients_df_features %>%
  group_by(feature, GeneSymbol) %>%
  summarise(appearance_count = n()) %>%
  mutate(percentage = (appearance_count / length(part_rand_seed)) * 100)

### Generate plots
## ROC Plots
training_roc_plot <- plot_roc_curves(roc_data_tra_list, median_model_idx, "Training ROC Curves")
validation_roc_plot <- plot_roc_curves(roc_data_val_list, median_model_idx, "Validation ROC Curves")
holdout_roc_plot <- plot_roc_curves(roc_data_hol_list, median_model_idx, "Holdout ROC Curves")
# Arrange the plots as subplots
combined_roc_plot <- subplot(training_roc_plot, validation_roc_plot, holdout_roc_plot, nrows = 1, titleX = TRUE, titleY = TRUE)

## Coefficient plots
# Get the features in the median model
median_model_features <- results_list[[median_model_idx]]$coefficients$feature
# Drop rows where the feature column contains the word "intercept" anywhere and print warning
intercept_rows <- grepl("intercept", median_model_features, ignore.case = TRUE)
median_model_features <- median_model_features[!intercept_rows]
median_model_features <- factor(median_model_features)

box_plot <- create_boxplot(coefficients_df_features, median_model_features, feature_counts, percentage_cutoff = 20)
bar_plot <- create_vertical_bar_plot(feature_counts, median_model_features, 20)

## Display plots
box_plot
bar_plot
combined_roc_plot

# Print confusion matrices, this needs more work
#print_confusion_matrix(results_list, median_model_idx, norm_counts_data_training, conditions_vector_training, "Training")
#print_confusion_matrix(results_list, median_model_idx, norm_counts_data_validation, conditions_vector_validation, "Validation")
#print_confusion_matrix(results_list, median_model_idx, norm_counts_data_holdout, conditions_vector_holdout, "Holdout")

### Save key results and files to output folder
# Define variable to post-pend to output files that incorporates key parameters
interaction_string <- paste(substr(interaction_params, 1, 4), collapse = "-")

file_version <- paste0("_num_runs_", number_of_runs, "_nfolds_", nfolds_val, "_num_top_genes_", num_top_genes, "_seeds_", length(part_rand_seed),"_params",interaction_string)
# Create output folder if it does not exist, based on file version and date
output_folder_ver <- paste0(output_folder, file_version, "_", format(Sys.Date(), "%Y%m%d"), "_", format(Sys.time(), "%H%M%S"))# if it doesn't exist, create it
if (!dir.exists(output_folder_ver)) {
  dir.create(output_folder_ver)
}
# Save key results and files to output folder
write.csv(metrics_df, file = paste0(output_folder_ver, "/metrics_", file_version, ".csv"), row.names = FALSE)
write.csv(coefficients_df, file = paste0(output_folder_ver, "/coefficients_", file_version, ".csv"), row.names = FALSE)
saveRDS(results_list[[best_model_idx]]$best_model, file = paste0(output_folder_ver, "/best_model_", file_version, ".rds"))
saveRDS(results_list[[median_model_idx]]$best_model, file = paste0(output_folder_ver, "/median_model_", file_version, ".rds"))
saveRDS(training_roc_plot, file = paste0(output_folder_ver, "/training_roc_plot_", file_version, ".rds"))
saveRDS(validation_roc_plot, file = paste0(output_folder_ver, "/validation_roc_plot_", file_version, ".rds"))
saveRDS(holdout_roc_plot, file = paste0(output_folder_ver, "/holdout_roc_plot_", file_version, ".rds"))
saveRDS(box_plot, file = paste0(output_folder_ver, "/box_plot_", file_version, ".rds"))
saveRDS(bar_plot, file = paste0(output_folder_ver, "/bar_plot_", file_version, ".rds"))
htmlwidgets::saveWidget(training_roc_plot, paste0(output_folder_ver, "/training_roc_curves_", file_version, ".html"))
htmlwidgets::saveWidget(validation_roc_plot, paste0(output_folder_ver, "/validation_roc_curves_", file_version, ".html"))
htmlwidgets::saveWidget(holdout_roc_plot, paste0(output_folder_ver, "/holdout_roc_curves_", file_version, ".html"))
### Save all plots as SVG files
# Save all plots as SVG files using kaleido
plotly::save_image(training_roc_plot, file = paste0(output_folder_ver, "/training_roc_curves_", file_version, ".svg"))
plotly::save_image(validation_roc_plot, file = paste0(output_folder_ver, "/validation_roc_curves_", file_version, ".svg"))
plotly::save_image(holdout_roc_plot, file = paste0(output_folder_ver, "/holdout_roc_curves_", file_version, ".svg"))
#plotly::save_image(bar_plot, file = paste0(output_folder_ver, "/bar_plot_", file_version, ".svg"))
#plotly::save_image(box_plot, file = paste0(output_folder_ver, "/box_plot_", file_version, ".svg"))

# Dump all configurations and parameters, along with date to config.txt
config_content <- capture.output({
  cat("Configuration and Parameters\n")
  cat("Counts file: ", counts_name, "\n")
  cat("Normalized counts file: ", normcounts_name, "\n")
  cat("Metadata file: ", meta_name, "\n")
  cat("Training fraction: ", training_fraction, "\n")
  cat("Validation fraction: ", validation_fraction, "\n")
  cat("Holdout fraction: ", holdout_fraction, "\n")
  cat("Design formula: ", as.character(design_formula), "\n")
  cat("Cutoff type: ", cutoff_type, "\n")
  cat("Number of top genes: ", num_top_genes, "\n")
  cat("Alpha values: ", alpha_values, "\n")
  cat("Number of runs: ", number_of_runs, "\n")
  cat("Number of folds for validation: ", nfolds_val, "\n")
  cat("Feature occurrence thresholds: ", feat_occ_thresholds, "\n")
  cat("Random seed: ", random_seed, "\n")
  cat("Output folder: ", output_folder, "\n")
  cat("File version: ", file_version, "\n")
  cat("Output folder version: ", output_folder_ver, "\n")
  cat("Start time: ", start_time, "\n")
  cat("End time: ", end_time, "\n")
})

# Write the captured output to the file
writeLines(config_content, con = paste0(output_folder_ver, "/config.txt"))

### TO DO:
# Fix export of bar and boxplots as svg
# Output confusion matrices