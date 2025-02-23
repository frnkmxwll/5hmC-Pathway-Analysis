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
library(nestedcv)

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
save_final_files <- TRUE
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
rep_val = 3
feat_occ_thresholds <- seq(0.7, 1, by = 0.05) #c(0.7,0.75,0.8,0.85,0.9,0.95,0.98,1)#
alpha_values <- seq(0, 1, by = 0.1)
standardize_bit = FALSE # should caret center and standardize before training models?
number_of_runs <- 50
part_rand_seed = seq(1,80, by=1)
plan(multisession) #sequential, multisession, multicore
n_inner_folds = 5
n_outer_folds = 5

rapid = FALSE
if(rapid){
  plan(sequential)
  num_raw_genes = 1000
  num_top_genes = 100
  feat_occ_thresholds <- seq(0.95, 1, by = 0.05) #c(0.7,0.75,0.8,0.85,0.9,0.95,0.98,1)#
  alpha_values <- seq(0.9, 1, by = 0.1)
  standardize_bit = FALSE # should caret center and standardize before training models?
  number_of_runs <- 10
  part_rand_seed = seq(1,3, by=1)
  save_final_files = FALSE
}

### Define Functions

### Main Script
#results <- initialize_results_dfs()

# Load and preprocess data
data_prepared <- load_and_prepare_data(counts_name, normcounts_name, meta_name, debug_data_summary, debug_progress, rapid, num_raw_genes)
# Add stratum collum based on interaction parameters
data_prepared$meta$stratum <- do.call(interaction, data_prepared$meta[interaction_params])
conditions <- define_conditions(data_prepared$meta)

#plan(list(multisession, sequential))#list(tweak(multisession, workers = 10), tweak(multisession, workers = 9),tweak(multisession, workers = 4)))

# print number of cores available
print(paste("Number of cores available: ", availableCores()))



start_time <- Sys.time()
print(paste("Started at: ",start_time))
source("~/5hmC-Pathway-Analysis/R Code/repeated_glmnet_funcs_v1.R")
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
    cat(paste("Starting on random seed:",random_seed))
    set.seed(random_seed)
    
    # Load and preprocess data
    partitions <- partition_data(data_prepared$meta, interaction_params, training_fraction, validation_fraction, holdout_fraction, random_seed, debug_progress)
    ordered_filtered_counts <- order_and_filter_data(partitions$training_meta, partitions$validation_meta, partitions$holdout_meta, data_prepared$counts_data, min_count, min_samples, debug_progress)
    normalized_counts <- normalize_counts(ordered_filtered_counts$training, ordered_filtered_counts$validation, ordered_filtered_counts$holdout, combat_norm, debug_progress)

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
    result <- train_nestedcv_model(
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
      root_folder = root_folder,
      random_seed = random_seed,
      n_outer_folds = n_outer_folds,
      n_inner_folds = n_inner_folds
    )
    
    if (is.null(result)) {
      cat("nestedcv result is NULL. Skipping evaluation.\n")
      return(NULL)
    }

    p_outer(message = sprintf("Processed Splits: %d/%d", iteration_counter, length(part_rand_seed)))
    return(result)
  }, future.packages = c("caret", "dplyr", "limma", "DESeq2", "sva", "doParallel", "foreach", "progress", "glmnet", "ROCR", "future.apply", "BiocParallel"), future.seed = TRUE)
})

end_time <- Sys.time()
print(paste("Total Operation Time: ",end_time - start_time))

### Prepare results for visualization
# Drop failed result from results
# Assuming results_list is your list
results_list <- lapply(results_list, function(item) {
  # Check if the item is not a function
  if (!is.function(item)) {
    return(item)
  } else {
    return(NULL)
  }
})

# Remove NULL elements from the list
results_list <- results_list[!sapply(results_list, is.null)]

## Find median model
results_df = do.call(rbind, lapply(results_list, flatten_result))
auc_tra_outters <- results_df$auc_tra_outter
# Calculate the median auc_val
custom_median <- function(x) {
  n <- length(x)
  if (n %% 2 == 1) {
    return(median(x))
  } else {
    sorted_x <- sort(x)
    middle_index <- n / 2
    return(sorted_x[middle_index + 1])
  }
}
median_auc_tra_outter <- custom_median(auc_tra_outters)

# Round values to 3 significant digits
rounded_auc_vals <- signif(auc_tra_outters, 3)
rounded_median_auc_tra_outter <- signif(median_auc_tra_outter, 3)

# Find indices of lists with the rounded median auc_val
median_model_idx <- which(rounded_auc_vals == rounded_median_auc_tra_outter)
# Extract lists with the median auc_val
median_lists <- results_list[median_model_idx]

# Find the number of features in each list with the median auc_val
num_features_vals <- sapply(median_lists, function(x) x$num_features)

# Find the index of the list with the smallest num_features among those with the median auc_val
index_min_num_features_in_median <- which.min(num_features_vals)

# Find the corresponding index in the original results_list
median_model_idx <- median_model_idx[index_min_num_features_in_median]

# Extract the list with the median auc_val and the smallest num_features
best_result <- results_list[[median_model_idx]]

best_training_innner_roc <- best_result$roc_data_tra_inner
best_training_outter_roc <- best_result$roc_data_tra_outter

### INCORPORATE COEFFS
# Pull coefficients of all models into single df
# Initialize an empty data frame for coefficients
coefficients_df <- data.frame()

# Iterate over best result for each random seed, to build coefficients table
for (item in results_list) {
  # Extract the best_coefficients
  best_coefficients <- rownames(item$coefficients)[-1]
  
  # Convert the coefficients to a data frame and add model_id
  coefficients_df_item <- data.frame(
    feature = rownames(item$coefficients)[-1],
    coefficient = item$coefficients$coef[-1],
    model_id = item$part_seed
  )
  
  # Bind the current coefficients to the overall coefficients_df
  coefficients_df <- rbind(coefficients_df, coefficients_df_item)
}


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
source("~/5hmC-Pathway-Analysis/R Code/repeated_glmnet_funcs_v1.R")
all_training_inner_roc <- sapply(results_list, function(x) x$roc_data_tra_inner) # inner
all_training_outter_roc <- sapply(results_list, function(x) x$roc_data_tra_outter) # outter
all_validation_roc <- sapply(results_list, function(x) x$roc_data_val) # outter
training_inner_roc_plot <- plot_nested_roc_curves(all_training_inner_roc, "Training ROC Curves")
training_outter_roc_plot <- plot_nested_roc_curves(all_training_outter_roc, "Validation ROC Curves")
validation_roc_plot <- plot_rocr_roc_curves(all_validation_roc, "Validation ROC Curves")

# Arrange the plots as subplots
combined_roc_plot <- subplot(training_inner_roc_plot,training_outter_roc_plot, validation_roc_plot, nrows = 1, titleX = TRUE, titleY = TRUE)

## Coefficient plots
# Get the features in the median model
median_model_features <- rownames(best_result$coefficients)[-1]
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

file_version <- paste0("REP_NSTEDCV_num_runs_", number_of_runs, "_nfolds_", nfolds_val, "_num_top_genes_", num_top_genes, "_seeds_", length(part_rand_seed),"_params",interaction_string)
# Create output folder if it does not exist, based on file version and date
output_folder_ver <- paste0(output_folder, file_version, "_", format(Sys.Date(), "%Y%m%d"), "_", format(Sys.time(), "%H%M%S"))# if it doesn't exist, create it
if (!dir.exists(output_folder_ver)) {
  dir.create(output_folder_ver)
}
# Save key results and files to output folder
write.csv(results_df, file = paste0(output_folder_ver, "/resultsdf_", file_version, ".csv"), row.names = FALSE)
write.csv(coefficients_df, file = paste0(output_folder_ver, "/coefficients_", file_version, ".csv"), row.names = FALSE)
saveRDS(best_result$repeated_fit, file = paste0(output_folder_ver, "/median_model_", file_version, ".rds"))
saveRDS(training_inner_roc_plot, file = paste0(output_folder_ver, "/training_inner_roc_plot", file_version, ".rds"))
saveRDS(training_outter_roc_plot, file = paste0(output_folder_ver, "/training_outter_roc_plot", file_version, ".rds"))
saveRDS(validation_roc_plot, file = paste0(output_folder_ver, "/validation_roc_plot", file_version, ".rds"))
htmlwidgets::saveWidget(training_inner_roc_plot, paste0(output_folder_ver, "/training_inner_roc_curves_", file_version, ".html"))
htmlwidgets::saveWidget(training_outter_roc_plot, paste0(output_folder_ver, "/training_outter_roc_curves_", file_version, ".html"))
htmlwidgets::saveWidget(validation_roc_plot, paste0(output_folder_ver, "/validation_roc_curves_", file_version, ".html"))
htmlwidgets::saveWidget(combined_roc_plot, paste0(output_folder_ver, "/combined_roc_curves_", file_version, ".html"))

### Save all plots as SVG files
# Save all plots as SVG files using kaleido
plotly::save_image(training_inner_roc_plot, file = paste0(output_folder_ver, "/training_inner_roc_curves_", file_version, ".svg"))
plotly::save_image(all_training_outter_roc, file = paste0(output_folder_ver, "/training_outter_roc_curves_", file_version, ".svg"))
plotly::save_image(validation_roc_plot, file = paste0(output_folder_ver, "/validation_roc_curves_", file_version, ".svg"))
plotly::save_image(training_inner_roc_plot, file = paste0(output_folder_ver, "/training_inner_roc_curves_", file_version, ".jpg"))
plotly::save_image(all_training_outter_roc, file = paste0(output_folder_ver, "/training_outter_roc_curves_", file_version, ".jpg"))
plotly::save_image(validation_roc_plot, file = paste0(output_folder_ver, "/validation_roc_plot_roc_curves_", file_version, ".jpg"))
#plotly::save_image(box_plot, file = paste0(output_folder_ver, "/box_plot", file_version, ".jpg"))
#plotly::save_image(bar_plot, file = paste0(output_folder_ver, "/bar_plot", file_version, ".jpg"))
#plotly::save_image(box_plot, file = paste0(output_folder_ver, "/box_plot", file_version, ".svg"))
#plotly::save_image(bar_plot, file = paste0(output_folder_ver, "/bar_plot", file_version, ".svg"))
#htmlwidgets::saveWidget(box_plot, paste0(box_plot, "/combined_roc_curves_", file_version, ".html"))
#htmlwidgets::saveWidget(bar_plot, paste0(bar_plot, "/combined_roc_curves_", file_version, ".html"))

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